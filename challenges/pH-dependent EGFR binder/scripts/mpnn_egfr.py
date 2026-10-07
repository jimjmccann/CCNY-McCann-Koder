#!/usr/bin/env python3
"""Stage 3 — SolubleMPNN sequence design on RFd3 EGFR binder backbones.

THE SETTLED FLAGS, implemented here and nowhere else:
  model_type=soluble_mpnn, checkpoint `solublempnn_v_48_020.pt`
  --omit_AA C        ⛔ an `allowed_aas` list that excludes P and PERMITTS C is a trap: it yields
                     2-3 free thiols per design that never pair, with partners 12-15 A apart.
                     Omit cysteine explicitly instead.
  P PERMITTED        ⛔ do NOT omit: the design needs crisp helix caps and turns to expose backbone
                     NH/C=O, and proline is the capping residue. Omitting it fights the spec.
  --bias_AA D/E      RFd3's own sequence head makes NEAR-NEUTRAL binders: the charge gate passed
                     only 3 of 8 smoke designs and ALL 5 failures were the ACIDITY bound. When only
                     a handful of designs can be tested, culling ~62% of the pool for a property
                     that can simply be requested is not worth it. ⚠️ THE MAGNITUDE IS NOT SETTLED
                     — tune on a small batch, measure NCPR, then commit. That is what `--bias` and
                     `seq_qc_egfr.py` are for.

⛔⛔ THE TRAP THIS SCRIPT EXISTS TO AVOID: **MPNN WILL HAPPILY REDESIGN THE GUIDEPOSTS.**
The whole point of the RFd3 stage is a pinned Asp carboxylate on H433 (P1) and, from rung 1, a
pinned Tyr OH on Q435p. Those are ordinary chain-A residues in the output -- nothing marks them --
so an unguarded MPNN run silently mutates them and the entire generation step is wasted.
⇒ this script LOCATES them per design by their own contact geometry and passes `--fixed_residues`.
⭐ The model chooses the guidepost's sequence position itself (68 distinct positions over 2-70 in
the n=400 baseline), so the index CANNOT be hardcoded -- it must be found per structure.

⛔ The target (chain B) is kept in the input as CONTEXT but is never designed
(`--chains_to_design A`), so MPNN sees the interface it is designing against.

⛔⛔ TWO REQUIREMENTS FOR A SPHERE-CUT TARGET, as opposed to a contiguous crop.
 1. **NO CONSTANT OFFSET EXISTS.** An anchor index written as `409 - 369` is valid only for a
    CONTIGUOUS crop starting at WT 370. A sphere target is 15 (site A) / 26 (site B) SEGMENTS and
    RFd3 renumbers the target 1..N, so there is no constant offset: WT H409 is output B141, while
    `409-369` = B40 is a CYS. ⇒ pass `--wt2out <bundle>/WT2OUT_site?.tsv`, the map the bundle ships
    (`make_wt2out_map.py`, verified against real RFd3 output 231/231 and 230/230).
    ⛔ `--legacy-crop-offset` reproduces the old constant-offset behaviour and exists ONLY to
    re-read output from superseded crop-era runs. Do not use it on sphere output.
 2. **`--no-guideposts` FOR HOTSPOTS-ONLY RUNS.** A run carrying select_hotspots and no guideposts
    at all gives the guidepost search nothing to find, so it raises on every design. With
    `--no-guideposts` nothing is passed to `--fixed_residues` and the whole binder is designed,
    which is correct when nothing was pinned.
    ⛔ The ANCHOR TYPE IS STILL ASSERTED at the mapped index: that one check is the cheapest
    target-identity test available, and a pipeline without it cannot detect a mis-cropped target.
"""
import argparse, glob, gzip, json, os, subprocess, sys
import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure.io.pdb import PDBFile

import wt2out as _wt2out

# \u26d4\u26d4 PARTNER RESIDUES ARE NAMED IN **WT** NUMBERING AND TRANSLATED THROUGH THE BUNDLE'S MAP.
#    They used to be written as `409 - 369`, i.e. a CONSTANT offset baked into the source. That is
#    valid only for a contiguous crop; a sphere cut has no constant offset (see the module header).
# \u26d4 EVERY partner residue is STILL asserted BY TYPE after translation. Reading site B's anchor
#    at site A's index raises "chain B40 is not HIS" -- the guard working as intended, and the
#    reason this is parameterised rather than hardcoded.
LEGACY_CROP_START = {"A": 370, "B": 4}   # \u26d4 the VOIDED crops. --legacy-crop-offset only.
D_P1, D_R1 = 2.570, 2.650 # the built constraint distances, A
TOL = 0.25                # a guidepost that has MOVED is not a guidepost; refuse to guess

SITES = {
 # anchor: (WT resnum, residue type, the atoms the guidepost points at)
 # gp:     (guidepost residue type, its atoms, the BUILT distance in A)
 "A": dict(anchor=(409, "HIS", ["ND1", "NE2"]),            gp=("ASP", ["OD1", "OD2"], 2.570),
           rung1=(411, "GLN", ["OE1"], "TYR", ["OH"], 2.650)),
 # \u2b50 site B's 2.807 A is MEASURED on the landed designs (exact in 4,666 of 4,704), not assumed.
 "B": dict(anchor=( 29, "ARG", ["NE", "NH1", "NH2"]),      gp=("ASP", ["OD1", "OD2"], 2.807),
           rung1=None),
}


class _LegacyOffsetMap:
    """\u26d4 The OLD constant-offset behaviour, for re-reading the VOIDED crop campaigns only."""
    def __init__(self, site):
        self.off = LEGACY_CROP_START[site] - 1
    def out(self, wt):
        return wt - self.off
    def describe(self):
        return f"\u26d4 LEGACY CONSTANT OFFSET out = WT - {self.off} (crop-era; INVALID for a sphere cut)"
    def check_target(self, out_res, path="?"):
        # \u26d4\u26d4 A CONSTANT OFFSET CANNOT BE CHECKED: it carries no residue types, so there is
        #    nothing to compare the design's target against. That is precisely why site B's old
        #    offset passed every assertion while being wrong -- see wt2out.WT2Out.check_target.
        raise SystemExit("\u26d4\u26d4 --legacy-crop-offset cannot be verified against the design's "
                         "target: a constant offset carries no residue types. It is accepted ONLY "
                         "with --i-am-rereading-a-voided-crop-campaign, which skips the check.")


def load(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        st = get_structure(CIFFile.read(fh), model=1)
    return st[st.element != "H"]


def nearest_res(A, resname, atoms, ref, target=None):
    """chain-A residue of `resname` whose `atoms` sit at `target` A from `ref`; (resid, dist).

    \u26d4\u26d4 WITH `target` THIS PICKS THE RESIDUE CLOSEST TO THE PINNED DISTANCE, NOT THE
       NEAREST ONE, AND THE DIFFERENCE IS LOAD-BEARING. The guidepost is pinned at an exact
       distance (site A 2.570 A, site B 2.807 A), but MPNN-free RFd3 output also contains
       DESIGNED Asp residues, and one of those can sit CLOSER to the anchor than the guidepost
       does. MEASURED: site B 38 of 4,704 designs (0.81%) place a designed Asp nearer R29 than
       the pin -- one at 2.46 A against the pin's 2.807 -- and site A 1 of 10,416 (0.01%).
       A nearest-wins rule therefore returns the WRONG residue in those designs, and because
       this function's answer becomes `--fixed_residues`, the consequence is that MPNN fixes a
       decoy and REDESIGNS THE ACTUAL GUIDEPOST -- the exact failure this script exists to
       prevent, reached by the guard itself. Distance-to-pin has no such failure mode.
    """
    best = (None, 1e9, 1e9)                      # (resid, dist, score)
    for r in np.unique(A.res_id):
        s = A[A.res_id == r]
        if s.res_name[0] != resname:
            continue
        p = s[np.isin(s.atom_name, atoms)].coord
        if not len(p):
            continue
        d = float(np.linalg.norm(p[:, None, :] - ref[None, :, :], axis=-1).min())
        score = abs(d - target) if target is not None else d
        if score < best[2]:
            best = (int(r), d, score)
    return best[0], best[1]


def guideposts(st, path, want_tyr, site="A", wmap=None, want_gp=True):
    """Locate the pinned guideposts; with want_gp=False only ASSERT the anchor and fix nothing.

    ⛔ `wmap` translates WT -> output chain-B index. There is no fallback: a sphere-cut target has
       no constant offset, and guessing one reads a different residue at full confidence.
    ⭐ want_gp=False is the HOTSPOTS-ONLY arm (the whole 2026-10-02 sphere campaign). Nothing was
       pinned, so nothing may be fixed -- but the anchor type check below still runs, because that
       assertion is the cheapest proof that this design was built against the target we think.
    """
    cfg = SITES[site]
    A, T = st[st.chain_id == "A"], st[st.chain_id == "B"]
    a_wt, a_type, a_atoms = cfg["anchor"]
    gp_type, gp_atoms, gp_d = cfg["gp"]
    a_id = wmap.out(a_wt)
    anc = T[T.res_id == a_id]
    if not len(anc) or anc.res_name[0] != a_type:
        got = anc.res_name[0] if len(anc) else "ABSENT"
        raise SystemExit(f"⛔ {os.path.basename(path)}: site {site} expects {a_type} at chain "
                         f"B{a_id} (WT {a_wt}), found {got} -- wrong target, wrong site, or a map "
                         f"that does not describe this design's cut ({wmap.describe()})")
    fixed, report = [], {}
    if not want_gp:
        # ⛔ NOT a bypass of the guidepost search: there is no guidepost to find. Returning an
        #    empty --fixed_residues is the CORRECT instruction for a hotspots-only backbone.
        report["anchor_ok"] = (a_wt, a_id, a_type)
        report["guideposts"] = "NONE -- hotspots-only arm"
        return fixed, report
    r, d = nearest_res(A, gp_type, gp_atoms, anc[np.isin(anc.atom_name, a_atoms)].coord,
                       target=gp_d)
    if r is None or abs(d - gp_d) > TOL:
        raise SystemExit(f"⛔ {os.path.basename(path)}: site {site} guidepost {gp_type} not found "
                         f"at {gp_d} A (best {r}, {d:.2f} A)")
    fixed.append(f"A{r}"); report["p1_asp"] = (r, round(d, 3))
    if want_tyr:
        if cfg["rung1"] is None:
            raise SystemExit(f"⛔ site {site} has no rung-1 island -- pass --no-tyr")
        q_wt, q_type, q_atoms, t_type, t_atoms, t_d = cfg["rung1"]
        q_id = wmap.out(q_wt)
        gln = T[T.res_id == q_id]
        if not len(gln) or gln.res_name[0] != q_type:
            raise SystemExit(f"⛔ {os.path.basename(path)}: chain B{q_id} (WT {q_wt}) is not {q_type}")
        r2, d2 = nearest_res(A, t_type, t_atoms, gln[np.isin(gln.atom_name, q_atoms)].coord,
                             target=t_d)
        if r2 is None or abs(d2 - t_d) > TOL:
            raise SystemExit(f"⛔ {os.path.basename(path)}: rung-1 {t_type} not found at {t_d} A "
                             f"(best {r2}, {d2})")
        fixed.append(f"A{r2}"); report["r1_tyr"] = (r2, round(d2, 3))
    return fixed, report


LMPNN = os.environ.get("LIGANDMPNN_DIR", "/PATH/TO/LigandMPNN")
PY = os.environ.get("MPNN_PY", "/PATH/TO/mpnn-python")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", help="directory of RFd3 *.cif.gz")
    ap.add_argument("outdir")
    ap.add_argument("--omit-aa", default="C",
                    help="residues MPNN may not use. Default C (free thiols: the cell-free "
                         "lysate is not reliably oxidising). \u2b50 CH also omits HISTIDINE, which "
                         "the build spec requires at ZERO -- our own His hands the pH setpoint "
                         "back to the target. \u26d4 --omit_AA C ALONE DOES NOT DO THAT: measured "
                         "26.15% of 24,400 designs carried >=1 designed His.")
    ap.add_argument("--bias", type=float, default=0.0, help="--bias_AA magnitude on D and E")
    ap.add_argument("--bias-raw", default="", help="raw --bias_AA string, e.g. 'A:-1.0'; overrides --bias")
    ap.add_argument("--site", choices=("A", "B"), default="A",
                    help="A = switch lobe (HIS409 anchor, Asp+Tyr islands); "
                         "B = domain I (ARG29 anchor, Asp only -- use with --no-tyr)")
    ap.add_argument("--wt2out", default=None,
                    help="the bundle's WT2OUT_site?.tsv (make_wt2out_map.py). ⛔ REQUIRED for "
                         "sphere-cut targets: they are multi-segment and RFd3 renumbers 1..N, so "
                         "NO constant offset exists. Falls back to $WT2OUT_MAP / $WT2OUT_MAP_<SITE>.")
    ap.add_argument("--legacy-crop-offset", action="store_true",
                    help="⛔ use the retired constant offset (site A WT-369, site B WT-3). ONLY to "
                         "re-read the VOIDED crop campaigns; wrong for every sphere design.")
    ap.add_argument("--i-am-rereading-a-voided-crop-campaign", action="store_true",
                    dest="skip_target_check",
                    help="⛔ skip the full-target-sequence check. The ONLY legitimate use is "
                         "re-reading voided crop output with --legacy-crop-offset.")
    ap.add_argument("--no-guideposts", action="store_true",
                    help="HOTSPOTS-ONLY arms (the 2026-10-02 sphere campaign): nothing was pinned, "
                         "so fix nothing and design the whole binder. The anchor type is still "
                         "asserted. ⛔ Without this the guidepost search raises on every design.")
    ap.add_argument("--seqs", type=int, default=8, help="sequences per backbone (batch_size)")
    ap.add_argument("--temperature", type=float, default=0.1)
    ap.add_argument("--no-tyr", action="store_true", help="baseline one-island designs (no rung-1 Tyr)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--run", action="store_true", help="actually run; otherwise print the command")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    pdbdir = os.path.join(a.outdir, "pdb"); os.makedirs(pdbdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.indir, "*.cif.gz")))
    if a.limit:
        files = files[:a.limit]
    if not files:
        raise SystemExit(f"⛔ no *.cif.gz in {a.indir}")

    if a.legacy_crop_offset:
        wmap = _LegacyOffsetMap(a.site)
        print("⛔⛔ --legacy-crop-offset: using the RETIRED constant offset. This is correct ONLY "
              "for the voided contiguous crops, NEVER for a sphere cut.")
    else:
        wmap = _wt2out.load(a.wt2out, site=a.site)
    print("WT->output map:", wmap.describe())
    if a.no_guideposts:
        print("⭐ --no-guideposts: hotspots-only arm. --fixed_residues will be EMPTY (nothing was "
              "pinned, so nothing may be fixed); the anchor type is still asserted per design.")

    paths, fixed_multi, rep = {}, {}, {}
    checked = False
    for f in files:
        st = load(f)
        # ⛔⛔ PROVE THE MAP DESCRIBES THIS DESIGN'S TARGET, by EVERY residue type, not just the
        #    anchor. Site B's retired offset passes an anchor-only check by luck (all 8 footprint
        #    residues sit in the sphere's first segment). Once per run: the target is held rigid
        #    and identical in every design of a shard -- MEASURED 1799/1799 atoms in 8/8.
        if not checked and not a.skip_target_check:
            T = st[st.chain_id == "B"]
            seen, pairs = set(), []
            for rid, rn in zip(T.res_id, T.res_name):
                if int(rid) not in seen:
                    seen.add(int(rid)); pairs.append((int(rid), rn))
            n = wmap.check_target(pairs, os.path.basename(f))
            print(f"✅ map verified against the design's target: {n}/{n} residue types identical")
            checked = True
        fx, r = guideposts(st, f, want_tyr=(not a.no_tyr) and not a.no_guideposts,
                           site=a.site, wmap=wmap, want_gp=not a.no_guideposts)
        out = os.path.join(pdbdir, os.path.basename(f).replace(".cif.gz", ".pdb"))
        pf = PDBFile(); pf.set_structure(st); pf.write(out)
        paths[out] = ""; fixed_multi[out] = " ".join(fx); rep[os.path.basename(out)] = r

    pj = os.path.join(a.outdir, "pdb_path_multi.json")
    fj = os.path.join(a.outdir, "fixed_residues_multi.json")
    json.dump(paths, open(pj, "w"), indent=1)
    json.dump(fixed_multi, open(fj, "w"), indent=1)
    json.dump(rep, open(os.path.join(a.outdir, "guideposts.json"), "w"), indent=1)

    nfix = {len(v.split()) for v in fixed_multi.values()}
    if a.no_guideposts:
        print(f"⭐ {len(paths)} backbones; anchor asserted in ALL of them; NO fixed residues "
              f"(hotspots-only arm) -- the whole binder is designed.")
        assert nfix == {0}, f"⛔ --no-guideposts but fixed residues appeared: {sorted(nfix)}"
    else:
        print(f"⭐ {len(paths)} backbones; guideposts located in ALL of them; "
              f"fixed residues per design = {sorted(nfix)}")
        if len(nfix) != 1:
            print("⛔ INCONSISTENT guidepost count across designs -- inspect guideposts.json "
                  "before running")

    cmd = [PY, os.path.join(LMPNN, "run.py"),
           "--model_type", "soluble_mpnn",
           "--checkpoint_soluble_mpnn", os.path.join(LMPNN, "model_params", "solublempnn_v_48_020.pt"),
           "--pdb_path_multi", pj,
           "--chains_to_design", "A",
           "--omit_AA", a.omit_aa,
           "--out_folder", os.path.join(a.outdir, "mpnn"),
           "--batch_size", str(a.seqs),
           "--number_of_batches", "1",
           "--temperature", str(a.temperature),
           "--save_stats", "1",
           "--seed", "37"]
    # ⛔ Pass --fixed_residues_multi ONLY when something is actually fixed. A file of empty
    #    strings is not obviously harmless to LigandMPNN's parser, and "fix nothing" is already
    #    the default; the json is still written, as the record of what was found.
    if not a.no_guideposts:
        cmd[cmd.index("--pdb_path_multi") + 2:cmd.index("--pdb_path_multi") + 2] = \
            ["--fixed_residues_multi", fj]
    # ⛔ ALANINE. soluble_mpnn cannot use the tetrad's fix: run.py:68-74 FORCES
    # ligand_mpnn_use_side_chain_context=0 for every model_type except ligand_mpnn, so that knob
    # is already applied and is NOT the cause of Ala excess here. The equivalent lever for this
    # model is a NEGATIVE --bias_AA on A (plus --temperature). Pass it with --bias-raw.
    if a.bias_raw:
        cmd += ["--bias_AA", a.bias_raw]
    elif a.bias:
        cmd += ["--bias_AA", f"D:{a.bias},E:{a.bias}"]
    print("\nCOMMAND:\n  " + " ".join(cmd) + "\n")
    if not a.run:
        print("(dry run -- pass --run to execute)"); return
    r = subprocess.run(cmd, cwd=LMPNN)
    print(f"rc={r.returncode}")
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
