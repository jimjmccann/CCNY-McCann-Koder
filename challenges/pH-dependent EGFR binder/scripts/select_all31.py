"""Site-A selection screen, length-agnostic (reads every shard, so output is not length-locked).

================================================================================================
⛔⛔ THREE CHANGES, each required before this screen is valid on a sphere-cut target.

1. ⛔ **`W2O = -369` IS GONE.** It was valid only for the contiguous 106 aa crop starting at
   WT 370. The site-A sphere is **231 residues in 15 SEGMENTS** and RFd3 renumbers the target
   1..N, so no constant offset exists. ⇒ `--wt2out <bundle>/WT2OUT_siteA.tsv`, the map the bundle
   ships (`make_wt2out_map.py`). ⭐ MEASURED: WT H409 is output **B141**; `409-369` = B40 is a
   **CYS**. Site A at least failed LOUDLY on the type assertion -- unlike site B, which passed by
   luck. The map is now checked against the design's target by EVERY residue type (231/231).

2. ⛔⛔ **THE GUIDEPOST EXCISION IS DELETED, NOT BYPASSED.** The old `:21-29` took the chain-A ASP
   whose carboxylate O sat NEAREST the H409 ring N -- with no pin window at all, so it was pure
   nearest-Asp -- and then `Aex = A[A.res_id != gp]` REMOVED it from the binder before measuring
   footprint contact. These arms are **hotspots-only: there is no guidepost**, so that deleted a
   real, designed Asp from the binder. ⛔ It fires exactly on the designs that place an Asp at the
   H409 imidazole, i.e. on the interaction the arm is FOR.
   ⭐ `--guidepost-pin <A>` restores PIN-based identification (2.570 A for the P1 island) for a
   future guidepost arm. There is deliberately NO fallback: a guidepost off its pin is not one.
   ⇒ `asp_h409` / `asp_h409_resid` now REPORT the nearest binder Asp instead of excising it.

3. ⚠️ **`bb_min` / `cterm` ARE MEASURED, NOT GATED.** Both are minima over ALL of chain B, which
   went 106 -> 231 residues, so each can only FALL for a purely geometric reason. The stored
   thresholds (`bb_min >= 2.0`, `cterm >= 5.0`) were calibrated on the CROP, where C-terminal
   clearance did essentially all the filtering (78.0% pass, measured on the site-A backbone screen).
   ⛔ They are **not comparable to the stored site-A baselines**. `target_nres`/`target_natoms`
   are emitted beside them so no later comparison can be made without noticing.
   ⇒ **Set the cuts from shard00's measured distribution.**
================================================================================================
"""
import argparse, glob, gzip, csv, os, numpy as np
from biotite.structure import annotate_sse
from biotite.structure.io.pdbx import CIFFile, get_structure
import wt2out as _wt2out

CUT = 4.5; BB = {"N", "CA", "C", "O"}
ANCHOR_WT, ANCHOR_TYPE = 409, "HIS"
RING_N = ["ND1", "NE2"]
GP_TOL = 0.05
# ⭐ WT (mature EGFR) numbering -- the project's ONLY residue scheme since 2026-10-02.
FP11 = {380: "PHE", 382: "LEU", 384: "GLN", 408: "GLN", 409: "HIS", 410: "GLY",
        411: "GLN", 412: "PHE", 417: "VAL", 438: "ILE", 465: "LYS"}
FP9 = {k: v for k, v in FP11.items() if k not in (380, 382)}   # drop F380/L382

# ⛔⛔ EPITOPE MUST BE ASSIGNED, NOT ASSUMED. The two footprints are 11.5 A apart, and MEASURED
#    2026-10-02: **each sphere contains the OTHER site's footprint almost entirely** -- the site-A
#    sphere holds all 8 of site B's residues (WT 19-30 in its segment 0, Y50 in segment 4), and the
#    site-B sphere holds 6 of site A's (WT 380, 408-412, 438, 465; only Q384 is absent). So a
#    site-A-arm design CAN land on site B, and this screen would score it as a low-coverage FAILURE
#    rather than a mis-targeted design. `cov_other` + `epitope` make that visible instead.
#    ⚠️ Residues outside this cut are SKIPPED, not fatal: `cov_other` is out of `n_other`, which is
#    printed, so the two sites' cross-coverage columns are not interchangeable.
FP_OTHER = {19: "THR", 20: "PHE", 21: "GLU", 25: "LEU",
            28: "GLN", 29: "ARG", 30: "MET", 50: "TYR"}      # site B's footprint
OTHER_NAME = "B"
CROP_ERA_STACK = "bb_min>=2.0 and helix>=20 and cterm>=5.0  (crop-era; NOT comparable)"

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default="all31", help="dir holding shard*/ of *.cif.gz")
ap.add_argument("out", nargs="?", default="all31_selection.csv")
ap.add_argument("--wt2out", default=None,
                help="the bundle's WT2OUT_siteA.tsv. ⛔ REQUIRED: there is no constant offset for "
                     "a 15-segment sphere cut.")
ap.add_argument("--guidepost-pin", type=float, default=None,
                help="restore PIN-based guidepost identification at this distance in A (2.570 for "
                     "the P1 island). ⛔ No nearest-Asp fallback exists any more.")
a = ap.parse_args()

wmap = _wt2out.load(a.wt2out, site="A")
print("WT->output map:", wmap.describe())
anchor_out = wmap.out(ANCHOR_WT)
print(f"anchor: WT {ANCHOR_TYPE}{ANCHOR_WT} -> output chain B{anchor_out} "
      f"(⛔ the retired constant gave B{ANCHOR_WT - 369}, which is a CYS)")
if a.guidepost_pin is None:
    print("⭐ hotspots-only mode: NO guidepost excision. The nearest binder Asp to H409 is "
          "REPORTED (asp_h409), never removed from the binder.")
else:
    print(f"⚠️ guidepost arm: identifying the guidepost by the pin {a.guidepost_pin} "
          f"+/- {GP_TOL} A ONLY -- no fallback.")
print(f"⚠️ bb_min / cterm are MEASURED, NOT GATED. Crop-era stack, reference only: "
      f"{CROP_ERA_STACK}")
_n_other = sum(1 for k in FP_OTHER if k in wmap.w2o)
print(f"⛔ cross-site: {_n_other} of {len(FP_OTHER)} site-{OTHER_NAME} footprint residues are ALSO "
      f"in this sphere ⇒ `cov_other` out of {_n_other}, and `epitope` says which site a design "
      f"actually engages. A design with epitope=OTHER is mis-targeted, NOT a coverage failure.")

rows, skipped = [], 0
files = sorted(glob.glob(f"{a.root}/shard*/*.cif.gz"))
print(f"designs: {len(files)}", flush=True)
checked = False

for i, f in enumerate(files):
    try:
        with gzip.open(f, "rt") as fh:
            st = get_structure(CIFFile.read(fh), model=1)
        st = st[st.element != "H"]
        A = st[st.chain_id == "A"]; B = st[st.chain_id == "B"]
        ids = np.unique(A.res_id); L = len(ids)

        # ⛔⛔ ONCE PER RUN: prove the map describes THIS target, by every residue type.
        if not checked:
            seen, pairs = set(), []
            for rid, rn in zip(B.res_id, B.res_name):
                if int(rid) not in seen:
                    seen.add(int(rid)); pairs.append((int(rid), rn))
            n = wmap.check_target(pairs, os.path.basename(f))
            print(f"✅ map verified against the design's target: {n}/{n} residue types identical",
                  flush=True)
            checked = True

        anc = B[(B.res_id == anchor_out) & (np.isin(B.atom_name, RING_N))]
        assert len(anc) and anc.res_name[0] == ANCHOR_TYPE, \
            f"anchor out B{anchor_out} (WT {ANCHOR_WT}) is not {ANCHOR_TYPE}"
        ringN = anc.coord

        cands = []
        for r in ids:
            s = A[A.res_id == r]
            if s.res_name[0] != "ASP":
                continue
            od = s[np.isin(s.atom_name, ["OD1", "OD2"])].coord
            if len(od):
                cands.append((int(r),
                              float(np.min(np.linalg.norm(od[:, None] - ringN[None], axis=-1)))))
        asp_id, asp_d = min(cands, key=lambda c: c[1]) if cands else (None, float("nan"))

        gp = None
        if a.guidepost_pin is not None:
            pinned = [c for c in cands if abs(c[1] - a.guidepost_pin) <= GP_TOL]
            if pinned:
                gp = min(pinned, key=lambda c: abs(c[1] - a.guidepost_pin))[0]
        Aex = A[A.res_id != gp] if gp is not None else A

        c11 = c9 = 0
        per = {}
        for p, typ in FP11.items():
            o = wmap.out(p)
            s = B[B.res_id == o]
            assert len(s) and s.res_name[0] == typ, \
                f"WT {p} -> out {o} is {s.res_name[0] if len(s) else 'ABSENT'}, want {typ}"
            d = float(np.min(np.linalg.norm(Aex.coord[:, None] - s.coord[None], axis=-1)))
            per[f"d{p}"] = round(d, 2)
            if d <= CUT:
                c11 += 1
                if p in FP9:
                    c9 += 1

        # ⭐ cross-site coverage: assign the epitope rather than assume it (see FP_OTHER above)
        c_oth = 0
        for p_o, typ_o in FP_OTHER.items():
            o_o = wmap.w2o.get(p_o)
            if o_o is None:
                continue
            s_o = B[B.res_id == o_o]
            if not len(s_o) or s_o.res_name[0] != typ_o:
                continue
            if float(np.min(np.linalg.norm(Aex.coord[:, None] - s_o.coord[None], axis=-1))) <= CUT:
                c_oth += 1
        bb = float(np.min(np.linalg.norm(
            A[np.isin(A.atom_name, list(BB))].coord[:, None] - B.coord[None], axis=-1)))
        sse = annotate_sse(A); hx = 100.0 * float((sse == "a").sum()) / max(len(sse), 1)
        last = int(np.max(ids))
        ct = float(np.min(np.linalg.norm(
            A[A.res_id == last].coord[:, None] - B.coord[None], axis=-1)))

        rows.append(dict(f=os.path.basename(f), shard=f.split("/")[-2], L=L,
                         asp_h409_resid=asp_id, asp_h409=round(asp_d, 3), gp=gp,
                         cov11=c11, cov9=c9, cov_other=c_oth,
                         epitope=("this" if c11 > c_oth else "OTHER" if c_oth > c11 else "tie"),
                         bb_min=round(bb, 2), helix=round(hx, 1),
                         cterm=round(ct, 2),
                         target_nres=int(len(np.unique(B.res_id))), target_natoms=int(len(B)),
                         **per))
    except Exception as e:
        skipped += 1
        print(f"  SKIP {f}: {e}", flush=True)
    if (i + 1) % 1000 == 0:
        print(f"  {i+1}/{len(files)}", flush=True)

if not rows:
    raise SystemExit("⛔ no designs scored")
with open(a.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"wrote {a.out}  n={len(rows)}  skipped={skipped}")

for col in ("bb_min", "cterm", "helix", "asp_h409"):
    v = np.array([r[col] for r in rows], dtype=float); v = v[~np.isnan(v)]
    if len(v):
        print(f"  {col:9s} n={len(v)} min={v.min():.2f} p5={np.percentile(v,5):.2f} "
              f"median={np.median(v):.2f} max={v.max():.2f}")
print("⚠️ SET bb_min / cterm CUTS FROM THE ABOVE, not from the crop-era 2.0 / 5.0.")
