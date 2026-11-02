#!/usr/bin/env python3
"""j895 -- build the RFd3 input PDB + spec for each core option (5XOX-scaffolded TfR1).

⛔ NOTHING HERE LAUNCHES OR SPENDS. It writes inputs and specs and asserts them. Validation
   through the library's own code path is a separate step (DesignInputSpecification.build).

WHY A NEW BUILDER AND NOT THE SINGLE-RESIDUE GUIDEPOST BUILDER
------------------------------
That builder is hardwired to SINGLE-RESIDUE indexed guideposts (one chain each) and asserts the binder
lands in a 45-90 aa band. These arms hold MULTI-RESIDUE secondary-structure segments and run
120-160 aa. ⛔ it is left untouched as the working reference instance; its conventions are
reused here deliberately and named at each point.

THE CONSTRUCT -- AND HOW IT DIFFERS FROM THE EARLIER CORE TRANSPLANT
--------------------------------------------------------------------
In the earlier core-transplant attempt the core was never meant to be connected: the chain was
asked to fold AROUND it and the core was to be backfilled afterwards, so it was present but was
not part of the generated backbone path.

⇒ That core was a SPATIAL OBSTACLE delivered through `unindex`, deliberately OFF the backbone
path, with the contig never naming it. ⛔ Its 0/144 continuity result therefore does not bear on
this construct. Here the held elements ARE the backbone path: every one is named IN THE CONTIG, in
sequence order, as an ordinary motif segment. That is what RFd3 is built for.

⭐ What DOES carry over from that work: ⛔⛔ a "did the held part land?" number is not a result.
Held coordinates read ~100% in place BY CONSTRUCTION, so element placement RMSD measures
constraint compliance, not design. ⛔ Never report it as a success metric.

THREE CONSTRUCT DECISIONS MADE HERE, STATED BECAUSE THEY ARE CHOICES
--------------------------------------------------------------------
D1. **Held elements are INDEXED motif segments at full backbone**, named in the contig in
    sequence order. ⛔ Not `unindex`, ⛔ not CA/CB-only -- the two bundled changes the earlier
    unindexed CA/CB-only control arm confounded (see ../methods/mining_negative_results.md).
D2. **Element SEQUENCE is released** (`select_unfixed_sequence`), so SolubleMPNN still designs
    the scaffold surface and a natural protein's sequence is not imported wholesale. What the
    construct keeps is the BACKBONE of the chosen strands and helices, not their sequence.
    ⚠️ The alternative -- holding 5XOX's own sequence to buy foldability -- is NOT taken here.
    ⚠️⚠️ CONSEQUENCE, stated because D3 below can read as if it survived this: a pin that sits
    INSIDE a held element is inside a sequence-released segment, so its identity is NOT held
    downstream. Only the bare guideposts are identity-constrained end to end.
D3. **A pin inside a held element keeps the SCAFFOLD's residue identity** (5XOX gives LYS at
    32/38/180/211 and HIS at 34/155), not the design's ARG. Reason: the constellation was matched
    on TYPED functional groups -- j851 treats an ARG guanidinium and a LYS NZ as the same CATION
    -- so the match that selected 5XOX already allowed either. ⚠️ It is still a real change at
    the cation positions: a guanidinium is planar and bidentate, a LYS NZ is a point charge.
    ⛔ Recorded per arm in the provenance sidecar.
    A pin OUTSIDE every held element stays a single-residue guidepost from the DESIGN input, at
    its spec position, with its own residue identity and `select_fixed_atoms` -- that builder's channel,
    unchanged.

CONVENTIONS INHERITED FROM THE EARLIER SPEC BUILDER, each load-bearing
----------------------------------------------------------------------
* target block `,/0,A193-298,/0,A332-379` with EACH `/0` its own comma-separated token.
  ⛔ `A6-157/0` registers NO break and silently changes the output layout. Segments are read from
  the design input's own chain A runs, never typed in.
* `ori_token` = the guidepost centroid (branch 1 of `set_com`, highest precedence), ⛔ not
  `infer_ori_strategy="com"` which is 94.8% target and carries no site information.
* `select_fixed_atoms` per residue type: HIS `ND1,NE2,CG` / ARG `NE,NH1,NH2` / LYS `NZ`.
* gap sizing is HELICAL, because our designs come out all-helical and a helix advances ~1.5 A per
  residue against an extended chain's 3.3: `gap_min = ceil(d/RISE_HELIX)`, asserted to be at or
  above the extended floor `ceil(d/RISE_EXT)`. ⇒ a RANGE, never an exact gap.
* the resolved config goes to a `_provenance.json` SIDECAR, never into the spec (the model is
  extra="forbid").

Local CPU only; no GPU and no accelerator.
"""
import argparse
import itertools
import json
import math
import os
import re

import gemmi
import numpy as np

STRICT = re.compile(r"([A-Za-z])(\d+)(?:-(\d+))?\Z")
TARGET_CHAIN = "A"
CORE_CHAIN = "H"            # all held 5XOX elements, native numbering
RISE_EXT = 3.3
RISE_HELIX = 1.5
SLACK = 3
FIXED_ATOMS = {"HIS": "ND1,NE2,CG", "ARG": "NE,NH1,NH2", "LYS": "NZ"}
ARG_AT = ("NE", "NH1", "NH2")
HIS_AT = ("ND1", "NE2", "CG")


class Check:
    def __init__(self):
        self.n, self.fails = 0, []

    def __call__(self, ok, label, detail=""):
        self.n += 1
        print("  %s  %-58s %s" % ("PASS" if ok else "FAIL", label, detail))
        if not ok:
            self.fails.append(label)

    def ok(self):
        print("== %d checks, %d FAIL ==" % (self.n, len(self.fails)))
        for f in self.fails:
            print("   FAILED: %s" % f)
        return not self.fails


def centroid(resname, at):
    if resname == "ARG":
        return np.mean([at[a] for a in ARG_AT], axis=0)
    if resname == "LYS":
        return np.asarray(at["NZ"], dtype=float)
    if resname == "HIS":
        return np.mean([at[a] for a in HIS_AT], axis=0)
    raise ValueError(resname)


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T, Pc, Qc


def rot_about_axis(axis, ang):
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * (K @ K)


def runs(nums):
    out = []
    for n in sorted(set(nums)):
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return [tuple(r) for r in out]


def read_design(path):
    st = gemmi.read_structure(path)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    chains = {}
    for ch in st[0]:
        rl = []
        for r in ch:
            at = {a.name: np.array([a.pos.x, a.pos.y, a.pos.z]) for a in r}
            rl.append((r.name, r.seqid.num, at))
        chains[ch.name] = rl
    return chains


def read_scaffold(path, chain):
    if path.endswith(".gz"):
        import gzip
        with gzip.open(path, "rt") as fh:
            doc = gemmi.cif.read_string(fh.read())
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        st = gemmi.read_structure(path)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    for ch in st[0]:
        if ch.name != chain:
            continue
        out = {}
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info is None or not info.is_amino_acid():
                continue
            at = {a.name: np.array([a.pos.x, a.pos.y, a.pos.z]) for a in r}
            if "CA" in at:
                out[r.seqid.num] = (r.name, at)
        return out
    raise SystemExit("chain %r not in %s" % (chain, path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--options", required=True, help="j894 options JSON")
    ap.add_argument("--scaffold", required=True)
    ap.add_argument("--scaffold-chain", default="B")
    ap.add_argument("--design-pdb", required=True)
    ap.add_argument("--spec", required=True, help="the existing arm spec, for select_fixed_atoms")
    ap.add_argument("--hits", required=True, action="append")
    ap.add_argument("--hit-chain", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rise-helix", type=float, default=RISE_HELIX)
    ap.add_argument("--rise-ext", type=float, default=RISE_EXT)
    ap.add_argument("--slack", type=int, default=SLACK)
    a = ap.parse_args()

    ck = Check()
    os.makedirs(a.outdir, exist_ok=True)
    opts = json.load(open(a.options))
    ELEMENTS = {int(k): v for k, v in opts["elements"].items()}

    # ---- pin order READ from the hit list, never typed in (the j894 bug)
    row = None
    for t in a.hits:
        for ln in list(open(t))[1:]:
            f = ln.rstrip("\n").split("\t")
            if f[1] == a.hit_chain:
                row = f
    if row is None:
        raise SystemExit("chain %r not in any --hits file" % a.hit_chain)
    PIN_RESI = [int(v) for v in row[4:10]]
    PIN_TYPE = row[10:16]

    spec0 = json.load(open(a.spec))
    spec0 = spec0[list(spec0)[0]]
    sel0 = spec0["select_fixed_atoms"]
    PIN_LABEL = sorted(sel0)                       # B1..G1 = the query order
    qtype = ["IMID" if "ND1" in sel0[k] else "CATION" for k in PIN_LABEL]
    ck(len(PIN_LABEL) == len(PIN_RESI) == 6, "six guideposts and six pins",
       "%s <- %s" % (PIN_LABEL, PIN_RESI))

    design = read_design(a.design_pdb)
    scaf = read_scaffold(a.scaffold, a.scaffold_chain)

    # guidepost functional centroids and their residues, from the DESIGN input
    gp = {}
    qcen = []
    for i, key in enumerate(PIN_LABEL):
        want = [x.strip() for x in sel0[key].split(",")]
        hit = None
        for rn, num, at in design.get(key[0], []):
            if all(w in at for w in want):
                hit = (rn, num, at)
        if hit is None:
            raise SystemExit("guidepost %s not found in design input" % key)
        gp[key] = hit
        qcen.append(np.mean([hit[2][w] for w in want], axis=0))
    qcen = np.asarray(qcen, dtype=np.float64)

    # ---- THE GUARD j894 was missing: pin type must agree with its guidepost, same index
    bad = []
    for i, (p, ty) in enumerate(zip(PIN_RESI, PIN_TYPE)):
        rn = scaf[p][0]
        want = {"CATION": {"ARG", "LYS"}, "IMID": {"HIS"}}[qtype[i]]
        if ty != qtype[i] or rn not in want:
            bad.append("%s(%s) <- %d %s/%s" % (PIN_LABEL[i], qtype[i], p, rn, ty))
    ck(not bad, "pin/guidepost residue-type agreement (6/6)", "; ".join(bad) or "all agree")

    target = design[TARGET_CHAIN]
    tgt_runs = runs([n for _, n, _ in target])
    ck(len(tgt_runs) >= 1, "target chain A runs read from the input",
       ",".join("A%d-%d" % r for r in tgt_runs))
    ori = [round(float(v), 2) for v in qcen.mean(axis=0)]

    built = []
    for oi, o in enumerate(opts["options"], 1):
        subset = [int(x) for x in o["subset"]]
        held = sorted(set(itertools.chain.from_iterable(
            [n for n in range(ELEMENTS[e][1], ELEMENTS[e][2] + 1) if n in scaf]
            for e in subset)))
        contained = [p for p in PIN_RESI if p in held]
        ci = [PIN_RESI.index(p) for p in contained]
        P = np.asarray([centroid(*scaf[p]) for p in contained])
        Q = qcen[ci]
        R, Pc, Qc = kabsch(P, Q)
        heldxyz = np.vstack([np.array(list(scaf[n][1].values())) for n in held])
        targxyz = np.asarray([at[k] for _, _, at in target for k in at])

        def clash_of(Rm):
            X = (heldxyz - Pc) @ Rm.T + Qc
            d = np.linalg.norm(X[:, None, :] - targxyz[None, :, :], axis=-1)
            return int((d < 2.5).sum())

        mode = "kabsch"
        if len(contained) == 2:                    # one rotational DOF survives
            axis = Q[1] - Q[0]
            best = None
            for deg in range(0, 360, 5):
                Rm = rot_about_axis(axis, math.radians(deg)) @ R
                c = clash_of(Rm)
                if best is None or c < best[1]:
                    best = (Rm, c, deg)
            R, mode = best[0], "sweep%d" % best[2]
        rms = float(np.sqrt((((P - Pc) @ R.T + Qc - Q) ** 2).sum(1).mean()))
        ck(abs(rms - float(o["pin_rmsd"])) < 2e-3,
           "opt%-2d placement reproduces j894 pin fit" % oi,
           "%.4f vs %s" % (rms, o["pin_rmsd"]))

        def xf(v):
            return (v - Pc) @ R.T + Qc

        # ---- segments in SEQUENCE order: held elements + bare guideposts
        segs = []
        for e in subset:
            lo, hi = ELEMENTS[e][1], ELEMENTS[e][2]
            segs.append(dict(key=lo, kind="elem", elem=e, lo=lo, hi=hi,
                             token="%s%d-%d" % (CORE_CHAIN, lo, hi), n=hi - lo + 1))
        for i, p in enumerate(PIN_RESI):
            if p not in held:
                lab = PIN_LABEL[i]
                segs.append(dict(key=p, kind="pin", pin=p, label=lab,
                                 token=lab, n=1))
        segs.sort(key=lambda s: s["key"])

        def anchor_ca(s, end):
            if s["kind"] == "pin":
                return gp[s["label"]][2]["CA"]
            return xf(scaf[s["hi"] if end else s["lo"]][1]["CA"])

        toks, gaps = [], []
        for i, s in enumerate(segs):
            toks.append(s["token"])
            if i + 1 < len(segs):
                d = float(np.linalg.norm(anchor_ca(s, True) - anchor_ca(segs[i + 1], False)))
                gmin = max(2, int(math.ceil(d / a.rise_helix)))
                floor = max(2, int(math.ceil(d / a.rise_ext)))
                if gmin < floor:
                    raise SystemExit("opt%d gap below the extended floor" % oi)
                toks.append("%d-%d" % (gmin, gmin + a.slack))
                gaps.append(dict(d=round(d, 2), lo=gmin, hi=gmin + a.slack, floor=floor))
        contig = ",".join(toks) + "," + ",".join(
            "/0,A%d-%d" % r for r in tgt_runs)
        for t in toks:
            if "-" in t and not t[0].isalpha():
                continue
            ck(bool(STRICT.match(t)), "opt%-2d token %-12s matches the strict contig dialect" % (oi, t))
        ck(contig.count("/0") == len(tgt_runs),
           "opt%-2d each /0 is its own comma token" % oi,
           "%d breaks, %d target runs" % (contig.count("/0"), len(tgt_runs)))

        lo_len = sum(s["n"] for s in segs) + sum(g["lo"] for g in gaps)
        hi_len = sum(s["n"] for s in segs) + sum(g["hi"] for g in gaps)

        # ---- the input PDB
        pdb = os.path.join(a.outdir, "j895_opt%02d_input.pdb" % oi)
        ser = 1
        with open(pdb, "w") as fh:
            for rn, num, at in target:
                for an, xyz in at.items():
                    fh.write(_atom(ser, an, rn, TARGET_CHAIN, num, xyz)); ser += 1
            for n in held:
                rn, at = scaf[n]
                for an, xyz in at.items():
                    fh.write(_atom(ser, an, rn, CORE_CHAIN, n, xf(xyz))); ser += 1
            for s in segs:
                if s["kind"] != "pin":
                    continue
                rn, num, at = gp[s["label"]]
                for an, xyz in at.items():
                    fh.write(_atom(ser, an, rn, s["label"][0], num, xyz)); ser += 1
            fh.write("END\n")

        sfa = {}
        for s in segs:
            if s["kind"] == "pin":
                rn = gp[s["label"]][0]
                sfa[s["label"]] = FIXED_ATOMS[rn]
        arm = {
            "dialect": 2,
            "input": os.path.basename(pdb),
            "contig": contig,
            "is_non_loopy": True,
            "select_fixed_atoms": sfa,
            # ⛔⛔ PER-ELEMENT RANGES, never one span. Chain H holds ONLY the held element
            # residues, so a span H<first>-<last> names residues that DO NOT EXIST in the
            # atom array and the library refuses: "[component=H41] Residue H41 not found".
            # That is the library's own guard catching it locally, for free.
            "select_unfixed_sequence": ",".join(
                "%s%d-%d" % (CORE_CHAIN, ELEMENTS[e][1], ELEMENTS[e][2])
                for e in subset),
            "ori_token": ori,
        }
        name = "tfr1_core_opt%02d" % oi
        with open(os.path.join(a.outdir, "%s_spec.json" % name), "w") as fh:
            json.dump({name: arm}, fh, indent=1)
        with open(os.path.join(a.outdir, "%s_provenance.json" % name), "w") as fh:
            json.dump(dict(option=oi, elements=subset, contained_pins=contained,
                           element_ranges={str(e): [ELEMENTS[e][1], ELEMENTS[e][2]]
                                           for e in subset},
                           element_lengths={str(e): len([n for n in range(
                               ELEMENTS[e][1], ELEMENTS[e][2] + 1) if n in scaf])
                               for e in subset},
                           bare_guideposts=[s["label"] for s in segs if s["kind"] == "pin"],
                           scaffold_identity_at_contained_pins={
                               str(p): scaf[p][0] for p in contained},
                           design_identity_at_those_positions={
                               str(p): gp[PIN_LABEL[PIN_RESI.index(p)]][0] for p in contained},
                           placement_mode=mode, pin_rmsd=round(rms, 4),
                           clash_atoms_vs_target=clash_of(R),
                           len_min=lo_len, len_max=hi_len, gaps=gaps,
                           rise_helix=a.rise_helix, rise_ext=a.rise_ext, slack=a.slack,
                           n_segments=len(segs), j894_total=o["total"],
                           control=bool(o.get("control"))), fh, indent=1)
        built.append((oi, name, len(segs), lo_len, hi_len, clash_of(R), mode, contained))
        print("  opt%-2d %-22s %2d seg  len %3d-%3d  clash %4d  pins-in-structure %s"
              % (oi, name, len(segs), lo_len, hi_len, clash_of(R), contained))

    print("\n== built %d arms into %s ==" % (len(built), a.outdir))
    print("⛔ NOT validated through the library yet, NOT staged, NOTHING launched.")
    ck.ok()


def _atom(ser, an, rn, ch, num, xyz):
    el = an[0] if not an[0].isdigit() else (an[1] if len(an) > 1 else "C")
    return ("ATOM  %5d %-4s%4s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (ser, an, rn, ch, num, xyz[0], xyz[1], xyz[2], el))


main()
