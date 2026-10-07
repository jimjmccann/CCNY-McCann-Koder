#!/usr/bin/env python3
"""SELF-CONSISTENCY for ecto621 folds: design (sphere target) vs Boltz (621-mer target).

⛔⛔ WHY A SECOND SCRIPT AND NOT A PATCH TO `selfconsistency_egfr.py`.
That script refuses these folds with "target residues differ: design 231 vs pred 621 (shared 231)".
It is RIGHT to refuse. RFd3 renumbers the sphere target **1..N**; the assay construct is numbered
**1..621 (WT)**. The integer overlap 1..231 is a COINCIDENCE: design residue 1 is WT 8 (site A),
fold residue 1 is WT 1. Relaxing the count check would superpose MISMATCHED RESIDUES and return a
confident, meaningless `pose_rmsd`. ⇒ the target leg must go through the bundle's WT->out map, and
residue IDENTITY is asserted after mapping, never just the count: an anchor-only check
cannot catch a wrong map.

TWO NUMBERS, same definitions as `selfconsistency_egfr.py` so the values read across:
  sc_rmsd    binder-onto-binder Ca RMSD (Kabsch on the binder ALONE)
             => FOLD fidelity, independent of docking.  NISE gate: < 1.5 A
             (Fry, Slaw & Polizzi, Nature 656:237, 2026 -- the first filter of their 100,000 -> 4 funnel)
  pose_rmsd  superpose on the TARGET (mapped Ca), then binder Ca RMSD with NO binder re-alignment
             => PLACEMENT fidelity; includes fold error, so pose_rmsd >= sc_rmsd always.
  => sc small + pose large = folded right, docked wrong
     sc large              = did not build the designed shape

⛔ CUT ON MEDIANS PER BACKBONE, NEVER ON COUNTS -- never rank on partial data: a truncated or
uneven set of draws reads too high and has already reversed a winner on this project once.
2026-10-03.
"""
import argparse, csv, glob, gzip, os, re, sys
import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure.io.pdb import PDBFile
import biotite.structure as struc


def kabsch_rmsd(P, Q):
    """RMSD of P onto Q after optimal superposition. Both (N,3)."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    return float(np.sqrt(((Pc @ R - Qc) ** 2).sum(1).mean())), R, P.mean(0), Q.mean(0)


def ca_map(aa, chain):
    """{res_id: (coord, resname)} for CA atoms of one chain."""
    s = aa[(aa.chain_id == chain) & (aa.atom_name == "CA")]
    return {int(r): (c, str(n)) for r, c, n in zip(s.res_id, s.coord, s.res_name)}


def load_design(path):
    with gzip.open(path, "rt") as fh:
        aa = get_structure(CIFFile.read(fh), model=1)
    aa = aa[struc.filter_amino_acids(aa)]
    return ca_map(aa, "A"), ca_map(aa, "B")


def load_pred(path):
    aa = PDBFile.read(path).get_structure(model=1)
    aa = aa[struc.filter_amino_acids(aa)]
    return ca_map(aa, "A"), ca_map(aa, "B")


def read_wt2out(path):
    """WT2OUT_site?.tsv -> {out_resnum: (wt_resnum, resname)}"""
    out = {}
    with open(path) as fh:
        # ⛔ the map ships with a 3-line '#' preamble; DictReader would take it as the header.
        lines = [ln for ln in fh if not ln.startswith("#")]
    rd = csv.DictReader(lines, delimiter="\t")
    need = {"wt_resnum", "out_resnum", "resname"}
    if not need <= set(rd.fieldnames or []):
        sys.exit(f"⛔ {path}: expected columns {sorted(need)}, got {rd.fieldnames}")
    for r in rd:
        out[int(r["out_resnum"])] = (int(r["wt_resnum"]), r["resname"].strip().upper())
    if not out:
        sys.exit(f"⛔ empty map {path}")
    return out


ap = argparse.ArgumentParser()
ap.add_argument("--design-dir", required=True)
ap.add_argument("--pred-dir", required=True)
ap.add_argument("--wt2out", required=True)
ap.add_argument("--gate", type=float, default=1.5)
ap.add_argument("--out", required=True)
ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--site", choices=("A", "B", "BR"), default=None,
                help="⭐ 'BR' added 2026-10-04 for the egfrBR3 BRIDGE arm, whose folds are named "
                     "'BR_<design>_id<N>_model_<M>'. The assert-then-strip logic below is generic, "
                     "so BR needs no new code path -- but ⛔ the BRIDGE TAKES ITS OWN MAP: pass "
                     "WT2OUT_bridge2.tsv (the 118 folded backbones came from stage_br2, 218->118). "
                     "bridge2 and bridge6 maps DIFFER (measured), and the identity assertion below "
                     "is what catches the wrong one. "
                     "⛔ REQUIRED for the p14 campaign, whose prediction names are "
                     "'<site>_<design>_id<N>_model_<M>'. Design names are NOT UNIQUE ACROSS "
                     "SITES (arm_shard01_shard01_9_model_4 exists at BOTH A and B), so the "
                     "site prefix is ASSERTED here, not silently stripped: a fold whose prefix "
                     "is not this site is REFUSED, never scored against the other site's "
                     "backbone. Omit for the ecto621 naming '<design>_s<NN>_model_<N>'.")
a = ap.parse_args()

W2O = read_wt2out(a.wt2out)
designs = {re.sub(r"\.cif\.gz$", "", os.path.basename(p)): p
           for p in glob.glob(os.path.join(a.design_dir, "*.cif.gz"))}
pdbs = sorted(glob.glob(os.path.join(a.pred_dir, "**", "*.pdb"), recursive=True))
if a.limit:
    pdbs = pdbs[:a.limit]
if not designs or not pdbs:
    sys.exit(f"⛔ designs={len(designs)} pdbs={len(pdbs)}")
print(f"{len(designs)} designs, {len(pdbs)} models, map {len(W2O)} residues", file=sys.stderr)

cache, rows, skipped = {}, [], []
checked_identity = False
for i, p in enumerate(pdbs):
    b = os.path.basename(p)[:-4]
    # ecto621 naming: <design>_s<NN>_model_<N>.  p14 naming: <site>_<design>_id<N>_model_<M>.
    stem = re.sub(r"_(?:s|id)\d+_model_\d+$", "", b)
    m = re.search(r"_(?:s|id)(\d+)_model_(\d+)$", b)
    if a.site is not None:
        # ⛔ ASSERT the site, do not merely strip it. Scoring site-A folds against site-B
        # backbones is silent and returns confident nonsense -- the 330-not-332 defect.
        if not stem.startswith(a.site + "_"):
            skipped.append((b, f"site prefix is not {a.site}_ -- REFUSED")); continue
        stem = stem[len(a.site) + 1:]
    if stem not in designs or not m:
        skipped.append((b, "no matching design")); continue
    sid, mid = m.group(1), m.group(2)
    if stem not in cache:
        cache[stem] = load_design(designs[stem])
    (dA, dB) = cache[stem]
    pA, pB = load_pred(p)

    # ---- binder leg: residue ids must correspond 1..L in both ----
    shared = sorted(set(dA) & set(pA))
    if len(shared) != len(dA) or len(shared) != len(pA):
        skipped.append((b, f"binder differs: design {len(dA)} vs pred {len(pA)}")); continue
    D = np.array([dA[r][0] for r in shared])
    P = np.array([pA[r][0] for r in shared])
    sc, _, _, _ = kabsch_rmsd(P, D)

    # ---- target leg: design out_resnum -> WT -> fold resnum, IDENTITY ASSERTED ----
    tgt_pairs, mism = [], 0
    for out_r, (coord, rn) in dB.items():
        wt = W2O.get(out_r)
        if wt is None:
            continue
        wt_r, wt_name = wt
        q = pB.get(wt_r)
        if q is None:
            continue
        if q[1] != rn:
            mism += 1
            continue
        tgt_pairs.append((coord, q[0]))
    if mism or len(tgt_pairs) < 50:
        skipped.append((b, f"target map: {len(tgt_pairs)} paired, {mism} identity mismatches")); continue
    if not checked_identity:
        print(f"✅ target map verified by residue identity: {len(tgt_pairs)}/{len(dB)} paired, "
              f"0 mismatches", file=sys.stderr)
        checked_identity = True
    TD = np.array([x[0] for x in tgt_pairs])   # design target
    TP = np.array([x[1] for x in tgt_pairs])   # pred target
    tgt_fit, R, muP, muD = kabsch_rmsd(TP, TD)  # rotation taking PRED -> DESIGN frame
    Pmoved = (P - muP) @ R + muD               # binder carried by the TARGET superposition
    pose = float(np.sqrt(((Pmoved - D) ** 2).sum(1).mean()))

    # ⛔⛔ tgt_fit IS THE TRIPWIRE FOR pose_rmsd AND MUST BE READ FIRST.
    # The design target is the crystal-derived (TETHERED) sphere; the fold target is Boltz's own
    # 621-mer. EGFR's tethered<->extended change is large, so if Boltz built a different global
    # conformation the two frames do not correspond and pose_rmsd measures THAT, not docking.
    #   tgt_fit small + pose large -> the binder really did dock elsewhere
    #   tgt_fit large              -> pose_rmsd is UNINTERPRETABLE; say so, do not rank on it
    rows.append(dict(design=stem, seq=int(sid), model=int(mid), n_ca=len(shared),
                     sc_rmsd=round(sc, 3), pose_rmsd=round(pose, 3),
                     tgt_fit_rmsd=round(tgt_fit, 3), n_target_paired=len(tgt_pairs)))
    if (i + 1) % 500 == 0:
        print(f"  {i+1}/{len(pdbs)}", file=sys.stderr)

if not rows:
    for s in skipped[:6]:
        print(f"    {s[0]}: {s[1]}", file=sys.stderr)
    sys.exit("⛔ nothing comparable")

with open(a.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

sc = np.array([r["sc_rmsd"] for r in rows])
po = np.array([r["pose_rmsd"] for r in rows])
print(f"\nwrote {a.out}  n={len(rows)}  skipped={len(skipped)}", file=sys.stderr)
if skipped:
    from collections import Counter
    for reason, k in Counter(re.sub(r"\d+", "N", s[1]) for s in skipped).most_common(4):
        print(f"  skip: {reason} x{k}", file=sys.stderr)
print(f"  sc_rmsd   min {sc.min():.2f}  p25 {np.percentile(sc,25):.2f}  median {np.median(sc):.2f}"
      f"  p75 {np.percentile(sc,75):.2f}  max {sc.max():.2f}", file=sys.stderr)
print(f"  pose_rmsd min {po.min():.2f}  median {np.median(po):.2f}  max {po.max():.2f}", file=sys.stderr)
tf = np.array([r["tgt_fit_rmsd"] for r in rows])
print(f"  ⭐ tgt_fit_rmsd min {tf.min():.2f}  median {np.median(tf):.2f}  max {tf.max():.2f}"
      f"   <- READ THIS BEFORE pose_rmsd", file=sys.stderr)
if np.median(tf) > 2.0:
    print(f"  ⛔⛔ TARGET FRAMES DO NOT CORRESPOND (median {np.median(tf):.2f} A). Boltz built a "
          f"different global conformation than the crystal-derived sphere (EGFR tethered<->extended). "
          f"pose_rmsd is UNINTERPRETABLE; sc_rmsd is UNAFFECTED (binder-only).", file=sys.stderr)
print(f"  ⛔ pose >= sc violated in {(po < sc - 1e-6).sum()} rows (must be 0)", file=sys.stderr)
print(f"  models passing sc_rmsd < {a.gate}: {(sc < a.gate).sum()}/{len(sc)} "
      f"= {100*(sc < a.gate).mean():.1f}%", file=sys.stderr)
