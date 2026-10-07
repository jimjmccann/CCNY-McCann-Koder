#!/usr/bin/env python3
"""BRIDGE selection screen — ONE binder, BOTH epitopes. The combined criteria, in executable form.

⭐⭐ WHY A SEPARATE SCORER AND NOT select_all31.py + select_siteB.py RUN TWICE. Those two ask
"does this binder cover ITS site", and for a bridge the answer must be "does it cover BOTH, with
one rigid body". Running them separately would also use two different WT->output maps, while a
bridge design has exactly ONE target (the 14 A union sphere, 224 res / 18 segs).

⛔⛔ THE CONJUNCTION IS A VERY LOW-PROBABILITY EVENT AND THE CUT IS GRADED FOR THAT REASON.
MEASURED on the single-site shards (shard10 each):
    site A full chain  12/336 = 3.57%      site B full chain (incl. glycan gate) 20/336 = 5.95%
    naive product if independent           0.21% = ~11 designs per 5,040-design arm
⭐ And the empirical check is harsher than the product suggests: across 672 single-site designs,
   **ZERO** achieved full coverage of their own site AND >=4 residues on the other. Site-A designs
   scored 0 on the other site 80.4% of the time; site-B designs 64.6%.
⛔ So ANDing `cov11==11` with `cov7==7` would very likely return nothing, and a screen that
   returns nothing teaches nothing. ⇒ the primary cut here is **GRADED total coverage plus BOTH
   ANCHORS**, with the all-or-nothing conjunction reported alongside as a stretch column:
     cov_total   of 18 (11 site-A + 7 site-B)        <- graded, the primary
     anchors_ok  H409 AND R29 both within 4.5 A      <- non-negotiable: a bridge that misses an
                                                        anchor is two half-interfaces, not a bridge
     cov11==11 and cov7==7                           <- the stretch goal, reported not required
⚠️ Independence is the WRONG model in both directions, so the ~11/arm figure is an order of
   magnitude to beat, not a projection: a bridge must satisfy both at once with one body (harder),
   but it is AIMED at both, so its per-site coverage is not drawn from the single-site distribution
   at all (easier). ⇒ SET THE CUT FROM THE FIRST SHARD'S OWN DISTRIBUTION, which this prints.

⛔ M30 IS EXCLUDED from the site-B footprint (not a key contact residue). It was
   reached by 0 of 336 site-B designs, and r(distance from M30, contact rate) = +0.967 across the
   other seven -- it is not part of the same surface patch.
⛔ bb_min / cterm are MEASURED, NOT GATED here, for the same reason as the single-site scorers:
   they are minima over a 224-residue target and are not comparable to the crop-era 2.0 / 5.0.
"""
import argparse, glob, gzip, csv, os, numpy as np
from biotite.structure import annotate_sse
from biotite.structure.io.pdbx import CIFFile, get_structure
import wt2out as _wt2out

CUT = 4.5
BB = {"N", "CA", "C", "O"}
FP_A = {380: "PHE", 382: "LEU", 384: "GLN", 408: "GLN", 409: "HIS", 410: "GLY",
        411: "GLN", 412: "PHE", 417: "VAL", 438: "ILE", 465: "LYS"}
FP_B = {19: "THR", 20: "PHE", 21: "GLU", 25: "LEU", 28: "GLN", 29: "ARG", 50: "TYR"}
ANCHORS = {409: ("HIS", ["ND1", "NE2"]), 29: ("ARG", ["NE", "NH1", "NH2"])}

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=".", help="dir holding shard*/ of *.cif.gz")
ap.add_argument("out", nargs="?", default="bridge_selection.csv")
ap.add_argument("--wt2out", required=True, help="the bridge bundle's WT2OUT_bridge?.tsv")
a = ap.parse_args()

wmap = _wt2out.load(a.wt2out)
print("WT->output map:", wmap.describe())
missing = [w for w in list(FP_A) + list(FP_B) if w not in wmap.w2o]
if missing:
    raise SystemExit(f"⛔ footprint residues absent from this cut: {missing}")
print(f"✅ all {len(FP_A)}+{len(FP_B)} = {len(FP_A)+len(FP_B)} footprint residues are in the map")
for w, (t, _) in ANCHORS.items():
    print(f"   anchor WT {t}{w} -> output chain B{wmap.out(w)}")
print("⚠️ bb_min / cterm MEASURED, NOT GATED (minima over a 224-residue target).")

rows, skipped, checked = [], 0, False
files = sorted(glob.glob(f"{a.root}/shard*/*.cif.gz") + glob.glob(f"{a.root}/*.cif.gz"))
print(f"designs: {len(files)}", flush=True)

for i, f in enumerate(files):
    try:
        with gzip.open(f, "rt") as fh:
            st = get_structure(CIFFile.read(fh), model=1)
        st = st[st.element != "H"]
        A = st[st.chain_id == "A"]; B = st[st.chain_id == "B"]
        ids = np.unique(A.res_id); L = len(ids)
        if not checked:
            seen, pairs = set(), []
            for rid, rn in zip(B.res_id, B.res_name):
                if int(rid) not in seen:
                    seen.add(int(rid)); pairs.append((int(rid), rn))
            n = wmap.check_target(pairs, os.path.basename(f))
            print(f"✅ map verified against the design's target: {n}/{n} types identical", flush=True)
            checked = True
        per, cov = {}, {}
        for tag, FP in (("A", FP_A), ("B", FP_B)):
            c = 0
            for w, typ in FP.items():
                o = wmap.out(w)
                s = B[B.res_id == o]
                assert len(s) and s.res_name[0] == typ, \
                    f"WT {w} -> out {o} is {s.res_name[0] if len(s) else 'ABSENT'}, want {typ}"
                d = float(np.min(np.linalg.norm(A.coord[:, None] - s.coord[None], axis=-1)))
                per[f"d{w}"] = round(d, 2)
                c += d <= CUT
            cov[tag] = c
        # both anchors, by their functional atoms
        anc = {}
        for w, (typ, atoms) in ANCHORS.items():
            s = B[(B.res_id == wmap.out(w)) & (np.isin(B.atom_name, atoms))]
            anc[w] = float(np.min(np.linalg.norm(A.coord[:, None] - s.coord[None], axis=-1))) \
                if len(s) else float("nan")
        anchors_ok = int(anc[409] <= CUT and anc[29] <= CUT)
        # nearest binder Asp carboxylate to each anchor -- now informative (no guideposts)
        asp = {}
        for w, (typ, atoms) in ANCHORS.items():
            ref = B[(B.res_id == wmap.out(w)) & (np.isin(B.atom_name, atoms))].coord
            best = float("nan")
            for r in ids:
                s = A[A.res_id == r]
                if s.res_name[0] != "ASP":
                    continue
                od = s[np.isin(s.atom_name, ["OD1", "OD2"])].coord
                if len(od):
                    d = float(np.min(np.linalg.norm(od[:, None] - ref[None], axis=-1)))
                    if not (best == best) or d < best:
                        best = d
            asp[w] = round(best, 2) if best == best else ""
        bb = float(np.min(np.linalg.norm(
            A[np.isin(A.atom_name, list(BB))].coord[:, None] - B.coord[None], axis=-1)))
        sse = annotate_sse(A); hx = 100.0 * float((sse == "a").sum()) / max(len(sse), 1)
        last = int(np.max(ids))
        ct = float(np.min(np.linalg.norm(
            A[A.res_id == last].coord[:, None] - B.coord[None], axis=-1)))
        rows.append(dict(f=os.path.basename(f), shard=f.split("/")[-2], L=L,
                         cov_total=cov["A"] + cov["B"], covA=cov["A"], covB=cov["B"],
                         anchors_ok=anchors_ok, d_H409=round(anc[409], 2), d_R29=round(anc[29], 2),
                         asp_H409=asp[409], asp_R29=asp[29],
                         full_both=int(cov["A"] == len(FP_A) and cov["B"] == len(FP_B)),
                         bb_min=round(bb, 2), helix=round(hx, 1), cterm=round(ct, 2),
                         target_nres=int(len(np.unique(B.res_id))), target_natoms=int(len(B)),
                         **per))
    except Exception as e:
        skipped += 1
        print(f"  SKIP {os.path.basename(f)}: {e}", flush=True)
    if (i + 1) % 500 == 0:
        print(f"  {i+1}/{len(files)}", flush=True)

if not rows:
    raise SystemExit("⛔ no designs scored")
with open(a.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"wrote {a.out}  n={len(rows)}  skipped={skipped}\n")

n = len(rows)
from collections import Counter
print("  cov_total (of 18) distribution:")
c = Counter(r["cov_total"] for r in rows)
for k in sorted(c, reverse=True):
    print(f"     {k:2d}  {c[k]:4d}  {100*c[k]/n:5.1f}%")
print(f"\n  both anchors within {CUT} A : {sum(r['anchors_ok'] for r in rows)}/{n} "
      f"({100*sum(r['anchors_ok'] for r in rows)/n:.1f}%)")
print(f"  FULL both sites (11 and 7)  : {sum(r['full_both'] for r in rows)}/{n} "
      f"  <- the stretch goal")
print(f"  covA==11                    : {sum(1 for r in rows if r['covA']==11)}/{n}")
print(f"  covB==7                     : {sum(1 for r in rows if r['covB']==7)}/{n}")
for col in ("bb_min", "cterm", "helix", "d_H409", "d_R29"):
    v = np.array([r[col] for r in rows], dtype=float); v = v[~np.isnan(v)]
    if len(v):
        print(f"  {col:9s} min {v.min():6.2f}  p5 {np.percentile(v,5):6.2f}  "
              f"median {np.median(v):6.2f}  max {v.max():6.2f}")
print("\n⚠️ SET THE cov_total CUT FROM THE ABOVE, not from site A's 11/11 or site B's FP7.")
