"""Epitope-contact scoring for the egfrBR3 BRIDGE arm -- BOTH sites, one binder.

A bridge design is supposed to span site A and site B with a single chain (closest approach
11.46 A, max span 45.59 A; measured 2026-10-02). So unlike contacts_a21.py
(site A only) this scores EVERY draw against BOTH footprints and reports the joint term.

⛔ IDENTICAL to contacts_a21.py / contacts_all_shards.py in every respect that affects a number,
   so BR3 numbers stay comparable to the site-A and P14 rounds:
  - CUT = 5.0 A, min over ALL atom pairs
  - chain A = binder, chain B = target
  - FOOT["A"] primary = the same 11 residues; FOOT["B"] primary = FP7 as SELECTED
  - numbering ASSERTED BY RESIDUE IDENTITY, not range (the folds are against full ecto621,
    chain B spans 1-621, so a range check is the wrong test)

⛔ THREE SITE-B FOOTPRINTS ARE REPORTED AND NONE IS TAKEN AS THE GATE. The repo is genuinely
   ambiguous and the differences are large, so the choice is made in the writeup, not here:
     n_B7   FP7, AS SELECTED   19 20 21 25 29 30 50      (the site-B footprint in setup_rfd3_sphere_campaign.py)
     n_B8   FP8, canonical 8   + Q28                     (the STRETCH diagnostic)
     n_B7p  FP7', drop M30     19 20 21 25 28 29 50
   ⛔⛔ M30 IS KNOWN UNSATISFIABLE: contacted by 0 of 336 site-B designs, and
   r(distance from M30, contact rate) over the other seven = +0.967 (measured), so it was
   DROPPED. A bridge scored on any footprint CONTAINING M30 is being charged for a residue nothing
   reaches -- which is why n_B7p exists and why n_B7/n_B8 are not the whole story.

   ⇒ `full_both` is deliberately the PAIR (n_A == 11) AND (n_B7p == 7): full site A, plus site B
   on the only footprint that drops the unsatisfiable residue. `full_both_fp7` uses FP7 as
   selected, for continuity with the published 0/10,080 figure.

usage: contacts_bridge.py <land_dir> <out.csv>
"""
import glob, os, sys, csv
import numpy as np
from biotite.structure.io.pdb import PDBFile
import biotite.structure as struc

CUT = 5.0
FOOT_A = [380, 382, 384, 408, 409, 410, 411, 412, 417, 438, 465]
FOOT_B7 = [19, 20, 21, 25, 29, 30, 50]
FOOT_B8 = [19, 20, 21, 25, 28, 29, 30, 50]
FOOT_B7P = [19, 20, 21, 25, 28, 29, 50]
IDENT = {380: "PHE", 409: "HIS", 412: "PHE", 438: "ILE", 417: "VAL",
         19: "THR", 20: "PHE", 21: "GLU", 25: "LEU",
         28: "GLN", 29: "ARG", 30: "MET", 50: "TYR"}


def main(land, outp):
    files = sorted(glob.glob(os.path.join(land, "**", "BR_*_model_*.pdb"), recursive=True))
    print(f"PDBs found: {len(files)}")
    if not files:
        print("⛔ no BR_*.pdb under that land dir"); sys.exit(2)
    rows, checked, refused = [], False, []
    for i, f in enumerate(files):
        bn = os.path.basename(f)[:-4]
        stem, mid = bn.rsplit("_model_", 1)
        try:
            design, sid = stem[3:].rsplit("_id", 1)      # strip the "BR_" site prefix
        except ValueError:
            refused.append((bn, "unparseable name")); continue
        aa = PDBFile.read(f).get_structure(model=1)
        aa = aa[struc.filter_amino_acids(aa)]
        b, t = aa[aa.chain_id == "A"], aa[aa.chain_id == "B"]
        if len(b) == 0 or len(t) == 0:
            refused.append((bn, "empty chain")); continue
        if not checked:
            obs = {}
            for rid, rn in zip(t.res_id, t.res_name):
                obs.setdefault(int(rid), str(rn))
            bad = {k: (v, obs.get(k, "ABSENT")) for k, v in IDENT.items() if obs.get(k) != v}
            if bad:
                print("⛔ NUMBERING NOT WT -- refusing:", bad); sys.exit(2)
            missing = [r for r in sorted(set(FOOT_A + FOOT_B8)) if r not in obs]
            if missing:
                print(f"⛔ footprint residues absent: {missing}"); sys.exit(2)
            print(f"✅ numbering asserted by identity on {len(IDENT)} residues across BOTH sites; "
                  f"all {len(set(FOOT_A + FOOT_B8))} footprint residues present")
            checked = True
        d = np.linalg.norm(t.coord[:, None, :] - b.coord[None, :, :], axis=2)
        hit = set(int(r) for r in t.res_id[d.min(1) < CUT])
        nA = sum(1 for r in FOOT_A if r in hit)
        nB7 = sum(1 for r in FOOT_B7 if r in hit)
        nB8 = sum(1 for r in FOOT_B8 if r in hit)
        nB7p = sum(1 for r in FOOT_B7P if r in hit)
        rows.append(dict(design=design, seq=sid, model=mid,
                         n_A=nA, n_A_total=len(FOOT_A),
                         n_B7=nB7, n_B7_total=len(FOOT_B7),
                         n_B8=nB8, n_B8_total=len(FOOT_B8),
                         n_B7p=nB7p, n_B7p_total=len(FOOT_B7P),
                         full_A=int(nA == len(FOOT_A)),
                         full_both=int(nA == len(FOOT_A) and nB7p == len(FOOT_B7P)),
                         full_both_fp7=int(nA == len(FOOT_A) and nB7 == len(FOOT_B7)),
                         m30=int(30 in hit),
                         n_target_contacts=len(hit)))
        if (i + 1) % 200 == 0:
            print(f"  ...{i + 1}/{len(files)}", flush=True)

    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote {outp}: {len(rows)} draws")
    if refused:
        print(f"⛔ REFUSED {len(refused)}:")
        for bn, why in refused[:20]:
            print(f"   {bn}: {why}")

    # ---- summary. Counts are over DRAWS; the independent unit is the SEQUENCE,
    #      so per-design bests are printed too and are what a ranking should use.
    n = len(rows)
    print(f"\n=== DRAWS (n={n}) ===")
    for k, tot in (("n_A", len(FOOT_A)), ("n_B7", len(FOOT_B7)),
                   ("n_B8", len(FOOT_B8)), ("n_B7p", len(FOOT_B7P))):
        v = np.array([r[k] for r in rows])
        print(f"  {k:5s} mean {v.mean():5.2f}/{tot}  median {np.median(v):4.1f}  "
              f"max {v.max()}  full {int((v == tot).sum())} ({100.0 * (v == tot).mean():.2f}%)")
    print(f"  M30 contacted in {sum(r['m30'] for r in rows)} of {n} draws")
    print(f"  full_both      (A 11/11 AND B 7/7 FP7') : {sum(r['full_both'] for r in rows)}/{n}")
    print(f"  full_both_fp7  (A 11/11 AND B 7/7 FP7 ) : {sum(r['full_both_fp7'] for r in rows)}/{n}")

    byd = {}
    for r in rows:
        byd.setdefault(r["design"] + "_id" + r["seq"], []).append(r)
    print(f"\n=== SEQUENCES (n={len(byd)}, {n / max(len(byd), 1):.1f} draws each) ===")
    bA = sum(1 for v in byd.values() if max(x["n_A"] for x in v) == len(FOOT_A))
    bB = sum(1 for v in byd.values() if max(x["n_B7p"] for x in v) == len(FOOT_B7P))
    bb = sum(1 for v in byd.values() if any(x["full_both"] for x in v))
    print(f"  seqs with >=1 draw at full site A (11/11) : {bA}/{len(byd)}")
    print(f"  seqs with >=1 draw at full site B (FP7')  : {bB}/{len(byd)}")
    print(f"  seqs with >=1 draw spanning BOTH          : {bb}/{len(byd)}")
    top = sorted(byd.items(), key=lambda kv: -max(x["n_A"] + x["n_B7p"] for x in kv[1]))[:10]
    print("\n  top 10 by best single-draw (n_A + n_B7p):")
    for name, v in top:
        best = max(v, key=lambda x: x["n_A"] + x["n_B7p"])
        print(f"    {name:48s} A {best['n_A']:2d}/11  B {best['n_B7p']}/7  "
              f"tgt {best['n_target_contacts']:3d}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2])
