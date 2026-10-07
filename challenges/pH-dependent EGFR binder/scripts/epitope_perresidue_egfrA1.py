"""PER-RESIDUE epitope-contact breakdown for egfrA1 (site A).

NEW FILE 2026-10-05. WHY IT EXISTS: `contacts_a21.py` computes the contacted-residue set
`hit` and then records only its SIZE (`n_epitope`). So every downstream table -- and
therefore `full_epitope_frac`, the metric the whole campaign ranks on -- can say that a
draw missed the ceiling but CANNOT say WHICH of the 11 footprint residues was lost.

That matters because the pooled c2+c90 per-draw read (2026-10-05) showed the misses are
overwhelmingly `n_epitope == 10`, i.e. ONE residue short, not a collapse. "Recover one
specific residue" is a targetable design objective; "fef < 1" is not. This script
recovers that identity from the structures already on disk. Read-only, no new compute.

⛔ The contact geometry here is COPIED from contacts_a21.py, not reinvented: same CUT,
same PRIMARY list, same chain convention (A = binder, B = target), same
`min over binder atoms < CUT` test, and the same numbering assertion by residue identity.
If that script's convention changes, this one must change with it or the two will
disagree silently. ⛔ Do NOT "improve" the cutoff here -- a different number makes these
counts incomparable with every table the campaign has already produced.

usage: epitope_perresidue_egfrA1.py <land_dir> [<land_dir> ...]
       (each <land_dir> is searched recursively for A_*_model_*.pdb)
"""
import glob
import os
import sys
import collections

import numpy as np
from biotite.structure.io.pdb import PDBFile
import biotite.structure as struc

# ⛔ COPIED from contacts_a21.py -- keep in lockstep, do not retune.
CUT = 5.0
PRIMARY = [380, 382, 384, 408, 409, 410, 411, 412, 417, 438, 465]
IDENT = {380: "PHE", 409: "HIS", 412: "PHE", 438: "ILE", 417: "VAL"}


def scan(land, checked_flag):
    files = sorted(glob.glob(os.path.join(land, "**", "A_*_model_*.pdb"), recursive=True))
    rows = []
    for f in files:
        bn = os.path.basename(f)[:-4]
        stem, mid = bn.rsplit("_model_", 1)
        design, sid = stem[2:].rsplit("_id", 1)
        aa = PDBFile.read(f).get_structure(model=1)
        aa = aa[struc.filter_amino_acids(aa)]
        b, t = aa[aa.chain_id == "A"], aa[aa.chain_id == "B"]
        if len(b) == 0 or len(t) == 0:
            continue
        if not checked_flag[0]:
            obs = {}
            for rid, rn in zip(t.res_id, t.res_name):
                obs.setdefault(int(rid), str(rn))
            bad = {k: (v, obs.get(k, "ABSENT")) for k, v in IDENT.items() if obs.get(k) != v}
            if bad:
                sys.exit(f"⛔ NUMBERING NOT WT -- refusing: {bad}")
            missing = [r for r in PRIMARY if r not in obs]
            if missing:
                sys.exit(f"⛔ footprint residues absent: {missing}")
            print(f"✅ numbering asserted by identity on {len(IDENT)} residues "
                  f"(same assertion as contacts_a21.py)")
            checked_flag[0] = True
        d = np.linalg.norm(t.coord[:, None, :] - b.coord[None, :, :], axis=2)
        hit = set(int(r) for r in t.res_id[d.min(1) < CUT])
        rows.append((design, sid, mid, [r for r in PRIMARY if r in hit]))
    return rows


def main(lands):
    checked = [False]
    rows = []
    for land in lands:
        got = scan(land, checked)
        print(f"  {land}: {len(got)} structures")
        rows += got
    if not rows:
        sys.exit("⛔ no structures found -- check the land dir(s)")

    n = len(rows)
    print(f"\ntotal structures: {n}")

    # How often is each footprint residue contacted?
    seen = collections.Counter()
    for _, _, _, got in rows:
        for r in got:
            seen[r] += 1
    print("\nPER-RESIDUE contact frequency over ALL draws "
          "(the column fef cannot show):")
    for r in PRIMARY:
        frac = seen[r] / n
        print(f"   {r:>4}  {seen[r]:>4}/{n}  {frac*100:5.1f}%  {'#' * int(round(frac * 50))}")

    # The decisive view: among draws that MISS the ceiling, which residue is absent?
    misses = [(d, s, m, got) for d, s, m, got in rows if len(got) < len(PRIMARY)]
    print(f"\ndraws BELOW the ceiling: {len(misses)} of {n} "
          f"({100*len(misses)/n:.1f}%)")
    if misses:
        lost = collections.Counter()
        for _, _, _, got in misses:
            for r in PRIMARY:
                if r not in got:
                    lost[r] += 1
        print("\n⭐ WHICH residue is LOST, counted only over those misses:")
        for r, c in lost.most_common():
            print(f"   {r:>4}  lost in {c:>3} of {len(misses)} misses  "
                  f"({100*c/len(misses):5.1f}%)  {'#' * int(round(50*c/len(misses)))}")
        bylen = collections.Counter(len(got) for _, _, _, got in misses)
        print("\n   how far short those misses fall:")
        for k in sorted(bylen, reverse=True):
            print(f"     n_epitope {k:>2}: {bylen[k]:>3} draws")
        # a miss of exactly one residue is the targetable case
        one = [(d, s, m, got) for d, s, m, got in misses
               if len(got) == len(PRIMARY) - 1]
        if one:
            lost1 = collections.Counter()
            for _, _, _, got in one:
                for r in PRIMARY:
                    if r not in got:
                        lost1[r] += 1
            print(f"\n⭐⭐ ONE-RESIDUE misses only ({len(one)} draws) -- the targetable set:")
            for r, c in lost1.most_common():
                print(f"     {r:>4}  {c:>3} of {len(one)}  ({100*c/len(one):5.1f}%)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
