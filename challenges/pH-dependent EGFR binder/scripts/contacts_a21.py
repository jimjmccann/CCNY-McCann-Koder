"""Epitope-contact scoring for the egfrA21 21-DRAW arm (site A only).

Why a separate script: contacts_all_shards.py requires a per-shard self-consistency CSV
(sc_shard??_A.csv) that egfrA21 has never been scored for. The settled gate -- epitope
CONTACT is PRIMARY, pose_rmsd is NEVER a cut -- needs only the PDB, so this scores contacts
alone and emits no pose/sc column rather than inventing one.

⛔ IDENTICAL to contacts_all_shards.py in every respect that affects the number:
  - FOOT["A"] primary = the same 11 residues
  - CUT = 5.0 A, min over all atom pairs
  - chain A = binder, chain B = target
  - numbering ASSERTED BY RESIDUE IDENTITY (not range): the folds are against the full
    ecto621 target, so a range check is the wrong test; the numbering is asserted by residue identity instead.
  Any divergence here would make the 21-draw and 3-draw numbers incomparable, which is
  the whole point of the comparison.

⭐⭐ 2026-10-05, `decisions/0047`: this script already computed the contacted-residue set
   `hit` and then threw the IDENTITIES away, recording only its SIZE. `0047` scores a draw on
   `min(n_epitope/11, c380)`, so the F380 bit is now load-bearing and must come out of the SAME
   geometry pass -- deriving it anywhere else would let the two disagree silently.
   ⭐ TWO COLUMNS ADDED, NOTHING CHANGED: `c380` (1/0) and `epitope_hits` (the hit residues,
   ";"-joined). `n_epitope`, `n_epitope_total` and `n_target_contacts` are byte-identical to what
   this script produced before, so every landed table still reproduces.
   ⛔ `c380` is NOT a new measurement -- it is `380 in hit`, the same 5.0 A min-over-atom-pairs
   test as the other ten residues. Do not re-derive it with a different cutoff.

usage: contacts_a21.py <land_dir> <out.csv>
"""
import glob, os, sys, csv, collections
import numpy as np
from biotite.structure.io.pdb import PDBFile
import biotite.structure as struc

CUT = 5.0
PRIMARY = [380, 382, 384, 408, 409, 410, 411, 412, 417, 438, 465]
F380 = 380                      # decisions/0047 -- the floor residue, asserted PHE in IDENT below
IDENT = {380: "PHE", 409: "HIS", 412: "PHE", 438: "ILE", 417: "VAL"}


def main(land, outp):
    files = sorted(glob.glob(os.path.join(land, "**", "A_*_model_*.pdb"), recursive=True))
    print(f"PDBs found: {len(files)}")
    rows, checked, refused = [], False, []
    for f in files:
        bn = os.path.basename(f)[:-4]
        stem, mid = bn.rsplit("_model_", 1)
        design, sid = stem[2:].rsplit("_id", 1)
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
            missing = [r for r in PRIMARY if r not in obs]
            if missing:
                print(f"⛔ footprint residues absent: {missing}"); sys.exit(2)
            print(f"✅ numbering asserted by identity on {len(IDENT)} residues; "
                  f"all {len(PRIMARY)} footprint residues present")
            checked = True
        d = np.linalg.norm(t.coord[:, None, :] - b.coord[None, :, :], axis=2)
        hit = set(int(r) for r in t.res_id[d.min(1) < CUT])
        got = [r for r in PRIMARY if r in hit]
        rows.append(dict(design=design, seq=sid, model=mid,
                         n_epitope=len(got),
                         n_epitope_total=len(PRIMARY),
                         n_target_contacts=len(hit),
                         # ⛔ decisions/0047: the F380 bit, from THIS pass. See the docstring.
                         c380=int(F380 in hit),
                         epitope_hits=";".join(str(r) for r in got)))
    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"scored: {len(rows)}")
    for k, n in sorted(collections.Counter((r["design"], r["seq"]) for r in rows).items()):
        print(f"  {k[0]} id{k[1]}: {n} draws")
    if refused:
        print("⛔ REFUSED:", refused)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
