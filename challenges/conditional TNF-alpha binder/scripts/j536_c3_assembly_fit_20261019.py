#!/usr/bin/env python3
"""j536 -- does the predicted C3 binder assemble as DESIGNED? The ASSEMBLY fit, not the protomer.

⛔⛔ WHY `j521`'s PER-COPY RMSD CANNOT ANSWER THIS. `j521` superposes each predicted binder copy
on its own designed chain **independently** (`:295-304`, `kabsch_rmsd(P, D)` per copy). Kabsch is
a free rigid motion, so each copy is moved onto its design chain on its own. ⇒ that measures the
**PROTOMER FOLD** and is **completely blind to the arrangement**: three perfectly folded protomers
stacked in the wrong twist score exactly as well as three in the designed trimer. -- *check the ASSEMBLY, never the angles.*

⭐ This script superposes **all three binder chains AT ONCE** as one rigid body. The residual then
contains the arrangement, which is the thing the C3 design actually specified.

⛔⛔ CHAIN PERMUTATION IS MANDATORY HERE (`J517`). The three binder chains are IDENTICAL in
sequence and Boltz assigns chain IDs per run with no guarantee of order, so a positional A->A,
B->B, C->C pairing is a GUESS. `J517` cost 34 of 88 C1 designs to exactly that assumption. All
**6** permutations are enumerated and the best is reported, with `frame` naming the winner.
⚠️ Per-copy RMSD is permutation-INSENSITIVE (all three design chains are the same shape, and
Kabsch absorbs the position), which is why `j521` got away without it and this script cannot.

⭐ READERS AND MATH ARE IMPORTED UNDER `j521`'s OWN SHA1 PINS -- never re-implemented, so these
numbers are comparable to `j533`/`j535` by construction.

INTERPRETATION KEY, written before the run:
  assembly_rmsd ~= percopy_rmsd  => protomers ARE arranged as designed; the fold is the only error
  assembly_rmsd >> percopy_rmsd  => protomers recovered individually, ARRANGEMENT WRONG
  both large                     => the protomer fold itself is not recovered; arrangement moot
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", required=True, help="the predictions/ dir, NOT out_<arm>")
    ap.add_argument("--design-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--binder-chains", default="A,B,C")
    ap.add_argument("--design-binder-chains", default="A,B,C")
    ap.add_argument("--percopy-tsv", default="",
                    help="j521 output for the SAME arm; merges percopy_rmsd for the comparison")
    a = ap.parse_args()

    if Path(a.out).exists():
        sys.exit(f"⛔ exists, refusing to overwrite: {a.out}  (use a fresh --out for each run)")

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "j521_mod", HERE / "j521_c3_boltz_selfconsistency_20261018.py")
    j521 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(j521)
    j497, sha497 = j521.load_j497(False)
    kabsch_rmsd, tm_score, sha465 = j497.load_j465_math(False)
    print(f"  pinned readers: j497 {sha497[:8]}  j465 {sha465[:8]}")

    bch = [c for c in a.binder_chains.split(",") if c]
    dch = [c for c in a.design_binder_chains.split(",") if c]
    if len(bch) != len(dch):
        sys.exit(f"⛔ {bch} vs {dch}: chain lists must be the same length")
    print(f"RESOLVED LAYOUT  predicted binder={','.join(bch)}  designed binder={','.join(dch)}")

    percopy = {}
    if a.percopy_tsv:
        for r in csv.DictReader(open(a.percopy_tsv), delimiter="\t"):
            percopy[r["design"]] = r

    import itertools
    perms = list(itertools.permutations(range(len(bch))))
    rows, notes = [], []
    for sub in sorted(Path(a.pred_dir).iterdir()):
        if not sub.is_dir():
            continue
        pdbs = sorted(sub.glob("*_model_0.pdb"))
        if not pdbs:
            notes.append(f"{sub.name}: no *_model_0.pdb")
            continue
        name = sub.name
        tag = j521_tag(name)
        dpath = Path(a.design_dir) / f"{tag}.cif.gz"
        if not dpath.exists():
            notes.append(f"{name}: no design file {dpath}")
            continue
        try:
            pred = j497.read_atoms(pdbs[0], bch)
            des = j497.read_atoms(dpath, dch)
        except Exception as exc:                                  # noqa: BLE001
            notes.append(f"{name}: read failed -- {exc}")
            continue

        # ⭐ one rigid body: concatenate the chains in the permuted order, then ONE Kabsch.
        nres = min(min(len(pred[c]["ca"]) for c in bch),
                   min(len(des[c]["ca"]) for c in dch))
        if nres < 10:
            notes.append(f"{name}: only {nres} paired CA per chain")
            continue
        D = np.vstack([des[c]["ca"][:nres] for c in dch])
        best = None
        for pm in perms:
            P = np.vstack([pred[bch[i]]["ca"][:nres] for i in pm])
            r = float(kabsch_rmsd(P, D))
            if best is None or r < best[0]:
                best = (r, float(tm_score(P, D)), "".join(bch[i] for i in pm))
        row = dict(design=name, backbone=tag, n_res_per_chain=nres,
                   n_res_assembly=nres * len(bch),
                   assembly_rmsd=round(best[0], 3), assembly_tm=round(best[1], 4),
                   frame=best[2], n_perms=len(perms))
        pc = percopy.get(name)
        if pc:
            # ⛔ j521's `_best` is np.max and `_worst` np.min REGARDLESS of direction, so for a
            # DISTANCE `_best` is the WORST copy. Named here by what the number IS. (J535)
            row["percopy_rmsd_worstcopy"] = pc.get("rmsd_binder_copy_best", "")
            row["percopy_rmsd_bestcopy"] = pc.get("rmsd_binder_copy_worst", "")
            row["percopy_tm_bestcopy"] = pc.get("tm_binder_copy_best", "")
            try:
                row["assembly_minus_percopy"] = round(
                    best[0] - float(pc["rmsd_binder_copy_best"]), 3)
            except (KeyError, TypeError, ValueError):
                pass
        rows.append(row)

    if not rows:
        sys.exit("⛔ 0 rows -- --pred-dir is probably out_<arm> and not its predictions/ dir")
    cols = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    cyc = {"".join(bch[k:] + bch[:k]) for k in range(len(bch))}
    nnc = sum(1 for r in rows if r["frame"] not in cyc)
    print(f"wrote {a.out}  n={len(rows)}")
    print(f"  ⭐ {nnc}/{len(rows)} best fit needed a NON-CYCLIC chain pairing "
          f"-- positional A->A would have mis-scored those (J517)")
    for note in notes[:10]:
        print(f"  ⚠️ {note}")
    return 0


def j521_tag(name):
    """Design-file stem for a prediction dir name. ⛔ j521's EXACT regex (`:186`), not a lookalike:
    the arm suffix is optional there and omitting it would silently fail to find the design file."""
    import re
    return re.sub(r"_seq\d+(_bound|_apo|_trimer|_holo)?$", "", name)


if __name__ == "__main__":
    sys.exit(main())
