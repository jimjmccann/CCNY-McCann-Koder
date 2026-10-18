#!/usr/bin/env python3
"""j409_path_or_fold_20261015.py -- is the binder a FOLD, or a PATH from pin to pin? Local, no spend.

⭐⭐ WHY: the PI, 2026-10-15, looking at the v38 `.pse`: *"They still look artificial. They're mostly
helical. They're just going from path to path just like we had the problem with the c3 rfd3."*
That is the same observation he made on 2026-10-10 — *"more globules connected by a linked tether
than you really see in normal proteins"* — which is what `min_cut` was built for.
This turns the impression into a number with named controls.

**`min_cut`** = the fewest long-range CA-CA contacts (<8 Å, |i-j|>4) crossing ANY interior cut of
the chain. A single folded domain keeps many crossing everywhere; a **path** — or two lobes on a
tether — has a point almost nothing crosses.

⛔⛔ **`min_cut` IS MARGIN-SENSITIVE AND PUBLISHED REFERENCE NUMBERS ARE MARGIN-SPECIFIC**
(`j153`:40-43 records m=10→26/10, m=12→29/15, m=13→30/17, m=15→40/25, m=25→47/54). ⇒ **the controls
are RECOMPUTED here at the same margin as the designs**, never quoted from another campaign across a
different setting.

Also reported, because `min_cut` alone cannot distinguish "compact fold" from "compact blob":
  * **LR fraction** — share of CA-CA contacts (<8 Å) that are long-range (|i-j| > 12). A helical
    hairpin path is dominated by local contacts; a real domain is not.
  * **Rg / Rg_globular** — Rg against 2.2·N^0.38 Å, the empirical globular expectation. ≫1 = extended.
⛔ None of these is an affinity statement, and none is sufficient alone.

⛔ A SCORER MUST NEVER NAME A CHAIN for our designs — roles are
DETECTED. Chains ARE named for the reference controls, deliberately, and each is printed.

usage: j409_path_or_fold_20261015.py --out <fresh.tsv> '<glob>' ['<glob>' ...] [--margin 13]
"""

import argparse
import glob as globmod
import os
import statistics
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from j393_continuity_gate_20261015 import classify_chains, read_cif_atoms

CUT_A, SEP, LR_SEP = 8.0, 4, 12

# ⛔ Chains NAMED on purpose: reference structures, deliberately selected. (label, path, chain)
CONTROLS = [
    ("LCB1 de novo helical minibinder", "status/competition2/data/J130_inputs_20261010_7jzu.pdb", "A"),
    ("SARS-CoV-2 RBD (natural)", "status/competition2/data/J130_inputs_20261010_7jzu.pdb", "B"),
    ("TNFa protomer (natural, our target)",
     "status/competition2/rfd3_tnf/1tnf_repaired_D143_20261007.pdb", "A"),
]


def metrics(X, margin):
    n = len(X)
    if n < 2 * margin + 2:
        return None
    d = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    ii, jj = np.triu_indices(n, k=SEP + 1)
    keep = d[ii, jj] < CUT_A
    ii, jj = ii[keep], jj[keep]
    cuts = np.asarray([int(((ii < k) & (jj >= k)).sum()) for k in range(margin, n - margin)])
    n_contacts = len(ii)
    n_lr = int((np.abs(jj - ii) > LR_SEP).sum())
    rg = float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))
    return dict(n_res=n, min_cut=int(cuts.min()), med_cut=int(np.median(cuts)),
                lr_frac=(n_lr / n_contacts if n_contacts else 0.0),
                rg=rg, rg_ratio=rg / (2.2 * n ** 0.38))


def ca_from_pdb(path, chain):
    X = []
    for L in open(path):
        if L.startswith("ATOM") and L[21] == chain and L[12:16].strip() == "CA":
            X.append([float(L[30:38]), float(L[38:46]), float(L[46:54])])
    return np.asarray(X)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--margin", type=int, default=13)
    ap.add_argument("globs", nargs="+")
    a = ap.parse_args(argv)
    if os.path.exists(a.out):
        sys.exit(f"⛔ {a.out} exists -- FRESH --out per experiment")
    files = sorted({p for g in a.globs for p in globmod.glob(g)})
    if not files:
        sys.exit("⛔ no files matched -- an empty glob is NOT a result")

    print("=== j409 path-or-fold ===")
    print(f"  files {len(files)} | margin {a.margin} | contact <{CUT_A} Å, |i-j|>{SEP}, "
          f"long-range |i-j|>{LR_SEP}")
    print(f"  ⛔ controls RECOMPUTED at THIS margin, not quoted from memory.\n")

    print(f"  {'REFERENCE CONTROL (chain NAMED on purpose)':<40} {'n':>4} {'min_cut':>8} "
          f"{'med_cut':>8} {'LR frac':>8} {'Rg/glob':>8}")
    ctrl = {}
    for label, path, ch in CONTROLS:
        if not os.path.exists(path):
            print(f"  ⚠️ MISSING {path}")
            continue
        m = metrics(ca_from_pdb(path, ch), a.margin)
        ctrl[label] = m
        print(f"  {label+' ('+ch+')':<40} {m['n_res']:4d} {m['min_cut']:8d} {m['med_cut']:8d} "
              f"{m['lr_frac']:8.3f} {m['rg_ratio']:8.2f}")

    rows = []
    for p in files:
        cols, atoms = read_cif_atoms(p)
        ch = "auth_asym_id" if "auth_asym_id" in cols else "label_asym_id"
        sq = "auth_seq_id" if "auth_seq_id" in cols else "label_seq_id"
        fx = "is_motif_atom_with_fixed_seq" if "is_motif_atom_with_fixed_seq" in cols else None
        binder, target, _, _ = classify_chains(cols, atoms, ch, sq, fx)
        if len(binder) != 1:
            continue
        b = binder[0]
        ca = sorted(((int(r[sq]), (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])))
                     for r in atoms if r[ch] == b and r["label_atom_id"] == "CA"))
        m = metrics(np.asarray([c for _, c in ca]), a.margin)
        if m:
            m["design"] = os.path.basename(p).replace(".cif.gz", "")
            m["shard"] = next((t for t in m["design"].split("_") if t.startswith("shard")), "?")
            rows.append(m)

    keys = ["design", "shard", "n_res", "min_cut", "med_cut", "lr_frac", "rg", "rg_ratio"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k])
                               for k in keys) + "\n")

    print(f"\n  OUR DESIGNS, n={len(rows)}        min / median / max")
    for k in ("min_cut", "med_cut", "lr_frac", "rg_ratio"):
        v = [r[k] for r in rows]
        print(f"    {k:<12} {min(v):8.3f} / {statistics.median(v):8.3f} / {max(v):8.3f}")
    print("\n  per arm (median min_cut / median LR frac):")
    for sh in sorted({r["shard"] for r in rows}):
        g = [r for r in rows if r["shard"] == sh]
        print(f"    {sh}: {statistics.median(r['min_cut'] for r in g):6.1f} / "
              f"{statistics.median(r['lr_frac'] for r in g):.3f}   (n={len(g)})")

    lcb1 = next((v for k, v in ctrl.items() if k.startswith("LCB1")), None)
    if lcb1:
        bar = lcb1["min_cut"]
        ge = sum(1 for r in rows if r["min_cut"] >= bar)
        z = sum(1 for r in rows if r["min_cut"] == 0)
        print(f"\n  ⇒ AGAINST THE ONLY LIKE-FOR-LIKE CONTROL (LCB1, a de novo helical minibinder,"
              f" min_cut {bar} at this margin):")
        print(f"     designs at or above LCB1 : {ge}/{len(rows)} ({100.0*ge/len(rows):.1f}%)")
        print(f"     designs with min_cut == 0 (chain in two pieces): {z}/{len(rows)}")
        print(f"     LCB1 LR fraction {lcb1['lr_frac']:.3f} vs our median "
              f"{statistics.median(r['lr_frac'] for r in rows):.3f}")
    print(f"\n  ⇒ wrote {a.out}")
    print("  ⛔ A fold test is NOT an affinity test, and min_cut is necessary, "
          "not sufficient.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
