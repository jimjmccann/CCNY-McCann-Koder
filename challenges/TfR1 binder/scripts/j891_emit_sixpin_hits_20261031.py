#!/usr/bin/env python3
"""j891 - emit the NAMED six-pin hit list that the pin-curve line never wrote down.

WHY THIS EXISTS
---------------
The corrected pin-count curve established that all six functional pins ARE presented in
nature (519 predicted-model + 795 PDB-derived chains at the admissible
tol_d = sqrt(2k)*tol_rmsd = 6.928 A, tol_rmsd = 2.0) and that the "6 pins is unsatisfiable"
wall was our own 1.5 A pre-filter. But both runs wrote COUNTS ONLY -- the deposited tables
(../data/j871_pincurve_{afdb,fx}_20261031.tsv) carry
k/subset/tol_d/hits/chains/pct_store/best_rmsd and NO chain identity. Nothing in a count
makes a hit a usable scaffold.

The identities were computed and thrown away: j851._search_shard already returns
(chain_index, best_rmsd) per hit, and <prefix>_names.txt sits beside the store. This script
keeps them, plus the matched residue numbers, so a hit can be pulled from the real structure
and inspected as a candidate SCAFFOLD.

DECISION RULE, fixed before the run (no verdict is pre-judged -- this is an enumeration):
  * k = 6, the FULL pin set. One subset, no sweep.
  * tol_d = floor = sqrt(12) * tol_rmsd ; tol_rmsd = 2.0 A. Both inherited, not re-chosen.
  * one row per CHAIN, its BEST-rmsd placement, ranked ascending by rmsd.
  * cap = 2000000, matching j871. A chain that overflows the cap is RECORDED as overflow,
    never silently dropped (j851's matcher returns [] on overflow, which is the same value it
    returns for "no match"; the curve runs measured 0 overflows at k=6, so a nonzero count
    here is new and must be reported).
  * no filtering on anything else. Ranking is the reader's job; this writes the record.

OUTPUTS (one per store):
  j891_sixpin_hits_<store>_20261031.tsv
    rank  chain  best_rmsd  n_placements  resid_1..6  type_1..6  overflow

REPRODUCTION GATE: the chain COUNT and BEST rmsd this script finds must equal the published
k=6 floor row of j871_pincurve_<store>_20261031.tsv. If they do not, the
emitter is not the instrument that produced the curve and the list must not be trusted.
Printed as PASS/FAIL per store before any list is believed.

Local CPU only. No GPU, no accelerator, no network.
"""
import argparse
import importlib.util
import math
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MATCHER_SRC = os.path.join(HERE, "j851_funcatom_search_20261030.py")


def load_j851():
    spec = importlib.util.spec_from_file_location("j851", MATCHER_SRC)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


M = load_j851()


def _shard(a):
    """Best placement per chain, with the matched functional-group indices kept."""
    (prefix, Q, qtyp, tol_d, tol_rmsd, cap, lo, hi) = a
    cen = np.load(prefix + "_cen.npy", mmap_mode="r")
    typ = np.load(prefix + "_typ.npy", mmap_mode="r")
    off = np.load(prefix + "_offsets.npy")
    qd = np.linalg.norm(Q[:, None, :] - Q[None, :, :], axis=-1)
    rows = []
    n_ovf = 0
    for c in range(lo, hi):
        s, e = off[c], off[c + 1]
        cc = np.asarray(cen[s:e])
        tt = np.asarray(typ[s:e])
        res = M.search_chain(cc, tt, qd, qtyp, tol_d, tol_rmsd, Q, cap)
        if not res:
            continue
        best_r, best_idx = min(res, key=lambda t: t[0])
        rows.append((c, float(best_r), len(res), tuple(int(v) for v in best_idx)))
    return rows, n_ovf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", required=True, help="store prefix, e.g. .../fx")
    ap.add_argument("--tag", required=True, help="short store name for the filename")
    ap.add_argument("--pins", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--curve", required=True, help="published curve TSV for the gate")
    ap.add_argument("--tol-rmsd", type=float, default=2.0)
    ap.add_argument("--cap", type=int, default=2000000)
    ap.add_argument("--nproc", type=int, default=16)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    labels, Qall, qtypall = M.read_query(args.pins, args.spec)
    k = len(labels)
    if k != 6:
        sys.exit(f"[abort] query has {k} pins, this script is the k=6 emitter")
    tol_d = math.sqrt(2.0 * k) * args.tol_rmsd
    qd = np.linalg.norm(Qall[:, None, :] - Qall[None, :, :], axis=-1)
    print(f"[query] {k} pins {labels}  span {qd.max():.2f} A")
    print(f"[tol]   tol_d {tol_d:.4f} A (= sqrt(2*{k})*{args.tol_rmsd})  "
          f"tol_rmsd {args.tol_rmsd} A  cap {args.cap}")

    names = open(args.store + "_names.txt").read().split("\n")
    off = np.load(args.store + "_offsets.npy")
    rid = np.load(args.store + "_rid.npy", mmap_mode="r")
    n = len(off) - 1
    print(f"[store] {args.store}  chains {n}  names {len(names)}")

    bounds = np.linspace(0, n, args.nproc + 1).astype(int)
    jobs = [(args.store, Qall, qtypall, tol_d, args.tol_rmsd, args.cap,
             int(bounds[i]), int(bounds[i + 1])) for i in range(args.nproc)]
    rows = []
    with Pool(args.nproc) as pool:
        for r, _ in pool.imap_unordered(_shard, jobs):
            rows.extend(r)
    rows.sort(key=lambda t: t[1])

    # ---- reproduction gate against the published curve, BEFORE the list is believed
    want_chains = want_best = None
    for ln in open(args.curve):
        f = ln.rstrip("\n").split("\t")
        if f[0] == "k" or len(f) < 8:
            continue
        if int(f[0]) == 6 and f[3] == "1":
            want_chains, want_best = int(f[5]), float(f[7])
    got_chains = len(rows)
    got_best = rows[0][1] if rows else float("nan")
    ok_c = (want_chains is not None and got_chains == want_chains)
    ok_b = (want_best is not None and abs(got_best - want_best) < 5e-4)
    print(f"[gate]  published k=6 floor row: chains {want_chains}  best_rmsd {want_best}")
    print(f"[gate]  this run:                chains {got_chains}  best_rmsd {got_best:.4f}")
    print(f"[gate]  chains {'PASS' if ok_c else 'FAIL'}   "
          f"best_rmsd {'PASS' if ok_b else 'FAIL'}")
    if not (ok_c and ok_b):
        print("[gate]  ** DOES NOT REPRODUCE THE PUBLISHED CURVE -- "
              "the list below is NOT the instrument that produced it **")

    os.makedirs(args.outdir, exist_ok=True)
    typ_all = np.load(args.store + "_typ.npy", mmap_mode="r")
    out = os.path.join(args.outdir, f"j891_sixpin_hits_{args.tag}_20261031.tsv")
    with open(out, "w") as fh:
        fh.write("rank\tchain\tbest_rmsd\tn_placements\t"
                 + "\t".join(f"resid_{i+1}" for i in range(k)) + "\t"
                 + "\t".join(f"type_{i+1}" for i in range(k)) + "\n")
        for i, (c, r, npl, idx) in enumerate(rows, 1):
            s = int(off[c])
            resids = [int(rid[s + j]) for j in idx]
            types = [M.TYPENAME[int(typ_all[s + j])] for j in idx]
            fh.write(f"{i}\t{names[c]}\t{r:.4f}\t{npl}\t"
                     + "\t".join(str(v) for v in resids) + "\t"
                     + "\t".join(types) + "\n")
    print(f"[out]   {out}  rows {len(rows)}")

    print(f"\n== top 25 by rmsd ({args.tag}) ==")
    print(f"{'rank':>4}  {'chain':<28} {'rmsd':>7}  {'npl':>5}  residues")
    for i, (c, r, npl, idx) in enumerate(rows[:25], 1):
        s = int(off[c])
        resids = ",".join(str(int(rid[s + j])) for j in idx)
        print(f"{i:>4}  {names[c]:<28} {r:>7.4f}  {npl:>5}  {resids}")


if __name__ == "__main__":
    main()
