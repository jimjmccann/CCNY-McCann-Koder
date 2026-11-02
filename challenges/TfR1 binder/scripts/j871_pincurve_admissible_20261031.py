#!/usr/bin/env python3
"""j871 - re-derive the pin-count curve at an ADMISSIBLE tol_d.

The decision rule, the floor derivation and the five predictions below were all fixed and
recorded BEFORE this file existed. See ../methods/structural_mining.md.

Why: an independent check showed j851's tol_d = 1.5 A pre-filter sits below the admissible
floor sqrt(2k)*tol_rmsd (6.928 A at k=6, r=2.0), so the whole published curve -- every row,
not just k=6 -- is a LOWER BOUND, and "6 pins is unsatisfiable" may be our own artefact.

Two things this driver does that j851 cannot:
  1. sweeps PIN SUBSETS (j851's CLI has no subset flag, and no 42-subset driver was ever
     committed), and
  2. COUNTS CAP OVERFLOWS. j851's matcher returns [] when a candidate set exceeds --cap,
     which is the same value as "no match" -- a silent drop that gets MORE likely exactly
     as tol_d is raised. A re-derived curve without an overflow count is not evidence.

j851 itself is NOT edited: it is imported and stays the reference instance. The vendored
search_chain_instr below is j851's own code plus an overflow counter, and --selfcheck
proves the two agree cell-for-cell before any curve is reported.
"""
import argparse
import importlib.util
import itertools
import math
import os
import time
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


# --------------------------------------------------------------- instrumented matcher
# j851.search_chain verbatim, with (a) an overflow flag instead of a silent [] and
# (b) counts-only accumulation so a high-tol_d cell cannot exhaust memory.
# Returns (n_hits, best_rmsd_or_None, overflowed).

def search_chain_instr(cen, typ, qd, qtyp, tol_d, tol_rmsd, Q, cap):
    L = len(cen)
    k = len(qd)
    if L < k:
        return 0, None, False
    D = np.sqrt(((cen[:, None, :] - cen[None, :, :]) ** 2).sum(-1))
    tok = [np.flatnonzero(typ == qtyp[i]) for i in range(k)]
    if any(len(t) == 0 for t in tok):
        return 0, None, False

    shells, counts = {}, {}
    for a, b in itertools.combinations(range(k), 2):
        Mm = np.abs(D - qd[a, b]) <= tol_d
        np.fill_diagonal(Mm, False)
        Ma = np.zeros_like(Mm)
        Ma[np.ix_(tok[a], tok[b])] = Mm[np.ix_(tok[a], tok[b])]
        if not Ma.any():
            return 0, None, False
        shells[(a, b)] = Ma
        shells[(b, a)] = Ma.T
        counts[(a, b)] = counts[(b, a)] = int(Ma.sum())

    first = min(itertools.combinations(range(k), 2), key=lambda p: counts[p])
    order = [first[0], first[1]]
    rest = [i for i in range(k) if i not in order]
    while rest:
        nxt = min(rest, key=lambda c: sum(counts[(c, p)] for p in order))
        order.append(nxt)
        rest.remove(nxt)

    A = np.argwhere(shells[(order[0], order[1])])
    if len(A) == 0:
        return 0, None, False
    if len(A) > cap:
        return 0, None, True                      # <-- j851 returned [] here
    for t in range(2, k):
        iq = order[t]
        mask = np.ones((L, len(A)), dtype=bool)
        for pos, p in enumerate(order[:t]):
            mask &= shells[(iq, p)][:, A[:, pos]]
        for pos in range(t):
            mask[A[:, pos], np.arange(len(A))] = False
        rows, cols = np.nonzero(mask)
        if len(rows) == 0:
            return 0, None, False
        if len(rows) > cap:
            return 0, None, True                  # <-- j851 returned [] here
        A = np.concatenate([A[cols], rows[:, None]], axis=1)

    inv = np.argsort(np.asarray(order))
    A = A[:, inv]
    nh = 0
    best = None
    for s in range(0, len(A), 20000):
        blk = A[s:s + 20000]
        r = M.kabsch_rmsd_batch(cen[blk].astype(np.float64), Q)
        ok = r[r <= tol_rmsd]
        nh += int(ok.size)
        if ok.size:
            b = float(ok.min())
            best = b if best is None else min(best, b)
    return nh, best, False


# --------------------------------------------------------------- cells

def floor_tol_d(k, tol_rmsd):
    """Admissible pre-filter floor: max |delta d_ij| = sqrt(2k)*r."""
    return math.sqrt(2.0 * k) * tol_rmsd


def build_cells(nq, kmin, kmax, tol_rmsd, tol_ds):
    cells = []
    for k in range(kmin, kmax + 1):
        for sub in itertools.combinations(range(nq), k):
            for td in tol_ds:
                t = floor_tol_d(k, tol_rmsd) if td == "floor" else float(td)
                cells.append((sub, round(t, 4), td == "floor"))
    return cells


def _shard(a):
    (prefix, cells, Qall, qtypall, tol_rmsd, cap, lo, hi) = a
    cen = np.load(prefix + "_cen.npy", mmap_mode="r")
    typ = np.load(prefix + "_typ.npy", mmap_mode="r")
    off = np.load(prefix + "_offsets.npy")
    acc = {i: [0, 0, None, 0] for i in range(len(cells))}   # hits, chains, best, ovf
    pre = []
    for sub, td, _ in cells:
        Q = Qall[list(sub)]
        qt = qtypall[list(sub)]
        qd = np.sqrt(((Q[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
        pre.append((Q, qt, qd, td))
    for c in range(lo, hi):
        s, e = off[c], off[c + 1]
        cc = np.asarray(cen[s:e])
        tt = np.asarray(typ[s:e])
        for i, (Q, qt, qd, td) in enumerate(pre):
            nh, best, ovf = search_chain_instr(cc, tt, qd, qt, td, tol_rmsd, Q, cap)
            r = acc[i]
            if ovf:
                r[3] += 1
            if nh:
                r[0] += nh
                r[1] += 1
                r[2] = best if r[2] is None else min(r[2], best)
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", required=True)
    ap.add_argument("--pins", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol-rmsd", type=float, default=2.0)
    ap.add_argument("--tol-d", nargs="+", default=["1.5", "floor"])
    ap.add_argument("--kmin", type=int, default=3)
    ap.add_argument("--kmax", type=int, default=6)
    ap.add_argument("--cap", type=int, default=2000000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--nproc", type=int, default=16)
    ap.add_argument("--selfcheck", type=int, default=0,
                    help="chains to cross-check vendored vs j851.search_chain, then exit")
    args = ap.parse_args()

    if os.path.exists(args.out):
        raise SystemExit(f"[abort] --out exists, refusing to append: {args.out}")

    labels, Qall, qtypall = M.read_query(args.pins, args.spec)
    qd_full = np.sqrt(((Qall[:, None, :] - Qall[None, :, :]) ** 2).sum(-1))
    print(f"[query] {len(labels)} pins {labels}  span {qd_full.max():.2f} A")
    print("[floor] admissible tol_d = sqrt(2k)*tol_rmsd:")
    for k in range(args.kmin, args.kmax + 1):
        print(f"    k={k}  {floor_tol_d(k, args.tol_rmsd):.3f} A")

    off = np.load(args.store + "_offsets.npy")
    n = len(off) - 1
    hi = n if not args.limit else min(n, args.limit)

    # ---- faithfulness gate: vendored copy must equal j851's own matcher
    if args.selfcheck:
        cen = np.load(args.store + "_cen.npy", mmap_mode="r")
        typ = np.load(args.store + "_typ.npy", mmap_mode="r")
        bad = 0
        for sub in list(itertools.combinations(range(len(labels)), 5))[:3] + \
                   [tuple(range(len(labels)))]:
            Q = Qall[list(sub)]; qt = qtypall[list(sub)]
            qd = np.sqrt(((Q[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
            for td in (1.5, floor_tol_d(len(sub), args.tol_rmsd)):
                for c in range(min(args.selfcheck, n)):
                    s, e = off[c], off[c + 1]
                    cc = np.asarray(cen[s:e]); tt = np.asarray(typ[s:e])
                    ref = M.search_chain(cc, tt, qd, qt, td, args.tol_rmsd, Q, args.cap)
                    nh, best, ovf = search_chain_instr(cc, tt, qd, qt, td,
                                                       args.tol_rmsd, Q, args.cap)
                    if ovf:
                        continue          # the one case where they legitimately differ
                    if len(ref) != nh:
                        bad += 1
                        print(f"  MISMATCH sub={sub} td={td:.3f} chain={c} "
                              f"j851={len(ref)} vendored={nh}")
                    elif ref and abs(min(r for r, _ in ref) - best) > 1e-9:
                        bad += 1
                        print(f"  RMSD MISMATCH sub={sub} td={td:.3f} chain={c}")
        print(f"[selfcheck] {'PASS' if not bad else 'FAIL'}  mismatches {bad}")
        return

    cells = build_cells(len(labels), args.kmin, args.kmax, args.tol_rmsd, args.tol_d)
    print(f"[cells] {len(cells)} = subsets x tol_d   store chains {n}, scanning {hi}",
          flush=True)

    bounds = np.linspace(0, hi, args.nproc + 1).astype(int)
    jobs = [(args.store, cells, Qall, qtypall, args.tol_rmsd, args.cap,
             int(bounds[i]), int(bounds[i + 1])) for i in range(args.nproc)]
    t0 = time.time()
    tot = {i: [0, 0, None, 0] for i in range(len(cells))}
    with Pool(args.nproc) as pool:
        for acc in pool.imap_unordered(_shard, jobs):
            for i, (h, ch, b, o) in acc.items():
                t = tot[i]
                t[0] += h; t[1] += ch; t[3] += o
                if b is not None:
                    t[2] = b if t[2] is None else min(t[2], b)
    el = time.time() - t0

    with open(args.out, "w") as fh:
        fh.write("k\tsubset\ttol_d\tis_floor\thits\tchains\tpct_store\t"
                 "best_rmsd\tcap_overflow_chains\n")
        for i, (sub, td, isf) in enumerate(cells):
            h, ch, b, o = tot[i]
            fh.write(f"{len(sub)}\t{','.join(labels[j].split(':')[0] for j in sub)}\t"
                     f"{td:.4f}\t{int(isf)}\t{h}\t{ch}\t{100.0 * ch / max(1, hi):.4f}\t"
                     f"{'' if b is None else f'{b:.4f}'}\t{o}\n")
    print(f"[done] {el:.0f}s  wrote {args.out}")

    print("\n== best subset per (k, tol_d), ranked by chains ==")
    for k in range(args.kmin, args.kmax + 1):
        for isf in (0, 1):
            sel = [(i, c) for i, c in enumerate(cells)
                   if len(c[0]) == k and int(c[2]) == isf]
            if not sel:
                continue
            i, c = max(sel, key=lambda x: tot[x[0]][1])
            h, ch, b, o = tot[i]
            tag = "FLOOR" if isf else "1.5  "
            print(f"  k={k} {tag} tol_d {c[1]:6.3f}  best "
                  f"{','.join(str(j) for j in c[0]):>12}  hits {h:>9}  chains {ch:>6}"
                  f"  ({100.0 * ch / max(1, hi):6.2f}%)  best_rmsd "
                  f"{'--' if b is None else f'{b:.3f}'}  ovf {o}")


if __name__ == "__main__":
    main()
