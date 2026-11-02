#!/usr/bin/env python3
"""j871b - does the earlier positive control's NEGATIVE half survive an admissible tol_d?

That control took one real chain's OWN 6 functional groups as the query and reported 0 hits
in 20,000 OTHER chains at tol_d 1.5 / tol_rmsd 2.0, concluding that a 6-group typed
constellation spanning >30 A is "essentially unique in nature". It was measured at the SAME
inadmissible tol_d as the zero it was written to explain, so it is in scope here.

Stronger than the original in one respect: the original used n=1 control chain. This samples
several, so the answer is not a single draw.

Imports j851's matcher via j871's vendored instrumented copy. Nothing is edited.
"""
import argparse
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


J = load("j871", "j871_pincurve_admissible_20261031.py")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--tol-rmsd", type=float, default=2.0)
    ap.add_argument("--tol-d", nargs="+", default=["1.5", "2.5", "floor"])
    ap.add_argument("--scan", type=int, default=20000, help="other chains to search")
    ap.add_argument("--controls", type=int, default=8, help="control chains to sample")
    ap.add_argument("--min-span", type=float, default=30.0)
    ap.add_argument("--cap", type=int, default=2000000)
    args = ap.parse_args()

    if os.path.exists(args.out):
        raise SystemExit(f"[abort] --out exists, refusing to append: {args.out}")

    cen = np.load(args.store + "_cen.npy", mmap_mode="r")
    typ = np.load(args.store + "_typ.npy", mmap_mode="r")
    off = np.load(args.store + "_offsets.npy")
    names = open(args.store + "_names.txt").read().split("\n")
    n = len(off) - 1

    # ---- pick control chains: the FIRST ones that have >= k groups spanning >= min_span,
    # deterministically, so the choice is not a tuned one.
    picks = []
    for c in range(n):
        s, e = off[c], off[c + 1]
        if e - s < args.k:
            continue
        cc = np.asarray(cen[s:e])
        tt = np.asarray(typ[s:e])
        D = np.sqrt(((cc[:, None, :] - cc[None, :, :]) ** 2).sum(-1))
        i, j = np.unravel_index(np.argmax(D), D.shape)
        if D[i, j] < args.min_span:
            continue
        # take the k groups greedily farthest-apart, to mimic a >30 A spread query
        idx = [int(i), int(j)]
        while len(idx) < args.k:
            rest = [q for q in range(len(cc)) if q not in idx]
            if not rest:
                break
            nxt = max(rest, key=lambda q: min(D[q, p] for p in idx))
            idx.append(int(nxt))
        if len(idx) < args.k:
            continue
        Q = cc[idx]
        qt = tt[idx]
        qd = np.sqrt(((Q[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
        picks.append((c, names[c], Q.astype(np.float64), qt, qd, float(qd.max())))
        if len(picks) >= args.controls:
            break

    print(f"[control] {len(picks)} control chains, k={args.k}, "
          f"scanning {min(args.scan, n)} others each")
    tol_ds = [J.floor_tol_d(args.k, args.tol_rmsd) if t == "floor" else float(t)
              for t in args.tol_d]
    print("[control] tol_d values: " + ", ".join(f"{t:.3f}" for t in tol_ds))

    with open(args.out, "w") as fh:
        fh.write("control_chain\tspan\tn_cation\tn_imid\ttol_d\tis_floor\t"
                 "self_rmsd\tother_hits\tother_chains\tbest_other_rmsd\tovf\n")
        for c, nm, Q, qt, qd, span in picks:
            # matches j851_funcatom_search_20261030.py -- CATION, IMID = 0, 1
            nca, nim = int((qt == 0).sum()), int((qt == 1).sum())
            # self-match check at the tightest tol_d
            s, e = off[c], off[c + 1]
            selfr = J.search_chain_instr(np.asarray(cen[s:e]), np.asarray(typ[s:e]),
                                         qd, qt, max(tol_ds), args.tol_rmsd, Q, args.cap)
            for td, raw in zip(tol_ds, args.tol_d):
                hits = chains = ovf = 0
                best = None
                for o in range(min(args.scan, n)):
                    if o == c:
                        continue
                    s2, e2 = off[o], off[o + 1]
                    nh, b, ov = J.search_chain_instr(
                        np.asarray(cen[s2:e2]), np.asarray(typ[s2:e2]),
                        qd, qt, td, args.tol_rmsd, Q, args.cap)
                    if ov:
                        ovf += 1
                    if nh:
                        hits += nh
                        chains += 1
                        best = b if best is None else min(best, b)
                fh.write(f"{nm}\t{span:.2f}\t{nca}\t{nim}\t"
                         f"{td:.4f}\t{int(raw=='floor')}\t"
                         f"{'' if selfr[1] is None else f'{selfr[1]:.4f}'}\t"
                         f"{hits}\t{chains}\t"
                         f"{'' if best is None else f'{best:.4f}'}\t{ovf}\n")
                print(f"  {nm:<22} span {span:5.1f}  tol_d {td:6.3f}"
                      f"{' FLOOR' if raw == 'floor' else '      '}  "
                      f"self {'--' if selfr[1] is None else f'{selfr[1]:.4f}'}  "
                      f"others: hits {hits:>7} chains {chains:>5}  "
                      f"best {'--' if best is None else f'{best:.3f}'}  ovf {ovf}",
                      flush=True)
    print(f"[done] wrote {args.out}")


if __name__ == "__main__":
    main()
