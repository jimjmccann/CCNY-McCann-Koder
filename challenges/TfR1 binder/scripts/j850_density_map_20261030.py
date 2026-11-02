#!/usr/bin/env python3
"""j850 — average ATOMIC DENSITY around the fixed positions, shape before secondary structure.

The question this answers: once an exact set of fixed positions is chosen, what is the average
atomic density around them -- the average SHAPE -- without looking at secondary structure at
all?

Method. Every match from the constellation search is superposed into OUR frame using the pin
correspondence (Kabsch, proper rotations only). The matched chain's CA atoms are then accumulated
into a voxel grid centred on the pin centroid. Averaging over matches gives the mean occupancy
natural scaffolds place at each point around our site -- a SHAPE, with no secondary-structure
assignment anywhere in it.

The grid is reported three ways:
  - a radial profile from the pin centroid,
  - split by half-space: TARGET side vs FREE side, the axis being pin-centroid -> target-centroid,
  - a PyMOL-loadable pseudo-atom PDB of the occupied voxels, b-factor = mean occupancy.

Local CPU only. No GPU and no accelerator.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import j848_cacb_pilot_20261030 as J  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", required=True)
    ap.add_argument("--pins", required=True)
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--limit", type=int, default=20000)
    ap.add_argument("--tol-d", type=float, default=1.0)
    ap.add_argument("--tol-rmsd", type=float, default=2.5)
    ap.add_argument("--gate-ca", type=float, default=1.5)
    ap.add_argument("--gate-ang", type=float, default=90.0)
    ap.add_argument("--voxel", type=float, default=2.0)
    ap.add_argument("--half-box", type=float, default=32.0)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    from scipy.spatial import cKDTree
    ks, (qca, qcb) = J.read_query(args.pins, args.target_chain)
    tgt = J.load_target(args.pins, args.target_chain)
    tree = cKDTree(tgt)
    cen_q = qca.mean(0)
    tdir = tgt.mean(0) - cen_q
    tdir /= np.linalg.norm(tdir)
    print(f"[j850] pins {len(ks)}  target heavy atoms {len(tgt)}")
    print(f"[j850] target direction from pin centroid: "
          f"[{tdir[0]:.3f} {tdir[1]:.3f} {tdir[2]:.3f}]")

    CA = np.load(args.store + "_ca.npy", mmap_mode="r")
    CB = np.load(args.store + "_cb.npy", mmap_mode="r")
    off = np.load(args.store + "_offsets.npy")

    hb, vx = args.half_box, args.voxel
    n = int(2 * hb / vx)
    grid = np.zeros((n, n, n), dtype=np.float32)
    n_match = 0
    all_r, all_proj = [], []

    for c in range(min(args.limit, len(off) - 1)):
        s, e = off[c], off[c + 1]
        Xca = np.asarray(CA[s:e])
        res = J.search_chain_cacb(Xca, np.asarray(CB[s:e]), qca, qcb,
                                  args.tol_d, args.tol_rmsd, 2000000, tree, 4.0, 8.0)
        for h in res:
            if h[1] > args.gate_ca or h[2] > args.gate_ang:
                continue
            pins_idx = np.asarray(h[5])
            P = Xca[pins_idx].astype(np.float64)
            R = J.rot_batch(P[None, :, :], qca)[0]
            Xt = (Xca.astype(np.float64) - P.mean(0)) @ R.T + cen_q
            d = Xt - cen_q
            all_r.append(np.linalg.norm(d, axis=1))
            all_proj.append(d @ tdir)
            ij = np.floor((d + hb) / vx).astype(int)
            ok = np.all((ij >= 0) & (ij < n), axis=1)
            ij = ij[ok]
            np.add.at(grid, (ij[:, 0], ij[:, 1], ij[:, 2]), 1.0)
            n_match += 1

    print(f"[j850] matches accumulated: {n_match}")
    if n_match == 0:
        print("[j850] nothing to map")
        return
    grid /= n_match
    r = np.concatenate(all_r)
    proj = np.concatenate(all_proj)

    print(f"\n=== RADIAL PROFILE: mean CA atoms per match in each shell ===")
    print(f"{'shell (A)':>12} {'mean CA/match':>14} {'cumulative':>12}")
    edges = [0, 5, 10, 15, 20, 25, 30, 40, 60]
    cum = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        v = float(((r >= a) & (r < b)).sum()) / n_match
        cum += v
        print(f"{a:>5d}-{b:<6d} {v:>14.2f} {cum:>12.2f}")

    print(f"\n=== HALF-SPACE SPLIT along pin-centroid -> target-centroid ===")
    nt = float((proj > 0).sum()) / n_match
    nf = float((proj <= 0).sum()) / n_match
    print(f"  TARGET side (toward TfR1) : {nt:>8.2f} CA per match  ({100 * nt / (nt + nf):.1f}%)")
    print(f"  FREE side   (away)        : {nf:>8.2f} CA per match  ({100 * nf / (nt + nf):.1f}%)")
    print(f"  ratio target/free         : {nt / max(nf, 1e-9):>8.2f}")
    print(f"  projection along axis (A) : min {proj.min():.1f}  median {np.median(proj):.1f}  "
          f"max {proj.max():.1f}")

    occ = np.argwhere(grid > 0)
    vals = grid[grid > 0]
    keep = vals >= np.percentile(vals, 50)
    occ, vals = occ[keep], vals[keep]
    outp = args.out_prefix + "_density.pdb"
    with open(outp, "w") as fh:
        fh.write(f"REMARK  j850 mean CA occupancy, {n_match} matches, voxel {vx} A\n")
        fh.write(f"REMARK  b-factor = mean CA per match per voxel (x1000)\n")
        for i, (ijk, v) in enumerate(zip(occ, vals)):
            p = (ijk + 0.5) * vx - hb + cen_q
            fh.write(f"HETATM{i % 99999 + 1:5d}  O   HOH X{i % 9999 + 1:4d}    "
                     f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00{min(999.99, v * 1000):6.2f}"
                     f"           O\n")
    np.save(args.out_prefix + "_grid.npy", grid)
    print(f"\n[j850] wrote {outp} ({len(occ)} voxels at/above median occupancy)")
    print(f"[j850] wrote {args.out_prefix}_grid.npy")


if __name__ == "__main__":
    main()
