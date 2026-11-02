#!/usr/bin/env python3
"""j855 — average atomic density around the fixed positions, masked to the FREE volume,
with a rotation null.

The original question: the average atomic density relative to the fixed positions -- the
average SHAPE, without looking at secondary structure at all.

The clash filter was retired: it could not separate real matches from randomly oriented ones,
because the pin centroid sits 1.91 A inside the target envelope, so everything overlaps the
target at every orientation. Two consequences for the density:

  1. The target-occupied volume must be MASKED OUT. Mass that a matched scaffold places where
     TfR1 already is tells us nothing about what a binder could build.
  2. The density needs a null, or it is uninterpretable; the retired clash filter is the proof.

THE NULL. For each hit, the actual superposed coordinates are rotated by a random rotation
ABOUT THE PIN CENTROID. That preserves every atom's distance from the centroid exactly -- so
the radial profile of null and actual are identical by construction -- and destroys only the
ANGULAR information. Any directional structure left in the actual density is therefore real,
and not an artefact of chain size, standoff or the radial distribution.

Local CPU only. No GPU and no accelerator.
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import j851_funcatom_search_20261030 as J      # noqa: E402
import j853_funcatom_density_20261030 as K     # noqa: E402


def rand_rot(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", required=True)
    ap.add_argument("--pins", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--subset", required=True)
    ap.add_argument("--afdb-tars", default="")
    ap.add_argument("--cif-dirs", default="")
    ap.add_argument("--tol-d", type=float, default=1.5)
    ap.add_argument("--tol-rmsd", type=float, default=2.0)
    ap.add_argument("--voxel", type=float, default=2.0)
    ap.add_argument("--radius", type=float, default=25.0)
    ap.add_argument("--free-cut", type=float, default=4.0)
    ap.add_argument("--n-null", type=int, default=50)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    from scipy.spatial import cKDTree
    rng = np.random.default_rng(11)

    labels, Q, qtyp = J.read_query(args.pins, args.spec)
    keys = [l.split(":")[0] for l in labels]
    sub = [keys.index(w.strip()) for w in args.subset.split(",")]
    Qs, qts = Q[sub], qtyp[sub]
    qd = np.sqrt(((Qs[:, None, :] - Qs[None, :, :]) ** 2).sum(-1))
    cen = Qs.mean(0)

    tgt = np.asarray([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                      for l in open(args.pins)
                      if l.startswith(("ATOM", "HETATM")) and l[21] == "A"
                      and l[76:78].strip() != "H"], dtype=np.float64)
    tree = cKDTree(tgt)
    tdir = tgt.mean(0) - cen
    tdir /= np.linalg.norm(tdir)
    print(f"[j855] pins {[keys[i] for i in sub]}  target atoms {len(tgt)}  "
          f"centroid->target {tree.query(cen)[0]:.2f} A")

    # ---- free-volume voxel mask
    R, vx = args.radius, args.voxel
    n = int(2 * R / vx)
    ax = (np.arange(n) + 0.5) * vx - R
    GX, GY, GZ = np.meshgrid(ax, ax, ax, indexing="ij")
    cpts = np.stack([GX.ravel(), GY.ravel(), GZ.ravel()], 1)
    inside = np.linalg.norm(cpts, axis=1) <= R
    dfree, _ = tree.query(cpts + cen, k=1)
    freemask = (dfree > args.free_cut) & inside
    FREE = freemask.reshape(n, n, n)
    print(f"[j855] grid {n}^3 at {vx} A; in-sphere {int(inside.sum())}, "
          f"FREE {int(freemask.sum())} ({100 * freemask.sum() / inside.sum():.1f}% of sphere)")

    # ---- collect hits
    cenS = np.load(args.store + "_cen.npy", mmap_mode="r")
    typS = np.load(args.store + "_typ.npy", mmap_mode="r")
    ridS = np.load(args.store + "_rid.npy", mmap_mode="r")
    off = np.load(args.store + "_offsets.npy")
    names = open(args.store + "_names.txt").read().split("\n")
    hits = []
    for c in range(len(off) - 1):
        s, e = off[c], off[c + 1]
        res = J.search_chain(np.asarray(cenS[s:e]), np.asarray(typS[s:e]),
                             qd, qts, args.tol_d, args.tol_rmsd, Qs, 2000000)
        if res:
            r, idxs = min(res)
            hits.append((r, names[c], tuple(int(ridS[s + i]) for i in idxs)))
    hits.sort()
    print(f"[j855] {len(hits)} chains to superpose")

    afdb_idx = K.index_afdb(args.afdb_tars) if args.afdb_tars else None
    cif_dirs = [d for d in args.cif_dirs.split(",") if d]

    gA = np.zeros((n, n, n), dtype=np.float64)
    gN = np.zeros((n, n, n), dtype=np.float64)
    nA = nN = 0
    projA, projN = [], []

    def deposit(grid, X):
        d = X - cen
        ij = np.floor((d + R) / vx).astype(int)
        m = np.all((ij >= 0) & (ij < n), axis=1)
        ij = ij[m]
        np.add.at(grid, (ij[:, 0], ij[:, 1], ij[:, 2]), 1.0)

    for r, nm, resnums in hits:
        src, chain_id = nm.rsplit("_", 1)
        got = K.load_chain(src, chain_id, afdb_idx, cif_dirs)
        if got is None:
            continue
        atoms, ca, _ = got
        P, ok = [], True
        for q in resnums:
            if q not in atoms:
                ok = False
                break
            fc = K.func_centroid(*atoms[q])
            if fc is None:
                ok = False
                break
            P.append(fc)
        if not ok:
            continue
        P = np.asarray(P)
        Rot = K.rot_for(P, Qs)
        Xt = (ca - P.mean(0)) @ Rot.T + cen
        deposit(gA, Xt)
        nA += 1
        dd = Xt - cen
        keep = np.linalg.norm(dd, axis=1) <= R
        # only count atoms that land in FREE voxels for the half-space statistic
        ij = np.floor((dd[keep] + R) / vx).astype(int)
        m = np.all((ij >= 0) & (ij < n), axis=1)
        ij = ij[m]
        fm = FREE[ij[:, 0], ij[:, 1], ij[:, 2]]
        projA.append((dd[keep][m][fm]) @ tdir)
        # ---- NULL: rotate the SUPERPOSED coords about the pin centroid
        for _ in range(args.n_null):
            Yt = (Xt - cen) @ rand_rot(rng).T + cen
            deposit(gN, Yt)
            nN += 1
            ddy = Yt - cen
            ky = np.linalg.norm(ddy, axis=1) <= R
            ijy = np.floor((ddy[ky] + R) / vx).astype(int)
            my = np.all((ijy >= 0) & (ijy < n), axis=1)
            ijy = ijy[my]
            fmy = FREE[ijy[:, 0], ijy[:, 1], ijy[:, 2]]
            projN.append((ddy[ky][my][fmy]) @ tdir)

    print(f"[j855] actual placements {nA}, null placements {nN}")
    if nA == 0:
        return
    gA /= nA
    gN /= max(1, nN)

    # ---- radial profile check: must be identical by construction
    rad = np.linalg.norm(cpts, axis=1).reshape(n, n, n)
    print("\n=== SANITY: radial profile, actual vs null (must agree -- the null preserves it) ===")
    print(f"{'shell':>9} {'actual':>9} {'null':>9}")
    for a, b in zip([0, 5, 10, 15, 20], [5, 10, 15, 20, 25]):
        sh = (rad >= a) & (rad < b)
        print(f"{a:>3}-{b:<5} {gA[sh].sum():>9.2f} {gN[sh].sum():>9.2f}")

    print("\n=== THE MEASUREMENT: mean CA per chain, FREE volume only ===")
    print(f"{'shell':>9} {'actual':>9} {'null':>9} {'ratio':>7}")
    for a, b in zip([0, 5, 10, 15, 20], [5, 10, 15, 20, 25]):
        sh = (rad >= a) & (rad < b) & FREE
        va, vn = gA[sh].sum(), gN[sh].sum()
        print(f"{a:>3}-{b:<5} {va:>9.2f} {vn:>9.2f} {va / max(vn, 1e-9):>7.2f}")
    fa, fn = gA[FREE].sum(), gN[FREE].sum()
    oa, on = gA[~FREE & (rad <= R)].sum(), gN[~FREE & (rad <= R)].sum()
    print(f"\n  FREE  total: actual {fa:.2f}  null {fn:.2f}  ratio {fa / max(fn, 1e-9):.3f}")
    print(f"  TARGET-occupied: actual {oa:.2f}  null {on:.2f}  ratio {oa / max(on, 1e-9):.3f}")

    pa = np.concatenate(projA) if projA else np.asarray([])
    pn = np.concatenate(projN) if projN else np.asarray([])
    if len(pa) and len(pn):
        print("\n=== HALF-SPACE, FREE volume only, along pin-centroid -> target ===")
        for lab, p in (("actual", pa), ("null  ", pn)):
            t = float((p > 0).mean())
            print(f"  {lab}: target side {100 * t:.1f}%   free side {100 * (1 - t):.1f}%")

    enr = np.zeros_like(gA)
    ok = (gN > 0) & FREE
    enr[ok] = gA[ok] / gN[ok]
    np.save(args.out_prefix + "_actual.npy", gA)
    np.save(args.out_prefix + "_null.npy", gN)
    np.save(args.out_prefix + "_enrich.npy", enr)
    sel = np.argwhere(ok & (gA > np.percentile(gA[ok], 90)))
    with open(args.out_prefix + "_enriched.pdb", "w") as fh:
        fh.write(f"REMARK  j855 enrichment actual/null, FREE voxels, {nA} chains\n")
        for i, ijk in enumerate(sel):
            p = (ijk + 0.5) * vx - R + cen
            v = enr[tuple(ijk)]
            fh.write(f"HETATM{i % 99999 + 1:5d}  O   HOH X{i % 9999 + 1:4d}    "
                     f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00{min(999.99, v * 100):6.2f}"
                     f"           O\n")
    print(f"\n[j855] wrote {args.out_prefix}_enriched.pdb ({len(sel)} voxels), "
          f"_actual/_null/_enrich .npy")
    if len(enr[ok]):
        e = enr[ok]
        print(f"[j855] enrichment over FREE voxels: median {np.median(e):.2f}  "
              f"p90 {np.percentile(e, 90):.2f}  max {e.max():.2f}")


if __name__ == "__main__":
    main()
