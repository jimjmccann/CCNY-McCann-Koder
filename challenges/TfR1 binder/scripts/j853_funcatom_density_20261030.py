#!/usr/bin/env python3
"""j853 — clash filter and average atomic density for the FUNCTIONAL-ATOM hits.

j850 did this on the CA representation and its conclusions are withdrawn: the CA/CB
coordinates of the pins are a single synthetic rotamer from our own builder, so the
superposition frame carried no site information. j851 fixed the query by matching functional
groups instead. This re-does the two downstream measurements on that valid frame.

For each hit the matched residues are recovered by residue number (stored at build time),
their functional atoms are re-read from the SOURCE structure, the whole chain is superposed
onto our pins through that correspondence, and then:

  - the chain is tested for clash against TfR1 (the filter as originally proposed: keep only
    scaffolds that present the constellation AND go the other way),
  - its CA are accumulated into a voxel grid in our pin frame -- shape, with no
    secondary-structure assignment anywhere.

Local CPU only. No GPU and no accelerator.
"""
import argparse
import glob
import gzip
import itertools
import os
import sys
import tarfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import j851_funcatom_search_20261030 as J  # noqa: E402

ARG_AT, HIS_AT = J.ARG_AT, J.HIS_AT


def rot_for(P, Q):
    """Proper rotation taking centred P onto centred Q."""
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, 1.0, d])
    return (U @ D @ Vt).T


def index_afdb(tar_dir):
    """member name -> tar path, scanned once."""
    idx = {}
    for t in sorted(glob.glob(os.path.join(tar_dir, "*.tar"))):
        with tarfile.open(t) as tf:
            for nm in tf.getnames():
                if nm.endswith(".cif.gz"):
                    idx[os.path.basename(nm).split(".")[0]] = (t, nm)
    return idx


def load_chain(src, chain_id, afdb_idx, cif_dirs):
    """-> (resnum -> {atom: xyz}, CA array (n,3), resnum array). None on failure."""
    import gemmi
    st = None
    if afdb_idx is not None and src in afdb_idx:
        t, nm = afdb_idx[src]
        with tarfile.open(t) as tf:
            raw = gzip.decompress(tf.extractfile(nm).read())
        doc = gemmi.cif.read_string(raw.decode())
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        for d in cif_dirs:
            p = os.path.join(d, src + ".cif.gz")
            if os.path.exists(p):
                st = gemmi.read_structure(p)
                break
    if st is None or len(st) == 0:
        return None
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    for ch in st[0]:
        if ch.name != chain_id:
            continue
        atoms, ca, rn = {}, [], []
        for res in ch:
            a = {x.name: (x.pos.x, x.pos.y, x.pos.z) for x in res}
            atoms[res.seqid.num] = (res.name, a)
            if "CA" in a:
                ca.append(a["CA"])
                rn.append(res.seqid.num)
        return atoms, np.asarray(ca, dtype=np.float64), np.asarray(rn)
    return None


def func_centroid(resname, atoms):
    if resname == "ARG" and all(a in atoms for a in ARG_AT):
        return np.mean([atoms[a] for a in ARG_AT], axis=0)
    if resname == "LYS" and "NZ" in atoms:
        return np.asarray(atoms["NZ"], dtype=float)
    if resname == "HIS" and all(a in atoms for a in HIS_AT):
        return np.mean([atoms[a] for a in HIS_AT], axis=0)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", required=True)
    ap.add_argument("--pins", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--subset", required=True, help="comma list of pin keys, e.g. B1,D1,E1,F1,G1")
    ap.add_argument("--tol-d", type=float, default=1.5)
    ap.add_argument("--tol-rmsd", type=float, default=2.0)
    ap.add_argument("--afdb-tars", default="")
    ap.add_argument("--cif-dirs", default="")
    ap.add_argument("--clash-cut", type=float, default=4.0)
    ap.add_argument("--voxel", type=float, default=2.0)
    ap.add_argument("--half-box", type=float, default=32.0)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--max-hits", type=int, default=0)
    args = ap.parse_args()

    from scipy.spatial import cKDTree
    labels, Q, qtyp = J.read_query(args.pins, args.spec)
    keys = [l.split(":")[0] for l in labels]
    want = [k.strip() for k in args.subset.split(",")]
    sub = [keys.index(w) for w in want]
    Qs, qts = Q[sub], qtyp[sub]
    qd = np.sqrt(((Qs[:, None, :] - Qs[None, :, :]) ** 2).sum(-1))
    cen_q = Qs.mean(0)

    tgt = J.read_query.__globals__  # unused; target read below
    tgt = []
    for ln in open(args.pins):
        if ln.startswith(("ATOM", "HETATM")) and ln[21] == "A" and ln[76:78].strip() != "H":
            tgt.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    tgt = np.asarray(tgt, dtype=np.float64)
    tree = cKDTree(tgt)
    tdir = tgt.mean(0) - cen_q
    tdir /= np.linalg.norm(tdir)
    print(f"[j853] subset {want}  target heavy atoms {len(tgt)}")

    cen = np.load(args.store + "_cen.npy", mmap_mode="r")
    typ = np.load(args.store + "_typ.npy", mmap_mode="r")
    rid = np.load(args.store + "_rid.npy", mmap_mode="r")
    off = np.load(args.store + "_offsets.npy")
    names = open(args.store + "_names.txt").read().split("\n")

    hits = []
    for c in range(len(off) - 1):
        s, e = off[c], off[c + 1]
        res = J.search_chain(np.asarray(cen[s:e]), np.asarray(typ[s:e]),
                             qd, qts, args.tol_d, args.tol_rmsd, Qs, 2000000)
        for r, idxs in res:
            hits.append((r, c, tuple(int(rid[s + i]) for i in idxs)))
    hits.sort()
    if args.max_hits:
        hits = hits[:args.max_hits]
    print(f"[j853] {len(hits)} hits to process")

    afdb_idx = index_afdb(args.afdb_tars) if args.afdb_tars else None
    if afdb_idx:
        print(f"[j853] indexed {len(afdb_idx)} AFDB members")
    cif_dirs = [d for d in args.cif_dirs.split(",") if d]

    hb, vx = args.half_box, args.voxel
    n = int(2 * hb / vx)
    grid = np.zeros((n, n, n), dtype=np.float32)
    rows = []
    allr, allproj = [], []
    nproc = 0
    for r, c, resnums in hits:
        nm = names[c]
        src, chain_id = nm.rsplit("_", 1)
        got = load_chain(src, chain_id, afdb_idx, cif_dirs)
        if got is None:
            continue
        atoms, ca, rn = got
        P = []
        ok = True
        for q in resnums:
            if q not in atoms:
                ok = False
                break
            fc = func_centroid(*atoms[q])
            if fc is None:
                ok = False
                break
            P.append(fc)
        if not ok or len(P) != len(sub):
            continue
        P = np.asarray(P)
        R = rot_for(P, Qs)
        Xt = (ca - P.mean(0)) @ R.T + cen_q
        d, _ = tree.query(Xt, k=1)
        nclash = int((d <= args.clash_cut).sum())
        rows.append((nm, r, nclash, len(ca), float(d.min())))
        dd = Xt - cen_q
        allr.append(np.linalg.norm(dd, axis=1))
        allproj.append(dd @ tdir)
        ij = np.floor((dd + hb) / vx).astype(int)
        m = np.all((ij >= 0) & (ij < n), axis=1)
        ij = ij[m]
        np.add.at(grid, (ij[:, 0], ij[:, 1], ij[:, 2]), 1.0)
        nproc += 1

    print(f"[j853] superposed {nproc} chains")
    if nproc == 0:
        return
    grid /= nproc
    nc = np.asarray([x[2] for x in rows])
    print(f"\n=== CLASH AGAINST TfR1 (cut {args.clash_cut} A, CA vs target heavy atoms) ===")
    print(f"  clash-free chains      : {int((nc == 0).sum())} of {nproc}")
    print(f"  <=5 clashing CA        : {int((nc <= 5).sum())}")
    print(f"  <=10 clashing CA       : {int((nc <= 10).sum())}")
    print(f"  clashing CA  min {nc.min()}  median {int(np.median(nc))}  max {nc.max()}")

    r = np.concatenate(allr)
    proj = np.concatenate(allproj)
    print(f"\n=== RADIAL PROFILE (mean CA per chain per shell) ===")
    cum = 0.0
    for a, b in zip([0, 5, 10, 15, 20, 25, 30, 40], [5, 10, 15, 20, 25, 30, 40, 60]):
        v = float(((r >= a) & (r < b)).sum()) / nproc
        cum += v
        print(f"  {a:>3}-{b:<3} {v:>8.2f}   cum {cum:>8.2f}")
    nt = float((proj > 0).sum()) / nproc
    nf = float((proj <= 0).sum()) / nproc
    print(f"\n=== HALF-SPACE along pin-centroid -> target-centroid ===")
    print(f"  TARGET side {nt:>8.2f} CA/chain ({100 * nt / (nt + nf):.1f}%)")
    print(f"  FREE   side {nf:>8.2f} CA/chain ({100 * nf / (nt + nf):.1f}%)")
    print(f"  ratio       {nt / max(nf, 1e-9):>8.2f}")

    with open(args.out_prefix + "_hits.tsv", "w") as fh:
        fh.write("chain\trmsd\tn_clash_CA\tchain_len\tmin_dist_target\n")
        for nm, rr, ncl, L, md in sorted(rows, key=lambda x: (x[2], x[1])):
            fh.write(f"{nm}\t{rr:.3f}\t{ncl}\t{L}\t{md:.2f}\n")
    occ = np.argwhere(grid > 0)
    vals = grid[grid > 0]
    keep = vals >= np.percentile(vals, 50)
    occ, vals = occ[keep], vals[keep]
    with open(args.out_prefix + "_density.pdb", "w") as fh:
        fh.write(f"REMARK  j853 mean CA occupancy, {nproc} chains, voxel {vx} A\n")
        for i, (ijk, v) in enumerate(zip(occ, vals)):
            p = (ijk + 0.5) * vx - hb + cen_q
            fh.write(f"HETATM{i % 99999 + 1:5d}  O   HOH X{i % 9999 + 1:4d}    "
                     f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00"
                     f"{min(999.99, v * 1000):6.2f}           O\n")
    np.save(args.out_prefix + "_grid.npy", grid)
    print(f"\n[j853] wrote {args.out_prefix}_hits.tsv and _density.pdb ({len(occ)} voxels)")


if __name__ == "__main__":
    main()
