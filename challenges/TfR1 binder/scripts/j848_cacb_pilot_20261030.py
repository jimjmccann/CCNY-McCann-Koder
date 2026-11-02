#!/usr/bin/env python3
"""j848 — PILOT: does adding CB restore selectivity to the constellation search?

The CA-only 5-point constellation search was measured to have essentially NO selectivity:
the pin set matched 73% of chains, and three control constellations drawn from real PDB
chains matched at the same rate or higher. The proposal tested here: add the beta carbon,
since CA->CB establishes orientation without needing full side chains.

This pilot answers ONE question, on real structures ALREADY ON DISK, before the full archive
mirror finishes: does the CA+CB representation separate the constellation from matched
controls, where CA alone did not?

Representation: each pin contributes TWO points, CA and CB, so a 5-pin query is a 10-point
constellation. Candidates are enumerated on the CA sub-constellation (reusing j846's batched
matcher) and then gated on the full 10-point Kabsch RMSD.

Gly has no CB and can never match a pin. That is correct, not a bug.

Local CPU only. No GPU and no accelerator.
"""
import argparse
import glob
import itertools
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from j846_pdb_constellation_search_20261030 import kabsch_rmsd_batch  # noqa: E402


# ----------------------------------------------------------------- extract

def extract_one(path):
    """-> list of (chain_label, ca (n,3) float32, cb (n,3) float32). Gly dropped."""
    import gemmi
    try:
        st = gemmi.read_structure(path)
    except Exception:
        return []
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    out = []
    if len(st) == 0:
        return []
    model = st[0]
    base = os.path.basename(path).split(".")[0]
    for chain in model:
        ca, cb = [], []
        for res in chain:
            a_ca = res.find_atom("CA", "*")
            a_cb = res.find_atom("CB", "*")
            if a_ca is None or a_cb is None:
                continue
            ca.append([a_ca.pos.x, a_ca.pos.y, a_ca.pos.z])
            cb.append([a_cb.pos.x, a_cb.pos.y, a_cb.pos.z])
        if len(ca) >= 5:
            out.append((f"{base}_{chain.name}",
                        np.asarray(ca, dtype=np.float32),
                        np.asarray(cb, dtype=np.float32)))
    return out


def build(dirs, out_prefix, nproc):
    files = []
    for d in dirs:
        files.extend(sorted(glob.glob(os.path.join(d, "*.cif.gz"))))
    print(f"[build] {len(files)} structure files", flush=True)
    t0 = time.time()
    CA, CB, off, names = [], [], [0], []
    cur = 0
    with Pool(nproc) as pool:
        for i, res in enumerate(pool.imap_unordered(extract_one, files, chunksize=32)):
            for nm, ca, cb in res:
                CA.append(ca)
                CB.append(cb)
                cur += len(ca)
                off.append(cur)
                names.append(nm)
            if (i + 1) % 5000 == 0:
                print(f"  .. {i + 1}/{len(files)}  chains {len(names)}  "
                      f"{time.time() - t0:.0f}s", flush=True)
    np.save(out_prefix + "_ca.npy", np.concatenate(CA))
    np.save(out_prefix + "_cb.npy", np.concatenate(CB))
    np.save(out_prefix + "_offsets.npy", np.asarray(off, dtype=np.int64))
    open(out_prefix + "_names.txt", "w").write("\n".join(names))
    print(f"[build] chains {len(names)}  residues {cur}  {time.time() - t0:.0f}s", flush=True)


# ----------------------------------------------------------------- query

def read_query(path, target_chain="A"):
    """CA and CB of every non-target chain. Chain is NAMED, never inferred."""
    ca, cb = {}, {}
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            ch, an = ln[21], ln[12:16].strip()
            if ch == target_chain:
                continue
            p = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            if an == "CA":
                ca[ch] = p
            elif an == "CB":
                cb[ch] = p
    ks = sorted(k for k in ca if k in cb)
    return ks, (np.asarray([ca[k] for k in ks], dtype=np.float64),
                np.asarray([cb[k] for k in ks], dtype=np.float64))


# ----------------------------------------------------------------- search

def search_chain_cacb(Xca, Xcb, qca, qcb, tol_d, tol_rmsd, cap,
                      tgt_tree=None, clash_cut=4.0, away_cut=8.0):
    """CA enumeration then full CA+CB 10-point gate. Returns list of (rmsd10, rmsd_ca, idx)."""
    L = len(Xca)
    k = len(qca)
    if L < k:
        return []
    qd = np.sqrt(((qca[:, None, :] - qca[None, :, :]) ** 2).sum(-1))
    D = np.sqrt(((Xca[:, None, :] - Xca[None, :, :]) ** 2).sum(-1))

    idx = list(range(k))
    shells, counts = {}, {}
    for a, b in itertools.combinations(idx, 2):
        M = np.abs(D - qd[a, b]) <= tol_d
        np.fill_diagonal(M, False)
        if not M.any():
            return []
        shells[(a, b)] = shells[(b, a)] = M
        counts[(a, b)] = counts[(b, a)] = int(M.sum())

    first = min(itertools.combinations(idx, 2), key=lambda p: counts[p])
    order = [first[0], first[1]]
    rest = [i for i in idx if i not in order]
    while rest:
        nxt = min(rest, key=lambda c: sum(counts[(c, p)] for p in order))
        order.append(nxt)
        rest.remove(nxt)

    A = np.argwhere(shells[(order[0], order[1])])
    if len(A) == 0 or len(A) > cap:
        return []
    for t in range(2, k):
        iq = order[t]
        mask = np.ones((L, len(A)), dtype=bool)
        for pos, p in enumerate(order[:t]):
            mask &= shells[(iq, p)][:, A[:, pos]]
        for pos in range(t):
            mask[A[:, pos], np.arange(len(A))] = False
        rows, cols = np.nonzero(mask)
        if len(rows) == 0 or len(rows) > cap:
            return []
        A = np.concatenate([A[cols], rows[:, None]], axis=1)

    # reorder columns back to pin order, then gate
    inv = np.argsort(np.asarray(order))
    A = A[:, inv]
    out = []
    CH = 20000
    Q10 = np.concatenate([qca, qcb], axis=0)
    # query CA->CB unit vectors, in the query frame
    qv = qcb - qca
    qv /= np.linalg.norm(qv, axis=1, keepdims=True)
    for s in range(0, len(A), CH):
        blk = A[s:s + CH]
        Pca = Xca[blk].astype(np.float64)
        Pcb = Xcb[blk].astype(np.float64)
        r_ca = kabsch_rmsd_batch(Pca, qca)
        P10 = np.concatenate([Pca, Pcb], axis=1)
        r10 = kabsch_rmsd_batch(P10, Q10)

        # ORIENTATION gate: superpose on CA ONLY, then compare CA->CB directions.
        # A 10-point positional RMSD lets a wrong CB direction hide behind good CA
        # positions; this does not. Returns the WORST of the k angles, in degrees.
        R = rot_batch(Pca, qca)
        v = Pcb - Pca
        v /= np.linalg.norm(v, axis=2, keepdims=True)
        vr = np.einsum("nij,nkj->nki", R, v)          # rotate into the query frame
        cos = np.clip((vr * qv[None, :, :]).sum(-1), -1.0, 1.0)
        ang = np.degrees(np.arccos(cos)).max(axis=1)

        keep = np.flatnonzero(r10 <= tol_rmsd)

        # ---- TARGET-CLASH FILTER (later retired; see ../methods/) -----------
        # Superpose the candidate's WHOLE chain into our frame using the pin
        # match, then reject it if its backbone runs into TfR1. What survives
        # is the subset that presents the constellation AND goes the other way,
        # which is the only subset that could ever be a binder here.
        if tgt_tree is not None and len(keep):
            n_clash = np.zeros(len(keep), dtype=np.int32)
            n_away = np.zeros(len(keep), dtype=np.int32)
            for jj, j in enumerate(keep):
                Rj = R[j]
                cen_p = Pca[j].mean(0)
                cen_q = qca.mean(0)
                Xt = (Xca.astype(np.float64) - cen_p) @ Rj.T + cen_q
                d, _ = tgt_tree.query(Xt, k=1)
                n_clash[jj] = int((d <= clash_cut).sum())
                n_away[jj] = int((d > away_cut).sum())
            for jj, j in enumerate(keep):
                out.append((float(r10[j]), float(r_ca[j]), float(ang[j]),
                            int(n_clash[jj]), int(n_away[jj]),
                            tuple(int(v2) for v2 in blk[j])))
        else:
            for j in keep:
                out.append((float(r10[j]), float(r_ca[j]), float(ang[j]), -1, -1,
                            tuple(int(v2) for v2 in blk[j])))
    return out


def rot_batch(P, Q):
    """Optimal proper rotation taking centred P onto centred Q, batched."""
    Pc = P - P.mean(axis=1, keepdims=True)
    Qc = Q - Q.mean(axis=0, keepdims=True)
    H = np.einsum("nki,kj->nij", Pc, Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("nij,njk->nik", U, Vt)))
    D = np.zeros((len(P), 3, 3))
    D[:, 0, 0] = D[:, 1, 1] = 1.0
    D[:, 2, 2] = d
    return np.einsum("nij,njk,nkl->nil", U, D, Vt).transpose(0, 2, 1)


def load_target(path, target_chain="A"):
    """Heavy atoms of the TARGET chain, in our input's own frame. Chain is NAMED."""
    P = []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")) and ln[21] == target_chain:
            if ln[76:78].strip() != "H":
                P.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.asarray(P, dtype=np.float64)


def run_query(args_tuple):
    (label, qca, qcb, prefix, tol_d, tol_rmsd, cap, lo, hi,
     tgt, clash_cut, away_cut) = args_tuple
    tree = None
    if tgt is not None and len(tgt):
        from scipy.spatial import cKDTree
        tree = cKDTree(tgt)
    CA = np.load(prefix + "_ca.npy", mmap_mode="r")
    CB = np.load(prefix + "_cb.npy", mmap_mode="r")
    off = np.load(prefix + "_offsets.npy")
    n_hit = n_chain = 0
    rows = []
    for c in range(lo, hi):
        s, e = off[c], off[c + 1]
        res = search_chain_cacb(np.asarray(CA[s:e]), np.asarray(CB[s:e]),
                                qca, qcb, tol_d, tol_rmsd, cap,
                                tree, clash_cut, away_cut)
        if res:
            n_chain += 1
            n_hit += len(res)
            rows.extend((r[0], r[1], r[2], r[3], r[4]) for r in res)
    return label, n_hit, n_chain, np.asarray(rows).reshape(-1, 5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", nargs="+", metavar="DIR")
    ap.add_argument("--out-prefix", default="")
    ap.add_argument("--search", metavar="PREFIX")
    ap.add_argument("--pins", default="")
    ap.add_argument("--controls", default="", help="comma list of control pin PDBs")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--tol-d", type=float, default=1.0)
    ap.add_argument("--tol-rmsd", type=float, default=2.5)
    ap.add_argument("--cap", type=int, default=2000000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--nproc", type=int, default=16)
    ap.add_argument("--target-pdb", default="", help="PDB whose target chain to clash against")
    ap.add_argument("--clash-cut", type=float, default=4.0)
    ap.add_argument("--away-cut", type=float, default=8.0)
    args = ap.parse_args()

    if args.build:
        build(args.build, args.out_prefix, args.nproc)
        return

    off = np.load(args.search + "_offsets.npy")
    n_chains = len(off) - 1
    hi = n_chains if not args.limit else min(n_chains, args.limit)
    print(f"[pilot] store chains {n_chains}, scanning {hi}", flush=True)

    queries = []
    ks, (qca, qcb) = read_query(args.pins, args.target_chain)
    print(f"[query] ours pins {ks}", flush=True)
    queries.append(("ours", qca, qcb))
    for i, p in enumerate([x for x in args.controls.split(",") if x]):
        kc, (cca, ccb) = read_query(p, args.target_chain)
        queries.append((f"ctrl{i + 1}", cca, ccb))

    TGT = None
    if args.target_pdb:
        TGT = load_target(args.target_pdb, args.target_chain)
        print(f"[target] {len(TGT)} heavy atoms from chain {args.target_chain} "
              f"of {os.path.basename(args.target_pdb)}", flush=True)

    # shard each query across processes
    jobs = []
    nsh = max(1, args.nproc // max(1, len(queries)))
    for label, a, b in queries:
        bounds = np.linspace(0, hi, nsh + 1).astype(int)
        for j in range(nsh):
            jobs.append((label, a, b, args.search, args.tol_d, args.tol_rmsd,
                         args.cap, int(bounds[j]), int(bounds[j + 1]),
                         TGT, args.clash_cut, args.away_cut))
    t0 = time.time()
    agg = {}
    with Pool(args.nproc) as pool:
        for label, nh, nc, rr in pool.imap_unordered(run_query, jobs):
            h, c, arr = agg.get(label, (0, 0, []))
            agg[label] = (h + nh, c + nc, arr + [rr])
    print(f"[pilot] {time.time() - t0:.0f}s\n", flush=True)

    allr = {l: (np.concatenate(v[2]) if len(v[2]) else np.zeros((0, 3)))
            for l, v in agg.items()}

    print("=== A. 10-point POSITIONAL rmsd (CA+CB folded together) ===")
    print(f"{'tol':>9} " + " ".join(f"{l:>14}" for l, _, _ in queries))
    for t in [2.5, 2.0, 1.5, 1.0, 0.75, 0.5]:
        cells = [f"{int((allr[l][:, 0] <= t).sum()):>14d}" for l, _, _ in queries]
        print(f"{t:>9.2f} " + " ".join(cells))

    print("\n=== B. ORIENTATION gate: CA rmsd <= 1.0 A AND worst CA->CB angle <= X deg ===")
    print(f"{'max_ang':>9} " + " ".join(f"{l:>14}" for l, _, _ in queries))
    for a in [180, 90, 60, 45, 30, 20, 15, 10]:
        cells = []
        for l, _, _ in queries:
            d = allr[l]
            cells.append(f"{int(((d[:, 1] <= 1.0) & (d[:, 2] <= a)).sum()):>14d}")
        print(f"{a:>9d} " + " ".join(cells))

    if args.target_pdb:
        print("\n=== C. TARGET-CLASH FILTER: CA rmsd <= 1.5 and worst angle <= 90 deg ===")
        print(f"{'':>26} " + " ".join(f"{l:>14}" for l, _, _ in queries))
        base, clean, usable = {}, {}, {}
        for l, _, _ in queries:
            d = allr[l]
            m = (d[:, 1] <= 1.5) & (d[:, 2] <= 90)
            base[l] = int(m.sum())
            clean[l] = int((m & (d[:, 3] == 0)).sum())
            usable[l] = int((m & (d[:, 3] == 0) & (d[:, 4] >= 20)).sum())
        print(f"{'matches (pre-clash)':>26} " + " ".join(f"{base[l]:>14d}" for l, _, _ in queries))
        print(f"{'NO clash with target':>26} " + " ".join(f"{clean[l]:>14d}" for l, _, _ in queries))
        print(f"{'+ >=20 res in free space':>26} " + " ".join(f"{usable[l]:>14d}" for l, _, _ in queries))
        for l, _, _ in queries:
            if base[l]:
                print(f"    {l}: clash-free {100.0 * clean[l] / base[l]:.2f}% of matches")

    print()
    for l, _, _ in queries:
        d = allr[l]
        h, c, _ = agg.get(l, (0, 0, []))
        best = f"{d[:, 0].min():.3f}" if len(d) else "n/a"
        bang = f"{d[:, 2].min():.1f}" if len(d) else "n/a"
        print(f"{l:>6}: hits {h:>8}  chains_with_hit {c:>6}  "
              f"best_rmsd10 {best}  best_worst_angle {bang} deg")


if __name__ == "__main__":
    main()
