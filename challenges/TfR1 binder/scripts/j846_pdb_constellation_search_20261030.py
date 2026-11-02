#!/usr/bin/env python3
"""j846 — mine the PDB for backbone scaffolds matching our fixed-residue constellation.

This is pass A of the mining line: match the guidepost CA positions only, with residue
identity IGNORED, against every chain in the local foldseek PDB100 set. The match criterion
is positional and not compositional. Pass A answers the SCAFFOLD question only. It cannot
speak to chemistry -- a CA match leaves the functional group free to rotate.

Two modes:

  --build    parse the extracted CA multimodel PDB into a compact numpy store (one-time)
  --search   run the constellation search and emit one row per hit

MATCHER. Per chain: build the 10 distance-shell boolean matrices, then grow partial assignments
in BATCH with numpy, never a Python loop over candidate pairs. Start from the shell with the
FEWEST candidates in that chain and add the next query point by the same criterion.

  v2, 2026-10-30: the v1 ordering ("tightest pair first") was measured at >300 ms/chain, i.e.
  >27 h single-core. It was backwards -- in a dense protein the SHORT-distance shell is the most
  heavily populated (every residue has near neighbours), so anchoring on 4.57 A maximised the
  branching factor.

Local CPU only. No GPU, no accelerator, no network.
"""
import argparse
import itertools
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
# The guidepost constellation file is a design input and is NOT distributed with this
# directory (see ../methods/README.md, "Reproducibility"). Point --pins at your own.
DEFAULT_PINS = os.path.join(REPO, "inputs", "guidepost_constellation.pdb")


# ----------------------------------------------------------------- build

def build(src, out_prefix):
    """Stream a multimodel CA-only PDB into coords.npy + offsets.npy + names.txt."""
    t0 = time.time()
    coords, offsets, names = [], [0], []
    cur = 0
    name = None
    with open(src, "rb") as fh:
        for raw in fh:
            tag = raw[:6]
            if tag == b"ATOM  ":
                coords.append((float(raw[30:38]), float(raw[38:46]), float(raw[46:54])))
                cur += 1
            elif tag == b"TITLE ":
                if raw[8:10] != b" 2":
                    name = raw[10:].decode("utf-8", "replace").strip().split()[0]
            elif tag == b"ENDMDL":
                offsets.append(cur)
                names.append(name or "?")
                name = None
    xyz = np.asarray(coords, dtype=np.float32)
    off = np.asarray(offsets, dtype=np.int64)
    np.save(out_prefix + "_coords.npy", xyz)
    np.save(out_prefix + "_offsets.npy", off)
    with open(out_prefix + "_names.txt", "w") as fh:
        fh.write("\n".join(names))
    print(f"[build] chains {len(names)}  CA atoms {len(xyz)}  "
          f"{time.time() - t0:.1f} s  -> {out_prefix}_*", flush=True)


# ----------------------------------------------------------------- query

def read_pins(path, target_chain="A"):
    """CA of every non-target chain, in chain order. Chain is NAMED, never inferred."""
    pins = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "CA":
                ch = ln[21]
                if ch != target_chain:
                    pins[ch] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    ks = sorted(pins)
    return ks, np.asarray([pins[k] for k in ks], dtype=np.float64)


# ----------------------------------------------------------------- batched Kabsch

def kabsch_rmsd_batch(P, Q):
    """RMSD after optimal superposition for a batch P (N,k,3) onto Q (k,3).

    Proper rotations only -- proteins are chiral, so a reflection is not a match.
    """
    Pc = P - P.mean(axis=1, keepdims=True)
    Qc = Q - Q.mean(axis=0, keepdims=True)
    H = np.einsum("nki,kj->nij", Pc, Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("nij,njk->nik", U, Vt)))
    S = S.copy()
    S[:, -1] *= d
    num = (Pc ** 2).sum(axis=(1, 2)) + (Qc ** 2).sum() - 2.0 * S.sum(axis=1)
    return np.sqrt(np.maximum(num, 0.0) / P.shape[1])


# ----------------------------------------------------------------- matcher

def search_chain(X, qd, idx, tol_d, tol_rmsd, Q, cap):
    """All matches of the query subset in one chain.

    Returns (list_of_(rmsd, residue_indices), overflowed_bool).
    """
    L = len(X)
    k = len(idx)
    if L < k:
        return [], False

    Xf = X.astype(np.float32)
    D = np.sqrt(((Xf[:, None, :] - Xf[None, :, :]) ** 2).sum(-1))

    # --- cheap rejections before any enumeration
    need = [qd[a, b] for a, b in itertools.combinations(idx, 2)]
    if D.max() < max(need) - tol_d:
        return [], False
    shells = {}
    counts = {}
    for a, b in itertools.combinations(idx, 2):
        M = np.abs(D - qd[a, b]) <= tol_d
        np.fill_diagonal(M, False)
        if not M.any():
            return [], False          # a required separation does not occur at all
        shells[(a, b)] = shells[(b, a)] = M
        counts[(a, b)] = counts[(b, a)] = int(M.sum())

    # --- order: start from the least-populated shell, then add the point whose shells
    #     against the already-placed set are least populated
    first = min(itertools.combinations(idx, 2), key=lambda p: counts[p])
    order = [first[0], first[1]]
    rest = [i for i in idx if i not in order]
    while rest:
        nxt = min(rest, key=lambda c: sum(counts[(c, p)] for p in order))
        order.append(nxt)
        rest.remove(nxt)

    # --- grow partial assignments in batch
    A = np.argwhere(shells[(order[0], order[1])])          # (m, 2)
    if len(A) == 0:
        return [], False
    if len(A) > cap:
        return [], True

    for t in range(2, k):
        iq = order[t]
        # mask: (L, n_partial) -- candidate row c is admissible for partial column n
        mask = np.ones((L, len(A)), dtype=bool)
        for pos, p in enumerate(order[:t]):
            mask &= shells[(iq, p)][:, A[:, pos]]
        # forbid reusing a residue already in the partial assignment
        for pos in range(t):
            mask[A[:, pos], np.arange(len(A))] = False
        rows, cols = np.nonzero(mask)
        if len(rows) == 0:
            return [], False
        if len(rows) > cap:
            return [], True
        A = np.concatenate([A[cols], rows[:, None]], axis=1)

    # --- final geometric gate
    P = Xf[A].astype(np.float64)                            # (n, k, 3)
    out = []
    CH = 20000
    qsub = Q[order]
    for s in range(0, len(P), CH):
        r = kabsch_rmsd_batch(P[s:s + CH], qsub)
        keep = np.flatnonzero(r <= tol_rmsd)
        for j in keep:
            out.append((float(r[j]), tuple(int(v) for v in A[s + j])))
    return out, False


# ----------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", nargs=2, metavar=("SRC_PDB", "OUT_PREFIX"))
    ap.add_argument("--search", metavar="STORE_PREFIX")
    ap.add_argument("--pins", default=DEFAULT_PINS)
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--subset", default="all", help="'all' or comma list of pin letters")
    ap.add_argument("--tol-d", type=float, default=1.0)
    ap.add_argument("--tol-rmsd", type=float, default=1.5)
    ap.add_argument("--limit", type=int, default=0, help="chains to scan (timing / sharding)")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--cap", type=int, default=2000000)
    ap.add_argument("--out", default="")
    ap.add_argument("--progress", type=int, default=0)
    args = ap.parse_args()

    if args.build:
        build(*args.build)
        return
    if not args.search:
        ap.error("need --build or --search")

    ks, Q = read_pins(args.pins, args.target_chain)
    qd = np.sqrt(((Q[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
    if args.subset == "all":
        idx = list(range(len(ks)))
    else:
        want = [s.strip() for s in args.subset.split(",")]
        idx = [ks.index(w) for w in want]
    print(f"[query] pins {ks}  subset {[ks[i] for i in idx]}  "
          f"tol_d {args.tol_d}  tol_rmsd {args.tol_rmsd}", flush=True)

    xyz = np.load(args.search + "_coords.npy", mmap_mode="r")
    off = np.load(args.search + "_offsets.npy")
    names = open(args.search + "_names.txt").read().split("\n")
    n_chains = len(off) - 1
    lo = args.start
    hi = n_chains if not args.limit else min(n_chains, lo + args.limit)

    fh = open(args.out, "w") if args.out else sys.stdout
    if args.out:
        fh.write("chain\trmsd\tresidue_indices\tchain_len\n")

    t0 = time.time()
    n_hit = n_ovf = n_chain_hit = 0
    for c in range(lo, hi):
        s, e = off[c], off[c + 1]
        res, ovf = search_chain(np.asarray(xyz[s:e]), qd, idx,
                                args.tol_d, args.tol_rmsd, Q, args.cap)
        if ovf:
            n_ovf += 1
            continue
        if res:
            n_chain_hit += 1
        for r, ridx in res:
            n_hit += 1
            if args.out:
                fh.write(f"{names[c]}\t{r:.3f}\t{','.join(map(str, ridx))}\t{e - s}\n")
        if args.progress and (c - lo + 1) % args.progress == 0:
            el = time.time() - t0
            print(f"  .. {c - lo + 1}/{hi - lo}  {el:.0f}s  "
                  f"{1000 * el / (c - lo + 1):.2f} ms/chain  hits {n_hit}", flush=True)
    dt = time.time() - t0
    if args.out:
        fh.close()
    done = hi - lo
    print(f"[search] chains {done}  chains_with_hit {n_chain_hit}  hits {n_hit}  "
          f"overflow {n_ovf}  {dt:.2f} s  {1000 * dt / max(1, done):.2f} ms/chain  "
          f"=> all {n_chains} ~ {dt / max(1, done) * n_chains / 60:.1f} min single-core",
          flush=True)


if __name__ == "__main__":
    main()
