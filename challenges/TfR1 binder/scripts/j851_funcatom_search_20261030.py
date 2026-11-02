#!/usr/bin/env python3
"""j851 — constellation search on FUNCTIONAL ATOMS, not on a synthetic backbone.

WHY THIS EXISTS. j846 then j848 searched on CA, then on CA+CB, and both returned ~0
informative signal. Two independent reviews found the cause: the guidepost CA/CB coordinates are
SYNTHETIC -- every Arg pin carries chi1 -176.4 / chi2 177.2 and every His pin chi1 -66.9 /
chi2 90.2, spread < 0.07 deg, i.e. one idealised rotamer emitted by our own builder. The
RFd3 spec pins only the FUNCTIONAL TIP ATOMS (select_fixed_atoms: NE,NH1,NH2 for Arg;
ND1,NE2,CG for His). The backbone was never part of the specification, so matching on it
was matching our own tooling's default. Measured consequence: superposing a guanidinium
exactly on real PDB Arg, two CA agree within 1.0 A only 4.5% of the time -- so a perfect
chemical match fails a CA gate essentially always.

WHAT THIS DOES INSTEAD. Each pin becomes a typed functional GROUP:
    CATION  Arg guanidinium (NE,NH1,NH2 centroid)  or  Lys NZ
    IMID    His imidazole   (ND1,NE2,CG centroid)
A chain's candidate points are its own functional groups -- roughly 30 per 250-residue
chain instead of 250 residues, which is ~5 orders of magnitude fewer 6-subsets. Type
compatibility is enforced, so a cation pin can only be served by a cation group. The
backbone is then free: where CA sits is an OUTPUT, not a constraint.

This is the "position not residue" query the line was built for -- position of the CHEMISTRY.

Local CPU only. No GPU and no accelerator.
"""
import argparse
import glob
import gzip
import io
import itertools
import os
import sys
import tarfile
import time
from multiprocessing import Pool

import numpy as np

CATION, IMID = 0, 1
TYPENAME = {CATION: "CATION", IMID: "IMID"}

ARG_AT = ("NE", "NH1", "NH2")
HIS_AT = ("ND1", "NE2", "CG")


# ----------------------------------------------------------------- extraction

def groups_from_residues(reslist):
    """reslist: (resname, resseq, {atomname: xyz}, plddt) -> arrays of groups."""
    cen, typ, ca, rid = [], [], [], []
    for resn, seq, at, _ in reslist:
        if "CA" not in at:
            continue
        if resn == "ARG" and all(a in at for a in ARG_AT):
            cen.append(np.mean([at[a] for a in ARG_AT], axis=0))
            typ.append(CATION)
        elif resn == "LYS" and "NZ" in at:
            cen.append(np.asarray(at["NZ"], dtype=float))
            typ.append(CATION)
        elif resn == "HIS" and all(a in at for a in HIS_AT):
            cen.append(np.mean([at[a] for a in HIS_AT], axis=0))
            typ.append(IMID)
        else:
            continue
        ca.append(at["CA"])
        rid.append(seq)
    if not cen:
        return None
    return (np.asarray(cen, dtype=np.float32), np.asarray(typ, dtype=np.int8),
            np.asarray(ca, dtype=np.float32), np.asarray(rid, dtype=np.int32))


def read_cif_gemmi(path_or_bytes, name, plddt_min):
    import gemmi
    if isinstance(path_or_bytes, bytes):
        # gemmi 0.7.5 has no read_structure_from_memory -- go via the cif document.
        doc = gemmi.cif.read_string(path_or_bytes.decode())
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        st = gemmi.read_structure(path_or_bytes)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    out = []
    if len(st) == 0:
        return out
    for chain in st[0]:
        reslist = []
        for res in chain:
            at = {}
            bs = []
            for a in res:
                at[a.name] = (a.pos.x, a.pos.y, a.pos.z)
                bs.append(a.b_iso)
            if plddt_min > 0 and bs and np.mean(bs) < plddt_min:
                continue          # AFDB: b-factor column IS pLDDT
            reslist.append((res.name, res.seqid.num, at, np.mean(bs) if bs else 0.0))
        g = groups_from_residues(reslist)
        if g is not None and len(g[0]) >= 4:
            out.append((f"{name}_{chain.name}",) + g)
    return out


def _one_file(a):
    path, plddt_min = a
    try:
        return read_cif_gemmi(path, os.path.basename(path).split(".")[0], plddt_min)
    except Exception:
        return []


def _one_tar(a):
    tpath, plddt_min = a
    out = []
    try:
        with tarfile.open(tpath) as tf:
            for m in tf:
                if not m.name.endswith(".cif.gz"):
                    continue
                f = tf.extractfile(m)
                if f is None:
                    continue
                raw = gzip.decompress(f.read())
                nm = os.path.basename(m.name).split(".")[0]
                try:
                    out.extend(read_cif_gemmi(raw, nm, plddt_min))
                except Exception:
                    pass
    except Exception:
        pass
    return out


def build(inputs, out_prefix, nproc, plddt_min):
    tars = [p for p in inputs if p.endswith(".tar")]
    files = []
    for d in inputs:
        if os.path.isdir(d):
            files.extend(sorted(glob.glob(os.path.join(d, "*.cif.gz"))))
    jobs = [(t, plddt_min) for t in tars]
    print(f"[build] {len(tars)} tar archives, {len(files)} loose files, "
          f"plddt_min={plddt_min}", flush=True)
    t0 = time.time()
    CEN, TYP, CA, RID, off, names = [], [], [], [], [0], []
    cur = 0

    def absorb(res):
        nonlocal cur
        for nm, cen, typ, ca, rid in res:
            CEN.append(cen); TYP.append(typ); CA.append(ca); RID.append(rid)
            cur += len(cen)
            off.append(cur)
            names.append(nm)

    with Pool(nproc) as pool:
        if jobs:
            for i, res in enumerate(pool.imap_unordered(_one_tar, jobs)):
                absorb(res)
                print(f"  .. tar {i + 1}/{len(jobs)}  chains {len(names)}  "
                      f"{time.time() - t0:.0f}s", flush=True)
        if files:
            for i, res in enumerate(pool.imap_unordered(
                    _one_file, [(f, plddt_min) for f in files], chunksize=32)):
                absorb(res)
                if (i + 1) % 5000 == 0:
                    print(f"  .. file {i + 1}/{len(files)}  chains {len(names)}  "
                          f"{time.time() - t0:.0f}s", flush=True)
    np.save(out_prefix + "_cen.npy", np.concatenate(CEN))
    np.save(out_prefix + "_typ.npy", np.concatenate(TYP))
    np.save(out_prefix + "_ca.npy", np.concatenate(CA))
    np.save(out_prefix + "_rid.npy", np.concatenate(RID))
    np.save(out_prefix + "_offsets.npy", np.asarray(off, dtype=np.int64))
    open(out_prefix + "_names.txt", "w").write("\n".join(names))
    print(f"[build] chains {len(names)}  functional groups {cur}  "
          f"({cur / max(1, len(names)):.1f} per chain)  {time.time() - t0:.0f}s", flush=True)


# ----------------------------------------------------------------- query

def read_query(pdb, spec_json):
    """Functional centroid + type per pin, from the SPEC's own select_fixed_atoms."""
    import json
    spec = json.load(open(spec_json))
    spec = spec[list(spec)[0]]
    sel = spec["select_fixed_atoms"]
    at = {}
    for ln in open(pdb):
        if ln.startswith(("ATOM", "HETATM")):
            ch = ln[21]
            at.setdefault(ch, {})[ln[12:16].strip()] = (
                float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    cen, typ, labels = [], [], []
    for key in sorted(sel):
        ch = key[0]
        atoms = [a.strip() for a in sel[key].split(",")]
        pts = [at[ch][a] for a in atoms]
        cen.append(np.mean(pts, axis=0))
        typ.append(IMID if "ND1" in atoms else CATION)
        labels.append(f"{key}:{TYPENAME[typ[-1]]}")
    return labels, np.asarray(cen, dtype=np.float64), np.asarray(typ, dtype=np.int8)


# ----------------------------------------------------------------- matcher

def kabsch_rmsd_batch(P, Q):
    Pc = P - P.mean(axis=1, keepdims=True)
    Qc = Q - Q.mean(axis=0, keepdims=True)
    H = np.einsum("nki,kj->nij", Pc, Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("nij,njk->nik", U, Vt)))
    S = S.copy()
    S[:, -1] *= d
    num = (Pc ** 2).sum(axis=(1, 2)) + (Qc ** 2).sum() - 2.0 * S.sum(axis=1)
    return np.sqrt(np.maximum(num, 0.0) / P.shape[1])


def search_chain(cen, typ, qd, qtyp, tol_d, tol_rmsd, Q, cap):
    L = len(cen)
    k = len(qd)
    if L < k:
        return []
    D = np.sqrt(((cen[:, None, :] - cen[None, :, :]) ** 2).sum(-1))
    tok = [np.flatnonzero(typ == qtyp[i]) for i in range(k)]
    if any(len(t) == 0 for t in tok):
        return []

    shells, counts = {}, {}
    for a, b in itertools.combinations(range(k), 2):
        M = np.abs(D - qd[a, b]) <= tol_d
        np.fill_diagonal(M, False)
        Ma = np.zeros_like(M)
        Ma[np.ix_(tok[a], tok[b])] = M[np.ix_(tok[a], tok[b])]
        if not Ma.any():
            return []
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

    inv = np.argsort(np.asarray(order))
    A = A[:, inv]
    out = []
    for s in range(0, len(A), 20000):
        blk = A[s:s + 20000]
        r = kabsch_rmsd_batch(cen[blk].astype(np.float64), Q)
        for j in np.flatnonzero(r <= tol_rmsd):
            out.append((float(r[j]), tuple(int(v) for v in blk[j])))
    return out


def _search_shard(a):
    (prefix, qd, qtyp, Q, tol_d, tol_rmsd, cap, lo, hi) = a
    cen = np.load(prefix + "_cen.npy", mmap_mode="r")
    typ = np.load(prefix + "_typ.npy", mmap_mode="r")
    off = np.load(prefix + "_offsets.npy")
    nh = nc = 0
    rm = []
    hits = []
    for c in range(lo, hi):
        s, e = off[c], off[c + 1]
        res = search_chain(np.asarray(cen[s:e]), np.asarray(typ[s:e]),
                           qd, qtyp, tol_d, tol_rmsd, Q, cap)
        if res:
            nc += 1
            nh += len(res)
            rm.extend(r for r, _ in res)
            hits.append((c, min(res)[0]))
    return nh, nc, np.asarray(rm), hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", nargs="+")
    ap.add_argument("--out-prefix", default="")
    ap.add_argument("--plddt-min", type=float, default=0.0)
    ap.add_argument("--search")
    ap.add_argument("--pins")
    ap.add_argument("--spec")
    ap.add_argument("--tol-d", type=float, default=1.5)
    ap.add_argument("--tol-rmsd", type=float, default=2.0)
    ap.add_argument("--cap", type=int, default=2000000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--nproc", type=int, default=16)
    ap.add_argument("--hits-out", default="")
    args = ap.parse_args()

    if args.build:
        build(args.build, args.out_prefix, args.nproc, args.plddt_min)
        return

    labels, Q, qtyp = read_query(args.pins, args.spec)
    qd = np.sqrt(((Q[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
    print(f"[query] {len(labels)} functional pins: {labels}")
    print("[query] pairwise functional-centroid distances (A):")
    for a, b in itertools.combinations(range(len(labels)), 2):
        print(f"    {labels[a]:>9} - {labels[b]:<9} {qd[a, b]:6.2f}")

    off = np.load(args.search + "_offsets.npy")
    n = len(off) - 1
    hi = n if not args.limit else min(n, args.limit)
    print(f"[search] store chains {n}, scanning {hi}", flush=True)

    bounds = np.linspace(0, hi, args.nproc + 1).astype(int)
    jobs = [(args.search, qd, qtyp, Q, args.tol_d, args.tol_rmsd, args.cap,
             int(bounds[i]), int(bounds[i + 1])) for i in range(args.nproc)]
    t0 = time.time()
    NH = NC = 0
    RM = []
    HITS = []
    with Pool(args.nproc) as pool:
        for nh, nc, rm, hits in pool.imap_unordered(_search_shard, jobs):
            NH += nh; NC += nc; RM.append(rm); HITS.extend(hits)
    RM = np.concatenate(RM) if RM else np.asarray([])
    print(f"[search] {time.time() - t0:.0f}s   hits {NH}   chains_with_hit {NC} / {hi}"
          f"  ({100.0 * NC / max(1, hi):.2f}%)")
    if len(RM):
        print(f"{'tol_rmsd':>9} {'hits':>10}")
        for t in [2.0, 1.5, 1.0, 0.75, 0.5, 0.35, 0.25]:
            print(f"{t:>9.2f} {int((RM <= t).sum()):>10d}")
        print(f"best rmsd {RM.min():.3f} A")
    if args.hits_out and HITS:
        names = open(args.search + "_names.txt").read().split("\n")
        with open(args.hits_out, "w") as fh:
            fh.write("chain\tbest_rmsd\n")
            for c, r in sorted(HITS, key=lambda x: x[1]):
                fh.write(f"{names[c]}\t{r:.3f}\n")
        print(f"[search] wrote {args.hits_out}")


if __name__ == "__main__":
    main()
