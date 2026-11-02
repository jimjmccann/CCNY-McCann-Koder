#!/usr/bin/env python3
"""j892 - THE FACE GATE: do a hit's six matched groups sit on one accessible face?

WHY
---
The hit list named the six-pin hits (795 PDB-derived + 519 predicted-model chains) but
flagged the gate it could not
answer: a constellation BURIED in a 500-aa core is not a binding site. Until that is
measured, a "hit" is a geometric coincidence and not a candidate scaffold. This decides it.

WHAT IS MEASURED, per hit chain, on the chain IN ISOLATION (a scaffold is used standalone)
------------------------------------------------------------------------------------------
A. EXPOSURE  -- relative solvent accessibility (RSA) of each of the 6 matched residues,
   from freesasa (the trusted instrument; no hand-rolled SASA). A = how many of 6 have
   RSA >= 0.20, the conventional "exposed" line.
B. CO-ACCESSIBILITY -- can a partner reach all six FROM ONE DIRECTION? For each of 642
   Fibonacci-sphere directions u, a group is "visible" along u if no chain heavy atom lies
   within 3.0 A of the outward ray between t = 1.0 and t = 15.0 A. Atoms of the group's OWN
   residue are excluded (self-occlusion is spurious); the other five pins' atoms are KEPT,
   because a partner really is blocked by them. B = max over u of the visible count.
C. descriptive -- chain length, best_rmsd (carried from the hit list), max pairwise angle between
   the six outward vectors (group centroid - chain centroid).

PRE-REGISTERED PASS RULE (fixed before the run; no tuning after seeing output)
-----------------------------------------------------------------------------
  CANDIDATE : A >= 5 and B == 6
  PARTIAL   : A >= 4 and B >= 5   (and not CANDIDATE)
  FAIL      : otherwise

PRE-REGISTERED PREDICTIONS -- recorded so the result can refute them
-------------------------------------------------------------------
  P1  overall CANDIDATE rate < 15 %.  ** If it is > 40 % that is GOOD news and must be
      reported as such, not quietly re-framed. **
  P2  CANDIDATE rate among 121-250 aa chains > CANDIDATE rate among > 500 aa chains.
  P3  LOOKUP_FAIL < 5 % of hits.

HARD CORRECTNESS ASSERTION (the residue-type verification guard)
------------------------------------------------------------------------
The residue found at a recorded residue number MUST match the recorded pin type:
CATION => ARG or LYS, IMID => HIS. A mismatch means the lookup landed on the wrong
residue, so the chain is recorded LOOKUP_FAIL and EXCLUDED from every statistic rather
than silently measured. ** If LOOKUP_FAIL >= 5 % the face numbers must NOT be published --
the identity mapping is broken and nothing downstream is trustworthy. **

Read path is j851's own (gemmi read_structure, remove_alternative_conformations,
remove_hydrogens; chain by auth name; residue by seqid.num) so the residue numbers in the
hit-list TSV mean here exactly what they meant when the store was built.

Local CPU only. No GPU, no accelerator, no network.
"""
import argparse
import gzip
import math
import os
import shutil
import tempfile

import numpy as np

CATION_RES = {"ARG", "LYS"}
IMID_RES = {"HIS"}
ARG_AT = ("NE", "NH1", "NH2")
HIS_AT = ("ND1", "NE2", "CG")
RSA_EXPOSED = 0.20
RAY_RADIUS = 3.0
RAY_TMIN, RAY_TMAX = 1.0, 15.0
NDIR = 642


def fib_sphere(n):
    i = np.arange(n, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    gold = math.pi * (1.0 + 5.0 ** 0.5)
    theta = gold * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


DIRS = fib_sphere(NDIR)


def split_chain(chain_full):
    """'4w5t_A' -> ('4w5t','A') ; 'AF-Q9I4A6-F1-model_v6_A' -> (...'_v6','A').
    rsplit, never split: an AFDB stem CONTAINS an underscore and split() would cut it
    in the wrong place and then look up the wrong chain."""
    return chain_full.rsplit("_", 1)


def mirror_path(mirror, stem, layout):
    if layout == "flat":
        return os.path.join(mirror, f"{stem}.cif.gz")
    pid = stem.lower()
    return os.path.join(mirror, pid[1:3], f"{pid}.cif.gz")


def read_chain(path, chain_name):
    """j851's read path. Returns (residues, ok) where residues is a list of
    (resname, seqid, {atom: xyz}) for amino acids of the named chain only."""
    import gemmi
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            doc = gemmi.cif.read_string(fh.read())
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        st = gemmi.read_structure(path)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        return None
    for chain in st[0]:
        if chain.name != chain_name:
            continue
        out = []
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            at = {a.name: (a.pos.x, a.pos.y, a.pos.z) for a in res}
            if "CA" not in at:
                continue
            out.append((res.name, res.seqid.num, at))
        return out
    return None


def group_centroid(resname, at):
    if resname == "ARG" and all(a in at for a in ARG_AT):
        return np.mean([at[a] for a in ARG_AT], axis=0)
    if resname == "LYS" and "NZ" in at:
        return np.asarray(at["NZ"], dtype=float)
    if resname == "HIS" and all(a in at for a in HIS_AT):
        return np.mean([at[a] for a in HIS_AT], axis=0)
    return None


def write_pdb(residues, tmp):
    """Isolated chain as PDB: chain 'A', sequential numbering (kills the >9999 and
    multi-char-chain format limits). Returns seqid -> new number map."""
    remap = {}
    with open(tmp, "w") as fh:
        n = 0
        serial = 1
        for resname, seqid, at in residues:
            n += 1
            remap[seqid] = n
            for aname, (x, y, z) in at.items():
                el = aname[0]
                if el.isdigit():
                    el = aname[1] if len(aname) > 1 else "C"
                fh.write(f"ATOM  {serial:>5} {aname:<4}{resname:>4} A{n:>4}    "
                         f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          "
                         f"{el:>2}\n")
                serial += 1
        fh.write("END\n")
    return remap


def rsa_of(tmp, remap, want):
    import freesasa
    st = freesasa.Structure(tmp)
    res = freesasa.calc(st)
    areas = res.residueAreas()
    out = {}
    ch = areas.get("A", {})
    for seqid in want:
        a = ch.get(str(remap[seqid]))
        out[seqid] = (None if a is None or a.relativeTotal is None
                      else float(a.relativeTotal))
    return out


def co_access(allxyz, owner, cents, pin_seqids):
    """max over directions of how many pin centroids are unobstructed."""
    best = 0
    masks = []
    for j, g in enumerate(cents):
        d = allxyz - g
        near = np.linalg.norm(d, axis=1) <= (RAY_TMAX + RAY_RADIUS)
        keep = near & (owner != pin_seqids[j])
        P = d[keep]
        if len(P) == 0:
            masks.append(np.ones(NDIR, dtype=bool))
            continue
        t = P @ DIRS.T                                  # (atoms, dirs)
        perp2 = (P * P).sum(1)[:, None] - t * t
        block = (t >= RAY_TMIN) & (t <= RAY_TMAX) & (perp2 <= RAY_RADIUS ** 2)
        masks.append(~block.any(axis=0))
    vis = np.stack(masks, axis=0).sum(axis=0)
    best = int(vis.max())
    return best, int(np.argmax(vis))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--mirror", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--layout", choices=("pdbmirror", "flat"), default="pdbmirror")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--poscontrol", action="store_true",
                    help="VACUITY GATE: replace the matched residues with the SIX "
                         "MOST-EXPOSED ARG/LYS/HIS of the same chain. If this does not "
                         "produce CANDIDATEs, the pass branch is unreachable and every "
                         "FAIL in the real run is an artefact, not a result.")
    args = ap.parse_args()

    rows = [l.rstrip("\n").split("\t") for l in open(args.hits)]
    hdr, rows = rows[0], rows[1:]
    if args.limit:
        rows = rows[:args.limit]
    print(f"[hits] {len(rows)} from {args.hits}")

    tmpdir = tempfile.mkdtemp(prefix="j892_")
    tmp = os.path.join(tmpdir, "chain.pdb")
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"j892_face_gate_{args.tag}_20261031.tsv")
    nfail = 0
    recs = []

    for i, r in enumerate(rows, 1):
        chain_full, rmsd = r[1], float(r[2])
        resids = [int(v) for v in r[4:10]]
        types = r[10:16]
        stem, chain_name = split_chain(chain_full)
        path = mirror_path(args.mirror, stem, args.layout)
        verdict = None
        if not os.path.exists(path):
            verdict = "LOOKUP_FAIL_NOFILE"
        else:
            residues = read_chain(path, chain_name)
            if not residues:
                verdict = "LOOKUP_FAIL_NOCHAIN"
        if verdict is None:
            by_id = {}
            for resname, seqid, at in residues:
                by_id.setdefault(seqid, (resname, at))
            cents, ok = [], True
            for seqid, ty in zip(resids, types):
                if seqid not in by_id:
                    ok = False
                    break
                resname, at = by_id[seqid]
                want = CATION_RES if ty == "CATION" else IMID_RES
                if resname not in want:
                    ok = False
                    break
                g = group_centroid(resname, at)
                if g is None:
                    ok = False
                    break
                cents.append(g)
            if not ok:
                verdict = "LOOKUP_FAIL_TYPE"
        if verdict is not None:
            nfail += 1
            recs.append((chain_full, rmsd, 0, verdict, -1, -1, -1.0, -1.0))
            continue

        cents = np.asarray(cents, dtype=np.float64)
        allxyz, owner = [], []
        for resname, seqid, at in residues:
            for aname, xyz in at.items():
                allxyz.append(xyz)
                owner.append(seqid)
        allxyz = np.asarray(allxyz, dtype=np.float64)
        owner = np.asarray(owner, dtype=np.int64)
        nres = len(residues)

        remap = write_pdb(residues, tmp)
        if args.poscontrol:
            pool = [(seqid, resname, at) for resname, seqid, at in residues
                    if resname in CATION_RES | IMID_RES
                    and group_centroid(resname, at) is not None]
            try:
                prsa = rsa_of(tmp, remap, [p[0] for p in pool])
            except Exception:
                nfail += 1
                recs.append((chain_full, rmsd, nres, "LOOKUP_FAIL_SASA",
                             -1, -1, -1.0, -1.0))
                continue
            pool = [p for p in pool if prsa.get(p[0]) is not None]
            pool.sort(key=lambda p: -prsa[p[0]])
            if len(pool) < 6:
                nfail += 1
                recs.append((chain_full, rmsd, nres, "LOOKUP_FAIL_POOL",
                             -1, -1, -1.0, -1.0))
                continue
            pick = pool[:6]
            resids = [p[0] for p in pick]
            cents = np.asarray([group_centroid(p[1], p[2]) for p in pick],
                               dtype=np.float64)
        try:
            rsa = rsa_of(tmp, remap, resids)
        except Exception:
            nfail += 1
            recs.append((chain_full, rmsd, nres, "LOOKUP_FAIL_SASA", -1, -1, -1.0, -1.0))
            continue
        A = sum(1 for s in resids if rsa.get(s) is not None and rsa[s] >= RSA_EXPOSED)
        rsa_min = min((rsa[s] for s in resids if rsa.get(s) is not None), default=-1.0)
        B, _ = co_access(allxyz, owner, cents, resids)

        com = allxyz.mean(axis=0)
        v = cents - com
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        cosm = np.clip(v @ v.T, -1, 1)
        maxang = float(np.degrees(np.arccos(cosm)).max())

        if A >= 5 and B == 6:
            verdict = "CANDIDATE"
        elif A >= 4 and B >= 5:
            verdict = "PARTIAL"
        else:
            verdict = "FAIL"
        recs.append((chain_full, rmsd, nres, verdict, A, B, rsa_min, maxang))
        if i % 100 == 0:
            print(f"  .. {i}/{len(rows)}")

    with open(out, "w") as fh:
        fh.write("chain\tbest_rmsd\tchain_len\tverdict\tA_exposed_of6\t"
                 "B_coaccess_of6\tmin_rsa\tmax_outward_angle_deg\n")
        for t in recs:
            fh.write(f"{t[0]}\t{t[1]:.4f}\t{t[2]}\t{t[3]}\t{t[4]}\t{t[5]}\t"
                     f"{t[6]:.4f}\t{t[7]:.2f}\n")
    shutil.rmtree(tmpdir, ignore_errors=True)

    n = len(recs)
    good = [t for t in recs if t[3] in ("CANDIDATE", "PARTIAL", "FAIL")]
    cand = [t for t in recs if t[3] == "CANDIDATE"]
    part = [t for t in recs if t[3] == "PARTIAL"]
    print(f"\n[out] {out}")
    print(f"== {args.tag}: {n} hits ==")
    print(f"  LOOKUP_FAIL (excluded) : {nfail}  ({100.0*nfail/n:.1f} %)   "
          f"{'** >=5 % -- DO NOT PUBLISH THE FACE NUMBERS **' if nfail/n >= 0.05 else 'P3 ok'}")
    if not good:
        return
    print(f"  measured               : {len(good)}")
    print(f"  CANDIDATE              : {len(cand)}  ({100.0*len(cand)/len(good):.1f} % of measured)")
    print(f"  PARTIAL                : {len(part)}  ({100.0*len(part)/len(good):.1f} %)")
    print(f"  FAIL                   : {len(good)-len(cand)-len(part)}")
    rate = len(cand) / len(good)
    print(f"  P1 (<15 % candidate)   : {'HELD' if rate < 0.15 else 'REFUTED'}"
          f"{'  ** >40 %: GOOD NEWS **' if rate > 0.40 else ''}")
    for lo, hi, lab in [(0, 120, '<=120'), (121, 250, '121-250'),
                        (251, 500, '251-500'), (501, 10**9, '>500')]:
        band = [t for t in good if lo <= t[2] <= hi]
        if band:
            c = sum(1 for t in band if t[3] == "CANDIDATE")
            print(f"    len {lab:>8}: {len(band):>4} measured, {c:>3} candidate "
                  f"({100.0*c/len(band):.1f} %)")
    if cand:
        print(f"\n== CANDIDATES, best rmsd first ==")
        print(f"  {'chain':<16} {'rmsd':>7} {'len':>6} {'A':>3} {'B':>3} {'minRSA':>7}")
        for t in sorted(cand, key=lambda x: x[1])[:40]:
            print(f"  {t[0]:<16} {t[1]:>7.4f} {t[2]:>6} {t[4]:>3} {t[5]:>3} {t[6]:>7.3f}")


if __name__ == "__main__":
    main()
