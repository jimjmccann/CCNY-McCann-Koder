#!/usr/bin/env python3
"""Measure LOCAL and GLOBAL identity for phmmer's hits, so the gate verdict
covers both tools.

phmmer's tabular output carries no identity column at all, so identity has to
be recomputed.  Every pair is realigned with Smith-Waterman (BLOSUM62,
gap -11/-1) and scored as:
    LOCAL  = nident / alnlen            (of the local alignment)
    GLOBAL = nident / full query length (same local alignment)
GLOBAL is deliberately derived from the LOCAL alignment: forcing a
Needleman-Wunsch alignment against a much longer target scatters identities
over the whole target and inflates nident/qlen (measured: 30.7% -> 40.0% for
id3017 vs 7udk_A), which is not what "sequence similarity" means.

Results are reported per significance tier, because phmmer calls repetitive
targets significant by accumulating many individually-insignificant repeat
matches -- so a full-sequence E-value alone over-calls homology.
"""
import csv
import os
from collections import defaultdict
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

REPO = Path("/PATH/TO/CHECKOUT")
DATA = REPO / "data" / "novelty_20261006"
SCRATCH = Path(os.environ["DESIGN_DATA_ROOT"])
WORK = SCRATCH / "novelty_search_20261006"
DBDIR = SCRATCH / "novelty_db_20261006"
DBS = {"swissprot": DBDIR / "uniprot_sprot.fasta",
       "pdbseqres": DBDIR / "pdb_seqres_prot.fasta"}
TOP_PH = 200        # phmmer hits per query per db to realign
SIG = 1e-3


def read_fasta(path):
    n, b = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    yield n, "".join(b)
                n, b = line[1:].rstrip("\n"), []
            else:
                b.append(line.strip())
    if n:
        yield n, "".join(b)


def norm_acc(name):
    p = name.split("|")
    return p[1] if len(p) >= 3 and p[0] in ("sp", "tr") else name


def main():
    order, tiers, qseqs = [], {}, {}
    for hdr, seq in read_fasta(DATA / "novelty_seq_query_20261006.fasta"):
        parts = hdr.split()
        order.append(parts[0])
        qseqs[parts[0]] = seq
        tiers[parts[0]] = next((p.split("=", 1)[1] for p in parts
                                if p.startswith("struct_tier=")), "?")

    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -11.0, -1.0
    al.mode = "local"

    rows = []
    for db in DBS:
        # ---- parse phmmer: best full E-value per (query, target) ---------
        best = defaultdict(dict)
        dompath = WORK / f"phmmer_{db}_domtbl_20261006.txt"
        with open(dompath) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                if len(f) < 22:
                    continue
                t, q = f[0], f[3]
                try:
                    fe, de = float(f[6]), float(f[12])
                except ValueError:
                    continue
                cur = best[q].get(t)
                # keep the full E-value and the BEST (smallest) domain E-value
                if cur is None or fe < cur[0]:
                    best[q][t] = (fe, de)
                elif de < cur[1]:
                    best[q][t] = (cur[0], de)

        want = set()
        shortlist = {}
        for q in order:
            hits = sorted(best[q].items(), key=lambda kv: kv[1][0])[:TOP_PH]
            shortlist[q] = hits
            want.update(t for t, _ in hits)

        tseq = {}
        for hdr, seq in read_fasta(DBS[db]):
            k = hdr.split()[0]
            if k in want:
                tseq[k] = seq
        print(f"[{db}] realigning {sum(len(v) for v in shortlist.values())} "
              f"pairs; resolved {len(tseq)}/{len(want)} target seqs")

        for q in order:
            qs = qseqs[q]
            # (local, global, nident, alnlen, target, fullE, domE)
            measured = []
            for t, (fe, de) in shortlist[q]:
                ts = tseq.get(t)
                if ts is None:
                    continue
                a = al.align(qs, ts)[0]
                x, y = str(a[0]), str(a[1])
                nid = sum(1 for p, r in zip(x, y) if p == r and p != "-")
                measured.append((nid / len(x), nid / len(qs), nid, len(x),
                                 t, fe, de))
            if not measured:
                continue

            def mx(key, filt=None):
                c = [m for m in measured if (filt is None or filt(m))]
                return max(c, key=lambda m: m[key]) if c else None

            sig_full = lambda m: m[5] <= SIG
            sig_dom = lambda m: m[6] <= SIG
            views = {
                "all_maxglobal": mx(1),
                "all_maxlocal": mx(0),
                "sigfull_maxglobal": mx(1, sig_full),
                "sigfull_maxlocal": mx(0, sig_full),
                "sigdom_maxglobal": mx(1, sig_dom),
                "sigdom_maxlocal": mx(0, sig_dom),
            }
            print(f"\n-- {q} ({tiers[q]}) [{db}]  "
                  f"{len(measured)} realigned")
            for tag, m in views.items():
                if m is None:
                    print(f"   {tag:20s}: none")
                    continue
                print(f"   {tag:20s}: {m[0]*100:5.1f}% local / "
                      f"{m[1]*100:5.1f}% global  ({m[2]}/{m[3]} aln, "
                      f"{m[2]}/{len(qs)} q)  {m[4]:10s} "
                      f"fullE={m[5]:.2g} domE={m[6]:.2g}")
                rows.append({"db": db, "query": q, "struct_tier": tiers[q],
                             "view": tag, "local_pct": round(m[0]*100, 1),
                             "global_pct": round(m[1]*100, 1), "nident": m[2],
                             "alnlen": m[3], "qlen": len(qs), "target": m[4],
                             "phmmer_full_evalue": f"{m[5]:.3g}",
                             "phmmer_best_dom_evalue": f"{m[6]:.3g}"})

    out = DATA / "novelty_seq_phmmer_identity_20261006.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
