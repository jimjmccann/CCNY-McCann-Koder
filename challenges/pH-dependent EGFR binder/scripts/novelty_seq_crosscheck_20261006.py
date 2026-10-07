#!/usr/bin/env python3
"""Cross-check the MMseqs2 sequence-novelty numbers two independent ways.

1. phmmer (--max, profile-based) vs the same two databases: does an
   INDEPENDENT tool find a significant homolog where MMseqs2 found none?
2. An independent Smith-Waterman realignment (Biopython PairwiseAligner,
   BLOSUM62, gap open -11 / extend -1 -- matching mmseqs' matrix and gap
   costs) of each query against the union of both tools' top hits, scoring
   LOCAL (nident/alnlen) and GLOBAL (nident/qlen) identity from scratch.

Purpose: the gate verdict rests on "max identity <= 30%", so the identity
numbers must not depend on one tool's conventions or one tool's aligner.
"""
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

REPO = Path("/PATH/TO/CHECKOUT")
DATA = REPO / "data" / "novelty_20261006"
QFA = DATA / "novelty_seq_query_20261006.fasta"
SCRATCH = Path(os.environ["DESIGN_DATA_ROOT"])
WORK = SCRATCH / "novelty_search_20261006"
DBDIR = SCRATCH / "novelty_db_20261006"

DBS = {"swissprot": DBDIR / "uniprot_sprot.fasta",
       "pdbseqres": DBDIR / "pdb_seqres_prot.fasta"}
TOP_N = 10          # top hits per query per tool to realign
SIG = 1e-3          # mmseqs' own default significance threshold


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def read_queries():
    order, tiers, seqs = [], {}, {}
    for hdr, seq in read_fasta(QFA):
        parts = hdr.split()
        q = parts[0]
        order.append(q)
        seqs[q] = seq
        tiers[q] = next((p.split("=", 1)[1] for p in parts
                         if p.startswith("struct_tier=")), "?")
    return order, tiers, seqs


def norm_acc(name):
    """Normalise a target id so mmseqs' and phmmer's names can be compared.
    mmseqs strips UniProt 'sp|ACC|ENTRY' down to ACC; phmmer keeps it whole."""
    m = re.match(r"^(?:sp|tr)\|([^|]+)\|", name)
    return m.group(1) if m else name


def parse_phmmer_domtbl(path):
    """-> {query: [(full_evalue, dom_ievalue, target, raw_target)]}"""
    out = defaultdict(list)
    if not path.exists():
        print(f"  [phmmer] MISSING {path}")
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            raw_t, q = f[0], f[3]
            try:
                full_e, dom_ie = float(f[6]), float(f[12])
            except ValueError:
                continue
            out[q].append((full_e, dom_ie, norm_acc(raw_t), raw_t))
    for q in out:
        out[q].sort(key=lambda r: r[0])
    return out


def parse_mmseqs(path):
    """-> {query: [row dicts]} sorted by descending bitscore"""
    out = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["evalue"] = float(r["evalue"])
            r["bits"] = float(r["bits"])
            r["nident"] = int(r["nident"])
            r["alnlen"] = int(r["alnlen"])
            r["qlen"] = int(r["qlen"])
            out[r["query"]].append(r)
    for q in out:
        out[q].sort(key=lambda r: -r["bits"])
    return out


def make_aligner():
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11.0      # mmseqs --gap-open aa:11
    al.extend_gap_score = -1.0     # mmseqs --gap-extend aa:1
    al.mode = "local"              # Smith-Waterman
    return al


def sw_identity(aligner, qseq, tseq):
    """-> (local_id, global_id, nident, alnlen) from an independent SW align."""
    try:
        aln = aligner.align(qseq, tseq)[0]
    except Exception as e:                      # pragma: no cover
        print(f"    [sw] align failed: {e}")
        return None
    a, b = str(aln[0]), str(aln[1])
    nid = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    alnlen = len(a)
    return nid / alnlen, nid / len(qseq), nid, alnlen


def main():
    order, tiers, qseqs = read_queries()
    aligner = make_aligner()

    ph = {db: parse_phmmer_domtbl(
              WORK / f"phmmer_{db}_domtbl_20261006.txt") for db in DBS}
    mm = {db: parse_mmseqs(
              DATA / f"mmseqs_{db}_hits_20261006.tsv") for db in DBS}

    # ---- which target sequences do we need? -----------------------------
    needed = {db: set() for db in DBS}
    for db in DBS:
        for q in order:
            for r in mm[db].get(q, [])[:TOP_N]:
                needed[db].add(r["target"])
            for _, _, acc, raw in ph[db].get(q, [])[:TOP_N]:
                needed[db].add(acc)

    # ---- load just those sequences --------------------------------------
    tseqs = {db: {} for db in DBS}
    for db, path in DBS.items():
        want = needed[db]
        for hdr, seq in read_fasta(path):
            first = hdr.split()[0]
            for key in (first, norm_acc(first)):
                if key in want:
                    tseqs[db][key] = seq
        missing = want - set(tseqs[db])
        print(f"  [seqs] {db}: resolved {len(tseqs[db])}/{len(want)} targets"
              + (f"  MISSING: {sorted(missing)[:5]}" if missing else ""))

    rows = []
    agree_delta = []
    for db in DBS:
        print(f"\n################ {db} -- cross-check ################")
        for q in order:
            qs = qseqs[q]
            mrows = mm[db].get(q, [])
            prows = ph[db].get(q, [])
            mm_sig = [r for r in mrows if r["evalue"] <= SIG]
            ph_sig = [r for r in prows if r[0] <= SIG]

            # independent SW over the union of both tools' top hits
            cands = {r["target"] for r in mrows[:TOP_N]}
            cands |= {acc for _, _, acc, _ in prows[:TOP_N]}
            best_loc = best_glob = None
            for t in sorted(cands):
                ts = tseqs[db].get(t)
                if ts is None:
                    continue
                res = sw_identity(aligner, qs, ts)
                if res is None:
                    continue
                loc, glob, nid, alen = res
                if best_loc is None or loc > best_loc[0]:
                    best_loc = (loc, glob, nid, alen, t)
                if best_glob is None or glob > best_glob[1]:
                    best_glob = (loc, glob, nid, alen, t)
                # agreement with mmseqs on the SAME target
                for r in mrows[:TOP_N]:
                    if r["target"] == t:
                        agree_delta.append(
                            abs(loc - r["nident"] / r["alnlen"]) * 100)
                        break

            mm_best_e = f"{mrows[0]['evalue']:.3g}" if mrows else "none"
            ph_best_e = f"{prows[0][0]:.3g}" if prows else "none"
            print(f"\n-- {q} ({tiers[q]})")
            print(f"   mmseqs: {len(mrows)} hits, best E={mm_best_e}, "
                  f"{len(mm_sig)} at E<={SIG}")
            print(f"   phmmer: {len(prows)} hits, best E={ph_best_e}, "
                  f"{len(ph_sig)} at E<={SIG}"
                  + (f"  top={prows[0][3][:40]}" if prows else ""))
            if best_loc:
                print(f"   indep SW max LOCAL : {best_loc[0]*100:.1f}% "
                      f"({best_loc[2]}/{best_loc[3]})  {best_loc[4]}")
            if best_glob:
                print(f"   indep SW max GLOBAL: {best_glob[1]*100:.1f}% "
                      f"({best_glob[2]}/{len(qs)})  {best_glob[4]}")
            rows.append({
                "db": db, "query": q, "struct_tier": tiers[q], "qlen": len(qs),
                "mmseqs_n_hits": len(mrows), "mmseqs_best_evalue": mm_best_e,
                "mmseqs_n_sig": len(mm_sig),
                "phmmer_n_hits": len(prows), "phmmer_best_evalue": ph_best_e,
                "phmmer_n_sig": len(ph_sig),
                "phmmer_best_target": prows[0][3] if prows else "",
                "sw_max_local_pct": round(best_loc[0]*100, 1) if best_loc else "",
                "sw_max_local_target": best_loc[4] if best_loc else "",
                "sw_max_global_pct": round(best_glob[1]*100, 1) if best_glob else "",
                "sw_max_global_target": best_glob[4] if best_glob else "",
            })

    out = DATA / "novelty_seq_crosscheck_20261006.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")
    if agree_delta:
        agree_delta.sort()
        n = len(agree_delta)
        print(f"\n[aligner agreement] |SW local id - mmseqs local id| over "
              f"{n} shared query-target pairs: "
              f"median {agree_delta[n//2]:.2f} pp, "
              f"p90 {agree_delta[int(n*0.9)]:.2f} pp, "
              f"max {agree_delta[-1]:.2f} pp")


if __name__ == "__main__":
    main()
