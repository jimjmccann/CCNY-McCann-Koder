#!/usr/bin/env python3
"""Analyse MMseqs2 novelty hits: LOCAL vs GLOBAL identity for the Proteinbase
Level-3/4 sequence gate (30% threshold).

LOCAL  identity = nident / alnlen   (this is mmseqs' native `fident` column)
GLOBAL identity = nident / qlen     (identical residues over the FULL query)

A hit can be 68% local and 30% global; the two answer different questions and
the gate verdict can differ between them.  Report both, always.
"""
import csv
import sys
from pathlib import Path

REPO = Path("/PATH/TO/CHECKOUT")
DATA = REPO / "data" / "novelty_20261006"
QFA = DATA / "novelty_seq_query_20261006.fasta"

# minimum alignment length, as a fraction of query length, for a hit to be
# considered a "comparable" alignment rather than a short spurious fragment.
COMPARABLE_COV = 0.50


def read_queries(path):
    order, tiers, seqs = [], {}, {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            parts = line[1:].split()
            name = parts[0]
            order.append(name)
            seqs[name] = ""
            tiers[name] = next(
                (p.split("=", 1)[1] for p in parts if p.startswith("struct_tier=")), "?"
            )
        elif name:
            seqs[name] += line.strip()
    return order, tiers, seqs


def load_hits(tsv):
    rows = []
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["fident"] = float(r["fident"])
            r["nident"] = int(r["nident"])
            r["alnlen"] = int(r["alnlen"])
            r["qlen"] = int(r["qlen"])
            r["tlen"] = int(r["tlen"])
            r["evalue"] = float(r["evalue"])
            r["bits"] = float(r["bits"])
            r["local_id"] = r["nident"] / r["alnlen"]
            r["global_id"] = r["nident"] / r["qlen"]
            r["qcov"] = r["alnlen"] / r["qlen"]
            rows.append(r)
    return rows


def short_name(theader, target):
    """Readable protein name from a SwissProt or pdb_seqres header."""
    h = theader.strip()
    if h.startswith("sp|") or "|" in h.split()[0]:
        rest = " ".join(h.split()[1:])
        for tag in (" OS=", " OX=", " GN=", " PE=", " SV="):
            if tag in rest:
                rest = rest.split(tag)[0]
        return rest or target
    # pdb_seqres: "1abc_A mol:protein length:123  DESCRIPTION"
    if "mol:" in h:
        tail = h.split("length:", 1)[-1]
        return tail.split(None, 1)[1].strip() if len(tail.split(None, 1)) > 1 else target
    return h


def best(rows, key, reverse=True, filt=None):
    cand = [r for r in rows if (filt is None or filt(r))]
    if not cand:
        return None
    return sorted(cand, key=lambda r: (r[key], -r["evalue"]), reverse=reverse)[0]


def fmt(r):
    if r is None:
        return "n/a"
    return (
        f"{r['local_id']*100:.1f}% local / {r['global_id']*100:.1f}% global "
        f"({r['nident']}/{r['alnlen']} aln, {r['alnlen']}/{r['qlen']} = "
        f"{r['qcov']*100:.0f}% qcov) {r['target']} E={r['evalue']:.2g}"
    )


def verify_from_alignments():
    """⛔ VERIFY, don't trust: recompute nident/alnlen straight from the qaln/taln
    alignment strings in the scratch copy and assert mmseqs' own columns agree.
    This is the check that would have caught the `-a`-missing bug, where mmseqs
    silently reported nident=0 / fident=0 for every hit."""
    import os

    scratch = os.environ.get("DESIGN_DATA_ROOT", "")
    base = Path(scratch) / "novelty_search_20261006"
    total = bad_nident = bad_alnlen = 0
    for name in ("swissprot", "pdbseqres"):
        f = base / f"mmseqs_{name}_hits_full_20261006.tsv"
        if not f.exists():
            print(f"  [verify] MISSING {f} -- cannot verify {name}")
            continue
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                qa, ta = r["qaln"], r["taln"]
                total += 1
                if len(qa) != len(ta):
                    bad_alnlen += 1
                    continue
                rec_nident = sum(1 for a, b in zip(qa, ta) if a == b and a != "-")
                if rec_nident != int(r["nident"]):
                    bad_nident += 1
                if len(qa) != int(r["alnlen"]):
                    bad_alnlen += 1
    print(f"  [verify] {total} hits re-derived from qaln/taln: "
          f"nident mismatches={bad_nident}, alnlen mismatches={bad_alnlen}")
    if total and (bad_nident or bad_alnlen):
        print("  [verify] ⛔ MISMATCH -- do not quote these numbers")
    elif total:
        print("  [verify] OK: mmseqs nident/alnlen reproduce exactly")
    return total, bad_nident, bad_alnlen


def main():
    order, tiers, seqs = read_queries(QFA)
    print("## verification pass")
    verify_from_alignments()
    dbs = {
        "swissprot": DATA / "mmseqs_swissprot_hits_20261006.tsv",
        "pdbseqres": DATA / "mmseqs_pdbseqres_hits_20261006.tsv",
    }
    allhits = {k: load_hits(v) for k, v in dbs.items()}

    out_rows = []
    for dbname, hits in allhits.items():
        byq = {}
        for r in hits:
            byq.setdefault(r["query"], []).append(r)
        for q in order:
            rows = byq.get(q, [])
            rec = {
                "db": dbname,
                "query": q,
                "qlen": len(seqs[q]),
                "struct_tier": tiers[q],
                "n_hits": len(rows),
            }
            # E-value tiers.  SIG=1e-3 is mmseqs' own default significance
            # threshold and is the tier the gate verdict is taken from; the
            # looser tiers are the no-threshold control that proves a "no
            # hits" answer is not just an insensitive search.
            SIG, MARG = 1e-3, 10.0
            sig = lambda r: r["evalue"] <= SIG
            marg = lambda r: r["evalue"] <= MARG
            cov = lambda r: r["qcov"] >= COMPARABLE_COV

            views = {
                # the hit mmseqs itself calls most significant
                "bestE": best(rows, "bits"),
                # ---- PRIMARY tier: E <= 1e-3 (mmseqs default significance) --
                "sig_maxglobal": best(rows, "global_id", filt=sig),
                "sig_maxlocal_cov": best(
                    rows, "local_id", filt=lambda r: sig(r) and cov(r)),
                "sig_maxlocal_raw": best(rows, "local_id", filt=sig),
                # ---- marginal tier: E <= 10 --------------------------------
                "marg_maxglobal": best(rows, "global_id", filt=marg),
                "marg_maxlocal_cov": best(
                    rows, "local_id", filt=lambda r: marg(r) and cov(r)),
                # ---- no-threshold control (E up to 10000) ------------------
                "ctl_maxglobal": best(rows, "global_id"),
                "ctl_maxlocal_cov": best(rows, "local_id", filt=cov),
                "ctl_maxlocal_raw": best(rows, "local_id"),
            }
            for tag, r in views.items():
                if r is None:
                    rec.update({f"{tag}_{k}": "" for k in
                                ("target", "name", "local", "global", "alnlen",
                                 "qcov", "evalue", "bits")})
                    continue
                rec[f"{tag}_target"] = r["target"]
                rec[f"{tag}_name"] = short_name(r["theader"], r["target"])
                rec[f"{tag}_local"] = round(r["local_id"] * 100, 1)
                rec[f"{tag}_global"] = round(r["global_id"] * 100, 1)
                rec[f"{tag}_alnlen"] = r["alnlen"]
                rec[f"{tag}_qcov"] = round(r["qcov"] * 100, 1)
                rec[f"{tag}_evalue"] = f"{r['evalue']:.3g}"
                rec[f"{tag}_bits"] = r["bits"]
            out_rows.append(rec)

    summ = DATA / "novelty_seq_summary_20261006.csv"
    with open(summ, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {summ}  ({len(out_rows)} rows)")

    # ---- top-5 per query per db, both conventions ------------------------
    top5 = DATA / "novelty_seq_top5_20261006.tsv"
    with open(top5, "w") as fh:
        fh.write("db\tquery\trank\tlocal_id_pct\tglobal_id_pct\tnident\talnlen\t"
                 "qlen\tqcov_pct\tevalue\tbits\ttarget\tprotein_name\n")
        for dbname, hits in allhits.items():
            byq = {}
            for r in hits:
                byq.setdefault(r["query"], []).append(r)
            for q in order:
                rows = sorted(byq.get(q, []), key=lambda r: -r["bits"])[:5]
                for i, r in enumerate(rows, 1):
                    fh.write(
                        f"{dbname}\t{q}\t{i}\t{r['local_id']*100:.1f}\t"
                        f"{r['global_id']*100:.1f}\t{r['nident']}\t{r['alnlen']}\t"
                        f"{r['qlen']}\t{r['qcov']*100:.1f}\t{r['evalue']:.3g}\t"
                        f"{r['bits']:.0f}\t{r['target']}\t"
                        f"{short_name(r['theader'], r['target'])}\n"
                    )
    print(f"wrote {top5}")

    # ---- console verdict -------------------------------------------------
    for dbname in dbs:
        print(f"\n################ {dbname} ################")
        for rec in [r for r in out_rows if r["db"] == dbname]:
            q, tier = rec["query"], rec["struct_tier"]
            print(f"\n-- {q}  (len {rec['qlen']}, struct_tier={tier}, "
                  f"{rec['n_hits']} permissive hits)")
            for tag, label in (
                ("bestE", "top bitscore hit            "),
                ("sig_maxglobal", "E<=1e-3  max GLOBAL        "),
                ("sig_maxlocal_cov", "E<=1e-3  max LOCAL qcov>=50%"),
                ("sig_maxlocal_raw", "E<=1e-3  max LOCAL any len "),
                ("marg_maxglobal", "E<=10    max GLOBAL        "),
                ("marg_maxlocal_cov", "E<=10    max LOCAL qcov>=50%"),
                ("ctl_maxglobal", "CONTROL  max GLOBAL        "),
                ("ctl_maxlocal_cov", "CONTROL  max LOCAL qcov>=50%"),
                ("ctl_maxlocal_raw", "CONTROL  max LOCAL any len "),
            ):
                if not rec.get(f"{tag}_target"):
                    print(f"   {label}: none")
                    continue
                print(f"   {label}: {rec[f'{tag}_local']}% local / "
                      f"{rec[f'{tag}_global']}% global  "
                      f"aln={rec[f'{tag}_alnlen']} qcov={rec[f'{tag}_qcov']}%  "
                      f"E={rec[f'{tag}_evalue']}  {rec[f'{tag}_target']}  "
                      f"{rec[f'{tag}_name'][:55]}")


if __name__ == "__main__":
    main()
