#!/usr/bin/env python3
"""Consolidate the MMseqs2 + phmmer sequence-novelty measurements into the
markdown tables used in SEQUENCE_NOVELTY_SWISSPROT_20261006.md.

Emits tables to stdout so every number in the report is machine-generated
rather than transcribed by hand.
"""
import csv
import os
from collections import defaultdict
from pathlib import Path

REPO = Path("/PATH/TO/CHECKOUT")
DATA = REPO / "data" / "novelty_20261006"
GATE = 30.0


def read_queries():
    order, tiers, qlen = [], {}, {}
    name = None
    for line in (DATA / "novelty_seq_query_20261006.fasta").read_text().splitlines():
        if line.startswith(">"):
            parts = line[1:].split()
            name = parts[0]
            order.append(name)
            qlen[name] = 0
            tiers[name] = next((p.split("=", 1)[1] for p in parts
                                if p.startswith("struct_tier=")), "?")
        elif name:
            qlen[name] += len(line.strip())
    return order, tiers, qlen


def load_mmseqs():
    """-> {(db,query): [rows]} with local/global identity attached."""
    out = defaultdict(list)
    for db in ("swissprot", "pdbseqres"):
        with open(DATA / f"mmseqs_{db}_hits_20261006.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                nid, aln, ql = int(r["nident"]), int(r["alnlen"]), int(r["qlen"])
                r.update(evalue=float(r["evalue"]), local=nid / aln * 100,
                         glob=nid / ql * 100, qcov=aln / ql, nid=nid,
                         aln=aln, db=db)
                out[(db, r["query"])].append(r)
    return out


def load_ph():
    out = defaultdict(dict)
    p = DATA / "novelty_seq_phmmer_identity_20261006.csv"
    for r in csv.DictReader(open(p)):
        out[(r["db"], r["query"])][r["view"]] = r
    return out


def mx(rows, key, filt=None):
    c = [r for r in rows if (filt is None or filt(r))]
    return max(c, key=lambda r: r[key]) if c else None


def main():
    order, tiers, qlen = read_queries()
    mm, ph = load_mmseqs(), load_ph()
    SIG = 1e-3

    # ================= TABLE 1: significant-hit verdict ==================
    print("### Table 1 — significant hits (E <= 1e-3) and the gate verdict\n")
    print("| design | len | struct_tier | MMseqs2 sig hits (SP / PDB) | "
          "phmmer domain-sig hits | max LOCAL id (sig) | max GLOBAL id (sig) | "
          "sequence axis |")
    print("|---|---|---|---|---|---|---|---|")
    verdict = {}
    for q in order:
        nsp = len([r for r in mm[("swissprot", q)] if r["evalue"] <= SIG])
        npd = len([r for r in mm[("pdbseqres", q)] if r["evalue"] <= SIG])
        # significant per mmseqs, both dbs pooled
        sig = [r for d in ("swissprot", "pdbseqres")
               for r in mm[(d, q)] if r["evalue"] <= SIG]
        # phmmer domain-level significant
        phd = [ph[(d, q)].get("sigdom_maxglobal")
               for d in ("swissprot", "pdbseqres")]
        phd = [r for r in phd if r]
        bl = mx(sig, "local", lambda r: r["qcov"] >= 0.5)
        bg = mx(sig, "glob")
        phl = max((float(r["local_pct"]) for r in phd), default=None)
        phg = max((float(r["global_pct"]) for r in phd), default=None)
        loc = max([v for v in (bl["local"] if bl else None, phl) if v is not None],
                  default=None)
        glb = max([v for v in (bg["glob"] if bg else None, phg) if v is not None],
                  default=None)
        has_sig = bool(sig) or bool(phd)
        if not has_sig:
            v = "**PASS** (no significant hit)"
        elif glb is not None and glb > GATE:
            v = f"**FAIL** ({glb:.1f}% > 30%)"
        else:
            v = f"**PASS** ({glb:.1f}% <= 30%)"
        verdict[q] = (has_sig, loc, glb, v)
        print(f"| `{q}` | {qlen[q]} | {tiers[q]} | {nsp} / {npd} | "
              f"{len(phd)} | {f'{loc:.1f}%' if loc is not None else 'n/a'} | "
              f"{f'{glb:.1f}%' if glb is not None else 'n/a'} | {v} |")

    # ============ TABLE 2: no-threshold control (worst case) =============
    print("\n### Table 2 — no-threshold CONTROL: the single worst hit at any "
          "E-value (up to E=10000)\n")
    print("| design | max LOCAL id (qcov>=50%) | its hit / E | "
          "max GLOBAL id | its hit / E | above 30% global? |")
    print("|---|---|---|---|---|---|")
    for q in order:
        allr = [r for d in ("swissprot", "pdbseqres") for r in mm[(d, q)]]
        phall = [ph[(d, q)].get("all_maxglobal")
                 for d in ("swissprot", "pdbseqres")]
        phall = [r for r in phall if r]
        bl = mx(allr, "local", lambda r: r["qcov"] >= 0.5)
        bg = mx(allr, "glob")
        # fold in phmmer's own worst case
        php = max(phall, key=lambda r: float(r["global_pct"])) if phall else None
        gl, gt, ge = (bg["glob"], bg["target"], f"{bg['evalue']:.2g}") if bg else (None, "", "")
        if php and float(php["global_pct"]) > (gl or 0):
            gl, gt, ge = (float(php["global_pct"]), php["target"],
                          f"{php['phmmer_full_evalue']} (ph)")
        flag = "yes" if gl is not None and gl > GATE else "no"
        print(f"| `{q}` | {bl['local']:.1f}% | {bl['target']} / "
              f"{bl['evalue']:.2g} | {gl:.1f}% | {gt} / {ge} | {flag} |")

    # ============ TABLE 3: top-5 significant hits per design =============
    print("\n### Table 3 — top hits per design (MMseqs2, by bitscore, both DBs)\n")
    print("| design | db | hit | protein | local id | global id | nident/alnlen | "
          "qcov | E-value |")
    print("|---|---|---|---|---|---|---|---|---|")
    with open(DATA / "novelty_seq_top5_20261006.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["rank"]) > 3:
                continue
            print(f"| `{r['query']}` | {r['db']} | {r['target']} | "
                  f"{r['protein_name'][:42]} | {r['local_id_pct']}% | "
                  f"{r['global_id_pct']}% | {r['nident']}/{r['alnlen']} | "
                  f"{r['qcov_pct']}% | {r['evalue']} |")

    # ================= summary for the PASS cohort =======================
    print("\n### PASS_L3 cohort summary\n")
    for q in order:
        if not tiers[q].startswith("PASS"):
            continue
        has_sig, loc, glb, v = verdict[q]
        allr = [r for d in ("swissprot", "pdbseqres") for r in mm[(d, q)]]
        phall = [ph[(d, q)].get("all_maxglobal") for d in ("swissprot", "pdbseqres")]
        phall = [r for r in phall if r]
        ctl_g = max([mx(allr, "glob")["glob"]]
                    + [float(r["global_pct"]) for r in phall])
        print(f"- `{q}` ({tiers[q]}): significant hits = "
              f"{'YES' if has_sig else 'NONE'}; "
              f"sig GLOBAL = {f'{glb:.1f}%' if glb is not None else 'n/a'}; "
              f"control worst GLOBAL = {ctl_g:.1f}%  -> {v}")


if __name__ == "__main__":
    main()
