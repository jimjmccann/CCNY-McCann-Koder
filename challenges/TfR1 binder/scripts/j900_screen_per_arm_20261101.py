#!/usr/bin/env python3
"""j900 -- aggregate the ten-arm 5XOX core screen PER CORE OPTION, from the per-design gate TSV.

Local, reads one TSV. ⛔ Runs no model and applies NO verdict column.

WHAT THIS IS FOR
----------------
400 designs, ten core options (opt01-opt10), 40 designs each. The question the screen was
launched to answer is WHICH CORE OPTION to take into depth. That is a per-arm question, and
The gate deliberately writes one row per DESIGN with no arm column, so the grouping happens here.

THREE RULES THIS OBEYS, EACH FOR A MEASURED REASON
--------------------------------------------------
1. ⛔⛔ n IS BATCHES, NOT DESIGNS. Each arm ran
   nbatch 5 x 8 models = 40 files, but one contig length is drawn PER BATCH, so the 8 models
   inside a batch share a draw and are NOT independent. Effective n per arm is 5. Every count
   below is therefore reported WITH its batch spread, and the per-batch table is printed so a
   40/40 that is really 5/5 cannot read as strong evidence.
2. ⛔ NO BINARY CONTINUITY CUT. The gate's own bar ladder on this cohort found the <=5/<=7/<=10 A
   counts DIFFER (351/353/371), so designs sit between "ideal" and "broken" and no cut is
   admissible, and continuity must be gated before any clash metric because a clash floor
   rewards a broken chain. max_step is reported as a
   DISTRIBUTION (median and the count under each bar), never as pass/fail.
3. ⛔ NEVER A MEAN OF FRACTIONS -- an averaged label pays partial credit. Arms are
   described by medians of continuous columns plus explicit counts; the conjunction over
   criteria is left to the reader and printed as its own column, not averaged in.

⛔ CONTINUITY AND FOLD TOPOLOGY ARE ORTHOGONAL: max_step and min_cut measure different
   failures and must not be collapsed. min_cut is the "is the chain ONE fold" number -- a low
   value means no crossing long-range contacts, i.e. two lobes on a tether; rg_ratio is the only
   length-safe compactness column, which matters here because the arms differ in length.

⚠️ CLASH IS REPORTED BUT NOT GATED. Raw generator output clashes normally; the documented raw
   bar is bb_min >= 2.0 and a tight
   gate belongs AFTER a relax. The clash columns are printed so the arms can be compared, and
   explicitly NOT used to rank.

usage: python3 -I j900_screen_per_arm_20261101.py ../data/j900_screen_gate_20261101.tsv [--out <tsv>]
"""
import argparse
import os
import re
import statistics as st
import sys

# design name: arm_shard00_tfr1_core_opt01_<batch>_model_<model>
NAME = re.compile(r"arm_shard(\d+)_tfr1_core_(opt\d+)_(\d+)_model_(\d+)$")

# bars printed for the continuity distribution; ⛔ re-measured per cohort, never inherited
BARS = (4.0, 4.1, 4.5, 5.0, 7.0, 10.0)


def med(v):
    """Median, rounded to 5 decimals.

    The median of an even-sized sample is the mean of the two middle values, which introduces
    binary floating-point noise in the last digits (3.9031 arriving as 3.9031000000000002). The
    inputs carry 4-5 significant decimals, so reporting more digits is noise, not precision.
    """
    if not v:
        return float("nan")
    return round(st.median(v), 5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    rows, hdr = [], None
    for ln in open(a.tsv):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            continue
        rows.append(dict(zip(hdr, f)))
    print("== j900 per-arm aggregation ==")
    print("   source : %s" % a.tsv)
    print("   rows   : %d" % len(rows))

    arms, unparsed = {}, []
    for r in rows:
        m = NAME.match(r["design"])
        if not m:
            unparsed.append(r["design"])
            continue
        _shard, opt, batch, _model = m.groups()
        arms.setdefault(opt, []).append((int(batch), r))
    # ⛔ an unparsed name is a silent row loss, so it is REFUSED, not skipped
    if unparsed:
        print("⛔⛔ %d design names did not parse, e.g. %s" % (len(unparsed), unparsed[:3]))
        sys.exit(1)
    print("   arms   : %d  (%s)" % (len(arms), ",".join(sorted(arms))))
    tot = sum(len(v) for v in arms.values())
    if tot != len(rows):
        print("⛔⛔ grouped %d of %d rows" % (tot, len(rows)))
        sys.exit(1)

    def col(items, name, cast=float):
        out = []
        for _b, r in items:
            try:
                out.append(cast(r[name]))
            except (KeyError, ValueError):
                pass
        return out

    out_rows = []
    print("\n" + "=" * 112)
    print("PER CORE OPTION. n=40 designs but n_eff=5 BATCHES -- the batch column is the honest n.")
    print("=" * 112)
    hdr_line = ("%-7s %4s %4s | %8s %8s | %7s %7s | %7s %7s | %7s %6s | %6s" %
                ("opt", "n", "nb", "medStep", "<=4.1", "medCut", "medLR", "medRgR",
                 "medBBmin", "bb>=2", "medPin", "nres"))
    print(hdr_line)
    print("-" * 112)
    for opt in sorted(arms):
        items = arms[opt]
        nb = len({b for b, _ in items})
        ms = col(items, "max_step")
        mc = col(items, "min_cut")
        lr = col(items, "lr_frac")
        rg = col(items, "rg_ratio")
        bb = col(items, "bb_min_to_target")
        pin = col(items, "pin_cb_min_to_target")
        nres = col(items, "n_res")
        n41 = sum(1 for v in ms if v <= 4.1)
        nbb2 = sum(1 for v in bb if v >= 2.0)
        print("%-7s %4d %4d | %8.3f %4d/%-3d | %7.1f %7.3f | %7.3f %7.3f | %4d/%-3d %6.3f | %6.0f"
              % (opt, len(items), nb, med(ms), n41, len(ms), med(mc), med(lr),
                 med(rg), med(bb), nbb2, len(bb), med(pin), med(nres)))
        row = dict(opt=opt, n=len(items), n_batches=nb,
                   med_max_step=med(ms), n_step_le_4p1=n41,
                   med_min_cut=med(mc), med_lr_frac=med(lr), med_rg_ratio=med(rg),
                   med_bb_min=med(bb), n_bb_ge_2p0=nbb2,
                   med_pin_cb_min=med(pin), med_n_res=med(nres))
        for bar in BARS:
            row["n_step_le_%s" % bar] = sum(1 for v in ms if v <= bar)
        out_rows.append(row)

    print("\n--- CONTINUITY DISTRIBUTION per arm (count with NO CA-CA step above each bar) ---")
    print("%-7s %s" % ("opt", " ".join("%7s" % ("<=%.1f" % b) for b in BARS)))
    for r in out_rows:
        print("%-7s %s" % (r["opt"],
              " ".join("%7d" % r["n_step_le_%s" % b] for b in BARS)))

    print("\n--- PER-BATCH continuity, the honest n (count of 8 models with max_step <= 4.1) ---")
    print("⛔ a 40/40 that is really 5 batches of 8/8 is FIVE observations, not forty.")
    print("%-7s %s" % ("opt", " ".join("%4s" % ("b%d" % i) for i in range(5))))
    for opt in sorted(arms):
        per = {}
        for b, r in arms[opt]:
            try:
                per.setdefault(b, []).append(float(r["max_step"]))
            except ValueError:
                pass
        cells = []
        for b in sorted(per):
            cells.append("%d/%d" % (sum(1 for v in per[b] if v <= 4.1), len(per[b])))
        print("%-7s %s" % (opt, " ".join("%4s" % c for c in cells)))

    if a.out:
        keys = ["opt", "n", "n_batches", "med_max_step", "n_step_le_4p1",
                "med_min_cut", "med_lr_frac", "med_rg_ratio", "med_bb_min",
                "n_bb_ge_2p0", "med_pin_cb_min", "med_n_res"] + \
               ["n_step_le_%s" % b for b in BARS]
        with open(a.out, "w") as fh:
            fh.write("# j900 per-arm aggregation of %s\n" % os.path.basename(a.tsv))
            fh.write("# n_eff is n_batches (8 models share one contig draw), NOT n.\n")
            fh.write("# ⛔ no verdict column, no cut applied; clash reported, NOT ranked on.\n")
            fh.write("\t".join(keys) + "\n")
            for r in out_rows:
                fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
        print("\n⇒ wrote %s" % a.out)
    print("\n⛔ No arm is declared a winner here. Ranking is the reader's job and the batch")
    print("   column is the n that any ranking must respect.")


main()
