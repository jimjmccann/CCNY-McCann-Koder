#!/usr/bin/env python3
"""j903 -- collect the SolubleMPNN output for the three leading arms into ONE table.

Local. ⛔ Designs nothing; it reads what the sequence-design driver already wrote.

⛔⛔ NEVER `cat *.fa | grep -c '^>'`. MEASURED 2026-11-01: that idiom reported 129 records where
the per-file sum is 144, undercounting by exactly (n_files - 1). LigandMPNN's .fa files have NO
trailing newline, so `cat` fuses each file's last SEQUENCE line onto the next file's '>' header
and the header stops being at line start. ⇒ every file is parsed on its own, and the per-arm
record count is asserted against n_backbones * (1 + n_seqs). A count is not a measurement.

⛔ THE FASTA HOLDS ALL CHAINS, ':'-separated. The DESIGNED chain is extracted by position and its
length is asserted against the designed chain's own length in the backbone; the target chains
are carried but never scored, and the chain is always named rather than assumed.

⛔ The first record of each file (no `id=`) is the INPUT sequence, not a design. It is kept,
labelled `is_input`, and excluded from every design statistic -- counting it as a design would
inflate each backbone's yield by 1/9.

usage: python3 -I j903_collect_mpnn_20261101.py --runs <dir> --arms opt05,opt01,opt04 --out <tsv>
"""
import argparse
import collections
import glob
import os
import re
import sys

HDR = re.compile(r"id=(\d+)")
CONF = re.compile(r"overall_confidence=([0-9.]+)")
REC = re.compile(r"seq_rec=([0-9.]+)")
NSEQ = 8


def read_fasta(path):
    """-> [(header, sequence)] parsed from ONE file. ⛔ never from a cat of many."""
    out, h, buf = [], None, []
    for ln in open(path):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if h is not None:
                out.append((h, "".join(buf)))
            h, buf = ln[1:], []
        elif ln:
            buf.append(ln)
    if h is not None:
        out.append((h, "".join(buf)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--arms", required=True)
    ap.add_argument("--dir-template", default="run_%s",
                    help="subdir name per arm; the bias rounds use '%%s' directly")
    ap.add_argument("--out", required=True)
    ap.add_argument("--qc-out", default="",
                    help="ALSO write a sequence-QC input TSV: header on line 1 (⛔ csv.DictReader "
                         "does not skip '#' comments), DESIGNED rows only (an is_input row has no "
                         "confidence and the QC reader refuses it), and pin columns from report.json")
    a = ap.parse_args()
    arms = a.arms.split(",")

    rows, fails = [], []
    print("== j903 collect ==")
    for arm in arms:
        seqd = os.path.join(a.runs, a.dir_template % arm, "mpnn", "seqs")
        files = sorted(glob.glob(os.path.join(seqd, "*")))
        if not files:
            fails.append("%s: no fasta in %s" % (arm, seqd))
            continue
        rep = {}
        rp = os.path.join(a.runs, a.dir_template % arm, "report.json")
        if os.path.exists(rp):
            import json as _j
            rep = _j.load(open(rp))
        n_rec = n_des = 0
        for f in files:
            recs = read_fasta(f)
            n_rec += len(recs)
            design = os.path.basename(f).rsplit(".", 1)[0]
            _r = rep.get(design + ".pdb", {})
            _pins = _r.get("pins", [])
            _pid = _r.get("pin_identities", [])
            for hdr, seq in recs:
                chains = seq.split(":")
                binder = chains[0]
                m = HDR.search(hdr)
                is_input = m is None
                if not is_input:
                    n_des += 1
                rows.append(dict(
                    arm=arm, design=design,
                    seq_id=("input" if is_input else m.group(1)),
                    is_input=int(is_input),
                    n_chains=len(chains), len_binder=len(binder),
                    overall_confidence=(CONF.search(hdr).group(1)
                                        if CONF.search(hdr) else ""),
                    seq_rec=(REC.search(hdr).group(1) if REC.search(hdr) else ""),
                    n_cys=binder.count("C"), n_his=binder.count("H"),
                    pin_positions=",".join(str(x) for x in _pins),
                    pin_identities="".join(_pid),
                    sequence=binder))
        want = len(files) * (1 + NSEQ)
        ok = n_rec == want and n_des == len(files) * NSEQ
        print("  %s  %-6s %2d backbones, %3d records (expect %d), %3d designed (expect %d)"
              % ("PASS" if ok else "FAIL", arm, len(files), n_rec, want,
                 n_des, len(files) * NSEQ))
        if not ok:
            fails.append("%s: record count mismatch" % arm)

    if fails:
        for f in fails:
            print("  ⛔ %s" % f)
        sys.exit(1)

    des = [r for r in rows if not r["is_input"]]
    ncys = sum(r["n_cys"] for r in des)
    print("\n  designed sequences : %d" % len(des))
    print("  ⛔ cysteines in designed binders: %d  (MUST be 0)" % ncys)
    if ncys:
        sys.exit("⛔⛔ a designed binder carries CYS -- omit_AA C did not hold.")
    lens = collections.Counter(r["len_binder"] for r in des)
    print("  binder length range: %d..%d aa" % (min(lens), max(lens)))
    print("  distinct sequences : %d of %d" % (len({r["sequence"] for r in des}), len(des)))
    for arm in arms:
        d = [r for r in des if r["arm"] == arm]
        if not d:
            continue
        conf = [float(r["overall_confidence"]) for r in d if r["overall_confidence"]]
        rec = [float(r["seq_rec"]) for r in d if r["seq_rec"]]
        conf.sort()
        rec.sort()
        print("  %-6s n=%3d  confidence med %.4f (%.4f-%.4f)  seq_rec med %.4f"
              % (arm, len(d), conf[len(conf) // 2], conf[0], conf[-1],
                 rec[len(rec) // 2]))

    keys = ["arm", "design", "seq_id", "is_input", "n_chains", "len_binder",
            "overall_confidence", "seq_rec", "n_cys", "n_his", "sequence"]
    with open(a.out, "w") as fh:
        fh.write("# j903 SolubleMPNN sequences, arms %s\n" % a.arms)
        fh.write("# is_input=1 is the backbone's OWN sequence, not a design -- exclude it.\n")
        fh.write("# sequence = the DESIGNED chain only; target chains dropped.\n")
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[k]) for k in keys) + "\n")
    print("\n  wrote %s (%d rows)" % (a.out, len(rows)))

    if a.qc_out:
        if os.path.exists(a.qc_out):
            sys.exit("⛔ %s already exists (use a fresh --qc-out per cohort)." % a.qc_out)
        qkeys = ["design", "arm", "seq_id", "length", "n_pins", "pin_positions",
                 "pin_identities", "n_cys", "n_his", "sequence"]
        with open(a.qc_out, "w") as fh:
            fh.write("\t".join(qkeys) + "\n")   # ⛔ header FIRST, no comment lines
            for r in des:                        # ⛔ designed rows only
                fh.write("\t".join(str({
                    "design": r["design"], "arm": r["arm"], "seq_id": r["seq_id"],
                    "length": r["len_binder"], "n_pins": len(r["pin_identities"]),
                    "pin_positions": r["pin_positions"],
                    "pin_identities": r["pin_identities"],
                    "n_cys": r["n_cys"], "n_his": r["n_his"],
                    "sequence": r["sequence"]}[k]) for k in qkeys) + "\n")
        print("  wrote %s (%d designed rows, QC-input format)" % (a.qc_out, len(des)))


main()
