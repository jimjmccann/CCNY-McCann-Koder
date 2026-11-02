#!/usr/bin/env python3
"""j953 -- measure the identity of the two CONTAINED constellation positions in every deposited
sequence, for both sequence cohorts.

WHY THIS EXISTS, AND WHY THE QC TABLES DO NOT ANSWER IT
-------------------------------------------------------
Four of the six constellation groups enter the design as single-residue guideposts: their
identity is fixed, and `pin_identity_held` in the sequence-QC tables checks exactly those four
(`n_pins = 4`, `pin_positions = 17,95,112,133`). The other two groups -- scaffold positions 34
and 38 -- are not guideposts at all. They sit INSIDE held element 3 (helix 33-40), whose
sequence the design run released so that the surface could be redesigned.

=> the QC column cannot see them, in EITHER cohort. A cohort can therefore report
`pin_identity_held = 1` for every sequence while having lost two of the six functional groups.
That is what happened to the first cohort, and this script is the measurement that shows it.

HOW THE TWO POSITIONS ARE LOCATED -- not hard-coded
---------------------------------------------------
Each backbone samples its own connector lengths, so scaffold position 34 lands at a different
index in every design. Element 3 is held, so its sequence `FHQFSKFY` is present verbatim and
identical in every design built on that backbone. The rows with `is_input = 1` carry the
backbone's own sequence; the motif is located in those, and the two positions are read at motif
offset 1 (scaffold 34) and offset 5 (scaffold 38). The motif is asserted to occur exactly once
in every input row, and the per-design offset is asserted to be reused for that design's
designed rows only.

usage: python3 -I j953_contained_pin_check_20261102.py <sequences.tsv> [<sequences.tsv> ...] --out <tsv>
"""
import argparse
import collections
import csv
import os
import sys

MOTIF = "FHQFSKFY"          # held element 3, scaffold residues 33-40
OFF_34 = 1                  # scaffold position 34 within the motif
OFF_38 = 5                  # scaffold position 38 within the motif
INTENT_34 = "H"             # the design specified an imidazole here
CATIONS = ("K", "R")        # the design specified a cation at 38; either satisfies the type


def read_rows(path):
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def offsets(rows, path):
    """motif offset per design, taken from that design's own is_input row."""
    off = {}
    for r in rows:
        if r["is_input"] != "1":
            continue
        s = r["sequence"]
        n = s.count(MOTIF)
        if n != 1:
            sys.exit("REFUSING: %s: input row %s carries the held motif %d times, expected 1"
                     % (os.path.basename(path), r["design"], n))
        off[r["design"]] = s.find(MOTIF)
    if not off:
        sys.exit("REFUSING: %s: no is_input=1 rows, so the held motif cannot be located"
                 % os.path.basename(path))
    return off


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tables", nargs="+")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit("REFUSING: %s already exists (use a fresh --out per run)" % a.out)

    out = []
    for path in a.tables:
        rows = read_rows(path)
        off = offsets(rows, path)
        des = [r for r in rows if r["is_input"] != "1"]
        c34 = collections.Counter()
        c38 = collections.Counter()
        for r in des:
            d = r["design"]
            if d not in off:
                sys.exit("REFUSING: %s: designed row %s/%s has no is_input row to locate the "
                         "motif" % (os.path.basename(path), d, r["seq_id"]))
            i = off[d]
            s = r["sequence"]
            aa34, aa38 = s[i + OFF_34], s[i + OFF_38]
            c34[aa34] += 1
            c38[aa38] += 1
            out.append({
                "cohort": os.path.basename(path),
                "design": d,
                "arm": r["arm"],
                "seq_id": r["seq_id"],
                "motif_offset": i,
                "aa_at_scaffold_34": aa34,
                "intent_34_imidazole_held": int(aa34 == INTENT_34),
                "aa_at_scaffold_38": aa38,
                "intent_38_cation_held": int(aa38 in CATIONS),
                "n_his_reported": r.get("n_his", ""),
            })
        n = len(des)
        print("%s: %d designed sequences" % (os.path.basename(path), n))
        print("   scaffold 34 (intent %s): held %d/%d (%.1f%%)   composition %s"
              % (INTENT_34, c34[INTENT_34], n, 100.0 * c34[INTENT_34] / n,
                 dict(c34.most_common())))
        nc = sum(c38[x] for x in CATIONS)
        print("   scaffold 38 (intent cation): held %d/%d (%.1f%%)   composition %s"
              % (nc, n, 100.0 * nc / n, dict(c38.most_common())))

    cols = ["cohort", "design", "arm", "seq_id", "motif_offset", "aa_at_scaffold_34",
            "intent_34_imidazole_held", "aa_at_scaffold_38", "intent_38_cation_held",
            "n_his_reported"]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        fh.write("# Identity at the two CONTAINED constellation positions (scaffold 34 and 38),\n")
        fh.write("# measured per deposited sequence. These two are NOT guideposts and are NOT\n")
        fh.write("# covered by pin_identity_held in the sequence-QC tables -- see the module\n")
        fh.write("# docstring of j953_contained_pin_check_20261102.py.\n")
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print("wrote %s (%d rows)" % (a.out, len(out)))


if __name__ == "__main__":
    main()
