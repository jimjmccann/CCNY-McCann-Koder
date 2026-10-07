#!/usr/bin/env python3
"""Write the public deposit's FASTA by READING the submitted CSV. No sequence is typed here.

WHY THIS FILE EXISTS (2026-10-06). An earlier staged FASTA, written at 10:35 EDT and not
deposited, held FOUR sequences including `id15`. The submission being sent
was built at 12:12 EDT and holds TEN, with `id15` DROPPED and `id86` ADDED
(by the submission CSV builder, which is internal to the campaign and is not
deposited). A deposit whose
sequence file disagrees with the submission is the one error in a deposit that cannot be
explained away, so the deposited FASTA is GENERATED from the submitted CSV and never edited.

Local, no compute, no network. Reads one CSV, writes one FASTA. Refuses an existing --out
unless --force is given (never silently overwrite).

Usage:
  python3 make_deposit_fasta.py --csv <submission csv> --out <fasta> [--force]
"""
import argparse
import csv
import os
import sys

STD = set("ACDEFGHIKLMNPQRSTVWY")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="the submitted CSV (name,sequence,molecule_class)")
    ap.add_argument("--out", required=True, help="FASTA to write")
    ap.add_argument("--force", action="store_true", help="allow overwriting --out")
    a = ap.parse_args()

    if os.path.exists(a.out) and not a.force:
        sys.exit(f"⛔ {a.out} exists. Pass --force to overwrite.")

    with open(a.csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"⛔ {a.csv} has no rows")
    for col in ("name", "sequence", "molecule_class"):
        if col not in rows[0]:
            sys.exit(f"⛔ {a.csv} has no '{col}' column -- is this the submission CSV?")

    # ⛔ Assert the organizer's own bounds again here. The CSV builder already does; a deposit
    #    that re-asserts cannot quietly ship a file the submission would have been rejected for.
    for r in rows:
        seq = r["sequence"].strip().upper()
        if not (10 <= len(seq) <= 250):
            sys.exit(f"⛔ {r['name']}: {len(seq)} aa outside the organizer's 10-250 bound")
        bad = sorted(set(seq) - STD)
        if bad:
            sys.exit(f"⛔ {r['name']}: non-standard residue(s) {bad}")
        r["sequence"] = seq
    if len(rows) > 20:
        sys.exit(f"⛔ {len(rows)} sequences exceeds the organizer's 20-per-submission cap")

    lines = [
        "; Challenge-1 submitted sequences, James McCann and the Koder Group, The City College of New York.",
        "; GENERATED from the submitted CSV by scripts/make_deposit_fasta.py -- do not hand-edit.",
        f"; source CSV: {os.path.basename(a.csv)}   sequences: {len(rows)}",
        "; Order is the ranking we believe, best first. Why these molecules: DESIGN_CHOICES.md.",
        "; ⛔ No value anywhere in this deposit is an affinity, and none of these has been tested.",
        "",
    ]
    for i, r in enumerate(rows, 1):
        lines.append(
            f">{r['name']} rank={i} len={len(r['sequence'])} class={r['molecule_class']}"
        )
        s = r["sequence"]
        for j in range(0, len(s), 60):
            lines.append(s[j:j + 60])

    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"wrote {a.out}  ({len(rows)} sequences)")
    for i, r in enumerate(rows, 1):
        print(f"   {i:2d}  {r['name']:22s} {len(r['sequence']):3d} aa")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
