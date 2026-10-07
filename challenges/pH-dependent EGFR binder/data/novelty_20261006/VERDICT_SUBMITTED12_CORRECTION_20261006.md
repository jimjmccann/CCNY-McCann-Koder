# Correction notice for `verdict_submitted12.json`

**`verdict_submitted12.json` is kept exactly as the program produced it and has not been edited.**
This note sits beside it because one of its fields has since been measured to be wrong, and a raw
output file is the wrong place to hide that.

## What the file is

Per-design verdicts from the **structural half** of the 2026-10-06 novelty measurement, computed on
**chain A extracted from a predicted binder-receptor complex**. Twelve entries. The field
`best_tm_cov70` is the best TM-score among foldseek hits covering more than 70% of the query, and
`high_struct` is simply `best_tm_cov70 >= 0.80`.

⚠️ **The `12` in the filename is the number of designs scored in that run, not a submitted-set
count.** Two of the twelve, `egfr-bridge-id1104` and `egfr-bridge-id47`, are **not** in the
submission; `SEQUENCE_SET_CORRECTION_20261006.md` is authoritative on which sequences are.

## The four `"high_struct": true` entries, and which one is wrong

| design | `best_tm_cov70` in this file | status |
|---|---|---|
| `egfr-siteA-id2015` | 0.8847 | ⛔ **DOES NOT REPRODUCE** — see below |
| `egfr-bridge-id1104` | 0.9751 | stands |
| `egfr-bridge-id47` | 0.9668 | stands |
| `egfr-bridge-id85` | 0.9773 | stands |

**`egfr-siteA-id2015`:** re-measured on the **monomer fold** — the object the organizers' pipeline
actually evaluates — the same design gives **0.590 median, 0.673 max, 0 of 63 draws at or above
0.80**. The 0.885 in this file is an artifact of measuring a chain pulled out of a complex. Full
measurement and its limits: `../../NOVELTY_CORRECTION_MONOMER_20261006.md`.

**The three bridge entries stand.** They were not re-measured on the monomer, and their high
structural similarity to the published consensus TPR proteins is the subject of
`../../accidental_ctpr_homology.txt`, which is independent of how the chain was extracted.

## Also note

Every `source` value in this file is a path on our own machines with the storage root replaced by
`$LS`. Those paths are provenance, not something a reader can resolve; the predicted structures
themselves are too large to deposit.
