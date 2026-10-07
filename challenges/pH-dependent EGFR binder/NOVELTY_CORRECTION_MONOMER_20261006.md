# Correction: the structural-novelty numbers were measured on the wrong object

**Measured 2026-10-06.** This file corrects two statements in
`NOVELTY_MEASURED_20261006.md` and one in `BRIDGE_MUTANT_FOLD_RESULT_20261006.md`. Those files are
left standing with correction banners rather than edited away, because a retraction is itself a
finding and the superseded numbers are what the earlier reasoning was built on.

**Nothing here is a binding measurement. No number in this directory is an affinity.**

## 1. What was wrong

Every structural-novelty number we had for the nine 75-80 aa single-site ("site A") designs came
from **chain A pulled out of a predicted binder-receptor complex**. The organizers' annotation
folds **the submitted sequence on its own**. Those are not the same object, and we had never
measured the one that is actually evaluated.

So we folded all eleven single-site molecules alone and re-ran the structural search.

## 2. What we measured

Eleven molecules x 3 seeds x 21 samples = **693 predicted monomer structures**, all 693 accounted
for by count rather than by exit code. Foldseek against the same local PDB structure set, same
convention as before: the best TM-score among hits covering more than 70% of the query; "high"
structural similarity is TM >= 0.80.

**0 of 693 draws reaches TM 0.80.** The single highest draw in the entire run is **0.738**.

| design | n | TM median | TM max | TM min | draws >= 0.80 | complex-derived prior | nearest (modal) |
|---|---|---|---|---|---|---|---|
| `egfr-siteA-id89` | 63 | 0.713 | 0.738 | 0.634 | **0/63** | 0.745 | `8jo3_A` |
| `cand-siteA-id2` | 63 | 0.644 | 0.734 | 0.633 | **0/63** | — | `8jo3_A` |
| `cand-siteA-id73` | 63 | 0.604 | 0.694 | 0.547 | **0/63** | — | `2crq_A` |
| `egfr-siteA-id3035` | 63 | 0.592 | 0.645 | 0.549 | **0/63** | 0.567 | `9c8w_A` |
| `egfr-siteA-id3036` | 63 | 0.592 | 0.668 | 0.516 | **0/63** | 0.646 | `2lse_A` |
| `egfr-siteA-id3017` | 63 | 0.591 | 0.643 | 0.547 | **0/63** | 0.745 | `9c8w_A` |
| `egfr-siteA-id3020` | 63 | 0.591 | 0.643 | 0.521 | **0/63** | 0.791 | `2lse_A` |
| `egfr-siteA-id2015` | 63 | 0.590 | 0.673 | 0.525 | **0/63** | **0.885 (flagged high)** | `2lse_A` |
| `egfr-siteA-id2048` | 63 | 0.588 | 0.672 | 0.509 | **0/63** | 0.757 | `2lse_A` |
| `egfr-siteA-id2046` | 63 | 0.587 | 0.669 | 0.517 | **0/63** | 0.682 | `2lse_A` |
| `egfr-siteA-id86` | 63 | 0.574 | 0.626 | 0.527 | **0/63** | 0.723 | `3key_A` |

Raw table: `data/siteA_monomer_20261006/monomer_foldseek_TM_20261006.csv`.

All eleven also **fold on their own** — they do not need the receptor to be structured. Every one
of the 63 draws per molecule is at or above 0.80 pLDDT; the weakest single draw in the whole 693
is 0.803 and the weakest molecule by median is `cand-siteA-id73` at 0.843.
Raw table: `data/siteA_monomer_20261006/monomer_plddt_20261006.csv`.

## 3. The three corrections

**(a) `NOVELTY_MEASURED_20261006.md` flags `egfr-siteA-id2015` at TM 0.885 as high structural
similarity, in its table of the twelve candidates. That does not reproduce.** On the monomer fold the
same design measures
0.590 median, 0.673 max, 0 of 63 draws at or above 0.80. **The high flag was an artifact of
measuring a chain extracted from a complex.** `id2015` is not in a different novelty class from
the rest of the set.

**(b) `BRIDGE_MUTANT_FOLD_RESULT_20261006.md` described the site-A designs as "TM 0.62-0.76" in its
§4. That range is not supported by either measurement** and should not be quoted. The
complex-derived table spans **0.567-0.885** — and silently omitted the one design flagged high.
The monomer medians span **0.574-0.713**. The *conclusion* that the single-site designs are
moderate rather than high on the structural half survives; the range used to state it did not
exist.

**(c) The nearest-neighbour references change completely when the object changes.** Not one
reference survives:

| | chain A extracted from a complex | monomer fold |
|---|---|---|
| nearest-neighbour set | `4n06`, `8j9b`, `8jo4`, `3t49`, `8ucf` | `2lse`, `9c8w`, `8jo3`, `2crq`, `3key` |

This is the part we think is worth other people's attention. A structural-novelty number is a
property of **the structure you submit for comparison**, not of the sequence. Six of nine designs
had been compared against the wrong reference, and nothing in the earlier run looked wrong.

## 4. What this does **not** change

- **The sequence half is untouched.** Zero statistically significant hits in Swiss-Prot 2026_03 or
  PDB seqres, by two independent search tools. See `SEQUENCE_NOVELTY_SWISSPROT_20261006.md`.
- **The three 145 aa bridge designs are untouched.** They were not re-measured here, and the CTPR
  homology described in `accidental_ctpr_homology.txt` stands exactly as written.
- **The `high_struct` flags in `data/novelty_20261006/verdict_submitted12.json` are not edited.**
  That file is a raw output of the complex-derived run and is kept as produced. Four of its entries
  carry `"high_struct": true`: `egfr-siteA-id2015` and the three bridges. **The `id2015` flag is
  the one corrected here. The three bridge flags stand.** See
  `data/novelty_20261006/VERDICT_SUBMITTED12_CORRECTION_20261006.md`.

## 5. Limits — read these before quoting any number above

- **Our TM convention is not necessarily theirs.** The "best TM among hits at >70% query coverage"
  rule, and the 0.80 threshold, are ours. A design at 0.738 is not far from 0.80 in anyone's units.
- **These are predicted monomers, not the organizers' own fold of them.** Closer to their pipeline
  than a complex-extracted chain. Still not their pipeline.
- **We searched a local PDB structure set only.** The annotation we are anticipating uses more
  structure databases than we did. More databases can only find more neighbours, never fewer, so
  every similarity figure here is a **lower bound** and every novelty figure an **upper bound**.
- **The sequence-novelty result is "no statistically significant hit", which is not the same claim
  as "measured low identity".** We searched **two** of the five databases the organizers use.
  Patent sequence collections are the unsearched category we would worry about first.
- **We make no claim about what the competition's annotation will report**, and nothing here has
  been tested at the bench.
