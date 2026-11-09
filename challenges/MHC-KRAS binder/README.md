# Challenge 5 — a binder for the KRAS G12D peptide on HLA-A\*11:01

De novo binders designed against a **peptide–MHC complex**: the 10-mer KRAS G12D neoepitope
presented by HLA-A\*11:01. The target is not a protein surface but a peptide lying in a groove, and
almost every residue of that groove is shared with every other peptide the same allele presents. So
the design problem and the *scoring* problem are the same problem: a binder that grips the MHC and
ignores the peptide looks excellent by most structural metrics and is useless.

⛔ **Every coordinate in this directory is a model. Nothing here was measured at a bench.**
⛔ **Seven sequences were entered. Four of them clear our robustness bar and three do not, and the
file says which is which.** Clearing a computational floor makes a design a candidate, never a
result.

## What was submitted

`submission_challenge5_20261109.csv` — 7 single chains, 66–82 aa, ranked.

| # | name | aa | seeds cleared | tier |
|---|---|---|---|---|
| 1 | `kras-mhc-mb-2-m3-r1` | 82 | 4 / 4 | **validated** |
| 2 | `kras-mhc-mb-11-m1-r2` | 75 | 4 / 4 | **validated** |
| 3 | `kras-mhc-mn-8-m2-r5` | 81 | 4 / 4 | **validated** |
| 4 | `kras-mhc-mb-6-m3-r4` | 82 | 4 / 4 | **validated** |
| 5 | `kras-mhc-mb-11-m2-r3` | 66 | 3 / 4 | ⛔ inside noise, not promoted |
| 6 | `kras-mhc-mb-19-m0-r2` | 76 | 3 / 4 | ⛔ inside noise, not promoted |
| 7 | `kras-mhc-mb-4-m0-r1` | 72 | 3 / 4 | ⛔ inside noise, not promoted |

⭐ **All 7 passed the organizers' automated novelty filter** (2026-11-09). ⛔⛔ That filter checks
only that a sequence is sufficiently distinct from known proteins to be worth testing — it says
**nothing** about binding, and the three designs that cleared only 3 of 4 seeds passed it just the
same. ⇒ **novelty and robustness are orthogonal here**, which is precisely why the tier of each
design is stated above rather than letting "submitted" stand in for "supported".

The four validated designs sit on **four distinct backbones** and come from **both** of the
no-motif design arms. `data/j1064_submission_provenance_20261109.tsv` maps every row back to its
internal design name, its sequence-design rank and the round that produced it.

⚠️ All seven fall in the 40–100 aa minibinder size class; **none** exceeds 100 aa, so only one of the
two judged size classes is occupied. That is a consequence of what passed, not a choice about what to
enter.

## The gate, and why it is a conjunction

Pre-registered before the scoring data existed. A design clears only if **both** hold:

1. its contacts **to the peptide** exceed the in-batch non-binder floor, **and**
2. its **binder↔peptide** interface confidence (a pairwise ipTM, not the global one) exceeds the
   same floor.

⛔ **Never either limb alone, and never a ratio of the two.** A ratio pays partial credit, and each
limb alone is defeated by a different failure: contacts alone reward a binder draped over the groove,
and interface confidence alone rewards a confident grip on the shared MHC surface.

### The floor is other people's binders

The non-binder floor is not a random baseline. It is **real, published binders for *other* pMHCs,
folded against our target** — molecules that are genuinely good at this class of job and have no
reason to recognise this peptide. They set peptide contacts of **9.6–12.6** and interface
confidence of **0.861–0.895**.

⇒ ⭐ **Across the whole screen the mean contact margin against that floor is negative at every
sequence-design rank.** The typical design touches the peptide *less* than a non-binder does by
accident. That one number is the most honest summary of the difficulty here.

The floor is replicated into every fold shard, so each shard is scored against its own floor rather
than a pooled one — the shards differ, and pooling them would smear the comparison.

## ⭐⭐ The finding worth carrying out of this challenge

**A design that clears the gate on one fold is close to a coin flip.** Seven designs were re-folded
from 5 draws to 20 draws across four independent seeds. **Three of the seven were demoted.**

And **nothing available at 5 draws predicted which way a design would fall:**

- ⛔ **not the contact margin** — it was *anti*-correlated here. The best margin in the cohort
  (**+6.6**) was demoted; the worst (**+0.2**) validated.
- ⛔ **not corroboration across sequences** — the design that failed had cleared four different
  designed sequences at four different seeds.
- ⛔ **not the per-draw spread.**

⇒ ⭐ **The re-draw is not a formality, it is the measurement.** The round that took the last three
designs to 20 draws cost **under one dollar** of rented GPU time in total, and it changed the verdict
on three of seven. ⚠️ We are deliberately not quoting a per-design figure: the folds are batched with
their controls, so dividing the round cost by the number of designs would invent a precision the
measurement does not have.

⚠️ Every demoted design still cleared **3 of 4** seeds and none collapsed to 0 of 4. The bar
separates **robust from marginal**, ⛔ not real from fake. That is why all seven were entered and why
the tier of each is stated rather than smoothed away.

### Why the best design failed, which is the mechanism

The design with the highest peptide contacts in the entire cohort passed the contact limb at **all
four** seeds — 15.6, 16.8, 17.2, 17.4 — and failed **only** the interface-confidence limb, at one
seed, where its value sat **inside the non-binder band**.

⇒ ⭐⭐ **Ranking on contacts alone would have made this the top design of the whole workstream.** It
makes plenty of peptide contacts and its interface is no more confident than a cross-paired
non-binder's. The conjunction earned its keep on the single design we had called strongest in
writing beforehand.

### A hypothesis we did **not** adopt

Sorted by interface confidence, the seven separate perfectly: every design that held scores ≥ 0.931,
every one demoted ≤ 0.916. The same ordering is present in the 5-draw data available at decision
time, while contact margin is fully interleaved there too.

⛔⛔ **This is post-hoc and is not used as a gate.** `n = 7`; a perfect split of this shape arises by
chance about 1 time in 35 before any penalty for the several statistics we inspected; it was not
pre-registered; and the separating gap at decision time is **0.0027**, far inside re-fold noise, so
any numeric cut there would be decided by noise. `methods/robustness_and_tiering.md` records it as a
prospective *ordering* test for a future cohort, with no threshold and explicitly not re-using these
seven.

## Layout

```
submission_challenge5_20261109.csv   the 7 entered sequences, ranked
data/                                scoring tables for every round, and submission provenance
figures/                             one rendered panel per submitted design, plus an overview
methods/                             how the designs were made, scored and tiered
scripts/                             the code that produced everything above
structures/                          PyMOL sessions for the submitted cohort
```

Start with `methods/README.md` for the pipeline and `methods/robustness_and_tiering.md` for what
"validated" is allowed to mean.
