# Robustness and tiering — what "validated" is allowed to mean

This is the part of the method we would most want another group to take.

## 1 · Draws, seeds, and why one fold is not a measurement

A co-fold is stochastic. A design scored on one fold of 5 draws has a number; it does not have a
measurement, because the within-fold spread understates how much the same molecule moves between
*independent* folds.

So every design that cleared the gate was **re-drawn**: re-folded to **20 draws across 4 independent
seeds**, each seed scored against **its own** floor. The seed is the shard — four shards, four seeds,
with the floor controls replicated into each.

## 2 · The tiers, pre-registered before any re-fold data existed

| tier | definition | may be called |
|---|---|---|
| **VALIDATED** | clears the conjunction at **4 / 4** seeds, 20 draws | validated |
| **CANDIDATE** | cleared at 5 draws only, never re-drawn | ⛔ must always be called **unvalidated** |
| **INSIDE NOISE** | clears **1–3 of 4** seeds at 20 draws | ⛔ **not promoted** |

⛔ A design that failed its re-draw is **not re-folded again**. Re-running a failed test hoping for a
pass converts a pre-registered gate into a best-of-N. In this campaign that rule is enforced by the
tooling rather than by restraint: the script that builds a re-draw cohort excludes any design that
already carries re-draw rows.

## 3 · ⭐⭐ The finding: a 5-draw clear is close to a coin flip

Seven designs were taken from 5 draws to 20. **Three were demoted.**

| round | held 4/4 | demoted 3/4 |
|---|---|---|
| first | 1 design | 1 design |
| second | 1 design | 1 design |
| third | 2 designs | 1 design |

And **nothing available at 5 draws predicted which way a design would fall:**

- ⛔ **Contact margin was *anti*-correlated.** The best 5-draw margin in the cohort (**+6.6**) was
  demoted. The worst (**+0.2**) validated.
- ⛔ **Corroboration across sequences did not help.** The design that failed had cleared **four
  different designed sequences at four different seeds** — which is evidence about the *backbone*,
  and it still failed.
- ⛔ **Per-draw spread did not help.**

⇒ ⭐ **The re-draw is not a formality; it is the measurement.** It is cheap — the round that re-drew
the last three designs cost under a dollar of rented GPU time — and it changed the verdict on three
of seven.

⚠️ **Every demoted design still cleared 3 of 4 seeds. None collapsed to 0 of 4.** ⇒ the bar separates
**robust from marginal**, ⛔ not real from fake. All seven were entered for this challenge, with the
tier of each stated; the three at 3/4 are not represented as validated anywhere.

## 4 · ⛔ A rank-1-only screen under-counts, by an amount we cannot estimate

Folding only the top-ranked designed sequence per backbone finds fewer usable backbones than folding
ranks 1–5. The **direction** has now replicated in three independent cohorts.

⛔⛔ **The magnitude has never been estimable in any of them.** In each cohort the implied factor has a
**denominator of 1 or 2** — it moves to a completely different number, or to infinity, if a single
design crosses the floor either way. ⇒ the factors our own tooling prints are explicitly **not
quotable** and no planning number should be taken from them. What is safe to carry is the direction:
**score more than one sequence per backbone.**

## 5 · ⚠️ A hypothesis recorded, deliberately NOT adopted

Sorted by binder↔peptide interface confidence, the seven re-drawn designs separate **perfectly**:
every design that held scores **≥ 0.931**, every one demoted **≤ 0.916**. Contact margin — the
*pre-registered* ranking statistic — does not separate them at all.

The ordering is also present in the **5-draw** data available at the moment a re-draw is chosen, so
it is not an artefact of the 20-draw numbers.

⛔⛔ **It is still not a gate, for four reasons stated so that a later reader cannot quietly promote
it:**

1. **`n = 7`.** A perfect split of this shape arises by chance about **1 time in 35** — before any
   penalty for the fact that several candidate statistics were inspected and the one that separated
   was reported.
2. ⛔ **It was not pre-registered.** What *was* pre-registered is contact margin, and the honest
   conclusion is that **the pre-registered ranking statistic is poor at this job** — ⛔ not that
   interface confidence is good at it.
3. ⛔⛔ **The separating gap at decision time is 0.0027**, which is far inside re-fold noise on this
   quantity. A numeric cut there would be decided by noise, and perfect separation across a gap that
   thin at `n = 7` is what overfitting looks like.
4. ⛔ A separation observed within the range of the data is not a lever until it **predicts**.

### ⭐ The prospective test, for whoever runs the next cohort

Stated here so it is falsifiable and so it cannot be softened later:

> Rank the candidates by **5-draw interface-confidence headroom above their own floor**. The designs
> that fail at 20 draws will be at the **bottom** of that ranking.

⛔ **No numeric cut is pre-committed** — the claim is about *ordering*, not a threshold.
⛔ **These seven may not be used as the test set.** They are the set the hypothesis was read off.
⚠️ With 3 failures among 7, a bottom-3 prediction is exactly right by chance about 1 time in 35, so a
single fresh cohort of this size is a real but **weak** test. **Two cohorts agreeing** is the first
point at which this would be worth acting on.
