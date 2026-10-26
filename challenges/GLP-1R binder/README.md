# GLP-1R — what the known binders told us, and what our designs did not do

James McCann and the Koder Group, The City College of New York.
Work dated 2026-10-16 .. 2026-10-24. Challenge 3 of the Anthropic × Adaptyv 2026 competition:
human GLP-1R, a binder that must also **activate** the receptor.

---

## ⛔ Read the limitations first. On this work they are most of the result.

1. **Nothing here has been tested at the bench. No number in this directory is an affinity.**
   Every coordinate is a model, and every confidence number is a structure predictor's output about
   its own prediction. None of it was measured on a protein.

2. ⛔⛔ **Half of what the challenge asks for is not addressed at all.** The challenge asks for a
   binder that **activates**. This directory contains **no activation measurement, no activation
   prediction and no activation score**, because — see `methods/README.md` §6 — there is no published
   method that predicts agonism for a new class B1 binder. The second objective is addressed only by
   a literature argument about architecture, and that argument says the architecture we chose is the
   wrong one for agonism.

3. ⛔ **This is not a design result. It is a calibration of other people's instruments on this
   target, using known binders.** Three separate comparisons, with three different evidence grades,
   which must not be blended into one claim:

   | | what it is | grade |
   |---|---|---|
   | **A** | two published de novo miniproteins with deposited cryo-EM structures, put through our own geometry programs | a **measurement**, on third-party experimental structures |
   | **B** | GLP-1 and liraglutide used as positive controls for a published confidence filter, against 960 folds of our own designs | **model on both sides** — a predictor's output either way |
   | **C** | what has been achieved computationally at this receptor, from the literature | **reading**, tagged source by source |

4. ⭐⭐ **The headline is about the filter, not about the designs.** A published five-criterion
   in-silico filter was run here with its own printed thresholds. MEASURED: its
   `complex pLDDT > 0.90` criterion **rejects the real GLP-1R agonists 27 of 27**, which they max out
   at **0.891**. ⇒ the conjunction as printed is **unusable on this target**, and any "none of our
   designs pass the filter" headline resting on it would be an artefact of a criterion that rejects
   the native ligand. ⛔ It is the pLDDT channel specifically that fails to transfer; the other two
   complex-side channels survive their controls cleanly.

5. ⛔⛔ **On the channels that do survive, our designs do not reach the native ligand.** MEASURED,
   against the worst single draw the native ligand produced under conditions **matched to ours**
   (no alignment, interface ipTM 0.9534, interface PAE minimum 0.607 Å): **0 of 960** of our folds
   reach it on ipTM and **0 of 960** on interface error. Against the broader control set's worst
   draw, **0 of 96 designs** reach it on both channels on any of their ten seeds.
   ⚠️ **But these are near misses, not different worlds** — our best ipTM is 0.9242 against 0.9457,
   our best interface error 0.745 Å against 0.754 Å. ⛔ Phrasing used elsewhere in this repository
   for a different challenge ("off by 5-10×") must **not** be carried over here.

6. ⛔ **Only three of the filter's five criteria are evaluated.** The other two — apo pLDDT and
   apo-to-complex Cα r.m.s.d. — need a second prediction of each binder **alone**, which does not
   exist for this cohort. ⚠️ That is stated rather than waved off: because the complex-side
   conjunction here is **not** a clean zero (2 of 960 folds clear the two surviving printed cuts),
   the apo half could in principle move a verdict.

7. ⛔ **There is no negative control.** The controls here are GLP-1 and liraglutide — both real
   agonists. They establish that the two surviving channels read **high** for a real binder on this
   receptor. They do **not** establish that those channels read **low** for a non-binder on it. ⇒ the
   discrimination is demonstrated in one direction only. Another challenge folder in this repository
   deposits both directions; this one cannot.

8. ⚠️ **The geometric reference is a proxy receptor, and it is labelled one throughout.** The two
   published GLP-1R binders have **no deposited bound structure** — checked. The two structures
   measured in §3 are of a *different* class B1 receptor from the same paper, lab and method. The
   strongest form of the geometric argument — *"here is what a working binder at **this** receptor
   looks like"* — cannot be made by anyone.

9. ⛔ **The predicted structures are not deposited** — 960 design folds plus the control folds. What
   is deposited is the scoring tables, which is what every number rests on. ⇒ the numbers here can
   be **read and audited as measurements of the predictions**, but they cannot be **re-measured from
   this directory** without re-running the prediction step. `methods/INPUTS.md` §5 and
   `methods/README.md` §7 state exactly what is and is not reproducible from what is here.

10. ⛔ **There is no figure.** Other challenge folders in this repository open with one; this
    material has none.

11. ⚠️ **No novelty verdict is reached, in either direction.** The three deposited fold tables do
    carry raw structural-search columns for each designed sequence — the best-matching database
    entry, coverage, identity, TM scores and an E-value — but they are deliberately **not** turned
    into a pass or fail, for two stated reasons in the tables' own header comments: which TM
    normalisation a 0.80 bar refers to was unresolved, and the database snapshot searched is dated,
    which makes any "nothing resembles this" reading a **lower bound** on resemblance rather than a
    measurement. ⛔ No novelty claim is made anywhere here. Another challenge folder in this
    repository does carry a full novelty measurement and a correction to it; this one does not.

---

## The one thing in here that is a measurement on an experimental structure

MEASURED with our own programs on the two deposited cryo-EM complexes, and with **the same code
paths** used on our own designs so the numbers are comparable by construction:

| | the two published minibinders | our 384 first-round designs |
|---|---|---|
| length | 65 and 66 | 75-105 |
| ⛔⛔ `min_cut` — long-range contacts crossing the cheapest interior cut | **13** and **13** | ⛔ **0, in 384 of 384** |
| `rg_ratio` — compactness against the globular expectation | **1.0424** and **1.0554** | median **1.5335** and **1.5522** |
| binder residues at the interface | 13 of 65 and 10 of 66 | — |
| contacting Cα deviation from a straight line | 4.35 Å and 4.91 Å | — |

⭐⭐ **A de novo miniprotein can grip a class B1 extracellular domain and still be a compact, closed,
foldable bundle.** At `rg_ratio` 1.04 these are more compact than a natural globular domain.
⇒ ⛔ **so our `min_cut` of 0 in 384 of 384 is a property of our designs, not of the site.** That is
the load-bearing inference in this directory, and it is an inference: it says the compact geometry is
reachable in this receptor family, so reaching it is not excluded by the target.

⛔ And the complement, MEASURED: pooling the three deposited fold tables (160 + 96 + 96 = 352),
**0 of 352** of our designed sequences folded alone reach a 0.8 structural match to their own
designed backbone. The ceiling is **0.6837**.

⛔⛔ **Neither of those is evidence about efficacy.** Both published designs are **antagonists**. See
`methods/README.md` §6, which is the part of this directory a reader should take most seriously.

---

## What is in this directory

| | |
|---|---|
| `README.md` | this file — the limitations, first |
| `methods/README.md` | the method and result in reader order, 9 sections. Every number marked **MEASURED** or **READ** |
| `methods/INPUTS.md` | the target, what must be obtained upstream, and what is deliberately not distributed |
| `methods/SCRIPTS.md` | what each deposited program does, what it expects, and what it refuses |
| `data/` | 8 scoring and gate tables, with `data/COHORTS.md` naming what each one is |
| `scripts/` | the 6 programs that produced every number here |
| `structures/` | one model, with `structures/STRUCTURES.md` stating what it is and how it is numbered |

⛔ There is no `figures/` directory, because there is no figure.

---

## A note on internal references

Several deposited programs carry identifiers like `J559`, `J637` or `J640`. Those point into this
group's internal note series — an append-only record where each entry states one settled question,
the measurement behind it, and the date. **They are not deposited.** They are kept rather than
stripped because they are the *provenance of a threshold*, and a threshold with its provenance
removed reads as an arbitrary constant. ⇒ read them as "this value was settled deliberately, on a
stated date, and the record exists." Where a sentence conveyed nothing at all without the reference,
the sentence was rewritten instead.

⚠️ **One qualification, so the convention is not oversold.** Where such an identifier referred to a
*section* of an undeposited note, the section number was removed and the substance kept, because a
section number that no reader can resolve is a dangling reference rather than provenance.

The marker glyphs (⛔ ⭐ ⚠️ ✅) are emphasis, carried over from the working documents as in the other
challenge folders. They rank nothing and encode nothing.
