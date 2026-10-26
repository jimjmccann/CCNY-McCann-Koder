# CCNY-McCann-Koder

Protein design work from James McCann and the Koder Group at **The City College of New York**, entered in
the [Anthropic × Adaptyv 2026 Protein Design Competition](https://proteinbase.com/competitions/anthropic-adaptyv-2026).

This repository holds the designs, the reasoning behind them, and the code used to produce them.

![The twenty submitted challenge-3 designs against the GLP-1 receptor. One panel per design, each in its own colour, all in a single shared orientation. The receptor is grey; the dotted horizontal layers mark the lipid bilayer boundaries, placed from OPM's curated assignment for PDB 6X18, with the extracellular side upward.](challenges/GLP-1R%20binder/figures/ccny_challenge3_submitted_designs_20261026.png)

*Challenge 3, the submitted twenty, against the GLP-1 receptor. One panel per design, each a
different colour, all in **one common orientation**. The receptor is grey; the dotted layers are the
**lipid bilayer**, taken from OPM's curated placement for 6X18 rather than drawn by eye, with the
**extracellular side up**. The designs cluster on the extracellular domain and the receptor mouth,
which is where they were built to bind. ⛔ Every coordinate shown is a model. Nothing here was
measured.*

## Status

**Challenge 3 — the GLP-1 receptor — has been submitted: 20 sequences, and all 20 passed the
competition's novelty filter.** They come from 20 distinct backbones across 11 generation arms, are
63–77 residues, and contain no cysteines.

How they were chosen is set out in
`challenges/GLP-1R binder/data/submission_20261026/SELECTION.md`. In short: 160 backbones were
scored on interface contacts, buried area per residue, interface predicted aligned error,
topological entanglement and co-fold confidence, ranked on the **minimum** of those percentile ranks
rather than the mean; the 40 highest were inspected individually; 11 were judged to engage the
receptor where the designs were built to reach, and the remaining 9 places were filled from the 120
not inspected, using a screen derived from the contrast between the two inspected groups.

⚠️ The honest reading of this set is weaker than the ranking makes it look, and the deposited files
say so:

- **The confidence axis is not a demonstrated predictor of binding.** On an earlier set of 248 of our
  designs with measured outcomes, the same panel reached a best AUROC of **0.601** with a confidence
  interval spanning 0.5. It was therefore given the least weight here, not the most.
- **The cohort's interfaces are predicted with low positional confidence throughout** — median
  interface predicted aligned error **18.0 Å**, best decile 11.8 Å. These numbers separate designs
  from each other; they certify none of them.
- **Entanglement and engagement are anti-correlated here, close to mutually exclusively.** Every
  uninspected design that beats the entanglement reference line set by two experimentally validated
  binders of a related receptor contacts **none** of the site residues the screen identifies. The
  submitted set prioritises engagement, and 9 of its 20 carry high linking numbers as a result. That
  is a property of the generation method and the clearest thing to fix next.
- **Nothing here addresses receptor activation.** This challenge asks for a binder that also
  activates; binding and activation are separate properties, and only the first was designed for.
- The receptor conformation is itself predicted and varies between runs by up to **9.5 Å**, so
  comparisons of where different designs sit on it are approximate.

Challenge 2 — the conditional TNF-α binder — was submitted: **20 sequences, in two arms, and all 20
passed novelty.** **17 are C1**, a single binder against the trimer; **3 are C3**, a *symmetric
homotrimeric* binder for a homotrimeric target. The result is negative and unvalidated, and we
submitted it anyway: our own predictor does not reproduce what the designs specify. For C1, **0 of
88** pass all five of the Baker *et al.* in silico criteria; for C3, **0 of 96** predicted trimers
are arranged as designed. That is one predictor's verdict, not a measurement — the test is
calibrated (a published de novo binder with measured affinity reproduces its own crystal pose to
**0.89 Å** through this identical path, and a negative control produces no confident interface) —
but a prediction is still not an affinity. Two of that challenge's three requirements, mouse
cross-reactivity and the pH switch, are not addressed. Start at
`challenges/conditional TNF-alpha binder/README.md`.

Challenge 1 has been submitted — some designs did not pass the novelty filter. **16 sequences were
submitted and 11 passed novelty**
(`challenges/pH-dependent EGFR binder/SEQUENCE_SET_CORRECTION_20261006.md` enumerates every set).
Deposited there: the sequence files, the submission file itself, a combined structure file holding
all eleven novelty-passed designs against the receptor with a per-design figure, the full method, the
scripts that produced and ranked the designs, and a written record of what the designs do **not**
address. Start at `challenges/pH-dependent EGFR binder/README.md`, which states the limitations
before anything else.

⚠️ Novelty is a statement about resemblance to known proteins, **not** about whether anything binds.

Nothing here has been tested at the bench. No number anywhere in this repository is an affinity.

## Novelty, and an accidental finding

While checking our designs against known proteins before submitting, we found that three of them are
close relatives of the published consensus tetratricopeptide repeat (CTPR) proteins — which are
themselves designed proteins. Nothing was copied, and no homolog was searched for at design time; the
resemblance was found by accident, by us, and is reported because it is useful to anyone running a
similar pipeline. The short note is
`challenges/pH-dependent EGFR binder/accidental_ctpr_homology.txt`; the measurements, scripts and raw
hit tables are indexed in `challenges/pH-dependent EGFR binder/NOVELTY_20261006.md`.

The follow-up is a negative result: we cut sequence identity by 30–45 points without touching a
single interface residue, and the fold did not move. Sequence novelty and structural novelty are
separable.

Our first structural-novelty figures for the single-site designs were measured on the wrong object — a
chain pulled out of a predicted complex, where the organizers fold the submitted sequence alone.
Re-measured on 693 predicted monomer structures, **0 of 693 draws reaches the high structural
similarity threshold**, one design previously flagged as high similarity is not, and no
nearest-neighbour reference survived the change. The superseded numbers are left in place with
correction banners rather than deleted, because what the earlier reasoning rested on is part of the
record. See `challenges/pH-dependent EGFR binder/NOVELTY_CORRECTION_MONOMER_20261006.md`.

These are our own measurements, with our own tools and thresholds, and they are a **lower bound** on
similarity: we searched two of the five databases the competition uses, and the sequence-novelty
result is "no statistically significant hit", which is not the same claim as "measured low identity".

## Competition context

The competition runs five weekly challenges from 2026-09-28 to 2026-10-31, with wet-lab validation by
Adaptyv Bio and results published openly on Proteinbase. Designs selected for testing have their
sequences, structures, methods and measurements published whether or not they work — negative results
included.

## Layout

One directory per weekly challenge. Each challenge gets its own target and its own approach, so
nothing is shared between them by default — no common `src/` that quietly couples one week's method to
another's. Directories are named after the target once it is announced.

```
challenges/
  pH-dependent EGFR binder/     designs, methods and notes for challenge 1
  conditional TNF-alpha binder/ designs, methods and notes for challenge 2
  GLP-1R binder/                designs, methods and notes for challenge 3
```

The directory name is the challenge title as the organizers published it, so it contains spaces: in a
URL they appear as `%20`, and a shell command that includes the path needs it quoted.

`PROVENANCE.tsv` lists every deposited file with its size and md5 as deposited.

## Scope and licence

This repository is narrower than the group's full research codebase. What is published for each
challenge is decided for that challenge.

There is no `LICENSE` file, so under default copyright all rights are reserved: the material may be
read, but not copied, modified, run or redistributed. Not everything needed to reproduce the pipeline
is distributable in any case — third-party model weights and third-party scoring code must be
obtained from their own upstream projects, and the full set of raw predicted structures is too large
to ship. See `challenges/pH-dependent EGFR binder/methods/INPUTS.md`.
