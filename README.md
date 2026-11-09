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

The competition is over. Where each challenge ended up:

| # | target | entered | novelty filter |
|---|---|---|---|
| 1 | pH-dependent EGFR binder | **16 sequences**, 2026-10-07 | **11 of 16 passed** |
| 2 | conditional TNF-α binder | **20 sequences**, in two arms | **all 20 passed** |
| 3 | GLP-1R binder | **20 sequences** | **all 20 passed** |
| 4 | TfR1 binder | ⛔ **nothing entered** — the design chain was unfinished | — |
| 5 | KRAS G12D / HLA-A\*11:01 | **7 sequences**, 4 validated and 3 not promoted | **all 7 passed** |

⛔⛔ **No wet-lab results have come back for any challenge.** Nothing in this repository has been
tested at a bench, and every structure in it is a model. The novelty filter above is the organizers'
automated check on how distinct a submitted sequence is; ⛔ it is not a measurement of whether
anything binds.

⚠️ Challenge 4 is deposited even though nothing was entered, because its negative results are the
useful part.

**Challenge 5 — the KRAS G12D peptide on HLA-A\*11:01 — 7 sequences were entered, of which 4 clear
our robustness bar and 3 do not, and the deposit says which is which.** The target is a
peptide in an MHC groove, so a binder that grips the shared MHC surface and ignores the peptide
scores well on most structural metrics and is worthless; the gate is built to refuse that molecule.
It is a conjunction, pre-registered: contacts **to the peptide** above a floor, **and** the
binder↔peptide interface confidence above the same floor. ⛔ Never either limb alone, never a ratio.

⭐ **That floor is other people's binders** — real published binders for *other* peptide–MHC
complexes, folded against our target. Against it, **the mean contact margin over the whole screen is
negative at every sequence-design rank**: the typical design touches the peptide *less* than one of
those cross-paired non-binders does by accident.

⭐⭐ **The transferable finding is about robustness, not about these molecules.** Seven designs that
cleared the gate on one fold were re-folded to 20 draws across 4 seeds, and **three of the seven were
demoted** — and nothing available at the one-fold stage predicted which. Contact margin, the
pre-registered ranking statistic, was *anti*-correlated: the best margin in the cohort was demoted
and the worst validated. Corroboration across four different designed sequences at four seeds did not
help either. ⇒ **the re-draw is not a formality, it is the measurement**, and it is cheap. The design
with the highest peptide contacts in the entire cohort failed on the interface limb alone, at one
seed, with a value inside the non-binder band — ranked on contacts alone it would have topped the
workstream. ⚠️ Every demoted design still cleared 3 of 4 seeds, so the bar separates **robust from
marginal**, ⛔ not real from fake.

### What was actually run, in order

1. **Target assembled** — HLA-A\*11:01 heavy chain, β2-microglobulin and the KRAS G12D 10-mer
   peptide, with the binder designed as a fourth chain of 65–85 residues.
2. **Backbone generation — RFdiffusion3** against the full four-chain complex, in a pre-registered
   **2 × 2**: held motif × hotspot breadth, with binder length and recycling identical in all four
   cells so a difference between cells is attributable. ⭐ Every design that survived came from a
   **no-motif** arm.
3. **Design-stage gates** — ⛔ backbone continuity **first** (a clash floor applied before it
   *rewards* a broken chain), then raw-clash survival on unrelaxed output.
4. **Sequence design — SolubleMPNN**, several sequences per backbone, ⛔ cysteine excluded outright.
   ⭐ Ranks 1–5 were carried forward, not rank 1 alone: folding only the top sequence per backbone
   under-counts usable backbones, which has now replicated in three cohorts.
5. **Co-folding — Boltz**, each designed sequence re-folded **as the full complex** with the target
   unchanged.
6. **Scoring** — peptide contacts, and the **binder↔peptide pairwise** interface confidence
   (⛔ not the global one, which a four-chain target dominates). Each is measured against **its own
   fold shard's floor**, set by real published binders for *other* peptide–MHC complexes folded
   against our target.
7. **The gate** — the pre-registered **conjunction**: both limbs above that floor. ⛔ Never one limb
   alone, never a ratio.
8. ⭐⭐ **Validation by re-draw** — every design that cleared was re-folded to **20 draws across 4
   independent seeds**, each seed scored against its own floor. **4 of 4 → validated; 1–3 of 4 → not
   promoted.** ⛔ A design that failed its re-draw was not re-folded again, because re-running a
   failed test until it passes converts a pre-registered gate into a best-of-N.
9. **Submission file built by derivation, not by hand** — sequences re-derived from the scored
   tables at each design's own clearing rank, each resolved to exactly one generation round by
   lookup (⛔ design names are **not** unique across rounds), then checked for chain lengths, absence
   of cysteine, the 10–250 aa window and uniqueness, and cross-checked by an independent validator.

⛔ **Nothing in Challenge 5 was tested at the bench, and every coordinate is a model.** See
`challenges/MHC-KRAS binder/` — start with its `README.md`, then
`methods/robustness_and_tiering.md`.

**Challenge 4 — the transferrin receptor (TfR1) — was worked but ⛔ nothing was submitted.** The
deadline passed with the design chain unfinished, and the deposit is published anyway because the
negative results in it are the useful part. ⛔ **No Challenge-4 sequence was entered into the
competition, and none of this was tested at the bench.** If it had been submitted, the deposited
sequence cohort is what would have gone in.

![Six generated TfR1 binder designs against one copy of the receptor surface, each design in its own colour in a shared orientation.](challenges/TfR1%20binder/figures/00_overview.png)

*Six designs from the ten-arm core-scaffolded screen, against one copy of the receptor surface, each
in its own colour in a common orientation. ⛔ Every coordinate shown is a model.*

Three lines were run. The honest reading of the first two is the negative one; the third is the only
one that reached its own test:

- **Structural mining — net yield of usable scaffolds: zero.** The site was defined as six typed
  side-chain functional groups at fixed relative positions, and two structure databases were
  searched for chains that already present that arrangement. **1,314** chains matched inside the
  RMSD bar; only **31 rows (30 distinct chains)** present all six groups on a reachable face; and
  all four survivors that were actually superimposed put their own body where the receptor is. The
  tight head of the match list is buried, so **fit and accessibility are opposed** in this problem.
- **Core-scaffolded redesign — the binding test was never run on this cohort.** Ten held-element
  options were built and screened, **400** designs generated, **39** passed continuity, clearance
  and fold integrity as a conjunction, and **312** sequences were designed on those 39 backbones.
  ⇒ ⭐ **CORRECTED 2026-11-02:** this entry previously said the independent co-fold of those 312 was
  never done. That was stale. The cohort **was** co-folded — **780** predictions, **156** designs ×
  5 draws, scored, with known natural binders folded alongside as controls. ⛔ **It still shows
  nothing about binding**, but for a different and more useful reason: the controls land **64–85 Å**
  from their own crystal poses, the same band as the designs, so **the gate cannot discriminate**.
  The screen measured whether the construct *clears* the receptor, never whether it *reaches* one.
- ⭐ **A generated N-terminal extension does reach the receptor — and it was predicted not to.**
  One cleared design was held almost entirely fixed (**143** of **170** residues) and the generator
  was asked for a 25-residue extension to span the remaining distance. **96/96** backbones came
  back continuous, **66** cleared the clearance conjunction, and the generated segment makes
  **28 to 190** contacts with the receptor under 4.5 Å, median **69**. The arm was launched
  expecting an extension that grew into solvent, so the prediction failed in the useful direction.
  ⭐ **And it survives the independent test the other two lines never reached**: the 20 sequences
  were co-folded against the whole 640-residue ectodomain, **100/100** predicted complexes, and the
  generated segment is in contact with the receptor in **77 of 100** draws.
  ⛔⛔ **A contact count on an unrelaxed model is not an affinity**, and ⛔ a defect in our own
  scoring code — which measured a 2-residue linker while claiming to measure the 27-residue
  extension — is reported in full in `challenges/TfR1 binder/methods/nterm_extension_arm.md`,
  because correcting it reordered the ranking completely while leaving the gate untouched.
- ⛔⛔ **And the same co-fold disqualifies the measurement used to rank that result — which is the
  most useful thing this challenge produced.** The design-frame contact count **does not predict**
  the folded one (backbones spanning **136 / 92 / 33** design-frame contacts give flat fold means of
  **56.5 / 62.8 / 56.9**); the metric's spread **within a single sequence across its own five
  draws** has median **130** against a between-backbone difference of **6.3**, so at five draws it
  is reading noise; and the predictor **does not reproduce the 143 residues the construct holds
  fixed**, returning them **6.7–16.2 Å** from the design backbone. ⇒ the reach is real, the ranking
  is not supported, and the honest version of both is deposited.

![Three generated TfR1 binder designs, each with its N-terminal extension in grey against the receptor surface in teal, in a shared orientation.](challenges/TfR1%20binder/figures/nterm_top3_corrected_20261102.png)

*The extension arm's corrected top three. Purple is the protein the generator was given and held
sequence-fixed; **grey is the N-terminal segment it built**; yellow is the one held guidepost; teal
is the receptor surface. ⚠️ Ranks 1 and 3 come from the same generation batch, so these are two
independent draws and not three. ⛔ Raw generator output, gated at a 2.0 Å minimum separation and
not relaxed — which is why the grey threads inside the receptor's surface envelope.*

⛔⛔ **One construct defect is reported in full rather than quietly fixed**, because it is the most
transferable thing in the deposit: four of the six specified functional groups were identity-fixed
during sequence design, but two sat inside a held element whose *sequence* was released to the
designer. Those two were redesigned. At the position where an imidazole was specified, the first
cohort came back with an imidazole in **0 of 312** sequences, and **190 of 312** of those positions
became an acid — while the quality-control column, which checks only the four fixed positions, read
"held" for every sequence. A corrected cohort with all six identity-locked holds **312 of 312** at
both positions; both cohorts are deposited, with the per-sequence measurement that distinguishes
them. See `challenges/TfR1 binder/methods/core_scaffolded_design.md` §7.5.

The transferable lessons from both lines, written to be useful without the target, are in
`challenges/TfR1 binder/FINDINGS.md`. Start at `challenges/TfR1 binder/README.md`, which states the
limitations first.

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

The competition ran five challenges, released weekly from 2026-09-28, with wet-lab validation by
Adaptyv Bio and results published openly on Proteinbase.

⚠️ **The schedule changed partway through, which is why the later challenges here carry November
dates.** After Challenge 1, entrants reported that one week per challenge was too short, and the
organizers extended each challenge's window from one week to two. The releases stayed weekly, so the
windows overlap and the fifth challenge ran into November. Our own working time was **two weeks on
Challenge 1 and one week on each of the four that followed**. Designs selected for testing have their
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
  TfR1 binder/                  designs, methods and notes for challenge 4 (⛔ not submitted)
  MHC-KRAS binder/              designs, methods and notes for challenge 5
```

⚠️ Challenge 5's target is a **complex**, not a single protein, and its published title contains a
solidus, which cannot appear in a path. The directory is therefore `MHC-KRAS binder` — the one place
in this repository where the directory name is not the organizers' title character-for-character.

The directory name is the challenge title as the organizers published it, so it contains spaces: in a
URL they appear as `%20`, and a shell command that includes the path needs it quoted.

`PROVENANCE.tsv` lists every deposited file with its size and md5 as deposited.

## Scope and licence

This repository is narrower than the group's full research codebase. What is published for each
challenge is decided for that challenge.

This repository is released under the **MIT License** — see `LICENSE`. You may use, copy, modify and
redistribute what is here, including commercially, provided the copyright notice and the licence text
travel with it. The licence covers the material **we** wrote: the code, the documentation, the
figures, and the sequences and coordinates we generated.

⛔ It cannot and does not license what is not ours to give. Not everything needed to reproduce the
pipeline is distributable in any case — third-party model weights and third-party scoring code must
be obtained from their own upstream projects under their own terms, and the full set of raw predicted
structures is too large to ship. See `challenges/pH-dependent EGFR binder/methods/INPUTS.md`.

⚠️ A licence is a copyright instrument and says nothing about patents: it grants no patent rights and
reserves none. MIT in particular states no express patent grant either way.
