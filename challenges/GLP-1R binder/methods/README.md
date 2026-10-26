# GLP-1R, the known-binder work: method and result, in reader order

Written to be read by someone **outside this group holding only this directory**: a sentence that is
true only because our internal repository exists is a defect, not a shortcut.

James McCann and the Koder Group, The City College of New York.
Work dated **2026-10-16 .. 2026-10-24**. Challenge 3 of the Anthropic × Adaptyv 2026 competition:
human GLP-1R, a binder that must also **activate** the receptor.

⚠️ **Read `../README.md` first.** It states the limitations, and on this work the limitations are
most of the result.

Statements below are marked **MEASURED** (recomputed from the deposited files when this was written)
or **READ** (quoted from a dated internal record that is not deposited). Where a number could not be
traced to a deposited file, this document says so rather than supplying a rationale.

The marker glyphs (⛔ ⭐ ⚠️) are emphasis carried over from the working documents, as in challenges
1 and 2. They rank nothing and encode nothing.

---

## 0. Scope, and what is NOT here

This file documents **three comparisons against known binders of GLP-1R** and nothing else. It does
**not** document the design campaign itself — the generation, the site selection, the sequence
design or the backbone gates — which is a much larger body of work and is not deposited here.

⛔ **The compute-orchestration and launch layer is deliberately absent**, as in challenges 1 and 2.
Nothing in this file requires it: every program here reads structures that were already on disk, and
no prediction or generation was run for any number below.

⛔ **The predicted complex structures themselves are not deposited** — 960 design folds plus the
control folds. What is deposited is the scoring tables, which is what every number below rests on.
⇒ Stated plainly: the numbers in this file can be **read and audited as measurements of the
predictions**, but they cannot be **re-measured from this deposit** without re-running the
prediction step. §7 says exactly what a reader can and cannot reproduce.

---

## 1. The target, and the two objectives

Human GLP-1R, UniProt **P43220**, residues **24-432 = 409 aa**. That span covers the extracellular
domain, the seven transmembrane helices and the first part of the intracellular tail. All of it is
the organisers' own specification, read off the challenge page.

⭐ **The published construct sequence is byte-identical to UniProt residues 24-432** — MEASURED
against the UniProt record: 409 aa, no mismatches. It carries **no visible engineered mutations**.

⚠️ **But the assay construct is not necessarily that protein.** The challenge page states that the
supplier includes **proprietary thermostabilising mutations without listing them**. In class B1
receptors such mutations characteristically lock a conformational state, and state selectivity is
exactly what this challenge is about. ⇒ every active-versus-inactive statement anywhere in this
deposit must be read as *measured on structures, not on the assay construct, whose stabilising
mutations are undisclosed*. ⛔ We do not guess the positions.

⚠️ **One substitution, in the structure rather than in the sequence, and it matters for anyone
reusing our inputs.** The reference crystal structure used for geometry carries **phenylalanine at
residue 260 where UniProt carries leucine**. The convention we work to, and which the organisers
stated, is **identity from the sequence database, geometry from the structure**; our structural input
keeps the structure's phenylalanine, which sits 29 Å from the nearest designed contact and inside the
membrane. ⛔ This is recorded so that a reader who diffs our input against UniProt is not left
wondering; it is not a defect and it is not open.

**The two objectives are assayed on two different preparations.** Binding is measured on the
purified, detergent-solubilised, thermostabilised 409-residue construct. Activation is measured in a
cell-based reporter assay on **confirmed binders only** — receptor expressed at the surface of living
human cells, with the G protein present, reading downstream signalling through a cAMP-response
element to a secreted reporter enzyme. ⇒ three consequences, stated so nobody rediscovers them:

1. A design can pass the first objective on one protein and fail the second on another.
2. Affinity is the **entry ticket** and explicitly **not sufficient**.
3. The readout reports **signalling**, not conformation ⇒ a design that binds and holds the receptor
   in its inactive state is **indistinguishable from a design that does nothing**. An antagonist and
   a dud look the same.

---

## 2. Why the known binders are the instrument, not a footnote

Three questions could not be answered from our own designs, and each has a published answer:

| question | the known binder that answers it |
|---|---|
| Can a de novo miniprotein grip a class B1 extracellular domain **and still be a compact foldable protein**? | the two published minibinders with deposited bound structures (§3) |
| Can our structure predictor **express a confident interface at all** on this receptor, under our fold conditions? | GLP-1 itself, and liraglutide (§4, §5) |
| Does binding this receptor's extracellular domain **activate** it? | every assayed extracellular-domain binder in the literature (§6) |

⇒ ⛔ Without §4's control, "none of our designs pass the filter" is uninterpretable: it is equally
consistent with bad designs and with a blind instrument. The whole point of running the controls is
that they can come out either way, and on one channel they came out against the filter.

---

## 3. The published minibinders, measured with our own instruments

### 3a. What the reference is, and that it is a proxy

Muratspahić / Feldman / Kim *et al.*, *Nature* **656**:1044-1053 (2026),
doi `10.1038/s41586-026-10656-8` — "De novo design of miniproteins targeting GPCRs". Free author
manuscript at `PMC12157396`; preprint bioRxiv `10.1101/2025.03.23.644666`. READ from the author
manuscript, not the typeset version.

That campaign produced two functional GLP-1R binders with measured affinities of **27 nM** and
**5.3 nM**. ⛔ **Neither has a deposited GLP-1R-bound structure** — checked. The structural work
covers two other receptors. So the structures measured here are the **closest published analogue**:

| structure | complex | resolution |
|---|---|---|
| **9MM5** | CGRP receptor + minibinder `dC2_049` | 3.2 Å |
| **9MNI** | CGRP receptor + minibinder `dC2_050` | 4.1 Å |

⭐ The CGRP receptor is **class B1, the same family as GLP-1R, with a homologous extracellular
domain**, and these come from the same paper, lab and method. ⚠️ **It is a proxy and is labelled one**
in this file, and in the docstrings of the deposited programs that measure and draw it.

⛔ **Chains are named, never inferred**, from the structures' own entity description records:
chain **A** = the de novo minibinder, **R** = the receptor, **E** = the accessory protein. Non-A/R/E
chains are removed before anything is measured. A hardcoded chain pair has silently changed a
conclusion in this group's work before.

### 3b. What is measured, and with what code

`../scripts/j559_baker_reference_20261019.py`. **The same code paths used on our own designs**, so the
numbers are comparable by construction rather than by assertion:

| quantity | definition |
|---|---|
| `min_cut` | the fewest long-range Cα-Cα contacts (< **8.0 Å**, \|i−j\| > **4**) crossing any interior cut of the chain, with a **13**-residue margin at each end |
| `med_cut` | the median of the same cut profile |
| `rg_ratio` | radius of gyration divided by the globular expectation **2.2·N^0.38 Å** |
| helix signature | median Cα(i,i+4) and Cα(i,i+5) distances |
| interface | binder heavy atoms within **4.5 Å** of receptor heavy atoms, both residue sets reported |
| interface line-fit | a straight line fitted to the **contacting** binder Cα positions; the reported number is the r.m.s. deviation from it |

⭐ **`min_cut` is the question "is this chain one fold or a path?"** A single folded domain keeps many
contacts crossing everywhere, so no cut is cheap. A path, or two lobes on a tether, has a point
almost nothing crosses. ⛔ It is **margin-sensitive**, so it is never quoted across settings; every
value in this file is at margin 13.

⭐ **The interface line-fit is the question "is this a groove-follower or a face?"** A chain threaded
along a groove hugs a line. One that presents a face does not.

⛔ Downloaded structures are untrusted data: they live in their own directory, the program lives
elsewhere, paths arrive as arguments, and it is run with the interpreter's isolated mode.

### 3c. ⇒ The result

**MEASURED**, recomputed for this document from `../data/j559_baker_reference_20261019.tsv`, and the
"ours" column from `../data/j557_topology_deep6_20261019.tsv` (288) + `../data/j557_topology_shallow4_20261019.tsv` (96):

| | `dC2_049` | `dC2_050` | our 384 first-round designs |
|---|---|---|---|
| length | **65** | **66** | 75-105 |
| ⛔⛔ `min_cut` | **13** | **13** | ⛔ **0 — in 384 of 384** |
| `med_cut` | 32 | 37 | median 16 (`deep6`) / 17 (`shallow4`), max **51** |
| ⛔ `rg_ratio` | ⭐ **1.0424** | ⭐ **1.0554** | median **1.534** / **1.552** |
| Cα(i,i+4) / (i,i+5) | 6.17 / 8.66 | 6.16 / 8.58 | 6.16 / 8.62 |
| binder residues at the interface | **13 of 65 = 20.0%** | **10 of 66 = 15.2%** | — |
| contacting Cα deviation from a straight line | **4.35 Å** | **4.91 Å** | — |

⭐⭐ **So a de novo miniprotein can bind a class B1 extracellular domain and be a compact, closed,
foldable bundle** — at `rg_ratio` 1.04 it is **more** compact than a natural globular domain, for
which the same instrument gives 1.138 (READ, from the same instrument run on a reference protomer in
the challenge-2 work; not recomputed here).

⇒ ⛔⛔ **Therefore our `min_cut` of 0 in 384 of 384, and our `rg_ratio` of 1.53, are a defect of our
designs and not a property of the site.** That is the load-bearing inference in this section, and it
is an inference, not a measurement: it says the compact geometry is reachable in this receptor
family, so reaching it is not excluded by the target.

⭐ **The two independent published designs converge on the same receptor patch** — receptor residues
**41, 72, 92, 94, 114, 124, 125** are shared between them. `dC2_049` adds 34, 38, 93, 95, 116, 121;
`dC2_050` adds 128, 129, 131. ⇒ their footprint spans receptor residues **34-131**, which is the
peptide-binding groove region of a homologous extracellular domain.

### 3d. ⛔ Two things that are NOT the difference, and both are worth knowing

1. ⛔ **Not helicity.** Cα(i,i+4)/(i,i+5) is 6.17/8.66 for theirs and 6.16/8.62 for ours —
   indistinguishable. Both are helical. ⇒ the difference is **topology**, how the helices pack
   against each other, not secondary structure.
2. ⛔ **Not length, and it runs the opposite way.** Theirs are **65-66** residues; ours are 75-105.
   ⇒ the **shorter** chains are the compact ones. Any argument that a fold can be bought with more
   residues argues against this observation.

### 3e. ⛔⛔ What §3 does NOT establish

Both published designs are **antagonists**, working by occluding the peptide site. §3 measures
*"can a compact de novo bundle grip a homologous class B1 extracellular domain and still be a
foldable protein?"* — a binding-and-foldability claim that would stand unchanged even if class B1
activation were entirely unrelated across receptors. ⛔ **Nothing in §3 is evidence about efficacy
and it must not be read that way.** See §6.

---

## 4. The native agonists as controls for a published confidence filter

### 4a. The filter, and only its published thresholds

A recent de novo enzyme-design paper's round-3 in-silico filter, READ from that paper's methods and
extended-data figure:

> complex pLDDT > 0.90 · interface ipTM > 0.90 · **min interface PAE < 1.2 Å** ·
> apo pLDDT > 0.90 · apo-to-complex Cα r.m.s.d. < 0.8 Å

⛔ **Only the published thresholds are used.** That paper's own code was not obtained and not run.

⛔ **"min PAE" is ambiguous** — the paper prints the criterion without defining the reduction over
the interface block. **Both readings are computed and neither is hidden**: the literal **minimum**
element over the binder × target block, and the **mean** over it (what the de novo-binder literature
usually means by an interface PAE). §4d settles which one is defensible here, on the controls.

⛔ **Only the complex-side half of the filter is evaluated.** The apo half needs a second prediction
of each binder alone, which does not exist for this cohort. ⚠️ That is stated rather than waved off:
because the complex-side conjunction here is **not** a clean zero, the apo half could in principle
change a verdict.

### 4b. The cohorts, and ⛔⛔ the one assertion the whole comparison rests on

`../scripts/j637_baker_filter_c3_20261024.py`.

| cohort | n | what it is | alignment |
|---|---|---|---|
| our designs | **960** folds = **96** designs × **10** seeds | binder + the whole continuous 409-residue receptor | **none**, on both chains |
| control, GLP-1 | **21** folds | GLP-1(7-36) amide + the same 409-residue receptor | on |
| control, GLP-1 | **3** folds | an earlier run of the same pair | on |
| control, liraglutide | **3** folds | liraglutide + the same receptor | on |

⛔⛔ **The receptor sequence is asserted byte-identical across the design and control cohorts, read
out of the input files rather than eyeballed.** If it differed, a confidence gap would not be
attributable to anything. The program refuses to score otherwise.

⛔ **n is 96 designs, not 960 folds.** The ten shards are ten seeds of the same 96 designs, so a
per-fold pass rate counts each design ten times. Both are printed, and the **per-design best-of-10**
is the generous reading and the one to quote against the designs.

⭐ **Pairwise ipTM, never the global number.** The global value averages all chain pairs; on a
multi-chain target that confounds the binder's interface with the target folding well. Here the
target is a single chain so the two nearly coincide, but the pairwise value is the one that means
*binder against target* and it is what is used.

### 4c. ⛔⛔ A real bug, caught and recorded — the chain order is not the same in every cohort

This is deposited as a finding, not buried. The first version of the program took protein-aligned-error
rows `[0:nb)` as the binder for **every** cohort. That is true of our design inputs, where the binder
is the first chain, and **false of both control inputs**, where the 409-residue receptor is first and
the peptide second. ⇒ for all 27 controls it measured an **intra-receptor** block — adjacent
residues, error about 0.3 Å — and reported it as a spectacular interface.

**MEASURED:** the correction roughly **doubled** the controls' interface error, from a median
minimum of 0.335 Å to **0.687 Å**. ⚠️ It changed **no pass/fail verdict**, but it changed the
*separation* from "about 2.2× better than our best" to "**8% better than our best**" — which is the
difference between the near-miss reading in §4d and a blowout.

⛔ **The guard that failed was a plausibility check `20 < nb < 120` on the inferred binder length.**
The control peptide is 31 residues, so an implausible layout produced a perfectly plausible length. A
range check cannot catch a transposition inside its own range.

⇒ Fixed by passing the block offsets **explicitly per cohort**, never inferring them from position,
plus the receptor-identity assertion above. The superseded output is kept in the private record as
evidence of the retraction and is **not deposited**.

### 4d. ⇒ The result, and ⛔⛔ the headline is about the filter

**MEASURED**, recomputed for this document from `../data/j637_baker_filter_c3_20261024_v2.tsv`
(min / median / max):

| channel | published cut | our 960 folds | ⭐ controls (27) | controls passing |
|---|---|---|---|---|
| complex pLDDT | > 0.90 | 0.659 / 0.737 / **0.787** | 0.878 / 0.886 / **0.891** | ⛔⛔ **0 / 27** |
| interface ipTM | > 0.90 | 0.082 / 0.542 / **0.924** | 0.946 / 0.958 / 0.966 | ✅ **27 / 27** |
| interface PAE, **min** | < 1.2 Å | **0.745** / 7.796 / 27.978 | 0.552 / 0.687 / 0.754 | ✅ **27 / 27** |
| interface PAE, **mean** | < 1.2 Å | 6.295 / 18.210 / 30.166 | 4.623 / 5.282 / 5.786 | ⛔ **0 / 27** |

**⛔⛔ Finding 1 — complex pLDDT > 0.90 does not transfer to this target.** Real GLP-1R agonists, in
the pose the predictor builds against the crystal structure, **max out at 0.891** and fail the
criterion **27 of 27**. ⇒ the three-criterion conjunction is **unusable here**, and any
"none of our designs pass the filter" headline built on it would be an artefact of a criterion that
rejects the native ligand. ⭐ It is the **pLDDT channel specifically** that does not transfer, not the
whole filter — the other two survive their controls cleanly.
**Why:** complex pLDDT is a mean over the whole complex. The published cut was calibrated on a
**180-200-residue** enzyme; this complex is about **470 residues**, of which **409** are a receptor
with long disordered stretches that no binder can rescue. The cut is target-size-gated.

**⛔ Finding 2 — the mean reading of the interface error is dead here.** It fails **27 of 27**
controls, so it cannot be the paper's intent on a target this size. ⇒ ⭐ **the literal minimum is the
only defensible reading**, which is also what a separate route concluded on the challenge-2 cohort.

**⭐⭐ Finding 3 — on the two channels that do survive their controls, the separation is near-total.**
Judged not against the paper's cut but against **the worst real binder in our own control set** — the
honest bar, because it is measured on this target with this predictor:

| | control worst | our folds reaching it | our best | per design, best of 10 |
|---|---|---|---|---|
| interface ipTM | 0.946 | **0 / 960** | 0.924 | **0 / 96** |
| interface PAE min | 0.754 Å | **1 / 960** | 0.745 Å | **1 / 96** |
| **both at once** | — | — | — | ⛔⛔ **0 / 96** |

⇒ **No design reaches the control band on both channels, on any of its ten seeds.**
⚠️ **The margins are near misses, not different worlds** — our best ipTM 0.924 against 0.946, our
best interface error 0.745 Å against 0.754 Å. ⛔ Phrasing used for the challenge-2 result ("off by
5-10×") must **not** be carried over here.

**⭐ Two folds clear the printed cuts on the two surviving channels. MEASURED, named, and not
promoted:**

| design | seed | ipTM | PAE min | complex pLDDT |
|---|---|---|---|---|
| `arm_shard08_v1_deep6_a207_5_model_2` | shard06 | 0.909 | 0.843 Å | 0.724 |
| `arm_shard06_v1_deep6_a207_4_model_2` | shard07 | 0.910 | 0.853 Å | 0.740 |

⛔ **2 of 960 folds, and one seed of ten each** — i.e. neither design does it reproducibly, and both
sit below every control draw. A binary gate at this n returns almost all zeros and ranks nothing.
⇒ they are named because naming them is how the claim stays checkable, **not** as a selection.

---

## 5. The confound, and the measurement that killed it

### 5a. The confound, stated before it was resolved

Every control in §4 was folded **with** a receptor sequence alignment; all 960 of our folds are
**single-sequence, no alignment on either chain**, on a byte-identical receptor. Confidence channels
are alignment-sensitive, so the §4d separation confounded design quality with the fold condition.
**MEASURED** at the time: receptor-chain `ptm` median **0.860** for our folds against **0.977** for
the controls — the *same* 409 residues, predicted far more confidently in the control. ⇒ ⛔ if the
alignment moves the receptor's own confidence by 0.12, it is not safe to read the ipTM and error gaps
as design quality.

⇒ §4d's separation was therefore recorded as **real as a measurement and not yet attributable**.
⛔ That is the state of the record at `J638`, and §5b supersedes it.

### 5b. ⇒ The matched control, and it was already on disk

`../scripts/j639_msa_effect_on_baker_channels_20261024.py`. An earlier run of the **same** ligand
against the **same** receptor exists in two arms that differ in **exactly one factor**: the receptor
alignment on, and the receptor alignment set to empty. 21 draws each.

⭐ **The alignment condition is verified in the runs' own bytes, not read off the input file** —
MEASURED: the alignment-on arm's processed-alignment directory holds **11** files, the
alignment-off arm's holds **0**. And the receptor and peptide sequences are **byte-identical across
the pair** (409 + 30 residues, 0 mismatches) ⇒ the arms differ in the alignment alone.

**MEASURED.** Rows 1-2 recomputed from `../data/j639_msa_effect_20261024.tsv`; rows 3-4 from
`../data/j637_baker_filter_c3_20261024_v2.tsv`, where they are the `ctrl_glp1_msa` and
`c3_design` cohorts:

| arm | n | ipTM, median (range) | interface PAE min, median (range) | complex pLDDT, median | clears **both** surviving cuts |
|---|---|---|---|---|---|
| GLP-1, template, alignment **ON** | 21 | 0.934 (0.924-0.944) | 0.753 Å (0.652-0.847) | 0.886 | ✅ **21 / 21** |
| ⭐ GLP-1, template, alignment **OFF** | 21 | **0.961** (0.953-0.974) | **0.576 Å** (0.485-0.607) | 0.876 | ✅ **21 / 21** |
| GLP-1, no template, alignment ON | 21 | 0.955 (0.946-0.964) | 0.688 Å (0.552-0.754) | 0.886 | ✅ **21 / 21** |
| **ours**, no template, alignment OFF | **960** | 0.542 (0.082-**0.924**) | 7.796 Å (**0.745**-27.978) | 0.737 | ⛔ **2 / 960** |

⛔⛔ **Finding 4 — the alignment does not help. It slightly hurts.** Single-factor deltas on the
native ligand, medians, everything else held:

| factor | Δ ipTM | Δ interface PAE min |
|---|---|---|
| alignment **OFF** versus ON (template held on) | **+0.027** (better) | **−0.177 Å** (better) |
| template **ON** versus off (alignment held on) | −0.020 (worse) | +0.065 Å (worse) |

⇒ ⛔ **The hypothesis that our designs look bad because they lack an alignment is refuted in the
direction opposite to the worry.** Removing the alignment made the native ligand score *better* on
both channels. The same goes for the template, the only other confound left standing: forcing the
receptor structure made the native ligand slightly *worse*. ⇒ neither factor inflates the controls.

**⇒ Finding 5 — the comparison that now means something.** Against the matched arm's **worst single
draw** (ipTM 0.953, interface error 0.607 Å):

| | our 960 folds reaching it |
|---|---|
| ipTM ≥ 0.953 | ⛔ **0 / 960** |
| interface PAE min ≤ 0.607 Å | ⛔ **0 / 960** |
| both | ⛔ **0 / 960** |

**Not one of 960 folds reaches the worst draw the native ligand produced under conditions matched to
ours.** Our best values — ipTM 0.924, interface error 0.745 Å — fall outside the **entire 63-draw**
range of the native ligand on both channels.

### 5c. ⚠️ The caveats, stated rather than buried

1. ⛔ **Both arms of the matched pair carry a forced structural template**; our 960 folds carry none.
   The template factor is measured directly above and it moves the native ligand the *wrong* way, so
   it cannot explain the gap — but there is **no control arm that is both untemplated and
   single-sequence**, so that one cell of the two-by-two is **interpolated from two measured
   single-factor deltas, not observed**.
2. ⚠️ The matched pair carries the structure's leucine at residue 260 while our folds carry the
   structure's phenylalanine (§1). Identical across the pair, so the within-pair delta is clean.
   ⛔ Not a proposal to change anything.
3. ⚠️ The peptide is 30 residues in the matched pair and 31 in one of the §4 control arms, and the
   alignment in the matched pair is the ligand-free one. Neither affects the within-pair delta.

### 5d. ⇒ What §4 and §5 together settle, and what they do not

✅ **Settled — the predictor, the alignment condition and the template are all exonerated.** On this
receptor, single-sequence and untemplated, the predictor gives a real GLP-1R agonist ipTM about 0.96
and an interface error about 0.58 Å. **The instrument can express a confident interface under exactly
our conditions**, so "the instrument is blind on this target" is no longer an available explanation.

⛔ **NOT settled, and these matter more than the above:**
- **That these designs cannot bind.** A prediction is not an affinity measurement. The filter's own
  source paper needed a 185-variant wet screen for its last improvement.
- ⛔⛔ **That complex pLDDT > 0.90 can be used here at all.** §4d Finding 1 stands: it rejects the
  native ligand 27 of 27 on the untemplated arm and 21 of 21 on the templated one. ⇒ the conjunction
  as printed remains unusable; only the two-channel form is defensible.
- **Why** the designs fail — interface chemistry, the pinned residue set, the backbone, or the
  receptor conformation they were designed against. Nothing here distinguishes those.
- **That the two surviving channels are correctly calibrated** for a binder task on this receptor.
  The controls show those channels *discriminate*; they do not show that 0.90 and 1.2 Å are the
  right numbers here.

---

## 6. ⛔⛔ The architecture result from the literature: every assayed binder at this domain is an antagonist

This section is **READ**, not measured, and every claim in the full prior-art record is tagged by
whether we read a primary source, a secondary one, or neither. The summary:

- **The published de novo campaign at this exact receptor** screened about **17,000** designs by
  yeast display against the **soluble extracellular domain**, expressed **49**, and obtained **2
  functional antagonists** (reporter IC50 61 ± 20 nM and 39 ± 18 nM; affinities 27 nM and 5.3 nM).
  ⛔ **No agonist was reported and none was tested.** The hotspot was in the extracellular domain and
  the stated rationale was that binding there "should induce steric hindrance and prevent peptide
  interaction" ⇒ **they were not trying for agonism.**
- That paper's **agonists** came from a different receptor class by a different route — designed
  against active-state structures of the transmembrane bundle.
- **A second, independent, state-specific design pipeline** aimed at the **active** structure of this
  same receptor also got **antagonists** (best IC50 874 nM, 4 of 10 peptides inhibitory); its
  agonists came from two other receptors, both a different class. ⇒ **two independent campaigns, two
  different methods, the same split.**
- **The one confirmed de novo class B1 agonist** is a **cyclic peptide**, in a preprint, n = 3
  synthesised, designed against a **different** receptor — activity at GLP-1R was cross-reactivity
  (best EC50 31.9 nM). ⛔ There is **no de novo designed protein agonist for any class B1 receptor.**
- ⭐⭐ **And the mechanism says why.** The family binds peptides in two steps: the peptide's
  C-terminal half binds the extracellular domain and supplies affinity and selectivity; its
  N-terminal half inserts into the transmembrane bundle and **activates**. N-terminal truncation
  converts class B1 agonists into competitive antagonists — a decades-old, family-wide demonstration
  that **extracellular-domain engagement buys affinity and nothing else**.
- ⛔ The one apparent counterexample is not one: in the paper reporting an "intrinsic agonist
  mechanism" for this receptor's extracellular domain, the activating constructs were **covalent
  fusions** through a linker or a membrane tether; a linker carrying **no ligand sequence at all**
  still gave about 30% activity, and **free truncated peptides supplied separately were not
  functional**. ⇒ tethering and proximity, not extracellular-domain binding.

⇒ ⛔⛔ **Therefore §3's reference structures are usable for GEOMETRY — how to grip a class B1
extracellular domain and still be a protein — and ⛔ NEVER for EFFICACY.** On efficacy they are
antagonists and demonstrate the opposite of what this challenge asks for.

⚠️ **And a found erratum, flagged because anyone citing that source will trip over it.** The cyclic-
peptide preprint's methods names its design target as this receptor at a given structure and chain.
We checked the structure database entry directly: that structure is a **different receptor**, and the
named chain is not a receptor at all. The paper's own appendix sequence and its results are
consistent with the other receptor, so the methods sentence is a wording error — but a reader taking
it at face value will believe this receptor was the design target. It was not.

---

## 7. What a reader can and cannot reproduce from this directory alone

**Can:**
- re-derive every number in §3c, §4d and §5b from the three deposited scoring tables;
- re-measure the two published complexes' geometry with the deposited program, after fetching
  those two structures from the structure database by accession -- they are named in §3a and
  are **not** deposited here, and the program reads the chain assignments out of their own
  entity records rather than being told them;
- confirm the 0-in-384 `min_cut` result and the `rg_ratio` distributions from the two deposited gate
  tables;
- confirm the 0-of-352 fold-back result by pooling the three deposited fold tables (160 + 96 + 96);
- check both readings of the ambiguous interface-error criterion, and the per-channel pass flags,
  because both columns are in the table.

**Cannot, without obtaining third-party components and re-running:**
- reproduce the 960 predicted complexes — the structures are not deposited, and the predictions
  would need the predictor, the receptor input and the seeds;
- reproduce the control folds, for the same reason;
- reproduce the 384 backbones — the generator is **unseeded**, so re-running gives a different
  cohort rather than this one. The artefact is the record.

⛔ **No model weights are distributed** and not everything needed to reproduce the pipeline is
distributable. The backbone generator, the sequence model, the structure predictor and the molecular
viewer must each be obtained from their own upstream projects.

---

## 8. The honest result

- ⭐⭐ **The published filter cannot be used as printed on this target, and the known binders are
  what proves it.** One of its three complex-side criteria rejects the native agonists 27 of 27.
  Running the controls is the only reason that is known rather than assumed.
- ⭐ **The instrument is not blind.** Single-sequence and untemplated — our own conditions — the
  predictor gives a real agonist ipTM about 0.96 and an interface error about 0.58 Å.
- ⛔⛔ **Our designs do not reach it.** 0 of 960 folds reach the native ligand's worst matched draw on
  either surviving channel; 0 of 96 designs reach the control band on both, on any seed.
- ⭐ **The compact-foldable geometry is reachable in this receptor family** — the published
  minibinders are closed bundles presenting a small face — ⛔ **and our 384 are not it**: `min_cut` 0
  in 384 of 384, and 0 of 352 folds of their own sequences reach a 0.8 match to their own backbone.
- ⛔⛔ **And the literature says the architecture itself is the harder problem**: every free binder at
  this receptor's extracellular domain that has ever been assayed is an antagonist, and no de novo
  designed protein agonist exists for any receptor in this family.

⇒ **This is a known-binder calibration of a published in-silico filter on a class B1 receptor, which
found the filter partly inapplicable and the designs short of the native ligand on the part that
applies, plus a measured geometric reference for what a working design of this kind looks like.
⛔ It is not a binder, it is not an activator, and nothing in it was measured at the bench.**

---

## 9. Provenance of this documentation

Written **2026-10-25**. Every count in §3c, §4d, §5b and §7 was recomputed from the files in
this directory, not copied from a lab note; where a figure comes from a dated
internal record that is not deposited it is marked READ.

⚠️ **One number in §3c is READ and was not recomputed here**: the reference protomer's `rg_ratio` of
1.138, which comes from the same instrument run on challenge-2 material.

⚠️ **The internal record for §4 and §5 is a chain of two notes, and the first is partly superseded by
the second.** §5a is the superseded state, kept because the retraction is itself a finding: the
confound was real, was stated before it was resolved, and was then resolved against the worry. ⛔ One
further correction belongs with it: the earlier note proposed running the matched control as new
work, when a matched arm had already been run and was sitting on disk unexamined.

Several deposited programs carry identifiers like `J559`, `J637` or `J640`. Those point into this
group's internal note series, an append-only record where each entry states one settled question, the
measurement behind it, and the date. **They are not deposited.** They are kept rather than stripped
because they are the *provenance of a threshold*, and a threshold with its provenance removed reads
as an arbitrary constant. ⇒ read them as "this value was settled deliberately, on a stated date, and
the record exists." Where a sentence conveyed nothing at all without the reference, the sentence was
rewritten instead.
