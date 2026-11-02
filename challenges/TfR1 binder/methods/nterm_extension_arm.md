# Giving a cleared scaffold something to reach with: a generated N-terminal extension

⛔ **Read this first.** This document reports a measurement defect that was found in this arm's own
scoring script *after* the arm's results had already been written up once, and that invalidated the
ranking the rest of the arm was built on. The defect, what it changed and what it did not change
are in §5 and §6. Nothing has been removed; the superseded ranking is deposited beside the
corrected one.

⛔ **Nothing in this document is a binding result.** Every coordinate is a generator's or a
structure predictor's model, and nothing here has been expressed, purified or assayed.

⇒ ⭐⭐ **The round's conclusion, stated up front so the rest can be read against it:** the generated
extension **does** reach the receptor, both in the generator's own frame (§6) and in an independent
co-fold (§8.2) — which is the opposite of what the arm predicted. ⛔⛔ **But the co-fold also shows
the design-frame engagement does not predict the folded engagement (§8.3), that this metric's
draw-to-draw noise is about twenty times the difference it would need to resolve (§8.4), and that
the predictor does not reproduce the 143 residues the construct holds fixed, returning them
6.7–16.2 Å away (§8.5).** ⇒ the positive result is real and the instrument used to rank it is not
good enough to rank with. Both halves are the finding.

---

## 1. The question this arm exists to answer

`methods/core_scaffolded_design.md` ends on a specific, narrow success: holding *part* of a mined
scaffold rather than the whole chain takes the collision with the receptor from 679 overlapping
atom pairs to zero. The redesign clears the receptor.

⛔ But clearing a receptor is not touching one. A design that no longer collides can simply be a
design that no longer reaches, and the ten-way screen that produced those scaffolds measured
**clearance**, never **engagement**. So the obvious next question is whether the cleared construct
can be given something to reach with.

This arm asks it directly. Take one design from that screen, keep almost all of it fixed, and ask
the generator for a **25-residue N-terminal extension** whose only job is to span the remaining
distance to the receptor.

⭐ **The arm was launched with a falsifiable prediction attached, and the prediction was that it
would fail.** The geometry of the kept scaffold puts its own N-terminus on the side *away* from the
receptor interface, so the most likely outcome was an extension that grows into solvent and never
makes contact — a clean negative that would close the line. §6 reports what was measured instead.

## 2. The construct

| | |
|---|---|
| source design | one design from the core screen's arm `opt01` (batch 3, model 3) |
| its first 3 residues | **dropped** — the generated extension replaces them |
| held from the source backbone | residues **4–145**, sequence-fixed (142 residues) |
| held guidepost | **1 arginine**, atoms `NE,NH1,NH2` held at their specified positions |
| generated | **25** residues of N-terminal extension **+ 2** residues of linker = **27** |
| final binder chain | **170 residues** |
| receptor held | **154 residues** in two segments, `A193–298` (106) and `A332–379` (48) |
| generator contig | `25-25,E1,2-2,H4-145,/0,A193-298,/0,A332-379` |

⭐ The two receptor segments are held as a **non-contiguous pair**, not as one crop — a sphere cut
around the designed site rather than a span. The two chain breaks in the contig are written with
explicit separators, because a separator written without its comma is silently swallowed and
changes the output layout without any error.

`[M]` Measured on the input before any generation: the held guanidinium sits **3.03 Å** from the
receptor and its backbone alpha-carbon **3.46 Å**. ⇒ the pin is already in contact; the extension
is not being asked to create the contact, only to carry a chain to it.

`[M]` Threading the design's own sequence onto the kept backbone required **93** substitutions
relative to the raw design backbone, which is recorded because it is the one place this construct
is not simply a subset of its parent.

## 3. Generation, and what the sample size actually is

**96 designs: 4 shards × 3 batches × 8 designs.** Settings were inherited from the screen and the
resolved configuration was read back from the run rather than assumed.

⛔⛔ **n is 12, not 96.** The generator produces a batch from one contig instantiation, so designs
within a batch are **clustered** — they share a length draw and much of their trajectory. Every
statistic below is reported over 96 rows because that is what was measured, but any claim about
*independent* draws has to be read against **12 batches**. The deposited table carries a `batch`
column precisely so a reader can do that arithmetic themselves.

⛔ **The runs are unseeded and not reproducible.** Re-running this contig does not return these 96
backbones. The deposited coordinates and tables are therefore the record of the run, not a recipe
that regenerates it.

## 4. The gate chain, and the order it runs in

⭐⭐ **Backbone continuity is gated first, before any clash or contact metric.** This order is not
stylistic. A raw generated backbone can break into pieces, and a chain that has broken apart scores
*better* on a minimum-distance clash floor than an intact one does, because the broken ends drift
apart. Gating clearance before continuity rewards exactly the failure it is supposed to catch.

The conjunction, in order:

| clause | what it means | `[M]` result over 96 |
|---|---|---|
| 1. no alpha-carbon step > 4.5 Å | the chain is one piece | ✅ **96 / 96** |
| 2. whole binder ≥ 2.0 Å from the receptor | raw, unrelaxed backbones clear a floor | **66 / 96** |
| 3. generated segment ≥ 2.0 Å from the receptor | the new segment specifically clears it | **66 / 96** |

⭐ **96 of 96 continuous** is the arm's most solid number and the one least likely to be an
artefact: the construct holds 143 of 170 residues, so the generator is solving a short
interpolation problem rather than a fold.

⛔ **Clause 3 can never bind, and it is reported so a reader does not credit it.** The generated
residues are a subset of the binder chain, so the minimum distance over the subset is always at
least the minimum over the whole chain. Clause 2 therefore **entails** clause 3 identically, for
any cohort. `[M]` Measured: the invariant holds in **96 / 96**, and **0 of 96** designs pass the
first two clauses and fail the third. ⚠️ A redundant clause in a conjunction is not harmless — it
reads as three independent tests when there are two.

⛔ **The 2.0 Å floor is a floor for raw generator output, not a packing criterion.** These
backbones have had no relaxation of any kind. A separation of 2.07 Å between a model's backbone and
a receptor's surface is not a measurement of good packing; it is the smallest value this gate was
willing to accept. Any tighter interpretation needs relaxation first.

## 5. ⛔⛔ The measurement defect, in full

The scoring script derived the set of **generated** residues from the generator's index map — the
structure that records where each input token ended up in the output.

⛔ **That map names only tokens, and a generated contig run has no token.** The script took the
lowest and highest mapped positions on the binder chain as the span of the chain, and subtracted
the held residues from it. Because every mapped position is by definition a *held* one, the lower
bound of that span was **the first held residue**, and every generated residue below it was
excluded before the subtraction ever ran.

`[M]` On a real design: the map covers binder positions **26–170**, while the structure has **170
residues starting at 1**. ⇒ residues **1–25** — the entire N-terminal extension this arm exists to
test — were never in the "generated" set, which collapsed to the **2-residue linker at 27–28**.

⇒ Both generated-segment columns, the contact count and the minimum distance, described a
**2-residue linker** and were labelled as describing a **27-residue extension**.

⭐ **The tell was in every output and was read as unremarkable:** a column reporting the size of the
generated set held the constant value **2** in all 96 rows, printed every time. ⛔ A column that
never varies is a question, not a fact.

**The fix** is to take the generated set from the structure rather than from the map — the residues
the binder chain actually has, minus the residues the specification holds — computed where the
structure is in hand, and raising rather than continuing if that set comes back empty. `[M]` The
reported size of the generated set goes from **2** to **27** in all 96 rows.

## 6. ⭐ What the correction changed, and what it did not

| | before | after |
|---|---|---|
| continuity | 96 / 96 | **unchanged, 96 / 96** |
| conjunction survivors | 66 / 96 | **unchanged, 66 / 96, 0 dropped and 0 admitted** |
| the ranking | — | ⛔ **completely different** |
| previous 1st place | 1st | **19th of 66** |
| previous 2nd place | 2nd | **60th of 66** |
| previous 3rd place | 3rd | 4th |

⭐ **The survivor set is identical for a provable reason, not a lucky one:** the only gate clause
the defective columns fed was clause 3, and §4 shows clause 3 cannot bind for any cohort. The gate
never depended on the broken measurement. The **ordering** did, entirely.

`[M]` **On the corrected measurement, over all 96 designs**, contacts made by the generated segment
with the receptor under 4.5 Å run **28 to 190**, median **69**. Over the 66 survivors: **28 to
165**, median **52**.

⇒ ⭐⭐ **The arm's prediction is NOT supported.** The extension does reach the receptor, in most
designs, across an order of magnitude of engagement. The expectation in §1 — that an extension off
the far side of the scaffold could only grow away from the interface — is not what was measured.

⛔⛔ **And engagement is not affinity.** A contact count says atoms are near atoms in a model. It
says nothing about whether the complex forms, or how tightly. n is 12 batches. Nothing in this
section has been folded independently, let alone tested.

⚠️ Two further cautions on the ranking, both structural: the corrected top three draw on only
**two** distinct batches, so three panels are not three independent results; and the contact count
is a sum over an atom-pair distance matrix, so it counts **pairs**, not distinct atoms.

## 7. ⛔ Sequence design, and a selection made on the superseded ranking

Sequences were designed on the generated backbones with a solubility-weighted message-passing
network, cysteine omitted from the alphabet throughout, the held guidepost's identity fixed, and a
diversity subselection taking **20 sequences** for independent folding.

⛔⛔ **Those 20 sequences sit on 3 backbones that were chosen by the ranking §5 invalidated.** Under
the corrected measurement, `[M]` those three backbones rank **4th, 19th and 60th of 66**, with
generated-segment contact counts of **136, 92 and 33**. The corrected 1st and 2nd place backbones
are **not** represented in the folded set at all.

⇒ This is stated rather than fixed, because the cohort was built, guard-checked and folded as it
stands, and re-selecting it after the fact would leave a published cohort whose provenance did not
match its selection rule. ⭐ One consequence is worth keeping: because the three backbones span
contact counts of 136 / 92 / 33, the folded set spans the corrected ordering instead of sitting at
the top of it, so it can be read as a test of whether that ordering predicts anything downstream —
which is a better experiment than a top-3 cohort would have been.

⚠️ 20 sequences over 3 backbones is **not** 20 independent design results.

## 8. The independent co-fold, and ⛔ what it cannot tell you

Each of the 20 sequences is folded against the **whole receptor ectodomain** rather than the
two-segment crop the design gate used — a 640-residue single chain, with **three N-linked sugars
bonded** at their sequons, 5 diffusion samples per sequence, single-sequence mode for both chains.

⭐ **Folding against the whole ectodomain, not the crop, is deliberate and was settled earlier in
this project:** a 106-residue segment cut out of the ectodomain does not fold the way it does in
the intact receptor, and a crop co-fold measured the receptor's own fold degrading from 0.95 Å to
7.21 Å against its deposited structure. The design gate used the crop; the final test does not.

⚠️ **Both chains are folded in single-sequence mode.** For the de novo binder that is correct — a
designed sequence has no evolutionary family and a profile built for it is noise. For a real
640-residue receptor it is the conservative choice rather than the best one, and an alignment-based
profile would normally improve the receptor's own fold. The setting is recorded here because it was
inherited from the builder rather than chosen for this run.

⛔⛔ **THE POSE READOUT ON THIS TARGET CANNOT DISCRIMINATE, AND THIS WAS ESTABLISHED FROM CONTROLS
BEFORE ANY DESIGN WAS SCORED.** An earlier round folded three positive controls through this same
instrument, with the criterion fixed in writing beforehand:

| control | what it is | fold vs. its own deposited chain | **pose vs. its own deposited complex** |
|---|---|---|---|
| Machupo virus GP1 | real, natural, **deposited complex** | ⭐ **0.59 Å** | ⛔ **64.68 Å** |
| parvovirus B19 VP1u | real, natural, **deposited complex** | 6.01 Å | ⛔ **84.71 Å** |
| a deposited **designed** binder | real, deposited | 14.50 Å | ⛔ **78.56 Å** |
| that round's 156 designs | ours | 14.32 Å median | 65.28 Å median |

⇒ ⭐⭐ **The first row ends the argument.** The predictor reproduces that protein's fold to 0.59 Å —
it knows exactly what the molecule looks like — and still places it **64.68 Å** from where the
crystal puts it on this receptor. ⛔ A pose metric that fails a perfectly-folded natural binder by
65 Å cannot be used to judge ours, and no threshold was set on it. That clause was **reported, not
gated**, per the pre-registration's own rule that a clause which cannot separate controls is never
quietly widened until something passes.

⇒ ⛔ **So this co-fold is not a screen and no design here is selected by it.** What it can report
is per-chain-pair interface confidence and the predictor's own fold agreement — ⛔ never a global
interface score, which is confounded when the target is multimeric. The results table is deposited
with that caveat attached to it rather than in a footnote.

### 8.1 What landed

`[M]` **20 sequences × 5 draws = 100 predicted complexes, 100/100 written, 0 failed examples.**
⛔ The count is the proof, not the exit status: this predictor exits **successfully** having written
nothing when it runs out of memory, and the run record's own completeness field counts *prediction
directories* (20), not structures (100). Both were checked, and the archive's checksum was matched
against the record before anything was read out of it.

### 8.2 ⭐⭐ The generated segment does engage in an independent fold — in most draws, not all

| `[M]` over 100 draws | |
|---|---|
| generated-segment contacts with the receptor under 4.5 Å | **0 to 223**, median **38.5** |
| whole-binder contacts | 0 to 771, median **374.5** |
| ⛔ draws where the generated segment touches nothing at all | **23 / 100** |
| draws where the **whole binder** touches nothing | 3 / 100 |
| binder–receptor interface confidence (`pair_chains_iptm`) | 0.131 to 0.805, median **0.282** |
| predicted complex confidence | 0.737 to 0.810 |

⇒ ⭐ The extension is not an artefact of the generator's own frame: an independent predictor, never
told where the binder was meant to sit, still puts the generated segment in contact with the
receptor in **77 of 100** draws. ⛔ **But in 23 it puts it nowhere near**, and the median engagement
is roughly half the design-frame figure.

### 8.3 ⛔⛔ The design-frame ranking does NOT transfer to the fold

This is what the accident in §7 made testable, and the answer is negative.

| backbone | design-frame generated contacts | `[M]` fold mean | `[M]` fold median | `[M]` mean `pair_chains_iptm` |
|---|---|---|---|---|
| ranked 60th of 66 | **33** | 56.9 | 38.0 | 0.334 |
| ranked 4th of 66 | **136** | 56.5 | 33.5 | 0.327 |
| ranked 19th of 66 | **92** | 62.8 | 40.0 | 0.347 |

⇒ ⛔⛔ **Across a 4× range of design-frame engagement, the folded engagement is flat** — the three
means span 56.5 to 62.8, and the backbone with **four times** the design-frame engagement has the
**lowest** fold median. ⇒ nothing here supports using the design-frame contact count to rank
candidates.
⚠️ **n is 3 backbones.** This is not a correlation and no statistic is quoted for it. It is a
negative result at small n, which is weaker than a refutation and much stronger than nothing.

### 8.4 ⛔⛔ And the reason ranking fails is measurable: the metric's noise swamps the effect

| `[M]` | |
|---|---|
| spread of generated-segment contacts **within one sequence**, across its own 5 draws | min 14, **median 130**, max 223 |
| the entire **between-backbone** difference in mean engagement | **6.3** |

⇒ ⭐⭐ **The draw-to-draw noise on this metric is about twenty times the difference it would have to
resolve.** ⛔ At five draws per sequence, this quantity cannot order anything, and any ranking built
on it is reading noise. That is a property of the measurement, not of the designs — and it is the
finding most worth carrying out of this round.

### 8.5 ⛔⛔ The predictor does not reproduce the scaffold it was given

The construct holds **143 of 170** residues sequence-fixed. The independent fold was checked against
the design backbone over exactly those held residues, within the binder chain, so no cross-frame
superposition is involved.

`[M]` **held-residue RMSD to the design backbone: 6.7 to 16.2 Å, median 12.0 Å.**

⇒ ⛔⛔ **This caveat sits underneath every number in §8.2 and §8.3.** The engagement measured in the
fold is the engagement of a structure that is **not** the design — the predictor rebuilt the held
protein somewhere else. ⚠️ Holding a backbone fixed during generation does not oblige an independent
predictor to agree that the sequence folds that way, and here it does not.
⚠️ `[M]` The most confident interface in the whole set (`pair_chains_iptm` **0.805**, 93 generated
contacts) is on a structure **14.9 Å** from its design. ⇒ confidence and fidelity to the design are
not the same axis, and the highest-confidence result is not the most faithful one.

### 8.6 ⚠️ Two measurement cautions, both from this run

- `[M]` **`pair_chains_iptm` barely tracks engagement**: Pearson r = **0.178** against
  generated-segment contacts over 100 draws (and r = 0.300 for held-residue RMSD, i.e. *worse*
  agreement with the design goes with *more* contacts). ⇒ the two readouts are close to
  independent, and neither can stand in for the other.
- ⛔⛔ **The first-listed draw of each sequence is the best of its five, not a sample.** This
  predictor **sorts** its diffusion samples by its own confidence. `[M]` Taking only that draw moves
  the median interface confidence from **0.282** to **0.559** and median engagement from 38.5 to
  60 — a difference produced entirely by selection. ⇒ the deposited table carries every draw, and
  the two summaries are reported apart.

### 8.7 ⭐ A free reproducibility replicate, and what it says

One of these 20 sequences was also folded locally, on different hardware, with the **same seed and
the same flags**, as the probe that preceded the run.

| `[M]` same sequence, same seed, first-listed draw | local 16 GB card | data-centre card |
|---|---|---|
| interface confidence | 0.650 | **0.623** |
| generated-segment contacts | 59 | **25** |
| held-residue RMSD to design | 7.3 Å | **9.2 Å** |

⇒ ⭐ **The confidence score reproduces across hardware to about 4%. The contact count does not — it
differs by 2.4×.** ⛔ A seed does not make this pipeline reproducible across accelerators, and the
quantity that moves most is the one §8.4 already shows is too noisy to rank on. ⇒ consistent
evidence, from a second direction, that the engagement count is a weak instrument.

## 9. What is not settled

1. ⛔ **No binding measurement of any kind exists** for this arm or any other in this directory.
2. ⛔ **Whether the corrected ordering predicts anything downstream.** §8.3 tested the design-frame
   contact count against the fold at **n = 3 backbones** and found it flat; §8.4 shows the metric's
   draw-to-draw noise is ~20× the difference it would need to resolve. ⇒ ⛔ the honest statement is
   that **this measurement cannot answer the question**, not that the ordering is refuted. Settling
   it needs many more draws per sequence, or a less noisy quantity, or both.
3. ⛔ **Nothing has been relaxed.** Every geometry here is raw generator output gated at a 2.0 Å
   floor, and §4's caution applies to every distance quoted.
4. ⛔⛔ **Whether the independent fold is testing the designed construct at all.** §8.5 measures the
   held 143 residues coming back **6.7–16.2 Å** from the design backbone. ⇒ the predictor rebuilt
   the protein that was supposed to be fixed, so the fold's engagement is a different structure's
   engagement. Nothing here separates "the extension works" from "a different structure's
   N-terminus happens to touch the receptor".
5. ⚠️ **Effective n is 12 batches for the backbones and 3 backbones for the sequences.**
6. ⛔ **The receptor was folded without an alignment profile** (§8).
7. ⛔ **Why 23 of 100 draws place the generated segment nowhere near the receptor** is unexplained,
   and nothing here distinguishes a real bimodality from sampling noise.
8. ⚠️ **The run is not reproducible, for two independent reasons**: the backbone generation is
   unseeded, and §8.7 shows that even a seeded fold of the same sequence gives a 2.4×-different
   contact count on different hardware.

---

## Files

| | |
|---|---|
| arm specification, provenance and input coordinates | `data/nterm_extension_arm_20261102/` |
| corrected scores, all 96 designs | `data/j948_nterm_score_CORRECTED_GENSET_20261102.tsv` |
| corrected ranking, 66 survivors | `data/j948_nterm_score_CORRECTED_GENSET_20261102_RANKED.tsv` |
| ⛔ superseded scores and ranking, kept | `data/j948_nterm_score_SUPERSEDED_20261102*.tsv` |
| the 20 designed sequences | `data/j949_diverse20_v2_20261102.tsv` |
| co-fold inputs, 20 | `data/nterm_ecto_fold_inputs_20261102/` |
| co-fold results | `data/j958_nterm_ecto_fold_scores_20261102.tsv` |
| scorer, builder, selector, figure scripts | `scripts/j945_*`, `j948_*`, `j949_*`, `j956_*`, `j957_*` |
| viewer session, corrected top 3 | `structures/j956_nterm_top3_corrected_20261102.pse` |
| panel figure | `figures/nterm_top3_corrected_20261102.png` |
