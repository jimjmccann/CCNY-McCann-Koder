# TfR1 — a designed binder that binds the wrong site, a structural-mining line that closed, and a scaffolded redesign

James McCann and the Koder Group, The City College of New York.
Work dated 2026-10-19 .. 2026-11-01. Challenge 4 of the Anthropic × Adaptyv 2026 competition:
human **transferrin receptor 1** (TfR1), the cell-surface receptor that imports iron.

---

## ⛔ Read the limitations first. On this work they are most of the result.

0. ⛔⛔ **Nothing from this challenge was submitted to the competition.** The deadline passed with
   the design chain unfinished — specifically, with the independent co-folding test still not run
   — and no TfR1 sequence was entered. This directory is published anyway because the negative
   results in it are the part worth having. The deposited sequence cohort is what would have been
   entered had the chain finished.

1. **Nothing here has been tested at the bench. No number in this directory is an affinity.**
   Every coordinate is a model or a database entry, and every confidence number is a structure
   predictor's statement about its own prediction.

2. ⛔⛔ **The first design round failed on targeting, and it failed cleanly.** Of **2,160**
   independently predicted binder–receptor complexes, **2,135 (98.84%)** reproduced **none** of the
   four designed contacts, and **0** reproduced three or four of them. The designs are not junk —
   they fold and they dock, in interfaces with interpenetrating side chains and no backbone clashes
   — they dock **somewhere else**: a single 36-residue patch, the same one in every generation arm,
   which is TfR1's own transferrin-binding site, **59–65 Å** from the site that was designed for.

3. ⛔⛔ **The score everyone ranks binder designs on prefers the wrong answer here.** Interface
   predicted TM-score correlates **+0.321** with being *further* from the designed site, and **0 of
   the top 20 predictions by that score** touch the designed site at all. Ranking this cohort on it
   would have selected hardest against the thing being made.

4. ⭐⭐ **The most likely cause is a specification error of ours: the binder was asked to be too
   short for the site it was asked to present.** Two independent routes agree. From the designs: a
   median-68-residue chain was asked to span **22.05 Å** when such a chain presents about
   **15.8 Å**. From the structural record (this directory's main content): of **795** PDB-derived
   chains that do present the six-group constellation specified here, **only 2 are under 120
   residues** and **91% are over 250**; across both databases searched, the shortest chain that
   presents it at all is **116 residues**. ⚠️ Marked **likely, not proven** — no run has tested it
   directly.

5. ⛔ **The structural-mining line documented here produced a statistic and no usable scaffold.**
   Net practical yield as drop-in scaffolds: **zero**. Three findings stack to close it, and they
   are set out in `methods/structural_mining.md` and `FINDINGS.md`. The negative results are
   the part of this directory worth reading.

6. ⭐ **A redesign that keeps only *part* of one of those scaffolds does clear the receptor, and a
   ten-way screen of which part separates cleanly.** This is the one positive structural result in
   the directory: holding subsets of secondary-structure elements, rather than a whole chain, takes
   the collision with the receptor from **679** heavy atoms to **zero**, and the arms that hold one
   particular large helix are the top four with **no overlap** against those that do not
   (`methods/core_scaffolded_design.md`). ⛔⛔ **It is still not a binding result.** 312 designed
   sequences exist and **none has been folded against the receptor or scored** — which is exactly
   the test the first round failed. ⚠️ And of the six functional groups the site specification
   names, only **four** remain identity-constrained through that construct; see §4 of that document.

7. ⛔ **A self-inflicted methodological error is corrected in place rather than removed.** An
   earlier conclusion that the six-group constellation is "unsatisfiable in nature, 0.00%" was an
   artefact of a pre-filter set tighter than the matching criterion algebraically allows. The
   corrected result is that roughly **1%** of natural chains present it. The superseded statement,
   the algebra that refutes it, and the re-derived numbers are all kept.

8. ⛔ **The mining line itself changed no design.** It characterised the instrument and the
   structural record; it moved nothing toward a better structure, and that is stated as a process
   outcome rather than softened (`FINDINGS.md`, §11). ⚠️ **Its closing result did then drive a
   reframing**, which is item 6 — so the line was not inert, but what it delivered was the argument
   for a different construct and not a structure.

9. ⚠️ **Two requirements of a real reagent are untouched here:** nothing in this directory
   addresses expression, solubility in practice, oligomeric state, or any functional consequence of
   binding TfR1. ⚠️ The sequence-design step reports composition, charge and isoelectric point,
   which are **predictions about** solubility and not measurements of it.

---

## The target, and the site

**TfR1** sits in the membrane as two identical subunits side by side. The site chosen was the
**apical patch** — a face on the top of the receptor — shifted toward the sub-patch that is
conserved between human and mouse. It was chosen because several natural ligands bind at that
patch, and because the conserved sub-patch carries chemistry a designer can aim at.

⚠️ **The tension in that choice was measured at the time and is not retrospective.** Every natural
ligand at the apical patch places most of its buried surface on residues that *differ* between
human and mouse. The conserved sub-patch is the narrow part of the target, not the comfortable part.

**The four receptor residues aimed at:** GLU294, LYS344, GLU350, LYS374.

**What the binder was asked to present:** six side-chain functional groups, held at fixed relative
positions — **four arginine guanidinium groups and two histidine imidazoles**. Only the functional
tip atoms were constrained (`NE,NH1,NH2` for arginine; `ND1,NE2,CG` for histidine); the backbone
was left free, deliberately, and is an output of the design rather than an input.

⭐ **That six-group constellation spans 32.55 Å**, and that number is the hinge of everything in
this directory. ⛔ No sub-constellation escapes it: the most compact five-group subset spans
**30.27 Å** and the most compact four-group subset **27.21 Å**, so dropping two of the six buys
**5.3 Å**. The site is intrinsically spread out.

---

## What was attempted

### Round 1 — de novo binders, and the targeting failure

| step | what it does | outcome |
|---|---|---|
| backbone generation (RFdiffusion3) | invents a backbone — a shape with no sequence yet — with the chosen contact points pinned in place | **640** designs in 4 arms |
| geometry quality control | discards shapes that clash, break, or misreport their crop | **577** survive |
| sequence design (SolubleMPNN) | writes an amino-acid sequence expected to fold into that shape, biased toward staying soluble | **4,616** = 577 × 8 |
| sequence quality control | composition, no cysteines, pins held | passes |
| independent co-folding (Boltz-2) | a separate structure predictor is given the binder sequence **and** the receptor together, and is **never told where the binder was meant to sit** | **432** sequences × 5 draws = **2,160** structures |

That last step is the test, and limitation 2 above is its result.

### Round 2 — scaffolded redesign, after the mining line closed

| step | what it does | outcome |
|---|---|---|
| scaffold decomposition | assign the chosen natural scaffold's secondary structure and enumerate which element subsets to hold | **10** core options, each sized to 120–160 residues |
| input build and validation | write the generator input and specification per option, then verify the five properties specific to this construct | 10 arms built; 5 properties checked directly, because two generic checks **did not apply** rather than passing |
| rate-measuring test run | one arm, to measure throughput and record the resolved configuration before committing to the screen | **16/16** designs; configuration identical to the screen that followed |
| backbone generation, all ten arms | the receptor present; everything between held elements generated from scratch | **400** designs, 40 per arm, 400/400 resolved |
| backbone gates | continuity **first**, then raw clearance, then fold integrity — as a conjunction | **39** designs pass all three; the arms separate |
| sequence design | write sequences for those 39 backbones, with a composition bias **refitted on this cohort** | **312** sequences; ⛔ the first cohort held 4 of the 6 specified groups, the corrected cohort holds 6 of 6 (§ below) |
| independent co-folding | — | ⛔ **NOT RUN.** This is the test, and it has not been done |

⇒ ⛔ **Round 2 stops one step short of its own test.** Everything in it is a statement about
backbone geometry or amino-acid composition.

---

## What was found, in the structural record

The question the mining line was built to answer: **do natural proteins already present this
constellation of functional groups, and if so, can one be used as a scaffold?**

⭐ **MEASURED, two databases, acceptance bar 2.0 Å throughout:**

| functional groups matched | PDB-derived store (85,699 chains) | predicted-model store (49,481 chains) |
|---:|---:|---:|
| 3 | 95.96% | 95.30% |
| 4 | 69.64% | 71.29% |
| 5 | 21.51% | 22.00% |
| **6 (all)** | **0.93% — 795 chains** | **1.05% — 519 chains** |

⛔⛔ **Carry the bar with every count above.** The best six-group match that exists anywhere in
either store is **1.068 Å** (PDB-derived) and **1.108 Å** (predicted models), against a 2.0 Å
acceptance bar. At a 1.0 Å bar **both stores return zero**. The geometry is satisfiable, but only
loosely.

Then three findings that stack, and together they close the line:

- ⛔ **Size and fit are opposed.** Nothing short presents the geometry well. The best-fitting chain
  at 250 residues or under is **1.4662 Å**, against the list's overall best of 1.0680 Å.
- ⛔⛔ **Accessibility and fit are opposed too.** Of **1,314** matching chains across both stores,
  **31 (2.4%)** present all six groups on a reachable face; **1,140 (86.8%)** present them where
  nothing can reach them. And the best fits are exactly the ones that fail: best match among
  failures **1.0680 Å**, best among survivors **1.5400 Å**. The whole tight head of the hit list is
  buried.
- ⛔⛔ **The survivors physically collide with the receptor.** Placed on the constellation, each of
  the four distinct small candidates puts **679 / 1,271 / 2,469 / 1,694** heavy atoms within 2.5 Å
  of the receptor crop, with closest approaches of **0.40 / 0.15 / 0.15 / 0.07 Å**. These are not
  grazing contacts; the candidate's body sits where the receptor is.

⭐⭐ **The honest reading.** The accessibility test asked whether a partner could *reach* the six
groups. It never asked whether the host protein's **body stays clear of the target**. That is a
strictly stronger requirement, and all four candidates fail it. ⇒ ⛔ **this hit list is not a
scaffold library, and scaffold grafting from it is not a shortcut to a better structure.**

⇒ ⭐ A reading consistent with both oppositions: a 32.55 Å constellation of four cations and two
imidazoles is, in nature, usually an **interior or active-site arrangement** rather than a surface
patch.

---

## What was done about it: keep part of a scaffold, generate the rest

If no natural chain can be used **whole**, the next question is whether **part** of one can. The
best-fitting survivor — 229 residues, matching at 1.5400 Å — was decomposed into its 16
secondary-structure elements, and subsets of those elements were held as a single rigid motif in
their native relative geometry while a backbone generator rebuilt everything between them **with the
receptor present**. Full method: `methods/core_scaffolded_design.md`.

⭐⭐ **MEASURED — the result that justifies the reframing:**

| placement of the scaffold on the constellation | heavy atoms within 2.5 Å of the receptor |
|---|---|
| the whole chain | **679** |
| the whole chain, both termini trimmed | **679** — trimming removes **none** |
| held element subsets, each re-placed | **0** for nine of the ten arms screened |

⇒ the collision is in the **middle** of the chain, so no trimming reaches it. Keeping elements and
discarding the rest is what takes it to zero.

⭐⭐ **And the screen separates the ten options.** 400 designs, 40 per option, ranked on a
**conjunction** of continuity, raw clearance and single-fold integrity — never a mean of fractions:

| arm | holds the large helix (residues 94–116)? | designs passing all three, of 40 |
|---|---|---:|
| best four | **yes** | 16, 12, 11, 6 |
| the other six | **no** | 5, 3, 3, 1, 0, 0 |

⇒ ⭐⭐ **the worst arm that holds that helix beats the best arm that does not, with no overlap
between the groups.** ⛔⛔ **How strongly to believe it:** this is 4 arms against 6, and the effective
sample size per arm is **5 batches, not 40 designs**, because one chain layout is drawn per batch.
Treat *"hold one large helix"* as strongly suggested and worth designing around, **not** as
established by one screen.

Two further structural readings, both stated as directions rather than rules:

- ⭐ **Fewer, larger elements connect; many small scattered ones fail — in two different ways.** One
  arm holding seven scattered elements passes continuity at 38/40 and then scores **0/40** on fold
  integrity: the chain threads between the fragments as two lobes on a tether instead of folding as
  one body. Another, holding three separate strands, fails the other way with the worst continuity in
  the cohort. ⇒ continuity and fold topology are **orthogonal** and must both be measured.
- ✅ **The deliberate control behaved as predicted.** The one arm built around the element known to
  sit inside the receptor's volume scores 0/40 with *perfect* continuity and a median backbone
  distance to the receptor of 0.379 Å. Its input carries 379 close atom pairs against the receptor, recorded before
  anything ran ⇒ it is a positive control that the clearance column reads what it claims to read,
  and that element is confirmed as the part to drop.

⭐ **The construct is also a large improvement on continuity** over an earlier attempt that delivered
the core as an obstacle the chain was never asked to connect: **7.2%** of designs here carry a
backbone break over 10 Å, against **69%** there.

**Sequence design** was then run on the 39 designs that passed all three criteria, giving **312**
sequences: 0 cysteines, all **identity-fixed** positions held, all distinct, 141–148 residues, 0 of
312 below the confidence health floor.

⛔⛔ **One correction, and it is a construct defect rather than a statistic.** "All fixed positions
held" means the **four** single-residue guideposts. The specification names **six** functional
groups; the other two sit inside a held element whose *sequence* was released to design, so they
were never identity-fixed. Measured per sequence: the position where an imidazole was specified
came back as one in **0 of 312** of the first cohort, and **190 of 312** of those positions became
an acid. A corrected cohort with all six positions identity-locked holds **312 of 312** at both,
and is deposited alongside the first. Full measurement and the reasoning:
`methods/core_scaffolded_design.md` §7.5, with the per-sequence table in
`data/j953_contained_pin_check_20261102.tsv`.

Two further findings from the cohort:

- ⚠️ **Sequence-design confidence does not separate the arms** — three arm medians inside 0.010 of
  each other with overlapping ranges. ⛔ It must not be quoted as if it confirmed the backbone
  ranking, and it is not evidence of binding.
- ⭐⭐ **A composition bias settled on an earlier cohort would have been wrong to inherit here**, and
  measuring this cohort instead showed why: its excess is **glycine** (13.18% against 7.07% natural),
  not the alanine the earlier bias targets, and two residues the earlier bias suppresses are already
  *below* natural here. Refitted on this cohort, glycine moves only 0.60 points under a large
  adjustment while alanine moves 1.71 under a smaller one ⇒ **the glycine demand is structural**, and
  a composition knob cannot argue these backbones out of it.

⛔⛔ **What none of this is.** No sequence has been co-folded against the receptor and none has been
scored for interface quality, so nothing above bears on whether any of these designs binds. The
clash figures are **raw** generator output, which clashes normally. The next real test is the co-fold
— the same test the first round failed.

---

## What is NOT settled

1. ⛔ **Whether a 1.1–1.9 Å constellation match is good enough to build on.** A design-and-relax
   pass can move a backbone, so a loose match is not self-evidently useless. No number is proposed
   here, in either direction.
2. ⛔ **Whether a candidate's accessible face is large enough to bury a partner at all**, and
   whether the six residues remain mutually compatible once a partner is docked. Unmeasured.
3. ⛔ **The four surviving candidates have not been looked at as structures** — only as geometry.
   Nothing is known here about their function, oligomeric state or expressibility.
4. ⛔ **Whether 2.0 Å is the right acceptance bar.** The work fixed the *admissibility of the
   pre-filter against* that bar; it did not examine the bar.
5. ⛔ **Protonated histidine is a cation and is not treated as one here.** Histidine enters the
   search as an imidazole only. A deliberate first-version choice, never validated.
6. ⛔ **Whether a length re-specification fixes the targeting failure.** The scaffolded redesign
   does raise the binder from a median 68 residues to 141–148, which is the direction the two
   independent routes in item 4 of the limitations list point. ⛔ But **no co-fold has been run on
   it**, so whether
   length was the cause remains untested.
7. ⛔⛔ **Whether a construct that constrains four of the six specified groups by identity, and the
   other two only by geometry, is an acceptable design.** It is a consequence of releasing the held
   element's sequence, it is measured, and it is **not** addressed as a decision anywhere in the
   material this directory was written from (`methods/core_scaffolded_design.md` §4).
8. ⛔ **Whether a larger helix is the mechanism, or merely correlated with it.** The screen shows the
   arms holding one 23-residue helix are the top four. It does not show *why*, and the obvious
   confound — those arms also hold **fewer, larger** segments overall — is not separated from the
   helix itself by any measurement here.
9. ⛔ **Nothing about the hardware requirement of the generative step is known.** The runs captured
   no memory figure, so whether these arms would fit a smaller accelerator is unmeasured.

---

## What is in this directory

| | |
|---|---|
| `README.md` | this file — the limitations, first |
| `methods/README.md` | reader-order guide to the method documents |
| `methods/structural_mining.md` | ⭐ the complete, self-contained account of the mining method: the constellation, the matching criterion, the tolerance algebra, both databases, the accessibility gate, the superposition and clash test |
| `methods/mining_negative_results.md` | the three measurements that were run, controlled, and returned nothing — the Cα search, the clash filter, the density map |
| `methods/core_scaffolded_design.md` | ⭐ the scaffolded redesign in full: the decomposition, the construct and its three decisions, the ten-arm screen, and the sequence-design step with its composition refit |
| `data/` | 16 result tables, 3 cohort lists, and the ten arms' complete input specifications and input coordinates, with `data/README.md` explaining every column and key |
| `scripts/` | the 24 programs behind every number here, with `scripts/README.md` saying what each does and the phase it belongs to |
| `structures/` | two viewer sessions: the four natural candidates on the constellation, and six generated designs from the screen |
| `figures/` | seven renders of the screen designs, with `figures/README.md` on the colouring and the three render defects that were fixed rather than shipped |

---

## A note on numbers in this directory

Every number here is marked by where it came from. **MEASURED** means a program in `scripts/` was
run on data in `data/` or on a database named in `methods/structural_mining.md`. **READ** means it
was taken from a published structure or paper. Nothing here is estimated, and no number was rounded
up from a smaller one.

The marker glyphs (⛔ ⭐ ⚠️ ✅) are emphasis, carried over from the working documents as in the
other challenge folders. They rank nothing and encode nothing.
