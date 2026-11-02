# Scaffolded design: keeping part of a real protein and generating the rest around it

This document covers the work that follows the structural-mining line in
`structural_mining.md`. Read that first — this one assumes its definitions (the constellation, the
six typed functional groups, the acceptance bar, the face gate, the superposition and clash test).

**The transition from mining to design in one sentence.** The mining line ended with four natural chains that present the
designed arrangement of functional groups on a reachable face, and all four put hundreds to
thousands of heavy atoms inside the receptor when actually placed; so none of them can be used
whole. The question this document answers is whether **part** of one of them can.

⛔ **Nothing here is a binding result.** No designed sequence in this document has been co-folded
against the receptor, scored for interface quality, or tested at the bench. What is measured here is
whether a generative backbone model, given a rigid fragment of a real protein and the receptor, can
build a chain that is **continuous, folds as one body, and stays out of the receptor**. That is a
prerequisite question, not the question.

---

## 1. The reframing, and the measurement that motivates it

The best-fitting candidate from the face gate is `5xox` chain B — 229 residues, numbered 5–233,
matching the six-group constellation at **1.5400 Å** against a 2.0 Å acceptance bar.

⭐ **MEASURED, and this is the whole argument for the reframing:**

| placement | heavy atoms of the candidate within 2.5 Å of the receptor crop |
|---|---|
| the **whole** chain, superimposed on the constellation | **679** atom pairs (closest approach 0.40 Å) |
| the whole chain with **both termini trimmed** | **679** — trimming removes **none** of them |
| **subsets of secondary-structure elements**, each re-placed | **0** for nine of the ten subsets screened. ⚠️ On the *first* enumeration of subsets (§4.1), where a closest-approach figure was also recorded, the best options reached zero clashing atoms at a closest approach of **2.54 Å** |

⇒ ⭐⭐ The collision is not in the tails. It is in the **middle** of the chain, so no amount of
trimming reaches it, and the only way to keep the useful part of this scaffold is to keep the
elements that do not collide and throw the rest away. A closest approach of 2.54 Å between two
chains is ordinary packing, not interpenetration.

⚠️ **Provenance note on the zero.** For the ten arms actually screened, the deposited sidecars record
the clashing-atom count (**0** for nine of ten) and **not** a closest-approach distance — that column
was only kept during the first enumeration. ⇒ the 2.54 Å figure describes element subsets of the same
scaffold on the same constellation, but not literally the ten arms in §5.

**What is kept, and what is not.** The elements keep 5XOX's **native relative geometry** — that
rigid arrangement is the entire reason to borrow from a real protein, because it supplies a span a
short de novo chain measurably cannot reach (`structural_mining.md` §6.1, "What the named list says: size and fit are opposed"). Everything
between the kept elements is thrown away and **generated from scratch with the receptor present**.

⚠️ **This is a different construct from an earlier attempt that failed**, and the difference is
worth stating because its result does not transfer. In that earlier attempt the core was delivered
as a **spatial obstacle** that the chain was asked to fold *around* and never to connect — it was
present in the box but was deliberately not part of the generated backbone path, and the contig
never named it. It produced **0 of 144** continuous designs. Here every kept element **is** named in
the contig, in sequence order, as an ordinary motif segment, which is what the generator is built
for. ⇒ the earlier 0-of-144 does not bear on this construct, and §5 below measures what this one
actually does.

---

## 2. The scaffold decomposed

Secondary structure was assigned on the chain itself (DSSP, simplified to contiguous H/E runs of at
least three residues), giving **16 elements — 11 helix, 5 strand**:

| element | type | residues | | element | type | residues |
|---|---|---|---|---|---|---|
| 1 | helix | 7–12 | | 9 | strand | 130–137 |
| 2 | strand | 23–30 | | 10 | helix | 140–169 |
| 3 | helix | 33–40 | | 11 | helix | 175–182 |
| 4 | helix | 48–64 | | 12 | helix | 187–195 |
| 5 | strand | 69–74 | | 13 | helix | 202–204 |
| 6 | strand | 78–83 | | 14 | helix | 207–210 |
| 7 | helix | 89–92 | | 15 | strand | 213–216 |
| 8 | **helix** | **94–116** | | 16 | helix | 229–232 |

**Two elements are excluded from the candidate pool, for measured reasons:**

- **Element 10** (helix 140–169, 30 residues) accounts for **337 of the 679** close atom pairs, closest
  approach 0.40 Å. It is also the **only** element that contains one of the six constellation
  residues other than the pair in element 3 — so excluding it caps how much of the constellation can
  be anchored in borrowed structure. It is kept in **exactly one arm, as a deliberate positive
  control** (§5.3).
- **Element 12** (helix 187–195) accounts for another **105** close atom pairs and contains no
  constellation residue at all. Excluded outright.

**Where the six constellation residues sit in this scaffold**, which determines what can be
anchored:

| constellation residue in the scaffold | inside an element? |
|---|---|
| 34, 38 | **yes** — both inside element 3 (helix 33–40) |
| 155 | **yes** — inside element 10, the excluded collider |
| 32, 211, 180 | **no** — loop residues, in no element at all |

⇒ ⛔ **At most 3 of the 6 can be anchored in real structure, and the third costs the collision.**
Every arm is therefore required to hold element 3, and the other constellation positions enter as
single-residue guideposts at their specified positions.

**Why a bare guidepost is sound here.** Each of the six residues, placed **alone**, is essentially
collision-free against the receptor crop: 0 atoms within 2.5 Å for four of them, 3 atoms for one and
1 atom for another. ⇒ what collides is never the functional group; it is the helix wrapped around
it.

---

## 3. Placement is re-solved for every option, and the degree-of-freedom count matters

The kept fragment's **position** is not inherited from the whole-chain superposition. It is set by
fitting the constellation residues that the fragment **contains** onto their specified positions, so
it changes with the element set:

| contained constellation residues | placement | why |
|---|---|---|
| **≥ 3** | Kabsch superposition | the rigid body is fully determined |
| **= 2** | Kabsch **plus a one-parameter sweep** | two points do not fix a rigid body: one rotation about the axis through the two centroids survives. It is swept in 5° steps and the angle that **minimises collision with the receptor** is kept |
| **≤ 1** | rejected | the placement would be arbitrary |

⛔⛔ **The resulting fit is over 2 or 3 points, not 6, and must never be quoted as "the scaffold fits
to 0.9 Å".** The reported figures — **0.9187 Å** for the two-point arms and **0.8644 Å** for the
three-point control — are **not comparable** to the 1.5400 Å six-group figure, and a two-point fit
being numerically smaller than a six-point fit is arithmetic, not an improvement. The bare
guideposts sit at their specified positions **exactly, by construction**, so the designed site is
still exactly specified; what is unproven is whether the generator can build linkers that reach
them.

⚠️ **A fitting bug was found here and it is recorded because of how it was caught, not because of
what it cost.** A first version of the option designer hard-coded the pin order as the residue
numbers **sorted numerically**, with a comment asserting that this was the query order. It is not:
the query order is the order of the hit list's `resid_1..6` columns. Three of the six therefore
fitted onto the wrong target positions. ⛔ The symptom that exposed it was an **impossibility**: a
four-point least-squares fit was reported as *worse* than the six-point fit of the same chain, which
cannot happen. Two of the three mismatches were also **residue-type** mismatches — an imidazole
fitted onto a cation position and a cation onto an imidazole position — ⇒ ⭐⭐ the residue-type
assertion that `structural_mining.md` §8 describes as standard in this line, and which this one
program omitted, would have caught it on the first subset. The order is now **read from the hit
list** and never typed in, and the type assertion runs before any enumeration (6/6 agreement).

---

## 4. The construct, and three decisions that are choices rather than defaults

**D1 — the kept elements are indexed motif segments at full backbone**, named in the contig in
sequence order. ⛔ Not the off-path obstacle channel, and ⛔ not a reduced Cα/Cβ representation.
Those two changes were bundled together in the earlier failed attempt and this construct unbundles
them.

**D2 — the kept elements' SEQUENCE is released.** The backbone of the chosen strands and helices is
held; their amino-acid identities are not, so the sequence-design step still designs the scaffold
surface and a natural protein's sequence is not imported wholesale. ⚠️ The alternative — holding
5XOX's own sequence to buy foldability — was **not** taken.

**D3 — a constellation residue inside a kept element takes the scaffold's own residue identity** in
the generator's input file, not the design's. The scaffold supplies lysine at four of these
positions and histidine at two, where the design asked for arginine at the cation positions. This is
defensible for the match itself, because the search that selected this scaffold treated an arginine
guanidinium and a lysine `NZ` as the **same** cation type. ⚠️ It is nonetheless a real chemical
change: a guanidinium is planar and bidentate, a lysine `NZ` is a point charge.

⛔⛔ **A consequence of D2 and D3 together, which the deposited inputs make plain and which should
be read as a limitation of this construct.** D2 releases the sequence of the whole kept element, and
the two constellation residues that sit inside element 3 are inside that element. So their identity
is held in the generator's **input file** and then released downstream. MEASURED in the deposited
specifications (`../data/core_arm_specs_20261101/`): each arm names **four** positions in
`select_fixed_atoms` — three cation groups and one imidazole — and the sequence-released selection
covers element 3 in full. MEASURED in the sequence-design output
(`../data/j904_pin_positions_20261101.tsv`): **four** fixed positions per design, identities
**R, H, R, R**, in all 312 sequences.

⇒ ⛔ **of the six functional groups the mining line specified, four are identity-constrained end to
end and two are not.** The two that are not keep their **geometry** — they are part of a held
backbone segment — but nothing downstream requires them to be a cation and an imidazole. Whether
that is acceptable is a design question this work does not settle, and it is **not** stated as
settled anywhere in the material this document was written from.

### 4.1 Linker sizing, and a re-selection that changed which options were screened

Linkers are sized from the **placed** geometry as a range, never as an exact length: the minimum is
`ceil(Cα–Cα distance / rise)` and the arm declares `[min, min + slack]`.

⛔ **Native loop lengths are not the budget.** MEASURED: 5XOX's own path between two of the
constellation residues takes **178 residues to cover 24.2 Å**, so the native route is a gross
detour and sizing on it would inflate every arm.

⚠️⚠️ **The rise constant was changed, and it changed which options are in the screen.** The first
option set was sized on an **extended-chain** bound of ~3.5 Å per residue. Designs from this pipeline
come out **all-helical**, and a helix advances only ~1.5 Å per residue, so every linker minimum
roughly doubles. Re-sized helically — with each gap asserted to be **at or above** the extended
floor, so the bound is never loosened — the first option set totals **173–229** residues, outside the
120–160 band the construct was specified for. A **second selection** was therefore made under
helical sizing, holding **fewer and smaller** elements, and **that** is the set that was built,
tested and screened.

⇒ ⛔ **the ten arms in §5 are not the ten options of the first selection.** They hold different
element sets and carry different collision counts. ⚠️ The deposited sidecars
(`../data/core_arm_specs_20261101/*_provenance.json`) are the authoritative record of what was
actually screened, and every per-arm number in §5 is checkable against them. The earlier set is
**not** deposited here and its table is not reproduced, because it was superseded before anything
ran. ⛔ **What is unresolved about the re-selection:** no written record of it exists — it is
recoverable only by comparing the built arms against each other and reading the sizing constant in
the builder's header — so the reason the two sets differ is reconstructed here rather than cited.

⛔ **One contig-syntax trap is load-bearing and is documented in the builder:** a chain-break token
must carry its own comma. Written without one it registers **no break at all** and silently changes
the output layout — a failure with no error message.

---

## 5. The ten-arm screen

**What ran:** ten arms, one per core option, **40 designs each = 400**, the receptor present in every
case, generated as 5 independent batches of 8 models per arm.

⭐ **The resolved configuration was recorded rather than assumed**, and it is **identical across all
ten arms and the single-arm test run**: `step_scale = 3`, `gamma_0 = 0.2`, `n_recycle = 3`,
classifier-free guidance off. ⇒ nothing drifted between the test and the screen, which is the only
reason the two can be compared at all.

**Rate, from the single-arm test run:** **19.12 s per design**, i.e. **153 s per batch of 8**, for a
149-residue binder against the 154-residue receptor crop, on an 80 GB data-centre GPU
on a single data-center accelerator. The screen itself ran at 18.15–18.85 s per design.

⛔⛔ **No memory figure exists for any of these runs, and this is a gap rather than a reassurance.**
The runner never samples the device, so it is **unmeasured** whether these arms would fit a smaller
card. ⇒ no hardware recommendation below the card actually used can be made from this work, and the
rate above should be read as *"this is what it cost on this card"*, not as a sizing.

⚠️ **An earlier rate from this project does not compare** to the 153 s figure: it was a 55–82
residue binder on a different card, so length **and** hardware both differ. Two variables, one
comparison; quote them separately or not at all.

### 5.1 ⛔ Continuity was gated first, and no binary cut was applied

⛔⛔ **Backbone continuity is measured and gated before any clash or interface number**, because a
minimum-distance-to-target floor **rewards** a chain that has broken apart: the pieces can retreat
from the target in a way an intact chain cannot. A clash column read first would have ranked the
breaks at the top.

The gate re-measured its own bar ladder on this cohort and **refused** a binary continuity cut:

| no Cα–Cα step above | designs, of 400 |
|---|---|
| 4.0 Å | 303 (75.8%) |
| 4.1 Å | 346 (86.5%) |
| 4.5 Å | 350 (87.5%) |
| 5.0 Å | 351 (87.8%) |
| 7.0 Å | 353 (88.2%) |
| 10.0 Å | 371 (92.8%) |

The rule is that a cut is admissible only where the counts stop moving. ⛔ Here 351 / 353 / 371 at
5 / 7 / 10 Å **differ**, so designs genuinely sit *between* "ideal, about 3.9 Å" and "broken". ⇒ **no
verdict column exists in the deposited table** and 4.1 Å is used below as a **reporting** bar only,
never as a pass mark.

⭐ **Context worth keeping:** the earlier off-path core-transplant cohort had a Cα–Cα step over 10 Å
in **69%** of designs. Here it is **7.2%** (29 of 400). ⇒ on exactly the failure the reframing was
meant to fix, the motif-scaffolded construct is a large improvement.

### 5.2 ⭐⭐ The result: the arms separate, and the split tracks one element

The ranking below is a **conjunction** — a design must pass every criterion — and ⛔ never a mean of
fractions, because an averaged label pays partial credit for failing one test outright.

| arm | elements held | `max_step ≤ 4.1` | `bb_min ≥ 2.0` | both | + `min_cut ≥ 15` | median length |
|---|---|---:|---:|---:|---:|---:|
| **opt05** | 2, 3, **8**, 15 | 39/40 | 16/40 | 16/40 | ⭐⭐ **16/40** | 145 |
| **opt01** | 3, 4, **8**, 15 | **40/40** | 12/40 | 12/40 | ⭐ **12/40** | 145 |
| **opt04** | 3, 6, 7, **8**, 15 | 38/40 | 11/40 | 11/40 | **11/40** | 142 |
| opt09 | 3, **8**, 13 | **40/40** | 8/40 | 8/40 | 6/40 | 129 |
| opt07 | 3, 4, 9, 15 | 37/40 | 6/40 | 5/40 | 5/40 | 130 |
| opt02 | 2, 3, 4, 5, 15, 16 | 22/40 | 4/40 | 3/40 | 3/40 | 145 |
| opt06 | 2, 3, 4, 6, 15 | 37/40 | 3/40 | 3/40 | 3/40 | 142 |
| opt08 | 2, 3, 5, 6, 15, 16 | **15/40** | 3/40 | 1/40 | 1/40 | 147 |
| opt03 | 3, 4, 5, 13, 14, 15, 16 | 38/40 | 10/40 | 10/40 | ⛔ **0/40** | 145 |
| opt10 *(control)* | 3, **10**, 13 | **40/40** | ⛔ **0/40** | 0/40 | 0/40 | 134 |

`max_step` is the largest Cα–Cα step in the chain (continuity); `bb_min` is the smallest backbone
distance to the receptor; `min_cut` is the minimum graph cut of the long-range contact map — the
"is this chain **one** fold" number.

⭐⭐ **Element 8 is the large helix, residues 94–116, 23 residues.** Sorting the arms by whether they
hold it:

- **holds element 8:** opt05 **16**, opt01 **12**, opt04 **11**, opt09 **6**
- **does not hold it:** opt07 5, opt02 3, opt06 3, opt08 1, opt03 0, opt10 0

⇒ ⭐⭐ **the worst arm that holds the large helix (6/40) beats the best arm that does not (5/40).**
The four arms holding it are *exactly* the top four, with **no overlap** between the two groups.

⛔⛔ **How strongly to believe that, stated plainly.** This is **4 arms against 6**, and **n is
batches, not designs**: one contig length is drawn per batch and the 8 models in a batch share that
draw, so each arm is **5** independent observations, not 40. The per-batch breakdown is deposited
(`../data/j900_per_arm_20261101.tsv`) for exactly this reason. ⭐ The continuity signal is
nonetheless **batch-consistent** rather than one lucky draw: opt01, opt09 and opt10 are 8/8 in every
one of their five batches, while opt08 runs 4/8, 2/8, **0/8**, 3/8, 6/8. ⇒ treat *"hold one large
helix"* as **strongly suggested and worth designing around**, ⛔ not as established by one screen.

### 5.3 Two arms fail for two different reasons, and the control works

⛔⛔ **opt03 is the arm a continuity-only read would have passed.** It scores 38/40 on continuity and
10/40 on the raw clash bar — mid-table — but its median `min_cut` is **4.0**, and adding the
single-fold criterion takes it to **0/40**. A low minimum cut means **no crossing long-range
contacts**: two lobes on a tether rather than one fold. ⇒ ⭐ it is a direct demonstration that
**continuity and fold topology are orthogonal**. opt03 holds **seven** scattered elements, the most
of any arm, and the chain threads between them instead of folding as one body.

⛔ **opt08 fails the other way:** **15/40** continuity, the worst in the cohort, one batch at 0/8. It
holds **three separate strands** plus two further elements — the most fragmented strand set — and the
chain cannot connect them.

⇒ ⭐⭐ **stated as a direction and not a rule: fewer, larger elements connect; many small scattered
ones either break the chain (opt08) or split the fold (opt03).**

✅ **opt10, the element-10 control, scores 0/40, which is the control working.** Its continuity is
perfect (40/40, 8/8 in all five batches) and its median backbone distance to the receptor is
**0.379 Å** — it interpenetrates. That was **built in**: opt10 is the only arm whose **input** carries
a nonzero collision, **379 heavy atoms** against the receptor crop, recorded in its deposited sidecar
**before** anything ran. ⇒ it is a **positive control that the clash column is reading what it claims
to read**, not a candidate. Element 10 is confirmed as the part to drop.

### 5.4 ⛔ What the screen does not say

- ⛔ **It is not a binding result and not an interface result.** Nothing was co-folded, and no design
  was tested for whether it engages the intended epitope. The screen answers **one** question: which
  core option gives the generator a connectable, single-fold, non-interpenetrating backbone.
- ⛔ **The clash column is RAW generator output, which clashes normally.** `bb_min ≥ 2.0` is the
  documented **raw** floor, not a post-relaxation gate, and it is used **only in conjunction with
  continuity**, never alone.
- ⛔ **Guidepost reach is reported and not ranked on.** The deposited table carries a Cβ-based
  distance from each fixed position to the receptor; clearance must be gated on the pins' **heavy
  atoms**, and a Cβ-based gate has previously passed a 1.67 Å clash elsewhere in this project. So the
  column is printed and not used.
- ⛔ **A "did the held part land?" number is not a result.** Held coordinates read as in place by
  construction, so an element-placement RMSD measures constraint compliance, not design quality.
  It is used below only as a **consistency check on the construct**, never as a success metric.

---

## 6. What the pictures show, and how the held/built split was measured

The deposited session file and renders (`../structures/README.md`, `../figures/README.md`) colour
each design by three **disjoint** residue sets. Those sets could not be read off the output file or
off the contig, and the reason is worth recording:

- the output coordinate file flags only the **four bare guideposts** as sequence-fixed motif atoms,
  because the kept elements' sequence was deliberately released (D2) — so the kept elements are not
  flagged in the file at all;
- the contig cannot be used either, because linker lengths are **ranges drawn per batch**, so a kept
  element's output index differs from design to design.

⇒ the kept residues are identified **geometrically**, which re-verifies the construct as a side
effect. The origin token was measured to be consumed as a **pure translation** at exactly 0.0000 Å,
so a held residue's output Cα must sit at its input Cα minus that translation. ⭐ **MEASURED:** every
design matched its full expected count — 43 of 43, 52 of 52 and 45 of 45 for the three arms shown —
with a worst deviation of **0.0623 Å**. ⇒ the motif really is held, not merely requested.

⭐⭐ **The split is the point:** the generator **built 85 to 98 residues of each 141–145 residue
binder**, which is most of the chain. That grey region is the answer to the question the screen was
launched to settle.

---

## 7. Sequence design on the three leading arms

The cohort is the **three-way conjunction survivors** — designs that passed continuity **and** the
raw clash floor **and** the single-fold criterion, never any one of them alone:

| arm | backbones | batches represented | sequences |
|---|---|---|---|
| opt05 | 16 | 5 of 5 | 128 |
| opt01 | 12 | 3 of 5 | 96 |
| opt04 | 11 | 4 of 5 | 88 |
| **total** | **39** | — | **312** |

Sequence design is SolubleMPNN at temperature 0.1, 8 sequences per backbone, designing the binder
chain only.

### 7.1 ⭐ The guards, every one measured on the output

| | |
|---|---|
| ⛔⛔ cysteines in the designed chain | **0 of 312**. A free cysteine breaks expression; it is suppressed by an explicit omission |
| ⭐ fixed positions preserved | **all four held their identity in all 312** — identities `R, H, R, R` |
| sequence diversity | **312 distinct of 312**, no duplicates |
| binder length | **141–148** residues, 0 out of range |
| non-standard residues | **0** |
| target chains | untouched; the target keeps its own two cysteines, ⚠️ which is why a naive search for cysteine across a whole multi-chain record reads as a failure when nothing is wrong |

⭐ **The fixed positions were read per design from the specification and cross-checked against the
model's own motif annotation** — ⛔ never hard-coded, because each batch samples its own contig
length, so the fixed positions sit at different indices in every batch.

### 7.2 ⚠️ The finding: sequence-design confidence does not separate the arms

| arm | n | confidence (median) | range |
|---|---:|---:|---|
| opt05 | 128 | **0.4095** | 0.3619–0.4539 |
| opt01 | 96 | **0.4192** | 0.3677–0.4496 |
| opt04 | 88 | **0.4133** | 0.3876–0.4523 |

⇒ ⚠️ the three medians sit inside **0.010** of each other and their ranges overlap almost completely.
⛔ **Sequence-design confidence therefore carries no usable signal for choosing between the arms, and
must not be quoted as if it confirmed the backbone ranking.** The arm ordering rests on the backbone
gate in §5 and on nothing measured here.

⭐ That is a reasonable outcome rather than a failure: this score says how **designable a backbone
is**, which is a different question from which core the chain was built around. ⛔ And a high value
is **not** evidence of binding.

*(Medians over an even-sized set are reported as the upper of the two central values; recomputing
them as the mean of the two central values shifts each by ≤ 0.0003 and changes nothing.)*

### 7.3 ⭐⭐ The composition bias had to be refitted, and the settled values would have been wrong

An earlier, unrelated cohort in this project established a composition bias for this sequence-design
step — a per-amino-acid adjustment that pulls the designed composition toward natural frequencies.
That earlier cohort's pathology was **alanine**, and its fitted bias is strongest on alanine.

⛔ **The first run on this cohort applied no bias at all.** That was an oversight, not a decision,
and it is corrected below. ⭐⭐ **But inheriting the settled values would have been a second mistake**,
and the standing rule from that earlier work anticipates exactly this: *measure this cohort rather
than assume the last one's pathology.*

MEASURED over **43,696** free residues in 312 sequences, with the fixed positions excluded:

| residue | **this cohort** | the earlier cohort the bias was fitted on | natural |
|---|---:|---:|---:|
| **G** | ⛔ **13.18%** | not a flagged excess | 7.07 |
| **A** | **11.24%** | ⛔ **21.03%** | 8.25 |
| **L** | **8.72%** ⇒ *below natural* | 10.68% | 9.66 |
| **S** | **4.81%** ⇒ *below natural* | — | 6.56 |
| E | 12.56% | 16.47% | 6.75 |

⇒ ⭐⭐ **the settled bias is aimed at the wrong residues for these backbones.** It would push alanine
far past natural; it suppresses leucine and serine, which are already **below** natural here, i.e.
the wrong direction; and its glycine term is the **smallest** of the four even though **glycine is
this cohort's dominant excess**.

⭐ **Why the two cohorts differ, and it is the construct.** The earlier binders were **52–79 residue
de novo helical bundles**. These are **141–148 residues built around a held fragment of a real
protein**. The unbiased composition here is already **closer** to natural than the earlier cohort
reached after nine rounds of refitting. ⇒ the alanine pathology is largely absent, and quoting it as
a property of the pipeline would have been the error the standing rule warns against.

**The refit.** Settings held from the settled configuration: the whole charged set
(D, E, K, R) **exempt**, negative-only adjustments, temperature 0.1. ⛔ The bias is applied
**pre-temperature**, so at temperature 0.1 every value is multiplied by 10 and the per-round
corrections **compose additively**. Each round's bias was read back from its own run log rather than
reconstructed.

| round | cumulative pre-temperature bias | G % | A % |
|---|---|---:|---:|
| 0 | *(none)* | 13.18 | 11.24 |
| 1 | `G −0.062, A −0.031, V −0.012, P −0.019, T −0.002` | 13.08 | 10.89 |
| 2 | `G −0.124, A −0.059, V −0.023, P −0.038, T −0.005` | 12.93 | 10.47 |
| 3 | `G −0.185, A −0.084, V −0.034, P −0.057, T −0.008` | 12.83 | 10.15 |
| 4 | `G −0.246, A −0.106, V −0.045, P −0.076, T −0.012` | 12.72 | 9.88 |
| **5 ⭐ adopted** | **`G −0.307, A −0.126, V −0.055, P −0.094, T −0.016`** | **12.58** | **9.53** |

⛔ **Stopped at round 5, on the criterion fixed in the earlier work**: returns had fallen to about
0.4 composition points per round, and alanine at 9.53% is already past the practical floor measured
on a different target in this project. ⛔ **The adopted cohort is the round-5 measured output, never
the round-6 proposal** — that proposal has never been run and its composition is unknown.

⚠️ **The aggregate deviation score is deliberately not quoted as a single headline number here**,
because it is defined over the *designable* residue set — the summed absolute difference from natural
over the residues the run was not told to omit — and it therefore changes when the omission list
changes. ⛔ Re-measured over **the same deposited sequences** under the two omission lists used in
this work, the unbiased cohort scores **32.3** under one and **30.9** under the other, so a
before-and-after pairing that crosses the two is partly an artefact of the definition. ⇒ **the
per-residue percentages above are the comparable quantities**, and they were re-measured directly
from the deposited sequence tables. Any aggregate figure should be quoted only with its support
stated.

### 7.4 ⭐⭐ The finding: the glycine excess is structural, not a sampling artefact

MEASURED: **glycine moved 13.18% → 12.58% — only 0.60 points — for a cumulative −0.307 bias, which
at temperature 0.1 is an effective shift of −3.07 in the logit.** Over the same five rounds alanine
moved **1.71** points on a shift less than half as large.

⇒ ⭐ **these backbones genuinely demand glycine, and the bias cannot argue them out of it.** That is
the same mechanism the earlier work names for alanine — the backbone has already specified what goes
there, so it is the wrong target for a composition knob — but acting here on a **different residue**.

⚠️ **The leading explanation, stated as an explanation and not as a measured cause:** the chain is
mostly generated connector (85–98 of 141–148 residues, §6) threading between held core elements, and
turns need glycine. ⛔ The test would be glycine content against loop fraction per design, and **it
was not run**.

### 7.5 ⛔⛔ Omitting histidine protected the *fixed* imidazole and destroyed the *contained* one — the construct defect this cohort exposed

⛔⛔ **This section replaces an earlier reading that was too narrow, and the correction is the
finding.** The earlier text reported that omitting histidine was safe because an explicit
fixed-residue list takes precedence over the omission list. That precedence rule is real, and it
does protect a **fixed** histidine. ⛔ **It says nothing about a constellation histidine that was
never fixed** — and one of the six was not.

**The construct.** Four of the six specified functional groups enter as single-residue guideposts
with their identity fixed (`n_pins = 4`, `pin_identities = RHRR`). The other two — scaffold
positions **34** (an imidazole was specified) and **38** (a cation was specified) — sit inside
held element 3, the helix at 33–40, whose **sequence was deliberately released** so the surface
could be redesigned. ⇒ the design run was free to change them, and did. The quality-control
column `pin_identity_held` cannot see this: it checks the four guideposts, so it reads true for
every sequence in **every** cohort while two of six groups are missing.

**MEASURED**, per sequence, over all three deposited cohorts — the two contained positions located
from each design's **own** held-element sequence rather than hard-coded, because every backbone
samples its own connector lengths (`../scripts/j953_contained_pin_check_20261102.py`, deposited
row-by-row as `../data/j953_contained_pin_check_20261102.tsv`):

| cohort | scaffold 34 keeps an imidazole | scaffold 38 keeps a cation |
|---|---:|---:|
| without histidine omission (`j903`, superseded) | **43 / 312 (13.8%)** | 146 / 312 (46.8%) |
| **with histidine omission (`j904`)** | ⛔⛔ **0 / 312 (0.0%)** | 159 / 312 (51.0%) |
| **six positions locked (`j908`)** | ✅ **312 / 312 (100%)** | ✅ **312 / 312 (100%)** |

⛔⛔ **The zero is the point.** Omitting histidine from the designable alphabet did not merely fail
to protect position 34 — it removed the only residue that could satisfy it, and **190 of 312
(60.9%)** of those positions came back as aspartate or glutamate: **an acid where an imidazole was
specified.** Without the omission the same position was already wrong 86.2% of the time, so
⛔ **the omission list was never the real defect** — it converted a frequent failure into a
certain one.

⭐⭐ **The arithmetic was self-refuting before anything was re-measured.** The specified
constellation contains **two** histidines, so an intact design must carry `n_his >= 2`. The
`j904` cohort reports `n_his = 1` in **all 312** sequences. A cohort cannot hold two imidazoles and
contain one histidine. ⇒ the same table that was read as evidence of success was proof of the
loss, and the "stray histidine" counter had been **counting a specified constellation position as a
stray** — which is why the unbiased cohort's 43 correct sequences were scored as defective.

**The fix, and it is deposited.** Both contained positions are given explicit fixed identities in
sequence design, exactly like the four guideposts, mapped through each design's own index map rather
than by position. All six groups then survive into every sequence (`j908`, the table above, and
`n_his = 2` in 312 of 312). ⚠️ **One honest consequence:** the contained positions are held at the
**scaffold's own** identities, which gives **lysine** at 38 where the specification said arginine.
Both are cations and the constellation was matched on functional-group type, so the type is
satisfied — but it is a real chemistry change and is recorded rather than absorbed.

⚠️ Whether stray histidines *should* be suppressed at all is a judgement about the intended pH
behaviour of the binder and is **not** settled here.

### 7.6 Quality control on the adopted cohort

| hard check (a failure is a defect) | |
|---|---|
| fixed-position identity lost | ✅ **0** |
| cysteine present | ✅ **0** |
| non-standard residue | ✅ **0** |
| length out of range | ✅ **0** (141–148) |

| reported | value |
|---|---|
| design confidence | median **0.4113**, minimum 0.3635 ⇒ **0 of 312** below the 0.15 health floor |
| ⭐ confidence cost of the bias | **−0.0007**, paired per backbone (n = 39) — effectively free; the earlier cohort's bias cost −0.0122 |
| **isoelectric point** | median **5.03**; **pI ≥ 7 in 4.81%** (15 of 312) |
| net charge per residue at pH 7.4 | median −0.0376, acidic |
| N-X-[S/T] sequons | 40.06% carry at least one — ⛔ **inert here**: the expression system has no glycosylation machinery. Counted for the record only |

⭐⭐ **The isoelectric-point figure is the one that matters, because it is the defect the earlier
cohort had to be refitted for.** That run's first attempt exempted only the **basic** residues from
the bias, which left the acidic counterweight to the fixed arginines exposed and drove the fraction
with pI ≥ 7 from 9.86% to **56.26%** — into a band where expression is known to fail. Exempting the
**whole** charged set here keeps it at **4.81%**. ⇒ the corrected convention transferred and the
defect did not recur.

⚠️ **A carried limit on the bias, from the work that settled it:** on another target the same bias
was **null for design quality** (paired test, p = 0.657). ⇒ ⛔ it earns its place on **composition
plausibility only**. Do not expect it to buy foldability and do not justify it by diversity. The
same holds here: confidence moved −0.0007, i.e. nothing.

### 7.7 ⚠️ The superseded first cohort is kept on purpose

The unbiased first run is deposited alongside the adopted one
(`../data/j903_mpnn_sequences_unbiased_20261101.tsv`). ⛔ **It is superseded for any downstream use**
and must not be ranked or folded. ⭐ It is kept because it is the **round 0** the refit was fitted
against and the **baseline** the paired confidence comparison is computed against — and because
every guard result in §7.1, and the finding in §7.2 that confidence does not separate the arms, were
measured on it first and still hold on the biased cohort.

⚠️ **One check in the sequence-design driver is vacuous on this cohort and is labelled rather than
counted as a pass.** The driver reconciles a side-chain-repaired coordinate file against the raw
file, because an earlier cohort's files lost their motif annotation during repair. This cohort needed
**no repair at all** — every non-glycine residue already carried side-chain atoms — so the same file
was passed as both and the comparison compared a file with itself. ⛔ **An inapplicable check is not
a passing check.**

---

## 8. ⛔ What this whole line is NOT, and what the next real test is

1. ⛔⛔ **Not a binding result.** 312 sequences exist. **None has been folded against the receptor or
   scored.** Every number in §5 and §7 is about backbone geometry or amino-acid composition.
2. ⛔ **The next real test is a co-fold of each sequence against the receptor, asking whether the
   designed contacts are reproduced** — which is the test the previous design round failed, with
   2,135 of 2,160 predictions reproducing **zero** of four designed contacts
   (`../README.md`, item 2 of its limitations list). Nothing in this document predicts that outcome.
3. ⛔ **The arm choice is measured; the depth is not.** §5 says which core option to take further. It
   does not say that any of them produces a binder.
4. ⚠️ **Two of the six specified functional groups are not identity-constrained downstream** (§4) and
   that is a property of the construct, not an oversight in a single run.
5. ⚠️ **No expression, solubility, oligomeric-state or functional measurement of any kind** appears
   anywhere in this directory.

---

## A note on numbers in this document

**MEASURED** means a program in `../scripts/` was run and produced the value, on data deposited in
`../data/`. Per-arm figures were re-derived from the deposited 400-row table during preparation of
this document and agree with the source record cell for cell, including the bar ladder
(303 / 346 / 350 / 351 / 353 / 371), the conjunction counts (16 / 12 / 11 / 6 / 5 / 3 / 3 / 1 / 0 / 0)
and the 29-of-400 over-10 Å count. Composition percentages in §7.3 were likewise re-measured from
the deposited sequence tables.

⛔ **Where a source disagreed with itself or with a re-measurement, nothing was silently chosen.**
Three such cases arose and all three are stated where they occur rather than smoothed over: the
composition-deviation score and its support (§7.3), the confidence medians and the median convention
(§7.2), and the two option sets — the one enumerated and the one screened (§4.1).
