# Three measurements that were run properly and returned nothing

These are the controlled negatives from the structural-mining line. Each was built to answer a real
question, each was run with a null or a control, and each came back empty. They are kept because a
negative with a working control is a result, and because two of them retired an instrument that
would otherwise still be in use.

Read `structural_mining.md` first — it is the method these three sit inside.

---

## 1. ⛔⛔ A Cα-only constellation search is not selective

### The question

Before the functional-group representation existed (see `structural_mining.md` §2.2), the first
search matched **Cα positions only**: find five Cα atoms anywhere in any chain whose internal
geometry matches the designed five-point constellation within an RMSD cutoff, ignoring residue
identity entirely. The query at that stage was a five-point set with a ~25 Å span and two tight
pairs (5.47 Å and 4.57 Å Cα–Cα).

### The control, which is what makes this a result

⭐ A raw hit count means nothing on its own. Three **control constellations** were drawn from real
chains, each with a sorted distance vector matching the query's to within 2.0 Å. If the designed
constellation is in any way special, it should match **fewer** chains than a generic one of the same
shape.

### The store this ran on

A local **foldseek PDB100** database, snapshot dated 2025-01-01: **324,204** searchable chains
covering **151,880** distinct PDB entries, clustered at 100% sequence identity so redundant chains
are collapsed. ⛔ **It holds Cα coordinates, sequence and secondary structure only — no side chains,
no ligands and no cofactors**, which is exactly why a Cα search was the first thing attempted: it
needed no download. Extracting the whole chain set took about 4 minutes and parsing it to an array
store 50 s, giving 72,861,389 Cα atoms; chain lengths ran from 4 to 4,579, median 170.

⚠️ **A coverage count that was wrong first time**, kept because the shape of the error is common:
taking distinct identifiers from the database's own lookup table gave **224,856**, but only 324,204
of those keys exist in the searchable index, and joining the index to the names gives **151,880**
distinct entries. ⇒ the lookup lists chains the database does not hold. **A key count is not a
coverage count.**

⛔ Note also what the Cα-only content *cannot* support: a match fixes the backbone but leaves the
functional group free to rotate, so nothing measured in this section speaks to chemistry.

### The result

**MEASURED**, query versus the three controls over the first 2,000 chains of that store —
hits / distinct chains:

| RMSD cutoff (Å) | designed query | control 1 | control 2 | control 3 |
|---|---|---|---|---|
| 2.50 | 359,853 / 1,466 | 122,774 / 1,446 | 274,027 / 1,490 | 250,258 / 1,495 |
| 1.50 | 200,848 / 1,442 | 70,069 / 1,402 | 182,610 / 1,469 | 228,380 / 1,492 |
| 1.00 | 96,189 / 1,388 | 39,995 / 1,356 | 96,077 / 1,430 | 131,378 / 1,455 |
| 0.50 | 4,282 / 982 | 3,880 / 914 | 7,274 / 1,160 | 10,081 / 1,206 |
| 0.35 | 172 / 142 | 259 / 210 | 489 / 330 | 923 / 510 |
| 0.25 | 5 / **5** | 17 / 17 | 37 / 32 | 73 / 64 |

⇒ ⛔⛔ **At the loosest cutoff the designed constellation is found in 73% of all protein chains, and
at every cutoff it sits at or below the controls.** Best RMSD achieved: query 0.206 Å, controls
0.201 / 0.094 / 0.132 Å.

⇒ ⭐⭐ **A Cα-only, identity-free five-point constellation search has essentially no selectivity, and
tightening the RMSD cutoff does not rescue it** — the ordering never inverts in the query's favour.

### Why, and why the standard fix does not apply

The combinatorics force it. A 170-residue chain (the median chain length in that store) holds
C(170,5) ≈ **1.1 × 10⁹** five-subsets, and ten pairwise distance constraints at ±1 Å cannot cut that
to a small number.

The standard cure is to match **segments** rather than isolated points — requiring i−2…i+2 around
each matched position turns 5 points into 25 and buys the selectivity back. That is what a
substructure-search tool such as **MASTER** does, and it is normally the right instrument for a
query of this kind.

⚠️ **And a fold-level search is the wrong instrument for a different reason.** A whole-chain
structure search (foldseek was available locally and was smoke-tested here) matches a contiguous
chain by its backbone and structural-alphabet profile and answers *"what fold is this like"*. It has
no way to express a **disjoint** positional constraint, and a fold-level hit says nothing about what
sits away from the motif. ⭐ It is a useful instrument for *characterising* a scaffold once found,
not for finding one. For the record, that smoke test also confirmed what was expected of a de novo
bundle: a 65-residue designed binder chain searched against all 324,204 chains returned 11 hits, the
best at 21.8% identity and an E-value of 1.6 — i.e. no real natural homologue.

⛔⛔ **It cannot be used here, because this query has no backbone context to match.** The guideposts
are **isolated single residues** with no neighbours; there is no i±2 to constrain. ⇒ installing such
a tool would not have helped this query. It becomes worth it the moment a query has real backbone
segments.

⭐ **The useful part of this negative is the diagnosis**, which is what led to the functional-group
representation: if a search is unselective *and* its query has no backbone, the backbone was never
the thing being specified.

### ⚠️ An engineering note worth keeping

**MEASURED:** the first implementation ordered the search "tightest pair first" and ran at over
300 ms per chain. That is backwards. In a dense protein the **short**-distance shell is the most
heavily populated — every residue has near neighbours — so anchoring on the tightest pair
**maximises** the branching factor. Picking the **least-populated shell in that particular chain**
instead, and growing assignments in batched array operations, gave **21.9 ms per chain, a 14×
speed-up.** The same principle is used in the functional-group matcher.

---

## 2. ⛔⛔ The target-clash filter has no discriminating power and is retired

### The question

Each match can be superimposed into the target's own coordinate frame. It then seemed natural to ask
whether the matched chain **clashes with the target** there, and to use that as a filter: keep the
matches that do not.

### ⛔ A control that does not apply, stated before the result

⛔⛔ **The control constellations from §1 cannot be compared on this metric.** A control is a
constellation lifted from an unrelated chain; superimposing a candidate onto *it* lands the candidate
in a frame with **no defined relationship to the target**, so nothing can clash by construction. The
run duly reported the controls 100% clash-free and the designed query 0% — ⛔ **that contrast is an
artefact of the frame, not a finding, and must never be quoted as one.** A valid control for a clash
filter has to be defined **at the real site**.

### The two nulls that killed it

**Null 1 — random orientation.** Take the same matched chains, centred at the same place, and give
each **1,000 random orientations**:

| | clash min | median | mean | clash-free |
|---|---|---|---|---|
| actual (constellation-matched) | 52 | **129** | 132 | **0.0%** |
| null (random orientation) | 36 | ~116 | ~119 | **0.0%** |

**Null 2 — binder scale.** Restrict to atoms within radius *R* of the constellation centroid, i.e.
the volume a 60–90-residue binder could occupy:

| *R* (Å) | clash, matched | clash, random | clash-free (both) |
|---|---|---|---|
| 12 | 31 | 29 | 0.0% / 0.0% |
| 15 | 49 | 50 | 0.0% / 0.0% |
| 20 | 86 | 87 | 0.0% / 0.0% |
| 25 | 113 | 118 | 0.0% / 0.0% |

⇒ ⛔⛔ **Indistinguishable at every scale.** The filter measures "a 167–1,099-residue protein was
placed next to a 1,174-atom target", which is unavoidable at any orientation. ⇒ **it carries no
information about the constellation and must not be used as a gate or quoted as a finding.**

⚠️ Size is not the explanation either: the correlation between chain length and clash count is
**+0.188**, and the smallest quartile of matched chains (≤ 357 residues) still buries **≥ 67**
backbone atoms, median 115.

### ⚠️ A correction to this measurement's own numbers, kept rather than removed

⛔ The actual-arm statistics were first published as **min 52 / median 92 / mean 93**. Those were
computed on the **first 25 rows of a 58-row file that was sorted ascending by clash count.** Over
all 58 rows the correct values are **min 52 / median 129.0 / mean 131.86**.

⇒ ⛔⛔ **The corrected sign is the opposite of what was published.** The original text concluded that
matched placements clash *slightly less* than random (92 against 116). The correct comparison is
**129 against ~116: matched placements clash MORE than random.**

⇒ ⭐⭐ **And the corrected sign has a mechanism**, which is why it is the more believable number:
**65.5% (38 of 58)** of matched chains have their centre of mass on the **target side** of the
constellation centroid, median **+3.67 Å**. Matching a constellation of this kind pulls a scaffold
*into* the target rather than around it.

⛔ **What does not change: 0 of 58 clash-free in both the truncated and the full statistics.** The
filter's retirement stands; only its numbers were wrong. ⚠️ The null's figures are a separate
random-orientation run, so small differences between quoted null medians are sampling, not
disagreement.

### ⭐⭐ Why, in one number

**MEASURED:** the centroid of the six functional pins sits **1.91 Å** from the nearest target heavy
atom — i.e. **inside the target's van der Waals envelope.** For the five-pin subset it is **0.93 Å**,
deeper still.

⇒ ⭐⭐ **That single fact explains every clash result in this line.** Any protein whose functional
groups match the constellation has its body centred on the target surface, so roughly half of it
lands inside the target no matter how it is rotated. ⇒ **the pins do not sit on a face a binder
approaches from outside; they straddle the surface.**

⚠️ **This is also the earliest form of the body-clearance requirement** that `structural_mining.md`
§8.3 states properly. It was visible here as a nuisance in a filter, and it was not recognised as a
requirement until the four surviving candidates were actually placed.

### ⭐ What replaced the filter: how much room is actually there

**MEASURED**, 1.0 Å grid about the constellation centroid; a grid point is free if it is more than
4.0 Å from every target heavy atom; connected components by 3-D labelling:

| *R* (Å) | grid points | free | free % | largest contiguous free region |
|---|---|---|---|---|
| 10 | 4,169 | 1,132 | 27.2% | 1,129 of 1,132 (**100% of free**) |
| 12 | 7,153 | 2,508 | 35.1% | 2,507 of 2,508 (100%) |
| 15 | 14,147 | 6,154 | 43.5% | 6,150 of 6,154 (100%) |
| 20 | 33,401 | 17,733 | 53.1% | 17,724 of 17,733 (100%) |
| 25 | 65,267 | 40,941 | 62.7% | 40,930 of 40,941 (100%) |

⇒ ⭐⭐ **The free space is ONE contiguous region at every radius — never fragmented pockets.** There
is a coherent place for a binder to go.

⚠️ **Scale check, an inference and not a measurement:** a 90-residue all-helical binder occupies
roughly 12,000 Å³. The free volume reaches that between *R* = 15 Å (6,154 Å³) and *R* = 20 Å
(17,733 Å³) ⇒ a design has enough contiguous room only once it extends **17–20 Å** from the
constellation centroid. ⛔ Not measured against a real binder's occupied volume.

### ⛔ A claim that was made here and is withdrawn

An earlier reading of the clash and density results was that the constellation *"describes a pocket
interior"* — the geometry nature builds by surrounding a cavity. ⛔ **That reading is withdrawn,
twice over:** first because it rested on the synthetic-backbone representation (`structural_mining.md`
§2.2), and then because the clash metric it rested on is uninformative against both nulls above.

⛔ Also withdrawn: *"natural scaffolds cannot present this geometry."* They demonstrably can. They
simply do it with bodies far larger than a small binder.

---

## 3. ⛔⛔ The masked density map returns a null: matched scaffolds put no more mass in the usable volume than chance

### The question

The original motivation for mining was not the match count at all. It was: **what do natural
proteins build around this geometry?** The designs in hand were criticised for being too simple —
the only thing they do is connect one point to another, where natural proteins have evolved to do
more. If matched scaffolds shared some common structural element positioned away from the fixed
residues, that element would be design guidance.

So: superimpose every match into the target frame, and average the atomic density relative to the
fixed positions. No secondary-structure assignment — just where the mass is.

Two things had to be fixed first:

1. **The target-occupied volume must be masked out.** Mass placed where the receptor already is
   cannot be part of a binder.
2. ⭐⭐ **It needs a null**, because without one the measurement is uninterpretable — if most chains
   host the motif, a raw density map is a picture of the database, not of the site. §2 is the proof
   of how badly this goes wrong.

### ⭐⭐ The null, which is the good part of this design

**Rotate each actual superimposed chain by a random rotation about the constellation centroid.**

This preserves every atom's **distance** from the centroid exactly, so the radial profile of null and
actual are identical **by construction**, and destroys only the **angular** information. Any
directional structure surviving in the actual density is therefore real, and not an artefact of chain
size, standoff or radial distribution.

⭐ **Sanity check passed:** radial profiles agree shell by shell (0–5 Å: 3.75 against 3.82; 5–10:
19.89 / 19.91; 10–15: 54.59 / 54.69; 15–20: 77.12 / 77.13; 20–25: 83.45 / 83.19). The null behaves
exactly as designed.

### Setup

56 matched chains from the five-pin subset, **50 null rotations each** (2,800 null placements), 2 Å
voxels out to *R* = 25 Å, free = more than 4 Å from any target heavy atom ⇒ **5,145 of 8,217
in-sphere voxels free (62.6%)**.

### ⛔⛔ The result: no signal

**MEASURED**, mean backbone atoms per chain, **free volume only**:

| shell (Å) | actual | null | ratio |
|---|---|---|---|
| 0–5 | 0.00 | 0.00 | — |
| 5–10 | 4.46 | 4.95 | 0.90 |
| 10–15 | 25.30 | 26.01 | 0.97 |
| 15–20 | 42.87 | 46.43 | 0.92 |
| 20–25 | 56.68 | 61.13 | 0.93 |
| **free total** | **129.32** | **138.51** | **0.934** |
| target-occupied | 109.48 | 100.22 | **1.092** |

Half-space split, free volume only, along centroid → target: actual 25.7% / 74.3%; null 23.2% /
76.8%. Voxel-level enrichment (actual / null) over free voxels: median **0.78**, p90 **2.32**, max
10.53.

⇒ ⛔⛔ **Constellation-matched scaffolds put no more mass in the binder-accessible volume than
randomly oriented ones do — the ratio is 0.934, marginally LESS.** The half-space split is within
2.5 points of the null.

⇒ ⭐⭐ **There is no common structural element positioned away from the fixed residues.** That is the
direct answer to the question the whole line was built for, and it is negative.

⭐ The one consistent asymmetry runs the other way: the actual placements put **1.09×** the null's
mass in the **target-occupied** volume. Matching the constellation pulls a scaffold *into* the target
— the same geometry the 0.93–1.91 Å centroid depth explains, and the same effect as the corrected
clash sign in §2. ⇒ ⭐ **the depletion in the free volume and the enrichment in the occupied volume
are one effect, not two.**

### ⛔ Scope

**Does not rule out:** a signal at the four-pin subset (a much larger population); a signal in
secondary-structure **type** rather than mass placement, which is a different measurement; a signal
for a different arm of the design; a signal at a tighter RMSD.

**Does rule out, for this subset and sample:** that the 56 closest natural chains share a common
shape in the volume a binder could occupy. ⚠️ With p90 enrichment 2.32 and median 0.78 there is no
broad directional preference to find; individual high-enrichment voxels exist (max 10.53) but at
*n* = 56 are expected from noise and ⛔ **must not be mined after the fact.**

---

## 4. What the three negatives add up to

⭐⭐ **Correctly executed, with controls and nulls, structural mining did not supply the guidance it
was built to supply.** The Cα search had no selectivity; the clash filter had no discriminating
power; the density map had no signal against its own null.

⛔ What the line did produce is in `structural_mining.md`: a validated matcher, a corrected
tolerance algebra, a named hit list, and three stacked findings that close the scaffold-grafting
route. ⛔ What it did **not** produce is a single usable scaffold.
