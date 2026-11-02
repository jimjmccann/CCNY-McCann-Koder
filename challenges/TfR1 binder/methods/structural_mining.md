# Structural mining: searching the structural record for a scaffold that already presents a designed binding site

This is a complete, self-contained account of the method. It should be possible to reimplement
everything below from this text alone. The programs that were actually run are in `../scripts/`,
and the tables they produced are in `../data/`.

**Work dated 2026-10-30 .. 2026-10-31.** All of it is local CPU. No structure prediction, no model
weights, and no network access is needed for the search itself — only to obtain the databases
(§3).

---

## 1. The question

A designed binding site can be written down as a small set of side-chain **functional groups** held
at fixed relative positions in space — here, four arginine guanidinium groups and two histidine
imidazoles, six groups in all, spanning 32.55 Å at their widest pair.

Two things could then be true, and they have different consequences:

- **The geometry is impossible.** No protein presents those six groups in that arrangement, in
  which case the design is asking for something nature does not build and the specification should
  change.
- **The geometry exists.** Some protein already presents it, in which case that protein is a
  candidate **scaffold** — a backbone one could graft onto or redesign, instead of inventing one.

⇒ The search is for chains whose own functional groups reproduce the designed constellation. It is
**identity-agnostic about everything except functional-group type**: it does not care what fold the
host has, what it does, or what sequence it carries.

⛔ **A match is a geometric coincidence, not a scaffold.** Sections 6–8 are three successive gates
that test whether a match could actually be used, and all three are where the method's real content
lies.

---

## 2. Defining the constellation

### 2.1 What is pinned, and what is not

The design specification constrains **only the functional tip atoms** of each guidepost residue:

| group type | residue | atoms constrained | point used in the search |
|---|---|---|---|
| `CATION` | arginine | `NE`, `NH1`, `NH2` | centroid of those three atoms |
| `CATION` | lysine | `NZ` | the `NZ` atom itself |
| `IMID` | histidine | `ND1`, `NE2`, `CG` | centroid of those three atoms |

The **backbone is not constrained**, by design. Where the Cα sits is an output of the design
process, not an input to it.

### 2.2 ⛔⛔ Why matching on backbone atoms is the wrong instrument, and the measurement that proves it

The first two attempts at this search matched on **Cα**, then on **Cα + Cβ**. Both failed, and the
cause is worth stating plainly because it is a general trap.

**MEASURED:** the guidepost coordinates in the design input carry one single idealised rotamer.
Every arginine pin has χ1 = −176.4°, χ2 = 177.2°; every histidine pin has χ1 = −66.9°, χ2 = 90.2°;
the spread across pins is **less than 0.07°**. Those backbone positions were never specified by the
designer — they are whatever the builder program emitted. ⇒ **matching on them matches the
builder's default rotamer, not the design.**

**MEASURED, independently, over 13,058 arginines in real structures:** if you superimpose a
guanidinium group *exactly*, the two Cα atoms agree to within 1.0 Å only **4.5%** of the time. Side
chains rotate; a guanidinium moves several ångström under χ rotation with the Cα unmoved.

⇒ ⛔⛔ **A chemically perfect match fails a Cα gate essentially always.** The search must be run on
the functional groups and the backbone must be left free.

⭐ It is also far cheaper. A 250-residue chain holds roughly 30–35 functional groups instead of 250
residues, which shrinks the space of 6-element subsets by about **five orders of magnitude**. A
full-store search takes seconds rather than hours.

### 2.3 The query as data

The query is reduced to:

- a **type vector** — one of `CATION` / `IMID` per pin, in a fixed order;
- the **pairwise distance matrix** between the six points (15 distances);
- the **point coordinates**, kept for the final superposition.

⭐ Recording the pin **order** matters and is easy to get wrong. Every downstream table — matched
residue numbers, matched types, the superposition's point correspondence — is written in that one
order. Fix it once, deterministically (here: the sorted order of the specification's own pin keys),
and assert it everywhere.

---

## 3. The two databases

Two stores were built, from independently sourced material. Each store is a flat array of
functional-group centroids plus a type array, a Cα array and a residue-number array, sharded by
chain, with a sidecar file of chain names.

| | PDB-derived store | predicted-model store |
|---|---|---|
| source | 25,821 experimentally determined mmCIF entries from the **Protein Data Bank**, held locally as two previously mined subsets | 15 whole prokaryote proteomes from the **AlphaFold Protein Structure Database**, v6 |
| filter | none beyond what the subsets already applied | per-residue pLDDT ≥ 70 |
| **chains** | **85,699** | **49,481** |
| functional groups | 2,927,334 | 1,778,798 |
| groups per chain | 34.2 | 35.9 |
| build time | 64 s | 77 s |

⚠️ The predicted-model download was verified before use: 15 of 15 proteome archives present, every
archive size matching the distributor's own manifest exactly, 51,680 predicted structures. pLDDT is
carried in the B-factor column (checked on a sample chain: min 46.2, median 92.8, max 96.9). ⚠️ Each
such archive ships both `.cif.gz` and `.pdb.gz`, so roughly half of what is transferred is a format
this method does not use.

⚠️ **One library trap, recorded because it fails silently.** In the structure-parsing library used
here, reading a structure from an in-memory buffer is not available in the installed version; the
call must go through the CIF string reader and then build a structure from the block. The direct call
failed **inside a try/except** and produced a **zero-chain store** that looked like a successful
build. ⇒ assert a non-zero chain and atom count immediately after any store build.

### 3.1 ⛔⛔ The two stores are NOT fully independent, and the agreement between them must not be oversold

The PDB-derived store carries **crystallographic redundancy** — homologues, re-depositions,
multiple chains of one assembly. The predicted-model store is predicted monomers built from UniProt
sequences, many of which have PDB relatives. ⇒ the two stores sample **overlapping biology by
different routes**.

⇒ ⭐ Agreement between them is a **strong cross-check on the method** and a **weak one on nature's
actual frequencies.** ⛔ **A percentage from either store is not a frequency in nature**, and a
rank correlation between them is not evidence that it is.

⇒ ⛔⛔ **And the sharper form of the same point, which is the lesson worth carrying:** a
cross-database replication **cannot validate a parameter the two arms share.** Both arms ran at the
same inadmissible tolerance (§5) and therefore agreed *while both were wrong*. What caught that
error was algebra, not replication.

### 3.2 A data-hygiene defect found in the PDB-derived store, recorded rather than fixed

**MEASURED:** the 795-row six-group hit list from that store holds **786 distinct chain names** —
nine names appear twice, because the same chain entered from two overlapping mined subsets. ⇒ the
store counts **entries, not unique chains**.

⚠️ This does not invalidate the curves in §5: both curves count store entries the same way, and the
hit emitter reproduced them to the digit. ⛔ But "795 unique chains" is wrong. It is **795 store
entries over 786 distinct names.**

---

## 4. The matching criterion and the search algorithm

### 4.1 Acceptance

A candidate is a choice of *k* functional groups from one chain, one per query pin, with **matching
type** (a cation pin can only be served by a cation group). It is accepted if the optimal rigid-body
superposition of those *k* points onto the *k* query points has

```
RMSD  ≤  tol_rmsd        (2.0 Å throughout this work)
```

The superposition is a standard Kabsch alignment, with **proper rotations only** — the reflection
case is corrected by flipping the sign of the smallest singular value, not accepted. Reflections
would match a mirror-image arrangement, which no protein can present.

⛔⛔ **The acceptance bar travels with every count this method produces.** At 2.0 Å the six-group
constellation is presented by about 1% of chains; at 1.0 Å **both stores return zero**. A count
quoted without its bar is not a result.

### 4.2 The pre-filter, and why enumeration needs one

Checking every *k*-subset is hopeless by brute force: a 250-residue chain with 35 functional groups
holds C(35,6) ≈ 1.6 million six-subsets, and the whole store is tens of thousands of chains. So the
candidate set is grown incrementally under a **pairwise-distance pre-filter**: a pair of chain
groups may serve pins *(a,b)* only if

```
| d_chain(i,j)  −  d_query(a,b) |  ≤  tol_d
```

The implementation:

1. Compute the chain's full group–group distance matrix.
2. For each of the C(k,2) query pairs, build a boolean "shell" matrix of chain pairs that satisfy
   the pre-filter **and** have compatible types on both ends. If any shell is empty, the chain
   cannot match — return immediately.
3. ⭐ **Choose the growth order by population, not by distance.** Start from the **least-populated**
   shell *in that particular chain*, then repeatedly append the pin whose shells against the
   already-placed pins are smallest in total.
4. Grow partial assignments in batched array operations, masking out any chain group already used.
5. Superimpose the survivors in blocks and keep those with RMSD ≤ `tol_rmsd`.

⛔ **Step 3 is load-bearing and the naive alternative is 14× slower.** An earlier version anchored
on the *tightest* query pair first. That is backwards: in a dense protein the short-distance shell
is the most heavily populated — every residue has near neighbours — so anchoring there *maximises*
the branching factor. Measured: 300+ ms per chain, versus 21.9 ms per chain for the
least-populated-shell order.

### 4.3 ⛔ An overflow guard that must be instrumented, and why

The implementation aborts a chain whose intermediate candidate set exceeds a cap (default
2,000,000) by returning **an empty list** — which is the same value it returns for *no match*.

⇒ ⛔⛔ **A cap overflow and a true zero are indistinguishable in the output.** Any driver built on
this matcher must count overflows **separately from zeros** and report them per cell. A re-derived
count that does not report an overflow count is not admissible evidence.

⚠️ **MEASURED, and stated as a prediction that failed:** across all 84 cells of each store's curve
sweep and all 12 control cells, overflows were **zero**. The defect is therefore **latent, not
active** — it corrupted no published number. ⛔ It is still a real defect.

---

## 5. ⛔⛔ The tolerance-floor error: the pre-filter was set tighter than the criterion allows

This is the central methodological finding of the whole line, it was our own error, and it is
transferable to any substructure search of this shape.

### 5.1 What went wrong

`tol_d` was set to **1.5 Å** while `tol_rmsd` was **2.0 Å**. The published result at that setting
was that the six-group constellation is matched by **0 chains, 0.00%** — reported as *"the geometry
is unsatisfiable"*, a statement about nature.

⛔⛔ **It was a statement about the pre-filter.** `tol_d` gates only a pre-filter; the acceptance
test is applied afterwards. A pre-filter that is tighter than its own criterion **discards
candidates the criterion would have accepted.**

### 5.2 The algebra, derived

For a Kabsch superposition of *k* points with RMSD ≤ *r*, the sum of squared per-point deviations
satisfies

```
Σ δ²  ≤  k · r²
```

A pairwise distance between two of those points can deviate by at most `δ_i + δ_j`. Maximising
`δ_i + δ_j` subject to `δ_i² + δ_j² ≤ k·r²` gives `δ_i = δ_j = r·√(k/2)`, so

```
max | Δd_ij |  =  2 · r · √(k/2)  =  √(2k) · r
```

⇒ **the admissible floor is `tol_d ≥ √(2k) · tol_rmsd`, and it is per-k.**

| *k* | admissible `tol_d` floor at `tol_rmsd` = 2.0 Å | what was used |
|---|---|---|
| 3 | 4.899 Å | 1.5 Å ❌ |
| 4 | 5.657 Å | 1.5 Å ❌ |
| 5 | 6.325 Å | 1.5 Å ❌ |
| **6** | **6.928 Å** | **1.5 Å** ❌ |

⛔ **Every row of the published curve was below its own floor, not only the six-pin row.** ⛔ And
never a single constant across *k*: the floor scales as √(2k).

### 5.3 ⭐⭐ The invariant that makes the damage diagnosable without re-running anything

Because `tol_d` gates only a pre-filter, raising it can **only admit** candidates previously
discarded — it can never reject one it previously kept.

⇒ ⭐⭐ **Hit counts are monotone non-decreasing in `tol_d`.** Therefore **any count measured at a
too-tight `tol_d` is a lower bound**, and conclusions of the form *"this is common"* are safe while
conclusions of the form *"this is absent"* are the ones at risk. That was written down before the
corrected run and it held: 0 monotonicity violations in 42 of 42 subsets, in each store.

⭐ This is a cheap, general check on any search result — and it tells you in advance which of your
conclusions need re-litigating.

### 5.4 The corrected result

**MEASURED**, all 42 subsets of the six pins (C(6,3..6) = 20+15+6+1), each at both `tol_d = 1.5 Å`
and its per-*k* admissible floor, `tol_rmsd = 2.0 Å` throughout. Best subset per pin count:

| pins | at 1.5 Å (inadmissible) | at the floor — predicted-model store | at the floor — PDB-derived store | best RMSD anywhere |
|---:|---:|---:|---:|---:|
| 3 | 81.41% | **95.30%** (47,153) | **95.96%** (82,236) | 0.009–0.011 Å |
| 4 | 15.17% | **71.29%** (35,273) | **69.64%** (59,682) | 0.173 / 0.197 Å |
| 5 | 0.11% (56 chains) | **22.00%** (10,885) | **21.51%** (18,432) | 0.545 / 0.572 Å |
| **6** | **0.00% (0 chains)** | **1.05% (519)** | **0.93% (795)** | **1.108 / 1.068 Å** |

⇒ ⭐⭐⭐ **"The six-group geometry is unsatisfiable" is dead as a statement about nature.** The wall
was the filter.

⇒ ⭐⭐ **And the regime claim moved by exactly one pin.** Five pins went from 56 chains ("a
shortlist you can inspect by hand") to 10,885 chains — a population. Six pins is now the
hand-inspectable shortlist that five was claimed to be. Even the *worst* of the six five-pin subsets
rises from 18 chains to 3,246, so the shift is not an artefact of picking the best subset.

⇒ ⛔ **The decline is also much gentler than published.** 95.30% → 1.05% is under two orders of
magnitude and it is smooth; the original "sharp four-decade transition across one added constraint"
was an artefact of the same filter.

### 5.5 ⛔⛔ The constraint that stops this being good news

**MEASURED:** widening the pre-filter produces **more loose matches, never a better one.** Per pin
count, the best RMSD achieved is **identical** at 1.5 Å and at the floor. At six pins the best match
that exists anywhere is **1.068 Å** (PDB-derived) / **1.108 Å** (predicted models), and the whole
six-pin population sits between about 1.07 Å and the 2.0 Å bar.

⚠️ **Read that per pin count, which is what it measures.** At the level of an individual
(*k*, subset) cell the floor *is* better in 3 of 84 cells — never worse, identical in the other 79
— because at a too-tight `tol_d` a *k*-pin match is only *guaranteed* to pass if
`r ≤ tol_d/√(2k)`, so widening can reveal a better match. No per-*k* number and no curve value
moves.

⇒ ⛔ The design question is therefore no longer *"can nature present this?"* but *"is a ~1.1 Å
constellation match good enough to build on?"* — which this work does not answer and for which it
proposes no number.

### 5.6 The uniqueness argument that shared the defect and fell with it

An earlier positive control had argued that the zero was a property of the *regime*, not of our
design: taking a **real** chain's own six functional groups as the query also returned 0 hits in
20,000 other chains, so *"a six-group typed constellation spanning more than 30 Å is essentially
unique in nature."*

⛔⛔ **That control was measured at the same inadmissible tolerance as the zero it was written to
explain.** **MEASURED**, four real control constellations, each searched against 20,000 other
chains: all four reproduce 0 hits at 1.5 Å, and at the admissible floor they return **286 / 14 /
1,242 / 156** chains, best RMSD 1.090–1.386 Å.

⇒ ⛔ *"Essentially unique in nature"* is **refuted**, and the honest version is *"our filter hid
these constellations — for our design and for nature alike."*

⇒ ⭐⭐ **The cleanest statement of the defect: the rescue and the thing it rescued were one
artefact.** When an anomaly and the control that explains it away share a parameter, the control
proves nothing.

⚠️ **Scope on that control:** four chains, chosen deterministically (first *N* with span ≥ 30 Å),
and they fall into two store-adjacent pairs, so they are nearer *n* = 2 than *n* = 4. The direction
is not in doubt — every cell is non-zero by one to two orders of magnitude — but the **magnitudes
must not be averaged** across them.

### 5.7 Gates that any re-derivation of this kind needs

Three gates were run before any corrected number was believed, and all three are reusable:

1. ⭐ **Faithfulness.** The matcher could not be instrumented for overflows without modifying it,
   so its core routine was copied verbatim with a flag added, and the copy was then **proved
   equivalent to the original**: 4 pin subsets × 2 tolerances × 400 chains = **3,200 paired
   comparisons**, 0 mismatches in both hit count and best RMSD. The original was imported, never
   edited, and stays the reference instance.
2. ⭐⭐ **Reproduction.** At the old, inadmissible tolerance the new driver must reproduce the
   published numbers **exactly** — counts, best RMSD, and the identity of the best subset at each
   pin count. It did, to the digit, on all four rows in both stores. ⇒ the arithmetic and the
   matcher were never wrong; **one parameter was.** ⛔ Had this gate failed, nothing else would have
   been reported.
3. ⭐ **Monotonicity** (§5.3), checked on every cell rather than on the best one.

⚠️ **One trap inside the reproduction gate, caught before it fired.** The gate was written to
require a best RMSD of "0.5720 Å" because the published source printed **0.572**, extended to the
driver's four decimals. The true value is **0.5718**. Enforcing the invented fourth decimal would
have failed the gate on a rounding artefact and sent someone looking for a non-existent database
mismatch. ⇒ ⛔ **never state a reproduction gate to more precision than its source carries.**

---

## 6. Gate one: the hits must be NAMED

⛔⛔ **The curves in §5 were published as counts only.** The output tables carried
`k / subset / tol_d / hits / chains / pct_store / best_rmsd` and **no chain identity**, so not one
of the 519 + 795 matching chains was known by name.

The identities had been **computed and discarded**. The matcher already returns a chain index and a
best RMSD per hit, and a chain-name sidecar sits beside each store; the driver summed them into a
count.

**MEASURED:** emitting the named list took **7.1 s** on one store and **4.3 s** on the other.

⇒ ⛔⛔ **Two days of searching produced a statistic and zero candidate structures, and the fix cost
eleven seconds of local CPU once it was asked for.** That is a process failure, not a measurement
one, and it is in `../FINDINGS.md` as the second transferable lesson.

The emitter writes one row per chain at its **best-RMSD placement**, ranked ascending, carrying the
six matched residue numbers and their six types — which is what lets a hit be pulled out of the real
structure and inspected. ⛔ Its own gate: the chain count and best RMSD it finds must equal the
published six-pin floor row for that store, or it is not the instrument that produced the curve and
its list must not be trusted. **PASS** on both stores, to the digit (795 / 1.0680 Å and
519 / 1.1079 Å).

### 6.1 ⛔⛔ What the named list says: size and fit are opposed

**MEASURED**, chain size of the six-group hits, as a **true residue count read from each
structure**. ⛔ An earlier version of this table binned by *last residue number*, which is not a
length and is not a bound on one: a chain numbered 296–411 reads as 407 by that measure and is
116 residues long. The numbers below replace it; see §6.1.1.

| length | PDB-derived hits | % | predicted-model hits | % |
|---|---:|---:|---:|---:|
| **≤ 120 aa** | **2** | **0.3%** | 0 | 0.0% |
| 121–250 aa | 71 | 8.9% | 21 | 4.0% |
| 251–500 aa | 372 | 46.8% | 233 | 44.9% |
| **> 500 aa** | **350** | **44.0%** | 265 | 51.1% |

⇒ ⛔⛔ **Only 2 of the 795 PDB-derived chains presenting the six groups are under 120 residues, and
91% are over 250. Across both databases the shortest chain that presents the constellation at all
is 116 residues.** And the small ones are the **worst** fits: the best hit at ≤ 250 aa is
**1.4662 Å**, against the list's overall best of **1.0680 Å**.

#### 6.1.1 ⛔ The correction, recorded rather than quietly applied

The first pass measured chain length as the highest residue *number* among a chain's matched
functional groups. For a chain numbered from 1 that approximates the length; for a domain or a
fragment numbered in the context of a larger parent it does not, and the error has no fixed sign.
The replacement counts amino-acid residues in the chain directly. **The conclusion of this section
survives the correction** — two chains under 120 residues out of 795 is still essentially none, and
the opposition between size and fit is unchanged — but the specific claim *"not one is under 120"*
was false, and the band percentages all moved.

⇒ ⭐ **Size and fit are opposed across this hit set.** The cause is the span, and no sub-constellation
escapes it: six pins span 32.55 Å, the most compact five-pin subset 30.27 Å, the most compact
four-pin subset 27.21 Å. Dropping two pins buys 5.3 Å.

---

## 7. Gate two: the accessibility ("face") gate

⛔ **A constellation buried in a 500-residue core is not a binding site.** Until that is measured, a
hit remains a geometric coincidence.

Measured per hit chain, **on the chain in isolation**, because a scaffold would be used standalone.

### 7.1 Channel A — exposure

Relative solvent accessibility of each of the six matched residues, computed with **freesasa**
rather than a hand-rolled implementation. A residue counts as exposed at **RSA ≥ 0.20**, the
conventional line.

⚠️ Implementation detail that matters: the isolated chain is written out with its chain renamed to a
single character and its residues renumbered sequentially before the SASA call, which sidesteps the
PDB format's four-digit residue-number and single-character chain-name limits. Several hit chains
exceed both.

```
A  =  how many of the six matched residues reach RSA ≥ 0.20      (0..6)
```

### 7.2 Channel B — co-accessibility from one direction

Exposure one residue at a time is not enough: a partner has to reach **all six from one side**.

- Sample **642 directions** on the sphere (Fibonacci lattice).
- For a direction **u** and a matched group, cast the outward ray from the group along **u**. The
  group is **visible along u** if no chain heavy atom lies within **3.0 Å** of the ray between
  **1 Å and 15 Å** from the group.
- ⭐ Exclude the group's **own** residue from the occlusion test — self-occlusion is spurious. ⭐
  **Keep the other five pins' atoms**, because a partner really is blocked by them.

```
B  =  max over the 642 directions of (number of groups visible along that direction)    (0..6)
```

Also recorded, descriptively: chain length, the carried best RMSD, and the maximum pairwise angle
between the six outward vectors.

### 7.3 The decision rule, fixed before the run

```
CANDIDATE   A ≥ 5  and  B = 6
PARTIAL     A ≥ 4  and  B ≥ 5
FAIL        otherwise
```

⛔ Fixed in advance and not tuned after seeing output.

### 7.4 ⭐⭐ The two gates that had to pass before any verdict was believed

1. ⭐⭐ **Vacuity — a gate that can only return FAIL is not a gate.** Re-run in a positive-control
   mode that swaps the matched residues for the **six most-exposed ARG/LYS/HIS of the same chain**:
   **11 of 12 return CANDIDATE**, with A = 6, B = 6 and minimum RSA 0.675–0.969. ⇒ ⭐ **the pass
   branch is reachable**, so every FAIL below is a verdict rather than an artefact of an
   unsatisfiable rule.
2. ⭐⭐ **Identity — the residue-type verification guard.** Hard assertion: the residue found at a
   recorded residue number **must** match the recorded pin type (`CATION` ⇒ ARG or LYS, `IMID` ⇒
   HIS). A mismatch means the lookup landed on the wrong residue, and that chain is recorded as a
   lookup failure and excluded from every statistic rather than silently measured. Rule fixed in
   advance: at a 5% failure rate or more, the face numbers must not be published at all, because the
   identity mapping would be broken. **MEASURED: 0 lookup failures in 1,314.**

⚠️ The structures are read by the search code's own read path — **gemmi** `read_structure`,
alternative conformations and hydrogens removed, chain selected by **author** name, residue by
`seqid.num` — so a recorded residue number means the same thing in the real structure as it meant in
the store. ⛔ Chain-stem splitting
must split on the **last** underscore: a predicted-model identifier contains underscores, and
splitting on the first one looks up the wrong chain.

### 7.5 ⛔⛔ The result: accessibility and fit are opposed

**MEASURED**, both stores, all 1,314 hits:

| verdict | n | % |
|---|---:|---:|
| **CANDIDATE** | **31** | **2.4%** |
| `PARTIAL` | 143 | 10.9% |
| `FAIL` | **1,140** | **86.8%** |

⇒ ⛔⛔ **86.8% of the chains that present the six groups present them where nothing can reach
them.** The 31 survivors collapse to **19 distinct PDB entries plus 5 predicted models**.

**MEASURED, and this is the sharp part:** best RMSD among **failed** hits **1.0680 Å**; best RMSD
among **candidates** **1.5400 Å**.

⇒ ⛔⛔ **The whole tight head of the hit list — every best-fitting match — is buried or
unreachable.** No candidate fits better than 1.5400 Å, and candidates run out to 1.9964 Å against
the 2.0 Å bar.

⇒ ⭐⭐ **Fit and accessibility are opposed, the same shape as size versus fit.** A reading
consistent with both: a 32.55 Å constellation of four cations and two imidazoles is, in nature,
usually an **interior or active-site arrangement**, not a surface patch.

**By chain length:**

| length | measured | CANDIDATE | % | PARTIAL |
|---|---:|---:|---:|---:|
| ≤ 120 aa | 2 | 1 | 50.0% | 1 |
| **121–250 aa** | 92 | **6** | **6.5%** | 24 |
| 251–500 aa | 605 | 5 | 0.8% | 63 |
| > 500 aa | 615 | 19 | 3.1% | 55 |

⚠️ The prediction that smaller chains would be more accessible **held weakly and
non-monotonically**: 121–250 aa does have the highest rate (6.5% against 3.1% above 500 aa), ⛔ but
251–500 aa is the *lowest* at 0.8%, so this is not a clean trend and ⛔ **must not be quoted as
"smaller is more accessible"**. The ≤ 120 aa band is **n = 2** and its 50% is one chain, ⛔ not a
rate.

### 7.6 What survives: seven rows, four proteins

Every CANDIDATE at ≤ 250 residues — small enough to be a binder scaffold **and** on a reachable
face:

| store | chain | RMSD (Å) | length | A | B | min RSA |
|---|---|---:|---:|---:|---:|---:|
| PDB-derived | **`5xox_B`** | **1.5400** | 229 | 5 | 6 | 0.173 |
| PDB-derived | `5xox_E` | 1.5477 | 229 | 6 | 6 | 0.211 |
| PDB-derived | `5xox_A` | 1.7315 | 229 | 5 | 6 | 0.186 |
| PDB-derived | **`9gyt_N`** | 1.8251 | **116** | 6 | 6 | **0.512** |
| PDB-derived | `5xox_D` | 1.8331 | 229 | 6 | 6 | 0.218 |
| predicted model | `AF-P56144` | 1.9216 | 242 | 5 | 6 | 0.135 |
| PDB-derived | **`9t6h_b`** | 1.9375 | 222 | 6 | 6 | 0.362 |

⛔⛔ **Four of the seven rows are chains A/B/D/E of `5xox` — one protein in four copies.** ⇒ the
honest count is **four distinct proteins**: `5xox` (229 aa), `9gyt_N` (116 aa), `9t6h_b` (222 aa),
`AF-P56144` (242 aa). `9gyt_N` is the most exposed of the set and the only one near the original
design's own length regime.

⭐ The next tier — PARTIAL at ≤ 250 aa — is nine more distinct entries (`3ozf_C` 1.5240 ·
`7TUX_B` 1.6016 · `8i09_B` 1.7158 · `9kyl_B` 1.8304 · `2V5B_A` 1.9109 · `3PIW_A` 1.9552 ·
`4K2P_A` 1.9608 · `6b0p_D` 1.9794 · `3U8I_B` 1.9969), ⚠️ all with B = 5, i.e. **one group
unreachable from any single direction.**

---

## 8. Gate three: the superposition, and the body-clearance requirement

The four surviving candidates were superimposed onto the constellation **in the target's own
coordinate frame**, which is what finally tests whether a candidate could be used.

### 8.1 How the superposition is defined, and the gate on it

- **Kabsch on six points only** — each candidate's matched functional-group centroid onto the
  corresponding guidepost centroid, in the query's own fixed pin order (§2.3), which is the order
  the hit list's matched-residue columns were written in.
- Centroids computed by the **same recipe as the search**: arginine as the mean of `NE`, `NH1`,
  `NH2`; lysine as `NZ`; histidine as the mean of `ND1`, `NE2`, `CG`.
- ⭐⭐ **The gate:** the RMSD achieved must equal the best RMSD recorded for that chain in the hit
  list, to **1e-3**, or the builder **refuses to write output**. A superposition that does not
  reproduce the published number is a *different* alignment, and a picture built on it would show a
  measurement that was never made. **All four passed to the digit: 1.5400, 1.8251, 1.9375, 1.9216.**
- ⭐ The residue **name** at each guidepost and at each matched residue is asserted against the
  recorded pin type before anything is drawn, so a silent off-by-one cannot superimpose on the wrong
  atoms and look perfect.

### 8.2 ⛔⛔ The result: every candidate occupies the receptor's volume

**MEASURED**, close **atom pairs** between the placed candidate and the receptor crop within
**2.5 Å**. ⛔ **This is a count of PAIRS, not of distinct atoms:** the builder sums every entry of
the candidate-by-receptor distance matrix below the cut, so one atom with several close partners is
counted once per partner. The number of distinct candidate heavy atoms involved is correspondingly
smaller, and is not quoted here because it was not separately measured.

| candidate | RMSD (Å) | clashing atom pairs | closest approach (Å) |
|---|---:|---:|---:|
| `5xox_B` | 1.5400 | **679** | 0.40 |
| `9gyt_N` | 1.8251 | **1,271** | 0.15 |
| `9t6h_b` | 1.9375 | **2,469** | 0.15 |
| `AF-P56144` | 1.9216 | **1,694** | 0.07 |

⇒ ⛔⛔ **These are not grazing contacts.** Closest approaches of 0.07–0.40 Å mean atoms essentially
superimposed; hundreds to thousands of close atom pairs place the candidate's body inside the
receptor volume. The candidate's body
sits where the receptor is.

⚠️ **Not an artefact of which placement was drawn:** three of the four chains have exactly **one**
placement in the hit list, and only `9t6h_b` has two — and it is the worst-fitting candidate anyway.

### 8.3 ⭐⭐ The reading, and the requirement the method was missing

⇒ ⭐⭐ **The accessibility gate (§7) asked whether a partner can REACH the six groups. It never
asked whether the host protein's BODY stays out of the target's way.** That is a strictly stronger
requirement, and **all four candidates fail it.**

⇒ ⛔⛔ **Net practical yield of the whole line as drop-in scaffolds: zero.** Four loose candidates
that interpenetrate the target is not a campaign. ⇒ ⭐ scaffold grafting from this mining set is not
a shortcut to a better structure, which strengthens the independent conclusion that the live lever
is the **binder length specification**, not a borrowed backbone.

### 8.4 ⛔ What the session file is, and what it is not

A session file is deposited in `../structures/`. It is **a view, not evidence**: nothing in it
ranks, scores or re-measures anything. Six scenes — all four candidates on the constellation, the
constellation alone, and one scene per candidate carrying its own residue count, RMSD against the
2.0 Å bar, clashing-atom count and closest approach.

⛔ **A maintained geometry is never evidence of affinity.**

⭐ Three construction rules that are worth copying, because each of them silently destroys a session
file otherwise:

- No `__main__` guard. Under a headless viewer invocation the module is not `__main__`, and a
  guarded script exits successfully having written nothing.
- Every scene message must be **pure ASCII**. A non-ASCII emphasis mark corrupts the heap and the
  saved file crashes on load.
- The builder refuses to overwrite an existing output.

⭐ And the file is **verified by reloading it in a fresh viewer process** before being handed on — 17
checks, 0 failures. Because the scene-message layout inside a saved session is version-specific, the
verifier **walks the whole saved session and asserts every string in it is ASCII**, which is
stronger than checking scene messages alone: it also covers labels and object names, and it cannot
rot when the viewer moves its internal layout.

---

## 9. What this method does not do

1. ⛔ **It does not measure chemistry beyond group type and position.** A `CATION` pin is served
   equally by an arginine guanidinium and a lysine `NZ`; nothing checks hydrogen-bond geometry,
   protonation, or the local electrostatic environment.
2. ⛔ **Protonated histidine is a cation and is not admitted as one.** Histidine enters only as
   `IMID`. A deliberate first-version choice, never validated.
3. ⛔ **It does not fold anything.** A matched constellation is a geometric coincidence until
   something folds; nothing here does, and nothing here buys a design.
4. ⛔ **It says nothing about whether a candidate's face is large enough to bury a partner**, nor
   about expressibility, oligomeric state or function.
5. ⛔ **The acceptance bar is unexamined.** The work fixed the admissibility of the pre-filter
   *against* a 2.0 Å bar; it did not establish that 2.0 Å is the right bar.
6. ⛔ **Neither store's percentage is a frequency in nature** (§3.1).
