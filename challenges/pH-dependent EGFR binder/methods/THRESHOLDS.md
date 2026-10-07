# THRESHOLDS — every constant, traced to its source

Site-A EGFR binder campaign. Companion to `README.md` and `OBJECTIVES.md`.

**Provenance labels used below:**

| label | meaning |
|---|---|
| MEASURED | measured on real data, with the measurement stated |
| DECIDED | a human decision on the record, with the reasoning stated |
| LITERATURE | taken from a published source |
| CONVENTION | a conventional value with **no justification recorded anywhere in our files**. Stated as such rather than given an invented rationale |

⛔ Four constants below are CONVENTION. They are load-bearing and their origin is
genuinely not in the record. Do not supply a rationale for them; if you change
one, every table already produced becomes incomparable.

---

## T1. Epitope contact: `CUT = 5.0 Å`, min over all atom pairs

**Value:** 5.0 Å, minimum distance over **all** atom pairs, chain A = binder,
chain B = target.

**Provenance: CONVENTION.** Traced back through the scorer lineage — the current
scorer states it is identical to its predecessor in every respect that affects
the number, that predecessor generalises an earlier one-shard script, and **that
earliest script states no rationale at all**. It is a conventional
protein–protein contact cutoff. ⛔ Our files do not say where 5.0 came from.

**Why it must not be changed:** every table the campaign has produced used it.
A different number makes the counts incomparable. The per-residue breakdown
script copies the constant from the scorer rather than redefining it, and both
headers say a change to one must change the other.

**What IS traced:** the *numbering* convention. Residue identity is asserted, not
a number range — because a range check was tried first and refused all 13
shard-sites. The folds are against the full 621-residue target, so "outside
370–475" was an expectation being wrong, not the data. Identity also catches a
target-crop shift (which would move H409 off 409) where a range check cannot.

---

## T2. Fold gate: `sc_rmsd < 1.5 Å`

**Value:** 1.5 Å, binder-onto-binder Cα RMSD.

**Provenance: LITERATURE**, as recorded in the scorer's own header: the first
filter of a published 100,000 → 4 design funnel (Fry, Slaw & Polizzi,
*Nature* 656:237, 2026). ⚠️ Cited from our script header; **the paper was not
independently read for this deposit**.

**Where it is enforced:** the scorer's `--gate` default; the cycle scorer's
`GATE = 1.5`; the arm configuration's `gate_sc_rmsd = 1.5`. Three places, one
value — check all three if you change it.

**What it buys, MEASURED:** arm-wide the gate removes **53% of full-footprint
contact, 422 draws → 200**. 31 of 118 backbones fail the median gate, and 16 of
those 31 had ≥5 full-footprint draws.

---

## T3. Draws per sequence: `10` for this 75-residue construct

**Value:** 10 structure-prediction draws per sequence.

**Provenance: DECIDED from a measured cost/precision table**, and the decision is
explicitly about **how the number is chosen**, not about one campaign's value.

**The rule:** the draw count is a function of **construct size** and is justified
per campaign against the construct actually being folded — token count, ligands,
whether a forced template plus potentials is in play. ⛔ **Never inherited
because the flag happens to be there.**

**What was measured:**

| arm | draws | tokens | forced template | outcome |
|---|---|---|---|---|
| three gate arms | 21 | 149–175 | no | fine — but a **crop**, not the full target |
| four early arms | 5 | 127–176 | no | fine — also crops |
| the 3-draw arm | **3** | **696 + 10 NAG** | **yes** | ✅ 12 shards, 4,200 structures, **0 OOM** |
| the 21-draw arm | **21** | **696 + 10 NAG** | **yes** | ⛔ **5 of 6 workers produced NOTHING** |

⇒ Every prior success at a high draw count was on a small crop. **3** is the only
count ever proven on the full construct; **10 is an interpolation**, and it was
then separately proven for this construct: cycles 0–2 ran 10 draws with **0
out-of-memory events over 1,040 structures**.

**Precision it buys, MEASURED:** **±15 percentage points** on one sequence's
full-footprint rate — enough to separate a >50% sequence from a ~33% pack,
**not** enough to separate 43% from 33%. Downstream ranking must respect that.
A single design's fraction carries **s.e. ≈ 0.16**.

**Why it cannot be raised:** 10 draws already measured **97.6% GPU memory** on a
24 GB card for this construct, and out-of-memory exits rc=0 writing nothing. ⇒
Precision comes from **repeated** 10-draw runs, not a bigger number.

**Two standing rules attached:**
- An unmeasured draw count is launched on **one canary shard first**. A failed
  sample consumes its full wall time and yields nothing, so a wrong count costs
  the whole arm rather than a fraction.
- A draw count may **not** be raised on the strength of a different construct's
  success.

⛔ **A different arm's 3 does not transfer here, and this 10 does not transfer
there.** The two-site arm runs at 3 for a 130–150-residue binder.

---

## T4. Candidate floor: `gate_frac ≥ 0.5`

**Value:** a sequence is a valid candidate only if ≥50% of its draws pass T2.

**Provenance: CONVENTION.** Earliest occurrence is the burn-in table builder,
which **states the definition and not a justification**. ⛔ Our files do not say
where 0.5 came from.

**What IS on the record about it, and it matters more than the value:**
- ⛔ The floor is a **candidate filter, never a term summed into the label.**
- ⛔ **It applies to the REFERENCE as well as to every candidate.** MEASURED: one
  design topped the table at `label` 10.333 on a gate of **3 of 10** passing
  draws with median `sc_rmsd` **1.871 Å** — it fails the fold gate 7 draws in 10.
  Anchoring a cycle on a thin-data outlier flatters everything downstream.
- Below-floor sequences are **kept and flagged** (`passes_floor = 0`), not
  dropped silently.

---

## T5. Site-A footprint: 11 residues

**Value:** `380 382 384 408 409 410 411 412 417 438 465` (WT numbering).

**Provenance: DECIDED** — the switch-lobe footprint from the site-choice decision
record, plus **F380 added** by the decision that selected the generator's
conditioning. Explicitly **not** the retired earlier site and **not** a
pre-decision epitope list that also exists in the record.

**Identity assertions carried in the scorer**, so a numbering error cannot pass
silently: 380 PHE · 409 HIS · 412 PHE · 438 ILE · 417 VAL.

⇒ **11 is also the ceiling of `label`.** That is a property of the footprint
size, and it is why `label` saturated (`OBJECTIVES.md` §1).

⚠️ **A second site's footprint is genuinely ambiguous in our record** — 8
residues in one family of scripts, 7 as generated (one residue dropped), and the
dropped residue is described elsewhere as "site B's F380". The scorer reports
**both** and takes **neither** as the gate. Not site A's problem, recorded so
nobody resolves it by guessing.

---

## T6. F380 as the floor residue

**Value:** `c380 ∈ {0,1}` is a **hard** requirement inside the per-draw `min`.

**Provenance: MEASURED, then DECIDED.** The measurement that identified it:

- contacted in **63.8%** of 600 draws (two cycles + the replicate round) and
  **63.7%** of a later cycle's 380 fresh draws — **weakest of the 11 by ~13
  points**; the strongest, H409, is 98.8%;
- of draws missing the ceiling by **exactly one** residue, F380 is that residue
  in **89.8%** and **81.9%** of the two cohorts;
- the rate **did not move** when the mutation order went 1–3 → 4–6 ⇒
  **architectural**, not a property of the substitutions explored;
- the accepted set the optimizer walked contains **no position that reaches it**.

⛔ **`c380` is NOT a new measurement** — it is `380 in hit`, out of the same
5.0 Å min-over-atom-pairs pass as the other ten residues. ⛔ Do not re-derive it
with a different cutoff; the scorer emits it and the footprint count from the
**same** geometry pass so the two can never disagree.

**Why hard and not weighted:** see the resolution table in `OBJECTIVES.md` §2.
The hard floor gives 31 distinct values over 98 sequences with **zero** zeros,
against the previous metric's 11. The soft variants resolve slightly more and
buy it by weakening the requirement, with no mechanistic justification.

---

## T7. Target sphere radius: 22 Å (site A)

**Value:** 22.0 Å measured to the footprint ⇒ **231 residues in 15 segments**.

**Provenance: DECIDED** ("be generous"), with two stated constraints:
- chosen to **match the other site's 20 Å / 230 residues**, so both arms run at
  the same rate;
- chosen to stay under the generator's **384-residue training crop** with the
  binder included: 231 + 80 = 311.

⛔⛔ **384 is a TRAINING crop, not a memory limit.** A bigger card buys
memory, not validity — past 384 the model runs outside the size distribution it
was trained on.

⛔ **The radius is a GUESS and the exclusion gate (T8) is the only thing that
proves it was generous enough.** Trimming preserves coordinates, not physics: a
residue whose neighbours were cut is artificially exposed.

**Hotspot assertions:** every hotspot is asserted by residue **type** against the
input file and asserted to be inside the sphere, before anything is written —
because a numbering error in this family of scripts "produced a bogus null once
already".

---

## T8. Full-receptor exclusion gate

| parameter | default | provenance |
|---|---|---|
| `--clash` | **2.0 Å** | MEASURED — the voiding metric. The void was declared on binder heavy atoms within 2 Å of omitted receptor: **70%** of site-A and **100%** of site-B designs, medians 1.14 Å and 0.34 Å, n=300 each |
| `--max-fit-rmsd` | **0.5 Å** | MEASURED — the generator returns the target rigid at **0.0402 Å** over 978/978 atoms, so a fit materially above that means the match is wrong. The RMSD is itself a tripwire, reported per design and asserted |
| `--warn-glycan` | **4.0 Å** | DECIDED, and the reason is the point: the crystal resolves only **1–7 sugars per tree** while native trees are ~15–20 residues and sialylated, with the antennae absent from every structure. A design 3.5 Å from a **truncated** tree may be inside the real one ⇒ a glycan distance is a **FLOOR on the clash, never the clash itself**, so there is a WARN tier rather than a pass |
| `--sequon-report` | **15.0 Å** | CONVENTION — reporting radius for per-sequon columns. No justification recorded; it affects which columns appear, not any verdict |

⛔⛔ **The gate MUST superpose, and that is not a style choice.** MEASURED: the
generator recenters its output by **66.93 Å**, the same offset on every model. A
gate that diffs raw coordinates calls every design a catastrophic clash; one that
"passes" everything is equally wrong. Superpose on the residues the design and
the receptor **share**, then measure the binder against everything else.

---

## T9. Sequence-step flags

| flag | value | provenance |
|---|---|---|
| model | SolubleMPNN, `solublempnn_v_48_020.pt` | DECIDED |
| `--omit_AA` | **C** | DECIDED on MEASURED cost — an earlier arm of this lab permitted cysteine and measured **2–3 cysteines per design** in a protein whose wild type has **zero**, never pairing (partners **12–15 Å** apart). Free thiols: oxidation and aggregation liability with no compensating structure. It confounded that design pool badly enough that the question had to be deferred by a formal decision rather than resolved |
| proline | **permitted** | DECIDED — free backbone NH/C=O exist only at helix caps and in loops (mid-helix amides are consumed by the helix's own i→i+4 bonds), and proline is the capping residue. Omitting it forces caps out of residues that cap badly |
| `--bias_AA` | ⛔ **0 in production — the direction was MEASURED WRONG.** See the box below | DECIDED as a direction from an n=8 smoke test (3 of 8 passed the charge gate, **all 5 failures on the acidity bound**), with the magnitude deliberately left open. ⛔ **Then measured to zero at n=3,320** |
| `--chains_to_design` | **A** | the target is context, never designed, so the model sees the interface |
| `--temperature` | **0.1** | CONVENTION (script default). No justification recorded |
| `--seqs` per backbone | **8** (script default); **100** were taken from the submitted backbone | DECIDED per run, not a fixed constant |

### ⛔⛔ T9-bis. The acidic bias: DECIDED, then MEASURED TO ZERO

The decision left the **magnitude** unsettled on purpose — *"tune on a small
batch, measure NCPR, then commit"*. **MEASURED at n = 3,320 sequences (332
backbones × 10 sequences, both sites):**

| | NCPR median | in band [−0.13, −0.02] | too POSITIVE | too NEGATIVE |
|---|---|---|---|---|
| site A | −0.100 | 1133/1460 (77.6%) | 44 | **283** |
| site B | −0.0909 | 1387/1860 (74.6%) | 138 | **335** |

⇒ **The pool already overshoots on acidity ~6:1.** A D/E bias would push more
sequences out the **bottom** of the band. ⇒ ⛔ **`--bias_AA D/E` STAYS AT 0**, and
the n=8 smoke test *"was WRONG about the direction"*.

⭐ **Measure-then-commit earned its keep: a bias set from n=8 would have hurt.**

**Other measurements from the same n=3,320 pass, same source:**
- cysteine-free **3,320 / 3,320** ⇒ the free-thiol hole stays closed;
- **322 sequences (9.7%)** carry an N-X-[ST] sequon **in the binder** ⇒ filtered;
- ⇒ **fold set: 2,260 sequences over 317 of 332 backbones** (in-band +
  cysteine-free + binder-sequon-free). **The filters cost 15 backbones, not 15%.**
- the anchor residue was asserted in **all 332** backbones; rc=0 both sites.

---

## T10. Net charge band: `−0.13 ≤ NCPR(7.4) ≤ −0.02`

**Definition:** `NCPR = (#R + #K + f·#H − #D − #E) / N` over the **binder chain
only**, with **f = 0.5** for histidine at pH 7.4.

**Provenance: DERIVED from a DECIDED target.** The brief's band is net
**−0.071 per residue** with a cap at **−0.094**; the gate used in the sequence QC
is the band above. ⚠️ The exact arithmetic turning −0.071/−0.094 into
−0.13/−0.02 is **not written down in our files** — the band is wider than the
target on both sides. Treat −0.13/−0.02 as the implemented gate and
−0.071/−0.094 as the target it was built around; ⛔ do not present the band as
derived until someone reconstructs it.

**f = 0.5 is a CONVENTION chosen for consistency** — it matches the separate
charge-gate script's convention specifically so the two agree. Stated as such in
the script.

**What this gate cannot do:** the largest-positive-patch metric is **not** here
and cannot be — it needs per-atom solvent accessibility on a folded structure, so
it runs after the fold-back. ⚠️ And whether that patch value is used as a **gate**
or as a **rank with the worst decile cut** is deliberately **open** in the record
(its own author recommended the rank).

---

## T11. N-glycosylation sequon: `N-[^P][ST]`

**Provenance: LITERATURE / third-party practice**, recorded as a free sequence
filter applied as routine by a frontier lab, and explicitly **not covered by
omitting cysteine**. Same cost, same source. Implemented as a regex on the
binder sequence.

---

## T12. Shard and campaign sizes (site-A sphere arm)

| constant | value | provenance |
|---|---|---|
| binder lengths | **60, 65, 70, 75, 80** | DECIDED. The window for the other site was knocked down 10 because it is a smaller site |
| shards per length (`REPS`) | **3** | DECIDED, and **only safe because inference is unseeded** — verified in source: the seed defaults to None and is only written into output metadata; seeding appears only in training code. If it were seeded, three shards at one length would be byte-identical and the campaign would pay 3× for one shard |
| shards total | **15** | MEASURED — `arm_shard00..arm_shard14` on disk |
| designs per shard | **336** (`n_batches = 42`) | DECIDED against a measured rate: 68.0 min/shard at 374 tokens |
| designs per arm | **5,040** (1,008 per length) | arithmetic: 3 × 5 × 336 |
| submitted backbone | `arm_shard11_shard11_24_model_2`, **75 aa** | MEASURED — `arm_shard11`'s contig reads `75,...`, consistent with the optimizer's `seq_len = 75` |

⚠️ **Correction recorded:** an earlier internal inventory described this stage as
"31 shards, one fixed length each, 70–100 aa". That is a **different, earlier
shard builder** (31 shards at lengths 70..100) and **not** the generator of the
submitted backbone.

⚠️ **The per-shard time bound is explicitly arbitrary.** A ~1 h bound was raised
to ~2 h when the real risk was computed, on the record, *because* the original
was arbitrary — it exists to cap the work that is unprotected by checkpointing.
⛔ Re-derive it from your own costs; never inherit the number. A shard sized past
~2 h is not covered by that decision at all.

---

## T13. Proposer constants

| constant | value | provenance |
|---|---|---|
| `--min-hub-margin` | **3** | CONVENTION (script default). ⛔ No justification recorded for 3 specifically. What IS recorded: the MEASURED margin on the real table was **50 neighbours against 4**, i.e. nowhere near the threshold, and the script **refuses rather than guessing** if the margin is ever thin, because a wrong reference silently redefines every mutation in the cycle |
| `--replicate-top` | **~5** | DECIDED — refold the top ~5 every cycle; more is acceptable where it buys something, but **must not approach doubling the cycle**, which is a budget constraint |
| replicate anchor | one design with a **known LOW** prior value | DECIDED. Without an anchor a repeat round can only confirm champions; with one it measures regression to the mean |
| accepted-set test | beat the reference on **BOTH** `label` **and** `fef` | MEASURED necessity — `label` saturates and cannot discriminate; `fef` alone is noisier at 10 draws |
| one substitution per position | enforced | two accepted mutations at the same position are mutually exclusive |
| enumeration order | lowest untried order first, then lexicographic by position tuple | DECIDED for determinism — a cycle is then reproducible from its training table alone |
| enumeration size, as run | **88** = 50 (orders 2–3) + 38 (orders 4–6); **no order 7** | MEASURED. Exhausted |

---

## T14. Surrogate gate

**Value:** out-of-bag `RMSE / sd(y) + 1 s.e. < 1.0`.

**Provenance: DECIDED** on the principle that a ratio ≥ 1.0 is at or worse than
predicting the mean — so the gate value is not a tuning choice, it is the
definition of "has any skill at all".

⚠️ **The `+ 1 s.e.` term is a later correction**, after a ratio of
**0.999** printed as "1.00" and was read as a pass. The verdict flipped on the
third decimal.

**Measured values** are in `OBJECTIVES.md` §6. Both arms fail.

---

## T15. Backbone pass-rate assertion

| constant | value | provenance |
|---|---|---|
| `EXPECT_DRAWS` | **300** | arithmetic for that campaign: 100 sequences × 3 draws per backbone. ⛔ Re-derive it for yours — the point is that the total is **asserted**, not the number |
| cut on | **medians per backbone**, never counts | MEASURED necessity — a truncated or uneven draw set reads too high and has already reversed a winner on this project |
| thresholds reported | **several levels, not one** | DECIDED, stated in the script: "the threshold is a choice and should be made in the open" |

⛔ Backbones **span shards** (the fold set is length-sorted, shards are
sequential chunks), so a backbone's draws can split across two shards. Aggregate
across every shard and assert the total.

---

## T16. Two-site (bridge) constants — recorded for completeness, NOT site A

| constant | value | provenance |
|---|---|---|
| objective | `min(n_A/11, n_B7p/7)` if the draw folds, else 0, meaned over **every** draw | DECIDED, and the superseded proposal's three defects are MEASURED (`OBJECTIVES.md` §5b) |
| draws | **3** | MEASURED for a 130–150-residue binder. ⛔ Does not transfer to site A's 10, and site A's 10 does not transfer here |
| sphere radius | **14 Å** to the **union** of both footprints ⇒ 224 residues / 18 segments | DECIDED against the 384 budget: 224 + 150 = 374 |
| binder length | **130–150** | DECIDED. ⚠️ The single-site measurement "extra length buys helices, not interface" **does not apply** — a bridge must physically cover a measured **45.59 Å** max span, so length is load-bearing here |
| geometry, MEASURED | closest approach between sites **11.46 Å**; footprint centroid separation **25.25 Å**; max span **45.59 Å** | MEASURED on the tethered model before anything was built |
| zero-inflation, MEASURED | **21 of 50** sequences at exactly 0 on the chosen backbone; 7 of 50 on another | MEASURED ⇒ check it every cycle before trusting any model ranking |

---

## T17. Constants that are NOT thresholds but will bite

- **Chain convention:** binder = **chain A**, target = **chain B** in the
  generator's output. ⛔ **Opposite to the previous generator version.** Chain
  selection has silently changed a conclusion on this project more than once —
  name the chain every time.
- **Binder length lives in the contig string's leading integer**, not in a
  `length` field. The spec forbids extra keys, so a `length` key is either
  silently ignored or fails at load.
- **The WT→output map is piecewise over 15 segments.** No constant offset exists.
  MEASURED: WT H409 is output B141; the old offset formula gives B40, which is a
  **cysteine**. One site failed loudly on a type assertion; another's offset
  "worked" by luck because all of its footprint residues fell in the first
  segment. Verify the map by **residue type** (231/231 and 230/230), not by count
  — an anchor-only check cannot catch a wrong map.
- **The prediction input's alignment path must be absolute.** A relative path
  makes the predictor **skip the design and exit 0**.
- **A `templates:` block must be present** or the forced-template flag is
  silently ignored — the runner enables potentials by grepping for the block.
  MEASURED: target fit **17.5 Å versus 0.69 Å**.
- **The sequence fasta's first record is the INPUT sequence**, not a design — it
  has no `id=` field. Counting it inflates every statistic. MEASURED: the
  submitted backbone's fasta holds **101 records = 1 input + 100 designs**.
- **A merged public table's `length` column was the FULL construct** — binder
  plus the constant 621-residue target, verified for all 400 rows. Rank orderings
  are safe (Spearman 1.00); every absolute value is **+621 too large**. The median
  binder is **88.5** residues, not 709.
