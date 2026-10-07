# De novo EGFR binder — method, in reader order

James McCann and the Koder Group, The City College of New York.
Site-A binder campaign, design work **2026-09-28 .. 2026-10-06** (see the date note in
`../README.md`).

This directory holds the documentation needed to reproduce the **design criteria** —
what we optimised, what the numbers mean, where every threshold came from, and in
what order the steps run. It is written for someone outside the lab with the
scripts in hand and none of our context.

**Read these five files:**

| file | what it answers |
|---|---|
| `README.md` (this) | what the pipeline is, the run order, the scope and the honest result |
| `SCRIPTS.md` | per script: what it does, what it expects, what it refuses |
| `OBJECTIVES.md` | every metric defined, and why it is what it is |
| `THRESHOLDS.md` | every constant and cutoff, traced to its source |
| `INPUTS.md` | required inputs, which are site-specific, which are not distributed |

Two further files are the record of what we got wrong:
`ideas_that_didnt_work.txt` and `design_guidelines.txt`.

⚠️ **Nothing in this campaign has been tested experimentally.** Every number is
computational, on one target, with one generator and one structure predictor.

---

## 0. Scope, and what is NOT here

This documents **site A** — the domain-III face of the EGFR extracellular region.
Where a number is site-specific it is marked.

⚠️ **The two-site "bridge" arm IS part of the submitted path.** It was settled on
2026-10-06 and its sequence `id85` is **rank 1** of the deposited set.
The bridge's own constants are in `THRESHOLDS.md` T16 and its objective in
`OBJECTIVES.md`; the decision and the evidence are in `DESIGN_CHOICES.md` §3. The
**second site** remains parked and contributes no deposited sequence.
⇒ Read this file as the **site-A** method, and `DESIGN_CHOICES.md` for the arm
structure and which molecules were deposited from which arm.

⛔ **The launch and compute-orchestration layer is deliberately absent.** It
carries account plumbing and is not part of the method. Any structure predictor
and any batch system reproduces this pipeline; §5 says exactly what has to happen
on the compute side and why, without the machinery. That layer carries account
identifiers and cost figures, which are not published.

---

## 1. The target

Human EGFR extracellular region. The assay construct is **binder + a constant
621-residue target**, read out of the competition's own released data rather than
chosen by us — so no design is forced to dock at a face we picked, and a design's
apparent epitope is a measurement.

**Numbering is WT throughout** (mature EGFR = UniProt P00533 minus 24). Older lab
documents carry precursor numbers; those are converted, not reinterpreted.

The **site-A footprint is 11 residues**:

    380 382 384 408 409 410 411 412 417 438 465

Asserted by residue identity, not by number range, everywhere it is used
(`THRESHOLDS.md` T1). Among them, **F380 is the rate-limiting contact** — see
`OBJECTIVES.md` §3.

---

## 2. Stage 1 — target preparation: a SPHERE, with an exclusion gate

Generation needs a smaller target than 621 residues.

⛔ **Do not cut a contiguous window.** Our first round did and was **voided**:
70% of site-A and 100% of site-B designs sat within 2 Å of receptor that had been
deleted — backbone overlap, not contact. The window had been validated by
*inclusion* ("all 11 footprint residues present") and never by *exclusion*, and
every downstream stage then treated the window AS the target, so the omission was
unfalsifiable from inside the pipeline. Full account in
`ideas_that_didnt_work.txt` entry 1.

**What replaced it.** A **22 Å sphere** around the site-A footprint: 231
residues in **15 segments**. Contiguity was never a requirement — measured, the
generator holds a non-contiguous target rigid to 0.0402 Å over 978/978 atoms,
nothing invented in the gaps, no linking.

Two things make the sphere auditable, and both are load-bearing:

1. **A persisted WT→output residue map.** The generator renumbers the target
   1..N contiguous and discards both the original numbers and every segment gap.
   With 15 segments **no constant offset exists**. Three of our scorers had one
   hardcoded; one mapped the key histidine to a cysteine. The map is written at
   cut time, shipped with the target, loaded by every consumer, and verified
   against real generator output **by residue type, 231/231**.
2. **A full-receptor exclusion gate.** Superpose each design back onto the whole
   receptor and measure overlap with receptor the generation target did not
   contain. It **must superpose** — the generator recenters its output by
   66.93 Å — and because the target comes back rigid at ~0.04 Å, the
   superposition RMSD is itself a tripwire.

⛔ The sphere is an improvement, not a fix. Its radius is a guess; the gate is
the only thing that proves the guess was generous enough. Trimming preserves
coordinates, not physics: a residue whose neighbours were cut is artificially
exposed.

Scripts: `setup_rfd3_sphere_campaign.py`, `make_wt2out_map.py`, `wt2out.py`,
`full_receptor_gate.py`.

---

## 3. Stage 2 — backbone generation (RFdiffusion3)

**Site-A sphere campaign, as built and as run** (MEASURED from the shard set on
disk, 15 shards, `arm_shard00..arm_shard14`):

| | |
|---|---|
| target | 22 Å sphere, 231 residues, 15 segments |
| binder lengths | **60, 65, 70, 75, 80** — five values |
| shards per length | **3** (`REPS = 3`) ⇒ 15 shards |
| designs per shard | 336 (`n_batches = 42`) ⇒ **5,040 per arm, 1,008 per length** |
| conditioning | `select_hotspots` (atom-level) on 3 residues, `infer_ori_strategy: "hotspots"`, `is_non_loopy: true` |
| classifier-free guidance | **OFF** |
| hotspots | H409 (ND1, NE2) · F412 (ring carbons) · I438 (CG2, CD1) |

The backbone the whole optimizer walks from is **`arm_shard11_shard11_24_model_2`**.
MEASURED: `arm_shard11`'s contig is `75,...`, so that backbone is **75 residues**,
which matches the optimizer's own `seq_len = 75`.

⚠️ **Correction to an earlier internal inventory**, recorded so it does not
propagate: that inventory described stage 2 as "31 shards, one fixed length each,
70–100 aa". That describes a *different*, earlier shard builder
(`make_rfd3_shards.py`, 31 shards at lengths 70..100) which is **not** the
generator of the submitted backbone. The site-A sphere set is 15 shards at
60–80 by 5. Verified by reading the contigs out of the shard JSONs.

**Four things about the generator that cost us time and are worth stating:**

- ⛔ **The length lives in the contig string's leading integer**, not in a
  `length` field. The spec forbids extra keys, so a `length` key is either
  silently ignored or fails at load. The builder rewrites the contig and refuses
  if it cannot find that integer.
- ⛔ **The generator can be told what to TOUCH. It cannot be told what to AVOID.**
  There is no anti-hotspot / exclusion / avoidance conditioning. The previous
  version exposed one; it is gone. Everything expressed as "avoid" therefore
  lives in the filter ⇒ **the filter is where the design happens; it is not QC.**
- ⛔ **Burial / exposure conditioning is untrained for the protein-only case**
  (~0% training share, and the transform early-returns without a ligand on both
  pipeline branches). "Keep this residue at the rim" is not sayable.
- ⛔ **Three shards per length is only safe because inference is unseeded.**
  Verified in source: the seed defaults to None and is only written into output
  metadata. If it were seeded, three shards at one length would be byte-identical
  and the campaign would pay 3× for one shard.

Sampler knobs we turned on and then **refuted against a null**: classifier-free
guidance 2.76× the null, extra recycling 2.55× null-to-worse. Both OFF in
production. A spread hotspot set is also a real risk — orientation is inferred
from the hotspot *centroid*, and one hotspot 12.2 Å from the centroid of the rest
was reached by 0 of 336 designs (r = +0.967 between distance-from-centroid and
contact rate).

Scripts: `setup_rfd3_sphere_campaign.py` (campaign + arm JSONs),
`prod_rfd3_egfr_multi.sh` (the runner), `select_all31.py`,
`backbone_passrate.py`.

---

## 4. Stage 3 — sequence design (SolubleMPNN)

Model: **SolubleMPNN, checkpoint `solublempnn_v_48_020.pt`**.

⭐ **On the two names, because the deposit uses both and they are not a contradiction:**
**SolubleMPNN is the MODEL** — the weights in `solublempnn_v_48_020.pt`, trained to avoid the
surface hydrophobics that a membrane-protein-containing training set otherwise encourages.
**LigandMPNN is the CODEBASE** those weights are run through. So a command line that invokes
LigandMPNN while loading a SolubleMPNN checkpoint is doing exactly one thing, not two.

Three flags, each closing a measured hole:

| flag | setting | why |
|---|---|---|
| `--omit_AA C` | cysteine **omitted** | an earlier arm of this lab permitted it: 2–3 cysteines per design in a protein whose wild type has zero, never pairing (partners 12–15 Å apart) ⇒ free thiols, an oxidation and aggregation liability with no compensating structure |
| proline | **permitted** | free backbone NH/C=O exist only at helix caps and in loops; mid-helix amides are consumed by the helix's own i→i+4 bonds. Proline is the capping residue, so omitting it forces caps out of residues that cap badly |
| `--bias_AA` on D/E | ⛔ **0 — NOT APPLIED in production** | decided as a direction, then **measured to zero**. See the box below; this is the one flag where the small-batch evidence was wrong |
| `--temperature` | 0.1 | script default |

### ⛔⛔ The acidic bias was DECIDED and then MEASURED TO ZERO. Read this before copying the flag.

The decision said "bias acidic", on an n=8 smoke test in which a charge gate
passed 3 of 8 designs and **all 5 failures were on the acidity bound** — and it
deliberately left the **magnitude** open, to be tuned on a batch with net charge
measured before committing.

**MEASURED at n = 3,320 sequences (332 backbones × 10), the direction was
wrong:**

| | NCPR median | in band [−0.13, −0.02] | too POSITIVE | too NEGATIVE |
|---|---|---|---|---|
| site A | −0.100 | 1133/1460 (77.6%) | 44 | **283** |
| site B | −0.0909 | 1387/1860 (74.6%) | 138 | **335** |

⇒ **The pool already overshoots on acidity ~6:1.** A D/E bias would have pushed
*more* sequences out the **bottom** of the band. ⇒ **`--bias_AA D/E` stays at 0.**

⭐ **Measure-then-commit earned its keep: a bias set from n=8 would have hurt.**
⇒ And the general lesson is not "don't bias" — it is that a decision can be right
about *what to control* and wrong about *which way*, and only the batch
measurement separates those.

⛔ `--chains_to_design A` — the target is kept in the input as **context** and is
never designed, so the model sees the interface it is designing against.

⛔ **The sequence model will happily redesign a deliberately placed residue.**
Nothing in the output marks one. Where a design has a pinned residue, locate it
per structure by its own contact geometry and pass it as fixed — the generator
chooses its sequence position itself (68 distinct positions over a 2–70 window in
an n=400 baseline), so the index cannot be hardcoded. The site-A sphere arms are
hotspots-only and carry no pinned residue, so this runs with `--no-guideposts`.

**Sequence QC runs on the fasta alone**, immediately after:
guidepost identity held · cysteine count == 0 · N-glycosylation sequon
`N-[^P][ST]` · net charge at pH 7.4 inside **−0.13 ≤ NCPR ≤ −0.02**.

⛔ The filter **splits in two**, because its inputs do. Net charge and
surface-patch metrics need a *sequence*, so they cannot run on backbones:

- **backbone level:** hotspot engagement, interface contact/burial against the
  target, secondary-structure composition, radius of gyration for length sanity;
- **sequence level (only after MPNN):** the four checks above;
- **after the fold-back:** anything needing per-atom solvent accessibility.

⚠️ The fasta's **first record is the INPUT sequence**, not a design — it has no
`id=` field. Counting it inflates every statistic. MEASURED: the submitted
backbone's fasta has **101 records = 1 input + 100 designs**.

Scripts: `mpnn_egfr.py`, `seq_qc_egfr.py`, `make_w10_tsv.py`.

---

## 5. Stage 4 — fold prediction (Boltz-2, against the full 621-mer)

One prediction input per **sequence** (not per draw). **10 draws per sequence**
for this 75-residue binder — a measured number, not an inherited one; see
`THRESHOLDS.md` T3.

The construct, which took measurement to get right: chain A = binder, chain B =
the full 621-residue target with its own alignment, **10 NAG glycans**, glycan
bond constraints, and a **forced template**.

⛔ **Do not rebuild that construct.** We did once, lost the script, and rebuilding
it from scratch is the expensive way to make a silent mistake. The builder in use
takes a **real landed input file as a skeleton and changes exactly one line** —
chain A's sequence — so every construct detail carries over byte-for-byte. It
then **asserts** the skeleton rather than trusting it:

- one chain-A sequence line;
- an **absolute** path to the chain-B alignment. A relative path makes the
  predictor **skip the design and exit 0**;
- a `templates:` block. The runner enables potentials by *grepping for it*, and
  without potentials the forced-template flag is **silently ignored** — measured,
  target fit 17.5 Å versus 0.69 Å.

⭐ And **prove** the new inputs differ from a real landed one in exactly the
intended way. Done here for a 50-sequence round: every file differed in exactly
the chain-A sequence line, 50 of 50 — because a rebuilt construct would have
silently broken the pooled comparison.

### What has to happen on the compute side, stated without the machinery

1. **One unit of work per worker, one output prefix per worker.** Two workers
   sharing an output prefix means the first to finish writes a complete-looking
   manifest and a verifier can bless or terminate the wrong one.
2. ⛔⛔ **Verify the output COUNT per unit of work. Never the exit code.** The
   predictor catches its own out-of-memory condition, skips the batch, writes no
   structure and **returns success**. rc=0 with zero structures is a normal,
   silent outcome, and it is non-deterministic — a re-run of the same shard
   reported zero failures. One worker of six surviving is not evidence the
   configuration works.
3. **Count inside the archives, and count — do not test for existence.** Outputs
   are packed per run, so searching for loose structure files returns 0 on a
   perfectly good round. ⛔⛔ And on a genuinely failed run the working directory
   comes back **fully populated**, so a file-existence check reports success as
   well. MEASURED on a deliberate probe: the heaviest glycan mode OOMed at rc=0
   with **zero** structures written while the directory looked complete. **Only
   the output count against the input count catches it.**
4. **Keep the log in the archive**, not just the outputs. Five of six failed
   shards in one of our rounds kept a fatal-error file and no log, so their cause
   is not proven and never will be.
5. **A tool's own error text can be wrong.** One fatal-error file blamed a missing
   alignment; the alignment was present and the cause was memory. Separately,
   "empty MSA" lines referring to the *binder* chain are expected for a de novo
   design and are not a fault.
6. **Ephemeral workers do not stop themselves.** On finish the container re-runs its
   start command, fails, and restart-loops at full rate. Termination is part of
   the cost model — a cost figure that assumes workers stop is wrong.

Scripts: `build_ecto621_yamls.py` (the canonical builder, where the construct was
first got right), `build_yamls_from_skeleton.py` (the one in the loop),
`build_ecto621_bundle.py`, `fold_shard_pod.sh`.

---

## 6. Stage 5 — scoring. THREE PROGRAMS, IN THIS ORDER, AND THE ORDER MATTERS

    extract the archives
      -> contacts_a21.py              (epitope contact, per draw)
      -> selfconsistency_ecto621.py   (fold + placement fidelity, per draw)
      -> score_cycle_egfrA1.py        (pool draws to one row per sequence)
      [-> score_0047.py               (the current objective, separate pass)]

⛔⛔ **Gate on fold quality BEFORE ranking on contacts.** This is the most
expensive lesson in the project. A contact-only backbone ranking named six
winners and **four of the six failed the fold gate**; the top pick had a median
Cα fold error of **8.45 Å**, so every one of its 11 "full-footprint" draws was
contact against a different structure. Arm-wide the gate removes **53% of
full-footprint contact, 422 draws → 200**, and 16 of the 31 failing backbones had
≥5 full-footprint draws — the contact metric was *actively promoting* unbuildable
backbones.

**Two fidelity numbers, kept separate:**

- `sc_rmsd` — binder-onto-binder Cα RMSD (superposed on the binder alone) ⇒
  **fold** fidelity, independent of docking;
- `pose_rmsd` — superpose on the mapped target, then binder Cα RMSD with **no**
  binder re-alignment ⇒ **placement** fidelity. Includes fold error, so
  `pose_rmsd ≥ sc_rmsd` always.

Read them together: small `sc` + large `pose` = folded right, docked wrong.
Large `sc` = it did not build the designed shape and nothing else means anything.

⛔ **Cut on medians per group, never on counts.** A truncated or uneven draw set
reads too high, and this has already reversed a winner on this project once.
Groups span shards, so pool every draw and **assert the expected total**.

⛔ `pose_rmsd` is reported as a diagnostic and is **never a cut** (an 8× measured
swing); epitope contact is the primary quantity and the fold gate is the filter.

Metric definitions: `OBJECTIVES.md`. Thresholds: `THRESHOLDS.md`.

---

## 7. Stage 6 — proposing the next cycle

1. **Find the reference by a rule, not by judgement.** The reference is the
   **hub**: the training row with the most Hamming-1 neighbours, because the
   accepted set is only meaningful relative to the sequence whose single mutants
   were actually measured. MEASURED on the cycle-2 table: 50 neighbours against 4
   for the runner-up. With a thin margin the proposer **refuses** — a wrong
   reference silently redefines every mutation in the cycle.
2. **Apply the fold floor to the reference too.** One design once topped the table
   at 10.333 on 3 of 10 passing draws with a median fold error of 1.871 Å.
3. **Recompute the accepted set**: mutations that beat the reference on **both**
   ranking metrics. `label` alone can no longer discriminate (it saturates) and
   the fraction alone is noisier at 10 draws.
4. **Enumerate untried combinations**, lowest untried order first then
   lexicographic by position tuple — so a cycle is reproducible from its training
   table alone. One substitution per position (two accepted mutations at the same
   position are mutually exclusive).
5. **Reserve ~5 slots for re-folds** of the previous cycle's best (§8).
6. **Never re-fold anything already in the table**, except as a declared repeat.

⛔⛔ **There is no model in this loop, and that is MEASURED, not a preference.**
Out-of-bag RMSE/sd(y) was **1.14** and **1.04** on the two labels — at or worse
than predicting the mean — and UCB **did not beat random** (+0.394 ± 0.634,
t=+0.62, Fisher p=1.000). Trained on the saturated label it predicted 11.696 and
12.101 against a physical ceiling of 11. **Do not re-introduce acquisition-based
proposal without re-measuring OOB first**, and quote the ratio with its standard
error: a ratio of 0.999 once printed as "1.00" and was read as a pass. Require
**ratio + 1 s.e. < 1.0**.

Exhaustive coverage needs no ranking, so dropping the model cost nothing.

⚠️ **Enumeration ends.** It completed here at **88 combinations** (50 at orders
2–3 plus 38 at orders 4–6; no order 7 exists), and the last enumerated cycle beat
nothing — "cycle sequences beating the previous best: 0". After that the proposer
has no derivable content and the next cycle is a science decision.

Script: `propose_cycle_egfrA1.py` (and `build_combo_cycle_egfrA1.py`,
`build_sitescan_cycle_egfrA1.py`, `stage_egfrA1_cycle0.py`,
`build_burnin_egfrA1.py` for the hand-built cycles).

---

## 8. Replication — standing, not optional

⛔ **This campaign ran four cycles before anything was ever re-measured.** Every
design had exactly **one** measurement, so there was no noise estimate at all.
The two rows that looked like replicates are **byte-identical across three
successive training tables** — carried forward, never re-measured.

**The policy now:** every cycle refolds its **top ~5** designs with fresh draws,
**riding inside the cycle's own slots** — zero extra workers, and it controls for
batch effects because repeats and new designs are predicted in the same run.
More than five where it buys something; nowhere near doubling the cycle.

⭐ **Include an anchor whose prior value is known and LOW.** Without one, a repeat
round can only confirm champions; with one it measures regression to the mean,
which is the quantity that says how much of each cycle's gain was *selection*.

MEASURED, first in-cycle run, n=6 spanning the score range:

    mean delta -0.058    mean |delta| 0.158    max |delta| 0.300
    Pearson(prior, fresh) +0.556    slope +0.343

Both champions at 1.000 returned 0.700 and 0.800; the low anchor came **up**,
0.400 → 0.600. **Slope +0.343 is regression to the mean ⇒ part of every cycle's
apparent gain was selection, not improvement.**

⛔ **Repeats are NOT new data points.** Never merge one into the training table:
it puts one sequence in twice, double-weights it in anything fitted on the table,
corrupts hub detection, and can let a repeat become the next cycle's parent. They
are scored, compared against their prior values, and written to a separate file.
⛔ Identify them from their **declared file only, never by id prefix** — "id9" as
a prefix already matches ten genuine designs.

---

## 9. The honest result

- The **aggregate** progression is real: between-design variance **3.02×** and
  **2.42×** the binomial null (χ² p=2.5e-11, p=1.3e-06), and the population
  per-draw full-footprint rate went 24.8% → 53.1% (p=5.7e-18) between two cycles.
- **No individual design is resolvable.** Both ranking metrics saturated at their
  ceilings; two designs sit at 1.000 with 10/10 fold passes and are not
  separable. A single design's fraction carries s.e. ≈ 0.16 at 10 draws.
- **Part of the apparent gain is selection**, measured: slope +0.343 (§8).
- **Draws cannot be raised to fix it.** 10 draws already measured **97.6% GPU
  memory** on a 24 GB card for this construct, and out-of-memory exits rc=0.
  Precision comes from repeated 10-draw runs.
- The rate-limiting contact is **architectural**, not a property of the
  substitutions explored, and the accepted set the optimizer walked contains **no
  position that reaches it** — which is why three cycles of enumeration beat
  nothing. See `OBJECTIVES.md` §3.
- A saturation scan at the one position that *governs* that contact found **no
  substitution above the noise floor** (best 0.609 against its own parent's
  0.600, floor 0.158).

⇒ **This is a working closed-loop optimizer with a measured progression and a
measured noise floor. It is not a finished binder.**

---

## 10. Provenance of this documentation

Written 2026-10-05. Every threshold is traced in `THRESHOLDS.md`; where a number
could not be traced, this documentation says so explicitly rather than supplying
a rationale. Statements are marked MEASURED (checked against the
filesystem when this documentation was written), READ (quoted from a dated lab
document or decision record) or INFERRED. Decisions cite the lab's decision log by number.

⚠️ **About those `decisions/NNNN` citations, which you cannot follow.** Several deposited scripts
cite a decision by number — `decisions/0047`, `decisions/0035`, `decisions/0045` and so on. That
is our internal decision log, an append-only record where each entry states one settled question,
the measurement behind it, and the date; it is not deposited. The numbers are kept rather than
stripped because they are the **provenance of a threshold**, and a threshold with its provenance
removed reads as an arbitrary constant. ⇒ Read `decisions/NNNN` as "this value was settled
deliberately, on a stated date, and the record exists" — and read `THRESHOLDS.md`, which traces
every constant that governs a result and says so explicitly where a number could not be traced.
⛔ **Corrected:** an earlier draft of this paragraph said the same applied to "the handful of
`[[slug]]`-style references". **No such references remain in this deposit** — they were an internal
note index and were removed during the public-release pass. This sentence is kept only so a reader
who saw the earlier text knows they are gone, not missing.
