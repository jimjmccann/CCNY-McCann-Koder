# Conditional TNF-α binder — method, in reader order

James McCann and the Koder Group, The City College of New York.
Campaign `v46`, design work **2026-10-17 .. 2026-10-18**.

⚠️ **Read `../README.md` first.** It states the limitations, and on this challenge the limitations
*are* the result. This file documents how the numbers were produced, not that they are good.

**Read these three files:**

| file | what it answers |
|---|---|
| `README.md` (this) | what the pipeline is, the run order, every threshold, the honest result |
| `SCRIPTS.md` | per script: what it does, what it expects, what it refuses |
| `INPUTS.md` | required inputs, and which are not distributed |

Statements below are marked **MEASURED** (checked against the deposited files when this was
written) or **READ** (quoted from a dated internal record). Where a number could not be traced to
a deposited file, this document says so rather than supplying a rationale.

---

## 0. Scope, and what is NOT here

This section documents the **C1** (single-binder) arm against the human TNF-α trimer, campaign
`v46`. ⭐ A second, **C3-symmetric** arm ran in parallel and **is now deposited** — 12 `v50`
homotrimeric-binder backbones, their sequences, and their scoring tables are in
`../data/submission_20261019/`, and **3 of the 20 submitted sequences come from it**. See
`C3_ARM.md` for the symmetric method, which differs from everything below in three ways that
matter: a replicated motif, a multi-copy chain layout, and **conjunction** aggregation over the
three copies rather than a mean.

⛔ **The launch and compute-orchestration layer is deliberately absent**, as in challenge 1. It
carries account identifiers and cost figures, which are not published. §3 says exactly what has to
happen on the compute side, without the machinery.

⛔ **Two of the challenge's three requirements are not addressed.** The challenge asks for a binder
that holds human TNF-α at neutral pH and releases at acidic pH (ranked first), that also binds
mouse TNF-α (ranked second), and that binds at all (ranked third). This cohort carries **one**
pinned histidine as the pH trigger; the group's own pre-registered objective for this challenge
(2026-10-07) specified a per-draw conjunction over **three** histidines with a geometric switch
criterion. **The switch was never measured on these designs, and mouse cross-reactivity is not
addressed in this cohort.** That is a descope and it is stated rather than smoothed over.

---

## 1. The target

Human TNF-α, soluble trimer, UniProt **P01375-1**, target residues **77–233 (157 aa)**, structure
reference **1TNF** chains A–C. The recommended epitope is the receptor-binding site between two
adjacent protomers. Both sequences and the epitope recommendation are the organizers' own, read off
the challenge page rather than chosen by us.

⚠️ **1TNF carries a defect inside the epitope** and was repaired before use: it has **LEU at
residue 143 where the challenge's own sequence has ASP**, in all three chains. The repair and its
verification are documented in the group's challenge-2 input record; ⛔ **that record is not part of
this deposit**, and the repaired target coordinates are not deposited either (`INPUTS.md`). The
deposited design files contain the target as it was used, so the as-used target is recoverable from
them.

**Generation addresses the target as residues 6–157 of each chain** — a 152-residue span per
protomer, three protomers, 456 residues of target held rigid.

---

## 2. Stage 1 — backbone generation (RFdiffusion3), conditioned on GUIDEPOSTS

⭐ **This is the method that is new this round.** Challenge 1 conditioned the generator on
**hotspots** — "touch these target atoms" — and let the sequence model choose everything else.
Here the generator is instead given **our own placed residues** and builds a chain that contains
them: three arginines positioned to bridge **Glu107** of a protomer, and one histidine as the pH
trigger. The generator holds those residues as a rigid template and designs around them.

⛔ **Guideposts are recorded ONLY in a non-standard column of the generator's own mmCIF output** —
`_atom_site.is_motif_atom_with_fixed_seq`. There is no REMARK, HETATM or segID record of them.
Every downstream consumer therefore **reads the flags from that column and never infers them**, and
**the deposited design files retain the column**, so the pins are recoverable from this deposit
alone. MEASURED: every design carries **4 pins = 3×ARG + 1×HIS**, at **per-design positions** (for
example `ARG32, HIS75, ARG119, ARG169` on one backbone and `ARG36, HIS76, ARG120, ARG166` on
another) — ⛔ so a hardcoded pin index is always wrong.

**The campaign as built and as run** (READ from `scripts/j484_make_c1_core_yield_bundle_20261017.py`):

| | |
|---|---|
| contig | `25-35,E1,39-48,D1,39-47,E3,43-56,E2,20-30,/0,A6-157,/0,B6-157,/0,C6-157` |
| binder length band | **170–220** residues |
| pins | **4** (`E1,D1,E3,E2`), interleaved with flexible spans |
| orientation token | `[0.927, 0.807, 40.226]` |
| arms | **b1** core N=8 "compact" · **b2** core N=8 "near_hyd" · **b3** core N=12 "compact" |
| recycling | `n_recycle = 3` on every arm |
| shards per arm | **8**, 16 designs each ⇒ **128 per arm, 384 total** |

⛔ **The chain-break token `/0` must be preceded by a comma.** Written without one it registers no
break at all and silently changes the output layout — three target chains become one.

⛔ **The three arms are not a sweep over one variable**; b1 and b2 are two different 8-residue core
selections that measured as tied, and b3 is a free hypothesis riding in the same run at no extra
wall-time cost. ⇒ a difference between b1 and b2 is not attributable to core size.

### 2a. Four things about this generator that cost this project time

- ⛔ **The binder length lives in the contig string's spans**, not in a `length` field; the spec
  forbids extra keys, so a `length` key is either ignored or fails at load.
- ⛔ **It can be told what to TOUCH; it cannot be told what to AVOID.** There is no anti-hotspot or
  exclusion conditioning. Everything expressed as "avoid" therefore lives in the filter ⇒ **the
  filter is where the design happens; it is not QC.**
- ⛔ **Runs are unseeded and not reproducible.** MEASURED on an earlier cohort: 16 of 16 outputs
  distinct from the same input. ⇒ the artefact is the record, and re-running this bundle will not
  reproduce these 384 backbones.
- ⛔ **Designs are clustered by batch** — one contig length per batch — so the effective n for any
  length-dependent claim is **batches, not designs**.

---

## 3. What has to happen on the compute side, stated without the machinery

1. **One unit of work per worker, one output prefix per worker.** Two workers sharing an output
   prefix means the first to finish writes a complete-looking manifest and a verifier can bless the
   wrong one.
2. ⛔⛔ **Verify the output COUNT per unit of work. Never the exit code.** Both the generator and the
   structure predictor used here can catch their own out-of-memory condition, write no structure and
   **return success**. rc=0 with zero structures is a normal, silent outcome.
3. **Count inside the archives.** Outputs are packed per run, so searching for loose structure files
   returns 0 on a perfectly good round, and on a genuinely failed run the working directory comes
   back fully populated — so a file-existence check reports success as well. **Only the output count
   against the input count catches it.**
4. **Keep the log in the archive**, not just the outputs. A failed shard that kept only a fatal-error
   file and no log has a cause that will never be proven.
5. ⛔ **A pre-launch gate that runs before the thing it checks exists is not a gate.** One bundle in
   this project passed its own canary because the canary ran before the bundle was assembled; 12 of
   12 shards then died on arrival. A check that cannot fire reads exactly like a check that passed.
6. **Ephemeral workers do not stop themselves.** On finish the container re-runs its start command,
   fails, and restart-loops at full rate. Termination is part of the cost model.

---

## 4. Stage 2 — the two BACKBONE gates, and they run before anything else

⛔⛔ **Continuity is gated FIRST, before any clash or interface metric.** This is the most expensive
lesson in the TNF-α campaign. An earlier cohort was ranked on a backbone-clash floor; **69% of it
carried a chain break over 10 Å**, and the clash floor *rewarded* the break — a chain with half of
itself 56 Å away has fewer close contacts, so 4 of the top 5 by clash score were disconnected
chains. A clash metric applied to a broken chain is measuring the wrong object.

**Gate 1 — continuity** (`j393_continuity_gate_20261015.py`). Maximum Cα–Cα step along the chain.
Ideal trans Cα–Cα is 3.80–3.81 Å; the bar is **`max_step ≤ 4.2 Å`**, with a second
**catastrophic** bar at 10.0 Å reported separately.
⛔ Numbering gaps are **not** breaks: steps are measured only between residues whose numbers are
consecutive, and a numbering jump is reported in its own column.
⛔ The primary column is the **continuous** quantity `max_step`, not the pass/fail flag — a binary
gate at small n returns all-zeros and ranks nothing.
⛔ Chain roles are **detected**, never named, from the fraction of residues carrying the
guidepost flag (~100% for a held target, a few percent for a designed binder). A hardcoded chain
pair once compared a designed chain against the target and reported it as a self-interface, with no
error.

**Gate 2 — topology** (`j409_path_or_fold_20261015.py`). Is the chain **one fold**, or a **path**
from pin to pin? `min_cut` is the fewest long-range Cα–Cα contacts (< 8 Å, |i−j| > 4) crossing any
interior cut of the chain, with a **margin of 13** residues at each end. A single folded domain
keeps many contacts crossing everywhere; a path, or two lobes on a tether, has a point almost
nothing crosses. The bar is **`min_cut ≥ 19`**.
⛔ `min_cut` is margin-sensitive, so the reference values are recomputed at the same margin as the
designs and never quoted across settings.
Also reported, because `min_cut` alone cannot separate "compact fold" from "compact blob":
long-range contact fraction (|i−j| > 12), and Rg against the globular expectation 2.2·N^0.38 Å.

⛔ **Continuity and topology are ORTHOGONAL** — measured on an earlier cohort, overlap zero. Ranking
on either alone ranks the wrong axis. So the gate is the **conjunction**, and it is called JOINT.

**MEASURED, reproduced for this document from `../data/backbone_gates_20261018/`:**

| gate | count of 384 |
|---|---|
| continuous (`max_step ≤ 4.2 Å`) | **23** — by arm b1 15 · b2 7 · b3 1 |
| topology (`min_cut ≥ 19`) | **191** (cohort `min_cut` median 18.0) |
| ⭐ **JOINT (both)** | **11** — by arm **b1 5 · b2 6 · b3 0** |
| cohort `rg_ratio` median | 0.969 |

⛔ **Quote a cohort's own distribution, never another campaign's.** An earlier build of the figure
script printed a different cohort's `rg_ratio` median (1.26) beside these 11, which would have read
as these designs beating their cohort by a wide margin. `v46`'s own median is 0.969 and the 11 JOINT
sit at 0.944 ⇒ **`min_cut` is what separates them, not compactness.**

⛔ **The 11 JOINT are the top 2.9% of the run and are not a typical design.**

---

## 5. Stage 3 — sequence design (SolubleMPNN)

Model: **SolubleMPNN**, checkpoint `solublempnn_v_48_020.pt`, run through the LigandMPNN codebase.
⭐ On the two names: **SolubleMPNN is the model** (weights trained to avoid the surface hydrophobics
a membrane-protein-containing training set otherwise encourages); **LigandMPNN is the codebase**
those weights run through. A command line that invokes LigandMPNN while loading a SolubleMPNN
checkpoint is doing one thing, not two.

Settings, MEASURED from the deposited FASTA headers in `../data/mpnn_20261018/seqs/`:
`T = 0.1`, `seed = 37`, `batch_size = 10`, `number_of_batches = 1`,
`use_ligand_context = True`, `ligand_cutoff_distance = 8.0` ⇒ **10 sequences per backbone, 110
total**, one FASTA per backbone.

| flag | setting | why |
|---|---|---|
| `--omit_AA C` | cysteine **omitted** | ⛔⛔ non-negotiable: free thiols are an oxidation and aggregation liability and cysteines break the assay's expression system. MEASURED: **zero cysteines** in all 110 designed chains |
| proline | **permitted** | proline is the helix-capping residue; omitting it forces caps out of residues that cap badly |
| `--chains_to_design A` | target is **context** | the target is kept in the input and never designed, so the model sees the interface it designs against |
| `--fixed_residues_multi` | the 4 pins, **per design** | ⛔⛔ without it the sequence model silently overwrites the entire Arg/His motif |

⛔ **The sequence model will happily redesign a deliberately placed residue, and nothing in the
output marks one.** Pin positions are read **per design** from the guidepost flag column — never
hardcoded, because the generator chooses the sequence positions itself. A zero pin count refuses.
The per-design pin record is deposited as `../data/mpnn_20261018/fixed_residues_multi.json`.

⛔ **Check the cysteine count on the DESIGNED chain's field, not field 1.** The FASTA holds **all**
chains, colon-separated.

⚠️ **The first record of each FASTA is the INPUT sequence, not a design** — it has no `id=` field.
Counting it inflates every statistic. MEASURED across the deposit: **121 records = 11 inputs + 110
designs.** ⛔ And `cat *.fa | grep -c '^>'` reads **111**, not 121, because these files carry no
trailing newline so each of the 10 file joins swallows one `>`. **Count per file, never across a
`cat`.**

---

## 6. Stage 4 — fold prediction (Boltz-2) against the WHOLE trimer

One prediction input per **sequence**: chain **A** = the designed binder, chains **B, C, D** = the
three target protomers, 152 residues each.

⛔ **Chain naming matches the design file** (binder A, target B/C/D) so the self-consistency RMSD
needs no remapping.
⭐ **The target sequence is not retyped.** It is read from the design file's own target chain and
asserted byte-identical to the 152-residue chain used by an earlier, proven set of inputs — two
independent sources agreeing is what settles identity.
**MSA convention:** the binder gets an **empty** MSA (it is de novo; an MSA is meaningless and
sequence models do not transfer to de novo backbones); the target gets the real TNF-α alignment,
**2198 rows**.

Prediction settings, READ from `scripts/j496_run_stage2_after_stage1.sh`:
`--diffusion_samples 1 --recycling_steps 3 --seed 20261018 --no_kernels --output_format pdb
--write_full_pae`.

⛔ **One sequence is one draw here.** `--diffusion_samples 1` means there is **no within-sequence
replicate**, so no per-design error bar exists from this run. That is a known weakness: in this
group's work, re-folding measures about **1.56×** the noise that within-batch standard error
suggests, so the spread quoted anywhere in this deposit understates the true spread.

⛔ **Folds were chained, not parallelised** — each holds ~9.4 GiB on a 16 GiB card, and an
out-of-memory exit returns 0 with no structure written. Every stage therefore verifies an **output
count**, never an exit status.

⛔ **Probe one input before a batch.** A previous batch on this project hung for 8.5 hours at 0% GPU
utilisation. The deposited run is `probe (1) → stage 1 (10) → stage 2 (99)`, and the probe's own
scored row is deposited separately.

---

## 7. Stage 5 — scoring. ONE program, and what it refuses to do

`j497_c1_boltz_selfconsistency_20261018.py`, **103 columns**, one row per prediction.

⭐ **The shared geometry is IMPORTED, not copied.** `kabsch_rmsd` and `tm_score` come from the
sibling C3 scorer `j465_boltz_selfconsistency_20261017.py` **under a SHA1 pin, checked every run**.
The C1 and C3 numbers get compared to each other, so a silent copy that drifts is the failure mode:
if the sibling changes, the scorer **refuses** rather than scoring with different math.
⚠️ **The pin was re-pointed for this deposit** and the reason is recorded in the file: one absolute
path in the sibling became a placeholder, which changed its hash while leaving the two functions
byte-identical. The hash the published numbers were produced under is named in the comment.

**Target residues are addressed by SEQUENTIAL INDEX, never by residue number.** The predictor
renumbers the target 1..152 while the design may carry TNF-α numbering. The offset was verified at
**three** points on real output — K11→6, E107→102, E110→105 ⇒ `predicted = TNF-α − 5` — and the
scorer **asserts the residue it lands on is actually GLU** before measuring from it. A silent miss
would turn the bridge column into a distance to an arbitrary residue and read exactly like a result.
MEASURED independently for this document, straight out of the deposited design coordinates:
sequential index 102 is **GLU in 33 of 33** target chains (11 designs × 3 protomers).

**The bridge metric.** Minimum distance from a binder arginine's guanidinium nitrogen
(`NE`/`NH1`/`NH2`) to that target chain's Glu107 carboxylate oxygen (`OE1`/`OE2`), **heavy atoms,
never Cα** — the pins are whole residues and every atom is held, and Cα-gating has passed a 1.67 Å
clash on this project. Best of the three protomers. Reported **twice**: once over *any* arginine and
once restricted to the *pinned* arginine. ⛔ The pinned column is the honest one; see `../README.md`.

**Two fidelity numbers, kept separate:**
- `rmsd_binder` / `tm_binder` — binder superposed on itself ⇒ **fold** fidelity, independent of
  placement;
- `rmsd_binder_on_target_frame` — superpose on the **target only**, then measure the binder with
  **no** re-fit ⇒ **placement** fidelity.

Read together: good fold + large placement = folded right, docked wrong, which is a real category
here and not a near miss.

**Two confounds are excluded in-script rather than argued away:**
1. ⛔ **The target is a C3 trimer, so there are three symmetry-equivalent sites.** Under the identity
   chain mapping, a binder at an equivalent site reads as a ~40 Å displacement while being the same
   mode. The scorer tries **all three cyclic mappings** and keeps the best. MEASURED over the first
   round (n=11): identity mapping median **46.08 Å** → best-of-three median **42.44 Å**, and the
   winning frame is split `CDB` 5× / `BCD` 4× / `DBC` 2×. ⇒ it barely helps, and the displacement is
   **not** a symmetry artefact.
2. ⛔ **A bad frame would fake the entire number.** The **target-only RMSD is reported beside every
   placement**, and a row whose target frame superposes badly is flagged `FRAME_SUSPECT` and its
   placement columns **must not be read**. MEASURED over the first round: trustworthy frames
   superpose at **1.48–1.72 Å** (n=9) against suspect rows at **21.83–21.93 Å** (n=2), meaning the
   predictor built the trimer *assembly* differently for those. Over all 110 scored folds, **37
   (34%) are flagged**, and the flag is **per-sequence, not per-backbone**.

⛔⛔ **Per-pair confidence, never the global number.** With 1 binder + 3 protomers there are 6 chain
pairs and **3 are target–target**. `iptm`, `plddt` and `pae` are therefore all reported **per chain
pair and per chain**, never pooled. See `../README.md` §6 — a global `iptm` on this construct is a
measurement of TNF-α, not of our binder.

⛔ `n_target_chains_engaged` is **2 of 3** for most designs. For a single-chain binder at an
inter-protomer site that is expected, not a defect.

---

## 8. The honest result

- **A guidepost route works mechanically.** The generator accepts our own placed residues, holds
  them rigidly, builds a chain around them, and the pins survive into the output where a scorer can
  find them without being told where they are.
- **The backbone gates are strict and they bite.** 23 of 384 continuous, 11 of 384 JOINT — **2.9%**.
  Continuity and topology are orthogonal, so both are needed and neither is sufficient.
- ⛔⛔ **The designed interaction is not reproduced by the fold-back check in 0 of 110 cases**, and
  the binder lands on a different face of the target. `../README.md` has the numbers.
- ⭐ **A positive control was run and it passes.** A published de novo binder with measured affinity
  reproduces its own crystal pose to **0.89 Å** through this identical path, and a negative control
  fails. ⇒ the fold-back check has demonstrated sensitivity on a de novo interface, so "the
  instrument is blind" is no longer an available explanation for the 0 of 110. ⛔ That still does
  not make "the designs are wrong" a measurement — only a bench assay does.
- ⛔ **Backbone identity predicts only how well the predictor reproduces that backbone's own fold,
  and predicts nothing about binding.** Within-backbone spread across sequences is about twice the
  between-backbone spread on the binding metrics.
- ⭐ **The cohort is complete and balanced: 110 folds, exactly 10 per backbone across 11.**

⇒ **This is a working guidepost-conditioned generation and gating pipeline with a measured,
strict backbone yield, and an unconfirmed interface. It is not a binder. It was submitted anyway,
as 20 sequences across two arms, because an untested design that is documented is worth more to
this competition than a design we never entered.**

---

## 9. Provenance of this documentation

Written **2026-10-18**. Every count in §4, §5 and §7 was recomputed from the files deposited in
`../data/` when this was written; where a figure comes from a dated internal record that is not
deposited, it is marked READ. Internal note and decision identifiers that remain in the deposited
scripts are explained at the end of `../README.md`.
