# SCRIPTS — what each one does, what it expects, what it refuses

Site-A EGFR binder campaign. Companion to `README.md`.

**Read this with `README.md` §6 open: the ORDER matters.** Several of these
scripts will run happily in the wrong order and produce a confident wrong number.

Conventions used throughout the pipeline, stated once:

- **chain A = binder, chain B = target** in generator output and in every
  prediction. ⛔ Opposite to the previous generator version. Name the chain.
- **WT numbering** everywhere (mature EGFR = UniProt P00533 − 24).
- Residue identity is **asserted**; a number range is never trusted.
- Scripts **refuse rather than degrade**. A refusal is a result.

Interpreters: the scorers need `biotite`; the sequence step needs the
sequence-model environment. Two separate environments in our setup; neither path
is portable and both are environment variables or script constants you will have
to set. See `INPUTS.md` §4.

⚠️ **Some scripts name sibling scripts that are not deposited**, because they
state a constraint the deposited script depends on — the batch runners
(`prod_rfd3_egfr_multi.sh`, `fold_shard_pod.sh`), the cycle driver and its
backend, the submission-CSV builder, the deposit-tree builder, the site-B
counterparts of the site-A screens, and the surrogate-learnability check. Those
names will not resolve in this repository. The sections below describe what each of
them has to do, so the constraint is readable without the file.

---

# Stage 1 — target preparation

## `setup_rfd3_sphere_campaign.py`
**Does:** cuts the generation target as a sphere around each site's footprint and
writes the whole per-length, per-shard generator campaign (both sites, plus the
two-site arms).

**Expects:** `--outdir`. The receptor model and the generator's Python are script
constants. `--validate` checks each written spec against the **installed** schema
rather than the documentation.

**Emits:** one bundle per site — the sphere target PDB, the shard arm JSONs, and
the WT→output map beside them.

**Refuses:** if any hotspot is outside the sphere, or if a hotspot's residue type
does not match the input file. Both assertions run **before anything is
written**.

**Constants that matter:** site A radius 22 Å ⇒ 231 residues / 15 segments;
lengths 60–80 by 5; 3 shards per length; 336 designs per shard. See
`THRESHOLDS.md` T7, T12.

⛔ **One target PDB per bundle** — the runner counts the PDBs in its input
directory and dies if there is not exactly one. Hence one bundle per site.

---

## `make_wt2out_map.py`
**Does:** persists the piecewise WT→generator-output residue map for a sphere cut.

**Expects:** the sphere PDB, `--out <map>.tsv`, optionally
`--verify-cif <generator output .cif.gz built from THIS target>`.

**Why it exists:** the generator renumbers the target 1..N contiguous and
**discards both the original numbers and every segment gap**. With 15 segments
**no constant offset exists**. Nothing used to persist the real map — the builder
computed it and threw it away, so it was recoverable only by re-running the cut
at the identical radius.

**Verified, not assumed:** output chain B residue *i* (1-based, file order) **is**
input residue *i* (1-based, file order), checked by residue type, **231/231**.

⛔ Residues are read in **file order**, because that order is what the generator
renumbers.

---

## `wt2out.py`
**Does:** a module, not a script. Loads the map the bundle ships.

**Imported by:** the self-consistency scorer, the sequence step, the
full-receptor gate and all three backbone selectors.

**Why a module and not a constant:** every scorer on this project once carried a
constant offset. MEASURED: WT H409 is output B141, and the old formula gives B40,
which is a **cysteine**. One site failed loudly on a type assertion; another's
offset "worked" by luck because all of its footprint residues fell in the
sphere's first segment.

⛔ The map ships **inside the bundle, beside the target it describes**, so a
design can never be scored against a map for a different cut.

---

## `full_receptor_gate.py` — the exclusion check
**Does:** superposes each design back onto the **whole** receptor and measures
overlap with receptor the generation target did not contain.

**Expects:** `indir` (designs), `--receptor`, `--wt2out`, `--csv`. Optional
`--write-superposed`, `--glycan-chain`.

**Emits, per design:** `fit_rmsd` · `d_omitted` · `d_glycan` · `glyc_sequon` ·
per-sequon minima · `d_sphere` · `n_omit_2A` (the voiding metric) · `n_glyc_2A` ·
`verdict` ∈ {PASS, WARN_GLYCAN, CLASH_OMITTED, CLASH_GLYCAN}.

**Why it exists:** both arms of an earlier round were voided because every design
was built against a crop and **70% of site-A / 100% of site-B designs sat within
2 Å of receptor that had been deleted** — backbone overlap, not contact. The crop
was validated by **inclusion** and never by **exclusion**, and every downstream
stage treated the crop as the target, so the omission was unfalsifiable from
inside the pipeline. This is the missing exclusion check.

⛔⛔ **It must superpose.** The generator recenters its output by **66.93 Å**, the
same offset on every model. And because the target comes back rigid at
**0.0402 Å**, the fit RMSD is itself a tripwire — asserted against
`--max-fit-rmsd`.

⚠️ **The glycan tier is a WARN, not a clash, on purpose.** The crystal resolves
only 1–7 sugars per tree; native trees are ~15–20 residues and sialylated, with
the antennae absent from every structure. A glycan distance is a **floor on the
clash, never the clash itself**.

⭐ This is the single most reusable thing in the project. Any design campaign
against a trimmed target needs its equivalent.

---

# Stage 2 — backbone generation

## `prod_rfd3_egfr_multi.sh`
**Does:** the worker-side generator runner. One GPU, many shards,
uploading each shard as it finishes and **skipping shards already uploaded on
restart**.

**Expects:** worker-absolute paths under a `/workspace` root — correct for a
worker, not for a workstation.

**Note:** this arm ships with **no checkpoint/resume** — a recorded, deliberate
deviation, valid only while a shard stays under a bounded wall time, because a
lost shard is **re-rolled, not resumed**. ⛔ The bound is explicitly arbitrary and
exists to cap unprotected work; re-derive it from your own costs
(`THRESHOLDS.md` T12).

---

## `make_rfd3_shards.py`
**Does:** writes a shard set at one fixed binder length per shard, by cloning a
base arm JSON that already validates and changing exactly one thing.

⚠️ **This is NOT the generator of the submitted backbone.** It builds a 31-shard
set at lengths 70..100; the submitted backbone comes from the 15-shard sphere set
(`setup_rfd3_sphere_campaign.py`, lengths 60–80 by 5). Kept because the trap it
documents is real and general.

**The trap:** ⛔ **the binder length lives in the contig string's leading
integer**, not in a `length` field. The spec sets `extra="forbid"`, so a stray
key fails at load — and a `length` key alongside the contig would be silently
ignored by the thing that matters. The script rewrites the contig and **refuses
if it cannot find that leading integer**.

⭐ **Why one fixed length per shard:** it makes "is longer more novel" free to
measure — every design in an output archive has a known length, with no extra
bookkeeping. Design the sharding so a question you will ask later costs nothing.

⛔ **Not set here, CLI-only, not spec fields:** batch count, sampler parameters,
guidance. Production runs with guidance **OFF**: measured 2.76× the null, and
extra recycling 2.55× null-to-worse.

---

## `select_all31.py`
**Does:** the site-A backbone selection screen; length-agnostic, reads every
shard, so its output is not length-locked. Its CSV is what the sequence step was
pointed at.

**Expects:** `--wt2out <map>.tsv`. ⛔ **Required** — the old `W2O = -369` constant
is **gone**, not deprecated.

**Two deletions worth knowing about, because both were silent-wrong-answer bugs:**
1. The constant offset (above). MEASURED: the map is now checked against the
   design's target by **every** residue type, 231/231.
2. ⛔⛔ **A guidepost-excision step was DELETED, not bypassed.** It took the
   binder's nearest aspartate to the key histidine and **removed it from the
   binder** before measuring footprint contact. On hotspots-only arms there is no
   guidepost, so that deleted a real, designed residue — and it fired precisely
   on the designs that place an aspartate at the imidazole, i.e. on the
   interaction the arm existed to create. `--guidepost-pin` restores
   **pin-based** identification for a future pinned arm, with deliberately **no
   fallback**: a guidepost off its pin is not one.

⚠️ Two reported quantities are **measured, not gated**, and both are minima over
all of chain B — which grew from 106 to 231 residues, so each can only fall for a
purely geometric reason. Do not read a change in them as a change in the designs.

---

## `backbone_passrate.py`
**Does:** the **backbone-level** pass rate. Not per sequence: for each backbone,
pool **all** of its draws and report how many pass the contact gate.

**Expects:** `<contacts.csv> [out.csv]`.

⛔ **Backbones span shards** — the fold set is length-sorted and shards are
sequential chunks, so a backbone's draws can split across two shards. This
aggregates across every shard and **asserts the total**, because a backbone
silently scored on a partial draw set reads **too high**. That has already
reversed a winner on this project.

⛔ **Thresholds are reported at several levels, not one** — stated in the script:
"the threshold is a choice and should be made in the open".

⛔ Contact is primary; placement RMSD is a **diagnostic and never a cut**.

---

# Stage 3 — sequence design

## `mpnn_egfr.py`
**Does:** runs SolubleMPNN on the generator's binder backbones with the decided
flag set baked in.

**Expects:** `indir` (generator `*.cif.gz`), `outdir`, `--site {A,B}`,
`--wt2out <map>`, and `--run` to actually execute (otherwise it prints the
command — use that). `--no-guideposts` for hotspots-only arms. The sequence-model
directory and interpreter are environment-overridable with this-machine
defaults; **an outside reader must set both** (`INPUTS.md` §4).

**Flags it bakes in:** `--omit_AA C`, proline permitted, `--bias_AA` on D/E with
the **magnitude open**, `--chains_to_design A`. Provenance for all four:
`THRESHOLDS.md` T9.

⛔⛔ **The trap it exists to avoid: the sequence model will happily redesign a
deliberately placed residue.** Pinned residues are ordinary chain-A residues in
the output — nothing marks them — so an unguarded run silently mutates them and
every pinned design you paid for is wasted. This script **locates them per design
by their own contact geometry** and passes them as fixed. ⭐ MEASURED: the
generator chooses the pinned residue's sequence position itself, **68 distinct
positions over a 2–70 window** in an n=400 baseline, so the index **cannot** be
hardcoded.

⛔ The target is kept as **context** and never designed, so the model sees the
interface it is designing against.

⛔ **No constant offset any more.** The anchor index used to be `409 − 369`, valid
only for a contiguous crop. `--legacy-crop-offset` reproduces the old behaviour
and exists **only** to re-read the voided crop campaigns — and is gated behind a
second, deliberately embarrassing flag name. ⛔ Do not use it on sphere output.

---

## `seq_qc_egfr.py`
**Does:** the **sequence-level** half of the design filter, on fasta alone.

**Expects:** `seqdir` (the sequence output directory),
`--guideposts <guideposts.json from mpnn_egfr.py>`, `--csv`.

**Checks:** guidepost identity held (the real test that "fixed" worked) ·
cysteine count == 0 · N-glycosylation sequon `N-[^P][ST]` · net charge at pH 7.4
inside the band.

⛔ **The largest-positive-patch metric is NOT here and cannot be** — it needs
per-atom solvent accessibility on a folded structure, so it runs **after** the
fold-back. The filter splits because its inputs do.

⚠️ **It SKIPS fasta entry 0** — the sequence model writes the input sequence
first and it has no `id=` field. Counting it as a design inflates every
statistic. MEASURED: 101 records = 1 input + 100 designs.

⚠️ Histidine counts as **f = 0.5** at pH 7.4, chosen to match the separate
charge-gate script's convention **specifically so the two agree**.

---

## `make_w10_tsv.py`
**Does:** builds the 100-sequence TSV for one chosen backbone from the sequence
model's fasta.

**Expects:** `<fasta> <skeleton.yaml> <out.tsv>`.

⛔ **Sequence length is ASSERTED against the skeleton's chain A.** A mismatch
would mean a different backbone's sequences, and every downstream contact number
would be wrong.

⭐ It draws from the **same** fasta an earlier arm drew from (ids 1–100), so the
existing results map onto these sequences **by id with no re-derivation**.

---

# Stage 4 — fold prediction

## `build_ecto621_yamls.py`
**Does:** the canonical prediction-input builder — the designs against the **full
ectodomain**, ± glycans. This is where the whole construct was first got right.

⭐ **The target is not a model choice; it is read out of the assay data.** The
released table's full-sequence column is binder + a **constant 621-residue
suffix**, one variant across 400 of 402 designs. MEASURED: that 621-mer matches
all 604 residues of the structural model at the same WT index with **zero
mismatches**, so structure and assay construct are the same protein in the same
numbering.

⭐ Folding against the whole ectodomain fixes **two** problems at once: nothing is
deleted, **and** every design is free to choose its own face — so "engages the
footprint" stops being partly forced by the construct.

---

## `build_yamls_from_skeleton.py` — the one in the loop
**Does:** writes prediction inputs for new binder sequences by **reusing a real
landed input verbatim** and changing exactly one line: chain A's sequence.

**Expects:** `--skeleton <a landed input file>`, `--seqs <TSV: name TAB sequence>`,
`--outdir`, `--shards N`. `--require-template` defaults **on**.

**Why:** the construct — target chain, 10 NAG glycans, glycan bond constraints,
a forced template with its threshold, and the alignment path — was discovered by
measurement and the script that first built it was not kept. Rebuilding it from
scratch is the expensive way to make a silent mistake. Everything that made the
first campaign fold is carried **byte-for-byte**.

**Refuses:**
- no absolute chain-B alignment path in the skeleton — a relative path makes the
  predictor **skip the design and exit 0**;
- no `templates:` block — the runner enables potentials by **grepping for it**,
  and without potentials the forced-template flag is silently ignored. MEASURED:
  target fit **17.5 Å versus 0.69 Å**;
- more or fewer than one chain-A sequence line.

⭐ **Prove the output differs from a real landed input in exactly the intended
way.** Done for a 50-sequence round: every file differed in exactly the chain-A
sequence line, **50 of 50** — because a rebuilt construct would have silently
broken the pooled comparison.

---

## `build_ecto621_bundle.py`
**Does:** packs a campaign into the per-worker layout: one alignment, the runner,
one archive per worker, and a **per-shard draw-count file**.

**Expects:** `--yaml-root`, `--a3m`, `--arm`, `--outdir`, `--draws`. The
repository root is a script constant.

⛔ **One unit of work per worker, one output prefix per worker.** Two workers on
one prefix means the first to finish writes a complete-looking manifest and a
verifier can bless or terminate the wrong one.

⛔ It **asserts** that the inputs reference the absolute alignment path rather
than trusting them.

⛔ It ships **no templates directory** — the driver copies it in and
**md5-checks** it, because without it the predictor exits rc=0 having written
nothing.

⛔ MEASURED precedent carried in its header: glycosylated inputs exhausted a
16 GB card at ~1,260 tokens and **still exited rc=0** with "ran out of memory,
skipping batch" — a silent no-output failure. **Check the per-shard prediction
count against the input count on every glycosylated shard.**

⭐ The draw count is a **bundle-level file**, not a command-line flag baked into
the runner — which is what makes per-campaign draw counts possible at all.

---

## `fold_shard_pod.sh`
**Does:** folds one shard on the worker and uploads the structures.

⭐ A deliberate **copy**, not an edit, of an older runner that another campaign
depends on — and the copy's header lists its three differences: a configurable
draw count via the bundle file, a real target alignment, and the glycan handling.
⛔ Do not edit a runner a live campaign shares.

⚠️ Its fatal-error text still names a missing alignment for what is really a
partial-output condition. Left unfixed deliberately while the campaign was
running — **editing a running script does not fix it.**

---

# Stage 5 — scoring. RUN IN THIS ORDER

## 1. `contacts_a21.py`
**Does:** per-draw epitope contact.

**Expects:** `<land_dir> <out.csv>`. Searches recursively for `A_*_model_*.pdb`.

**Emits:** `n_epitope` · `n_epitope_total` · `n_target_contacts` · `c380` (1/0) ·
`epitope_hits` (the hit residues, `;`-joined).

**Definition:** the 11 footprint residues, **5.0 Å min over all atom pairs**,
chain A = binder / chain B = target, **numbering asserted by residue identity,
not range** — a range check was tried first and refused all 13 shard-sites,
because the folds are against the full 621-residue target.

⛔ **`c380` and `epitope_hits` were ADDED; nothing was changed.** `n_epitope`,
`n_epitope_total` and `n_target_contacts` are byte-identical to what this script
produced before, so **every landed table still reproduces**.

⛔ **`c380` is not a new measurement** — it is `380 in hit`, out of the same
geometry pass. Do not re-derive it with a different cutoff.

⚠️ A contacts file written before `c380` existed does **not** have the column, and
the only honest response is to re-run this script over the structures. That is
free and local.

---

## 2. `selfconsistency_ecto621.py`
**Does:** `sc_rmsd` (fold fidelity) and `pose_rmsd` (placement fidelity).

**Expects:** `--design-dir`, `--pred-dir`, `--wt2out`, `--out`, `--gate` (default
1.5), `--site {A,B,BR}`, `--limit`.

⛔⛔ **Why a second script and not a patch to the older one.** The older scorer
refuses these folds — *"target residues differ: design 231 vs pred 621 (shared
231)"* — and **it is right to refuse**. The generator renumbers the sphere target
1..N; the assay construct is numbered 1..621. The integer overlap 1..231 is a
**coincidence**: design residue 1 is WT 8, fold residue 1 is WT 1. Relaxing the
count check would superpose **mismatched residues** and return a confident,
meaningless placement number. ⇒ The target leg goes through the map, and residue
**identity** is asserted after mapping, never just the count — an anchor-only
check cannot catch a wrong map.

**Definitions kept identical to the older scorer so values read across:**
`sc_rmsd` = Kabsch on the binder **alone**; `pose_rmsd` = superpose on the mapped
target Cα, then binder Cα RMSD with **no** binder re-alignment ⇒ includes fold
error, so `pose_rmsd ≥ sc_rmsd` always.

⛔ **Cut on medians per backbone, never on counts.**

---

## 3. `score_cycle_egfrA1.py`
**Does:** pools draws to one row per sequence and writes the next cycle's
training table.

**Expects, positionally:**
`<base_training.csv> <contacts.csv> <sc.csv> <cycle.tsv> <cycle_map.csv> <out_training.csv> <out_parent.txt>`.
A replicate map, if any, is found **beside** `<cycle_map.csv>` as
`<cycle>_replicates.csv`.

**Emits:** `label` · `gate_frac` · `passes_floor` · `mean_nep_all` ·
`full_epitope_frac` · `median_sc` · `mean_tgt`. Definitions and provenance:
`OBJECTIVES.md` §1.

⛔ **The label definition is copied exactly** from the cycle-0 inline script and
the burn-in builder. Any drift silently changes what a surrogate is trained on
and the earlier cycles stop reproducing.

⛔ **The floor is a candidate filter, never a term summed into the label.**

⭐⭐ **The replicate hold-out.** A cycle may carry re-folds of earlier designs in
its own slots. Those are real measurements of sequences **already in the table**,
so they are **never merged** — merging would put one sequence in the table twice,
corrupt hub detection (which counts Hamming-1 neighbours), double-weight one
design in anything fitted on the table, and let a re-fold become the next cycle's
parent. They are scored, compared against their prior values, and written to
`<cycle>_replicates_SCORED.csv` with the test–retest read.

⛔ **Repeats are identified from that file ONLY, never from the id.** "id9" as a
prefix already matches ten genuine designs (id9, id90, id91..id99, id901, id902).
If the file is absent the script behaves **exactly** as it did before repeats
existed — which is what makes cycles 0–3 reproduce byte-for-byte.

⚠️ **It computes the older metric**, not the current objective. Every landed
training table predates the objective change. **Label the metric.**

---

## 4. `score_0047.py` — the current objective, as a separate pass
**Does:** implements the current objective verbatim (`OBJECTIVES.md` §2).

**Expects:** `--pair <contacts.csv> <sc.csv>` (repeatable, so several cycles pool
and a sequence id reused across cycles cannot collide — each pair is tagged by
its contacts directory name), `--tsv` to carry sequences through, `--out`,
`--alias-map`, `--augment <training.csv> --augment-out <out.csv>`.

⛔⛔ **It REFUSES rather than degrades.** `c380` must be present in the contacts
CSV. There is **no fallback, no inference of `c380` from `n_epitope`, and no
silent default** — that failure mode is exactly the one that publishes a stale
table at rc=0.

⛔ The mean is over **every** draw, not over gate-passers. This is the one place
it deliberately differs in **shape** from the older scorer. Do not "fix" them
into agreement.

---

# Stage 6 — proposing

## `propose_cycle_egfrA1.py`
**Does:** finds the reference, recomputes the accepted set, enumerates untried
combinations, reserves replicate slots.

**Expects:** `--training`, `--cycle`, `--cap`, `--out-tsv`, `--out-map`,
`--design`, `--min-hub-margin` (default 3), `--replicate-top` (default 5),
`--replicate-anchor` (default `auto`), `--out-replicates`.

**Refuses** on: a thin reference margin · a reference that fails the fold floor ·
an empty accepted set · a fully-enumerated set. ⭐ **Each refusal is a result, not
a failure.**

**The reference rule, and it is the one inference in the file made measurable:**
the accepted set is only meaningful relative to the sequence whose single mutants
were actually **measured** — the **hub**, the training row with the most
Hamming-1 neighbours. MEASURED on the cycle-2 table: **50** neighbours against
**4** for the runner-up. Not a tie, not a judgement call.

**The accepted-set test:** beat the reference on **BOTH** `label` **and**
`full_epitope_frac` — `label` saturates and can no longer discriminate, and the
fraction alone is noisier at 10 draws.

⛔ **One substitution per position.** Two accepted mutations at the same position
are mutually exclusive; a combination takes one or the other.

⛔ **Never re-fold** anything already in the table — except the declared repeats,
which are the one documented exception.

⭐ **Enumeration order is deterministic:** lowest untried order first, then
lexicographic by position tuple ⇒ a cycle is reproducible from its training table
alone.

⛔⛔ **No model in this loop, measured.** See `OBJECTIVES.md` §6. Do not
re-introduce acquisition-based proposal without re-measuring out-of-bag error
first.

---

## `build_combo_cycle_egfrA1.py`
**Does:** the hand-path builder that produced cycle 2. Enumerates **combinations**
of an accepted set rather than single mutants of one parent. The driver's proposer
was written later and reproduces its sequence set exactly.

**Why combinations:** cycle 1 ended with **seven** mutations that beat the parent
on both metrics and were **mutually indistinguishable** (pairwise Welch dmean
+0.00, Fisher p=0.628–1.000) while `label` had saturated. Picking one winner would
have been arbitrary, so **all** the indistinguishable top scorers became the
baseline — a decision on the record.

**The accepted set, as run:** `L16D A17M T19Q L49V L49K R55S S75A` — 7 mutations
over 6 positions ⇒ **88 combinations**. Cycle 2 took 50 (orders 2–3), cycle 3 the
remaining 38 (orders 4–6). **There is no order 7.**

---

## `build_sitescan_cycle_egfrA1.py`
**Does:** builds a **single-position saturation scan** instead of combinations.

**Why:** the combination walk is **complete, not stalled** — 88 enumerated, the
proposer returns 0 new at every order. ⛔ Re-running the combo builder produces
nothing, and that is not a bug to fix.

**And the deeper reason:** the accepted set **cannot reach** the one thing left to
optimise. MEASURED on the 98-sequence cohort, the current objective equals the
plain F380 contact rate exactly for **63 of 98** sequences and for **all ten**
re-measured top designs — among reliably-folding designs the objective **is** the
F380 rate, because the other ten footprint residues are already satisfied
whenever F380 is. One binder position effectively governs F380
(`OBJECTIVES.md` §3).

⛔ **Check for collisions before launching.** One substitution in the planned scan
already existed as an earlier design; including it would have collided and made
the scorer exit **after the money was spent**. Caught by hand **and** against the
scorer's own guard, and excluded.

---

## `build_burnin_egfrA1.py`
**Does:** turns 100 sequences × 10 already-folded draws into one row per sequence
— the campaign's starting table.

⭐ Keeps below-floor sequences with `passes_floor = 0` rather than **dropping data
silently**.

---

## `stage_egfrA1_cycle0.py`
**Does:** stages the four-candidate head-to-head at 10 draws each — two consensus
sequences plus two individual designs.

**Why it belongs in the record:** cycle 0 **refuted the consensus starting
sequence** and set the reference. A reader following the narrative needs it. A
measurement that refutes your starting point is the cheapest one in the campaign.

---

## `singlesite_surrogate.py`
**Does:** enumerates all 75 × 19 = **1,425** single mutants and scores them
exhaustively with the model ensemble. No genetic algorithm, so the ranking is
**deterministic**.

**Status:** ⛔ **out of the loop**, by measurement. It is also the file that
measured the surrogate out.

⚠️ **Packaging problem, flagged not fixed:** it inserts a path and imports a
module that is not part of this deposit. Everything is imported rather than
retyped, so drift shows up as a diff — but publishing this script means
publishing that module too, or stubbing the import.

---

⚠️ **Worker-absolute paths.** Several scripts here hard-code `/workspace/msa/target.a3m`
(`build_ecto621_bundle.py`, `build_ecto621_yamls.py`, `build_yamls_from_skeleton.py`). That path is
the layout **inside the compute worker's container**, not a path on any local machine, and the
bundle builder asserts on it deliberately so a YAML that would not resolve on the worker fails
early. ⇒ An outside reader re-running these must substitute their own worker path; the scripts are
deposited as they ran, not rewritten to be portable.

---

# Supporting analysis — not in the runnable chain

## `epitope_perresidue_egfrA1.py`
**Does:** per-residue epitope breakdown — **which** of the 11 residues a draw
lost.

**Why it exists:** the contacts scorer computes the hit set and records only its
**size**, so every downstream table could say a draw missed the ceiling but not
**which** residue was lost. MEASURED: misses are overwhelmingly **one residue
short**, not a collapse — and *"recover one specific residue"* is a targetable
design objective where *"the fraction is below 1"* is not. **This is the analysis
that identified F380.**

⛔ Contact geometry is **copied** from the scorer, not reinvented: same cutoff,
same residue list, same chain convention, same test, same identity assertion. If
one changes the other must change with it or the two disagree silently. ⛔ Do not
"improve" the cutoff here.

## `learnability_egfrA1.py`
**Does:** the **zero-compute pre-spend gate** — is the label learnable from the starting
sequences at all? Ridge + random forest on encoded sequence.

**Why:** if the label cannot be predicted from sequence at all, the loop stalls at
cycle 1 and you would only find out **after paying**.

⚠️ **It is a FLOOR, not the production surrogate.** A positive means "there is
signal a simple model can find". **A negative is the informative result.**

⚠️ And it is the **harder** problem in one respect: held-out starting sequences
sit ~21 mutations from the training set, while the optimizer's proposals sit 1
mutation from a parent with data nearby. ⛔ Do not read a weak r here as "the loop
cannot rank single mutants".

## `analyze_stacking_egfrA1.py`
**Does:** answers the cycle-3 fork — do beneficial single mutations stack? Prints
the comparison **and the criteria** and **deliberately does not decide**.

⭐ A script that lays out the criteria and leaves the call to a human is the right
shape for a fork. Worth copying as a pattern.

---

# The driver

## `autocycle.py` + `cycle_lib.py` + `egfr_backend.py`
**Does:** chains score → train → propose → stage → arm → launch → land across
cycles. The campaign plugs in as a **backend**: the cycle machinery is shared,
every per-step program is this campaign's own.

**What it carries that matters scientifically:**
- a **double-spend guard** that asks the filesystem whether a round already
  landed;
- `--resume` **idempotency** — which is what let a hand-built cycle be picked up
  without re-deriving it, because the propose step is skipped when the cycle's
  TSV already exists. ⭐ That is how to run any cycle the proposer cannot derive;
- three **refuse-gates**, each of which cost money or a wrong answer to learn:
  **per-shard output count** (the predictor's out-of-memory exits rc=0, so the
  count is the only proof), **templates-directory md5**, and a **file-set diff of
  the new bundle against the previous cycle's**;
- a **verified-terminator gate** before each launch — it refuses to launch
  without one;
- `--only score,train` so a landed round can be scored with **nothing proposed,
  staged or launched**, i.e. at zero cost.

⛔ **It never decides what a cycle designs.**

⚠️ **Known limits, stated because they are reproducibility-relevant:**
- it does **not** survive the controlling shell exiting; an unattended-launch mode removes
  the prompt, not the fragility;
- the backend hardcodes the campaign directory in **7 places** (MEASURED; an
  earlier estimate of "≥8" is withdrawn) plus **two
  interpreter paths** that an outside reader does not have. ⛔ **Parameterise, do
  not fork** — a forked copy drifts invisibly, a parameter shows up in a diff;
- a ledger line is written **unconditionally** before the launch branch, so every
  resume re-appends one. ⇒ **Count distinct workers for spend, never ledger
  lines.**

⛔ **The launch and termination layer is not documented here and is not proposed
for publication.** `README.md` §5 states what has to happen on the compute side,
and why, without the machinery.

---

# Deliberately omitted from the documented path

| what | why |
|---|---|
| the crop-era builders, stagers and length-band analyses | every fold they produced is from the voided round; superseded by the sphere path |
| the older self-consistency scorer | it **correctly refuses** these folds. It is kept in the lab record *as the reason the second scorer exists*; in a public tree it would read as the live scorer |
| the training-table merge script | ⛔ **that path never worked.** It pools sequences already in the base table and **skips new ones** — which is every sequence in a mutant cycle. The cycle scorer exists because of it |
| ⚠️ the two-site bridge family | **The bridge's `id85` is rank 1 of the deposited set, but its closed-loop optimizer was never built** — the arm was decided by two replicated depth rounds, not by cycling — so the proposer and driver for it do not exist. ⇒ the **scoring** chain it was decided on (`contacts_bridge.py`, the shared self-consistency scorer with its two-site option, `score_br50.py`, `bridge_depth2_analyze.py`) is reproducibility-relevant and is listed separately in `DESIGN_CHOICES.md` §3. Publishing it is a separate decision from publishing the site-A chain |
| the second site's family | parked; no sequence from it is in the submission |
| the pH-switch exploration | real upstream **reasoning**, not a step that emitted a sequence. A README paragraph, not code |
| the metric-calibration family against prior designs | it chose *which* metric to score on; it generated no sequence. ⛔ And one of its outputs **must be quoted stratified** — a trap for an outside reader (`OBJECTIVES.md` §5a) |
| the pinned-residue ("guidepost") fork | anchors were ruled **QC-only, not a design lever**. Not on the shipped path |
| figure and one-off analysis scripts | outputs, not the generation path |
| third-party vendored code | separate licences. ⛔ Cite, do not redistribute |
| every `*.pre-*` backup | pre-edit backups by convention. Never publish |
