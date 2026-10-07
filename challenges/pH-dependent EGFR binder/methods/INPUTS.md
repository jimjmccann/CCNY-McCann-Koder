# INPUTS — what the pipeline requires, and which parts are site-specific

Site-A EGFR binder campaign. Companion to `README.md` and `SCRIPTS.md`.

⛔ **None of the data files below are deposited in this repository.** They live on
a local bulk volume in the lab. ⚠️ **"the code tree" below always means the lab's own
working repository, NOT this deposit** — a file described as "mirrored in the code tree"
is mirrored *there*, and is still **not present here**. Nothing in this section is
retrievable from this archive. This file says exactly what each one is, which stage consumes it,
whether it is site-specific, and what an outside reader has to substitute. Sizes
and record counts were **MEASURED on 2026-10-05** so a substitute can be
sanity-checked against them.

---

## 1. The five external inputs the optimizer declares

The arm configuration declares exactly five, and **every one is REQUIRED to
exist** — the backend refuses rather than defaulting.

⛔⛔ **Every one of the five is SITE-SPECIFIC, and pointing any of them at
another site's file silently scores the WRONG THING rather than erroring.** The
live precedent in this lab: two residue maps for two different cuts differ, and
**only an identity assertion catches the swap**. That is why they are required
rather than defaulted.

| # | input | what it is | MEASURED size | site-specific? |
|---|---|---|---|---|
| 1 | **skeleton** | one real landed prediction input file, used verbatim as the construct template. `build_yamls_from_skeleton.py` changes only chain A's sequence | 1,882 B | ⚠️ yes — it encodes the binder length and the target construct |
| 2 | **a3m** | the chain-B (target) multiple-sequence alignment | 3.42 MB | no — it is the target's, shared across sites |
| 3 | **template_cif** | the forced structural template: the tethered full-length receptor | 319,611 B, md5 `15a67d7d49cb02ca45d05fae1573b96d` | no — same receptor |
| 4 | **design_dir** | the generator's output backbones to score against, as `.cif.gz` | directory, **10 files** | ⛔ **yes** |
| 5 | **wt2out** | the persisted WT→generator-output residue map for this cut | 3,531 B = 4 header lines + **231 mapped residues** | ⛔⛔ **yes — this is the one that silently reads the wrong residue** |

**The map file's own header states the three facts that matter**, verbatim:

    # WT -> RFd3 output residue map for egfr_siteA_sphere22A_WT.pdb
    # out_resnum = 1-based rank in file order; TARGET is chain B in RFd3 output.
    # BINDER is chain A. Piecewise over 15 segments: NO constant offset exists.

⭐ **The md5 of the template is CHECKED, not assumed.** The bundle builder ships
no templates directory; the driver copies it in and verifies the checksum —
because without the template the predictor exits **rc=0 having written nothing**.

---

## 2. Upstream data, by stage

### Stage 0 — the target

| input | what it is | notes |
|---|---|---|
| the receptor model | tethered full-length extracellular region, **WT-numbered** | a script constant, not a flag |
| the assay construct | binder + a **constant 621-residue** target, read out of the competition's released table — a single variant shared by 400 of 402 rows, so it is the assay construct **as published** | ⭐ MEASURED, and **not** whole-sequence verification: against the structural model at the same WT index, **604 of 604 modelled residues match, zero mismatches** — **97.3% of the target**. The remaining **17 residues are absent from the model** (WT 1-3, 295-299, 307, 577, 615-621) and are unverified by that comparison. Whether the published sequence matches the physical protein in the tube cannot be verified here and is not claimed |


⛔ **The target is read out of the assay data, not chosen.** That is what makes a
design's apparent epitope a measurement.

### Stage 1 — target preparation

| produced | what it is | MEASURED |
|---|---|---|
| the sphere target PDB | 22 Å sphere around the site-A footprint | 231 residues, **15 segments** |
| the WT→output map | input #5 above | 231 mapped residues |
| 15 shard arm JSONs | `arm_shard00 .. arm_shard14` | contigs 60/65/70/75/80, three shards per length |

### Stage 2 → 3 — the sequence ancestry of every submitted design

⭐⭐ **This is the single most important data provenance fact in the campaign:**

| input | what it is | MEASURED |
|---|---|---|
| the chosen backbone | `arm_shard11_shard11_24_model_2` | `arm_shard11`'s contig reads `75,...` ⇒ **75 residues**, matching the optimizer's `seq_len = 75` |
| **the sequence-model fasta** | **every submitted sequence descends from this one file** | 43,743 B, **101 records = 1 input sequence + 100 designs** |
| the same 100 as a TSV | the staged form | 11,392 B |

⚠️ **Record 0 is the INPUT sequence, not a design** — it has no `id=` field.
Counting it inflates every statistic, and the QC script skips it deliberately.

### Stage 5 → 6 — the campaign tables

| table | MEASURED size | note |
|---|---|---|
| `training_c1.csv` | 11,942 B | mirrored in the lab repository; ⛔ not deposited here |
| `training_c2.csv` | 17,895 B | mirrored in the lab repository; ⛔ not deposited here |
| `training_c3.csv` | 23,889 B | ⛔ **bulk volume only** |
| `training_c4.csv` | 28,435 B | 240 sequences |
| `training_c5.csv` | 30,575 B | 258 sequences, newest at the time of writing |
| `cycle90_map.csv` | 516 B | the first replicate round's map, `replicate_id → source_id` plus prior values |

⭐ **What IS in the lab's working repository** (⛔ **not in this deposit** — listed so an
outside reader knows the trail exists and what to ask for, not so they can open it here):
cycles 0–2 — `cycle0/1/2.tsv`, `cycle1/2_map.csv`, `contacts_c0/c1.csv`,
`sc_c0/c1.csv`, `training_c1/c2.csv`, `parent_c1/c2.txt`, the proposal CSVs — plus
the burn-in and single-backbone scoring tables.

⛔ **Cycle 3 onward and the replicate rounds exist only on the bulk volume.** That
is the reproducibility gap in the data, not in the code. ⛔ **Restated plainly for an
outside reader: none of the tables in this section ship with this archive.**

---

## 3. What is NOT distributable, and what to substitute

| not distributed | why | substitute |
|---|---|---|
| structure-prediction weights | third-party | obtain from the upstream project |
| sequence-model weights (`solublempnn_v_48_020.pt`) | third-party | obtain from the upstream project |
| generator weights | third-party, and shipped inside a container image in our setup | obtain from the upstream project |
| vendored third-party scoring code | separate licences. ⛔ **Cite, do not redistribute** | obtain from each upstream project |
| raw structures (hundreds of thousands) | size | regenerate, or ask |
| the launch / termination / account layer | ⛔ carries account plumbing and is not part of the method | any batch system. `methods/README.md` §5 states what must happen and why |

⛔ **FILENAMES ONLY for credentials.** Credentials in our setup are **files, not
commands**, and live outside the code tree. Nothing in this pipeline reads one
except the launch layer, which is not published.

---

## 4. Environment paths an outside reader must set

These are the places with this-machine defaults. All are either environment
variables or single named constants; none is buried.

| where | what | portable? |
|---|---|---|
| sequence step | the sequence-model directory and its interpreter | ⚠️ **environment-overridable**, but the defaults are this machine |
| bundle builder | the repository root | a named constant at the top of the file |
| backbone shard builder, campaign builder | the generator's interpreter | a named constant |
| single-mutant surrogate | a path inserted to import a module from **another campaign** | ⚠️ packaging problem, flagged not fixed |
| the backend | the campaign directory, hardcoded in **7 places**, plus **two separate interpreters** | ⛔ **parameterise, do not fork** |
| the worker runners | worker-absolute `/workspace` paths | ✅ correct for a worker |

⛔ **Two separate interpreters** are genuinely required in our setup: the scorers
need one environment (for `biotite`) and the sequence step needs another. An
outside reader needs both, or one environment carrying both dependency sets.

---

## 5. The bulk-storage rule, stated because it is a reproducibility hazard

Heavy intermediates go on a dedicated bulk volume, never the root filesystem.
This is not tidiness: a full root filesystem has, on this machine, **disabled the
tooling needed to fix a full root filesystem**, and separately caused a
verifying terminator to block silently for **6 h 21 min** while every liveness
check reported healthy.

⚠️ **The practical consequence for a reader of our record:** several of our own
"file is missing" conclusions were wrong because a relative path resolved under
the code tree instead of the bulk volume. ⛔ Check with an **absolute** path, and
**without** suppressing stderr. "No output" is not "gone".

---

## 6. Minimum set to reproduce one cycle end to end

Given the weights above and the environments of §4:

1. The receptor model, WT-numbered, and the assay construct's target sequence.
2. The sphere target PDB **and** its WT→output map — or run stage 1 to build both.
3. The chosen backbone's generator output (`.cif.gz`) — or run stage 2.
4. The sequence-model fasta, or run stage 3 on the backbone.
5. **One real landed prediction input file as the skeleton**, plus the target
   alignment and the forced template. ⛔ These three are the construct; do not
   rebuild them (`SCRIPTS.md`, `build_yamls_from_skeleton.py`).
6. A training table to propose from — or `build_burnin_egfrA1.py` on a folded
   starting set.
7. A structure predictor and somewhere to run it, with the output **count**
   verified per unit of work (`README.md` §5).

⛔ **Order matters** (`README.md` §6). ⛔ **Fresh output path per experiment** —
the scorers append, so pointing a run at an existing results file silently
**merges two campaigns**, with no error and no symptom until an analysis averages
across both.
