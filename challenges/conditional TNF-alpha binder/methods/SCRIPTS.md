# Scripts — what each one does, what it expects, what it refuses

Every script named in this table is present in `../scripts/`. ⛔ A citation to a program that is
not deposited is a dangling reference, not provenance — that was a defect in our challenge-1
deposit and it is not repeated here.

**Run order.** Stages 2–5 of `README.md` run left to right:

    j484 (build)
      -> [generation, off-site]
      -> j393 (continuity gate)   +   j409 (topology gate)      ... both, then intersect
      -> j168 (sequence design)
      -> j496_make_boltz_inputs   -> j496_run_stage2 (fold)
      -> j497 (score)             [j465 imported for geometry]
      -> j494 / j498 (figures, each with its own verifier)

⛔ **These are lab programs, not a packaged tool.** They take absolute paths, several carry
**placeholder** paths marked in the file, and they assume the directory layout of the campaign that
produced them. They are deposited so the method can be *read and checked*, not so it can be
re-run unchanged.

---

## Generation

### `j484_make_c1_core_yield_bundle_20261017.py`
Builds the `v46` generation bundle: three arms × 8 shards × 16 designs = **384**. Writes the contig,
the orientation token, the per-arm core selection and a human-readable note into the bundle, then
**gates its own output before it can be used**: it re-derives the arm tokens and refuses if the
downstream scorer's own regular expression would parse one differently; it refuses to overwrite an
existing bundle; it re-digests the runner against a pinned SHA1; and it prints the resolved
configuration rather than inheriting a default.
⛔ Writes a **new** directory; never edits a source bundle. Does not launch anything.
⚠️ `REPO` is a **placeholder** path.
⚠️ Its wall-time projection is retained; the rental-cost projection was removed for this deposit.

---

## The two backbone gates — ⛔ run continuity first

### `j393_continuity_gate_20261015.py`
Is the designed chain **one connected chain**? Reports the maximum Cα–Cα step along the chain, where
the worst step falls relative to the pinned span, and a separate count of numbering gaps.
- Bars: `--step-warn 4.2` (the continuity bar) and `--step-break 10.0` (catastrophic), both printed.
- ⛔ **Chain roles are detected, never named** — from the fraction of residues carrying the guidepost
  flag. The chains actually used are printed *and* written into the output, so a wrong assignment is
  visible rather than silent. It **skips** a file with no guidepost column rather than guessing.
- ⛔ **One row per (design, binder chain)**, never collapsed: a trimer with one broken copy is not
  "two thirds good".
- ⛔ The primary column is the continuous `max_step`, not the pass/fail flag.
- ⛔ Refuses an existing `--out` (these scripts **append**, so re-using a path silently merges two
  campaigns) and refuses an empty glob — an empty glob is not a result.
- Reports distances at **three** decimals, deliberately: a step of 4.204 Å prints as "4.20" at two
  decimals and then passes a filter written against the printed value.

### `j409_path_or_fold_20261015.py`
Is the binder **one fold**, or a **path** from pin to pin? Reports `min_cut` (fewest long-range
Cα–Cα contacts crossing any interior cut, `--margin 13`), the long-range contact fraction, and Rg
against the globular expectation.
- ⛔ `min_cut` is **margin-sensitive**: reference controls are recomputed at the same margin as the
  designs and never quoted across settings.
- ⛔ Roles are detected for our designs; chains **are** named for the reference controls,
  deliberately, and each is printed.
- ⛔ States in its own output that none of its numbers is an affinity statement and none is
  sufficient alone.
- Same `--out` and empty-glob refusals as `j393`.

---

## Sequence design

### `j168_mpnn_tnf.py`
SolubleMPNN driver for this target. The invocation itself is carried over from challenge 1's
deposited `../../pH-dependent EGFR binder/scripts/mpnn_egfr.py`, so the settled flags transfer; what
differs is that this target has no WT→output residue map and no single anchor type to assert.
Three things it guards, because getting them wrong is **silent**:
1. **Chain identity, per design, every chain.** Roles come from the guidepost-flag fraction, not from
   length. More than one binder chain **refuses** and points at the symmetric driver. An earlier
   version inspected only chains A and B, so a three-copy file matched both and sailed through —
   after which one copy was designed against two neighbours still carrying the generator's own
   mostly-alanine sequence, with no error raised.
2. **The pinned motif is passed as `--fixed_residues_multi`.** Without it the model silently
   overwrites the whole Arg/His motif. ⛔ Positions are read **per design** from the file; a zero
   count refuses.
3. **No cysteine in the designed chain**, checked after the run on the **designed chain's** FASTA
   field — the FASTA holds all chains colon-separated, so checking field 1 can check the target.
⚠️ `LIGANDMPNN_DIR` and `MPNN_PY` have **placeholder** defaults; set them by environment variable.

---

## Fold prediction

### `j496_make_boltz_inputs_20261018.py`
Builds one prediction input per designed sequence: designed binder as chain **A**, the whole trimer
as **B/C/D**, chain naming matched to the design file so no remap is needed downstream.
- ⭐ The target sequence is **read from the design file's own target chain** and asserted
  byte-identical to the 152-residue chain a previously proven input set used. It is never retyped.
- ⛔ The designed chain is pulled from the FASTA **by length**, and a record without `id=` is the
  input sequence, not a design, and is skipped.
- ⛔ Re-asserts **no cysteine** on the sequence that actually goes to the predictor — the last point
  before the fold.
⚠️ `REPO` is a **placeholder** path.

### `j496_run_stage2_after_stage1.sh`
Runs the remaining 99 folds, chained to start when the first stage's process exits.
- ⛔ **Chained, not parallel**: each fold holds ~9.4 GiB on a 16 GiB card.
- ⛔⛔ **Verifies an output COUNT at every stage, never an exit status** — an out-of-memory exit
  returns 0 having written nothing. If stage 1 is short it refuses to start stage 2.
- Carries the prediction settings used: `--diffusion_samples 1 --recycling_steps 3 --seed 20261018
  --no_kernels --output_format pdb --write_full_pae`.
⚠️ `D` and `BOLTZ` are **placeholder** paths.

---

## Scoring

### `j497_c1_boltz_selfconsistency_20261018.py`
The scorer behind every number in `../README.md`. **103 columns**, one row per prediction. Measures
per-chain and complex pLDDT, `ptm`/`iptm`/confidence, **per-chain-pair** `iptm` and **per-target-chain**
PAE (never pooled), binder Rg, per-protomer contact counts and chains engaged, the Arg→Glu107 and
Arg→Glu110 bridges **both over any arginine and restricted to the pinned one** with the responsible
residue number, fold fidelity (`rmsd_binder`, `tm_binder`), and placement fidelity
(`rmsd_binder_on_target_frame`) with all three cyclic chain mappings tried and a target-frame
control that **flags** a row rather than reporting an untrustworthy placement.
- ⭐ **Geometry is imported from `j465` under a SHA1 pin, not copied.** If the sibling drifts, this
  **refuses** rather than scoring with different math. ⚠️ The pin was re-pointed for this deposit;
  the file says why and names the hash the published numbers were produced under.
- Includes an **mmCIF reader that parses by declared column name, not by position** — the loop
  header is explicit, and a PDB-column parser reads a CIF as **zero atoms**, which looks exactly
  like an empty file.
- Reads guideposts only from the generator's own flag column; never infers them.
- Addresses target residues by **sequential index**, and asserts the residue is actually GLU before
  measuring from it.
- Echoes the resolved layout and configuration on every run.

### `j465_boltz_selfconsistency_20261017.py`
⛔ **Deposited as a dependency, not as a runnable part of this pipeline.** It is the
**C3-symmetric** sibling scorer; its own cohort and its own design file are *not* deposited, and
its `DESIGN` constant is a **placeholder**. It is here because `j497` imports `kabsch_rmsd` and
`tm_score` from it under the SHA1 pin, so removing it would make `j497` refuse to start. Read it for
the shared geometry; do not expect it to run against this deposit.

---

## Figures

Each builder has a **separate verifier**, and that separation is itself a finding: an in-process
`reinitialize()` + `load()` **segfaults** at this object count and takes the buffered stdout with
it, so a perfect session file prints nothing and reads as having done nothing. The verifier opens
the saved file in a **fresh process**, which is also exactly what opening it elsewhere does.

| script | what it builds / checks |
|---|---|
| `j494_v46_joint11_pse_20261018.py` | the 11 JOINT **designed** backbones in one session, 16 objects / 13 scenes. ⛔ Re-derives the JOINT gate in-script and **asserts** it reproduces 23 continuous and 11 JOINT with the 5/6/0 arm split — the figure cannot ship disagreeing with the gate. ⛔ Quotes this cohort's own distribution, never another campaign's. ⚠️ `REPO`, `CONT_TSV` and `PATH_TSV` are **placeholders**; the two gate tables are deposited in `../data/backbone_gates_20261018/`. |
| `j494_verify_pse_20261018.py` | opens that file cold and asserts the object and scene counts |
| `j498_pse_boltz_c1_20261018.py` | designed **and** predicted positions in one shared frame, for the scored folds. ⛔ The frame transform is applied in NumPy and written as **coordinates**, never handed to the viewer as a matrix (row-major vs column-major is a silent-wrong-answer risk), and each written file is **re-measured and asserted against the scorer's own placement number** before the viewer sees it. ⛔ Scene captions are asserted **ASCII** — a non-ASCII character corrupts the session and it segfaults on load. Scene calls are emitted inside a Python block because the viewer's own command parser treats a comma in a caption as the next argument, which silently loses scenes while still exiting 0. |
| `j498_verify_pse_20261018.py` | loads the saved file cold, recalls every scene, counts atoms in all three populations and **re-measures every placement out of the loaded session**. ⛔ Writes its result to a **file**, not stdout: `pymol -cq` echoes the script source, so a stdout grep for a sentinel fires on the echoed source whether that branch ran or not. |

⛔ **Neither session file (`.pse`) is deposited.** They are figures of the deposited coordinates and
tables, and both exceed what is useful to ship here; the builders and verifiers are deposited so the
figures can be rebuilt and, more importantly, so the assertions they make are auditable.
