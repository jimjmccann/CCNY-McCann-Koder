# ⛔⛔ Which cohort is which — read this before quoting any number from `data/`

This directory holds **three different cohorts of designs**. They are not interchangeable, and
**only one of them was submitted**. A count taken from the wrong one is arithmetically correct and
scientifically wrong, so every table below is labelled by the cohort it belongs to.

| directory | cohort | n sequences | **submitted?** |
|---|---|---|---|
| `mpnn_20261018/` · `scoring_20261018/` | `v46 joint11` (C1) | **110** (10 × 11 backbones) | ⛔ **NO — not one of these was submitted** |
| `submission_20261019/scoring/c1_*` | `v46 rebias` (C1) | **88** (8 × 11 backbones) | ✅ **17 of these** |
| `submission_20261019/scoring/c3_*` · `seqs_c3/` | `v50` (C3, symmetric) | **96** (8 × 12 backbones) | ✅ **3 of these** |

⭐ **The authority on what was submitted is `../submission_challenge2_20261019.csv`** (20 rows).
`submission_20261019/selected_20_with_scores.tsv` carries the same 20 with an explicit `arm` column.

## ⛔ The trap, stated plainly
The **110** cohort is the one described in the challenge `README.md` section *"The result"* and
throughout `methods/README.md`. It is the fullest characterisation we have, which is why it is
deposited — but **it is not the submission**. The submitted cohorts' own numbers are in the
challenge `README.md` section *"The submitted cohorts' result"* and in `../methods/C3_ARM.md`.

## File-by-file

**`mpnn_20261018/` — `v46 joint11`, ⛔ NOT SUBMITTED**
- `seqs/*.fa` — 11 files, 11 records each = **11 inputs + 110 designs**
- `winners.txt` — the 11 backbone filenames this cohort was designed on
- `fixed_residues_multi.json` — the pinned residues, 11 keys

**`scoring_20261018/` — `v46 joint11`, ⛔ NOT SUBMITTED**
- `j497_ALL110_selfconsistency_20261018.tsv` — **110** rows, the cohort's scored folds
- `manifest.tsv` — **110** rows, one per sequence (`tag`, `binder_aa`, `total_res`)
- `j496_stage1_confidence_…tsv` (11) · `j497_stage1_…tsv` (10) · `j497_probe_…tsv` (1) — the
  first-round subset and the single probe fold

**`submission_20261019/` — ✅ THE SUBMITTED SET and its scoring**
- `selected_20_with_scores.tsv` — the 20, with `arm`, backbone and rank value
- `scoring/c1_holo_selfconsistency_framefixed.tsv` (**88**) — ⭐ use the **`_v2`** columns; the
  non-`_v2` frame columns are superseded (34 of 88 were mis-scored before the correction)
- `scoring/c1_baker_five_criteria.tsv` (**88**)
- `scoring/c3_apo_…` · `c3_trimer_selfconsistency` · `c3_trimer_assembly_fit` (**96** each)
- `scoring/c3_holo_selfconsistency` · `…_framefixed` · `c3_holo_assembly_fit` (**24** each)
- `seqs_c3/*.fa` — 12 files, 9 records each = **12 inputs + 96 designs**
- `c3_design_plan.json` — 12 backbones, **12 distinct pin sets**

## ⛔⛔ Column direction in the C3 tables — `_best` is not always better
The per-copy aggregates in the `c3_*` tables are emitted as `_mean` / `_best` / `_worst`, where
`_best` is the **maximum** over the three binder copies and `_worst` the **minimum**, *regardless
of which direction is good for that column*. So:

| column family | higher is better? | the **conjunction** (all three copies) is |
|---|---|---|
| `plddt_*`, `tm_*`, `iptm_*` | yes | **`_worst`** |
| `bridge_copy_best_*`, `rmsd_binder_copy_*`, `pae_*` | ⛔ **no — lower is better** | ⛔ **`_best`** |

⇒ selection ranked C3 on **`bridge_copy_best_best`**, which is the *worst copy's best bridge*.
Full explanation in `../methods/C3_ARM.md` §3a.

## ⚠️ Reading these tables
Several `.tsv` files here begin with `#` comment lines carrying the cohort label (the same
convention `backbone_gates_20261018/j490_continuity_…tsv` already uses). ⛔ **Skip lines
starting with `#` before parsing** — `pandas.read_csv(..., sep='\t', comment='#')`, or filter
them out before `csv.DictReader`. The column header is the first non-`#` line.

## ⚠️ Not deposited
The predicted complex structures themselves (tens of thousands of files) are not deposited. Every
number quoted in the prose comes from the tables here, and the scripts in `../scripts/` that
produced them are deposited.
