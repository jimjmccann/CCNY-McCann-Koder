# Cohorts — what each deposited table is, and what it is NOT comparable to

⛔⛔ **Every number in every table here is a model's output, with one exception**:
`j559_baker_reference_20261019.tsv` is a geometric measurement of two third-party **experimental**
structures. Everything else is a structure predictor's or a generator's output about its own work.
⛔ Nothing in this directory was measured on a protein.

Eight tables, in three groups. Row counts below are data rows, excluding the `#` header comments.

---

## Group 1 — the measurement on experimental structures (strand A)

| file | rows | what it is |
|---|---|---|
| `j559_baker_reference_20261019.tsv` | **2** | the geometry of the two published de novo minibinders, measured with our own programs on their deposited cryo-EM complexes. 17 columns: length, `min_cut` / `med_cut`, radius of gyration and `rg_ratio`, the helix signature, the end-to-end distance, both interface residue sets, and the straight-line fit to the contacting binder Cα positions |

⭐ This table and the two in Group 3 are the two columns of `../methods/README.md` §3c, measured with
the **same code paths** so they are comparable by construction.

---

## Group 2 — the published filter and its controls (strand B)

| file | rows | what it is |
|---|---|---|
| `j637_baker_filter_c3_20261024_v2.tsv` | **987** | per-fold channel values and per-channel pass flags. **960** rows are our designs (96 designs × 10 seeds) and **27** are known-binder controls, split across three control cohorts named in the `cohort` column |
| `j639_msa_effect_20261024.tsv` | **42** | the matched single-factor control: the native ligand against the same receptor, **21** draws with a receptor alignment and **21** without, same template throughout |

⛔⛔ **Three things about these two tables that change how they read.**

1. **`iptm_pair`, not `iptm_global`.** Both columns are present. The global value averages every chain
   pair; the pairwise one is the number that means *binder against target*. On this construct the
   target is a single chain so the two nearly coincide — but the column to use is the pairwise one,
   and it is kept beside the global one so that is checkable rather than asserted.
2. **Both readings of the ambiguous interface-error criterion are present** — `pae_iface_min` (the
   literal minimum element over the binder × target block) and `pae_iface_mean` (the mean over it),
   with a pass flag for each. The source paper prints the criterion without defining the reduction.
   `../methods/README.md` §4d settles which is defensible here, **on the controls**.
3. ⛔⛔ **`n` is 96 designs, not 960 folds.** The ten shards are ten seeds of the same 96 designs, so
   a per-fold pass rate counts each design ten times. Rank on the per-design best-of-ten.

⚠️ **The `binder_first` column is not cosmetic.** It records which chain came first in that cohort's
own input, because it is **not the same across cohorts** — our design inputs are binder-first and
both control inputs are receptor-first. An earlier version of the scoring program inferred this from
position and measured the wrong block for all 27 controls. `../methods/README.md` §4c is the full
account. The column exists so the correction is visible in the data.

⚠️ **The 27 controls were folded WITH a receptor alignment and all 960 of our designs WITHOUT one.**
That is a confound, not a like-for-like comparison, and it is why `j639_msa_effect_20261024.tsv`
exists. ⛔ Do not quote the control-versus-design gap from the 987-row table alone.

---

## Group 3 — our own design cohorts, deposited as the denominator for two claims

⛔ These five tables are **not** known-binder data. They are our own designs, and they are here
because two claims in `../methods/README.md` are untraceable without them.

| file | rows | what it is |
|---|---|---|
| `j557_topology_deep6_20261019.tsv` | **288** | the topology and compactness gate on the `deep6` arm of the 384 first-round backbones |
| `j557_topology_shallow4_20261019.tsv` | **96** | the same gate on the `shallow4` arm. 288 + 96 = **384** |
| `glp1r_fold_novelty_sel20_20261016.tsv` | **160** | each designed sequence folded **alone**, with the fold-back column `selfconsist_tm_vs_backbone` |
| `j558_fold_unbiased_20261019.tsv` | **96** | the same instrument on the unbiased sequence arm |
| `j558_fold_biased_20261019.tsv` | **96** | the same instrument on the composition-biased sequence arm. 160 + 96 + 96 = **352** |

⭐ **The pooling is the point, and it is arithmetic a reader can check.** The `0 of 384` `min_cut`
result pools the two topology tables; the `0 of 352` fold-back result pools the three fold tables.
`../methods/README.md` quotes both denominators without naming the tables inline, so they are all
deposited together.

⛔⛔ **Within the topology tables, the `c_*` continuity columns are a PREFILTER and the tables say so
in their own header.** A design with `c_continuous=0` has a broken chain; its clearance (`p_*`) and
topology (`t_*`) numbers are **not comparable** to a continuous design's, because a detached element
scores *better* on clearance and zero on contact. ⇒ rank within `c_continuous=1` only.

⛔ **The fold tables' `plddt` is single-sequence pLDDT and is not comparable to any alignment-mode
figure**, as their header states. And pLDDT cannot call a conformational state.

⚠️ **These three fold tables carry the designed sequences themselves**, in the `sequence` column:
352 sequences, 75-104 residues, MEASURED cysteine-free in all 352.

⚠️ **They also carry raw structural-search output** — `best_target`, `fident`, `qcov`,
`qtmscore_FOLDED`, `evalue` — against a dated snapshot of the structure database. ⛔ No novelty
verdict is drawn from it, for the two reasons the tables' own header comments give. Treat those
columns as raw search output, not as a novelty result.

---

## What is NOT in this directory

| | why |
|---|---|
| the **960 predicted complexes** themselves | not distributed; `../methods/INPUTS.md` §5 states the consequence |
| the other 383 first-round backbones' coordinates | the gate measurements for all 384 **are** here |
| the superseded output of the corrected scoring program | it is the chain-order bug's output, and it measured an intra-receptor block for all 27 controls. It is kept in the private record as the retraction's evidence; the correction, the `_v2` table, is what is deposited, and the bug is described in `../methods/README.md` §4c |
| any apo (binder-alone, in-complex-numbering) prediction | ⛔ it does not exist for this cohort, which is why only three of the filter's five criteria are evaluated |
