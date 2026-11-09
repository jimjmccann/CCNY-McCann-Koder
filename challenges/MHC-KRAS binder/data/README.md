# Data — the scoring tables

Every number quoted in `../README.md` and `../methods/` comes from these files. Each round has a
short human-readable Markdown summary and a per-design table (the file whose name ends
*_per_design.tsv*) holding one row per design per fold.

| file | round |
|---|---|
| `t151_2x2_refold_conjunction_20261108.*` | first re-fold of the 2 × 2 cohort |
| `t152_seqdepth_20261108.*` | sequence depth: ranks 1–5 over the 2 × 2 backbones |
| `t153_mb_depth_raw_20261108.tsv` | raw per-structure scores for the MB-only round |
| `t154_mb_depth_20261109.*` | the MB-only round, 2 shards of new backbones |
| `t155_mb_redraw_20261109.*` | the final re-draw, taking the last candidates to 20 draws |
| `j1064_submission_provenance_20261109.tsv` | maps each submitted row to its design, rank and round |
| `j1063_picks_20261109.json` | the figure cohort input — ⭐ the only place in the figure pipeline where a number is written by hand, and every value in it is asserted against the scored tables before a panel is drawn |
| `j1063_c5_submission_7designs_20261109_build.json` | build manifest for `../structures/j1063_c5_submission_7designs_20261109.pse`: which structure file each panel came from, and the alignment atom counts |
| `j1063_labels_20261109.json` | every string stamped on every figure panel, with its measured width — the record behind the "0 overflowing labels" check |

## Columns in the per-design tables

| column | meaning |
|---|---|
| `tag` | the fold's own label. ⭐ A `REDRAW`-prefixed tag marks a **re-draw** row — part of the 4-seed, 20-draw test rather than the first 5-draw pass. |
| `design` | the internal design name, `arm_shardNN_<arm>_<backbone>_model_<n>` |
| `arm` | which of the four design arms produced the backbone |
| `shard` | the fold shard, which **is** the seed — four shards, four independent seeds |
| `rank` | which designed sequence for that backbone, by sequence-design confidence |
| `cpep` | ⭐ contacts between binder and **peptide** — the **first** limb of the gate |
| `dpep` | ⭐ the **binder↔peptide pairwise interface confidence** (0–1) — the **second** limb of the gate |
| `chla` | contacts between binder and the HLA surface, for contrast |
| `bpl` | the binder's own mean pLDDT (**0–100**, typically 88–97). ⛔ This is a confidence in the binder's *own fold*, ⛔ **not** an interface measure and ⛔ **not** part of the gate. |
| `margin` | `cpep` minus **this shard's own** non-binder floor |
| `imargin` | ⭐ `dpep` minus this shard's own floor — interface-confidence headroom |
| `conj_len` | `conj` **and** the binder length inside the 65–85 aa design window |
| `pass_c`, `pass_i` | the two limbs of the gate, separately |
| `conj` | ⭐ **the gate itself** — `pass_c AND pass_i`. This is the only column that decides anything. |

⚠️ **`dpep` and `bpl` are easy to confuse and are on completely different scales.** The gate uses
`dpep` (0–1). `bpl` is a pLDDT on 0–100 and decides nothing. An earlier draft of this very file
swapped them; it was caught by an independent pre-publication audit, and the mistake is recorded here
rather than quietly corrected, because anyone reading the tables can make it too.

⚠️ **The first round's table has a different schema.** `t151_2x2_refold_conjunction_20261108_per_design.tsv`
carries **13** columns, omitting `design` and `rank`; the other three carry **15**. ⇒ code that reads
all four must not assume a fixed header.

## ⛔ Two traps in reading these tables

1. ⛔⛔ **A design name is not unique across rounds.** `MB_nomotif_broad_11_model_1` exists under
   **both** `arm_shard03_` (the 2 × 2 round) and `arm_shard01_` (the MB-only round) — two different
   backbones behind the same stem — and two other stems each occur under three different prefixes.
   ⇒ **always match the full name including the `arm_shardNN_` prefix, and know which round's table
   you are in.** Matching by stem can silently give you the wrong molecule.
2. ⛔ **Margins are per-shard and are not comparable to a pooled floor.** `margin` and `imargin` are
   already differences against the floor measured **in that shard**, because the shards differ enough
   that pooling smears the comparison. Do not re-derive them against an average.

## ⛔ What a row is not

A row is a comparison between two sets of model outputs. ⛔ It is not a measurement of affinity, and
`conj = True` on a single row is not a result — see `../methods/robustness_and_tiering.md`, where
three of seven designs that cleared on one fold were demoted at four seeds.
