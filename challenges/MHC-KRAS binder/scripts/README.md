# Scripts

The code that produced everything in `../data`, `../structures` and `../submission_challenge5_20261109.csv`.
These are the working scripts, deposited as they ran, not tidied re-implementations.

## The chain, in pipeline order

| script | stage |
|---|---|
| `s112_build_designability_2x2_bundle.py` | builds the pre-registered 2 × 2 generation spec (held motif × hotspot breadth) and asserts every hotspot atom against the input |
| `s113_upload_designability_2x2_bundle.sh` | stages that bundle to object storage for the GPU workers |
| `s114_designability_2x2.py` | the design-stage gates: continuity first, then clash survival |
| `s115_select_refold_cohort.py` | samples the re-fold cohort **one design per diffusion batch**, seeded |
| `s116_make_2x2_cofold_inputs.py` | co-fold inputs for the 2 × 2 cohort |
| `s117_shard_2x2_cofold.py` | splits folds across shards and replicates the floor controls into each |
| `s118_score_2x2_cofold_per_shard.py` | scores each shard against **its own** floor |
| `s119_make_seqdepth_cofold_inputs.py` | co-fold inputs for sequence ranks beyond rank 1 |
| `s120_score_seqdepth.py` | scores the sequence-depth round |
| `s121_build_mb_depth_bundle.py` | the MB-only round's spec — **copies** the winning cell byte-for-byte and asserts its hash rather than regenerating it |
| `s122_upload_mb_depth_bundle.sh` | stages that bundle |
| `s123_make_mb_depth_cofold_inputs.py` | co-fold inputs for the MB round, and ⭐ derives re-draw targets so a re-draw re-folds **the rank that cleared**, not rank 1 |
| `s124_score_mb_depth.py` | scores the MB round |
| `s125_build_mb_depth_pse.py` | builds the PyMOL sessions in `../structures` |
| `j1064_build_submission_csv.py` | builds and validates the submission file |

Supporting, imported by the above: `s59_make_r4_cofold_inputs.py`, `s60_score_r4_cofold.py`,
`s71_make_5b_cofold_inputs.py`, `s75_score_scaffolded_arm.py`,
`s88_cont_pep_excluding_pins.py` (the consistent-atom-set contact metric),
`s89_pool_gate_plus_prod.py`.

## ⭐ How to read these scripts

They are written to **refuse rather than guess**, and most of their length is assertions. Several
assert the *text of a sibling script* before relying on its behaviour, so that if the sibling changes
the dependent script stops instead of silently diverging. `j1064` is the clearest example: it derives
the submitted cohort from the scoring tables, resolves each design to exactly one generation round by
lookup, refuses if a name resolves to two, and re-derives each sequence the same way the folded one
was obtained.

⛔ **Nothing here re-ranks or re-scores at submission time.** The submission builder reads the
already-scored tables; it cannot promote a design.

## ⛔ Notes for anyone trying to run them

- **Paths written `/PATH/TO/...` are placeholders, ⛔ not working defaults.** Object-store locations
  and credential file paths were replaced with placeholders for publication; the upload scripts will
  not run as deposited, by design.
- Third-party model weights and third-party scoring code are **not** redistributed here and must be
  obtained from their own upstream projects.
- `s125` deliberately has **no `__main__` guard**: under the headless viewer that runs it the module
  is not `__main__`, so a guarded script would exit successfully having written nothing.
- Comment markers like `[[some-slug]]` are pointers into the group's internal engineering-lessons
  notes. ⚠️ They are **not links and resolve to nothing here**; they are left in place because each
  one marks a specific trap the surrounding code exists to avoid.
- `[M]` in a comment marks a measured value and `[R]` a derived or estimated one.

## ⚠️ One deposited script carries a known false label

`s60_score_r4_cofold.py` hard-codes a provenance string naming an adopted cohort and a seed that are
**both wrong for the rounds deposited here**. It is included because later scripts import it, and the
label is called out here rather than quietly corrected, because the per-shard outputs it writes are
scratch and should not be cited. ⛔ The **sequences and scores** in `../data` are unaffected — those
are asserted per design against their own source. This is a false *label* in a tool, not an error in
the results.
