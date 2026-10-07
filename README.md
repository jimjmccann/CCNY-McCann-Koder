# CCNY-McCann-Koder

Protein design work from the McCann and Koder groups at **The City College of New York**, entered in
the [Anthropic × Adaptyv 2026 Protein Design Competition](https://proteinbase.com/competitions/anthropic-adaptyv-2026).

This repository holds the designs, the reasoning behind them, and the code used to produce them.

## Status

Challenge 1 is deposited: the submitted sequences, the submission file itself, four predicted
complexes, the full method, the scripts that produced and ranked the designs, and a written record of
what the designs do **not** address.

Start at `challenges/pH-dependent EGFR binder/README.md`, which states the limitations before
anything else. `challenges/pH-dependent EGFR binder/SEQUENCE_SET_CORRECTION_20261006.md` is
authoritative on which sequences make up the set, and the deposited
`challenges/pH-dependent EGFR binder/submission_challenge1_20261006.csv` is the authority on its own
contents.

Nothing here has been tested at the bench. No number anywhere in this repository is an affinity.

## Novelty, and an accidental finding

While checking our designs against known proteins before submitting, we found that three of them are
close relatives of the published consensus tetratricopeptide repeat (CTPR) proteins — which are
themselves designed proteins. Nothing was copied, and no homolog was searched for at design time; the
resemblance was found by accident, by us, and is reported because it is useful to anyone running a
similar pipeline. The short note is
`challenges/pH-dependent EGFR binder/accidental_ctpr_homology.txt`; the measurements, scripts and raw
hit tables are indexed in `challenges/pH-dependent EGFR binder/NOVELTY_20261006.md`.

The follow-up is a negative result: we cut sequence identity by 30–45 points without touching a
single interface residue, and the fold did not move. Sequence novelty and structural novelty are
separable.

Our first structural-novelty figures for the single-site designs were measured on the wrong object — a
chain pulled out of a predicted complex, where the organizers fold the submitted sequence alone.
Re-measured on 693 predicted monomer structures, **0 of 693 draws reaches the high structural
similarity threshold**, one design previously flagged as high similarity is not, and no
nearest-neighbour reference survived the change. The superseded numbers are left in place with
correction banners rather than deleted, because what the earlier reasoning rested on is part of the
record. See `challenges/pH-dependent EGFR binder/NOVELTY_CORRECTION_MONOMER_20261006.md`.

These are our own measurements, with our own tools and thresholds, and they are a **lower bound** on
similarity: we searched two of the five databases the competition uses, and the sequence-novelty
result is "no statistically significant hit", which is not the same claim as "measured low identity".

## Competition context

The competition runs five weekly challenges from 2026-09-28 to 2026-10-31, with wet-lab validation by
Adaptyv Bio and results published openly on Proteinbase. Designs selected for testing have their
sequences, structures, methods and measurements published whether or not they work — negative results
included.

## Layout

One directory per weekly challenge. Each challenge gets its own target and its own approach, so
nothing is shared between them by default — no common `src/` that quietly couples one week's method to
another's. Directories are named after the target once it is announced.

```
challenges/
  pH-dependent EGFR binder/     designs, methods and notes for challenge 1
```

The directory name is the challenge title as the organizers published it, so it contains spaces: in a
URL they appear as `%20`, and a shell command that includes the path needs it quoted.

`PROVENANCE.tsv` lists every deposited file with its size and md5 as deposited.

## Scope and licence

This repository is narrower than the group's full research codebase. What is published for each
challenge is decided for that challenge.

There is no `LICENSE` file, so under default copyright all rights are reserved: the material may be
read, but not copied, modified, run or redistributed. Not everything needed to reproduce the pipeline
is distributable in any case — third-party model weights and third-party scoring code must be
obtained from their own upstream projects, and the full set of raw predicted structures is too large
to ship. See `challenges/pH-dependent EGFR binder/methods/INPUTS.md`.
