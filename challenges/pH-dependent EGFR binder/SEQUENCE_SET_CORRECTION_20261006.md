# ⛔⛔ Which sequences were submitted — and how that set changed

Read this before any other sequence claim in this directory.

**Where this stands: we submitted designs, some did not pass the competition's novelty filter, and
we are waiting for results.**

⛔ **The final submission was 16 sequences, of which 11 passed the novelty filter.**
⭐ **Both files are deposited here**, so every count in this document is checkable:

| file | n | what it is |
|---|---|---|
| **`submission_challenge1_SUBMITTED_20261007.csv`** | **16** | ⭐ **the file actually submitted**, byte-for-byte (md5 `a2bd128ab18dc69a80495a997470657a`) |
| **`submission_challenge1_NOVELTY_PASSED_20261007.csv`** | **11** | the subset that passed the competition's novelty filter |
| `submission_challenge1_20261006.csv` | 10 | ⛔ **superseded.** The file prepared on 2026-10-06, an earlier state of the set — kept as the record of that day |
| `submitted_sequences.fasta` | 10 | generated from the 10-row file above, and therefore also superseded |

The set changed between 2026-10-06 and submission on 2026-10-07.

| set | n | members |
|---|---|---|
| **submitted** | **16** | the 11 below, plus `egfr-bridge-id85`, `egfr-bridge-id1104`, `egfr-bridge-id47`, `egfr-bridge-id85-mutB`, `egfr-siteA-id19sh09` |
| **passed novelty** | **11** | `egfr-siteA-id2015` `id86` `id89` `id3035` `id3036` `id2046` `id3020` `id2048` `id3017` `id2` `id73` |
| **did not pass** | **5** | all four bridges and `egfr-siteA-id19sh09` |

The 2026-10-06 ten-sequence file described below was a different, earlier set, and one sequence
that appears in earlier files than that was **deliberately dropped**.

## The submitted set

Authoritative file: `submission_challenge1_20261006.csv` — deposited here byte-for-byte as built.
⛔ **Nothing had been sent at the time this was written.** The CSV is the file prepared for
upload, deposited unaltered; it is not a record of a completed submission. ⚠️ Treat the CSV as the
authority on its own contents.
Its builder is internal to the campaign and is **not** deposited (it reads campaign directories
that are not part of this deposit), so the authority rests on the deposited CSV itself, not on a
script you cannot see. `submitted_sequences.fasta` in this directory is **generated from that
CSV** by `scripts/make_deposit_fasta.py`, which **is** deposited and re-asserts the organizer's
bounds (10–250 aa, ≤ 20 per submission, standard residues only) per sequence before writing —
so you can re-run it against the deposited CSV and check that the two cannot disagree.

| rank | name | aa | arm | the measurement behind it |
|---|---|---|---|---|
| 1 | `egfr-bridge-id85` | 145 | bridge | engages **both** footprints in the same draw, `full_both` 20/24 |
| 2 | `egfr-siteA-id2015` | 75 | site A | simultaneous occupancy **0.857** over 70 draws — the best measured in the programme |
| 3 | `egfr-siteA-id86` | 80 | site A | a **different backbone** (shard14_35); simultaneous occupancy 0.619, n = 21 |
| 4 | `egfr-siteA-id89` | 75 | site A | simultaneous occupancy 0.484, pooled n = 31 |
| 5 | `egfr-siteA-id3035` | 75 | site A | tied leader, 0.850 |
| 6 | `egfr-siteA-id3036` | 75 | site A | tied leader, 0.844 |
| 7 | `egfr-siteA-id2046` | 75 | site A | tied leader, 0.838 |
| 8 | `egfr-siteA-id3020` | 75 | site A | tied leader, 0.812 |
| 9 | `egfr-siteA-id2048` | 75 | site A | tied leader, 0.812 |
| 10 | `egfr-siteA-id3017` | 75 | site A | tied leader, 0.800 |

Ranks and notes above are quoted from the ranked plan in the undeposited submission builder
that produced the CSV. Lengths are measured from the CSV.

⚠️ **Ranks 5–10 are one molecule's neighbourhood, not six independent bets.** They are the tied
leaders of the same site-A family and differ from one another at a handful of positions. The
genuinely independent bets in the set are ranks 1–4: a 145-mer bridge, the best-measured site-A
sequence, a **different** site-A backbone, and a fourth site-A sequence. Treat the set as
**four independent designs plus a replicate band**, and read `DESIGN_CHOICES.md` §4 on why
farthest-point selection exhausts site A after three.

## ⛔ Two changes from the earlier files in this directory, and why

**1. `id15` was dropped.** `SUBMITTED_SEQUENCES.md` lists `id15` at rank 3 as a diversity
bet, as did an earlier four-sequence FASTA that is not deposited. It is **not in the submission.** It was dropped on the
measurement, not on taste — four independent signals, all pointing the same way:

- simultaneous site-A occupancy **0.178** in the 90-draw deep round;
- interface-screen score **0.000** in 10/10 draws;
- interface packing **0.000** in 10/10 draws, and backbone τ 4.29× / clashes 1.98× its own
  crystal control;
- ⛔ **2 of its 10 draws dock ~29 Å off-site entirely**, with the nearest receptor residue at
  position 103 — nowhere near the designed epitope.

⇒ A sequence measured as not binding the designed site is not diversity. ⛔ Where
`SUBMITTED_SEQUENCES.md` or `DESIGN_CHOICES.md` still describe a **four**-sequence set with
`id15` in it, **this file is the correction** and those passages describe the set as it stood
earlier the same day. They are kept rather than rewritten because the reasoning in them about
ranks 1, 2 and 4 is unchanged and still load-bearing.

**2. `id86` was added.** It is the only submitted site-A sequence on a **different backbone**
(shard14_35), which is the only backbone-level diversity in the site-A half of the set. Its
sequence is carried as a literal in the submission builder with its provenance recorded there,
because it exists only in that arm's input file and was not re-parsed.

## What this means for the structures

`structures/` holds one combined predicted complex containing **every** design that passed the
organisers' novelty check, each as its own chain, against the receptor — plus
single-design complexes for `id2015` and `id89`. `structures/STRUCTURES.md` says exactly which
prediction each chain and each file is, and how each frame was chosen. ⛔ Designs that did **not**
pass novelty, the bridge `id85` among them, have no deposited structure; do not read the absence of
a structure as a statement about a sequence.
