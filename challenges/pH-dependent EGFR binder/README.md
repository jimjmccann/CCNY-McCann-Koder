# pH-dependent EGFR binder

Anthropic × Adaptyv 2026, Challenge 1. Design work ran **2026-09-28 to 2026-10-06**; the
competition's own submission window closed **2026-10-04**. ⭐ Both dates appear in this deposit on
purpose and they are not the same thing — where a document gives a date range it is the **work**,
not the competition window.
The brief: design a de novo binder to the human EGFR extracellular region that binds at
**pH 6.5** and shows **no detectable binding at pH 7.4**, with mouse cross-reactivity.
De novo designs only.

**Target, per the organizers:** the **full ectodomain, tethered** — human Met 1 – Ser 645,
mouse Met 1 – Ser 647. Our design and scoring target is the mature-numbered **621-mer**
(645 − 24 signal residues = 621).

## ⛔⛔ What these designs do and do not address — read this first

| requirement | addressed? |
|---|---|
| binds the human EGFR ectodomain | ⭐ **this is what the campaign optimised** — epitope contact geometry and fold self-consistency on predicted structures |
| pH-dependent binding: **yes** at 6.5, **not detectable** at 7.4 | ⚠️ **addressed by construction in the design, NOT by a scoring term.** See below — this row is the one most easily misread in both directions |
| mouse cross-reactivity | ⛔ **NOT addressed.** No mouse term exists in the chain that scored or proposed these sequences |

### The pH requirement, stated precisely

**The pH dependence is built into the geometry of the interface, not bolted on afterwards.** The
epitope was chosen for it. The binder presents a **carboxylate to one of the receptor's own
histidines** (H409 in mature numbering), and the binder itself carries **zero designed
histidines** — so the titrating group whose pKa sets the pH response is the **target's**, which has
the same measured pKa in every design we submitted, rather than one of ours that we would have had
to predict. The one published pH-dependent result on this target has exactly this polarity: a
carboxylate opposite that histidine gave 13.3× preference at 6.5 over 7.4, while a histidine in
the same position did nothing.

**What we did not do, and will not claim:**

- ⛔ **No pH-dependence term appears anywhere in the chain that scored or ranked these sequences.**
  Measured 2026-10-06. The design intent is in the epitope and the chemistry; it was never
  expressed as a scoring objective, so **nothing here ranks designs by predicted pH selectivity**.
- ⚠️ The only pH that appears numerically in the deposited pipeline is a **net-charge composition
  filter at a single pH (7.4)**, used for solubility and composition. That is not a pH measurement
  and must not be read as one.
- ⚠️ **Our own honest expectation is modest**: of order **12–32×** preference for 6.5 over 7.4, not
  the 100–1000× that an unexamined framing would assume.
- ⛔ **The single largest risk, stated plainly:** no published protein-protein interface has a
  measured histidine pKa up-shift with the partner carboxylate visible in the structure. Our
  strongest calibration for the size of that shift is a **protein-folding** measurement transferred
  to a **binding** one. The design carries a second, independent feature that cannot invert
  precisely because of this.
- ⛔ **None of this has been tested.** No pH-dependent binding has been measured for any molecule
  in this directory.

⚠️ An earlier version of this README said these were "binder candidates, not pH-switch candidates".
**That framing was wrong and has been withdrawn**: it judged the designs against a *conformational*
switch, which was never the goal, and in doing so it understated the design. The brief asks for
pH-**sensitive** binding, which is what the interface chemistry above is for. The accurate
limitation is the narrower one stated in the bullets — no pH term in the scoring chain, and no
measurement.

## ⛔ What no number here is

⛔⛔ **Nothing in this directory is an affinity, and nothing here has been tested
experimentally.** Every value is epitope contact geometry or fold self-consistency computed on
**predicted** structures. `DESIGN_CHOICES.md` §2 states exactly what that does and does not
license, and §5 carries the per-sequence flags.

## What is here

| file | what it is |
|---|---|
| `submission_challenge1_SUBMITTED_20261007.csv` | ⭐ **the 16 sequences actually submitted**, byte-for-byte as uploaded |
| `submission_challenge1_NOVELTY_PASSED_20261007.csv` | the **11** of those 16 that passed the competition's novelty filter |
| `submitted_sequences.fasta` | the **ten** sequences of the 2026-10-06 file, generated from that CSV. ⛔ Superseded — see the two files above and `SEQUENCE_SET_CORRECTION_20261006.md` |
| `submission_challenge1_20261006.csv` | the submission file |
| `SEQUENCE_SET_CORRECTION_20261006.md` | ⛔ which sequences were submitted and how that set changed during the day — read before any other sequence claim |
| `SUBMITTED_SEQUENCES.md` | the deposit record, and what in it is **not** verified |
| `DESIGN_CHOICES.md` | **why** these molecules: the arm structure, the selection, the flags, the negatives |
| `structures/` | all eleven novelty-passed designs in one predicted complex with the receptor (`.pdb` + a PyMOL `.pse` whose eleven designs toggle independently + a one-panel-per-design figure), the same twelve molecules as separate files in `structures/per_molecule/`, the single-design complexes for `id2015` and `id89`, and the receptor model, with provenance |
| `methods/` | the pipeline, every metric defined, every threshold traced, every script described, every input listed |
| `scripts/` | the code that produced and ranked the designs |
| `ideas_that_didnt_work.txt` | what we tried that failed, with the measurement — including things that worked mechanically and were retracted on their meaning |
| `design_guidelines.txt` | the transferable lessons, deliberately **not** EGFR-specific |

## ⭐ Novelty measurement, and an accidental finding

**`accidental_ctpr_homology.txt`** — a short plain-text note, and the thing most worth reading
here if you read nothing else. While checking our own designs' novelty before submitting, we found
that our three 145 aa "bridge" designs are close structural and sequence relatives of the published
**consensus tetratricopeptide repeat (CTPR)** proteins — which are themselves **designed**
proteins, not natural ones. Nothing was copied from any existing protein and no homolog was
searched for at design time; we found it ourselves, by accident, and we are reporting it because it
is useful to anyone running a similar pipeline. The measured follow-up is that sequence novelty and
structural novelty are **separable**: we moved one 30–45 points and the other did not move at all.

**`NOVELTY_20261006.md`** — the index for that material: the reports, the scripts, and the raw hit
tables, per-residue identity maps, PyMOL session and renders under `data/`.

**`NOVELTY_CORRECTION_MONOMER_20261006.md`** — ⛔ **a correction to our own earlier numbers.** The
structural-novelty figures for the single-site designs were first measured on a chain pulled out of
a predicted complex; the organizers fold the submitted sequence **alone**. Re-measured on 693
monomer structures, **0 of 693 draws reaches the "high similarity" threshold**, one design that had
been flagged as high similarity is not, and **not one** nearest-neighbour reference survived the
change. The superseded numbers are left in place with banners rather than deleted.

⚠️ The novelty reports are **working lab reports**. They state their own limitations; read those
sections. In particular the sequence-novelty result is **"no statistically significant hit"**,
which is **not** the same claim as "measured low identity", and we searched **two** of the five
databases the organizers use — so every identity figure is a **lower bound**.

### ⚠️ These are our own novelty measurements, not the competition's

The competition's annotation pipeline runs automatically after upload and reports which sequences
passed its novelty filter. **Nothing in this directory states that outcome.** Everything above is
*our own* measurement, with our own tools and thresholds, anticipating theirs.

## ⛔ Reproducing what is here

⚠️ **Not everything needed to run the pipeline is distributable** — third-party model weights and
third-party scoring code must be obtained from their own upstream projects, and the full set of
raw predicted structures is too large to ship. See `methods/INPUTS.md`. Every required input is
**named with a substitute**, but a bare clone will not run end to end.

⚠️ **Several scripts here were redacted for publication** — machine-absolute paths replaced with
`/PATH/TO/...` placeholders. ⛔ A redacted path is a placeholder, not a working default: set it
before running.

⚠️ **LICENSE.** See the repository root. At the time this directory was prepared there was no
`LICENSE` file, which under default copyright means all rights reserved.
