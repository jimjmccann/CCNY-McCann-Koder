# Methods — how the Challenge 5 designs were made, scored and tiered

| file | what it covers |
|---|---|
| `design_and_scoring.md` | the target, the four design arms, sequence design, co-folding, and the scoring gate |
| `robustness_and_tiering.md` | seeds and draws, what "validated" is allowed to mean, and the re-draw finding |

## The target

A peptide–MHC complex. In the **co-folded complexes** — which is what the deposited sessions and
figures show — the four chains are:

| chain | what | residues |
|---|---|---|
| A | HLA-A\*11:01 heavy chain | 276 |
| B | β2-microglobulin | 100 |
| C | the KRAS G12D 10-mer peptide | 10 |
| D | **the designed binder** | 65–85 |

⛔⛔ **TWO CONVENTIONS EXIST IN THIS CAMPAIGN AND THEY DISAGREE. ALWAYS NAME THE FILE WITH THE
CHAIN.** `[M]` Measured on both file types:

| | binder | HLA | β2m | peptide |
|---|---|---|---|---|
| **co-fold output** (sessions, figures, scoring) | **D** | A (276) | B (**100**) | C (10) |
| **design-stage output** (sequence files, generator PDBs) | **A** | B (276) | C (**99**) | D (10) |

⇒ ⭐ Two things change between them: the **letters rotate**, and **β2m is 100 residues in the complex
but 99 in the design-stage files**, because the fold input omits its leading residue. Both numbers are
correct in their own frame — ⛔ neither is a typo, and code that mixes them silently compares the
binder to the peptide. `scripts/s123_make_mb_depth_cofold_inputs.py` compares against β2m **minus its
first residue** for exactly this reason, and `scripts/j1064_build_submission_csv.py` asserts the
design-stage lengths (276 / 99 / 10) rather than the complex ones.

⛔ **The scoring problem is the design problem.** Nearly the whole groove is shared with every other
peptide this allele presents, so a binder that grips the MHC and ignores the peptide scores well on
most structural metrics and is worthless. Every gate below is built to refuse that molecule.

## The pipeline, in order

1. **Backbone generation** — RFdiffusion3 against the full four-chain complex, on rented GPUs.
2. **Sequence design** — SolubleMPNN over each generated backbone, several sequences per backbone,
   ⛔ with cysteine excluded outright (`--omit_AA C`).
3. **Co-folding** — each designed sequence re-folded *as a complex* with the unchanged target, and
   scored only on the resulting binder↔peptide interface.
4. **Gating** — the pre-registered conjunction in `design_and_scoring.md`.
5. **Re-draw** — every design that cleared is re-folded across four independent seeds before it is
   allowed to be called anything. See `robustness_and_tiering.md`.

## ⭐ The one experiment that defined the generation stage

A **2 × 2, pre-registered, in which exactly two things vary and nothing else does**:
**held motif × hotspot breadth.** Binder length (65–85) and recycling (3) were identical in all four
cells, so a difference between cells is attributable.

| arm | held motif | hotspot breadth |
|---|---|---|
| `R_motif_narrow` | a CDR3β-derived motif held in the binder chain | narrow |
| `B_motif_broad` | the same held motif | broad |
| `M_nomotif_narrow` | ⛔ none — binder built from nothing | narrow |
| `MB_nomotif_broad` | ⛔ none | broad |

**Hotspot breadth** is how much of the peptide the generator is told to aim at:

- **narrow** — 2 atoms on a single peptide residue. This was the *entire* steering signal of the
  preceding campaign.
- **broad** — 6 atoms spread across four peptide residues (backbone O of two, three side-chain
  carbons of a third, backbone N of a fourth). ⛔ The broad set was **re-derived** against this
  target's own measured anchor slots rather than copied from the earlier arm, so it is auditable
  from `scripts/s112_build_designability_2x2_bundle.py` alone.

⇒ ⭐⭐ **Every design that survived to submission came from a *no-motif* arm.** Three of the four
validated designs are `MB_nomotif_broad` and one is `M_nomotif_narrow`; the two held-motif arms
contributed nothing that cleared at 20 draws. Grafting a known recognition motif into the binder did
not help here.

## ⛔ Two counting rules that change what the numbers mean

- **`n` is BATCHES, not designs.** RFdiffusion3 draws one contig length per diffusion batch, so
  designs inside a batch are **not independent**. Anywhere a count of designs is reported, the batch
  structure is reported with it, and the re-fold cohorts sample one design per batch.
- **A metric computed over a different atom set cannot compare arms.** The held-motif arms hold
  binder residues fixed; counting those held residues as achievements makes the arms incomparable.
  Peptide-contact counts are therefore computed over a consistent atom set
  (`scripts/s88_cont_pep_excluding_pins.py`).

## ⛔ What is not here, and cannot be

Third-party model weights and third-party scoring code are not ours to redistribute and must be
obtained from their own upstream projects. The full set of raw predicted structures is far too large
to ship; what is deposited is the scoring tables, the sessions for the submitted cohort, and the code
that produced both. Paths written as `/PATH/TO/...` in the scripts are placeholders, ⛔ **not working
defaults**.
