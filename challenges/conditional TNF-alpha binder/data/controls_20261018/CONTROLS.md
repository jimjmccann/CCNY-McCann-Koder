# The positive and negative controls — what was run, what it shows, and what it does NOT show

This directory holds the evidence for the claim made in the challenge `README.md` and in
`methods/README.md`: that the structure predictor used throughout this deposit is **not blind on a
de novo interface**. Everything quoted below is in `controls_scored.tsv` (2 rows).

## The question these controls answer
Our designs' interface is not reproduced by the predictor (0 of 110 on the `v46 joint11` cohort;
0 of 88 on the Baker five criteria for the submitted C1 cohort). There are two explanations and
they are not the same: **our designs are wrong**, or **the predictor cannot see a de novo interface
at all**. A control that is itself a *de novo* binder with *measured* affinity separates them.

## The positive control — LCB1
**LCB1** (Cao *et al.* 2020) is a computationally designed miniprotein binder against the
SARS-CoV-2 receptor-binding domain, with measured picomolar-range affinity and a solved complex,
**PDB `7JZU`**. ⭐ It is the right control precisely because it is **de novo** — it has no natural
homologs, so the predictor cannot be succeeding by recall of a protein family.

It was pushed through the **identical path** used for our own designs: binder sequence on an
**empty MSA**, no templates, same predictor, same settings.

| | value |
|---|---|
| binder length | 55 aa |
| CA RMSD of the predicted binder to its crystal pose | ⭐ **0.89 Å** |
| crystal epitope residues recovered | ⭐ **30 of 31** (recall 0.968) |
| interface ipTM | **0.955** |
| complex pLDDT | 0.946 |
| target superposition RMSD (the frame the binder RMSD is measured in) | **0.59 Å** |

⇒ on a de novo interface, through this path, the predictor places the binder essentially correctly.

## The negative control — TNFα protomer vs the same RBD
A well-folded natural protein with **no RBD affinity**, substituted for LCB1 against the same
target. It folds fine (pLDDT 0.831) and produces **no confident interface**: ipTM **0.161** against
0.955, complex iPDE **11.4** against 0.57.

### ⛔⛔ The negative does NOT fail on every metric — and this matters
| metric | positive (LCB1) | negative (TNFα) | discriminates? |
|---|---|---|---|
| interface ipTM | **0.955** | 0.161 | ✅ yes, strongly |
| complex iPDE | 0.57 | 11.4 | ✅ yes |
| **epitope enrichment** | 5.66× | ⛔ **5.78×** | ⛔ **NO — the non-binder scored HIGHER** |

⇒ **"the negative fails cleanly" is true only of the confidence metrics.** On epitope enrichment
the non-binder outscored the real binder, which is why **epitope enrichment was removed as a
discriminator** in this project and is not used to rank anything in this deposit. ⛔ Do not read
"fails cleanly" as "fails on everything".

## ⛔ What these controls do NOT establish
1. ⛔⛔ **`7JZU` was published in 2020 and may sit in the predictor's training data.** The control
   therefore cannot fully separate *generalisation* from *recall*. This is the single largest
   surviving caveat and it is not resolved by anything in this deposit.
2. ⛔ **n = 1 per arm.** Both arms are single folds with one seed; no spread is measured.
3. ⛔ **The arms are not length-matched** — 55 aa binder vs a 152 aa negative.
4. ⛔ A control establishes that the *instrument* responds. It says nothing about whether **our**
   designs bind. **No prediction, positive control or not, is an affinity measurement.**

## Files here
| file | what it is |
|---|---|
| `controls_scored.tsv` | both arms, 2 rows — the source of every number above |
| `pos_lcb1_denovo_bound.fasta` | positive input: LCB1 (55 aa) + RBD (193 aa), both `empty` MSA |
| `neg_tnfa_vs_rbd_bound.fasta` | negative input: TNFα protomer (152 aa) + the same RBD |
| `confidence_*_model_0.json` | the predictor's own confidence output, per arm |
| `../../structures/controls_20261018/*_model_0.pdb` | the two predicted complexes |
| `../../scripts/j505_make_denovo_control_inputs_20261018.py` | builds the inputs |
| `../../scripts/j506_score_denovo_control_20261018.py` | computes the RMSD and the epitope recovery |

⭐ `j506` carries a **built-in self-validation**: scoring the crystal against its own reference
epitope must return recall 1.000, and the script asserts it.

⚠️ **The crystal reference `7JZU` is not redistributed here.** It is an unmodified RCSB download —
fetch it with `curl https://files.rcsb.org/download/7JZU.pdb`. The reference epitope is the 31 RBD
residues within 5.0 Å (heavy atom) of LCB1 in that structure.
