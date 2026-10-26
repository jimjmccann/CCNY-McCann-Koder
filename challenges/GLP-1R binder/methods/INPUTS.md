# Inputs — what is required, what must be obtained upstream, and what is not distributed

Written to be read by someone **outside this group holding only this directory**. A sentence that is
true only because our internal repository exists is a defect, not a shortcut.

---

## 1. The target, and it is recoverable without us shipping anything

**Human GLP-1R.** UniProt **P43220**, residues **24-432 = 409 aa**. That span covers the
extracellular domain, the seven transmembrane helices and the first part of the intracellular tail.
All of it is the organisers' own specification, read off the challenge page.

⭐ **Nothing has to be taken on trust here**: the construct is a plain UniProt range, so a reader can
reconstruct the exact 409-residue sequence from a public record. The published construct sequence was
checked against UniProt residues 24-432 and is byte-identical, with **no visible engineered
mutations**.

⚠️ **But the assay construct is not necessarily that protein.** The challenge page states that the
supplier includes **proprietary thermostabilising mutations without listing them**. In class B1
receptors such mutations characteristically lock a conformational state, and state selectivity is
exactly what this challenge is about. ⇒ every active-versus-inactive statement anywhere in this
directory must be read as *measured on structures, not on the assay construct, whose stabilising
mutations are undisclosed*. ⛔ We do not guess the positions.

⚠️ **One substitution, in the structure rather than in the sequence.** The reference crystal
structure used for geometry carries **phenylalanine at residue 260 where UniProt carries leucine**.
The convention we work to, and which the organisers stated, is **identity from the sequence database,
geometry from the structure**; our structural inputs keep the structure's phenylalanine, which sits
about 29 Å from the nearest designed contact and inside the membrane. ⛔ This is recorded so a reader
who diffs our inputs against UniProt is not left wondering. It is a deliberate, settled deviation,
not a defect.

---

## 2. The two published complexes — ⛔ NOT deposited here, fetch them by accession

`../scripts/j559_baker_reference_20261019.py` measures two third-party cryo-EM structures. They are
**not** included in this directory. Fetch them from the structure database:

| accession | complex | resolution |
|---|---|---|
| **9MM5** | CGRP receptor + RAMP1 + de novo minibinder `dC2_049` | 3.2 Å |
| **9MNI** | CGRP receptor + RAMP1 + de novo minibinder `dC2_050` | 4.1 Å |

⭐ **The program does not need to be told their chain layout.** It reads the chain assignments out of
each file's own `_entity.pdbx_description` records and refuses if a named chain is absent, so the
measurement can be re-run on freshly fetched files without any information from us.

⇒ This matches how the other challenge folders in this repository handle a third-party structure: the
accession is given and the file is not redistributed.

The paper those structures come from is cited in full in `README.md` §3a.

---

## 3. Models and third-party code — obtain from upstream

⛔ **No model weights are distributed**, and not everything needed to reproduce the pipeline is
distributable.

| component | what was used | where it comes from |
|---|---|---|
| backbone generator | **RFdiffusion3** | its own upstream project |
| sequence model | **SolubleMPNN**, checkpoint `solublempnn_v_48_020.pt` | the LigandMPNN project's parameter set |
| sequence-model codebase | **LigandMPNN** | its own upstream project |
| structure predictor | **Boltz-2** (2.2.1) | its own upstream project |
| structural search / TM scoring | **Foldseek**, against a local snapshot of the structure database | its own upstream project |
| geometry and numerics in the deposited programs | **gemmi**, **NumPy**, **SciPy** | their own upstream projects |
| viewer-session builders | **PyMOL** | its own upstream project |

⭐ **SolubleMPNN is the model; LigandMPNN is the codebase.** A command line that invokes LigandMPNN
while loading a SolubleMPNN checkpoint is doing one thing, not two.

⚠️ **The structure-database snapshot used for structural search is a dated local copy**, so any
"nothing resembles this" statement derived from it is a **lower bound** on resemblance, not a
measurement of novelty against the current archive. ⛔ No novelty claim is made anywhere in this
directory in any case.

In-file paths in the deposited programs are **placeholders**, not working defaults. They are spelled
`/path/to/your/working/tree` and carry a comment saying so.

---

## 4. The sequence alignment condition — and it is the point of §5

⛔ **No alignment file is deposited**, and none is needed to read anything here.

All **960** folds of our own designs were run **single-sequence: no alignment on either chain**. The
control folds of the native ligand were run **with** a receptor alignment, which is a confound and is
treated as one. `../README.md` limitation 5 and `README.md` §5 carry the matched control that
resolves it.

⛔ The de novo binder chain is given an **empty** alignment deliberately in every condition: an
alignment is meaningless for a sequence with no natural relatives, and sequence models trained on
natural proteins do not transfer to de novo backbones.

---

## 5. ⛔ What is NOT distributed, and why

| | why |
|---|---|
| the **960 predicted complexes** and the control folds | too large to ship for a cohort this size. The per-fold scoring tables are deposited instead, and they are what every number here rests on |
| the two third-party cryo-EM structures | somebody else's deposit; §2 gives the accessions |
| the generation bundle and the launch layer | the compute-orchestration layer is deliberately absent, as in the other challenge folders. Nothing here requires it: every deposited program reads structures that were already on disk |
| the receptor structural template used by the control arms | third-party derived. Its accession is named in `../scripts/j639_msa_effect_on_baker_channels_20261024.py`, and our own 960 folds use **no** template at all |
| the other 383 of our 384 first-round backbones | the gate measurements for **all 384** are deposited in `../data/`, so the yield is checkable without the coordinates |
| the sequence-generation and fold-gate programs | ⚠️ **a real gap, named rather than hidden**: five of the eight deposited tables carry a header naming the program that wrote them, and those two programs are not here. See `SCRIPTS.md` §3 |
| the superseded output of a corrected program | a chain-order bug's output, kept in the private record as the retraction's evidence. The correction is what is deposited, and the bug is described in full in `README.md` §4c |

⇒ **Consequence, stated plainly:** the numbers in this directory can be **read and audited as
measurements of the predictions**, but they cannot be **re-measured from this directory** without
re-running the prediction step. `README.md` §7 says exactly which claims fall on which side of that
line.

---

## 6. One thing that is NOT reproducible even with every component

⛔ **The backbone generator is unseeded.** Re-running it gives a different cohort, not this one. ⇒ for
the 384 first-round backbones **the artefact is the record**; there is no seed that recovers them. The
sequence step and the prediction step are seeded and are reproducible given the same codebases and
checkpoints.
