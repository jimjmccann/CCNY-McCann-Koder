# Structures — what this file is, how it is numbered, and ⛔ the numbering is a trap

⛔⛔ **The one coordinate file here is a MODEL, produced by a generative model. It was not measured.**
Nothing in this directory is an experimental structure and nothing here has been tested at the bench.

---

## `ours_best_medcut51.pdb` — our highest-contact-density first-round design

Generator output, **unmodified** — the file as it came off the run.

MEASURED from the file, for this document:

| | |
|---|---|
| chains | **A** = the designed binder · **B** = the receptor target |
| atoms | **3,765**, all of them `ATOM` records; **0** heteroatoms |
| elements present | C, N, O, S only — ⛔ **no pseudo-atoms**, which this generator can emit |
| residue names | 20, all standard |
| chain A | **96** residues, numbered **1-96**, contiguous |
| chain B | **384** residues, numbered **1-384**, contiguous |
| non-coordinate records | one `HEADER` carrying the campaign-internal design name, and one placeholder `CRYST1`. No `REMARK`, no `SEQRES` |

**Why this design and not one of the other 383.** It has the highest `med_cut` — contact density
through the chain — of all 384 first-round backbones: **51**, against a median of 16 in one arm and
17 in the other. ⛔ **It is therefore the most favourable case in the cohort and not a typical
design.** It is deposited so the single most favourable example of the geometry criticised in
`../methods/README.md` §3c can be inspected, not as a representative sample. ⚠️ It still scores
`min_cut` **0**, like all 384.

---

## ⛔⛔ The trap: chain B is renumbered and the segment breaks are collapsed

This is the most important sentence in this file. **Chain B's residue numbers are NOT the receptor
construct's residue numbers, and the offset between them is not even constant along the chain.**

The generator was given a multi-segment crop of the 409-residue receptor construct and wrote its
output target chain as **1..384 contiguous**, discarding the breaks. MEASURED by aligning chain B's
own sequence back against the construct:

| chain B residues | = construct residues | = UniProt P43220 residues | length | offset to UniProt |
|---|---|---|---|---|
| 1-101 | 6-106 | **29-129** | 101 | +28 |
| 102-304 | 113-315 | **136-338** | 203 | +34 |
| 305-384 | 321-400 | **344-423** | 80 | +39 |

⇒ construct residues **107-112** and **316-320** are absent, their flanks are **covalently joined in
this file**, and the numbering runs straight through the joins. Within the middle segment there is
exactly one substitution, at chain B **226** = UniProt **260**, leucine to phenylalanine — the
deliberate, settled structure-versus-sequence deviation described in `../methods/INPUTS.md` §1.

### ⛔ Why this bites, and why a guard does not catch it

Selecting a named receptor site on chain B by its **construct** residue numbers **returns atoms** —
the numbers 1..384 all exist — so an is-the-selection-empty check passes and the wrong residues get
highlighted and labelled. ⇒ `../scripts/j559_build_baker_pse_20261019.py` **deliberately does not
store** the scene that did exactly that; the file carries the full reason at the point where the site
residue list is defined.

⇒ **Anyone reusing these coordinates must carry the table above with them.** A site, epitope or
contact claim stated in construct or UniProt numbering cannot be applied to chain B directly. ⛔ And
the general form of the lesson: a generator's output is a *product* of the target, never a *source*
for the target's numbering or sequence.

---

## What is NOT in this directory

| | why |
|---|---|
| the two published cryo-EM complexes (`9MM5`, `9MNI`) | somebody else's deposit. `../methods/INPUTS.md` §2 gives the accessions, and the measurement program reads their chain assignments out of their own entity records, so it runs on freshly fetched files |
| the **960 predicted complexes** | not distributed; the per-fold scoring tables are deposited instead, and `../methods/INPUTS.md` §5 states what that costs |
| the receptor on its own, in construct numbering | ⛔ not deposited, and the table above is why that matters: there is no file here in which chain B's numbering is the construct's |
| the other 383 first-round backbones | the gate measurements for all 384 are in `../data/` |
| the viewer session (`.pse`) files | the builders and the reload verifier are deposited in `../scripts/` |

---

## Counts stated here

**1 coordinate file in this directory**, 2 chains, 3,765 atoms, 480 residues in total (96 + 384).
⚠️ Any count in this file must be re-checked if a file is added to or removed from this directory.
