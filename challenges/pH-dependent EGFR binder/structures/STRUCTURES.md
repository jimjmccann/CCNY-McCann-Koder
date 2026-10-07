# The deposited structures — what each file is, and what it is not

⛔⛔ **Every structure here is a PREDICTION. None has been measured, and none is an affinity.**

## ⛔ Chains — never read one of these files without the chain

Measured on the files themselves:

| file | chain A | chain B |
|---|---|---|
| the three complexes | **the binder** (75, 80 or 145 residues) | **the receptor**, 621 residues |
| `EGFR_receptor_WTnumbering_tethered.pdb` | **the receptor only**, 604 residues, WT numbering 4–614 | — |

⛔ The receptor model is **WT-numbered** (4–614, 604 residues present). The complexes' chain B is
the **mature 621-mer** the assay construct uses. These are two numbering systems for the same
protein; `methods/INPUTS.md` states the relationship. Do not compare residue numbers across the
two files without going through that map.

## The files

| file | sequence | md5 of the source prediction |
|---|---|---|
| `id85_bridge_complex.pdb` | `egfr-bridge-id85`, rank 1, 145 aa | `278e92734ff451dc8ccbd9f651005742` |
| `id2015_siteA_complex.pdb` | `egfr-siteA-id2015`, rank 2, 75 aa | `5ae475747facabde43c358e37bf67e7e` |
| `id89_siteA_complex.pdb` | `egfr-siteA-id89`, rank 4, 75 aa | `b80926b9bbb5fb0e36c9029cf3879aa7` |
| `EGFR_receptor_WTnumbering_tethered.pdb` | the receptor model the campaign designed against | `ba9d8669e48deb5b9d61a1ace89e25b0` |

⛔ **There is no structure here for `egfr-siteA-id86` (rank 3) or for ranks 5–10.** Those
sequences were selected from scored tables, and a representative prediction was never extracted
for them. ⚠️ **Absence of a structure is not a statement about a sequence.**

## ⛔⛔ How each frame was chosen, and why that makes it flattering

Each complex is **one draw** out of that sequence's set, chosen as its **best by footprint
engagement** — site-A contacts descending, ties broken on the smaller binder-to-footprint
centroid distance (site A *and* the second footprint, for the bridge). All three selected frames
reach 11/11 on the site-A footprint, and `id85` additionally reaches 7/7 on the second footprint.

⇒ **The frame is the best case by construction. It is not the typical case.** For `id89` the gap
between this frame and its typical draw is the whole point: it is **bimodal** and docks off-site in
a fraction of draws. Read the tables in `DESIGN_CHOICES.md`, never the picture.

⚠️ An earlier selection of these frames was **wrong** and has been corrected: the picker sorted
each sequence's draws on an interface score that is `0.000` in every draw for `id89`, so for that
sequence it had nothing to order on and returned file order. The selection rule above is the fix.

## ⛔ What you cannot measure off these files

1. **Contacts.** The committed scoring tables measured contacts in **each prediction's own
   frame**. If you superpose these files onto a common receptor frame to view them together, the
   residual (~0.5–0.9 Å over ~600 Cα) can move a borderline contact by about 1 Å. Measure in the
   file's own frame, or quote the tables.
2. **Comparability.** `id85` is 145 aa against 75, and the bridge and site-A arms are scored
   against **different footprints**. Their sizes and contact counts are not comparable by eye.
3. **Confidence as affinity.** The interface screen's own output says *"SCREEN, NEVER RANK … on
   de novo designs confidence does not correlate with affinity at all."*

## Provenance

These four files were copied, unmodified, out of the campaign's prediction output. The md5 column
above is the md5 of the source file; `../../PROVENANCE.tsv` carries the md5 of every file as
deposited. The full set of raw predictions — tens of thousands of structures — is not
distributable at this size; `methods/INPUTS.md` names every input and its substitute.
