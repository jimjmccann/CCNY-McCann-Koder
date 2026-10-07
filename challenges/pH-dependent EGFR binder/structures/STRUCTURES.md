# The deposited structures — what each file is, and what it is not

⛔⛔ **Every structure here is a PREDICTION. None has been measured, and none is an affinity.**

## ⛔ Chains — never read one of these files without the chain

Measured on the files themselves:

| file | chain A | chain B |
|---|---|---|
| the two single-design complexes | **the binder** (75 residues) | **the receptor**, 621 residues |
| `EGFR_receptor_WTnumbering_tethered.pdb` | **the receptor only**, 604 residues, WT numbering 4–614 | — |

⛔ The receptor model is **WT-numbered** (4–614, 604 residues present). The complexes' chain B is
the **mature 621-mer** the assay construct uses. These are two numbering systems for the same
protein; `methods/INPUTS.md` states the relationship. Do not compare residue numbers between the
receptor model and any complex (chain B, or chain R of the combined file) without going through
that map.

## The files

| file | what it is | md5 of the source prediction |
|---|---|---|
| `CCNY_challenge1_11designs_20261007.pdb` | **all eleven novelty-passed designs**, one chain each, against one copy of the receptor | built, not copied — see below |
| `CCNY_challenge1_11designs_20261007.pse` | the same, as a PyMOL session: the receptor plus **eleven independently toggleable design objects**, sites highlighted | built, not copied |
| `CCNY_challenge1_grid_20261007.png` | the figure: **one panel per design**, the receptor and both sites in every panel, one shared camera | built, not copied |
| `CCNY_challenge1_overview_20261007.png` | all eleven overlaid in a single view — the same session, rendered without the grid | built, not copied |
| `per_molecule/` | **the same twelve molecules as twelve separate `.pdb` files** — see below | extracted from the combined file |
| `id2015_siteA_complex.pdb` | `egfr-siteA-id2015`, 75 aa, its own unfitted complex | `5ae475747facabde43c358e37bf67e7e` |
| `id89_siteA_complex.pdb` | `egfr-siteA-id89`, 75 aa | `b80926b9bbb5fb0e36c9029cf3879aa7` |
| `EGFR_receptor_WTnumbering_tethered.pdb` | the receptor model the campaign designed against | `ba9d8669e48deb5b9d61a1ace89e25b0` |

⭐ **Every design that passed the organisers' novelty check has coordinates here**, as a chain of
`CCNY_challenge1_11designs_20261007.pdb`. ⛔ Designs that did **not** pass novelty have
no deposited structure. ⚠️ **Absence of a structure is not a statement about a sequence.**

## The combined file — what is in it, and what each chain is

Built from the eleven per-design predicted complexes by a PyMOL script that is **not deposited**
(it carries machine-specific paths and reads inputs that are too large to ship); everything it did
that affects a coordinate is stated in this file. ⛔ The designs were matched to the submitted sequences **by exact chain-A sequence
equality**, never by the id in a filename: the campaign's internal ids are reused for different
molecules, and six of the eleven live in directories named after a different id.

| chain | contents |
|---|---|
| `R` | EGFR extracellular region, tethered, 621 residues. **PREDICTED.** It is the receptor copy from the `egfr-siteA-id2015` complex, so that one binder sits against it **exactly** as predicted (no fit). |
| `X Y Z W V U T S Q P O` | the eleven designs, in the order of the submitted CSV. **PREDICTED.** |

⛔ **There are no glycans in these files.** The receptor is the protein chain only. The campaign
itself folded against a glycosylated receptor — that is a separate matter, documented in
`methods/` and `DESIGN_CHOICES.md`, and nothing here should be read as a statement that the
receptor is unglycosylated in life.

⛔ **Every binder except the reference one was fitted into this frame** (0.33–1.01 Å RMSD on its own
receptor copy), so its contacts with chain R carry that fit error. The exact predicted geometry for
`id2015` and `id89` is in their own complex files here; for the other nine it is in the raw
predictions, which are not distributable at this size.

### `per_molecule/` — twelve files, so they load as twelve named objects

A single `.pdb` loads as **one** object, so the eleven designs arrive as *chains*, not as things you
can switch on and off by name. `per_molecule/` holds the same twelve molecules as twelve files:

| file | molecule | residues |
|---|---|---|
| `EGFR_receptor.pdb` | the receptor, tethered | 621 |
| `d01_id2015.pdb` … `d11_id73.pdb` | the eleven designs, in the order of the submitted CSV | 75, or 80 for `id86` and `id73` |

⭐ **The filenames are the `.pse` object names**, so a design has the same name whether you toggle it
in the session or load its file on its own. No hyphens: PyMOL parses `-` as an operator inside a
selection, so `egfr-siteA-id2015` would be an awkward object name.

⛔⛔ **These are EXTRACTS, not rebuilds.** Each file is a byte-for-byte copy of the ATOM records of
one chain of `CCNY_challenge1_11designs_20261007.pdb`, so **all twelve share one frame** and load
already superposed, matching the deposited figure. Verified per file: every ATOM record is identical
to its source chain, the twelve files together hold **11,490 atoms — exactly the combined file's
count**, and loading all twelve alongside the combined file gives a maximum coordinate difference of
**0.0 Å**. Rebuilding them from the original per-design complexes would have re-introduced each
complex's own frame and silently disagreed with the figure.

⛔ **The superposition caveat travels with every file**, in its own REMARK header: only
`egfr-siteA-id2015` sits against the receptor exactly as predicted — the receptor in this frame is
*that* complex's own copy. Every other design was brought in by superposing its own receptor copy
(0.33–1.01 Å), so its contacts with the receptor carry that fit's error.

⚠️ `per_molecule/` does not replace the combined file, which is what the grid figure is built from.

### The session, and how to use it

The `.pse` is an **overlay**: one receptor object, `EGFR_receptor`, and **eleven separate design
objects** named `d01_id2015` … `d11_id73` in the order of the submitted CSV. ⭐ They are deliberately
**not** grouped, so each design can be switched on and off on its own in the object panel. The named
selections are `site_A` and `site_B`.

Highlighted in the `.pse`: site A (dark red spheres), site B (blue spheres), the receptor as a
semi-transparent grey cartoon, and one distinct colour per design. Two distance dashes
only, both on the reference design, whose geometry is unfitted: `SER12 OG–HIS409 O` (2.77 Å — HIS409
is the site-A histidine the pH switch acts on) and `SER73 OG–GLN384 OE1` (2.95 Å). Each was required
to be a real polar contact in the **unfitted** complex **and** to still be drawn at that length here.
There are no text labels anywhere in the session.

⭐ `CCNY_challenge1_grid_20261007.png` is the same session rendered **one design per panel**, eleven
panels sharing a single camera with the receptor and both sites in every one, because eleven binders
at one site are an unreadable pile when overlaid. The per-panel receptor copies exist only in the
render; the session itself holds one.

## ⛔⛔ How each frame was chosen, and why that makes it flattering

Each complex is **one draw** out of that sequence's set, chosen as its **best by footprint
engagement** — site-A contacts descending, ties broken on the smaller binder-to-footprint
centroid distance. **Every selected frame reaches 11/11 on the site-A footprint** — which is a
statement about the frame that was picked, not about the design: the fraction of that design's draws
reaching 11/11 ranges from 11/69 to 156/180 across the set.

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
2. **Comparability.** In the combined file every binder but the reference one was brought into a
   common frame by superposing its own receptor copy (0.33–1.01 Å); its contacts with chain R
   therefore carry that fit error. Compare contact counts from the tables, never by eye off this file.
3. **Confidence as affinity.** The interface screen's own output says *"SCREEN, NEVER RANK … on
   de novo designs confidence does not correlate with affinity at all."*

## Provenance

The two single-design complexes and the receptor model were copied, unmodified, out of the
campaign's prediction output; the combined file was built from those same predictions, as
described above. The md5 column
above is the md5 of the source file; `../../PROVENANCE.tsv` carries the md5 of every file as
deposited. The full set of raw predictions — tens of thousands of structures — is not
distributable at this size; `methods/INPUTS.md` names every input and its substitute.
