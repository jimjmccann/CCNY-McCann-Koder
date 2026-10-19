# Inputs — what is required, and what is not distributed

Written to be read by someone **outside this group holding only this directory**. A sentence that is
true only because our internal repository exists is a defect, not a shortcut.

---

## 1. The target, and it is recoverable from this deposit

**Human TNF-α**, soluble trimer. UniProt **P01375-1**, target residues **77–233 = 157 aa**,
structure reference **1TNF** chains A–C. All of that is the organizers' own specification, read off
the challenge page.

The 157-residue human sequence as published by the organizers:

    VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTI
    SRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL

⭐ **The target as actually used is inside every deposited design file** — chains B, C and D, 152
residues each, renumbered 1..152. MEASURED for this document: that sequence is **byte-identical to
residues 6–157 of the organizers' own human sequence**, in all three chains. ⇒ nobody has to trust
a description of the target; it can be read out of the coordinates.

⚠️ **The reference structure was repaired before use, and the repair is already in the deposited
coordinates.** 1TNF carries **LEU at residue 143 where the challenge sequence carries ASP**, in all
three chains, *inside the epitope* — and the deposited mutation annotation for 1TNF is empty, so a
routine "is this wild type?" check passes it. The aspartate was transplanted from a wild-type
structure by a local four-atom backbone fit. ⛔ **The repair script, the repair log and the standalone
repaired target coordinates are NOT in this deposit.** The sequence assertion above is the check that
matters, and it can be re-run against the deposited files.

⚠️ **Mouse TNF-α** (UniProt **P06804**, 156 aa, reference **2TNF**) is part of the challenge and is
**not an input to this cohort**. No mouse structure, sequence or scoring appears here.

---

## 2. Models and third-party code — obtain from upstream

⛔ **No model weights are distributed**, and not everything needed to reproduce the pipeline is
distributable.

| component | what we used | where it comes from |
|---|---|---|
| backbone generator | **RFdiffusion3** | its own upstream project |
| sequence model | **SolubleMPNN**, checkpoint `solublempnn_v_48_020.pt` | the LigandMPNN project's parameter set |
| sequence-model codebase | **LigandMPNN** | its own upstream project |
| structure predictor | **Boltz-2** | its own upstream project |
| session figures | **PyMOL** | its own upstream project |

⭐ **SolubleMPNN is the model; LigandMPNN is the codebase.** A command line that invokes LigandMPNN
while loading a SolubleMPNN checkpoint is doing one thing, not two.

The scripts expect `LIGANDMPNN_DIR` and `MPNN_PY` as environment variables; their in-file defaults
are **placeholders**, not working paths.

---

## 3. The target alignment — not distributed

The predictor is given a real multiple-sequence alignment for the target chains: **2198 rows**. The
binder is given an **empty** alignment deliberately — it is a de novo sequence, an alignment is
meaningless for it, and sequence models trained on natural proteins do not transfer to de novo
backbones.

⛔ **The alignment file is not deposited.** It is third-party derived and regenerable from the target
sequence with any standard search.

---

## 4. Generation inputs that are NOT distributed

| | why |
|---|---|
| the generation bundle itself (`v46`) | contains the runner and launch configuration; the launch layer is deliberately excluded. `scripts/j484_…py` is deposited and shows exactly what the bundle contains |
| the four guidepost residues as a standalone coordinate file | ⭐ **not needed**: the pins are present *in* every deposited design file, flagged in the `_atom_site.is_motif_atom_with_fixed_seq` column, and their per-design residue numbers are in `../data/mpnn_20261018/fixed_residues_multi.json` and in the scoring tables' `pins` column |
| the 373 designs that did **not** pass the two backbone gates | the gate measurements for **all 384** are deposited in `../data/backbone_gates_20261018/`, so the yield is checkable without the coordinates |
| the generator's per-design statistics tensors | large binary artefacts, no bearing on any stated number |

---

## 5. Prediction outputs that are NOT distributed

⛔ **The predicted complex structures are not deposited.** All **110 of 110** are scored, at
roughly 35 MB for the first round alone. What is deposited is the
**full scoring table**, 103 columns per prediction, which is what every number in `../README.md`
rests on.

⇒ **Consequence, stated plainly:** the placement and bridge numbers in this deposit can be *read and
audited as measurements*, but they cannot be *re-measured from the deposit* without re-running the
prediction step. The design side **can** be re-measured, because the designed coordinates are here.

---

## 6. What a reader can and cannot reproduce from this directory alone

**Can:**
- read the pins out of the designed coordinates and confirm there are 4 per design (3 ARG + 1 HIS);
- confirm the designed Arg→Glu107 distance is 3.45 Å, by measuring it;
- confirm the target is the challenge's own sequence;
- re-derive both backbone gates and the JOINT count of 11 from the deposited gate tables;
- confirm every designed sequence is cysteine-free, and that there are 110 of them;
- read every per-chain-pair confidence number and check the global-`iptm` confound for themselves.

**Cannot, without obtaining third-party components and re-running:**
- reproduce the 384 backbones — ⛔ the generator is **unseeded**, so re-running the bundle gives a
  different cohort, not this one;
- reproduce the 110 sequences (seeded at 37, so reproducible given the same codebase and checkpoint);
- reproduce the predictions (seeded at 20261018, likewise).

---

## 7. One observed consistency worth recording

MEASURED across all 11 deposited FASTA files: the generator-written header field `num_res` is
**exactly 4 fewer** than the designed chain's actual length, in every one (for example 192 against
196). Four is the pin count in every design. ⇒ the two independent records of "how many residues
were held" agree. It is recorded because it is the kind of cheap cross-check that catches a
pin set silently failing to be passed through.
