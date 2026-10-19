# Structures — what each file is, which frame it is in, and whether it was modified

⛔⛔ **Every coordinate here is a MODEL, produced by a generative model. None of it was measured.**
Nothing in this directory is an experimental structure, and nothing here has been tested at the
bench.

---

## `designs_v46_joint11/` — the 11 designed backbones

Eleven gzipped mmCIF files, generator output from campaign `v46`, **unmodified** — these are the
files the generator wrote, byte-for-byte as they came off the run, which is why they still carry the
generator's own non-standard annotation column.

MEASURED from the files, for this document:

| | |
|---|---|
| chains | **A** = designed binder · **B, C, D** = the three target protomers |
| binder length | **187–204 residues**, 9 distinct values across the 11 |
| target | **152 residues per chain**, renumbered 1..152 |
| guideposts | **4 per design = 3×ARG + 1×HIS**, at per-design positions |
| designed Arg→Glu107 | **3.45 Å in 11 of 11** (pinned arginine guanidinium N to Glu107 carboxylate O, heavy atoms, best of three protomers) |

**Why these 11 and not the other 373.** They are the designs that passed *both* backbone gates —
continuity (`max_step ≤ 4.2 Å`) **and** topology (`min_cut ≥ 19`). ⛔ **That is the top 2.9% of the
run and is not a typical design.** The gate measurements for all 384 are in
`../data/backbone_gates_20261018/`, so the yield can be checked without the other coordinates.

### ⭐ The guideposts are IN these files, and that is deliberate

The generator records a held residue **only** in a non-standard mmCIF column,
`_atom_site.is_motif_atom_with_fixed_seq`. There is no REMARK, HETATM or segID record of it. These
files retain that column, so:

- the four pinned residues of each design can be read out of the deposit with no external
  information — their numbers differ per design, so they could not be hardcoded anyway;
- the designed Arg→Glu107 distance can be **re-measured**, not merely read from our table;
- any tool that strips non-standard columns, or converts these to PDB format, **loses the pin
  record**. Keep the mmCIF.

Per-design pin residue numbers are also written out in `../data/mpnn_20261018/fixed_residues_multi.json`
and in the `pins` column of the scoring tables, so they are recorded three independent ways.

### The frame

All 11 share a **common target frame** — the target was held rigid across the whole campaign. That
is what makes the designed binder positions directly comparable to one another, and it is what a
figure overlaying them relies on. ⚠️ The **predicted** structures are in their own frames and were
superposed onto the target before any placement number was taken; see `../methods/README.md` §7.

---

## What is NOT in this directory

| | why |
|---|---|
| the **predicted** complexes | ~35 MB for the first round alone, across the complete **110-fold** cohort. The full 103-column scoring table is deposited instead; `../methods/INPUTS.md` §5 states the consequence |
| the repaired target as a standalone file | the target as used is inside all 11 designs above and is byte-identical to residues 6–157 of the challenge's own human sequence — asserted, and re-checkable |
| the other 373 generated backbones | the gate tables cover all 384 |
| the session (`.pse`) figures | the builders and their verifiers are deposited in `../scripts/` |
| any mouse-target structure | ⛔ mouse cross-reactivity is not addressed in this cohort |

---

## Counts stated here

11 design files · 1 binder chain + 3 target chains each · 4 guideposts each. ⚠️ Any count in this
file must be re-checked whenever a file is added or removed from this directory.
