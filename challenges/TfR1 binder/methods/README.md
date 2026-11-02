# Methods — reader order

Three documents, and the order matters.

| | |
|---|---|
| **1.** `structural_mining.md` | ⭐ **Start here.** The complete, self-contained account of the structural-mining method: what it searches for, how the constellation is defined and why the backbone is deliberately not part of it, the matching criterion and the tolerance algebra that criterion implies, the two databases and why their agreement is weaker evidence than it looks, the accessibility gate, and the superposition and clash test that closed the line. Someone outside the group should be able to reimplement the whole thing from this text. |
| **2.** `mining_negative_results.md` | The three measurements that were run with controls and returned nothing: the Cα-only search that has no selectivity, the target-clash filter that was retired for having no discriminating power, and the density map that returned a null against its own rotation control. Read after (1) — it assumes that document's definitions. |
| **3.** `core_scaffolded_design.md` | ⭐ **The work that follows from (1).** The mining line ended with four natural chains that present the designed arrangement on a reachable face and whose bodies all collide with the receptor, so none can be used whole. This document is the test of whether **part** of one can: keep subsets of a real protein's secondary-structure elements as a rigid core, generate the rest around them with the receptor present, screen ten such subsets, and design sequences on the three that survive. Read after (1); it assumes the constellation, the acceptance bar and the clash result. |

For what each program does and the order they run in, see `../scripts/README.md`. For what each
result table contains, column by column, see `../data/README.md`. For the renders and the viewer
sessions, `../figures/README.md` and `../structures/README.md`.

---

## ⛔ What is NOT documented here

**The first design round is not written up in this directory.** Its backbone generation, sequence
design and independent co-folding step are summarised in the challenge front page (`../README.md`,
"What was attempted") and no further — and it is that round's targeting failure that the mining line
exists to explain.

⚠️ Document (3) **partly** closes this gap, because it documents a second, later design round in
full: the construct, the ten-arm screen, the backbone gates and the sequence-design step with its
composition refit. ⛔ What is still missing is the **first** round — the one whose 2,160 predicted
complexes reproduced none of the designed contacts — and the receptor preparation both rounds share.
Neither is described anywhere in this directory.

---

## ⛔ Reproducibility: what can and cannot be re-run from here

**Can be re-run from what is here:** no **analysis** program in `../scripts/` depends on a GPU, on
model weights, or on anything proprietary. Every curve, count and ranking in all three method
documents can be recomputed from `../data/` alone — including the whole ten-arm screen aggregation
and every composition figure in `core_scaffolded_design.md` §7.

⚠️ **The generative step is the exception and it cannot be re-run from here.** The backbone
generation that produced the 400 screened designs needs the generator, its weights and an
accelerator; the sequence-design step needs its own model weights. What is deposited instead is the
**complete input specification of all ten arms** (`../data/core_arm_specs_20261101/`) and the
**measured output of both steps** (`../data/`), so the construct is auditable and the results are
checkable even though the runs are not reproducible from this directory alone.

**Cannot be re-run from what is here, and must be obtained upstream:**

- **The structure databases.** Two stores, 85,699 and 49,481 chains, built from a local mirror of
  Protein Data Bank mmCIF entries and from 15 whole prokaryote proteomes of the AlphaFold Protein
  Structure Database (v6). Neither is distributable from here: the mmCIF mirror is tens of gigabytes
  and the predicted-structure archives are the distributor's to serve. `structural_mining.md` §3
  states exactly what each store was built from, with the filters applied, so an equivalent store can
  be rebuilt. The third-party libraries the programs rely on — a structure parser, a solvent-accessible
  surface implementation, an array library and a molecular viewer — must be installed from their own
  upstream projects.
- **The design input.** The six-guidepost specification and its coordinate file are not deposited
  here. `structural_mining.md` §2 gives the constellation's definition — group types, constrained
  atoms, and the spans (32.55 Å at six points, 30.27 Å and 27.21 Å for the most compact five- and
  four-point subsets) — which is what the search actually consumes.
- **The receptor coordinates** used for the clash and free-volume measurements, and for every arm in
  `core_scaffolded_design.md`. ⛔ How the receptor was cropped — which entry it came from, and how its
  two-subunit assembly was handled — is **not described anywhere in this directory**, and the clash,
  free-volume and arm-collision numbers all depend on that crop.
- **The generated coordinates.** The 400 screened designs and the 312 designed sequences' structures
  are not deposited; the sequences themselves are (`../data/`), and six of the 400 backbones are
  viewable in `../structures/`.
- **Model weights** for the backbone generator and the sequence-design model, which are their own
  projects' to distribute.

⇒ ⛔ So the numbers here can be **read and audited as measurements of the databases**, but the
database-building step cannot be reproduced from this directory without obtaining the upstream
material.

---

## A note on how numbers are marked

**MEASURED** means a program in `../scripts/` was run and produced the value. **READ** means it was
taken from a published structure or paper. Where a value is an inference from measurements rather
than a measurement, it says so in the sentence that carries it.

⛔ Two conventions are applied throughout and are worth stating once:

1. **The acceptance bar travels with every count.** Every hit count in these documents was measured
   at a 2.0 Å RMSD acceptance bar, and both databases return **zero** six-group matches at a 1.0 Å
   bar. A count quoted without its bar is not a result.
2. **Supersessions and retractions are kept, not deleted.** Where an earlier conclusion was wrong
   about a method or a measurement, the corrected conclusion is carried forward *and* the fact that
   there was a correction is stated, because the correction is itself a finding. What is not kept is
   the prose history of how the wording changed.
