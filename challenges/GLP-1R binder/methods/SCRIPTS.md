# Programs — what each one does, what it expects, and what it refuses

Six programs are deposited in `../scripts/`. Together they produced every number in this directory.

⛔ **A citation to a program that is not deposited is a dangling reference, not provenance.** §3 names
the two programs that are *not* here and says exactly which deposited tables depend on them.

**Run order.** Generation and prediction happen upstream and are not deposited; everything below
reads files that were already on disk.

    [generation, upstream]  ->  [the two gate programs, NOT deposited -- see section 3]
    [prediction, upstream]
      ->  j559_baker_reference      (measures the two published complexes)  -> data/j559_...tsv
      ->  j637_baker_filter_c3      (the filter: 960 design folds + 27 controls) -> data/j637_..._v2.tsv
      ->  j639_msa_effect_...       (the matched single-factor control) -> data/j639_...tsv
                                     [imports the scoring function from j637]
      ->  j559_build_baker_pse      (viewer session: the published-versus-ours comparison)
      ->  j641_build_baker_best_pse (viewer session: the best folds, held residues marked)
            ->  j641_verify_pse     (reloads it and recalls every scene)

⛔ **These are lab programs, not a packaged tool.** They take absolute paths, several carry
**placeholder** paths marked in the file, and they assume the directory layout of the campaign that
produced them. They are deposited so the method can be *read and checked*, not so it can be re-run
unchanged.

All six were re-parsed as valid Python after the editorial pass that prepared this directory.
In-file paths are **placeholders** spelled `/path/to/your/working/tree`, with a comment saying so;
none of them is a working default.

---

## 1. The measurement programs

### `j559_baker_reference_20261019.py` — the geometry of the two published minibinders

**What it does.** Measures `min_cut` / `med_cut`, `rg_ratio`, the helix signature, the interface
residue sets on both sides, and the straight-line fit to the contacting binder Cα positions, on the
two deposited cryo-EM complexes named in `INPUTS.md` §2. ⭐ **The same code paths used on our own
designs**, so the two columns of `README.md` §3c are comparable by construction rather than by
assertion.

**Expects.** `--cif <file> [<file> ...] --out <fresh.tsv>`, plus `--binder-chain` / `--receptor-chain`
(defaults `A` and `R`). Needs `gemmi`, NumPy and SciPy.

**Refuses.**
- ⛔ if `--out` already exists — it **appends**, so pointing it at an existing file would silently
  merge two measurement runs;
- ⛔ if a named chain is absent from a structure, printing the chains it did find.

**⛔ Chains are named, never inferred.** The program reads each file's own
`_entity.pdbx_description` records and prints which entity it treated as the binder, so a reader can
see it picked the right one. A hardcoded chain pair has silently changed a conclusion in this group's
work before.

⛔ The downloaded structures are treated as untrusted data: they live in their own directory, the
program lives elsewhere, paths arrive as arguments, and it is run with the interpreter's isolated
mode.

**Writes.** `../data/j559_baker_reference_20261019.tsv` (2 rows, 17 columns).

---

### `j637_baker_filter_c3_20261024.py` — the published filter, on 960 design folds and 27 controls

**What it does.** Applies a published five-criterion in-silico filter, **using only its printed
thresholds**, to 960 folds of our designs and 27 known-binder control folds. Computes **both**
readings of the paper's ambiguous "min PAE" criterion — the literal minimum element over the
binder × target block, and the mean over it — and hides neither. Reports per fold (n=960) and per
design best-of-ten-seeds (n=96).

**Expects.** `--out <fresh.tsv>`. Reads prediction outputs that are **not deposited**
(`INPUTS.md` §5); it does not run any prediction.

**Refuses.**
- ⛔ if `--out` exists;
- ⛔⛔ **if the receptor sequence is not byte-identical across the design and control cohorts.** This
  is the one assertion the whole comparison rests on, and it is read out of the input files rather
  than eyeballed. If the receptors differed, a confidence gap would not be attributable to anything;
- ⛔ if a prediction directory does not hold exactly one error matrix, or if the inferred binder
  length is implausible, or if the design and control counts are not 960 and 96 distinct.

**⛔⛔ It carries its own bug report.** The function that reduces the interface error block documents,
in its docstring, the chain-order transposition that an earlier version of this program had, what it
produced, and why a plausibility range check could not catch it. `README.md` §4c is the long form.
⇒ the block offsets are now passed in **explicitly per cohort** and never inferred from position.

**⭐ Pairwise ipTM, never the global number**, with the reason in a comment at the point of use.

**Writes.** `../data/j637_baker_filter_c3_20261024_v2.tsv` (987 rows: 960 + 27).

---

### `j639_msa_effect_on_baker_channels_20261024.py` — the matched single-factor control

**What it does.** Scores two arms of an earlier run of the **same** ligand against the **same**
receptor that differ in **exactly one factor**, the receptor sequence alignment, 21 draws each. The
difference between them *is* the alignment's effect on each filter channel, measured on the native
ligand with everything else held. This is what resolves the confound stated in `README.md` §5a.

**Expects.** `--out <fresh.tsv>`. Imports the scoring function from
`j637_baker_filter_c3_20261024.py`, which must sit beside it. Reads prediction outputs that are not
deposited.

**Refuses.**
- ⛔ if `--out` exists;
- ⛔ if an arm's prediction directory is empty, naming the pattern it looked for;
- ⛔ if a scored record's binder length is not the control peptide's length, which is the check that
  the chain layout was read correctly.

**⚠️ It states its own limits in its docstring, and they matter**: both arms carry a forced
structural template while our 960 folds carry none, so the arm is a matched control for the
*alignment question only* and its absolute values must not be read against our designs as though the
template were absent.

**Writes.** `../data/j639_msa_effect_20261024.tsv` (42 rows).

---

## 2. The viewer-session programs

### `j559_build_baker_pse_20261019.py` — the strand-A comparison as a navigable session

**What it does.** Builds a PyMOL session: the two published complexes, their measured contact
residues, the published binder alone, and one of our own designs in the same scene for contrast.
Receptor-contact residue lists come from the measurement program's output, not from prose.

**Expects.** `pymol -cq <this file>`, with the paths at the top set for your own tree. Needs the two
cryo-EM structures from `INPUTS.md` §2 and the model in `../structures/`.

**Refuses.** ⛔ **Every selection is asserted non-empty** and the program raises rather than saving a
scene that shows nothing.

⛔ **Two PyMOL hazards are handled in the file and both are real.** There is deliberately **no
`__main__` guard**, because `pymol -cq` on a guarded script exits 0 having written nothing; and every
stored scene message is **pure ASCII**, because a non-ASCII character in one corrupts the heap and
the saved session then crashes when it is opened.

⛔⛔ **One scene is deliberately not stored, and the reason is a trap worth reading.** A scene that
painted a named receptor site on our own design's target chain was removed rather than fixed: that
chain is renumbered by the generator with its segment breaks collapsed, so a selection by site
residue number **returns atoms** — which is exactly why an is-the-selection-empty guard cannot catch
it — but they are the **wrong residues**. The comment block in the file gives the measured
chain-to-construct mapping. See `../structures/STRUCTURES.md`.

### `j641_build_baker_best_pse_20261024.py` — the best-scoring folds, with the held residues marked

**What it does.** Builds a session of the best and worst-contrast design folds by the filter's two
surviving channels, with the seven held ("pinned") binder residues drawn and labelled by which atoms
of each were held, plus a native-ligand reference object so a reader can see where a real agonist
sits against the same receptor. All objects are superposed on the **receptor**, never on the binder,
because the binder pose is the thing being looked at.

**Expects.** `pymol -cq <this file>`. Reads prediction outputs that are not deposited.

**Refuses.** ⛔ if a fold file is missing; ⛔ if a receptor chain does not hold 409 Cα; ⛔ if a binder
chain falls outside the length range its generator specification predicts; ⛔ if any object does not
carry all seven held residues; ⛔ if a scene message is not pure ASCII.

⭐ **The held positions are derived, not remembered.** The program walks the generator contig
specification for that batch and derives the positions from it, then checks the derivation against an
independent prediction — that the specification implies binder lengths 55-65, which is what the fold
cohort measures. ⚠️ That specification file is **not deposited**; the contig string itself is
reproduced in the docstring, which is what the derivation actually rests on.

### `j641_verify_pse_20261024.py` — the reload check

**What it does.** Reloads a saved session, recalls **every** scene, and asserts that every design
object still carries all seven held residues.

**Why it exists.** ⛔ Writing a session file proves nothing; only a reload does. The failure it guards
against — a non-ASCII character in a stored message corrupting the heap so the file crashes on open —
is not hypothetical.

**Refuses.** ⛔ if any object is missing a held residue. Like its sibling it has **no `__main__`
guard**, for the same reason.

---

## 3. ⚠️ The two programs that are NOT deposited, and what depends on them

Named here rather than left for a reader to discover.

| program, as the tables' own headers name it | what it wrote | which deposited tables |
|---|---|---|
| `glp1r_fold_novelty_20261016.py` | folded each designed sequence **alone** and scored it against its own designed backbone | `glp1r_fold_novelty_sel20_20261016.tsv`, `j558_fold_unbiased_20261019.tsv`, `j558_fold_biased_20261019.tsv` — the three tables the **0 of 352** fold-back result pools |
| `j557_topology_gate_20261019.py` | the topology and compactness gate on the 384 first-round backbones | `j557_topology_deep6_20261019.tsv`, `j557_topology_shallow4_20261019.tsv` — the "our 384" column of `README.md` §3c |

⇒ **What this costs a reader.** Both results can be *read* out of the deposited tables, and the
pooling arithmetic checks (160 + 96 + 96 = 352; 288 + 96 = 384). What cannot be checked from this
directory is **how** either number was computed — the fold condition, the scoring definition and the
gate's own thresholds are described in `README.md` §3b but the code that implements them is not here.
⛔ By the convention stated at the top of this file, those two header citations are dangling
references.
