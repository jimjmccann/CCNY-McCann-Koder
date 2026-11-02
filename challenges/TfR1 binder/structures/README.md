# Structures

Two molecular-viewer session files, one per line of work. ⛔ Neither holds coordinates of a
designed binder in a deposited, reusable form — see the closing note.

| file | what it holds |
|---|---|
| `j893_candidate_scaffolds_on_constellation.pse` | the four natural candidate scaffolds the mining line ended on, superimposed on the designed constellation |
| `j901_screen_best_designs_20261101.pse` | six generated designs from the scaffolded-design screen, coloured by what was given and what was built |

---

## `j893_candidate_scaffolds_on_constellation.pse`

A molecular-viewer session (PyMOL 3.1 format, 2.7 MB) holding the **four distinct small candidate
scaffolds** that survived the accessibility gate, each superimposed on the six designed functional
groups in the receptor's own coordinate frame, with the receptor crop present.

⛔⛔ **This is a VIEW, not evidence.** Nothing in it ranks, scores or re-measures anything, and ⛔ **a
maintained geometry is never evidence of affinity.** It exists because the thing it shows is easier
to see than to read: all four candidates occupy the receptor's volume.

**Six scenes:**

| scene | what it shows |
|---|---|
| 1 | all four candidates on the constellation at once |
| 2 | the six designed groups alone, against the receptor crop — four arginines and two histidines |
| 3–6 | one candidate each, carrying that candidate's residue count, its constellation RMSD against the 2.0 Å bar, its clashing-atom count and its closest approach to the receptor |

**What is drawn, and how it was checked:**

- The superposition is Kabsch on **six points only** — each candidate's matched functional-group
  centroid onto the corresponding designed centroid, in the query's own fixed pin order.
- ⭐ The builder **refuses to save** unless the RMSD it achieves equals the value published for that
  chain to 1e-3. All four passed to the digit: `5xox_B` 1.5400, `9gyt_N` 1.8251, `9t6h_b` 1.9375,
  `AF-P56144` 1.9216 Å. A superposition that does not reproduce the published number is a different
  alignment, and a picture built on it would show a measurement nobody made.
- ⭐ The residue **name** at every designed position and at every matched position is asserted
  against the recorded group type before anything is drawn, so a silent off-by-one cannot superimpose
  on the wrong atoms and look perfect.
- ⭐ The file was **verified by reloading it in a fresh viewer process** before being handed on: 17
  checks, 0 failures, including an assertion that every string stored anywhere in the saved session
  is ASCII. A non-ASCII character in a stored scene message corrupts the heap and crashes the file on
  load, and the internal layout of scene data is version-specific, so walking the whole session is
  more durable than checking scene messages alone.

**The measurement the file surfaced**, close **atom pairs** between the placed candidate and the
receptor crop within 2.5 Å, with closest approach. ⛔ **Pairs, not distinct atoms** — see
`../methods/structural_mining.md` §8.2.

> ⚠️ **The scene captions stored inside `j893_candidate_scaffolds_on_constellation.pse` still
> call these numbers "heavy atoms".** That wording is wrong in the same way this table's was: the
> builder sums every entry of the candidate-by-receptor distance matrix below 2.5 Å, so the figure
> counts close atom **pairs**. The builder script was corrected, but the session file is deposited
> as it was generated rather than rebuilt, so the numbers in it are reproducible from the deposited
> scripts. ⛔ Read the captions' counts as atom pairs. The counts themselves are correct.

| candidate | length | RMSD (Å) | clashing atom pairs | closest approach (Å) |
|---|---:|---:|---:|---:|
| `5xox_B` | 229 | 1.5400 | 679 | 0.40 |
| `9gyt_N` | 116 | 1.8251 | 1,271 | 0.15 |
| `9t6h_b` | 222 | 1.9375 | 2,469 | 0.15 |
| `AF-P56144` | 242 | 1.9216 | 1,694 | 0.07 |

⚠️ **Not an artefact of which placement was drawn:** three of the four chains have exactly one
accepted placement, and only `9t6h_b` has two — and it is the worst-fitting candidate anyway.

---

## `j901_screen_best_designs_20261101.pse`

A molecular-viewer session (PyMOL 3.1 format, 4.4 MB) holding **six generated designs** from the
ten-arm scaffolded-design screen — three from the leading arm, two from the runner-up and one from
the third arm — in the receptor's own coordinate frame, with the receptor surface present. The method
is `../methods/core_scaffolded_design.md`; renders of this session are in `../figures/`.

⛔⛔ **This is a VIEW, not evidence.** Nothing in it is scored for binding: no co-fold and no
sequence design is represented here, and the clash numbers in its scene text are **raw** generator
output, which clashes normally.

**Six designs from six distinct batches**, which is deliberate: one contig length is drawn per batch,
so six designs from one batch would be six views of one draw.

**What the colouring shows** — three disjoint residue sets per design:

| colour | meaning |
|---|---|
| purple | the kept 5XOX core elements, i.e. what the generator was **given** |
| grey | what the generator **built** to connect them |
| yellow | the four fixed guidepost residues |
| teal | the receptor |

⭐ **The split is the result:** the generator built **85 to 98 residues** of each **141–145** residue
binder.

**How it was checked, and why each check is there:**

- ⭐ The held set is **measured geometrically**, not read from the file or the contig, and it could
  not have been: the output file flags only the four guideposts as sequence-fixed, and linker lengths
  are ranges drawn per batch so a held element's output index moves design to design
  (`../methods/core_scaffolded_design.md` §6). Every design matched its full expected residue count —
  43/43, 52/52, 45/45 — with a worst deviation of **0.0623 Å**.
- ⭐ **A common frame is asserted, not assumed.** Every design carries the same receptor crop, so the
  target-on-target RMSD to the first object must be about zero. It is **printed at 0.0000 Å for all
  five comparisons**, which is what makes the six directly comparable.
- ⭐ The file was **verified by reloading it in a fresh viewer process**: **62 checks, 0 failures** —
  6 objects, 7 scenes, per-object atom counts against the manifest, target present, zero
  pseudo-atoms, held and pin sets both selectable at full count, the two sets **disjoint**, held plus
  pins plus built covering each chain exactly, and the held core confirmed to carry a *different
  colour value* from the built region rather than merely having had a colour command issued at it.
- ⭐ Every stored scene string is **asserted ASCII before it is stored**, and the verifier re-checks
  this by recovering the stored text from the **saved bytes**, because this viewer version offers no
  way to read stored scene messages back through its interface. A non-ASCII character in a stored
  scene message corrupts the heap and the saved file crashes on load.
- ⭐ **Pseudo-atoms are excluded by element *and* name**, never by a text match on the file — a text
  match returns zero and reads as a false clean. The input cohort measured zero, and the verifier
  asserts zero per object anyway rather than depending on that staying true.

⚠️ **String-audited for release.** Unlike the file above, this one's full string table was dumped and
read: it contains object names, scene text (measurement sentences and the colour legend) and nothing
else — no filesystem path, no machine name, no working-directory string. ⚠️ The same audit has **not**
been done on the older session file above; its stored strings were verified to be ASCII and its scene
text was read, but its binary has not been searched exhaustively.

---

⛔ **No coordinates of our own designs are deposited in this directory in a reusable form.** Both
files are viewer sessions. No predicted binder–receptor complex from either design round is
deposited, and no designed sequence here has been co-folded at all.
