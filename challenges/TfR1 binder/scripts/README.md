# The programs

Thirty-three files. Every number in `../methods/` and `../data/` came from one of them.

⚠️ **How these differ from the prose.** The documents in `../methods/` were written for a reader
from outside the group. These programs were written to be run, and their headers are working notes:
they argue with earlier versions of themselves, name defects they deliberately did not fix, and say
why a check exists. That is on purpose — a threshold with its reasoning stripped reads as an
arbitrary constant — and it is why the headers are worth reading alongside the methods documents.

⚠️ **Prepared for release 2026-11-01, extended 2026-11-02.** Comments and docstrings were edited:
machine-specific paths, internal cross-references and personal attributions are gone. **Executable
code was not touched**, with **nine declared exceptions**, each flagged in the file where it occurs:

- three sets of input-path defaults replaced with paths relative to this directory or with marked
  placeholders, and one printed column label (the 2026-11-01 pass);
- ⭐ five **message strings** in the phase-5 programs (the 2026-11-02 pass): in (25) a refusal
  message, a check label and one `print`; in (27) an argument's help text and a refusal message.
  ⚠️ Each carried an internal reference that had to go. **Only the text inside those messages
  changed** — no control flow, no threshold, no number, no comparison.

Every other program's parsed syntax tree is **identical** before and after, and for the two
exceptions above the trees differ **only** in those five string constants. That is what makes the
claim checkable rather than asserted.

⚠️ **What that verification does and does not prove.** It proves no executable statement changed
other than the nine declared. It does **not** prove the programs still run, because they were not
re-executed after editing.

⛔ **Two programs have usage lines naming inputs that are not distributed here** — the guidepost
constellation file, the receptor coordinates, and the unpacked generator output. They are flagged in
the file where they occur and in `../methods/README.md`, "Reproducibility".

---

## Phases

The work splits into five phases. Phases 1 and 2 produced negative results and their instruments are
retired. Phase 3 is the mining method that closed the line. Phase 4 is the scaffolded-design line
that followed from it, and phase 5 is the extension arm that tests whether what phase 4 built can
reach the receptor as well as clear it.

### Phase 1 — the Cα era (retired: no selectivity)

| # | file | what it does |
|---|---|---|
| 1 | `j847_rsyncPDB_mmcif_20261030.sh` | Mirrors the public mmCIF archive by rsync, from the distributor's own template script with the local paths and one server target filled in. ⚠️ **Carries a defect banner in its own header: the log file is written inside the mirrored directory, which the rsync deletes, so a run destroys its own log while writing it.** The mirror itself was unaffected and was verified against the filesystem instead. The banner is deliberately left in place and the defect deliberately unfixed — this file is the record of a run that already happened. |
| 2 | `j846_pdb_constellation_search_20261030.py` | Pass A: the **Cα-only** constellation search, residue identity ignored, over a 324,204-chain Cα store. `--build` parses the extracted coordinates into a compact array store; `--search` runs the search and emits one row per hit. Produced the negative in `../methods/mining_negative_results.md` §1, and the batched least-populated-shell search order that the later matcher inherits. |
| 3 | `j848_cacb_pilot_20261030.py` | Pilot: does adding the **beta carbon** restore selectivity, since Cα→Cβ fixes an orientation without needing full side chains? Run on structures already on disk. It did not. |
| 4 | `j850_density_map_20261030.py` | Average atomic density around the fixed positions in the Cα frame — shape before secondary structure. ⛔ **Its conclusions are withdrawn**: the superposition frame was built on a synthetic backbone and carried no site information. Superseded by (7) and (8). |

### Phase 2 — the functional-atom representation, and two measurements that returned nothing

| # | file | what it does |
|---|---|---|
| 5 | `j851_funcatom_search_20261030.py` | ⭐ **The core instrument, and the reference implementation.** Builds a functional-group store from structures (one typed centroid per arginine, lysine or histidine) and searches it for the designed constellation, with group-type compatibility enforced and the backbone left free. Everything from (7) onward imports this file and **does not edit it** — it stays the reference instance against which later copies are checked. |
| 6 | `j853_funcatom_density_20261030.py` | Re-runs the clash filter and the density map on the valid functional-atom frame, replacing (4). The clash filter was retired on the strength of this plus the two nulls (`../methods/mining_negative_results.md` §2). |
| 7 | `j855_masked_density_20261030.py` | The density map restricted to the **free volume** (target-occupied space masked out), with a **rotation null about the constellation centroid** as its control. Returned a null: ratio 0.934, i.e. matched scaffolds place marginally *less* mass in the usable volume than randomly oriented ones. |

### Phase 3 — the corrected curve, the named hits, and the three gates

| # | file | what it does |
|---|---|---|
| 8 | `j871_pincurve_admissible_20261031.py` | ⭐⭐ Re-derives the whole pin-count curve at an **admissible** pre-filter tolerance, `tol_d ≥ √(2k)·tol_rmsd`. Sweeps all 42 subsets of the six pins at two tolerances each (84 cells), **counts cap overflows separately from true zeros** — which the reference matcher cannot do — and carries a `--selfcheck` mode that proves its vendored copy of the matcher equivalent to (5) over 3,200 paired comparisons before any curve is believed. Also refuses a pre-existing output file, because the writers append. |
| 9 | `j871b_control_uniqueness_20261031.py` | Re-runs the "is such a constellation unique in nature?" positive control at an admissible tolerance, over several control chains rather than one. The original control's negative half does not survive. |
| 10 | `j891_emit_sixpin_hits_20261031.py` | Emits the **named** six-group hit list that the curve runs never wrote down: one row per chain at its best-RMSD placement, with the six matched residue numbers and their types. Its own gate: the chain count and best RMSD must equal the published six-pin row for that store, or the list is not from the instrument that produced the curve and must not be trusted. |
| 11 | `j892_hit_face_gate_20261031.py` | ⭐ The **accessibility gate**: per-residue relative solvent accessibility from an established SASA library, plus co-accessibility over 642 sphere directions, against a decision rule fixed before the run. Carries `--poscontrol`, a vacuity control that substitutes the chain's six most-exposed cationic and imidazole residues — because a gate that can only return FAIL is not a gate. Asserts the residue found at each recorded number against the recorded group type, and excludes any mismatch from every statistic. |
| 12 | `j893_candidate_scaffold_pse_20261031.py` | Superimposes each surviving candidate onto the constellation (Kabsch on six points, in the query's own pin order), measures body clash against the receptor, and builds a six-scene session file. ⭐ **Refuses to write output unless the RMSD it achieves equals the published value for that chain to 1e-3**, and asserts the residue name at every guidepost and every matched position before drawing anything. |
| 13 | `j893_verify_pse_20261031.py` | ⭐ Verifies the session file **by reloading it in a fresh viewer process**. A builder that exits successfully is not a verification: a non-ASCII character stored in a scene message crashes the file on load, and a named selection that matched nothing saves happily and shows nothing. This walks the whole saved session, asserts every stored string is ASCII, and asserts every object non-empty with its atom count. |

### Phase 4 — the scaffolded-design line: keep part of a real protein, generate the rest

These are the programs behind `../methods/core_scaffolded_design.md`. ⚠️ Programs (14) to (16) build
and check the design inputs; the generative run itself is not in this directory, and (17) onward read
its output.

| # | file | what it does |
|---|---|---|
| 14 | `j894_core_option_designer_20261031.py` | ⭐ Enumerates **core options**: which subsets of the scaffold's secondary-structure elements to hold as one rigid motif. Assigns elements, re-solves the placement per option (Kabsch at ≥3 contained guideposts; Kabsch plus a swept rotation at exactly 2, because two points do not fix a rigid body; rejected at ≤1), measures the resulting collision with the receptor, and sizes every linker as a range. ⛔ Writes no design input and runs no model. ⚠️ Its header records the pin-ordering defect that invalidated its own first run, and how an *impossible* number — a 4-point fit worse than the 6-point fit — exposed it. |
| 15 | `j895_build_core_arms_20261101.py` | Builds the generator input coordinates and the specification for each option. ⭐ Its header states the **three construct decisions** as choices rather than defaults, including the one whose consequence is that two of the six specified groups are not identity-constrained downstream. ⛔ Not a copy of the earlier single-residue-guidepost builder: that one is left untouched as the working reference, and its conventions are reused here and named at each point. |
| 16 | `j897_validate_core_arms_20261101.py` | ⭐⭐ Verifies the five properties of **this** construct that the generic spec validator cannot. Two of that validator's checks do not fail here, they **do not apply** — one treats every fixed residue as a 3-atom guidepost, which is false when whole motif residues are held, and the other's arithmetic test is skipped rather than passed. ⛔ **An inapplicable check is not a passing check**, so both underlying properties are measured directly: guideposts identified by **atom-name set** and never by "has 3 fixed atoms", and the origin token confirmed consumed as a pure translation over the guidepost atoms only. |
| 17 | `j900_screen_per_arm_20261101.py` | Aggregates the 400-design screen **per core option**, which the per-design gate output deliberately cannot do. ⛔ Applies **no** verdict column: it re-measures its own bar ladder, finds no admissible binary continuity cut on this cohort, and reports distributions plus explicit counts. ⛔ Reports `n` as **batches**, not designs, and prints the per-batch table so a 40/40 that is really 5/5 cannot read as strong evidence. |
| 18 | `j901_prep_pse_inputs_20261101.py` | Selects six designs across six distinct batches and **measures**, per design, which residues were given, which are fixed guideposts and which the generator built. ⭐ The held set cannot be read from the output file or the contig and the header says why; it is identified geometrically and the count is asserted. ⛔ Checks for generator pseudo-atoms by **element and name**, never by a text match, which returns zero and reads as a false clean. |
| 19 | `j901_build_screen_pse_20261101.py` | Builds the six-design viewer session in `../structures/`. |
| 20 | `j901_verify_pse_20261101.py` | ⭐ Verifies it **by reloading the saved file** in a fresh viewer process — 62 checks, 0 failures, including that the held and built sets are disjoint and together cover each chain exactly, and that the held core carries a *different colour value* rather than merely having had a colour command issued at it. The ASCII check walks the **saved bytes**, because stored scene text cannot be read back through this viewer's interface. |
| 21 | `j902_render_screen_pse_20261101.py` | Ray-traces the session to the PNGs in `../figures/`. ⚠️ Its header records the clipping defect that rendered a black wedge over the binding interface, and labels the interior-colour setting as a rendering choice rather than data. |
| 22 | `j902_stamp_renders_20261101.py` | Crops the uniform white border and stamps the title band. ⭐ Fits every stamped line to the available width and **refuses** below a minimum size, because a truncated column name is a wrong label and not a cosmetic problem. |
| 23 | `j903_collect_mpnn_20261101.py` | Collects the sequence-design output into one table. ⛔⛔ Its header records the counting trap it exists to avoid: concatenating the output FASTA files and counting headers **undercounts by exactly one per file**, because those files carry no trailing newline, so each file's last sequence line fuses onto the next file's header. Every file is parsed on its own and each arm's record count is asserted. ⛔ Also handles the first record of each file, which is the **input** sequence and not a design. |
| 24 | `j953_contained_pin_check_20261102.py` | ⭐⭐ Measures the identity of the **two contained** constellation positions in every deposited sequence, for all three cohorts — the two that no quality-control column covers, because they are not guideposts. Locates them from each design's **own** held-element sequence rather than by position, since every backbone samples its own connector lengths, and **refuses** if that element does not occur exactly once in a backbone's input row. |

⛔ **The sequence-design driver itself is not in this directory.** It is an existing program used
**unmodified** for this cohort, and it is shared with other targets in the group's work, so it is not
part of this challenge directory. What it was run with is stated in
`../methods/core_scaffolded_design.md` §7.

### Phase 5 — the extension arm: hold almost all of it, generate the reach

These are the programs behind `../methods/nterm_extension_arm.md`.

| # | file | what it does |
|---|---|---|
| 25 | `j945_build_nterm_extension_arm_20261102.py` | Builds the input coordinates and specification for the extension construct: drops three N-terminal residues, holds 4–145 sequence-fixed, re-places the guidepost arginine and writes the contig. ⛔ Writes no design and runs no model. ⭐ Asserts the finished construct rather than trusting the edit — total length, the held span, the guidepost's atom set, and that the held protein's sequence is **not** released to the sequence-design step. |
| 26 | `j948_nterm_design_score_20261102.py` | ⭐⭐ The arm's scorer, and the one program here whose own defect is documented in the method text. It gates continuity **before** any clash metric, resolves chain roles from each design's **own** specification, drops generator pseudo-atoms by element **and** name, and measures the generated segment's engagement separately from the whole chain's. ⛔ Its header records why the standard gate-chain driver **refused** this arm — that driver assigns chain roles by fixed-sequence fraction, and under partial diffusion every chain reads as fully fixed — and ⭐ the **corrected** derivation of the generated set: the binder chain's own residues minus the residues the specification holds, taken from the **structure**, never from the generator's token map, which names only held positions. It raises rather than continuing if that set comes back empty. |
| 27 | `j949_select_diverse20_20261102.py` | Picks the 20 sequences for the co-fold by **farthest-point diversity** rather than by model confidence. ⛔ Its header records why: a previous round measured that this model's confidence does not even separate the design arms, so ranking on it would be ranking on a quantity already shown not to carry the signal. ⛔ Also records the FASTA counting trap — those files carry no trailing newline — and parses each file on its own. |
| 28 | `j956_nterm_top3_pse_20261102.py` | Builds the corrected top-three viewer session in `../structures/`. ⭐⭐ **It refuses to write the session if the contact count it measures from the coordinates disagrees with the ranked table** — a figure must not claim a number it cannot reproduce. Also asserts the shared receptor frame across designs, drops pseudo-atoms by element and name, and asserts every stored scene message is pure ASCII. |
| 29 | `j956_verify_pse_20261102.py` | ⭐ Verifies that session by **reloading the saved file in a fresh viewer process** — 36 checks, 0 failures: 143 held and 27 built per design, the two sets disjoint, the receptor present at its full residue count, no surviving pseudo-atoms. The ASCII sweep walks the **saved bytes**, not strings still held in memory, because a non-ASCII mark in stored scene text is only dangerous once it has been written. |
| 30 | `j957_probe_pse_20261102.py` | Reads the session and reports what is in it before anything is rendered, so the render step is configured from measurement rather than from assumption. |
| 31 | `j957_render_top3_pse_20261102.py` | Ray-traces the panels. ⚠️ Carries an **additive** surface-transparency override, left at its default and raised only for the shipped figure; the defect that forced it is recorded in `../figures/README.md`. |
| 32 | `j957_stamp_panels_20261102.py` | Assembles and labels the panel figure, including the caption line that states the structures are raw and unrelaxed. |
| 33 | `j957_verify_labels_20261102.py` | ⭐ Re-opens each **finished** PNG and re-measures every stamped string with a real ink box rather than a font advance width, then reports overflowing labels and images with ink touching an edge. ⚠️ It passed on the first figure, which is the point of the note in `../figures/README.md`: the defect that figure had was one no label check can see. |

---

## Running them

All but the generative step are local-CPU only, and none of the programs in this directory needs a
GPU, model weights, or network access at run time — only (1) fetches anything, and (16) needs the
generator's own environment because it validates through the library's own code path.

Each program carries its own usage in its header or `--help`. The two that need saying out loud:

- ⛔ **Every writer appends to its output file.** Pointing a run at an existing results file silently
  merges two runs. Use a fresh output path per run; (8) enforces this by refusing a path that
  already exists.
- ⛔ **(12), (13), (19), (20) and (21) run under a headless viewer invocation and have no
  `__main__` guard, on purpose.** Under that invocation the module is not `__main__`, so a guarded
  script exits successfully having written nothing.
- ⛔ **Three programs need an input this directory does not distribute** — (2) and (12) the guidepost
  constellation file, (18) the unpacked generator output. Each carries a flagged placeholder where
  the path goes. ⛔ (1) carries two more, and the defect note above its path settings should be read
  before any re-mirror.
