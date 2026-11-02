# The result tables

Two groups of files.

**§1–§3 — the structural-mining line.** Six tab-separated files, all with a header row. Every count
and every ranking in `../methods/structural_mining.md` and `../methods/mining_negative_results.md`
can be recomputed from these.

**§4–§6 — the scaffolded-design line** (`../methods/core_scaffolded_design.md`): the ten arms'
complete input specifications, the 400-design screen, and the designed sequences with their quality
control.

---

Two structure stores are reported separately throughout §1–§3 and are never pooled silently. In the
filenames:

- **`_fx_`** — the **PDB-derived store**, 85,699 chains built from 25,821 Protein Data Bank mmCIF
  entries held locally.
- **`_afdb_`** — the **predicted-model store**, 49,481 chains from 15 whole prokaryote proteomes of
  the AlphaFold Protein Structure Database (v6), filtered at per-residue pLDDT ≥ 70.

⛔⛔ **Every RMSD in every one of these files was measured against a 2.0 Å acceptance bar.** At a
1.0 Å bar both stores return **zero** six-group matches. A count taken from these files without its
bar is not a result.

⛔ **The stores count chain ENTRIES, not unique chain names.** The 795-row PDB-derived hit list holds
786 distinct names, because nine chains entered from two overlapping source subsets. See
`../methods/structural_mining.md` §3.2.

---

## 1. `j871_pincurve_fx_20261031.tsv` · `j871_pincurve_afdb_20261031.tsv`

**84 data rows each.** The pin-count curve: how many chains present *k* of the six designed
functional groups, for every subset of the six, at two pre-filter tolerances.

| column | meaning |
|---|---|
| `k` | how many of the six groups this row's subset asks for (3, 4, 5 or 6) |
| `subset` | which groups, by their identifiers in the design specification, comma-separated. `B1`, `D1`, `E1`, `G1` are cation groups; `C1` and `F1` are imidazoles |
| `tol_d` | the pairwise-distance **pre-filter** tolerance in Å used for this row |
| `is_floor` | `0` = the historical, **inadmissible** `tol_d` of 1.5 Å, kept so the published numbers can be reproduced. `1` = the admissible floor `√(2k)·tol_rmsd` for this `k`, i.e. 4.899 / 5.657 / 6.325 / 6.928 Å at `k` = 3/4/5/6. ⛔ **Only `is_floor = 1` rows are valid statements about the databases**; see `../methods/structural_mining.md` §5 |
| `hits` | accepted placements. One chain can contribute several, so this is always ≥ `chains` |
| `chains` | distinct store entries with at least one accepted placement |
| `pct_store` | `chains` as a percentage of the store (denominator 85,699 or 49,481) |
| `best_rmsd` | the smallest RMSD achieved in this cell, in Å. **Empty when `hits` is 0** |
| `cap_overflow_chains` | chains abandoned because an intermediate candidate set exceeded the implementation's cap. ⛔ **This column exists because the underlying matcher returns the same empty result for an overflow and for a true zero.** It is **0 in all 84 rows of both files**, so the defect is latent here and corrupted nothing — but a curve reported without this column is not auditable |

⭐ **The two rows that carry the central correction**, both files, `k = 6`:

| store | `tol_d` 1.5 Å (inadmissible) | `tol_d` 6.9282 Å (floor) |
|---|---|---|
| PDB-derived | 0 hits / 0 chains / 0.0000% / no RMSD | 1,004 hits / **795 chains** / **0.9277%** / **1.0680 Å** |
| predicted models | 0 hits / 0 chains / 0.0000% / no RMSD | 657 hits / **519 chains** / **1.0489%** / **1.1079 Å** |

⇒ the zero was the pre-filter, not the biology.

⚠️ **Monotonicity is the cheap self-check on this file:** for every subset, the `is_floor = 1` row
must have `hits` and `chains` greater than or equal to its `is_floor = 0` partner, because widening a
pre-filter can only admit candidates it previously discarded. 0 violations in 42 of 42 subsets in
each store.

⚠️ Neither file records how many worker processes the run used, so the two runs' wall times are **not
comparable**. The counts are unaffected — the store is sharded by chain and the shard results summed.

---

## 2. `j891_sixpin_hits_fx_20261031.tsv` · `j891_sixpin_hits_afdb_20261031.tsv`

**795 and 519 data rows.** The named six-group hit list — one row per chain, at that chain's
best-RMSD placement, ranked ascending by RMSD. This is the file that turns a count into something
that can be pulled out of a real structure and looked at.

| column | meaning |
|---|---|
| `rank` | 1 = best fit in this store |
| `chain` | the chain, as the store names it. PDB-derived entries are `<entry>_<chain>` with the **author** chain name, which is case-sensitive (`9t6h_b` is not `9t6h_B`). Predicted-model entries are the distributor's own accession stem plus `_<chain>`; ⛔ split such a name on its **last** underscore, because the stem itself contains underscores |
| `best_rmsd` | Å, from the superposition of the six matched group centroids onto the six designed ones |
| `n_placements` | how many distinct accepted placements this chain had. Only the best one is described by the residue columns |
| `resid_1` … `resid_6` | the residue sequence numbers of the six matched residues in that chain, **in the designed pins' own fixed order** — the same order as `type_1..6`, and the same correspondence the superposition uses |
| `type_1` … `type_6` | `CATION` (arginine guanidinium or lysine `NZ`) or `IMID` (histidine imidazole). ⛔ Constant down each file, because it is a property of the query, not of the hit: `CATION IMID CATION CATION IMID CATION` |

⭐ **These files are the evidence for the size-versus-fit opposition.** Only 2 of the 795
PDB-derived chains are under 120 residues; 91% are over 250; the shortest chain presenting the
constellation in either database is 116 residues; and the best hit at ≤ 250 residues is
1.4662 Å against the list's overall best of 1.0680 Å.

⛔ **Length here is a true residue count read from each structure**, taken from the `chain_len`
column of the face-gate files. ⛔ It is **not** the highest residue number in a chain — that
measure is not a length and not a bound on one, and an earlier version of these notes used it.

⛔ **A row is a geometric coincidence, not a scaffold.** The next two files are what test that.

---

## 3. `j892_face_gate_fx_20261031.tsv` · `j892_face_gate_afdb_20261031.tsv`

**795 and 519 data rows** — one per hit from the file above, same chain names, same carried RMSD.
The accessibility gate: are the six matched groups on one reachable face, or buried?

| column | meaning |
|---|---|
| `chain` | as above, so this file joins to the hit list on `chain` |
| `best_rmsd` | carried through unchanged from the hit list, for convenience |
| `chain_len` | residue count of the chain as read from the structure |
| `verdict` | `CANDIDATE`, `PARTIAL` or `FAIL`, from a rule fixed **before** the run: `CANDIDATE` = `A_exposed_of6` ≥ 5 **and** `B_coaccess_of6` = 6; `PARTIAL` = `A` ≥ 4 **and** `B` ≥ 5; `FAIL` otherwise |
| `A_exposed_of6` | how many of the six matched residues reach relative solvent accessibility ≥ 0.20, computed on the chain **in isolation** with `freesasa` |
| `B_coaccess_of6` | co-accessibility: the maximum, over 642 sphere directions, of how many of the six groups are visible along that one direction. A group is visible if no chain heavy atom lies within 3.0 Å of its outward ray between 1 and 15 Å. The group's **own** residue is excluded from the occlusion test; the other five groups' atoms are **kept**, because a partner really is blocked by them |
| `min_rsa` | the lowest relative solvent accessibility among the six matched residues — the limiting one |
| `max_outward_angle_deg` | largest pairwise angle between the six groups' outward vectors. Descriptive only; no threshold is applied to it |

**Totals across both files:** 31 `CANDIDATE` (26 + 5), 143 `PARTIAL` (93 + 50), 1,140 `FAIL`
(676 + 464) of 1,314. ⚠️ These are **row** counts; per the entry note above, the 26 PDB-derived
`CANDIDATE` rows are **25 distinct chains** (`9s2e_C` entered from two overlapping subsets), so the
headline set is **30 distinct chains**, which collapse to 19 PDB entries plus 5 predicted models.

⛔⛔ **The ordering of the `verdict` column against `best_rmsd` is the finding:** the best RMSD among
`FAIL` rows is **1.0680 Å** and the best among `CANDIDATE` rows is **1.5400 Å**. Sorting either file
by `best_rmsd` and reading the `verdict` column shows it directly — the whole tight head of the list
is buried.

⚠️ **What a `CANDIDATE` verdict does and does not mean.** It means a partner could reach all six
groups from one direction. ⛔ It does **not** mean the host protein's body stays clear of the target:
all four of the small surviving candidates interpenetrate the receptor by hundreds to thousands of
close **atom pairs** when actually placed (⛔ pairs, not distinct atoms). See `../methods/structural_mining.md` §8.

⚠️ There are **no lookup-failure rows** in either file. The gate asserts that the residue found at
each recorded number matches the recorded group type, and excludes any mismatch from every statistic;
it fired **0 times in 1,314**. A run of this kind with 5% or more failures should not be published at
all, because the identity mapping between the store and the real structures would be broken.

---
---

# The scaffolded-design line

Method: `../methods/core_scaffolded_design.md`. ⛔ **Nothing in §4–§6 is a binding result.** No
sequence in §6 has been folded against the receptor or scored.

## 4. `core_arm_specs_20261101/` — the ten arms, as actually built and screened

Twenty small JSON files, two per arm: `tfr1_core_opt<NN>_spec.json` and
`tfr1_core_opt<NN>_provenance.json`, for `NN` = 01 … 10.

⚠️⚠️ **These are the authoritative record of what was screened**, and they matter because an earlier
option set was superseded before anything ran. The first selection sized its linkers on an
extended-chain bound; re-sizing them for an all-helical chain pushed every total outside the
specified 120–160 residue band, so a second selection holding fewer and smaller elements was made and
**that** is what ran. ⛔ The earlier set is **not** deposited here. ⇒ if a figure elsewhere disagrees
with these files about which elements an arm holds, **these files govern**.

### The specification file

| key | meaning |
|---|---|
| `input` | the coordinate file the generator was given. ⛔ **Not deposited here** — it contains the receptor crop |
| `contig` | the generated layout, in sequence order: held element segments (`H<lo>-<hi>`), linker **ranges** (`<lo>-<hi>`), the four guidepost labels, then the target block. ⛔ Each chain-break token carries **its own comma**; written without one it registers no break and silently changes the output layout |
| `select_fixed_atoms` | the **four** guideposts whose functional tip atoms are constrained, with the atom names: `NE,NH1,NH2` for a cation group, `ND1,NE2,CG` for an imidazole |
| `select_unfixed_sequence` | the held element segments whose **sequence is released**. ⛔ Read this against `select_fixed_atoms`: the two constellation residues that fall inside a held element are inside a sequence-released segment, so only **four** of the six groups are identity-constrained downstream (`../methods/core_scaffolded_design.md` §4) |
| `ori_token` | the guidepost centroid, used as a pure translation. Measured consumed at exactly 0.0000 Å |

### The provenance sidecar

| key | meaning |
|---|---|
| `option`, `control` | the arm number, and whether it is the deliberate positive control |
| `elements` | which scaffold secondary-structure elements this arm holds. Ranges are in `../methods/core_scaffolded_design.md` §2 |
| `contained_pins`, `bare_guideposts` | constellation residues delivered **by the held backbone**, and those delivered as single-residue guideposts |
| `scaffold_identity_at_contained_pins` / `design_identity_at_those_positions` | ⚠️ where the two differ, the construct took the **scaffold's** identity. The search that selected this scaffold treated an arginine guanidinium and a lysine `NZ` as the same cation type, so the match allowed it — but a guanidinium is planar and bidentate and a lysine `NZ` is a point charge, and that is a real chemical change |
| `placement_mode` | `kabsch` where ≥3 guideposts are contained; `sweep<angle>` where exactly 2 are, naming the swept rotation kept |
| `pin_rmsd` | ⛔⛔ the fit over the **contained** guideposts only — **2 or 3 points, not 6.** It is **not** comparable to the 1.5400 Å six-group figure and must never be quoted as "the scaffold fits to 0.9 Å" |
| `clash_atoms_vs_target` | close **atom pairs** between the placed motif and the receptor crop within 2.5 Å, **before** anything was generated. ⛔ **Pairs, not distinct atoms** — the builder sums every entry of the motif-by-receptor distance matrix below the cut, so an atom with several close partners counts once per partner. The field name is kept as written because the deposited scripts read it. ⭐ **0 for nine of the ten arms; 379 for `opt10`** |
| `len_min`, `len_max` | the realised binder length band, from the linker ranges |
| `gaps` | one entry per linker: `d` the Cα–Cα distance to span, `lo`/`hi` the range requested, `floor` the extended-chain minimum. ⛔ Each helical `lo` is asserted at or **above** its extended floor, so the bound is never loosened |
| `rise_helix`, `rise_ext`, `slack`, `n_segments` | the sizing constants (1.5 and 3.3 Å per residue), the range width, and how many held segments the arm has |

⭐ **`n_segments` is the column to read against the screen result:** the arms that fail do so with
many small scattered segments, and the two distinct failure modes are in
`../methods/core_scaffolded_design.md` §5.3.

## 5. `j900_screen_gate_20261101.tsv` and `j900_per_arm_20261101.tsv` — the ten-arm screen

**400 data rows and 10 data rows.** 40 designs per arm, ten arms, with the receptor present.

⛔⛔ **No verdict column exists in either file, and that is deliberate.** The gate re-measured its own
bar ladder on this cohort and **refused** a binary continuity cut, because the counts at 5, 7 and
10 Å differ (351 / 353 / 371) — so designs genuinely sit between "ideal" and "broken". Every bar used
in the method document is a **reporting** bar. Anyone applying one is making a choice, not reading a
verdict.

### `j900_screen_gate_20261101.tsv` — one row per design, 35 columns

Its first three lines are `#` comments and carry the aggregation rule, the worse-direction of every
column, and the statement that every column is continuous. Read them; they are part of the table.

| column | meaning |
|---|---|
| `design` | the design identifier. The **arm** is the `opt<NN>` inside it |
| `chains_binder`, `chains_target` | which chains are which, **detected** and then recorded, never named by hand |
| `n_copies`, `n_copies_used` | copies of the design compared, and how many entered the aggregation |
| `max_step`, `median_step` | Cα–Cα step lengths. **`max_step` is the continuity column** — an ideal chain sits near 3.9 Å |
| `min_cut`, `med_cut` | minimum graph cut of the long-range contact map: the "is this chain **one** fold" number. A low value means no crossing long-range contacts, i.e. two lobes on a tether |
| `lr_frac` | fraction of contacts that are long-range |
| `rg`, `rg_ratio` | radius of gyration, and its ratio to the expectation for that length. ⛔ **`rg_ratio` is the only length-safe compactness column**, which matters because the arms differ in length |
| `bb_min_to_target`, `heavy_min_to_target`, `n_res_under_2p5` | clearance from the receptor. ⛔ **Raw generator output clashes normally**; `bb_min ≥ 2.0` is the documented **raw** floor, not a post-relaxation gate, and it is never used alone |
| `pin_cb_min_to_target`, `pin_cb_max_to_target` | guidepost reach, **Cβ-based**. ⛔ Printed and **not** ranked on: clearance must be gated on pin **heavy** atoms, and a Cβ-based gate has passed a 1.67 Å clash elsewhere in this project |
| `n_res`, `n_pins` | realised binder length, and fixed positions found. **4 in 360 rows, 3 in 40** — the 40 are the control arm, which has one fewer bare guidepost because one of its guideposts falls inside the element it holds |
| `*_spread` | max − min across copies of the same design. ⭐ **Nonzero means the copies disagree**, so a single-copy reading of that column would have been wrong |
| `n_numbering_gaps` | breaks in residue numbering |

### `j900_per_arm_20261101.tsv` — one row per arm

Medians of the continuous columns plus explicit counts: `n_bb_ge_2p0`, and
`n_step_le_4.0 … n_step_le_10.0`, which is the bar ladder per arm.

⛔⛔ **`n_batches` is the effective sample size, not `n`.** One contig length is drawn per batch and
the eight models in a batch share that draw, so each arm is **5** independent observations and not 40.
The file's own header comment says so. ⇒ ⛔ **a 40/40 is 5/5.** The per-batch consistency that the
method document reports — three arms at 8/8 in every batch, one running 4/8, 2/8, 0/8, 3/8, 6/8 — is
what makes the continuity signal more than one lucky draw.

⚠️ **The conjunction counts are not a column.** They are computed by intersecting criteria over the
per-design file, deliberately, because ⛔ **a mean of fractions pays partial credit** for failing one
test outright. The intersections quoted in the method document (16 / 12 / 11 / 6 / 5 / 3 / 3 / 1 / 0 / 0)
were re-derived from the 400-row file during preparation of this directory and agree exactly.

## 6. The designed sequences

### `j903_cohort_opt05.txt`, `…opt01.txt`, `…opt04.txt`

The **39** backbones sequence design was run on — 16, 12 and 11 — one identifier per line. These are
the **three-way conjunction survivors**: designs that passed continuity **and** the raw clash floor
**and** the single-fold criterion, never any one alone.

### `j904_mpnn_sequences_biased_20261101.tsv` — ⛔⛔ superseded for the designed site; kept as the record

> ⛔⛔ **SUPERSEDED by `j908_mpnn_sequences_sixpin_20261101.tsv` for anything that depends on the
> designed site.** This cohort holds four of the six specified functional groups. The two that sit
> inside held element 3 were released to sequence design, and position 34 — where an imidazole was
> specified — came back as an imidazole in **0 of 312** sequences. ⚠⚠ Its bulk composition,
> confidence, charge and cysteine statistics remain valid and are still the ones quoted; it is the
> **site chemistry** that does not hold. Full measurement: `../methods/core_scaffolded_design.md`
> §7.5 and `j953_contained_pin_check_20261102.tsv`.
>
> ⭐ **It is deliberately kept**, because the interaction that broke it — an omission list meeting a
> released held element — is the finding, and a cohort deleted is a lesson deleted.

**351 rows: 312 designed sequences plus 39 inputs.**

| column | meaning |
|---|---|
| `arm`, `design`, `seq_id` | which arm, which backbone, which of the 8 sequences for it |
| `is_input` | ⛔ **`1` means this row is the backbone's own input sequence, not a design.** It is kept and labelled, and excluded from every statistic — counting it would inflate each backbone's yield by one ninth |
| `n_chains`, `len_binder` | chains in the record, and the designed chain's length (141–148) |
| `overall_confidence`, `seq_rec` | the model's confidence, and recovery against the input sequence. ⛔ **Confidence does not separate the arms** — the three medians sit inside 0.010 of each other — and it is **not** evidence of binding |
| `n_cys` | ⭐ **0 in all 312 designed chains.** ⚠️ The *target* chains carry their own two cysteines, which is why a naive search across a whole multi-chain record reads as a failure when nothing is wrong |
| `n_his` | histidines in the designed chain. **1 in all 312.** ⛔⛔ **Read this as a defect, not a clean sheet:** the specification names **two** imidazoles, so an intact design must carry `n_his >= 2`. The single histidine is the *fixed* guidepost; the *contained* one at scaffold 34 is absent in every sequence. An earlier reading of this column treated the missing constellation histidine as a suppressed stray |
| `sequence` | the designed chain only; target chains dropped |

### `j903_mpnn_sequences_unbiased_20261101.tsv` — ⛔ superseded, kept on purpose

> ⛔⛔ **SUPERSEDED for any downstream use. Do not rank or fold from this file.** It was produced
> without the composition bias and without histidine omission. The measured cost: stray histidines in
> **200 of 312** sequences (64.1%), since reduced to **0 of 312**, and an excess of glycine and
> alanine since reduced.
>
> ⭐ **It is deliberately kept.** It is the unbiased **round 0** the bias refit was fitted against and
> the baseline the paired confidence comparison is computed against. ⭐⭐ Everything in §6's guard list
> — 0 cysteines, all fixed positions held, 312 distinct sequences, lengths in band — was measured on
> **this** file first and still holds on the adopted one, as does the finding that confidence does not
> separate the arms.

Same columns as the adopted file. ⚠️ Its `n_his` ranges over `{1,2,3,4,5,6}`; the adopted file's is
uniformly 1.

### `j904_seq_qc_biased_20261101.tsv` — quality control, 312 rows

| column | meaning |
|---|---|
| `design`, `arm`, `seq_id`, `length`, `sequence` | as above |
| `pin_identity_held` | ⭐ **true in all 312.** The four fixed positions kept their identities `R, H, R, R` |
| `n_cys`, `n_nonstandard` | ✅ **0** and **0** |
| `n_sequon`, `sequon_positions` | N-X-[S/T] motifs. 40.06% carry at least one — ⛔ **inert here**, because the expression system has no glycosylation machinery. Counted for the record only |
| `ncpr` | net charge per residue at pH 7.4. Median −0.0376, acidic |
| `pI` | ⭐⭐ **the column that matters**, because it is the defect an earlier cohort had to be refitted for. Median **5.03**; **pI ≥ 7 in 4.81%** (15 of 312). That earlier run's first attempt exempted only the basic residues from the bias and drove this fraction to 56.26%, into a band where expression is known to fail; exempting the whole charged set keeps it here |
| `longest_run` | longest single-residue run |

### `j904_pin_positions_20261101.tsv` — what was actually fixed, 312 rows

The file that answers "how much of the specified site is identity-constrained?" directly.

| column | meaning |
|---|---|
| `design`, `arm`, `seq_id`, `length`, `sequence` | as above |
| `n_pins` | ⭐ **4 in all 312**, not 6 |
| `pin_positions` | the fixed indices in the designed chain, **read per design** and ⛔ never hard-coded — each batch samples its own contig length, so the positions differ between batches |
| `pin_identities` | ⭐ **`RHRR` in all 312** — three cation groups and the one imidazole |
| `n_cys`, `n_his` | as above |

⛔⛔ **Read `n_pins = 4` against the six groups the mining line specified.** Four are
identity-constrained end to end. The other two sit inside a held backbone element whose sequence was
deliberately released, so they keep their **geometry** and not their chemistry. That is a property of
the construct, explained in `../methods/core_scaffolded_design.md` §4, and it is **not** stated as a
settled design decision anywhere in the material this directory was written from.

⚠️ **The aggregate composition-deviation score is not in these files and is not quoted as a headline
anywhere in this directory**, because it is defined over the *designable* residue set and so changes
when the set of omitted residues changes — which makes a before-and-after pairing that crosses two
different omission lists partly an artefact of the definition. The per-residue percentages in
`../methods/core_scaffolded_design.md` §7.3 are the comparable quantities, and they were re-measured
from these files.

---

## 7. The corrected cohort, and the measurement that justifies it

### `j908_mpnn_sequences_sixpin_20261101.tsv` — ⭐ the cohort that holds all six groups

**351 rows: 312 designed sequences plus 39 inputs.** Same columns as `j904`, produced with the same
settings and on the same 39 backbones, with one change: the **two contained** constellation
positions (scaffold 34 and 38) are given explicit fixed identities in sequence design, exactly like
the four guideposts, mapped through each design's own index map rather than by position.

| what | `j904` | `j908` |
|---|---:|---:|
| imidazole at scaffold 34 | ⛔ **0 / 312** | ✅ **312 / 312** |
| cation at scaffold 38 | 159 / 312 | ✅ **312 / 312** |
| `n_his` | 1 in all 312 | **2 in all 312**, as two specified imidazoles require |
| cysteines | 0 | ✅ **0** |

⚠️ **The contained positions are held at the scaffold's own identities**, which gives **lysine** at
38 where the specification said arginine. Both are cations and the constellation was matched on
functional-group *type*, so the type is satisfied — but it is a real chemistry change, recorded here
rather than absorbed.

⛔ `pin_identity_held` in the QC file still checks the **four** guideposts (`n_pins = 4`), in this
cohort as in the others. It is not evidence about the two contained positions in *either* direction
— that is what the next file is for.

### `j908_seq_qc_sixpin_20261101.tsv` · `j908_qc_input_sixpin_20261101.tsv`

Quality control and its input for the corrected cohort, same columns as the `j904` pair.

### `j953_contained_pin_check_20261102.tsv` — ⭐⭐ the two contained positions, per sequence

**936 rows: all 312 designed sequences of each of the three cohorts.** This is the file that makes
the site-chemistry claim checkable instead of asserted, because no QC column covers these two
positions.

| column | meaning |
|---|---|
| `cohort` | which sequence table the row was measured from |
| `design`, `arm`, `seq_id` | joins to that cohort's own table |
| `motif_offset` | where held element 3 (`FHQFSKFY`) begins in this design's chain. ⛔ **Read per design and never hard-coded** — each backbone samples its own connector lengths, so the element starts at **7 distinct offsets** across the 39 backbones |
| `aa_at_scaffold_34`, `intent_34_imidazole_held` | the residue found where an imidazole was specified, and whether it is one |
| `aa_at_scaffold_38`, `intent_38_cation_held` | the residue found where a cation was specified, and whether it is one (`K` or `R`) |
| `n_his_reported` | the `n_his` value that cohort's own table reports, so the `n_his >= 2` arithmetic can be checked against the identity columns |

Built by `../scripts/j953_contained_pin_check_20261102.py`, which refuses to run if the held motif
does not occur exactly once in a backbone's own input row.

## 8. A note on the `arm` column's prefixes

The `arm` column carries a short run label: bare `opt05` / `opt01` / `opt04` in the screen tables,
`r5_` prefixed in `j904`, `fx2_` prefixed in `j908`. These identify **which sequence-design run**
produced the rows — the production run and the six-position-locked re-run — and the three arms are
the same three core options throughout. ⛔ The prefix is a run label only; it carries no scoring or
selection meaning, and rows join across files on `design` and `seq_id`.

---

# The extension arm

## 9. `nterm_extension_arm_20261102/` — the arm, as actually built

Three files, the same shape as the ten screen arms in §4: the generator **specification**, a
**provenance sidecar**, and the **input coordinates** the specification names.

The specification's one load-bearing field is the contig,
`25-25,E1,2-2,H4-145,/0,A193-298,/0,A332-379`, which reads left to right as the output chain:
**25** generated residues, then **1** held guidepost (chain `E`, its `NE,NH1,NH2` held by
`select_fixed_atoms`), then **2** generated linker residues, then the source design's residues
**4–145** held sequence-fixed, then two chain breaks and the receptor's two held segments.

| provenance key | what it records |
|---|---|
| `design` | the source design's filename — one design from the screen's `opt01` arm |
| `ext`, `gap` | the generated lengths: 25 and 2 |
| `dropped_nterm` | the source design residues the extension replaces: 1, 2, 3 |
| `held_binder_residues` | `[4, 145]`, the held span |
| `total_length` | 170 |
| `threaded_substitutions_vs_design_backbone` | 93 — ⚠️ the one place this construct is not simply a subset of its parent |
| `measured.R1_guanidinium_to_target_A` | 3.03 Å, the held guidepost's distance to the receptor **before** any generation |
| `measured.rank1_draw_worst_target_facing_pin_A` | 11.31 Å — ⛔ in the source draw, the designed contacts are **not** all maintained |

## 10. `j948_nterm_score_CORRECTED_GENSET_20261102.tsv` — 96 designs, 26 columns

⛔⛔ **Read `../methods/nterm_extension_arm.md` §5 before using any column named `new_*`.** An
earlier scoring run measured those columns on a 2-residue linker instead of the 27-residue
generated segment. The files here are the corrected ones; the superseded pair is kept beside them
under `SUPERSEDED` names, with a banner in their own first lines.

| column | what it is |
|---|---|
| `design`, `model` | the design, and which model of its batch |
| ⚠️ `batch` | **the unit of independence.** 12 distinct values over 96 rows; designs sharing a batch are not independent draws |
| `binder_chain`, `target_chains` | resolved from each design's **own** specification, never assumed |
| `n_res_binder`, `n_gen_res` | 170 and 27 in every row. ⭐ `n_gen_res` read **2** in every row before the correction — the constant column that was the tell |
| `n_pseudo_dropped` | generator pseudo-atoms removed before measuring, matched on element **and** name |
| `max_step`, `median_step`, `n_break_4p5` | backbone continuity. `n_break_4p5 == 0` is gate clause 1 and holds in 96/96 |
| `max_ca_dev` | largest alpha-carbon deviation of the held residues from their input positions |
| `rg`, `rg_ratio` | radius of gyration, and that value relative to a length-matched expectation — ⭐ the ratio is the only compactness measure here that is length-safe |
| `bb_min_to_target` | whole binder's closest approach to the receptor. `>= 2.0` is gate clause 2 |
| `n_contact_4p5` | whole-binder contacts with the receptor under 4.5 Å |
| `intra_min_nonadj` | closest approach between non-adjacent residues of the binder itself |
| `new_min_dist` | ⛔ the **generated segment's** closest approach. Gate clause 3, which is **provably redundant**: the segment is a subset of the chain, so this value is always ≥ `bb_min_to_target`, and clause 2 entails it. Verified 96/96 |
| ⭐ `new_contacts_4p5` | the arm's actual answer: contacts the **generated** segment makes with the receptor under 4.5 Å. ⚠️ a sum over an atom-pair distance matrix ⇒ it counts **pairs**, not distinct atoms |
| `pin_min_to_target` | the held guidepost's closest approach |
| `rfd3_*` | ⚠️ the **generator's own** reported statistics, carried for comparison only. Not independent of it |

`…_RANKED.tsv` holds the 66 designs clearing the conjunction, ordered by `new_contacts_4p5`
descending then `bb_min_to_target` descending. ⛔ There is **no verdict column and no binary cut**
beyond the conjunction; the ordering is not a prediction of binding.

## 11. `j949_diverse20_v2_20261102.tsv` — the 20 sequences that were folded

| column | what it is |
|---|---|
| `design`, `seq_id` | the backbone and the sequence index on it |
| `sequence`, `length` | the designed binder chain, 170 residues, cysteine-free by construction |
| `overall_confidence`, `seq_rec` | ⚠️ the sequence-design model's own confidence, and its recovery of the input sequence at designed positions. ⛔ **Neither is a binding predictor** — a previous round measured that confidence does not even separate the design arms — so they were used for reporting and diversity, never as the selection rank |

⛔⛔ **These 20 sit on 3 backbones chosen by the ranking §10's correction invalidated.** Under the
corrected measurement those three rank **4th, 19th and 60th of 66**, with generated-segment contact
counts of **136, 92 and 33**. ⇒ the folded set **spans** the corrected ordering rather than topping
it. `../methods/nterm_extension_arm.md` §7 explains why that was left standing rather than reselected.

## 12. ⭐ `nterm_top3_coordinates_20261102/` — the three structures the figure is built from

The **corrected** top three backbones as the generator produced them, gzipped mmCIF plus each
design's own specification record. ⇒ ⭐ the contact counts in `…_RANKED.tsv`, the panel figure and
the viewer session can all be **recomputed from these three files**, which is the point of
depositing them.

⚠️ **These are raw generator output.** They carry the generator's own pseudo-atoms (named `V` plus
digits, element `V`), which any measurement must remove **by element and name** — a text match on
the name returns zero and reads as a false clean. No relaxation, no minimisation, no cleanup.

Two reproducibility facts read out of the specification records themselves rather than assumed:

- `[M]` **`seed` is `null` in all 96 designs.** ⇒ the generation really is unseeded, so re-running
  this contig does not return these backbones. The coordinates are the record, not a recipe.
- ⚠️ `[M]` **`ckpt_path` is `/weights/rfd3_latest.ckpt`** — a *moving* checkpoint, not a pinned
  version. ⇒ even with a seed, a later run of "the same" generator would not be guaranteed to be
  the same model. ⛔ Which version was in that file on the run date is not recorded anywhere, and
  that is a gap, not a detail.

## 13. `nterm_ecto_fold_inputs_20261102/` — the 20 co-fold inputs, verbatim

The exact 20 input files the co-fold consumed, one per sequence. Each names two protein chains —
`A`, the 170-residue designed binder, and `B`, the **640-residue receptor ectodomain** as one
contiguous chain — plus three `NAG` ligands bonded `ND2 → C1` at the receptor's sequons.

⚠️ Both chains carry `msa: empty`, i.e. **single-sequence mode**. For the de novo binder that is
correct; for the receptor it is the conservative choice and not the best one.

⭐ Every one of these files was machine-checked before it was folded: 2 protein blocks, 3 ligand
blocks all `NAG`, 3 bond constraints, chain ids `A,B,C,D,E` unique and single-character, an
identical 640-residue receptor across all 20, and — the check that matters most — **every bonded
position confirmed to be an asparagine**, because a sugar bonded to the wrong residue folds happily
and says nothing.

## 14. ⭐ `j958_nterm_ecto_fold_scores_20261102.tsv` — the co-fold, 100 rows

One row **per draw**: 20 sequences × 5 diffusion samples. ⛔ Read
`../methods/nterm_extension_arm.md` §8 before using it; three of its columns are easy to misread.

| column | what it is |
|---|---|
| `sequence`, `design`, `draw` | the sequence, its backbone, and which of the 5 draws |
| ⛔ `draw` | **`0` is the BEST of the five, not a sample** — the predictor sorts its diffusion samples by its own confidence. Taking only `draw == 0` moves median interface confidence from 0.282 to 0.559 **by selection alone** |
| ⭐ `pair_iptm_binder_target` | the **binder–receptor** interface confidence, the mean of the two directed values. ⭐ This is the number to read |
| `pair_iptm_a_to_b`, `pair_iptm_b_to_a` | the two directed values, kept because they are not equal |
| ⛔ `iptm_CONFOUNDED` | the **global** interface figure. Named to be awkward on purpose: it folds in the three NAG ligands, so it is not a binder–receptor score |
| `binder_ptm`, `target_ptm`, `ptm`, `complex_plddt`, `confidence_score` | the predictor's other statements about its own prediction |
| ⭐ `gen_contacts_4p5`, `gen_min_dist` | the **27 generated residues** `{1..25, 27, 28}` against the receptor. ⚠️ `gen_contacts_4p5` counts atom **pairs**. ⛔ Its spread **within one sequence across its own 5 draws** has median **130**, against a between-backbone difference in means of **6.3** ⇒ it cannot order anything at 5 draws |
| `all_contacts_4p5`, `all_min_dist` | the whole 170-residue binder against the receptor |
| `pin_contacts_4p5`, `pin_min_dist` | the held guidepost residue alone |
| ⛔⛔ `held_rmsd_to_design`, `n_held_compared` | the **143 held residues** against the design backbone they were supposed to reproduce, within the binder chain so no cross-frame superposition is involved. `[M]` **6.7 to 16.2 Å, median 12.0** ⇒ the predictor rebuilt the protein the construct holds fixed. This caveat sits under every other column |
| `design_frame_gen_contacts` | the same backbone's engagement in the **generator's** frame, carried on every row so the two can be compared without a join |

⛔⛔ **There is no pose column, and that is deliberate.** Positive controls run earlier through this
same instrument — including two natural, deposited complexes, one of whose folds the predictor
reproduces to 0.59 Å — all land **64–85 Å** from their own crystal poses. ⇒ no pose threshold is
admissible on this target, so none is offered.

Built by `../scripts/j958_score_ecto_fold_20261102.py`, which refuses to write if any sequence is
short of its 5 draws — ⛔ because this predictor exits **successfully** having written nothing when
it runs out of memory, so the structure count is the only proof a run finished.
