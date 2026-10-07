# Deposited sequences — the record

⛔⛔ **READ `SEQUENCE_SET_CORRECTION_20261006.md` FIRST.** This file's tables and counts describe
the **four**-sequence set as it stood earlier on 2026-10-06. `id15`, which appears below, was
**deliberately dropped**, and `id86` was added. Nothing below is deleted — it is the record of how
the set stood — but the correction file is the authority.

⛔ **The final submission was 16 sequences, 11 of which passed the novelty filter.** Neither the
four below nor the ten in `submission_challenge1_20261006.csv` is that set; both are earlier states
of it. ⇒ `SEQUENCE_SET_CORRECTION_20261006.md` enumerates all three.

James McCann and the Koder Group, The City College of New York. **2026-10-06.**

Four sequences, two backbones, against the **monomeric tethered EGFR extracellular
region**. ⛔ Superseded: see the banner above. The deposited sequence file
`submitted_sequences.fasta` holds the **ten** of the 2026-10-06 file — not the 16 finally
submitted.

Why the sequences are these: **`DESIGN_CHOICES.md`**.
How they were made: **`methods/README.md`** (pipeline), `methods/OBJECTIVES.md`
(every metric defined), `methods/THRESHOLDS.md` (every constant traced).

⚠️ **No number in this deposit is an affinity.** Every value is epitope contact
geometry or fold self-consistency, computed on predicted structures.
`DESIGN_CHOICES.md` §2 states exactly what that does and does not mean.

---

## The set

| rank | id | arm | aa | backbone | the one measurement behind it |
|---|---|---|---|---|---|
| 1 | `id85` | bridge | 145 | bridge backbone | both footprints in the **same draw**, 20/24 |
| 2 | `id2015` | site A | 75 | site-A backbone | site-A deep label **0.857** over 70 draws |
| 3 | `id15` | site A | 75 | site-A backbone | **25/75** residues from rank 2 (diversity) |
| 4 | `id89` | site A | 75 | site-A backbone | **17/75** residues from rank 2 (diversity) |

**Footprints.** Site A is 11 WT residues: `380 382 384 408 409 410 411 412 417
438 465`. The bridge's second footprint is 7 WT residues: `19 20 21 25 28 29 50`.
Both are asserted by residue **identity**, never by number range, everywhere they
are used.

⛔ **Ranks 3 and 4 carry flags, stated in `DESIGN_CHOICES.md` §5.** They are deposited as a
deliberate diversity bet: ranks 2–4 are **one backbone**, and ranks 3 and 4 are the only
sequence diversity in the set. Dropping them would leave one 145-mer and one 75-mer
scaffold.

---

## How the FASTA was built, and why that matters

⛔ **One sequence WAS typed into the builder, and the exception is the interesting
part.**

**Nine of the ten residue strings** are read out of the same scored table the ranking
was computed from. Before writing, the builder **asserts** each record's expected
length and both terminal residues. On a deposit a single wrong residue is invisible
and unrecoverable, so transcription was removed as a possible failure — for those nine.

⭐ **`id86` is a hard-coded literal in the submission builder, and that was the correct
call.** The scored tables **reuse design ids across different molecules**: the entry
`id86` in `training_c5.csv` is a **75-mer**, while the `id86` we submitted is **80 aa**
— a different protein with the same label. Sourcing it by id would have looked clean,
passed a naive length check against the wrong expectation, and **submitted the wrong
sequence**. Writing it as a literal, with its provenance recorded at the point of use,
was the safer of the two options rather than a shortcut. It is also the one record
whose correctness rests on that provenance note rather than on a table lookup.
⇒ `SEQUENCE_SET_CORRECTION_20261006.md` states the same thing.

⭐ It also **recomputes the site-A pairwise Hamming distances from what it wrote**
(25 / 26 / 17 of 75), so the file cannot silently disagree with the diversity
argument in `DESIGN_CHOICES.md` §4.

---

## ⚠️ Two things about this deposit that are NOT verified

1. ⚠️ **17 of the 621 target residues are unverified against our structural
   model, and whether the published target matches the physical protein cannot
   be checked here.** The 621-residue target is read out of the constant suffix
   of the competition's released `full_sequence` column, where it is a single
   variant shared by 400 of 402 rows — so it is the assay construct **as
   published**. Against our structural model it has been checked
   residue-for-residue at the same WT index: **604 of 604 modelled residues
   match, zero mismatches (97.3% of the target)**. The remaining **17 residues
   are not present in the model** (WT 1-3, 295-299, 307, 577, 615-621) and are
   therefore unverified by that comparison. Whether the published sequence
   matches the physical protein in the tube cannot be verified here and is not
   claimed.
   ⛔ **An earlier version of this item claimed the residue-level check had never
   been done and that only the length arithmetic (645 − 24 signal residues = 621)
   was ever verified. That was false:** the residue-level check exists, is in the
   deposit and reproduces. What that earlier wording got right is that it asked a
   different question — match against the *assay construct* rather than against our
   *model* — and that distinction is preserved above.
2. ⛔ **Two novelty checks HAVE been run, and their headline is the opposite of what an
   earlier version of this item implied.**
   - **Sequence half:** MMseqs2 and phmmer, independently, against Swiss-Prot 2026_03
     and PDB seqres. The nine **submitted** site-A designs returned **zero statistically significant
     hits** by either tool.
   - **Structural half:** a foldseek search against a local PDB structure set, first on
     a chain extracted from a predicted complex and then re-measured on 693 predicted
     **monomer** structures.
   - ⛔ **It also found something we did not expect:** the three 145 aa bridge designs
     are close relatives of the published consensus TPR proteins — which are themselves
     designed proteins. See `accidental_ctpr_homology.txt`.
   ⇒ `NOVELTY_20261006.md` is the index; `NOVELTY_CORRECTION_MONOMER_20261006.md`
   carries the corrections to our own first numbers.

   ⚠️ **Three caveats that travel with that result and must not be dropped:**
   (a) "**No statistically significant hit**" is **not** the same claim as "measured low
   identity" — it means nothing rose above the search's significance threshold.
   (b) We searched **two** of the five databases the competition's annotation uses.
   Additional databases can only raise the maximum similarity found, never lower it, so
   **every identity figure is a lower bound and every novelty figure an upper bound**.
   **Patent sequence collections are the unsearched category we would worry about
   first.**
   (c) We make **no claim about what the competition's own annotation will report**.

   ⚠️ The original point in this item still stands and is unaffected: if novelty is
   assessed at the **backbone** level, the two backbones are the relevant unit, and our
   searches were run per sequence rather than per backbone.

---

## Reproducing a score from the deposit

The per-stage scripts, their inputs, and what each refuses are documented in
`methods/SCRIPTS.md`. ⚠️ **Not everything needed is distributable**: third-party
model weights (structure generator, sequence model, structure predictor) and
third-party scoring code must be obtained from their own upstream projects, and
the raw predicted structures are too large to ship. `methods/INPUTS.md` lists
every required input, marks which are site-specific, and names what substitutes for
each one that is not distributed.
