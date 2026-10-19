# The C3 arm — a symmetric homotrimeric binder for a homotrimeric target

James McCann and the Koder Group, The City College of New York. Campaign `v50`, 2026-10-18..19.
**3 of the 20 submitted sequences come from this arm.** It is documented separately from
`README.md` because the symmetry changes the method at every stage, not just the input.

## 1. Why build a C3 binder at all

TNF-α is a **homotrimer**. A conventional binder engages one protomer and ignores the other two.
A **C3-symmetric** binder is one designed protomer in three copies, arranged on the same
three-fold axis as the target, so the three symmetry-equivalent sites are engaged at once. The
hypothesis is avidity: three copies of a weak interface may hold where one does not.

⛔ **We do not claim it works.** §5 is the result, and it is negative.

## 2. What changes, relative to the C1 arm

| | C1 | **C3** |
|---|---|---|
| binder chains | 1 (`A`) | **3 (`A,B,C`)**, tied — one sequence, three copies |
| target chains | 3 (`B,C,D`) | 3 (`a,b,c`) |
| chain pairs scored | 3 | **15, of which only 9 are binder↔target** |
| designed guideposts | 4 per design | replicated into **every** copy |
| per-copy aggregation | n/a | ⭐ **conjunction, never a mean** |

### 2a. ⛔⛔ The motif must be replicated into every entity
RFdiffusion3 treats a single-copy guidepost in a symmetric spec as **multiplicity 1** and the run
dies on arrival. Every guidepost residue is replicated into all three binder copies, so the
designed salt bridge is specified three times, once per protomer.

### 2b. ⛔⛔ Per-backbone pins, never a cohort-wide list
Each of the 12 backbones has its **own** guidepost residue numbers — `c3_design_plan.json` carries
**12 backbones with 12 distinct pin sets**. A cohort-wide `--pin-resnums` list would be wrong for
11 of the 12 and would silently measure the wrong residues.

### 2c. ⛔⛔ Sequence design must tie the copies, not copy the target
The three binder chains are designed as one tied sequence (`--chains_to_design A,B,C`,
SolubleMPNN, `--omit_AA C`, T = 0.1, 8 sequences per backbone ⇒ **96**). The deposited FASTA holds
**all six chains** per record — the three binder copies **and** the three target protomers. ⛔ Read
chain `A` only when computing binder composition; reading the whole record silently averages the
binder with TNF-α.

## 3. Scoring — `j521_c3_boltz_selfconsistency_20261018.py`

⛔⛔ **A global `iptm` is meaningless here.** With 6 chains there are 15 pairs; **3 are
binder↔binder** (our own designed trimer interface) and **3 are target↔target** (TNF-α's own
crystal trimer, which the predictor rebuilds near-perfectly from a deep alignment). A global
`iptm` on this system is mostly a measurement of TNF-α folding. `j521` extracts
`pair_chains_iptm` **per pair** and reports binder↔target separately from binder↔binder.

⭐ `j521` imports `kabsch_rmsd` and `tm_score` from `j465` under the **same SHA1 pin** `j497` uses,
so C1 and C3 numbers are comparable by construction rather than by assertion.

### 3a. ⭐ The conjunction, and a labelling trap in our own tables
A C3 binder may grip **at most 2 of 3** protomers, so a mean over copies pays partial credit for a
design where one copy does all the work. Every per-copy quantity is reported three ways.

⛔⛔ **`j521`'s `_best`/`_worst` suffixes are `np.max`/`np.min` regardless of direction.** That is
correct for pLDDT / TM / ipTM, where higher is better, and **inverted** for every **distance**
column — `bridge_copy_best_*`, `rmsd_binder_copy_*`, `pae_*` — where lower is better. ⇒ **for a
distance, the conjunction (the worst copy) is the `_best` column.** Verified: `_best > _worst` in
**96 of 96** rows. Selection used `bridge_copy_best_best`, the worst copy's best bridge. ⚠️ A
reader using these tables must apply the same correction; the numbers are right and only the
suffixes mislead.

### 3b. ⛔ Chain permutation — the binder chains too
The three binder copies are identical in sequence and the predictor assigns chain IDs with no
guaranteed order, so a positional A→A pairing is a guess. `j532` enumerates **all 6** permutations
for the target frame and `j536` does the same for the binder assembly. **59 of 96** C3
trimer folds fit best on a **non-cyclic** pairing — a positional assumption would have mis-scored
most of the cohort. The same defect, found in the C1 arm first, had wrongly voided **34 of 88** C1
designs; all 34 were recovered, and the correction reproduces the original exactly where it was
right (**54/54 unchanged at 0.0000 Å**).

### 3c. ⭐ `j536` — the assembly fit, which per-copy RMSD cannot see
`j521` superposes each copy on its own design chain **independently**. Kabsch is a free rigid
motion, so that measures the **protomer fold** and is blind to the **arrangement**: three
perfectly folded protomers in the wrong twist score identically to the designed trimer. `j536`
fits **all three binder chains at once as one rigid body**, so the residual contains the
arrangement — which is the thing a C3 design actually specifies.

## 4. The unbound reference, and one wrong question we asked first

The C3 binder is designed as an **obligate homotrimer**, so a lone protomer is a state it was
never meant to occupy. We folded one anyway (`c3_apo`) and it scored badly — apo pLDDT median
**0.428**, apo→design TM median **0.285**. ⛔ That is the **wrong question**: for an obligate
oligomer a lone protomer may be molten by construction. The relevant unbound state is the
assembled binder with no target (`c3_trimer`).
⛔⛔ **C1's apo pLDDT must never be compared to C3's** — same column name, different physical
question. C1 apo is a monomer that *is* the whole binder; C3 apo is a fragment of one.

## 5. The result — negative, and clearly so

The right question was then asked, and it gives a **worse** answer than the wrong one:

| | pLDDT median | TM-to-design median |
|---|---|---|
| `c3_apo` (one protomer alone) | 0.428 | 0.285 |
| **`c3_trimer`** (the designed assembly) | **0.326** | **0.242** |

⇒ assembling the three copies **did not rescue the fold**. TM > 0.7: **0 of 96**.

**The assembly is not recovered at all.** Assembly RMSD median **23.5 Å** (min 12.1); **0 of 96**
come within 1 Å of their own per-copy fit, so **no design has its protomers arranged as designed**,
and the arrangement adds ~6 Å on top of an already-unrecovered protomer.

**Against the target** (`c3_holo`, 24 folds): complex pLDDT median 0.548. The designed Arg→Glu107
bridge, as the conjunction over copies: median **13.8 Å**, best **2.02 Å**. **9 of 24** have all
three copies within 9 Å and **3 of 24** within 4 Å. ⚠️ Engagement is **all-or-nothing** — every
design is 0 or 3 copies, never 1 or 2 — which is what the symmetry predicts and is therefore not
independent evidence of binding.

⛔ **A caveat we can measure, and one we cannot.** Measurable: ranking C3 on the bridge is
ranking inside structures whose fold is not recovered, so a short bridge may be an accident of a
collapsed chain lying against the target. Across all 24, the rank correlation between assembly
RMSD and bridge is only **−0.165** — weak, true at the extremes only. Not measurable here:
⛔⛔ **we have no positive control for this predictor on a de novo symmetric homotrimer with no
alignment.** Without one, TM 0.24 cannot be split between *the designs fail* and *the predictor
cannot score this class of object at all*. That control is the first thing we would run next.

### 5a. ⚠️ A composition result we could not explain
The C3 designed sequences are compositionally unlike the C1 ones: Gly+Ser **0.321** vs **0.175**,
hydrophobic **0.255** vs **0.352**, top-3-residue fraction **0.455** vs **0.328** (natural globular
proteins run ≈ 0.12 / 0.32 / 0.22). A sequence that low in hydrophobic content has little to build
a core with, which would explain the folding result.
⛔ **The obvious cause is refuted.** We expected the symmetric backbones to specify a smaller
buried core; they do not — RFd3 marks **20.0%** of C3 positions as core against **20.6%** for C1.
⇒ the mechanism is **open**. Our current suspicion, untested, is the three-way tying: a residue
must satisfy three environments at once, which may push the optimiser toward small flexible
residues. We state it as a suspicion, not a finding.

### 5b. What was checked rather than assumed
- Provenance: every submitted C3 sequence **matches a record in the deposited FASTAs exactly**,
  which a reader can check here. ⚠️ That the *folded* chain matched was checked against the
  prediction files, which are **not deposited** (size) — that step is not reproducible from
  this repository alone.
- `c3_design_plan.json` resolves **12 backbones, 12 distinct pin sets**.
- Every submitted C3 sequence is **cysteine-free** and its three binder chains are **identical**,
  so each submitted chain genuinely represents the homotrimer that was designed.
- Predictions used **one diffusion draw per sequence**, so no spread is measured and none is
  claimed.

## 6. Files

```
structures/designs_v50_c3/                    12 designed backbones, .cif.gz, pin column intact
data/submission_20261019/c3_design_plan.json  12 backbones, 12 distinct pin sets
data/submission_20261019/seqs_c3/             12 FASTAs, 8 designs each (96); all 6 chains per record
data/submission_20261019/scoring/             c3_apo / c3_trimer / c3_holo, + both assembly fits
scripts/j521_c3_boltz_selfconsistency_20261018.py   the C3 scorer
scripts/j532_frame_permutation_fix_20261018.py      6-permutation target frame
scripts/j536_c3_assembly_fit_20261019.py            the assembly fit
scripts/j534_select_submission_c1_c3_20261018.py    the two-arm selector
figures/c3_designed_backbones_20261019.png          all 12 designs
figures/c3_predicted_complexes_20261019.png         the 3 submitted, predicted
```

⛔ The predicted complex structures themselves are not deposited (size). The scoring tables carry
every number quoted above, and the scripts that produced them are here.
