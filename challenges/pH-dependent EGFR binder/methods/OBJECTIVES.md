# OBJECTIVES AND METRICS — every definition, and why it is what it is

Site-A EGFR binder campaign. Companion to `README.md`; constants are traced in
`THRESHOLDS.md`.

⛔ **Two objectives coexist in the record and they rank differently.** Kendall
tau between them is **0.707** over 142 designs / 1,520 draws, and one design sits
at rank 63 under the older metric and rank 7 under the newer one. **Label the
metric on every number you quote.** §1 is the metric every landed table was
scored on; §2 is the metric in force at the end of the campaign.

---

## 0. The per-draw primitives

Everything is built from three per-draw quantities. A "draw" is one structure
prediction of one sequence.

| symbol | definition |
|---|---|
| `n_epitope` | how many of the 11 site-A footprint residues the binder contacts: **min distance over all atom pairs < 5.0 Å**, chain A = binder, chain B = target. Numbering **asserted by residue identity**, not by range |
| `c380` | **1** if that same contact test puts F380 in the contacted set, else **0**. ⛔ Not a new measurement — it is `380 in hit`, out of the *same* geometry pass, so the two can never disagree |
| `sc_rmsd` | binder-onto-binder Cα RMSD, superposed on the binder alone ⇒ **fold** fidelity |
| `pose_rmsd` | superpose on the mapped target, then binder Cα RMSD with **no** binder re-alignment ⇒ **placement** fidelity. `pose_rmsd ≥ sc_rmsd` always |

**"The draw folds"** means `sc_rmsd < 1.5 Å` everywhere in this document.

⛔ `pose_rmsd` is a **diagnostic and never a cut** (ruled on after an 8× measured
swing). Epitope contact is the primary quantity; the fold gate is the filter.

⚠️ The contact geometry is **copied, not re-derived**, between the scorer and the
per-residue breakdown script, and both say so in their headers. Two independent
derivations of "the same" geometry disagree silently. A different cutoff makes
every table incomparable with every table already produced.

---

## 1. The metric the landed tables were scored on

Per sequence, pooling its draws:

    label             = mean n_epitope over GATE-PASSING draws (sc_rmsd < 1.5)
    gate_frac         = n_gated / n_draws
    passes_floor      = gate_frac >= 0.5
    mean_nep_all      = mean n_epitope over ALL draws
    full_epitope_frac = fraction of ALL draws with n_epitope == 11      ("fef")
    median_sc         = median sc_rmsd over ALL draws
    mean_tgt          = mean n_target_contacts over ALL draws

⛔ **`passes_floor` is a CANDIDATE FILTER, never a term summed into the label.**
Mixing a quality floor into the score is how a score stops meaning anything.

⛔ **The label definition is copied verbatim** from the cycle-0 inline script and
the burn-in builder. Any drift silently changes what the surrogate is trained on,
and the earlier cycles would stop reproducing.

### Why this metric stopped working — both halves saturated

- **`label` has a hard ceiling of 11**, the footprint size. MEASURED at cycle 1:
  the top three sat at **identical means**, dmean +0.00, Fisher p=0.628–1.000.
- **`fef` then hit its ceiling of 1.000** on two designs, both passing the fold
  gate 10/10, not separable by anything available.
- A single design's `fef` carries **s.e. ≈ 0.16** at 10 draws
  (√(p(1−p)/10) at p=0.5), so 1.000 vs 1.000 vs 0.900 is **inside the noise**.
- **Draws cannot be raised to fix it**: 10 draws already measured **97.6% GPU
  memory** on a 24 GB card for this 75-residue construct, and out-of-memory exits
  rc=0 writing nothing. ⇒ Precision must come from **repeated** 10-draw runs.

⇒ Saturation is foreseeable from the definition. Check for a ceiling before
building the loop.

---

## 2. The current objective — a per-draw floor that REQUIRES F380

    c380           = 1 if the draw contacts F380, else 0
    per-draw score = min(n_epitope / 11, c380)   if the draw folds (sc_rmsd < 1.5)
                   = 0                             otherwise
    label_A        = mean of that over EVERY draw the predictor produced

For the two-site bridge the same term nests: `min(n_A/11, c380, n_B7p/7)`.

### Why F380 goes INSIDE the min

The requirement is *"hold the footprint **and** F380 in the same draw"*. Putting
`c380` inside the `min` enforces that by construction. Nothing is averaged
per-residue or per-region first.

⛔ **Never a sum and never a per-residue average.** MEASURED on real draws:

- one design holds every one of the 11 residues in **≥50% of its draws** and
  covers all 11 together in **1 draw in 10**;
- on the bridge, one sequence scores **0.857** if each site is averaged
  separately — near the top of the table — and has **ZERO** both-site draws; the
  per-draw min gives it 0.571;
- Kendall tau between the two rankings is **+0.783**, so this changes the
  **order**, not just the scale.

⇒ Marginal coverage overstates same-draw coverage. Averaging marginals rewards a
design that never once does the thing we want.

### Why the mean is over EVERY draw, including misfolds

A misfolded or mispositioned draw contributes a real **0**. *Scoring zero is
different from not reporting a score*; a misfolded protein is still a data point.
This is the one place this objective deliberately differs in **shape** from
`label` in §1, which averages over gate-passers only. ⛔ Do not reconcile them.

The only admissible missing value is a draw where the **predictor produced no
structure at all**.

### Why HARD and not a weighted F380 term — and this reversed the first call

The expectation going in was a **soft** floor, on the argument that zeroing ~36% of
draws would recreate a resolution collapse. **That argument was wrong, and the
measurement is why.** Over 98 sequences with per-draw data:

| per-draw term | distinct values / 98 seqs | zeros | median | max |
|---|---|---|---|---|
| `fef` (all-11, §1's metric) | **11** | 4.1% | 0.500 | 1.000 |
| `n_epitope/11`, no F380 requirement | 34 | 0 | 0.936 | 1.000 |
| **`min(n_epitope/11, c380)` — adopted** | **31** | **0** | **0.600** | 1.000 |
| soft floor 0.5 on F380 | 36 | 0 | 0.755 | 1.000 |
| soft floor 0.7 on F380 | 56 | 0 | 0.850 | 1.000 |

⭐ **The hard floor collapses nothing.** Because the term is averaged over 10
draws, a sequence almost never misses F380 in *all* of them — so there are **no
zeros at all**, it yields **31 distinct values against `fef`'s 11** (~3× the
resolution), and the median sits at 0.600 rather than `n_epitope/11`'s
near-saturated 0.936. Hard is simultaneously the cleanest statement of the
requirement *and* the better-resolved metric. The soft variants resolve slightly
more and buy it by weakening the requirement, with no mechanistic justification.

⚠️ **ERRATUM carried from the decision record, not corrected away:** the decision
document's own evidence table was computed **without** the misfold-zero rule,
while the objective text includes it. On the same 98-sequence cohort the
no-misfold-zero variant reproduces that table exactly; the objective **as
written** gives 30 distinct values and median 0.509, because 99 of 980 draws
(10.1%) fail the fold gate. The decision is unchanged and the conclusion holds
either way — the hard floor still collapses nothing. **The implementation follows
the objective text (misfold = 0).**

⚠️ **SCOPE: this table is measured at 10 draws.** The two-site arm runs at 3,
where 21 of 50 sequences score exactly 0. ⛔ **Do not assume the no-zeros property
transfers to a 3-draw arm** — it depends on averaging over enough draws that an
all-miss sequence is rare.

### What this objective does and does not fix

⭐ At the **top** of the table it reduces **exactly** to the F380 contact rate —
63 of 98 sequences, and all ten re-measured top designs. ⇒ **It does not
unsaturate the leaders. Only replication does** (`README.md` §8).

⭐ In the **middle** of the table it does real work: pooled re-score of 142
designs / 1,520 draws gives **53 distinct values against `fef`'s 14**, Kendall
tau 0.707, and `fef` had buried one design at **rank 63 where this objective puts
it 7th** — a design that had been picked independently on other grounds.

---

## 3. Why F380, and what kind of claim that is

F380 is the **rate-limiting** contact of the 11, not merely a necessary one.
MEASURED:

- contacted in **63.8%** of 600 draws from two cycles plus the replicate round,
  and **63.7%** of a later cycle's 380 fresh draws. **Weakest of the 11 by ~13
  points** — the strongest, H409, is 98.8%;
- of draws missing the ceiling by **exactly one** residue, F380 is that residue
  in **89.8%** and **81.9%** of the two cohorts;
- ⭐ the rate **did not move** when the campaign's mutation order went 1–3 → 4–6
  ⇒ it is **architectural**, not a property of the substitutions explored;
- ⛔ the accepted set the optimizer walked contains **no position that reaches
  F380** — which is why cycles 1–3 enumerated all 88 combinations and the last
  beat nothing.

**One binder position effectively governs it:** within 5 Å of F380 in 248 of the
310 draws that contact it (**80%**); the next position is second at 7.7%.

Causal evidence that it governs: substituting that position to glycine takes the
score **0.40 → 0.00** while **folding better** than its parent (median fold 0.587
vs 0.612, gate 10/10, floor PASS). 18 of 19 substitutions there are unfolded.

⛔ **Governing is not the same as improvable.** A saturation scan at that position
found **no substitution above the noise floor**: best 0.609 against its own
parent's 0.600, floor 0.158, and the leading champion at 1.000. The identity at
that position was not the lever; the pose was.

⚠️ The 380-reaching geometry is measured on **one backbone family's ~500 draws**.

---

## 4. Backbone-level metrics

Backbones are not scored like sequences.

    pool EVERY draw of a backbone across shards, ASSERT the total,
    then report the pass rate at SEVERAL thresholds, not one.

⛔ Backbones **span shards** (the fold set is length-sorted and shards are
sequential chunks), so a backbone's draws can split across two shards. A backbone
silently scored on a partial draw set reads **too high**, and that has already
reversed a winner on this project.

⛔ Thresholds are reported at several levels rather than one, **because the
threshold is a choice and should be made in the open.**

⛔⛔ **Gate on fold quality before ranking on contacts.** A contact-only ranking
named six backbones and **four of the six fail the fold gate**; the top pick had
a median Cα fold error of **8.45 Å**. Arm-wide the gate removes **53% of
full-footprint contact: 422 draws → 200**, and 16 of the 31 failing backbones had
≥5 full-footprint draws — the contact metric was actively promoting unbuildable
backbones.

⚠️ **Sequence-model confidence is a BETWEEN-backbone signal, not a within-backbone
one.** MEASURED: rho = −0.299 against median fold error across 1,082 sequences
campaign-wide, and **absent** within a single backbone at n=60. A between-group
predictor used as a within-group ranker is a new, untested claim.

⚠️ And no draw count settles the backbone question: within-backbone sequence
variance is **10×** the between-backbone variance, and the leading site-A
backbones are statistically inseparable.

---

## 5. Metrics we tried and demoted or retracted

### 5a. Footprint engagement vs measured binding — RETRACTED as an epitope claim

Folding 402 prior competition designs with published outcomes against our face
gave **AUROC 0.633** (CI excluding 0.5, p=0.0055), beating all 21 confidence
metrics available. The **pre-registered decoy control** — the same designs folded
against a face almost nobody targeted — scored **higher**, 0.656; paired
difference **−0.022, 95% CI [−0.135, +0.094]**.

Stratifying then revealed the pool is a mixture: 99 cysteine-rich ligand-mimic
designs (binding rate 29.3%) and 149 genuinely designed interfaces (11.4%). The
mimics engage the control face more and ours less, so pooling inflated the
control and suppressed ours — **Simpson's paradox**. In the designed stratum ours
survives a length control (**0.715, p=0.0045**) and the control face does not
(**0.557, p=0.45**), with a quartile dose–response of 2.7% → 5.4% → 15.8% →
21.6%.

⛔ **But the direct contrast still fails:** paired difference within the designed
stratum **+0.087, 95% CI [−0.097, +0.264]**. *"A passes its own test and B fails
its own"* is weaker than *"A beats B"*. The epitope-specific claim stays
retracted; the AUROC itself reproduces byte-for-byte. What was withdrawn is what
it **meant**.

**Three standing rules from it:**
1. Quote that metric **stratified or not at all** — state the cysteine stratum
   and whether length was held.
2. It is kept under the weaker name *"docks tightly to the face it is handed"*.
   **It ranks; it sets no absolute threshold.**
3. ⛔ **It does not become a design objective until binder length is regressed
   out.** Length alone scores 0.575 and correlates with engagement (+0.249 ours,
   +0.379 control); length is a variable we set at design time, so optimising
   engagement would partly optimise "make it longer" — Goodhart, not affinity.
   The stratified amendment makes this worse, not better: length is an
   **independent** predictor in that stratum (0.700 p=0.0055 raw; 0.676 p=0.0157
   after regressing out engagement).

### 5b. Gate on one site, grade on the other — REPLACED

Proposed for the bridge: keep draws that fold **and** reach ≥9/11 at site A, then
average site-B contact over those. Three measured defects:

1. ⛔ **It selects on the outcome.** corr(n_A, n_B) = **+0.590** (p=1.9e-15) over
   150 fold-passing draws: mean site-B is **4.15** when n_A ≥ 9 and **2.81** when
   not. Dropping non-qualifying draws preferentially deletes the **low** site-B
   measurements. The label became an average over a sample chosen by the thing
   being averaged.
2. ⛔ **The inflation lands where it does damage.** One sequence scored **7.00**
   from a single qualifying draw while its three draws were [7, 6, 0] — true mean
   4.33, inflation **+2.67**. 10 of 36 labels rested on one draw, 4 of the top 10.
   Kendall tau between biased and fixed rules is **+0.714** and that sequence
   falls from 1st to outside the top 5.
3. ⛔ **The gate discarded data it should have scored.** 14 of 50 sequences got no
   label at all, including one that reached site B 7/7 in a draw.

⭐ Under the replacement, the separate "≥50% of draws must qualify" floor becomes
**redundant**: the 12 sequences it would have rejected score exactly 0.00 anyway.
Gate, floor and label collapse into one honest number.

### 5c. Counting the draws that satisfy everything — ARITHMETICALLY IMPOSSIBLE

| candidate objective | usable seqs | distinct values | at zero |
|---|---|---|---|
| gate-and-grade (5b) | 36 | 17 | — |
| mean site-B over folded draws | 50 | 20 | — |
| **count of strict both-site draws** | 50 | **1** (all zero) | 50 |
| count, relaxed | 50 | **4** | 28 |
| **per-draw min floor — adopted** | **50** | **33** | **7** |

A count over 3 draws can take at most **4** values. That is a hard ceiling on
resolution, not a tuning problem.

⇒ **Tabulate distinct values / zeros / median / max for every candidate objective
on real data before adopting it.** A metric with fewer distinct values than you
have sequences cannot rank them.

### 5d. The rare true criterion — WATCHED, never regressed on

Where the thing actually wanted is too rare to optimise, optimise the graded
proxy and keep the rare criterion as its own column.

⛔ And know they can **disagree about which backbone to take**. MEASURED: one
backbone scored far better on the graded label (**0.464 vs 0.195**,
Mann-Whitney p=7.2e-06) and had **never once** satisfied both sites in 150 draws,
while the other was the only one that ever had — **8/109 vs 0/804** across the
rest, Fisher p=3.3e-08. The graded label rewards moderate contact across many
draws; the criterion rewards rare complete events. That choice was made by a
human, on the record, and not by the label.

---

## 6. Surrogate gate — the one number that keeps a model out of the loop

    GATE:  out-of-bag RMSE / sd(y)  +  1 s.e.  <  1.0

A ratio ≥ 1.0 means the model is at or worse than predicting the mean.

MEASURED, site A: **1.14** on `label`, **1.04** on `fef`. UCB did not beat random
(+0.394 ± 0.634, t=+0.62, Fisher p=1.000). Trained on the saturated label it
predicted 11.696 and 12.101 against a ceiling of 11.

MEASURED, bridge: fails at n=50 (1.056) and n=60 (**1.028 ± 0.019**, Spearman
+0.136, AUROC 0.579). The learning curve, by subsampling, 3 repeats per size:

| n | OOB RMSE/sd | sd |
|---|---|---|
| 20 | 1.077 | 0.080 |
| 30 | 1.031 | 0.122 |
| 40 | 1.064 | 0.133 |
| 50 | 1.077 | 0.024 |
| 60 | 1.028 | — |

**FLAT** across a 3× range of n, and n=60 sits inside n=30's scatter. ⇒ The
apparent gain from pooling an earlier campaign was noise, and **more sequences is
not a path to a working surrogate**; the label at 3 draws carries too little
information per sequence. An earlier conclusion that "n is the constraint" is
withdrawn by that curve. A large part of the label is simply whether the
predictor folds the design at all (`sc_rmsd` rho −0.506), and predicting that
from sequence is the hard part.

⚠️ **Check zero-inflation every cycle before trusting any model ranking** — 21 of
50 sequences at exactly 0 on one backbone, 7 of 50 on another.

⚠️ **The ratio needs its standard error.** A value of 0.999 printed as "1.00" was
once read as a pass; the verdict flipped on the third decimal.

**A free pre-spend check:** ask whether the label is learnable at all with a ridge
+ random-forest floor on data already held, before paying for any prediction. A
negative there is the informative result — the loop would stall and you would
learn it after paying. ⚠️ With its own caveat recorded: held-out burn-in
sequences sit ~21 mutations from the training set while the optimizer's proposals
sit 1 mutation from a parent with data nearby, so local ranking can be easier
than the global test. Do not read a weak global r as "the loop cannot rank single
mutants".

---

## 7. Noise floor — the number every comparison must clear

MEASURED, first in-cycle replicate run, n=6 spanning the score range:

    replicate  source   prior_fef  fresh_fef   delta
    id4901     id3017       1.000      0.700  -0.300
    id4902     id3020       1.000      0.800  -0.200
    id4903     id2015       0.950      0.800  -0.150
    id4904     id2007       0.900      1.000  +0.100
    id4905     id3035       0.900      0.900  +0.000
    id4906     id41         0.400      0.600  +0.200   <- the LOW anchor

    n=6   mean delta -0.058   mean |delta| 0.158   max |delta| 0.300
    Pearson(prior, fresh) +0.556   slope +0.343

⇒ **Any difference smaller than ~0.158 is not a difference.** The cycle's best
new design beat its own parent by 0.009.

⇒ **Slope +0.343 is regression to the mean**: part of every cycle's apparent gain
was **selection**. The low anchor is what makes that visible — a champions-only
repeat round would have shown the fall and not explained it.

⛔ **Apply the floor to the baseline as well.** The honest baseline is the hub
design at **10.05 on 20 draws, fef 0.40** — not its single best reading.
