# Design choices — what we deposited, and why

⛔⛔ **READ `SEQUENCE_SET_CORRECTION_20261006.md` FIRST.** Where this file speaks of **four**
deposited sequences including `id15` — §5, and the closing summary — it describes the set as it
stood earlier on 2026-10-06. **Ten** were submitted, `id15` was **deliberately dropped** and
`id86` added. The reasoning in those sections is unchanged and still the reasoning behind the
submitted set; only the count and the membership moved.

James McCann and the Koder Group, The City College of New York.
De novo binders against the monomeric tethered EGFR extracellular region.

`methods/README.md` is the **pipeline** — the eight stages in run order, and every
threshold traced to its source. This file is the **decisions**: why the deposited
molecules are these molecules and not others, what the numbers behind them do and
do not mean, and what we got wrong.

Read `methods/README.md` first if you want to reproduce the pipeline. Read this
file first if you want to judge the designs.

⚠️ **Nothing here has been tested experimentally.** Every number is computational,
from one structure generator and one structure predictor, on one target.

---

## 1. Three architectures were attempted. Two produced a molecule.

The target's two candidate faces are **11.46 Å** apart at closest approach, with
footprint centroids **25.25 Å** apart and a maximum span of **45.59 Å** (measured
on the tethered model before anything was built). That geometry is what the three
architectures are about.

| architecture | the idea | state |
|---|---|---|
| **Site A** | one binder on the domain-III face, 11-residue footprint | ⭐ **built, optimised through six cycles, plateaued** |
| **Bridge** | one chain long enough to engage **both** faces at once | ⭐⭐ **built, 145 aa, and it is the strongest thing we have** |
| **Loop insert** | take a site-A binder and splice a short element that reaches the second face | ⛔ **the element was never generated.** Not "untested" — it does not exist |

⛔ **The loop insert is a generation gap, and we say so rather than omit it.** The
splice anchor can only reach three residues of the second face, and the published
geometry for it is **one-way Euclidean reach**, while an insert is a *loop* and
spends **contour length**. A lower bound on the residues needed per target residue
(so an upper bound on feasibility) puts the two residues worth having at **5
residues of insert each**, robust over clash tolerance 3.0–4.0 Å, glycans in or
out, and a stricter contact proxy. That is the next thing to build; it is not a
stretch goal, and it is not in this deposit.

⚠️ A side-by-side **two-domain** construct was also screened and is geometrically
feasible — of 27,156 candidate pairs, **5,193 are connectable** and 294 need
**≤ 4 residues** of linker, with three of four folded second-face backbones
reachable (best linker **3**). ⛔ **No fusion has ever been folded as one chain**,
so there is no two-domain sequence here either — only the geometry saying one
could be built.

---

## 2. ⛔⛔ WHAT OUR LABELS MEASURE — read this before reading any number

Every ranking in this programme is built from **two** quantities, both computed on
predicted structures:

1. **Epitope contact geometry** — for each predicted draw, how many footprint
   residues the binder touches within 5 Å, asserted by residue **identity**, never
   by number range.
2. **Fold self-consistency** — `sc_rmsd`, the binder-onto-binder Cα RMSD against
   the designed backbone, used as a **gate**, plus `pose_rmsd` (placement) as a
   diagnostic that is never a cut.

⛔⛔ **Neither is affinity, and nothing in this programme predicts affinity.** The
labels say *"in this predicted structure, the binder is in contact with the
intended residues, and it built the shape it was designed to build."* They cannot
say how tightly it binds, or whether it binds at all.

This is not a hedge — it is the measured state of the art. No method predicts de
novo binding free energy better than about 1.9 kcal/mol RMSE, which is roughly a
25-fold error in dissociation constant, and **on de novo designs predictor
confidence does not correlate with affinity at all.**

⇒ Three consequences we hold ourselves to:

- ⛔ **A contact count is not a ranking of goodness.** Ranking backbones on
  interface contact alone named six winners of which **four failed the fold
  gate**; the top pick had a median Cα fold error of **8.45 Å**, so each of its
  eleven "full-footprint" draws was contact against a *different* structure.
  Arm-wide, gating on fold removes **53% of full-footprint contact (422 draws →
  200)**, and 16 of the 31 failing backbones had ≥5 full-footprint draws — the
  contact metric was *actively promoting* unbuildable backbones. **Gate on fold
  before ranking on contact.**
- ⛔ **A large contact count can be the bad sign.** One deposited sequence makes
  the **most** contacts of its cohort (mean 107, up to 172 in one draw) and scores
  **zero** on every interface-packing metric. A high count with no packing is what
  a spread, non-specific interface looks like — and it is precisely the pattern a
  footprint **count** cannot see.
- ⛔ **Both of our ranking metrics saturated.** At the end, two designs sat at the
  metric ceiling with a perfect fold rate and were **not separable**: a single
  design's full-footprint fraction carries s.e. ≈ **0.16** at ten draws. We could
  not raise draws to fix it — ten draws already measured **97.6%** of a 24 GB
  card's memory on this construct, and the predictor's out-of-memory path **exits
  successfully having written nothing**. Precision had to come from *repeating*
  the ten-draw measurement, not from enlarging it.

---

## 3. Why the bridge sequence is `id85`

The bridge objective is a **per-draw floor**: `min(n_A/11, n_B7′/7)` if the draw
folds, else **0**, meaned over **every** draw the predictor produced.

Three choices inside that one line, each of which we got wrong first:

1. ⛔ **A misfolded or mispositioned draw scores 0. It is never dropped.** Scoring
   zero and not reporting a score are different things; a draw that folded is a
   data point whatever it says. Dropping non-qualifying draws preferentially
   deletes the *low* second-site values — measured, mean second-site contact is
   **4.15** when the first site reaches ≥9/11 and **2.81** when it does not — so
   conditional averaging inflates exactly the designs that rarely qualify. One
   design scored a perfect conditional mean from a **single** qualifying draw while
   its three draws were [7, 6, 0].
2. ⛔ **The floor, not the average of the two sites.** Averaging each requirement
   separately rewards a design that never once satisfies both at the same time: one
   design ranked near the top on the separated average and had **zero** both-site
   draws. The `min` is same-draw by construction.
3. ⛔ **Not a count of draws that satisfy everything.** At three draws a count can
   take at most **four** values, which collapsed a 50-sequence table to a single
   non-zero row. Measured side by side: 50 distinct values on the mean, **1** on
   the count.

**The result, over 24 fresh draws at depth:**

| id | deep label | both footprints in the same draw | fold rate | verdict |
|---|---|---|---|---|
| ⭐ **`id85`** | **0.973** | **20/24 = 0.833** | 1.000 | the bridge's sequence |
| `id47` | 0.687 | 16/36 = 0.444 | 1.000 | loses by **9.30 s.e.** |
| `id76` | 0.368 | 9/36 = 0.250 | 0.861 | loses by **11.61 s.e.** |

Three things make that comparison trustworthy, and one weakens it:

- ⭐ **An in-batch anchor with a known prior value was carried in the same round.**
  It returned **+0.037**, i.e. **1.75 s.e.** of the *combined* uncertainty of the
  two estimates — within tolerance ⇒ no batch shift, the round is comparable to
  the one before it.
  ⚠️ We first computed that against the new round's s.e. **alone**, got 3.55 s.e.,
  and declared the round void. That was wrong: the prior value is itself an
  estimate, and ignoring its uncertainty makes any well-replicated new round look
  like a batch shift. The corrected arithmetic is the one above.
- ⭐ **The backbone was separable before any sequence was scored.** Of eight
  candidate bridge backbones, this one produced both-site draws at **8/109**
  against **0/804** for the rest, Fisher **p = 3.3e-08**.
- ⭐ **An unrelated measurement family agrees.** An interface-packing panel
  (LIS/iLIS, ipSAE, actifpTM, pDockQ, interface pLDDT — a completely different
  axis from contact counts) puts them in the **identical** order,
  `id85 > id47 >> id76`, with `id76` below the one published threshold in that
  family that applies. Two unrelated methods, same answer.
- ⚠️ **One anchor cannot separate a batch shift from a change in that one
  sequence.** The previous round used two, precisely so they could move in
  opposite directions. This round's control is weaker by construction, and the
  anchor sat **+0.037 high** — below threshold, not zero, and in the direction
  that would flatter it. It does not decide anything here only because the margins
  being judged are **8–16× larger**.

⇒ **`id85` is deposited. The bridge decision is closed**, and we would not reopen
it on the evidence available without an assay.

---

## 4. Why site A contributes **one** molecule and two diversity bets

Six cycles of a closed-loop optimizer ended with **seven designs tied at the top**
and no statistical test separating them. The obvious reading is "our metric is too
noisy." ⛔ **It is the wrong reading, and one cheap measurement settles it.**

All seven are **75 aa**, so Hamming distance is exact. Measured:

> **The seven "tied" leaders differ by 1–4 of 75 residues (mean 2.2).**

They are **one molecule with point mutations**. Nothing separates them because
**there is nothing there to separate** — not because the measurement is too coarse.

**Farthest-point selection** over the whole site-A pool, seeded at the best design:

| added | minimum Hamming distance to the set so far |
|---|---|
| `id2015` | — (seed; highest site-A deep label, **0.857** over 70 draws) |
| ⭐ `id15` | **25 / 75** |
| ⭐ `id89` | **17 / 75** |
| next candidate | ⛔ **4 / 75** — diversity is exhausted after three |

⇒ ⭐⭐ **Site A offers three distinct molecules, not ten.** Each further slot spent
on the leading family buys **≤ 4 residues** of difference, which is not diversity.
All three sit on **one** backbone.

### Why the loop stopped rather than running another cycle

Every one of these is measured, and each independently removes a reason to
continue:

- **Enumeration is complete.** The combination walk over the accepted mutation set
  finished at **88** combinations (50 at orders 2–3 plus 38 at orders 4–6; no order
  7 exists), and the last enumerated cycle produced **zero** sequences beating the
  previous best.
- **Both ranking metrics saturated** at their ceilings (§2).
- **A saturation scan at the one position that governs the rate-limiting contact
  found nothing**: best **0.609** against its own parent's **0.600**, against a
  noise floor of **0.158**, and **0 of 42** scan positions beat the parent.
- ⛔ **There is no model in the loop, and that is measured, not a preference.**
  Out-of-bag RMSE/sd(y) was **1.14** and **1.04** on the two labels — at or worse
  than predicting the mean — and upper-confidence-bound acquisition **did not beat
  random** (+0.394 ± 0.634, t = +0.62, Fisher p = 1.000). Trained on the saturated
  label it predicted **11.696** and **12.101** against a physical ceiling of **11**.
  ⇒ Do not reintroduce acquisition-based proposal without re-measuring
  out-of-bag error first, and require **ratio + 1 s.e. < 1.0** — a ratio of 0.999
  once printed as "1.00" and was read as a pass.
- **The limiting contact is architectural.** One footprint residue is contacted in
  only **63.8%** of 600 draws and **63.7%** of a later 380 fresh draws — weakest of
  the eleven by ~13 points, where the strongest is 98.8% — and of draws missing the
  ceiling by **exactly one** residue, it is that residue. The accepted set the
  optimizer walked contains **no position that reaches it.**

⇒ **A further site-A cycle needs a new premise, not another iteration.** That is
why the deposit is three site-A molecules and not a seventh cycle.

---

## 5. The four deposited sequences, with the flags on two of them

⛔ **SUPERSEDED, 2026-10-06.** Ten sequences were submitted, not four. `id15`, ranked 3 below, was **deliberately
dropped** before submission, and `id86` was added. See
`SEQUENCE_SET_CORRECTION_20261006.md`, which is authoritative and gives the submitted ten in
order. The per-design reasoning and the flags below are unchanged and are why the surviving
designs were chosen.

Ranked as deposited. ⛔ **The flags are not hidden and are not softened.**

| # | id | arm | aa | the measurement behind it | flags |
|---|---|---|---|---|---|
| 1 | **`id85`** | bridge | 145 | both footprints in the **same draw** in 20/24, fold rate 1.000 | — |
| 2 | **`id2015`** | site A | 75 | highest site-A deep label **0.857** over **70** draws, fold rate 1.000 | — |
| 3 | **`id15`** | site A | 75 | **25/75** distinct from #2 — the real diversity axis | ⛔ **three independent flags**, below |
| 4 | **`id89`** | site A | 75 | **17/75** distinct from #2 | ⛔ interface screen; ⚠️ only 10 draws |

### ⛔⛔ `id15` is flagged by three unrelated kinds of measurement

| axis | `id15` | the leading family |
|---|---|---|
| strict fold-**and**-full-footprint rate | **0.178** | 0.80–0.86 |
| interface packing (iLIS) | **0.000 in 10/10 draws** — computed, not missing | 0.486 |
| backbone geometry vs the crystal target in the **same file** | ⛔ **τ 4.29×**, heavy-atom clashes **1.98×** — the only candidate above its own target on **both** | 0.13× / 1.80× |

⚠️ Its averaged label (0.779) sat mid-pack and its label **stability** was the best
of the eight (sd 0.038). **Those were measuring consistency, not quality.** This is
the single most useful cautionary result in the deposit: a low-variance mid-pack
average survived two of our own ranking metrics and failed three others.

⚠️ Honest limits on that verdict: all three axes read the **same** predicted
structures, so a systematic predictor artefact would move them together; and the τ
ratio rests on **10** models.

### ⭐ `id89` is a different kind of risk, and that is why it is still here

`id89`'s backbone is **better than its own crystal target** on both τ (0.85×) and
clashes (0.54×), with **zero** buried unsatisfied charges. Its interface score of
0.000 is therefore about **where and how it sits on the target**, not about whether
the chain is built properly. ⇒ a *placement* risk, not a *foldability* risk — a
better diversity bet than `id15` if only one of the two is kept.

### ⭐ `id2015` is the cleanest sequence in the programme on every axis measured

Best or tied-best on backbone τ (**0.13×** its own target), buried unsatisfied
charges (**0.02 per 100 aa** — one across 70 models) and exposed apolar area
(**4.32 Å² per residue**), on the largest model count of any candidate (70), and
first on the interface panel (iLIS **0.486**, ipSAE **0.598**, actifpTM **0.940**),
stably so across all seven replicates (0.414–0.515).

⚠️ ipSAE **0.598** is the closest anything here comes to the published **0.61**
threshold. ⛔ That threshold was calibrated on a different predictor and on natural
protein-protein discovery, so it **must not be lifted**: this is stated as
proximity, **not** as a pass.

⚠️ **One comparison here has no precedent in our record.** The interface panel is
the first measurement that places a site-A design above the bridge on any axis;
every contact-based comparison kept them in separate tables because they are
scored against different footprints. ⇒ treat the `id2015`-versus-`id85` **order**
as untested. The finding is that **both clear the field**, not which of the two is
first.

### ⛔ What the deposited set is NOT diverse in

`id2015`, `id15` and `id89` are **one backbone**; `id85` is a second. The deposit
is **two backbones**, four sequences. If a novelty or diversity criterion is
applied at the backbone level, that is the number that matters.

---

## 6. ⚠️ Two caveats that are properties of the measurement, not of the designs

1. **The interface panel discards, it does not rank.** Its own documentation puts
   this target near the **bottom of fifteen** for filterability (average precision
   ~0.2–0.25 against a 0.57 median), so the right expectation is
   **3–6× enrichment, not separation.** And `pDockQ` does not discriminate on this
   cohort at all (0.248–0.394, with the worst sequence *above* two better ones) —
   ⛔ quoting it in support of any pick would be selecting the one metric that
   disagrees.
2. **The post-fold geometry panel is relative, not calibrated.** It is **not** a
   MolProbity clashscore and **not** a MolProbity Ramachandran test — the hydrogen
   placement and contour grids for those were unavailable. It is a heavy-atom clash
   count and a coarse favoured-region test, each defined explicitly in the script,
   and each reported **as a ratio to the crystal-derived target in the same file,
   scored by the same code**. ⛔ No calibrated threshold for any of it exists, and
   inventing one is how a screen quietly becomes a fake gate. Hydrogens were not
   placed, so the buried-unsatisfied-charge counts are a **floor**, not a count.

---

## 7. The honest result, in one place

- The **aggregate** progression is real: between-design variance **3.02×** and
  **2.42×** the binomial null (χ² p = 2.5e-11 and 1.3e-06), and the population
  per-draw full-footprint rate went **24.8% → 53.1%** (p = 5.7e-18) between two
  cycles.
- ⛔ **No individual design is resolvable** on our own metrics (§2).
- ⛔ **Part of the apparent per-cycle gain was selection, not improvement** —
  measured. Re-folding six designs spanning the score range gave a fresh-on-prior
  slope of **+0.343**: both champions at 1.000 returned 0.700 and 0.800, and a
  deliberately included **low** anchor came **up**, 0.400 → 0.600. Without a low
  anchor this is invisible; re-measuring only champions can confirm them and can
  never measure regression to the mean.
- ⛔ **The campaign ran four cycles before anything was ever re-measured.** The two
  rows that looked like replicates were **byte-identical across three successive
  training tables** — carried forward, never re-measured. Replication is now
  standing: the top ~5 of every cycle are re-folded with fresh draws **inside the
  cycle's own slots**, which costs no extra compute and controls for batch effects
  because repeats and new designs are predicted in the same run.
- ⛔ **The rate-limiting contact is architectural** and the explored substitutions
  do not reach it (§4).
- ⛔ **An entire earlier round was voided** because the generation target was a
  contiguous crop: 70% of first-face and 100% of second-face designs sat within
  2 Å of receptor that had been deleted — backbone overlap, not contact. The
  window had been validated by **inclusion** ("all footprint residues present") and
  never by **exclusion**, and every downstream stage then treated the window **as**
  the target, so the omission was unfalsifiable from inside the pipeline. The
  replacement is a sphere cut plus a **full-receptor exclusion gate**; full account
  in `ideas_that_didnt_work.txt` entry 1, and it is the most generally useful
  lesson in this deposit.

⇒ **This is a working closed-loop design pipeline with a measured progression, a
measured noise floor, and four deposited molecules. It is not a finished binder,
and no number here is an affinity.**

---

## 8. What would settle it

In order of how much they would change the picture, cheapest first:

1. **Calibrate the expression-associated geometry panel against measured
   outcomes.** A set of 402 prior designs with measured experimental results is
   held and has **never** been used to calibrate any of it. No new compute. ⛔ Not
   done.
2. **Run the full-receptor exclusion gate on the folded complexes** of the
   deposited sequences. It has been run on the site-A *backbones* (146/146 pass)
   and never on these complexes. It also names which glycosylation sequon a design
   sits in, and three sequons lie within 1.5 Å of the site-A target. ⚠️ A distance
   to a modelled glycan is a **floor**: the crystal resolves 1–7 sugars where a
   native tree is ~15–20 residues. No new compute. ⛔ Not done.
3. **Re-fold `id89` at depth** so it stands on the same footing as the others
   (10 draws against 70 and 90). This costs compute.
4. **Generate the loop insert** (§1). It does not exist.
5. ⭐⭐ **The assay.** Nothing computational settles the ranking; the only
   measurement that does is the experiment.
