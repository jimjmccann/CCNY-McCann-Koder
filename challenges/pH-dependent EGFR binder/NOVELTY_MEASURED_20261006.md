# Novelty measured on the actual submission candidates

Measured 2026-10-06, on local CPU.

> **Corrected on the monomer fold — read `NOVELTY_CORRECTION_MONOMER_20261006.md` before quoting any
> structural number below.** Every structural figure in this file was measured on **chain A extracted
> from a predicted complex**. The organizers fold the submitted sequence **alone**. Re-measured on 693
> monomer structures, **0 of 693 draws reaches TM 0.80** and the single highest draw is **0.738**.
> Specifically, the **`egfr-siteA-id2015` row below, flagged `YES` HIGH at TM 0.885, does not
> reproduce** (monomer: 0.590 median, 0/63 draws ≥ 0.80) — that flag was an extraction artifact. And
> **not one** of the nearest-neighbour PDB codes in the tables below survives the switch to the
> monomer fold. The file is kept as written, because its numbers are what the reasoning of the time
> was built on.

This was the first evaluation of either half of the Level-4 novelty criterion on the molecules
actually being submitted. Both halves are measurable locally: the structural half always was, and the
sequence half is too, once the scored object is the Boltz-folded chain A rather than the RFd3
generator backbone.

## 0. The headline: there are two novelty classes, and the set mixes them

- **The nine 75–80 aa site-A binders pass the sequence half** (20.7–28.2% identity) and fail only the
  structural half, placing them at Level 2–3 by the framework's own worked example.
- **The three 145 aa bridges fail both halves, and fail them badly:** sequence identity
  **59.0–64.2%** (the bar is ≤ 30%) against **high** structural similarity (TM **0.967–0.977** over
  more than 70% of the design; high is TM ≥ 0.80). They sit one axis away from **Level 1 (known)**,
  which is > 70% identity plus at-least-moderate structure.

## 1. Method, and the two traps avoided

`scripts/novelty_chainA_foldseek.py`. `foldseek easy-search` against a local foldseek PDB database,
`-e 10 --max-seqs 2000`.

1. **Chain A only.** Chain B is EGFR; including it matches EGFR itself and reads as "not novel" for
   the wrong reason. An earlier run established this.
2. **Structures matched by sequence, never by id.** Design ids are reused across campaigns — `id86`
   and `id47` each appear in four different arms, and an `id2015` directory exists under a *bridge*
   arm holding a different 145 aa molecule. 42,850 prediction directories were indexed by their
   chain-A sequence; each submitted sequence was matched exactly, and the written chain-A PDB was
   **re-read and asserted equal to the submitted sequence** before use.

An earlier run scored RFd3 generator backbones, whose sequence is ~52% alanine and is not what gets
submitted. These are Boltz folds of the **submitted** sequence, so `fident` is the real identity of
the real sequence. Identity is quoted at **coverage > 0.50**, so a 24-residue alignment cannot
inflate it.

## 2. The twelve candidates measured

| design | len | struct (TM@cov>.7) | HIGH? | seq-id @cov>.5 | seq half | nearest |
|---|---|---|---|---|---|---|
| `egfr-siteA-id3036` | 75 | 0.646 | no | **20.7%** | ✅ pass | 4n06 |
| `egfr-siteA-id2015` | 75 | **0.885** RETRACTED | YES RETRACTED | 22.5% | ✅ pass | 4n06 |
| `egfr-siteA-id3017` | 75 | 0.745 | no | 23.8% | ✅ pass | 8j9b |
| `egfr-siteA-id2046` | 75 | 0.682 | no | 24.0% | ✅ pass | 4n06 |
| `egfr-siteA-id2048` | 75 | 0.757 | no | 24.0% | ✅ pass | 8j9b |
| `egfr-siteA-id3020` | 75 | 0.791 | no | 25.3% | ✅ pass | 8j9b |
| `egfr-siteA-id89` | 75 | 0.745 | no | 25.6% | ✅ pass | 8jo4 |
| `egfr-siteA-id3035` | 75 | **0.567** | no | 25.9% | ✅ pass | 3t49 |
| `egfr-siteA-id86` | 80 | 0.723 | no | 28.2% | ✅ pass | 8ucf |
| `egfr-bridge-id47` | 145 | **0.967** | **YES** | **59.0%** | FAIL | 7obi |
| `egfr-bridge-id1104` | 145 | **0.975** | **YES** | **59.5%** | FAIL | 5a01 |
| `egfr-bridge-id85` | 145 | **0.977** | **YES** | **64.2%** | FAIL | 7obi |

**The `id2015` row above is retracted** (monomer fold): 0.885/HIGH does not reproduce — monomer 0.590
median, 0/63 draws ≥ 0.80 — and its `4n06` nearest neighbour becomes `2lse`. The row is left in place
because later reasoning was built on it.

**Structural half: 12/12 defeated** (TM ≥ 0.50 over more than 70% of the design), consistent with an
earlier 95.1% on n = 400. This is a property of the all-α fold class, not of these designs.
**`id3035` is the structurally least-defeated molecule in the set** (0.567, barely over the moderate
line) and `id3036` the most sequence-novel (20.7%).

## 3. The bridge scaffold is a consensus TPR repeat protein — the nearest neighbours agree

`id85`'s top sequence hits, all with TM 0.95–1.00: **2avp** 67.6% · **1na3** 64.2% · **2wqh** 61.0% ·
**3kd7** 60.2% · **2hyz** 55.3% · **7obi** · **5a01**. The submitted bridge sequence carries the
canonical repeat motif
`…DPNNAEAYYNRGNVYAFAGKYEEAIKDYEKAIKLDPNFAEAYYNLGDTYAEAGKYEEAIEYYEKAIKL…`, a
**tetratricopeptide-repeat consensus**. ⚠️ These are PDB codes and a motif; the published titles of
those entries were not independently confirmed.

**Why this is expected rather than a bug:** the TPR consensus *is* the energetically optimal sequence
for a TPR backbone, so a sequence designer run on a repeat backbone converges toward it. It remains a
novelty liability, and it is measured rather than inferred.

All three bridges share the backbone `arm_shard11_shard11_41_model_5` — **one scaffold, one novelty
verdict.** Adding bridge molecules cannot diversify away from it.

## 4. Four further candidates, where the binding-evidence ordering inverts on novelty

| candidate | len | struct | HIGH? | seq-id | seq half | novelty verdict |
|---|---|---|---|---|---|---|
| `cand-siteA-id2` | 75 | 0.714 | no | **27.5%** | ✅ pass | **best profile measured** |
| `cand-siteA-id73` | 80 | 0.751 | no | 29.5% | ✅ pass | passes, close to the line |
| `cand-siteA-id19-sh09` | 75 | **0.953** | **YES** | **37.1%** | FAIL | fails both halves |
| `cand-bridge-id76` | 145 | **0.964** | **YES** | **58.9%** | FAIL | the TPR problem |

The candidates ranked best on binding evidence are the novelty-worst, and the ones ranked weakest are
the novelty-best. `id19` had been favoured as the only candidate on a fourth backbone and as the
lowest-identity molecule in the set at 60.0% — but that 60.0% is identity **to our own other
designs**, a *diversity* number. Its identity to the **PDB** is 37.1%, which fails the gate actually
applied. Two different quantities were doing each other's work.

## 5. What this does not show — read before acting

1. **The bar is novelty ≥ 3/4, not 4/4.** Failing Level 4 is **not** automatic rejection. Levels 2
   and 3 are not precisely defined in the published framework; only L4 (≤ 30% *and* less than
   moderate) and L1 (> 70% *and* at least moderate) are. **No claim is made here that any design will
   be rejected.** The sound claim is ordinal: the bridges fail both halves where the site-A binders
   fail one, so the bridges are strictly worse and sit nearer L1.
2. **One database.** The competition annotates against **pdb100 + afdb50 + cath50**; this searched a
   local `pdb` only. More databases can only find more hits, so every number here is an **upper bound
   on novelty**, never a floor. `afdb50` especially will contain more helical-bundle neighbours.
3. **Different pipeline.** Theirs is ESMFold, then FoldSeek and TM-align with a **consensus of three
   domain predictors**; ours is Boltz chain A, then foldseek alone. Domain-level coverage could
   differ.
4. `fident` is identity over the aligned region, not global; MMseqs2 as they run it may differ.
5. The annotation pipeline is **not necessarily** the same thing as the selection criteria, and
   novelty — a computed label — and the written rationale — a scored selection input — are two
   different quantities.

## 6. What the measurement implied for the set

The result quantified a **binding-versus-novelty trade** on both axes for the first time, and it cut
against adding more bridges: they are the worst novelty class in the set and they add no scaffold
diversity, since they share `id85`'s backbone. Against that, the bridge is the only architecture that
engages **both** footprints, and `id85` is the measured champion on binding evidence — so dropping
the bridges entirely would trade a real measurement for an unvalidated novelty mapping.

The molecules that pass the sequence half and are not in the submitted set are `id2` and `id73`;
their binding evidence is weaker. `id19`'s case rested on backbone diversity, and it fails both
novelty halves (37.1% identity, TM 0.953).
