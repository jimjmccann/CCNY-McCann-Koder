# Conditional TNF-α binder — a negative, unvalidated result

James McCann and the Koder Group, The City College of New York.
Design work 2026-10-07 .. 2026-10-18. Challenge 2 of the Anthropic × Adaptyv 2026 competition:
*"Design a conditional TNF-α binder."*

---

## ⛔ Read the limitations first. They are the result.

1. **Nothing here has been tested at the bench. No number in this directory is an affinity.**
   Every coordinate is a model — either generator output or a structure-predictor output — and is
   labelled as such. None was measured.

2. **Nothing was submitted for this challenge.** There is no submitted sequence set, no
   novelty-passed set and no competition result to report. The 110 sequences deposited here are
   candidates that were never entered.

3. ⛔⛔ **The designed interaction is not reproduced by our own structure predictor.** The whole
   affinity hypothesis in these designs is one salt bridge: a pinned arginine on the binder to
   **Glu107** of a TNF-α protomer. The designs specify it at **3.45 Å in 69 of 69** scored
   predictions' input structures. Re-folding the designed sequences against the whole trimer with
   Boltz-2 reproduces it in **0 of 69** — median **21.5 Å**, nothing under 4 Å, and only 4 of 69
   under 9 Å. The binder also lands **34–70 Å (median 43 Å) from the site it was designed for**,
   on a trimer only ~50–60 Å across, so it is a different face of the target and not a near miss.
   §"The result" has the numbers and how they were measured.

4. ⛔ **But "the designs are wrong" is NOT established, and we are not claiming it.** This is one
   predictor's verdict. **No positive control exists** — we have no case where this predictor
   confirmed one of our own designed interfaces, so its sensitivity on a de novo interface is
   **unmeasured**. The same predictor is already on record in this group's work misplacing exactly
   this kind of feature (right fold, wrong side-chain rotamer, bond absent) and building a
   different target conformation than the one asked for. "Our predictor cannot see our interfaces"
   is as consistent with these numbers as "our interfaces are not there." Separating those two
   requires the positive control, which was not run.

5. ⚠️ **The scored set is incomplete: 69 of 110 folds.** The remaining 41 were still running when
   this was written, and the per-backbone counts are unbalanced (4–9 of 10). A truncated run has
   reversed a ranking in this group's work before, so no per-backbone claim here should be treated
   as settled.

6. ⛔ **A confidence score we previously quoted was confounded, and the retraction is part of the
   record.** On a trimeric target, the global `iptm` averages 6 chain pairs of which **3 are
   target–target** — TNF-α's own crystal trimer, which the predictor rebuilds near-perfectly from
   a deep alignment. Global `iptm` median **0.691** reads as a good interface; the binder-to-target
   pairs alone are **0.240** while target-to-target is **0.965**. ⇒ `iptm ≈ 0.7` on this construct
   means *TNF-α folded*, not *our binder bound*. An earlier internal reading of these same files
   was corrected on this basis. The deposited scoring tables carry the per-pair columns so the
   confound is checkable rather than asserted.

7. ⛔ **Two of the three things the challenge asked for are not addressed here.** The challenge
   asks for a binder that (a) binds human TNF-α, (b) **also binds mouse** TNF-α, and (c) **holds at
   neutral pH and releases at acidic pH** — and ranks the pH switch first, mouse second, affinity
   third. This cohort carries **one pinned histidine** per design as a pH trigger; our own
   pre-registered objective for this challenge specified a conjunction over **three** histidines.
   The switch was **never measured** on these designs, and **mouse cross-reactivity is not
   addressed in this cohort at all**. We state the descope rather than quietly re-scoping the
   challenge to the part we worked on.

8. ⚠️ **No novelty check was run.** In challenge 1 the organizers' novelty filter removed 5 of our
   16 submitted sequences. Nothing in this directory has been checked against known proteins, so
   these sequences carry no novelty claim in either direction.

9. ⛔ **The compute-orchestration and launch layer is deliberately absent**, as it was for
   challenge 1. It carries account identifiers and cost figures, which are not published. Any
   batch system reproduces what is here; `methods/README.md` says what has to happen on the
   compute side without the machinery.

---

## What this is

A de novo binder campaign against the soluble human TNF-α trimer, targeting the receptor-binding
site between two adjacent protomers. The generator is **RFdiffusion3**, conditioned on a small
rigid motif ("guideposts") rather than on hotspot residues: three arginines placed to bridge
Glu107 of a protomer, plus one histidine as the intended pH trigger. Sequences come from
**SolubleMPNN**; the fold-back check is **Boltz-2** against the whole 3×152-residue trimer.

This round's new method, relative to challenge 1, is the **guidepost** route — conditioning the
generator on our own placed residues and then gating the resulting backbones on whether they are
*one connected chain* and *one fold*, before anything else is measured. Challenge 1's deposit
holds the hotspot-conditioned route, the target-sphere preparation and the full-receptor exclusion
gate; those are not duplicated here.

## What is deposited

| | |
|---|---|
| `structures/designs_v46_joint11/` | the **11 designed backbones** that passed both backbone gates, as generator output, with the guidepost flags intact |
| `data/mpnn_20261018/` | **110 sequences** (10 per backbone), the pinned-residue record, and the run settings |
| `data/scoring_20261018/` | the **scored fold tables**, 103 columns, including every per-chain-pair confidence column |
| `scripts/` | the 12 programs that built, gated, sequenced, folded and scored the above |
| `methods/` | the method in reader order, per-script notes, and what is not distributed |

⛔ **Not deposited:** the predicted complex structures themselves (35 MB+ for an incomplete set),
the raw generator output beyond the 11 that passed (384 designs), model weights, and the
third-party alignment. `methods/INPUTS.md` says what to obtain and from where.

## The result

All figures below were **recomputed from the deposited tables** for this README, not copied from a
lab note. Reproduce them from `data/scoring_20261018/`.

**The pinned bridge** — minimum distance from the *pinned* arginine's guanidinium nitrogen to that
protomer's Glu107 carboxylate oxygen, heavy atoms, best of the three protomers:

| set | n | design | prediction | prediction under 4 Å |
|---|---|---|---|---|
| first round + probe | 11 | 3.45 Å (3.45–3.45) | median **19.08 Å** (4.77–32.76) | **0 / 11** |
| all folds scored so far | 69 | 3.45 Å (3.45–3.45) | median **21.52 Å** (4.03–52.40) | **0 / 69** |

⭐ **The design side being exactly 3.45 Å everywhere is a consistency check, not a result.** The
guidepost is a rigid template transplant, so the design file essentially restates the request —
and the deposited tables show it: the designed distance is identical to the two decimals printed
across **11 different backbones of 9 different lengths**, which is what a held template looks like
and not a property any of those backbones achieved. That is precisely what makes the prediction
column a fair comparison: the *ask* was identical in all of them. ⛔ So 3.45 Å is not evidence of
affinity.

⛔ **Restricting to the pinned arginine matters.** Measuring *any* arginine gives a median of
**13.97 Å** instead of 19.08 Å on the first round — so where a bridge-like distance does appear, it
is usually a *different*, non-pinned arginine that the sequence model happened to place, not the
designed guidepost. The designed mechanism is not what occasionally forms, and measuring "any
arginine" would have overstated the result. Both columns are in the tables.

**Placement, fold, and confidence** (69 scored folds):

- binder displacement from its designed position, superposing on the target only and not refitting
  the binder: median **43.4 Å**, range 34.2–70.3 Å over the **45 rows whose target frame is
  trustworthy**;
- **24 of 69 rows (35%) are flagged `FRAME_SUSPECT`** — the predictor built the TNF-α *trimer
  assembly* itself differently (target RMSD ~21.8 Å rather than ~1.5 Å). Their placement columns
  must not be read. This is a fact about the predictor, not about the designs, and it is
  per-sequence rather than per-backbone;
- `tm_binder` (does the predictor reproduce the backbone's own fold) median **0.493**, range
  0.219–0.947. ⇒ **"right fold, wrong place" is a real and separate category here.** Two designs
  reproduce their fold well and are still 37–47 Å out of position, so a fold-only self-consistency
  gate would have passed them and said nothing about binding;
- binder-chain-alone `plddt` median **0.590** against a `complex_plddt` of **0.792** — the complex
  figure is the target's number with ours diluted into it;
- `pae` binder→target **17.1–26.2 Å**. A confident interface cannot carry 20 Å of expected
  positional error. This is an independent channel agreeing with the per-pair `iptm`, and it is why
  we treat the global `iptm` as the outlier rather than the per-pair numbers.

**Backbone choice did not predict binding.** Shuffling the backbone labels 4000 times against the
between-backbone spread of per-backbone medians separates **only `tm_binder`** (SD of medians
0.125, p = 0.007). Placement (p = 0.26), interface `iptm` (p = 0.26) and the bridge (p = 0.37) are
**flat**. ⚠️ Four metrics were tested, so the one surviving p-value needs reading with that in
mind, and it says only *"the predictor reproduces this backbone's own fold"* — not *"this backbone
binds."* The four best bridges sit on four *different* backbones. ⇒ on this evidence, sequence
choice dominates the binding metrics and backbone choice does not; picking a backbone on these
numbers would be picking noise.

**On a deliberately lenient three-way gate** — bridge < 9 Å and `tm` > 0.6 and placement < 40 Å,
trustworthy frames only — **exactly 1 of 69 passes**: bridge 7.1 Å, `tm` 0.753, placement 35.7 Å.
Its interface `iptm` of 0.167 sits in the bottom eighth of the cohort (rank 8 of 69 from the
bottom; the cohort minimum is 0.141), and 35.7 Å is still a different face of the target. ⇒ it is
the least-bad fold, not a hit.

## What would decide it

A **positive control**: push a known TNF-α binder's sequence through the identical path and see
whether the predictor reproduces *its* interface. Without that, the sensitivity of this test is
unmeasured and these numbers cannot distinguish a failed design from a blind instrument. It is
cheap and local and was not run before this deposit was assembled.

## A note on internal references

Several deposited scripts carry identifiers like `J490`, `ADR 0030` or `decisions/0055`. Those point
into this group's internal note series and decision log — an append-only record where each entry
states one settled question, the measurement behind it, and the date. **They are not deposited.**
They are kept rather than stripped because they are the *provenance of a threshold*, and a threshold
with its provenance removed reads as an arbitrary constant. ⇒ Read them as "this value was settled
deliberately, on a stated date, and the record exists." Where a sentence conveyed nothing at all
without the reference, the sentence was rewritten instead.

⚠️ **Two qualifications, so the convention is not oversold.** Some of these identifiers name a
*program and its outputs* rather than a written note — `J465`, `J494`, `J497` and `J498` are of that
kind, and for `J494` and `J498` no note was ever written, so there is nothing further to read even
internally. And where the identifier refers to something that *is* deposited, it resolves inside
this directory: `J497`'s record is the scoring tables in `data/scoring_20261018/`, and `J490`'s is
the gate tables in `data/backbone_gates_20261018/`.

The marker glyphs (⛔ ⭐ ⚠️) are emphasis, carried over from the working documents as in challenge 1.
They rank nothing and encode nothing.
