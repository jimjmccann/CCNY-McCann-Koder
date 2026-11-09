#!/usr/bin/env python3
"""s115 -- select the REFOLD cohort from the 2x2, one design per diffusion batch, at RANDOM.

⛔ READ-ONLY. Writes only a cohort list. Nothing launches, spends, or edits a table.

⭐⭐ THE SELECTION RULE, DECLARED BEFORE THE REFOLD DATA EXISTS:
  ⑴ **Conjunction survivors only** -- `continuous AND pass_raw_clash AND one_fold`, `s51`'s own
     gate as applied by `s75` in `t149`. ⛔ A clashing or broken backbone cannot test engagement.
  ⑵ **ONE design per diffusion batch.** `n` is BATCHES ([[rfd3-designs-are-clustered-by-batch]]):
     all 4 designs in a batch share one sampled contig, so taking 2 from a batch would be
     pseudo-replication -- the error that nearly shipped a 3.4x enrichment resting on one backbone.
  ⑶ **At RANDOM within the batch, seeded.** ⛔⛔ NOT "the best" by any design-stage score. Every
     available quality criterion is correlated with the ENDPOINT: picking the largest
     `bb_min_to_target` selects designs sitting FURTHER from the target, which biases peptide
     contacts DOWN, and differentially if the arms' distributions differ; picking on `cont_pep`
     selects the endpoint itself. ⭐ Random is the only rule provably unbiased with respect to every
     endpoint, and it is applied IDENTICALLY to all four arms.
     ⚠️ The cost is accepted and stated: this cohort tests the ARMS, it is ⛔ not a hand-picked
     candidate shortlist, and the best design in the round is probably not in it.

⚠️ THE COHORT WILL NOT BE PERFECTLY BALANCED, and that is a RESULT, not a defect. A batch with zero
  survivors cannot contribute, and `[M]` the motif arms lose ~28% of designs at the raw-clash gate
  against ~10% for the no-motif arms (`t150` §3). ⇒ the realised n per arm is reported here and must
  be carried into the refold comparison rather than assumed to be 20/20/20/20.

Usage:
  python3 -I s115_select_refold_cohort.py <t149.tsv> <out_cohort.txt>
"""
import collections
import csv
import random
import sys

SEED = 20261108

if len(sys.argv) != 3:
    sys.exit(__doc__)
REF, OUT = sys.argv[1], sys.argv[2]

rows = []
with open(REF) as fh:
    for r in csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"):
        rows.append(r)
if not rows:
    sys.exit("⛔ no rows read -- refusing to write an empty cohort")

surv = [r for r in rows
        if int(r["continuous"]) and int(r["pass_raw_clash"]) and int(r["one_fold"])]

by_batch = collections.defaultdict(list)
for r in surv:
    by_batch[(r["arm"], r["batch"])].append(r)

rnd = random.Random(SEED)
picked = []
for key in sorted(by_batch):
    cand = sorted(by_batch[key], key=lambda x: x["design"])     # deterministic order first
    picked.append(rnd.choice(cand))

print(f"[M] {len(rows)} designs, {len(surv)} conjunction survivors, "
      f"{len(by_batch)} batches with >=1 survivor, {len(picked)} selected (1 per batch)")
print()

arms = collections.defaultdict(lambda: [0, 0, 0])
nb = collections.defaultdict(set)
for r in rows:
    nb[r["arm"]].add(r["batch"])
for r in surv:
    arms[r["arm"]][1] += 1
for r in rows:
    arms[r["arm"]][0] += 1
sel = collections.Counter(r["arm"] for r in picked)

print(f"{'arm':<22} {'designs':>8} {'survivors':>10} {'batches':>8} {'batches w/ surv':>16} "
      f"{'SELECTED':>9}")
for arm in sorted(arms):
    bws = len({b for (a, b) in by_batch if a == arm})
    print(f"{arm:<22} {arms[arm][0]:>8} {arms[arm][1]:>10} {len(nb[arm]):>8} {bws:>16} "
          f"{sel[arm]:>9}")

print()
lost = {arm: len(nb[arm]) - len({b for (a, b) in by_batch if a == arm}) for arm in sorted(arms)}
print(f"[M] batches contributing NOTHING (all 4 designs failed the gate): "
      + ", ".join(f"{a} {n}" for a, n in lost.items()))
print(f"[M] cohort is {'BALANCED' if len(set(sel.values())) == 1 else 'NOT balanced'} across arms "
      f"-- carry the realised n into the comparison")

with open(OUT, "w") as fh:
    for r in sorted(picked, key=lambda x: x["design"]):
        fh.write(r["design"] + "\n")
print()
print(f"wrote {OUT} -- {len(picked)} design basenames, seed {SEED}")
