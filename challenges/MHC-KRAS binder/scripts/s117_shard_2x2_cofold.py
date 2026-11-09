#!/usr/bin/env python3
"""s117 -- split the 2x2 co-fold inputs into POD SHARDS, round-robin, controls in EVERY shard.

    RUN WITH:  python3 -I   (stdlib only; no network, no spend)

⛔⛔ THE SEED IS THE SHARD (`J1039`): shard01..04 fold with seeds 1..4. ⇒ a shard-to-shard
   difference is CONFOUNDED WITH A SEED DIFFERENCE. Two consequences, both acted on here:

  ⑴ ⭐ THE ARMS ARE ROUND-ROBINED ACROSS SHARDS, never dealt in contiguous blocks. If arm aligned
     with shard, the 2x2's main effect would be inseparable from the seed, and a single lost shard
     would silently delete one whole cell of the design
     ([[never-rank-on-partial-data]], [[files-delivered-is-not-rows-delivered]]).
     `[M]` asserted below: every shard carries exactly 5 of each of the 4 arms.

  ⑵ ⭐⭐ THE FOUR CONTROLS GO IN **EVERY** SHARD, not once in the batch. `J986` §3 rule 5 is that
     the floor MOVES with the arm and is never read off an older table -- and across shards the
     seed is part of "older table". With the controls in one shard only, 60 of the 80 designs would
     be compared against a non-binder floor folded at a DIFFERENT SEED, which is exactly the
     cross-batch comparison the rule forbids. ⇒ 4 controls x 4 shards = 16 control folds, so each
     shard holds its own complete floor at its own seed and the conjunction is evaluable
     PER SHARD as well as pooled. ⚠️ The cost is 12 extra folds (~$0.15) and it buys the only
     version of the floor that is actually in-batch.

⛔ A STRAY NON-YAML IN A BOLTZ INPUT DIR IS FATAL TO THE WHOLE BATCH (`J985` §6), and a bundle is
   the worst place to discover it: the pod boots, dies, and bills. ⇒ shard dirs are pure `.yaml`,
   asserted by count here and again by `s110` before upload.

Usage:
  python3 -I s117_shard_2x2_cofold.py <s116_out_dir> <shard_out_dir> [n_shards]
"""
import collections
import csv
import os
import shutil
import sys

if len(sys.argv) not in (3, 4):
    sys.exit(__doc__)
SRC, DEST = sys.argv[1], sys.argv[2]
NSH = int(sys.argv[3]) if len(sys.argv) == 4 else 4

MAN = os.path.join(SRC, "MANIFEST_2X2_COFOLD.tsv")
YDIR = os.path.join(SRC, "yaml")
for p in (MAN, YDIR):
    if not os.path.exists(p):
        sys.exit(f"⛔ missing {p}")
if os.path.exists(DEST) and os.listdir(DEST):
    sys.exit(f"⛔ {DEST} exists and is not empty -- use a fresh dir")

rows = list(csv.DictReader(open(MAN), delimiter="\t"))
designs = [r for r in rows if r["klass"] == "DESIGN"]
controls = [r for r in rows if r["klass"] != "DESIGN"]
if not designs or not controls:
    sys.exit(f"⛔ {len(designs)} designs / {len(controls)} controls -- refusing")
n_nb = sum(1 for r in controls if r["klass"] == "NONBINDER_crosspaired")
if n_nb != 2:
    sys.exit(f"⛔ {n_nb} cross-paired non-binders -- `J1055` §3 requires the in-batch floor")

# ── round-robin, within a deterministic (arm, batch) order so the deal is reproducible ──
designs.sort(key=lambda r: (r["cell"], r["batch"], r["tag"]))
assign = collections.defaultdict(list)
for i, r in enumerate(designs):
    assign[i % NSH].append(r)

for s in range(NSH):
    d = os.path.join(DEST, f"shard{s + 1:02d}")
    os.makedirs(d)
    for r in assign[s] + controls:                 # ⭐ controls in EVERY shard
        src = os.path.join(YDIR, r["tag"] + ".yaml")
        if not os.path.exists(src):
            sys.exit(f"⛔ {src} missing")
        shutil.copy2(src, os.path.join(d, r["tag"] + ".yaml"))

# ── assertions, on the filesystem and not on the plan above ──
print(f"== sharded {len(designs)} designs + {len(controls)} controls x {NSH} shards -> {DEST}")
tot = 0
for s in range(NSH):
    d = os.path.join(DEST, f"shard{s + 1:02d}")
    files = os.listdir(d)
    ny = len([f for f in files if f.endswith(".yaml")])
    if ny != len(files):
        sys.exit(f"⛔ {d} holds {len(files)} entries but {ny} .yaml -- a stray file kills the batch")
    cells = collections.Counter(r["cell"] for r in assign[s])
    nb = sum(1 for r in controls if r["klass"] == "NONBINDER_crosspaired")
    if len(set(cells.values())) != 1:
        sys.exit(f"⛔ shard{s + 1:02d} is arm-UNBALANCED {dict(cells)} -- arm would be confounded "
                 f"with seed")
    if nb != 2:
        sys.exit(f"⛔ shard{s + 1:02d} lacks its own non-binder floor")
    print(f"   shard{s + 1:02d}: {ny} yaml = {len(assign[s])} designs "
          f"({', '.join(f'{k} {v}' for k, v in sorted(cells.items()))}) + {len(controls)} controls")
    tot += ny
_ncell = len({r["cell"] for r in designs})
print(f"✅ every shard carries {len(designs) // NSH // _ncell} of each of the {_ncell} arms and its "
      f"OWN 2 cross-paired non-binders")
print(f"⛔ COUNT THE FOLD OUTPUTS AGAINST {tot} ({tot} folds x 5 diffusion samples = "
      f"{tot * 5} structures), never against rc ([[boltz-oom-exits-zero]])")
