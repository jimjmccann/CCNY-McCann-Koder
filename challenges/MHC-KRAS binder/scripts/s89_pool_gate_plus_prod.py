#!/usr/bin/env python3
"""s89 -- pool a GATE prefix with its PRODUCTION prefix for the two 2026-11-05 arms, and ASSERT the
shard/design census before anything is scored.

⛔ Extracts to a FRESH tree under `$LAB_SCRATCH`. ⛔ The tarballs under `<VERIFY_DIR>/` are NEVER
   touched -- they are the record ([[rfd3-runs-are-unseeded-and-not-reproducible]]: RFd3 is unseeded,
   so the ARTEFACT is the only record of a run).

⛔⛔ WHY THIS IS NOT `s61`, AND WHY REUSING `s61` WOULD HAVE BEEN WRONG IN THE WORST WAY.
   `s61` pools the MATRIX arm, and two of its deliberate decisions INVERT here:

   ⑴ ⛔⛔ **`s61` EXCLUDES THE GATE PREFIX BY DESIGN** -- correctly, because there
      `gate_20261103_1400` held **shard00 only** while `prod_`/`more5_` ALSO covered shard00, so
      pooling it gave cell 0 nine batches against eight elsewhere and re-created the imbalance a
      re-roll had been bought to remove.
      ⇒ ⛔ **HERE THE GATE *IS* CELL 0 AND MUST BE INCLUDED.** The gate ran shard **00** and
      production ran shards **01-07** -- DISJOINT. Dropping the gate drops a whole cell and leaves
      exactly **280**, a number that looks like the expected answer
      ([[a-clean-empty-report-is-a-failure-mode]]), and reproduces the `J996` §4 failure where
      scoring one prefix silently lost shard03 and reported a clean 200.
   ⑵ ⛔ **`s61` RENUMBERS BATCHES; THIS MUST NOT.** `s61`'s collision was real because both its
      prefixes ran the SAME shards, each numbering batches from 0, so one basename named two
      different designs ([[pod-batch-indices-collide-across-prefixes]]). Here the shard ids are
      disjoint, so basenames are already unique and renumbering would CORRUPT the batch index that
      `n` is counted on. ⇒ ⛔ no renumbering, and the no-collision claim is **ASSERTED, never
      assumed** -- the script refuses if any basename appears in both prefixes.
   ⑶ `s61` also hard-codes 6 shards and REFUSES a prefix that does not cover all of them, so it
      would refuse both of these arms outright.

⇒ ⭐ `s75` then scores the pooled tree UNMODIFIED. Nothing here re-implements a gate or a metric.

⛔⛔ COUNT THE DESIGNS OFF THE FILESYSTEM, NEVER OFF A REPORT LINE. Designs are `.cif.gz` INSIDE
   the per-shard tarballs, so a `find -name '*.pdb'` returns 0 on a GOOD round
   ([[pod-verify-holds-pdbs-inside-tarballs]]).

Usage:
  python3 -I s89_pool_gate_plus_prod.py <arm_root> <dest_dir> <expect_shards> <expect_designs>
    e.g. ... <PROD_DIR>  $LAB_SCRATCH/c5_nterm_pooled_20261105  8  320
"""
import collections
import hashlib
import os
import sys
import tarfile

if len(sys.argv) != 5:
    sys.exit(__doc__)
ROOT, DEST = sys.argv[1], sys.argv[2]
EXP_SHARDS, EXP_DESIGNS = int(sys.argv[3]), int(sys.argv[4])

ROOT = os.path.expanduser(ROOT)
DEST = os.path.expanduser(DEST)
if not os.path.isdir(ROOT):
    sys.exit(f"⛔ arm root not found: {ROOT}")

# ── 1 · census the tarballs, per prefix, WITHOUT extracting ───────────────────────────────────
prefixes = sorted(d for d in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, d)))
if not prefixes:
    sys.exit(f"⛔ no prefix directories under {ROOT}")

per_prefix = {}          # prefix -> {shard: [member names]}
print(f"== census of {ROOT}")
for p in prefixes:
    shards = {}
    pdir = os.path.join(ROOT, p)
    for sd in sorted(os.listdir(pdir)):
        tgz = os.path.join(pdir, sd, f"{sd}.tar.gz")
        if not os.path.exists(tgz):
            continue
        with tarfile.open(tgz) as tf:
            members = [m.name for m in tf.getmembers() if m.name.endswith(".cif.gz")]
        shards[sd] = members
        print(f"   {p:<34} {sd}  {len(members):>4} designs")
    per_prefix[p] = shards

all_shards = {}
for p, shards in per_prefix.items():
    for sd, members in shards.items():
        all_shards.setdefault(sd, []).append(p)

print()
dupes = {sd: ps for sd, ps in all_shards.items() if len(ps) > 1}
if dupes:
    print("⛔⛔ A SHARD APPEARS IN MORE THAN ONE PREFIX -- this is the `s61` situation, and the")
    print("   batch indices collide. ⛔ REFUSING: renumbering is required and this script does not")
    print("   renumber (see its docstring ⑵).")
    for sd, ps in dupes.items():
        print(f"   {sd}: {ps}")
    sys.exit(1)
print(f"✅ shard ids are DISJOINT across the {len(prefixes)} prefixes "
      f"⇒ no batch-index collision, no renumbering needed")

# ── 2 · basename collision, asserted on BYTES not on the claim above ──────────────────────────
basenames = collections.defaultdict(list)
for p, shards in per_prefix.items():
    for sd, members in shards.items():
        for m in members:
            basenames[os.path.basename(m)].append(f"{p}/{sd}")
collide = {b: v for b, v in basenames.items() if len(v) > 1}
if collide:
    print(f"⛔⛔ {len(collide)} BASENAME(S) APPEAR TWICE -- a flat extract would overwrite them.")
    for b, v in list(collide.items())[:5]:
        print(f"   {b}: {v}")
    sys.exit(1)
print(f"✅ all {len(basenames)} design basenames are UNIQUE across prefixes (asserted, not assumed)")

# ── 3 · the census gate, BEFORE extracting ────────────────────────────────────────────────────
n_shards, n_designs = len(all_shards), len(basenames)
print()
print(f"`[M]` shards {n_shards} (expect {EXP_SHARDS})   designs {n_designs} (expect {EXP_DESIGNS})")
if n_shards != EXP_SHARDS or n_designs != EXP_DESIGNS:
    print("⛔⛔ CENSUS MISMATCH -- refusing to extract or score. A scorer pointed at a short census")
    print("   reports a clean number for an incomplete arm (`J996` §4).")
    sys.exit(1)
print("✅ census gate PASSED")

# ── 4 · extract to a fresh tree + mapping ─────────────────────────────────────────────────────
if os.path.exists(DEST) and os.listdir(DEST):
    sys.exit(f"⛔ {DEST} exists and is not empty -- use a fresh dir "
             f"(⛔ never merge two campaigns into one tree)")
os.makedirs(DEST, exist_ok=True)
rows = [("design", "prefix", "shard", "md5")]
n = 0
for p, shards in per_prefix.items():
    pdir = os.path.join(ROOT, p)
    for sd in sorted(shards):
        tgz = os.path.join(pdir, sd, f"{sd}.tar.gz")
        with tarfile.open(tgz) as tf:
            for m in tf.getmembers():
                if not m.name.endswith((".cif.gz", ".json")):
                    continue
                data = tf.extractfile(m).read()
                base = os.path.basename(m.name)
                with open(os.path.join(DEST, base), "wb") as fh:
                    fh.write(data)
                if m.name.endswith(".cif.gz"):
                    rows.append((base, p, sd, hashlib.md5(data).hexdigest()))
                    n += 1
with open(os.path.join(DEST, "00_PREFIX_MAPPING.tsv"), "w") as fh:
    for r in rows:
        fh.write("\t".join(r) + "\n")

on_disk = len([f for f in os.listdir(DEST) if f.endswith(".cif.gz")])
print()
print(f"extracted {n} designs -> {DEST}")
print(f"`[M]` COUNTED ON THE FILESYSTEM AFTER EXTRACTION: {on_disk} `.cif.gz`")
if on_disk != EXP_DESIGNS:
    print(f"⛔⛔ {on_disk} on disk against {EXP_DESIGNS} expected -- a flat extract LOST files.")
    sys.exit(1)
print("✅ post-extraction count matches the census ⇒ nothing was overwritten")
print(f"wrote 00_PREFIX_MAPPING.tsv ({len(rows)-1} rows), originals under {ROOT} UNTOUCHED")
print()
print("⇒ next: score with `s75` (UNMODIFIED) over "
      f"'{DEST}/*.cif.gz'")
