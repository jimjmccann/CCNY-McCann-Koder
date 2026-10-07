#!/usr/bin/env python3
"""Pack an ecto621 fold campaign into a bundle, in the layout the batch runner expects.

Bundle layout (copied from an internal bundle builder that is not part of this deposit):
  bundle_<arm>/msa/target.a3m      <- chain B MSA; YAMLs must reference /workspace/msa/target.a3m
  bundle_<arm>/fold_shard_pod.sh   <- the fold runner (not part of this deposit)
  bundle_<arm>/shardNN.tar.gz      <- one tar per worker
  bundle_<arm>/drawsNN             <- draw count per shard

⛔ ONE SHARD PER WORKER, ONE OUTPUT PREFIX PER WORKER. Two workers sharing one output prefix means
the first to finish writes a complete-looking result set and both are then treated as done.
⛔ The YAMLs are written by build_ecto621_yamls.py with --msa /workspace/msa/target.a3m; this script
ASSERTS that rather than trusting it, because a missing MSA makes boltz SKIP the design and exit 0.
⛔ Glycosylated inputs OOM'd a 16 GB GPU at ~1,260 tokens and still exited rc=0 with
"ran out of memory, skipping batch" -- a silent no-output failure. Check the per-shard prediction
count against the input count on EVERY glyco shard; do not read rc=0 as success.

2026-10-02.
"""
import argparse, glob, hashlib, os, shutil, subprocess, sys

REPO = "/PATH/TO/CHECKOUT"
RUNNER = f"{REPO}/<worker-runner-not-distributed>"

ap = argparse.ArgumentParser()
ap.add_argument("--yaml-root", required=True, help="dir containing shardNN/ of YAMLs")
ap.add_argument("--a3m", required=True)
ap.add_argument("--arm", required=True, help="e.g. egfrEctoA (aglyco) / egfrEctoG (glyco)")
ap.add_argument("--outdir", required=True, help="bundle parent; must not already hold this bundle")
ap.add_argument("--draws", type=int, default=3)
a = ap.parse_args()

for p in (a.yaml_root, a.a3m, RUNNER):
    if not os.path.exists(p):
        sys.exit(f"⛔ missing {p}")

shards = sorted(os.path.basename(d) for d in glob.glob(os.path.join(a.yaml_root, "shard*")))
if not shards:
    sys.exit(f"⛔ no shard*/ under {a.yaml_root}")
yamls = {s: sorted(glob.glob(os.path.join(a.yaml_root, s, "*.yaml"))) for s in shards}
total = sum(len(v) for v in yamls.values())
if total == 0:
    sys.exit("⛔ zero YAMLs")

# ---- assert the worker-absolute MSA path in EVERY yaml -------------------------------------
bad = [y for s in shards for y in yamls[s] if "/workspace/msa/target.a3m" not in open(y).read()]
if bad:
    sys.exit(f"⛔ {len(bad)} YAML(s) do not reference /workspace/msa/target.a3m, e.g. "
             f"{os.path.basename(bad[0])} -- a missing MSA makes boltz SKIP and exit 0")
print(f"✅ all {total:,} YAMLs carry the worker-absolute MSA path")

# ---- assert the a3m query matches the chain-B sequence in the YAMLs ---------------------
def a3m_query(path):
    for ln in open(path):
        ln = ln.rstrip("\n")
        if ln.startswith("#") or ln.startswith(">") or not ln:
            continue
        return ln
    sys.exit(f"⛔ no query sequence in {path}")

q = a3m_query(a.a3m)
probe = yamls[shards[0]][0]
tgt = None
lines = open(probe).read().splitlines()
for i, ln in enumerate(lines):
    if ln.strip() == "id: B":
        for j in range(i, min(i + 4, len(lines))):
            if lines[j].strip().startswith("sequence:"):
                tgt = lines[j].split("sequence:", 1)[1].strip()
if tgt is None:
    sys.exit(f"⛔ could not read chain B sequence from {probe}")
if q != tgt:
    sys.exit(f"⛔⛔ MSA/TARGET MISMATCH -- refusing to bundle.\n"
             f"   a3m query  {len(q)} aa: {q[:50]}...\n"
             f"   YAML chain B {len(tgt)} aa: {tgt[:50]}...\n"
             f"   boltz would fold with an MSA describing a DIFFERENT protein.")
print(f"✅ a3m query matches the YAML chain-B target exactly ({len(q)} aa)")

bundle = os.path.join(a.outdir, f"bundle_{a.arm}")
if os.path.exists(bundle):
    sys.exit(f"⛔ {bundle} exists -- fresh path per experiment; move it by hand")
os.makedirs(os.path.join(bundle, "msa"))
shutil.copy(a.a3m, os.path.join(bundle, "msa", "target.a3m"))
if os.path.exists(RUNNER):
    shutil.copy(RUNNER, bundle)
else:
    print(f"⚠️ no runner at {RUNNER}: the worker-side runner is not part of this deposit.")
    print("   Put your own runner in the bundle before submitting the shards.")
for s in shards:
    nn = s.replace("shard", "")
    with open(os.path.join(bundle, f"draws{nn}"), "w") as fh:
        fh.write(f"{a.draws}\n")
    subprocess.run(["tar", "-C", a.yaml_root, "-czf",
                    os.path.join(bundle, f"{s}.tar.gz"), s], check=True)

runner_path = os.path.join(bundle, "fold_shard_pod.sh")
runner_txt = open(runner_path).read() if os.path.exists(runner_path) else ""
print(f"\n⭐ BUNDLE READY: {bundle}")
print(f"   {len(shards)} shard(s) = {len(shards)} worker(s), {total:,} YAMLs, {a.draws} draws each")
print("   per-shard: " + ", ".join(f"{s}={len(yamls[s])}" for s in shards))
if runner_txt:
    print(f"   runner sha1 {hashlib.sha1(runner_txt.encode()).hexdigest()}")
    print(f"   '--no_kernels' in runner: {'--no_kernels' in runner_txt}   "
          f"(⚠️ PRESENT is the intended setting -- removing it has been observed to break the fold)")
else:
    print("   ⚠️ no runner in the bundle, so the sha1 / flag check is skipped.")
print(f"\nUPLOAD: copy {bundle} to wherever the workers read bundles from.")
print("SUBMIT: one job per shard on your batch system, each writing to its own output prefix.")
print(f"        shards: {','.join(s.replace('shard','') for s in shards)}")
