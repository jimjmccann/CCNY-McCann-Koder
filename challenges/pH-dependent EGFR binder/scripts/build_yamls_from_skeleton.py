#!/usr/bin/env python3
"""Write Boltz-2 YAMLs for new BINDER sequences by reusing a LANDED campaign's YAML verbatim.

⭐ WHY THIS AND NOT ANOTHER BUILDER. The egfrP14 YAMLs (ecto621 chain B + 10 NAG + the N-glycan
bond constraints + the FORCED template block) were written inline and the script was not kept.
Rewriting that construct from scratch is the expensive way to make a silent mistake: the glycan
bond list, `template_id: Axp`, `force: true` + `threshold: 2.0` and the worker-absolute MSA path all
had to be discovered by measurement (2026-10-03). So this script takes a
REAL P14 yaml as the skeleton and changes exactly ONE line: chain A's sequence. Everything that
made that campaign fold against a correctly-built target is carried over byte-for-byte.

⛔ IT ASSERTS the skeleton instead of trusting it: one chain-A sequence line, a worker-absolute
chain-B MSA (a relative path makes boltz SKIP the design and exit 0), and a `templates:` block
(the runner turns --use_potentials on by grepping for it, and without potentials `force: true` is
silently ignored -- measured tgt_fit 17.5 A vs 0.69 A).

2026-10-03.
"""
import argparse, os, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("--skeleton", required=True, help="a landed P14 yaml to copy")
ap.add_argument("--seqs", required=True, help="TSV: name<TAB>binder_sequence")
ap.add_argument("--outdir", required=True, help="yaml root; shardNN/ written under it")
ap.add_argument("--shards", type=int, required=True)
ap.add_argument("--require-template", action="store_true", default=True)
a = ap.parse_args()

AA = set("ACDEFGHIKLMNPQRSTVWY")
sk = open(a.skeleton).read()

# ---- assert the skeleton is the construct we think it is --------------------------------
if "msa: /workspace/msa/target.a3m" not in sk:
    sys.exit("⛔ skeleton has no worker-absolute chain-B MSA -- boltz would SKIP and exit 0")
if a.require_template and not re.search(r"^templates:", sk, re.M):
    sys.exit("⛔ skeleton carries no templates: block -- the runner would fold without "
             "--use_potentials and the forced template would be silently ignored")
if "force: true" not in sk:
    sys.exit("⛔ skeleton template is not forced (force: true) -- a soft template does nothing")

lines = sk.splitlines()
ida = [i for i, ln in enumerate(lines) if ln.strip() == "id: A"]
if len(ida) != 1:
    sys.exit(f"⛔ expected exactly one 'id: A' in the skeleton, found {len(ida)}")
seq_i = None
for j in range(ida[0], min(ida[0] + 4, len(lines))):
    if lines[j].strip().startswith("sequence:"):
        seq_i = j
        break
if seq_i is None:
    sys.exit("⛔ could not find chain A's sequence line")
indent = lines[seq_i][:len(lines[seq_i]) - len(lines[seq_i].lstrip())]
print(f"✅ skeleton OK: chain A sequence at line {seq_i + 1} "
      f"({len(lines[seq_i].split('sequence:', 1)[1].strip())} aa), forced template, absolute MSA")

rows = []
for ln in open(a.seqs):
    ln = ln.rstrip("\n")
    if not ln:
        continue
    name, seq = ln.split("\t")
    seq = seq.strip().upper()
    if set(seq) - AA:
        sys.exit(f"⛔ {name}: non-standard residues {sorted(set(seq) - AA)}")
    rows.append((name, seq))
if len({n for n, _ in rows}) != len(rows):
    sys.exit("⛔ duplicate names in --seqs")
print(f"✅ {len(rows):,} sequences, {len({s for _, s in rows}):,} distinct, "
       f"len {min(len(s) for _, s in rows)}-{max(len(s) for _, s in rows)}")

if os.path.exists(a.outdir):
    sys.exit(f"⛔ {a.outdir} exists -- fresh path per experiment")

# round-robin, so a shard is never one backbone's worth of a single length
per = [[] for _ in range(a.shards)]
for k, r in enumerate(rows):
    per[k % a.shards].append(r)

for s, chunk in enumerate(per):
    d = os.path.join(a.outdir, f"shard{s:02d}")
    os.makedirs(d)
    for name, seq in chunk:
        out = list(lines)
        out[seq_i] = f"{indent}sequence: {seq}"
        with open(os.path.join(d, f"{name}.yaml"), "w") as fh:
            fh.write("\n".join(out) + "\n")

print(f"⭐ {a.outdir}: {a.shards} shard(s), " +
      ", ".join(f"shard{s:02d}={len(c)}" for s, c in enumerate(per)))
