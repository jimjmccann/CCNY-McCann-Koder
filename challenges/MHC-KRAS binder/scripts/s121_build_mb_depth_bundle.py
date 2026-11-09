#!/usr/bin/env python3
"""s121 -- build the MB-depth RFd3 bundle (J1058): 2 shards, both `MB_nomotif_broad`.

Usage: python3 -I s121_build_mb_depth_bundle.py <outdir>

⛔⛔ THIS IS A DEPTH ROUND IN ONE CELL, NOT A NEW EXPERIMENT. Nothing varies between the two
shards and nothing varies from the 2x2's winning arm. Per [[never-inherit-a-default-record-the-
resolved-config]]: CUT ARMS, NEVER n.

⭐⭐ WHY THE SPEC IS **COPIED**, NOT REGENERATED. `s112` built the 2x2's four specs from literals,
and re-running that construction here would re-derive `MB_nomotif_broad` from the same literals --
which is *probably* byte-identical and cannot be *proven* so without comparing hashes anyway. So
this script copies `rfd3_c5_designability_2x2_v1/shards/arm_shard03.json` byte-for-byte and asserts
its sha1 against the value `ARMS.tsv` recorded for the arm that actually ran. ⇒ the new round's
specs are the SAME BYTES the measured cell ran, not a reconstruction that resembles them.
⛔ The atom-existence assertions are still run independently against the input PDB, because a
sha1 match proves the spec is unchanged, NOT that the input still satisfies it
[[grep-before-citing-is-not-sufficient]].

⛔ RFd3 is UNSEEDED [[rfd3-runs-are-unseeded-and-not-reproducible]] ⇒ two byte-identical shards are
two INDEPENDENT draws. That is the point of the round; the artefact is the record.
⛔ Batch indices restart per shard and COLLIDE [[pod-batch-indices-collide-across-prefixes]] ⇒
every downstream extract goes to its OWN root, md5 BEFORE extracting.
⛔ RFd3 output here is `.cif.gz`, NOT `.pdb` -- a `find -name '*.pdb'` returns 0 on a GOOD round.
"""
import hashlib
import json
import os
import shutil
import sys

SRC = "/PATH/TO/bundles"
TWOBYTWO = os.path.join(SRC, "rfd3_c5_designability_2x2_v1")
MATRIX = os.path.join(SRC, "rfd3_pmhc_matrix_v1")

# ---- the three values this round inherits from the MEASURED cell, by hash ----
RUNNER_SHA1 = "623b3c69924b4e0afb245330d5eb2b6f39df10de"
MB_SPEC_SHA1 = "23d9c25cd0689ebe2e68a9d280f9d807ba9e3787"   # ARMS.tsv row 03
ARM_NAME = "MB_nomotif_broad"
N_SHARDS = 2
INPUT_PDB = "arm3_hotspot_only_input.pdb"

OUT = sys.argv[1]


def sha1(p):
    h = hashlib.sha1()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def pdb_atoms(p):
    """{(chain, resnum): set(atom names)} over ATOM records only."""
    d = {}
    with open(p) as fh:
        for ln in fh:
            if ln.startswith("ATOM"):
                d.setdefault((ln[21], int(ln[22:26])), set()).add(ln[12:16].strip())
    return d


os.makedirs(os.path.join(OUT, "shards"), exist_ok=True)

# ---- 1. the source spec, verified by hash BEFORE it is copied anywhere ----
src_spec = os.path.join(TWOBYTWO, "shards", "arm_shard03.json")
got = sha1(src_spec)
if got != MB_SPEC_SHA1:
    sys.exit("REFUSING: %s sha1 %s != the ARMS.tsv value %s for arm %s. ⛔ The cell that was "
             "MEASURED is not the cell on disk -- do not build a depth round on it."
             % (src_spec, got, MB_SPEC_SHA1, ARM_NAME))
spec = json.load(open(src_spec))
if list(spec) != [ARM_NAME]:
    sys.exit("REFUSING: %s holds arms %s, expected exactly [%r]" % (src_spec, list(spec), ARM_NAME))
print("✅ source spec sha1 %s matches ARMS.tsv row 03 (%s)" % (got, ARM_NAME))

body = spec[ARM_NAME]
if body["input"] != INPUT_PDB:
    sys.exit("REFUSING: spec names input %r, expected %r" % (body["input"], INPUT_PDB))
if body.get("select_fixed_atoms", {}).get("E101") or any(
        k.startswith("E") for k in body.get("select_fixed_atoms", {})):
    sys.exit("REFUSING: the MB cell must hold NO E-chain motif, but select_fixed_atoms has one. "
             "⛔ The held motif is refuted on two independent axes; this round cannot rescue it.")
print("✅ no held motif in the spec (select_fixed_atoms: %s)" % sorted(body["select_fixed_atoms"]))

# ---- 2. the runner and the MANDATORY engine patch, byte-for-byte, sha1-verified ----
for fn in ("prod_rfd3_multi.sh", "patch_held_motif_guard.py"):
    shutil.copy2(os.path.join(TWOBYTWO, fn), os.path.join(OUT, fn))
got = sha1(os.path.join(OUT, "prod_rfd3_multi.sh"))
if got != RUNNER_SHA1:
    sys.exit("REFUSING: runner sha1 %s != expected %s. ⛔ RFD3_SCRIPT_SHA1 must match or the pod "
             "burns an 1800 s sleep at full billing rate." % (got, RUNNER_SHA1))
print("✅ runner sha1 %s matches the canonical value" % got)
print("✅ patch_held_motif_guard.py shipped unchanged (mandatory; the runner dies if it refuses)")

# ---- 3. the input PDB, and the assertions a sha1 match does NOT cover ----
shutil.copy2(os.path.join(MATRIX, INPUT_PDB), os.path.join(OUT, INPUT_PDB))
at = pdb_atoms(os.path.join(OUT, INPUT_PDB))
chains = sorted({c for c, _ in at})
if chains != ["A", "B", "C"]:
    sys.exit("REFUSING: %s has chains %s, expected ['A','B','C']" % (INPUT_PDB, chains))
pep = sorted(n for c, n in at if c == "C")
if pep != list(range(1, 11)):
    sys.exit("REFUSING: %s chain C is %s, expected residues 1-10" % (INPUT_PDB, pep))
print("  %-34s %5d residues  chains %s  chain C = the 10-mer peptide"
      % (INPUT_PDB, len(at), chains))

# ⛔⛔ `select_hotspots` is the ONLY placement signal when infer_ori_strategy is "hotspots" -- a
#    silently-dropped atom changes the ask, and a missing atom is a pod-side failure at full rate.
hs = body["select_hotspots"]
for res, atoms in hs.items():
    key = (res[0], int(res[1:]))
    if key not in at:
        sys.exit("REFUSING: %s has no residue %s" % (INPUT_PDB, res))
    for a in atoms.split(","):
        if a not in at[key]:
            sys.exit("REFUSING: %s %s has no atom %s -- present: %s"
                     % (INPUT_PDB, res, a, sorted(at[key])))
nat = sum(len(v.split(",")) for v in hs.values())
print("✅ all %d BROAD hotspot atoms exist in this round's own input PDB (%s)"
      % (nat, ";".join("%s:%s" % kv for kv in sorted(hs.items()))))

# ⛔ every fixed target range must exist too, or the pod holds atoms that are not there.
for rng in body["select_fixed_atoms"]:
    ch, span = rng[0], rng[1:]
    lo, hi = (int(x) for x in span.split("-"))
    missing = [n for n in range(lo, hi + 1) if (ch, n) not in at]
    if missing:
        sys.exit("REFUSING: select_fixed_atoms %s names %d residues absent from %s (first: %s)"
                 % (rng, len(missing), INPUT_PDB, missing[:5]))
print("✅ every select_fixed_atoms range is fully present in the input")

# ---- 4. N identical shards, each asserted byte-identical to the measured cell ----
rows = []
for sh in range(N_SHARDS):
    p = os.path.join(OUT, "shards", "arm_shard%02d.json" % sh)
    shutil.copy2(src_spec, p)
    s1 = sha1(p)
    if s1 != MB_SPEC_SHA1:
        sys.exit("REFUSING: written shard %02d sha1 %s != %s" % (sh, s1, MB_SPEC_SHA1))
    with open(os.path.join(OUT, "shards", "arm_shard%02d.flags" % sh), "w") as fh:
        fh.write("inference_sampler.n_recycle=3\n")
    rows.append((sh, s1))

with open(os.path.join(OUT, "ARMS.tsv"), "w") as fh:
    fh.write("shard\tarm\tmotif\thotspots\tinput\tcontig\tspec_sha1\n")
    for sh, s1 in rows:
        fh.write("%02d\t%s\t%s\t%s\t%s\t%s\t%s\n"
                 % (sh, ARM_NAME, "NONE",
                    ";".join("%s:%s" % kv for kv in sorted(hs.items())),
                    INPUT_PDB, body["contig"], s1))

print()
print("=== %d shards written to %s ===" % (len(rows), OUT))
for sh, s1 in rows:
    print("  shard %02d  %-20s  spec sha1 %s" % (sh, ARM_NAME, s1))
print()
print("⭐ both shards are the SAME BYTES as the 2x2 arm that was measured; RFd3 is unseeded, so")
print("   they are two INDEPENDENT draws of 20 batches each = 40 new MB batches.")
print("⛔ n is BATCHES, not designs: 20 nbatch x diffusion_batch_size 4 = 80 designs/shard.")
print("⛔ NOTHING IS LAUNCHED OR UPLOADED BY THIS SCRIPT.")
