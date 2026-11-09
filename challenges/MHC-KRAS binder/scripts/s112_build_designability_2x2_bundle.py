#!/usr/bin/env python3
"""s112 -- build the RFd3 bundle for the option-③ designability 2x2 (J1055).

Usage: python3 -I s112_build_designability_2x2_bundle.py <outdir>

⛔⛔ THE 2x2 IS held-motif x hotspot-breadth, AND NOTHING ELSE VARIES. Pre-registered in
`J1055` §4. Length is held at L65_85 in all four arms; `n_recycle`, `dialect`,
`infer_ori_strategy`, `is_non_loopy` are identical; guidepost radius and design symmetry are
held at their single absent/default value because `[M]` neither setting exists anywhere in the
previous config.

⛔⛔ WHY THE NO-MOTIF ARMS ARE **NOT** A VERBATIM COPY OF THE MATRIX `arm3` SPEC, EVEN THOUGH
   `arm3` IS THE PROVEN-GOOD SHAPE. `rfd3_pmhc_matrix_v1/shards/arm_shard05.json` carries
   ⛔ NO `select_fixed_atoms` AT ALL -- not even the target chains -- while every cdr3b cell sets
   `A1-276/B1-99/C1-10: ALL`. Copying `arm3` verbatim would therefore change TWO things at once
   (the held E motif AND whether the target is explicitly held) and the 2x2 would no longer
   isolate one variable [[test-both-branches-of-a-guard]].
   ⇒ ⭐ M/MB = R/B **minus only the `E...` lines**, target fixing retained. Both halves of that
   combination are individually proven (target-fixed is the cdr3b shape; motif-free is arm3's),
   ⛔ but the COMBINATION has not run before.
   ⚠️ IF A NO-MOTIF SHARD FAILS THE POD'S OWN PREFLIGHT, the known-good fallback is `arm3`'s exact
   spec (drop `select_fixed_atoms` entirely) -- ⛔ and that fallback REINTRODUCES the confound, so
   it must be recorded as a deviation, never taken silently.

⛔⛔ THE BROAD HOTSPOT SET IS RE-DERIVED AGAINST `t142`, ⛔ NOT COPIED FROM `arm3`. `arm3` chose from
   the `t44`/`t92` maps, where `[R]` `J1047`:124-125 measured that **28 of 60 entries were
   artefacts**. `arm3`'s `C6:CB`, `C6:OD1` and `C7:CA` are ABSENT from `t142`'s usable-32 and are
   DROPPED here. Every atom below is asserted present in `t142`'s usable list AND in the input PDB
   [[grep-before-citing-is-not-sufficient]].
"""
import hashlib
import json
import os
import shutil
import sys

SRC = "/PATH/TO/bundles"
CDR3B = os.path.join(SRC, "rfd3_cdr3b_v1")
MATRIX = os.path.join(SRC, "rfd3_pmhc_matrix_v1")
RUNNER_SHA1 = "623b3c69924b4e0afb245330d5eb2b6f39df10de"

OUT = sys.argv[1]

# ---- the two hotspot sets ----
NARROW = {"C8": "CG1,CG2"}                 # the previous campaign's ENTIRE steering signal
BROAD = {"C6": "O", "C7": "O", "C8": "CB,CG1,CG2", "C9": "N"}

# ⛔ every BROAD atom, with its `t142` slot count, so the choice is auditable from this file alone.
#   `[M]` read off tables/t142_site_map_donor_and_rim_fixed_20261108.txt:252-284.
T142_USABLE_PEPTIDE = {
    ("C7", "O"): 1454, ("C4", "N"): 405, ("C9", "N"): 350, ("C6", "O"): 312,
    ("C8", "CB"): 298, ("C8", "CG2"): 277, ("C8", "CG1"): 203, ("C5", "N"): 145,
    ("C5", "CB"): 130,
}

TARGET_FIXED = {"A1-276": "ALL", "B1-99": "ALL", "C1-10": "ALL"}
MOTIF_FIXED = {"E96-103": "N,CA,C,O", "E101": "CG,OD1,ND2"}

ARMS = [
    # (shard, name, input pdb, contig, has_motif, hotspots)
    (0, "R_motif_narrow", "cdr3b_trim_pmhc_input.pdb",
     "28-38,E96-103,29-39,/0,A1-276,/0,B1-99,/0,C1-10", True, NARROW),
    (1, "B_motif_broad", "cdr3b_trim_pmhc_input.pdb",
     "28-38,E96-103,29-39,/0,A1-276,/0,B1-99,/0,C1-10", True, BROAD),
    (2, "M_nomotif_narrow", "arm3_hotspot_only_input.pdb",
     "65-85,/0,A1-276,/0,B1-99,/0,C1-10", False, NARROW),
    (3, "MB_nomotif_broad", "arm3_hotspot_only_input.pdb",
     "65-85,/0,A1-276,/0,B1-99,/0,C1-10", False, BROAD),
]


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


# ---- assert the BROAD set is a subset of t142's usable peptide anchors ----
for res, atoms in BROAD.items():
    for a in atoms.split(","):
        if (res, a) not in T142_USABLE_PEPTIDE:
            sys.exit("REFUSING: %s:%s is NOT in t142's usable peptide anchors -- the broad set must "
                     "be re-derived against the CORRECTED map, never a superseded one" % (res, a))
print("✅ broad set is a subset of t142's usable peptide anchors (%d atoms, slots %s)"
      % (sum(len(v.split(',')) for v in BROAD.values()),
         [T142_USABLE_PEPTIDE[(r, a)] for r, v in BROAD.items() for a in v.split(",")]))

os.makedirs(os.path.join(OUT, "shards"), exist_ok=True)

# ---- copy the runner and the MANDATORY engine patch, byte-for-byte, sha1-verified ----
for fn, src in (("prod_rfd3_multi.sh", CDR3B), ("patch_held_motif_guard.py", CDR3B)):
    shutil.copy2(os.path.join(src, fn), os.path.join(OUT, fn))
got = sha1(os.path.join(OUT, "prod_rfd3_multi.sh"))
if got != RUNNER_SHA1:
    sys.exit("REFUSING: runner sha1 %s != expected %s. ⛔ RFD3_SCRIPT_SHA1 must match or the pod "
             "burns an 1800 s sleep at full billing rate." % (got, RUNNER_SHA1))
print("✅ runner sha1 %s matches the canonical value" % got)
# ⛔ the patch is NOT optional and the runner dies if it refuses; M/MB hold no motif so it has
#   nothing to guard, but it is shipped UNCHANGED because the runner's own guard is not ours to
#   branch on.
print("✅ patch_held_motif_guard.py shipped unchanged (mandatory; runner dies if it refuses)")

# ---- copy both input PDBs and validate them ----
inputs = {}
for fn, src in (("cdr3b_trim_pmhc_input.pdb", CDR3B), ("arm3_hotspot_only_input.pdb", MATRIX)):
    shutil.copy2(os.path.join(src, fn), os.path.join(OUT, fn))
    at = pdb_atoms(os.path.join(OUT, fn))
    chains = sorted({c for c, _ in at})
    inputs[fn] = at
    print("  %-34s %5d residues  chains %s" % (fn, len(at), chains))

want = {"cdr3b_trim_pmhc_input.pdb": ["A", "B", "C", "E"],
        "arm3_hotspot_only_input.pdb": ["A", "B", "C"]}
for fn, w in want.items():
    got_ch = sorted({c for c, _ in inputs[fn]})
    if got_ch != w:
        sys.exit("REFUSING: %s has chains %s, expected %s" % (fn, got_ch, w))
# ⛔ the peptide must be the 10-mer in BOTH inputs or the two halves of the 2x2 are not the same
#   target [[always-name-the-pdb-chain]].
for fn in inputs:
    pep = sorted(n for c, n in inputs[fn] if c == "C")
    if pep != list(range(1, 11)):
        sys.exit("REFUSING: %s chain C is %s, expected residues 1-10" % (fn, pep))
print("✅ both inputs carry chain C residues 1-10 (the same 10-mer peptide)")

specs = []
for sh, name, pdb, contig, has_motif, hs in ARMS:
    at = inputs[pdb]
    # ⛔⛔ ASSERT EVERY HOTSPOT ATOM EXISTS IN THIS ARM'S OWN INPUT PDB. A missing atom is a pod-side
    #   failure at full billing rate, and `select_hotspots` is the ONLY placement signal when
    #   `infer_ori_strategy` is "hotspots" -- a silently-dropped atom changes the ask.
    for res, atoms in hs.items():
        rn = int(res[1:])
        key = (res[0], rn)
        if key not in at:
            sys.exit("REFUSING: %s has no residue %s (arm %s)" % (pdb, res, name))
        for a in atoms.split(","):
            if a not in at[key]:
                sys.exit("REFUSING: %s %s has no atom %s (arm %s) -- present: %s"
                         % (pdb, res, a, name, sorted(at[key])))
    fixed = dict(TARGET_FIXED)
    if has_motif:
        fixed.update(MOTIF_FIXED)
        for res in ("E101",):
            if (res[0], int(res[1:])) not in at:
                sys.exit("REFUSING: %s has no %s but arm %s pins it" % (pdb, res, name))
    spec = {name: {"dialect": 2,
                   "infer_ori_strategy": "hotspots",
                   "input": pdb,
                   "contig": contig,
                   "select_fixed_atoms": fixed,
                   "select_hotspots": hs,
                   "is_non_loopy": True}}
    p = os.path.join(OUT, "shards", "arm_shard%02d.json" % sh)
    with open(p, "w") as fh:
        json.dump(spec, fh, indent=2)
        fh.write("\n")
    # the per-shard Hydra sidecar, identical on every shard exactly as the previous round
    with open(os.path.join(OUT, "shards", "arm_shard%02d.flags" % sh), "w") as fh:
        fh.write("inference_sampler.n_recycle=3\n")
    specs.append((sh, name, pdb, contig, has_motif, hs, sha1(p)))

with open(os.path.join(OUT, "ARMS.tsv"), "w") as fh:
    fh.write("shard\tarm\tmotif\thotspots\tinput\tcontig\tspec_sha1\n")
    for sh, name, pdb, contig, has_motif, hs, s1 in specs:
        fh.write("%02d\t%s\t%s\t%s\t%s\t%s\t%s\n"
                 % (sh, name, "held_E96-103+E101" if has_motif else "NONE",
                    ";".join("%s:%s" % kv for kv in sorted(hs.items())), pdb, contig, s1))

print()
print("=== %d shards written to %s ===" % (len(specs), OUT))
print("  %-5s %-20s %-18s %s" % ("shard", "arm", "motif", "hotspots"))
for sh, name, pdb, contig, has_motif, hs, s1 in specs:
    print("  %-5s %-20s %-18s %s"
          % ("%02d" % sh, name, "held" if has_motif else "NONE",
             ";".join("%s:%s" % kv for kv in sorted(hs.items()))))
print()
print("⛔ the 2x2 is held-motif x hotspot-breadth. Length L65_85 in all four. n_recycle=3 in all four.")
print("⛔ NOTHING IS LAUNCHED OR UPLOADED BY THIS SCRIPT.")
