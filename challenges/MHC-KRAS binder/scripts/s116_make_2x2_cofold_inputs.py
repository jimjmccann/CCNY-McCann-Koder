#!/usr/bin/env python3
"""s116 -- emit the Boltz co-fold inputs for the designability-2x2 refold cohort, plus the controls.

    RUN WITH:  python3 -I   (stdlib only; no network, no spend)

⭐⭐ WHY THIS IS NEITHER `s59` NOR `s71`, each checked against its OWN guards rather than remembered
   ([[check-a-bundle-against-its-runners-own-guards]]):
     - ⛔ `s59` HARD-EXITS on this cohort twice: `s59`:114 demands exactly **128** gated designs
       (we have 260) and `s59`:124 demands exactly **6** cells (we have **4**). It also re-does its
       own pre-registered top-N-per-cell ranking on `cont_pep`, which would ⛔⛔ SELECT ON THE
       ENDPOINT -- the cohort for this round was already chosen by `s115` at RANDOM, one per batch,
       precisely so the refold comparison is not selected on the thing it measures.
     - ⛔ `s71` is driven by a `--picks` TSV expressing SEVERAL sequences per backbone (22 from 6).
       Our shape is ONE sequence per backbone over 80 backbones, which it cannot express.
   ⇒ a session that ran either verbatim would meet `FATAL` **after** the pod is up and billing.

⭐ WHAT IS REUSED, IMPORTED NOT COPIED, because comparability depends on it being the same code:
   `s71.lift_targets` (targets lifted from `s35`'s SOURCE by `ast`, ⛔ never from a model output
   [[never-source-target-sequence-from-model-output]]), `s71.pdb_chains`, `s71.write_yaml`
   (`msa: empty` on EVERY chain -- `J986` §3 rule 1, `[M]` 2.3x better class separation).
   ⭐ `s71` is importable: it guards its entry with `if __name__ == "__main__"`, so importing runs
   nothing. `s59` does NOT (it executes at module level and would hard-exit), so the two pieces
   needed from it -- `best_sequence` and the CONTROLS block -- are copied and ASSERTED against
   `s59`'s source at run time.

⛔ THE CHAIN REMAP IS NOT COSMETIC [[always-name-the-pdb-chain]]:
     design PDB / MPNN FASTA order:  A=binder  B=HLA  C=b2m  D=peptide
     fold YAML order (the study's):  A=HLA     B=b2m  C=peptide  D=binder
   ⇒ binder<->peptide is `pair_chains_iptm(3,2)` in the fold output, and ⛔ `J1054` §2 measured that
   this pair ipTM is CONFOUNDED BY THE SHARED ANCHOR CHAIN: C and D both hang off A, so a binder
   glued anywhere on HLA scores high without touching the peptide (Spearman +0.681 to HLA contacts
   against +0.351 to peptide). ⇒ it is ⛔ never read alone; the endpoint is the CONJUNCTION.

⛔⛔ THE FOUR CONTROLS ARE REBUILT IN THIS BATCH (`J986` §3 rule 5 -- the floor MOVES with the arm
   and is never read off an older table). Two POSITIVE crystal complexes and ⭐⭐ TWO CROSS-PAIRED
   NON-BINDERS: a real published binder for a DIFFERENT pMHC, folded against OUR target. That is
   `J1055` §3's pre-registered requirement, and its absence is exactly why `J1054`'s cohort had no
   internal null and `s60` could not be used on it.

⚠️ `OUR_B2M` keeps the initiator Met in the YAML although the design PDBs' b2m starts one residue
   later (`OUR_B2M[1:]`). That asymmetry is `s59`'s and is KEPT DELIBERATELY: the measured G9 floor
   was built that way and changing it would break comparability
   ([[comparability-to-a-broken-baseline-is-not-a-virtue]] cuts the other way only when the baseline
   is BROKEN; here it is merely asymmetric, and consistently so).

⛔ WHAT THIS FOLD CAN AND CANNOT SAY:
   - it ⛔ CANNOT rank peptide-IDENTITY selectivity (rung ⑶) -- `J989` §1: no bound-state interface
     metric from any model ranks WT-vs-G12D here, and the crystals invert the same way;
   - the designed hydrogen bond is a HYPOTHESIS this fold TESTS, ⛔ not delivered geometry
     [[rfd3-guidepost-functional-group-is-not-bonded-on]];
   - ⭐ UNLIKE 5B, this cohort was NOT selected on `cont_pep`, so a high folded peptide-contact
     count IS informative here rather than circular -- that is the whole point of `s115`'s random
     pick;
   - ⇒ what it tests is whether the SEQUENCE rebuilds the designed POSE
     [[design-vs-prediction-selfconsistency]].

⛔ COUNT THE OUTPUTS, NEVER THE EXIT CODE: a Boltz CUDA OOM exits **rc=0** having written nothing
   [[boltz-oom-exits-zero]]. The manifest row count printed at the end is what the fold outputs are
   checked against.

Usage:
  python3 -I s116_make_2x2_cofold_inputs.py --cohort <s115 cohort.txt> --t149 <t149.tsv>
          --mpnn <s53 outdir> --s35 <s35 py> --out <fresh dir>
"""
import argparse
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
S71 = os.path.join(HERE, "s71_make_5b_cofold_inputs.py")
S59 = os.path.join(HERE, "s59_make_r4_cofold_inputs.py")


def load_s71():
    """s71 guards its entry point, so importing it executes nothing."""
    spec = importlib.util.spec_from_file_location("s71_mod", S71)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for need in ("lift_targets", "pdb_chains", "write_yaml", "TARGET_KEYS"):
        if not hasattr(mod, need):
            sys.exit(f"⛔ s71 lacks `{need}` -- it has changed shape; re-audit before use")
    return mod


s71 = load_s71()

# ⛔ copied from `s59`:79-107 and ASSERTED against its source. Record 0 of an MPNN fasta is the
#   INPUT sequence and carries no `id=`; taking it would fold the poly-Ala backbone.
_S59_SRC = open(S59).read()
for _needle in ('def best_sequence(fa):', '", id=" not in hdr', 'overall_confidence='):
    if _needle not in _S59_SRC:
        sys.exit(f"⛔ REFUSING: `{_needle}` is absent from s59's source -- `best_sequence` would no "
                 f"longer be the selection the measured floors were built with, and the drift "
                 f"would be silent.")
for _needle in ('NONBINDER_crosspaired', 'POSITIVE_crystal_2.05A', 'POSITIVE_crystal_2.27A'):
    if _needle not in _S59_SRC:
        sys.exit(f"⛔ REFUSING: control `{_needle}` is absent from s59's source -- the floor would "
                 f"not be the one `s60` scores against.")


def best_sequence(fa):
    best, hdr = None, None
    for line in open(fa):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            hdr = line
            continue
        if hdr is None or ", id=" not in hdr:
            continue
        conf = None
        for field in hdr.lstrip(">").split(","):
            if "overall_confidence=" in field:
                conf = float(field.split("=")[1])
        if conf is None:
            sys.exit(f"FATAL: design record with no overall_confidence in {fa}")
        if best is None or conf > best[0]:
            best = (conf, line, hdr)
    if best is None:
        sys.exit(f"FATAL: no design records (id=) in {fa}")
    return best


ap = argparse.ArgumentParser()
ap.add_argument("--cohort", required=True)
ap.add_argument("--t149", required=True)
ap.add_argument("--mpnn", required=True)
ap.add_argument("--s35", required=True)
ap.add_argument("--out", required=True)
A = ap.parse_args()

if os.path.exists(A.out) and os.listdir(A.out):
    sys.exit(f"⛔ {A.out} exists and is not empty -- use a fresh dir (STOP S5: scripts APPEND, and "
             f"merging two fold campaigns into one input dir is unrecoverable)")

C = s71.lift_targets(A.s35)

ref = {}
with open(A.t149) as fh:
    for r in csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"):
        ref[r["design"]] = r

cohort = [l.strip() for l in open(A.cohort) if l.strip()]
if not cohort:
    sys.exit("⛔ empty cohort -- refusing to build an empty fold batch")
missing = [d for d in cohort if d not in ref]
if missing:
    sys.exit(f"⛔ {len(missing)} cohort design(s) absent from t149, e.g. {missing[:3]}")

os.makedirs(os.path.join(A.out, "yaml"))   # ⛔ YAMLs in their own dir: a stray non-YAML file inside
                                           # a Boltz input dir killed a launch (`J986` §6).
man = [["tag", "klass", "cell", "design", "batch", "mpnn_id", "mpnn_conf", "cont_pep", "cont_hla",
        "bb_min_to_target", "binder_len", "n_pins", "n_residues_total"]]

for design in sorted(cohort):
    r = ref[design]
    pdb = os.path.join(A.mpnn, "pdb", design + ".pdb")
    fa = os.path.join(A.mpnn, "mpnn", "seqs", design + ".fa")
    for p in (pdb, fa):
        if not os.path.exists(p):
            sys.exit(f"FATAL {design}: missing {p}")
    order, lens = s71.pdb_chains(pdb)
    conf, seq, hdr = best_sequence(fa)
    parts = seq.split(":")

    # ⛔ FATAL CHECKS -- refuse rather than fold the wrong thing.
    if len(parts) != len(order):
        sys.exit(f"FATAL {design}: {len(parts)} fasta parts vs {len(order)} pdb chains")
    if [len(p) for p in parts] != [lens[c] for c in order]:
        sys.exit(f"FATAL {design}: part lengths {[len(p) for p in parts]} vs chains "
                 f"{[lens[c] for c in order]}")
    binder = parts[0]                     # designed chain = A, first by appearance
    if len(binder) != int(r["binder_len"]):
        sys.exit(f"FATAL {design}: binder {len(binder)} vs t149 binder_len {r['binder_len']}")
    if "C" in binder:
        sys.exit(f"FATAL {design}: cysteine in a competition design "
                 f"([[no-cysteines-in-competition-designs]])")
    # ⛔ the design must have been built against the target we are about to fold against
    if parts[1] != C["OUR_HLA"]:
        sys.exit(f"FATAL {design}: design's HLA chain != s35 OUR_HLA")
    if parts[2] != C["OUR_B2M"][1:]:
        sys.exit(f"FATAL {design}: design's b2m chain != OUR_B2M[1:]")
    if parts[3] != C["G12D"]:
        sys.exit(f"FATAL {design}: design's peptide != G12D")

    tag = "D_" + design.replace("arm_shard", "s")
    tot = s71.write_yaml(A.out, tag, [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]),
                                      ("C", C["G12D"]), ("D", binder)])
    mid = [f for f in hdr.lstrip(">").split(",") if " id=" in f or f.strip().startswith("id=")]
    man.append([tag, "DESIGN", r["arm"], design, r["batch"], mid[0].strip() if mid else "?",
                f"{conf:.4f}", r["cont_pep"], r["cont_hla"], r["bb_min_to_target"],
                str(len(binder)), r["n_pins"], str(tot)])

# ---- the four controls, rebuilt in THIS batch (`J986` §3 rule 5) ----
CONTROLS = [
    ("B4_9MIN_own_ss", "POSITIVE_crystal_2.05A",
     [("A", C["MIN_HLA"]), ("B", C["MIN_B2M"]), ("C", C["MIN_PEP"]), ("D", C["MIN_BINDER"])]),
    ("B5_9O5S_own_ss", "POSITIVE_crystal_2.27A",
     [("A", C["O5S_HLA"]), ("B", C["O5S_B2M"]), ("C", C["O5S_PEP"]), ("D", C["O5S_BINDER"])]),
    ("N1_9MINbinder_x_ours_ss", "NONBINDER_crosspaired",
     [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]), ("C", C["G12D"]), ("D", C["MIN_BINDER"])]),
    ("N2_9O5Sbinder_x_ours_ss", "NONBINDER_crosspaired",
     [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]), ("C", C["G12D"]), ("D", C["O5S_BINDER"])]),
]
for tag, klass, chains in CONTROLS:
    tot = s71.write_yaml(A.out, tag, chains)
    man.append([tag, klass, "-", "-", "-", "-", "-", "-", "-", "-",
                str(len(chains[3][1])), "-", str(tot)])

with open(os.path.join(A.out, "MANIFEST_2X2_COFOLD.tsv"), "w") as fh:
    fh.write("\n".join("\t".join(r) for r in man) + "\n")

n_yaml = len([f for f in os.listdir(os.path.join(A.out, "yaml")) if f.endswith(".yaml")])
if n_yaml != len(man) - 1:
    sys.exit(f"⛔ {n_yaml} yaml written against {len(man) - 1} manifest rows -- refusing")
n_nb = sum(1 for r in man[1:] if r[1] == "NONBINDER_crosspaired")
if n_nb != 2:
    sys.exit(f"⛔ {n_nb} cross-paired non-binders -- `J1055` §3 requires the in-batch floor; "
             f"without it there is NO conjunction to evaluate and the round cannot be read")

cells = {}
for r in man[1:]:
    if r[1] == "DESIGN":
        cells[r[2]] = cells.get(r[2], 0) + 1
print(f"✅ {n_yaml} YAMLs -> {A.out}/yaml   ({len(man) - 1} manifest rows)")
print(f"   designs {sum(cells.values())} over {len(cells)} cells: "
      + ", ".join(f"{k} {v}" for k, v in sorted(cells.items())))
print(f"   controls: 2 POSITIVE crystals + 2 NONBINDER_crosspaired (the in-batch floor)")
print(f"   residues per fold: {min(int(r[-1]) for r in man[1:])}-{max(int(r[-1]) for r in man[1:])}")
print(f"⛔ COUNT THE FOLD OUTPUTS AGAINST {len(man) - 1}, never against rc "
      f"([[boltz-oom-exits-zero]])")
