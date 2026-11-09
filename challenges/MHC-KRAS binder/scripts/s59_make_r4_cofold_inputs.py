#!/usr/bin/env python3
"""Emit the Boltz co-fold inputs for the ADOPTED r4 MPNN cohort, plus the controls.

This is the fold step J993 §6 item 1 owes. It takes the adopted bias cohort
(`$LAB_SCRATCH/c5_mpnn_20261103/r4`, bias A:-0.495,E:-0.354,L:-0.092 under
`--hold D,K,R --negative-only --omit_AA C`) and asks what Boltz builds from SEQUENCE ALONE.

⛔⛔ WHAT THIS FOLD CAN AND CANNOT SAY
  - It CANNOT rank step-1 selectivity. `J989` §1: no bound-state interface metric from ANY model
    ranks WT-vs-G12D here, and the CRYSTALS invert the same way. ⛔ Never use these numbers for it.
  - The designed hydrogen bond is a HYPOTHESIS this fold TESTS, not delivered geometry
    (`J992` §2 -- the pin weld is intact in 17 % of pins against a matched 100 % control).
    Boltz takes sequence only and builds its own side chains.
  - ⛔ The selection below ranks on `cont_pep`, a DESIGN-STAGE backbone contact count. So a high
    folded peptide-contact count is NOT independent evidence that the designs engage the peptide --
    they were chosen for that. What the fold tests is whether the SEQUENCE rebuilds the pose.

⭐ OPERATING RULES APPLIED, from `J986` §3:
  1. Single-sequence, never MSA -- `msa: empty` on EVERY chain, both classes. [M] 2.3x better
     class separation, and the MSA's extra confidence lands where the model is wrong.
  5. Cross-paired non-binders in the SAME batch -- the floor moves with the arm, so it is rebuilt
     here rather than read off `t62`.
  7. Count the output files -- a Boltz CUDA OOM exits rc=0 having written nothing
     ([[boltz-oom-exits-zero]]). The manifest row count is the number to check against.

⛔ CHAIN IDS ARE REMAPPED AND THIS IS NOT COSMETIC ([[always-name-the-pdb-chain]]):
      design PDB / MPNN FASTA order:  A=binder  B=HLA  C=b2m  D=peptide
      fold YAML order (study's own):  A=HLA     B=b2m  C=peptide  D=binder
  The YAML follows the STUDY's convention so `s36`/`s37` and the `t62` floor apply unchanged.

⛔ Target sequences are NOT read from any model output ([[never-source-target-sequence-from-model-output]]).
  They are lifted from `s35_make_g9_inputs.py` by parsing its source, so that file stays the single
  source of truth and there is no transcription copy to drift. `s35` took them from SEQRES.
  ⚠️ `OUR_B2M` carries the initiator Met; the design PDBs' b2m chain starts one residue later
  (`OUR_B2M[1:]`, from 7OW6). The full Met form is used here DELIBERATELY, because that is what the
  G9 batch folded -- changing it would break comparability with the measured floor.

Usage: python3 -I s59_make_r4_cofold_inputs.py <r4_cohort_dir> <t72_tsv> <s35_py> <out_dir>
"""
import ast
import collections
import csv
import os
import sys

R4, T72, S35, OUT = sys.argv[1:5]

N_PER_CELL = 3        # [M] 18 designs + 4 controls = 22 folds; ~140 s/fold => ~51 min, inside the
                      # A4000's "minutes to ~1 hour" rule ([[a4000-short-jobs-only]]).

# ---- STOP S5: a script that appends or overwrites silently merges two campaigns ----
if os.path.exists(OUT):
    sys.exit(f"FATAL: {OUT} exists. Use a fresh --out per experiment (STOP S5).")

# ---- the target, from s35's source (never from a model output) ----
C = {}
for node in ast.parse(open(S35).read()).body:
    if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
        try:
            C[node.targets[0].id] = ast.literal_eval(node.value)
        except Exception:
            pass
for k in ("OUR_HLA", "OUR_B2M", "G12D", "MIN_HLA", "MIN_B2M", "MIN_PEP", "MIN_BINDER",
          "O5S_HLA", "O5S_B2M", "O5S_PEP", "O5S_BINDER"):
    if k not in C or not isinstance(C[k], str) or not C[k]:
        sys.exit(f"FATAL: could not lift {k} from {S35}")


def pdb_chains(path):
    """(chain_order, {chain: n_residues}) by first appearance."""
    seen = collections.OrderedDict()
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                seen.setdefault(line[21], set()).add(line[22:27])
    return list(seen.keys()), {c: len(v) for c, v in seen.items()}


def best_sequence(fa):
    """The highest-overall_confidence DESIGN from an MPNN fasta.

    ⛔ Record 0 is the INPUT sequence and carries no `id=` -- it must be skipped, exactly as
    `j745_tfr1_composition_bias_20261027.py`:78 skips it. Taking it would fold the poly-Ala
    backbone instead of a design.
    """
    best = None
    hdr = None
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


# ---- the gated universe and the PRE-REGISTERED selection ----
rows = [r for r in csv.DictReader(
    [l for l in open(T72) if not l.startswith("#")], delimiter="\t")]
gated = [r for r in rows if r["continuous"] == "1" and r["pass_raw_clash"] == "1"]
if len(gated) != 128:
    sys.exit(f"FATAL: expected the 128-design r4 cohort, the gates give {len(gated)}")

# ⭐ PRE-REGISTERED, stated before any fold runs: within each of the 6 cells rank on the
# peptide-STRATIFIED contact count (⛔ never pooled -- `J986` §3 rule 3: the highest HLA-contact
# count in the whole G9 table belongs to a NON-BINDER), tie-break on target clearance, then on name
# for determinism. Take the top N_PER_CELL. Balanced by construction, so no cell can dominate.
by_cell = collections.defaultdict(list)
for r in gated:
    by_cell[r["arm"]].append(r)
if len(by_cell) != 6:
    sys.exit(f"FATAL: expected 6 cells, found {sorted(by_cell)}")

picked = []
for cell in sorted(by_cell):
    ranked = sorted(by_cell[cell],
                    key=lambda r: (-int(r["cont_pep"]), -float(r["bb_min_to_target"]), r["design"]))
    picked.extend((cell, r) for r in ranked[:N_PER_CELL])
if len(picked) != 6 * N_PER_CELL:
    sys.exit(f"FATAL: picked {len(picked)}, expected {6 * N_PER_CELL}")

os.makedirs(os.path.join(OUT, "yaml"))   # ⛔ YAMLs in their own dir: a stray non-YAML file inside
                                         # a Boltz input dir killed a launch (`J986` §6).
man = [["tag", "klass", "cell", "design", "mpnn_id", "mpnn_conf", "cont_pep", "cont_hla",
        "bb_min_to_target", "binder_len", "n_pins", "n_residues_total"]]


def write_yaml(tag, chains):
    """chains = [(chain_id, sequence)]; every chain single-sequence (`J986` §3 rule 1)."""
    body = "".join(
        f"  - protein:\n      id: {cid}\n      sequence: {seq}\n      msa: empty\n"
        for cid, seq in chains)
    with open(os.path.join(OUT, "yaml", tag + ".yaml"), "w") as fh:
        fh.write("version: 1\nsequences:\n" + body)
    return sum(len(s) for _, s in chains)


for cell, r in picked:
    design = r["design"]
    order, lens = pdb_chains(os.path.join(R4, "pdb", design + ".pdb"))
    conf, seq, hdr = best_sequence(os.path.join(R4, "mpnn", "seqs", design + ".fa"))
    parts = seq.split(":")

    # ⛔ FATAL CHECKS -- refuse rather than fold the wrong thing.
    if len(parts) != len(order):
        sys.exit(f"FATAL {design}: {len(parts)} fasta parts vs {len(order)} pdb chains")
    if [len(p) for p in parts] != [lens[c] for c in order]:
        sys.exit(f"FATAL {design}: part lengths {[len(p) for p in parts]} vs chains "
                 f"{[lens[c] for c in order]}")
    binder = parts[0]                     # designed_chain = A, first by appearance
    if len(binder) != int(r["binder_len"]):
        sys.exit(f"FATAL {design}: binder {len(binder)} vs t72 binder_len {r['binder_len']}")
    if "C" in binder:
        sys.exit(f"FATAL {design}: cysteine in a competition design "
                 f"([[no-cysteines-in-competition-designs]])")
    # the design must have been built against the target we are about to fold against
    if parts[1] != C["OUR_HLA"]:
        sys.exit(f"FATAL {design}: design's HLA chain != s35 OUR_HLA")
    if parts[2] != C["OUR_B2M"][1:]:
        sys.exit(f"FATAL {design}: design's b2m chain != OUR_B2M[1:]")
    if parts[3] != C["G12D"]:
        sys.exit(f"FATAL {design}: design's peptide != G12D")

    tag = "D_" + design.replace("arm_shard", "s")
    tot = write_yaml(tag, [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]),
                           ("C", C["G12D"]), ("D", binder)])
    mid = [f for f in hdr.lstrip(">").split(",") if " id=" in f or f.strip().startswith("id=")]
    man.append([tag, "DESIGN", cell, design, mid[0].strip() if mid else "?", f"{conf:.4f}",
                r["cont_pep"], r["cont_hla"], r["bb_min_to_target"], str(len(binder)),
                r["n_pins"], str(tot)])

# ---- controls, rebuilt in THIS batch (`J986` §3 rule 5) ----
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
    tot = write_yaml(tag, chains)
    man.append([tag, klass, "-", "-", "-", "-", "-", "-", "-",
                str(len(chains[3][1])), "-", str(tot)])

with open(os.path.join(OUT, "MANIFEST_R4_COFOLD.tsv"), "w") as fh:
    fh.write("\n".join("\t".join(r) for r in man) + "\n")

n_yaml = len([f for f in os.listdir(os.path.join(OUT, "yaml")) if f.endswith(".yaml")])
if n_yaml != len(man) - 1:
    sys.exit(f"FATAL: {n_yaml} yaml files vs {len(man)-1} manifest rows")

print(f"wrote {n_yaml} YAMLs -> {os.path.join(OUT, 'yaml')}")
print(f"manifest ({len(man)-1} rows) -> {os.path.join(OUT, 'MANIFEST_R4_COFOLD.tsv')}")
print(f"\n⛔ EXPECT EXACTLY {n_yaml} prediction dirs. Fewer = a silent OOM (rc=0, wrote nothing).\n")
cells = collections.Counter(r[2] for r in man[1:] if r[1] == "DESIGN")
for cell in sorted(cells):
    print(f"  {cell:<30} {cells[cell]} designs")
print(f"  {'CONTROLS':<30} {sum(1 for r in man[1:] if r[1] != 'DESIGN')} "
      f"(2 crystal-true binders, 2 cross-paired non-binders)")
tot = [int(r[11]) for r in man[1:]]
print(f"\n  residues per fold: min {min(tot)}  max {max(tot)}")
