#!/usr/bin/env python3
"""s71 -- emit the Step 5B Boltz co-fold inputs for the CDR3β cohort, plus the controls.

    RUN WITH:  python3 -I   (stdlib only; no network, no spend)

⭐⭐ WHY THIS IS NOT `s59`. `CMP-11_00_PLAN.md` §5B says to run `s59_make_r4_cofold_inputs.py`.
   `[M-me]` **`s59` CANNOT BUILD THIS BATCH AND WILL HARD-EXIT**, for three independent reasons
   ([[check-a-bundle-against-its-runners-own-guards]]):
     1. `s59:114` -- `if len(gated) != 128: sys.exit("FATAL: expected the 128-design r4 cohort")`.
        Our cohort is `[M-me]` 6 backbones, not 128 designs.
     2. `s59:124` -- `if len(by_cell) != 6` over the `arm` column. Our 6 backbones span `[M-me]`
        THREE cells, not six.
     3. `s59:79 best_sequence()` returns **ONE** sequence per design (the highest
        `overall_confidence`). Step 5B as ruled folds **22 sequences from 6 backbones**, i.e.
        several sequences per backbone -- a shape `s59` cannot express at all.
   ⇒ ⛔ a session that follows the plan verbatim meets `FATAL` **after** the pod is up and billing.

⭐ WHAT IS REUSED FROM `s59` DELIBERATELY, because it is proven and comparability depends on it:
   - the chain REMAP, which is ⛔ not cosmetic [[always-name-the-pdb-chain]]:
        design PDB / MPNN FASTA order:  A=binder  B=HLA  C=b2m  D=peptide
        fold YAML order (the study's):  A=HLA     B=b2m  C=peptide  D=binder
   - `msa: empty` on EVERY chain (`J986` §3 rule 1: `[M]` 2.3x better class separation, and the
     MSA's extra confidence lands where the model is wrong);
   - the four CONTROLS rebuilt in THIS batch (`J986` §3 rule 5 -- the floor moves with the arm,
     so it is never read off an older table);
   - target sequences lifted from `s35_make_g9_inputs.py`'s SOURCE by `ast`, so there is no
     transcription copy to drift [[never-source-target-sequence-from-model-output]];
   - ⚠️ `OUR_B2M` keeps the initiator Met in the YAML although the design PDBs' b2m starts one
     residue later (`OUR_B2M[1:]`). That asymmetry is `s59`'s and is KEPT, because the measured
     G9 floor was built that way and changing it would break comparability.

⛔⛔ WHAT THIS FOLD CAN AND CANNOT SAY -- unchanged from `s59`, restated so it cannot be lost:
   - It ⛔ CANNOT rank peptide-IDENTITY selectivity (rung ⑶). `J989` §1: no bound-state interface
     metric from ANY model ranks WT-vs-G12D here, and the CRYSTALS invert the same way.
   - The designed hydrogen bond is a HYPOTHESIS this fold TESTS, ⛔ not delivered geometry: Boltz
     takes sequence only and builds its own side chains, and the guidepost weld is broken in 96/96
     pinned designs [[rfd3-guidepost-functional-group-is-not-bonded-on]].
   - The 22 were ranked partly on design-stage quantities, so a high folded peptide-contact count
     is ⛔ NOT independent evidence of peptide engagement -- they were selected for it.
   - ⇒ what 5B actually tests is whether the SEQUENCE rebuilds the designed POSE
     [[design-vs-prediction-selfconsistency]].

⛔ COUNT THE OUTPUTS, NEVER THE EXIT CODE: a Boltz CUDA OOM exits **rc=0** having written nothing
   [[boltz-oom-exits-zero]]. The manifest row count printed at the end is the number to check the
   fold outputs against.

Usage: s71_make_5b_cofold_inputs.py --picks <t134 tsv> --run <5A outdir> --s35 <s35 py>
                                    --out <fresh dir>
"""
import argparse
import ast
import collections
import csv
import os
import sys

TARGET_KEYS = ("OUR_HLA", "OUR_B2M", "G12D", "MIN_HLA", "MIN_B2M", "MIN_PEP", "MIN_BINDER",
               "O5S_HLA", "O5S_B2M", "O5S_PEP", "O5S_BINDER")


def lift_targets(s35):
    """Target sequences from s35's SOURCE. ⛔ never from a model output."""
    C = {}
    for node in ast.parse(open(s35).read()).body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
            try:
                C[node.targets[0].id] = ast.literal_eval(node.value)
            except Exception:
                pass
    for k in TARGET_KEYS:
        if k not in C or not isinstance(C[k], str) or not C[k]:
            sys.exit("⛔ could not lift %s from %s. Refusing." % (k, s35))
    return C


def pdb_chains(path):
    """(chain_order, {chain: n_residues}) by first appearance."""
    seen = collections.OrderedDict()
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                seen.setdefault(line[21], set()).add(line[22:27])
    return list(seen.keys()), {c: len(v) for c, v in seen.items()}


def write_yaml(out, tag, chains):
    body = "".join(
        "  - protein:\n      id: %s\n      sequence: %s\n      msa: empty\n" % (cid, seq)
        for cid, seq in chains)
    with open(os.path.join(out, "yaml", tag + ".yaml"), "w") as fh:
        fh.write("version: 1\nsequences:\n" + body)
    return sum(len(s) for _, s in chains)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--picks", required=True, help="t134 from s70: the 22 ruled sequences")
    ap.add_argument("--run", required=True, help="the Step 5A outdir (pdb/ + report.json)")
    ap.add_argument("--s35", required=True, help="s35_make_g9_inputs.py -- the target source")
    ap.add_argument("--out", required=True, help="⛔ fresh dir; refused if it exists (STOP S5)")
    a = ap.parse_args()

    # ⛔ STOP S5: a script that appends or overwrites silently merges two campaigns.
    if os.path.exists(a.out):
        sys.exit("⛔ %s exists. Use a fresh --out per experiment (STOP S5)." % a.out)

    C = lift_targets(a.s35)
    rows = [r for r in csv.DictReader(
        [l for l in open(a.picks) if not l.startswith("#")], delimiter="\t")]
    if not rows:
        sys.exit("⛔ %s has no rows. Refusing." % a.picks)

    # ⛔ YAMLs in their own dir: a stray non-YAML file inside a Boltz input dir killed a
    # launch (`J986` §6).
    os.makedirs(os.path.join(a.out, "yaml"))

    man = [["tag", "klass", "backbone", "seq_id", "mpnn_conf", "bsa_pep", "bsa_nopep",
            "binder_len", "n_residues_total"]]
    seen_tags = set()

    for r in rows:
        design, sid, binder = r["backbone"], r["seq_id"], r["sequence"]

        # the design must have been built against the TARGET we are about to fold against
        order, lens = pdb_chains(os.path.join(a.run, "pdb", design + ".pdb"))
        if order[:4] != ["A", "B", "C", "D"]:
            sys.exit("⛔ %s: chain order %s is not A,B,C,D. Refusing." % (design, order))
        if len(binder) != lens["A"]:
            sys.exit("⛔ %s id=%s: sequence %d vs chain A %d residues. Refusing."
                     % (design, sid, len(binder), lens["A"]))
        if len(binder) != int(r["seq_len"]):
            sys.exit("⛔ %s id=%s: seq_len column disagrees with the sequence. Refusing."
                     % (design, sid))
        if "C" in binder:
            sys.exit("⛔ %s id=%s: cysteine in a competition design "
                     "([[no-cysteines-in-competition-designs]]). Refusing." % (design, sid))
        if lens["B"] != len(C["OUR_HLA"]):
            sys.exit("⛔ %s: HLA chain %d vs s35 OUR_HLA %d. Refusing."
                     % (design, lens["B"], len(C["OUR_HLA"])))
        if lens["C"] != len(C["OUR_B2M"]) - 1:
            sys.exit("⛔ %s: b2m chain %d vs OUR_B2M[1:] %d. Refusing."
                     % (design, lens["C"], len(C["OUR_B2M"]) - 1))
        if lens["D"] != len(C["G12D"]):
            sys.exit("⛔ %s: peptide chain %d vs G12D %d. Refusing."
                     % (design, lens["D"], len(C["G12D"])))

        tag = "D_" + design.replace("arm_shard", "s") + "_id" + str(sid)
        if tag in seen_tags:
            sys.exit("⛔ duplicate tag %s. Refusing to overwrite a fold input." % tag)
        seen_tags.add(tag)

        tot = write_yaml(a.out, tag, [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]),
                                      ("C", C["G12D"]), ("D", binder)])
        man.append([tag, "DESIGN", design, sid, r["overall_confidence"],
                    r["bsa_pep"], r["bsa_nopep"], str(len(binder)), str(tot)])

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
        tot = write_yaml(a.out, tag, chains)
        man.append([tag, klass, "-", "-", "-", "-", "-", str(len(chains[3][1])), str(tot)])

    with open(os.path.join(a.out, "manifest.tsv"), "w") as fh:
        for row in man:
            fh.write("\t".join(row) + "\n")

    nyaml = len(os.listdir(os.path.join(a.out, "yaml")))
    toks = [int(row[-1]) for row in man[1:]]
    print("# s71 — Step 5B co-fold inputs")
    print("   designs   : %d" % (len(man) - 1 - len(CONTROLS)))
    print("   controls  : %d  (2 positive crystals + 2 cross-paired non-binders)"
          % len(CONTROLS))
    print("   YAMLs     : %d   manifest rows: %d" % (nyaml, len(man) - 1))
    print("   tokens/fold: min %d  median %d  max %d"
          % (min(toks), sorted(toks)[len(toks) // 2], max(toks)))
    print()
    print("   ⛔ COUNT THE FOLD OUTPUTS AGAINST %d PER SEED — a Boltz OOM exits rc=0 having"
          % (len(man) - 1))
    print("      written nothing [[boltz-oom-exits-zero]].")
    print("   ⛔ 2 seeds as ruled ⇒ %d folds total." % (2 * (len(man) - 1)))
    print("   wrote %s" % a.out)


if __name__ == "__main__":
    main()
