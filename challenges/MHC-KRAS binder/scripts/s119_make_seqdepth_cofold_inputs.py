#!/usr/bin/env python3
"""s119 -- co-fold inputs for the SEQUENCE-DEPTH round: MPNN ranks 2-5 on the no-motif backbones,
plus the two surviving candidates re-folded for extra draws.

    RUN WITH:  python3 -I   (stdlib only; no network, no spend)

⭐⭐ THE QUESTION THIS ROUND ASKS, and why it is the cheapest one left.
   `[M]` `t151`: 2 of 80 designs cleared the pre-registered conjunction. `[M]` But every one of those
   80 folds used **MPNN rank 1 only** (`s116` takes `best_sequence`, the highest
   `overall_confidence`), so **4 of 5 sequences per backbone were never tested** -- 560 of 640
   sequences. The `J1028`-era question *"it may be the sequence, not the backbone"* is therefore
   still open on this cohort, and testing it costs fold time only: ⛔ no new design round, no RFd3.

⛔ SCOPE, STATED SO IT CANNOT BE MISREAD: ranks 2-5 are folded for the **no-motif arms ONLY**
   (`M_nomotif_narrow`, `MB_nomotif_broad`). The held-motif arms are refuted on TWO independent axes
   -- design-stage survival and contact margin (`t150`, `t151`) -- so they are deliberately not
   deepened. ⇒ ⛔ THIS ROUND CANNOT RESCUE THE HELD MOTIF, and a null here says nothing about it.

⭐ THE CANDIDATE RE-DRAW. The 2 survivors rest on **5 diffusion draws each**, and `[M]` re-fold noise
   runs ~**1.56x** the within-batch s.e. ([[within-batch-se-understates-refold-noise]]) ⇒ a +0.6
   margin is well inside noise and even +6.0 is one cohort of draws. Their rank-1 sequences are
   re-emitted here and, because `s117` copies the control set into EVERY shard, they are folded at
   **all 4 seeds** ⇒ 20 draws each instead of 5. ⛔ They are tagged `REDRAW_` and carry klass
   `DESIGN_REDRAW` so they are NEVER pooled with the fresh rank-2-5 designs.

⛔ EVERY GUARD FROM `s116` IS KEPT, because the comparison to `t151` depends on the inputs being
   built the same way: chain remap (design A=binder,B=HLA,C=b2m,D=peptide -> yaml A=HLA,B=b2m,
   C=peptide,D=binder), `msa: empty` on every chain, targets lifted from `s35`'s SOURCE by `ast`,
   no cysteine, and the design's OWN target chains asserted equal to `s35`'s.

⚠️ The manifest is written as `MANIFEST_2X2_COFOLD.tsv` because `s117` reads that exact name. The
   name is `s117`'s contract, ⛔ not a claim that this is the 2x2 round.

Usage:
  python3 -I s119_make_seqdepth_cofold_inputs.py --t149 <t149.tsv> --t151 <t151_per_design.tsv>
          --cohort <s115 cohort.txt> --mpnn <s53 outdir> --s35 <s35 py> --out <fresh dir>
          [--ranks 4]
"""
import argparse
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
S71 = os.path.join(HERE, "s71_make_5b_cofold_inputs.py")
S59 = os.path.join(HERE, "s59_make_r4_cofold_inputs.py")

spec = importlib.util.spec_from_file_location("s71_mod", S71)
s71 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s71)
for need in ("lift_targets", "pdb_chains", "write_yaml"):
    if not hasattr(s71, need):
        sys.exit(f"⛔ s71 lacks `{need}` -- re-audit")

_S59 = open(S59).read()
for _n in ('NONBINDER_crosspaired', 'POSITIVE_crystal_2.05A', 'POSITIVE_crystal_2.27A'):
    if _n not in _S59:
        sys.exit(f"⛔ REFUSING: control `{_n}` absent from s59's source -- the floor would not be "
                 f"the one s60 scores against.")

ap = argparse.ArgumentParser()
ap.add_argument("--t149", required=True)
ap.add_argument("--t151", required=True)
ap.add_argument("--cohort", required=True)
ap.add_argument("--mpnn", required=True)
ap.add_argument("--s35", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--ranks", type=int, default=4, help="how many sequences AFTER rank 1 to fold")
A = ap.parse_args()

if os.path.exists(A.out) and os.listdir(A.out):
    sys.exit(f"⛔ {A.out} exists and is not empty -- use a fresh dir (STOP S5)")

C = s71.lift_targets(A.s35)

ref = {r["design"]: r for r in csv.DictReader(
    [l for l in open(A.t149) if not l.startswith("#")], delimiter="\t")}
prev = {r["tag"]: r for r in csv.DictReader(open(A.t151), delimiter="\t")}
cohort = [l.strip() for l in open(A.cohort) if l.strip()]

# ⛔ no-motif arms only, and the membership test is the ARM NAME from t149, never a filename guess
nomotif = [d for d in sorted(cohort)
           if not ref[d]["arm"].endswith(("_motif_narrow", "_motif_broad"))]
if len(nomotif) != 40:
    sys.exit(f"⛔ {len(nomotif)} no-motif backbones, expected 40 -- refusing (the cohort was "
             f"balanced 20/20/20/20, so a different number means the arm test is wrong)")

# ⭐ the survivors, read from t151 rather than retyped
winners = sorted(t for t, r in prev.items() if r["conj"] == "True")
if not winners:
    sys.exit("⛔ no conjunction survivors in t151 -- nothing to re-draw; re-check the input")


def ranked_sequences(fa):
    """All design records (id=), highest overall_confidence first. ⛔ Record 0 is the INPUT
    sequence and carries no `id=`; folding it would fold the poly-Ala backbone."""
    out, hdr = [], None
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
        for f in hdr.lstrip(">").split(","):
            if "overall_confidence=" in f:
                conf = float(f.split("=")[1])
        if conf is None:
            sys.exit(f"FATAL: design record with no overall_confidence in {fa}")
        out.append((conf, line, hdr))
    out.sort(key=lambda x: -x[0])
    return out


def checked_binder(design, seq):
    """s116's FATAL checks, unchanged."""
    parts = seq.split(":")
    order, lens = s71.pdb_chains(os.path.join(A.mpnn, "pdb", design + ".pdb"))
    if len(parts) != len(order):
        sys.exit(f"FATAL {design}: {len(parts)} fasta parts vs {len(order)} pdb chains")
    if [len(p) for p in parts] != [lens[c] for c in order]:
        sys.exit(f"FATAL {design}: part lengths differ from pdb chains")
    b = parts[0]
    if len(b) != int(ref[design]["binder_len"]):
        sys.exit(f"FATAL {design}: binder {len(b)} vs t149 {ref[design]['binder_len']}")
    if "C" in b:
        sys.exit(f"FATAL {design}: cysteine in a competition design")
    if parts[1] != C["OUR_HLA"] or parts[2] != C["OUR_B2M"][1:] or parts[3] != C["G12D"]:
        sys.exit(f"FATAL {design}: design's own target chains != s35 source")
    return b


os.makedirs(os.path.join(A.out, "yaml"))
man = [["tag", "klass", "cell", "design", "batch", "mpnn_rank", "mpnn_conf", "binder_len",
        "n_residues_total"]]

n_new = 0
for design in nomotif:
    recs = ranked_sequences(os.path.join(A.mpnn, "mpnn", "seqs", design + ".fa"))
    if len(recs) < A.ranks + 1:
        sys.exit(f"FATAL {design}: {len(recs)} sequences, need {A.ranks + 1}")
    for rank, (conf, seq, hdr) in enumerate(recs[1:A.ranks + 1], start=2):
        b = checked_binder(design, seq)
        tag = f"R{rank}_" + design.replace("arm_shard", "s")
        tot = s71.write_yaml(A.out, tag, [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]),
                                          ("C", C["G12D"]), ("D", b)])
        man.append([tag, "DESIGN", ref[design]["arm"], design, ref[design]["batch"], str(rank),
                    f"{conf:.4f}", str(len(b)), str(tot)])
        n_new += 1

# ⭐ the re-draws go in as CONTROLS so s117 copies them into every shard -> all 4 seeds
redraw = []
for tag151 in winners:
    design = tag151[2:].replace("s0", "arm_shard0", 1) if tag151.startswith("D_") else None
    cand = [d for d in cohort if tag151 == "D_" + d.replace("arm_shard", "s")]
    if len(cand) != 1:
        sys.exit(f"⛔ could not map t151 tag {tag151!r} back to exactly one cohort design "
                 f"(got {cand}) -- refusing to guess")
    design = cand[0]
    recs = ranked_sequences(os.path.join(A.mpnn, "mpnn", "seqs", design + ".fa"))
    conf, seq, hdr = recs[0]
    b = checked_binder(design, seq)
    tag = "REDRAW_" + design.replace("arm_shard", "s")
    tot = s71.write_yaml(A.out, tag, [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]),
                                      ("C", C["G12D"]), ("D", b)])
    redraw.append([tag, "DESIGN_REDRAW", ref[design]["arm"], design, ref[design]["batch"], "1",
                   f"{conf:.4f}", str(len(b)), str(tot)])

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
    man.append([tag, klass, "-", "-", "-", "-", "-", str(len(chains[3][1])), str(tot)])
man.extend(redraw)

with open(os.path.join(A.out, "MANIFEST_2X2_COFOLD.tsv"), "w") as fh:
    fh.write("\n".join("\t".join(r) for r in man) + "\n")

n_yaml = len([f for f in os.listdir(os.path.join(A.out, "yaml")) if f.endswith(".yaml")])
if n_yaml != len(man) - 1:
    sys.exit(f"⛔ {n_yaml} yaml against {len(man) - 1} manifest rows -- refusing")
if sum(1 for r in man[1:] if r[1] == "NONBINDER_crosspaired") != 2:
    sys.exit("⛔ the in-batch floor is missing -- without it there is NO conjunction to evaluate")

cells = {}
for r in man[1:]:
    if r[1] == "DESIGN":
        cells[r[2]] = cells.get(r[2], 0) + 1
print(f"✅ {n_yaml} YAMLs -> {A.out}/yaml   ({len(man) - 1} manifest rows)")
print(f"   fresh designs {n_new} (ranks 2-{A.ranks + 1}) over {len(cells)} no-motif arms: "
      + ", ".join(f"{k} {v}" for k, v in sorted(cells.items())))
print(f"   re-draws {len(redraw)}: " + ", ".join(r[0] for r in redraw))
print(f"   + 2 POSITIVE crystals + 2 NONBINDER_crosspaired (the in-batch floor)")
print(f"⛔ s117 copies every non-DESIGN row into EVERY shard ⇒ re-draws and controls are folded at "
      f"all 4 seeds. Folds = {n_new} + {(len(redraw) + 4) * 4} = {n_new + (len(redraw) + 4) * 4}")
