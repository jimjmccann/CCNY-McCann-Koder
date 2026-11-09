#!/usr/bin/env python3
"""s123 -- co-fold inputs for the MB-DEPTH round (`J1058`): MPNN ranks 1-5 on the NEW MB backbones,
plus a re-draw of the two designs `t152` rescued but could only afford 5 draws for.

    RUN WITH:  python3 -I   (stdlib only; no network, no spend)

⭐⭐ WHY THIS IS A SIBLING OF `s119` AND NOT AN EDIT TO IT. `s119` produced `t152`; editing it would
   make `t152` unreproducible. It also cannot be pointed at this round even if one wanted to:
   `[M]` `s119:85-87` hard-asserts **exactly 40** no-motif backbones across TWO arms, and
   `s119:91-92` requires conjunction survivors in a `t151`-shaped file. Same choice, same reason, as
   `s120` being a sibling of `s118`.

⛔⛔ THE TWO THINGS THIS CHANGES FROM `s119`, BOTH OF WHICH WOULD BE SILENT IF GOT WRONG:
 ⑴ **RANKS 1-5, NOT 2-5.** `s119`'s `--ranks` means *"how many sequences AFTER rank 1"* because
    `t151` had already folded rank 1 for its cohort. ⛔ For this round **no rank has been folded**,
    so rank 1 must be INCLUDED or the best sequence per backbone is simply never tested.
 ⑵ **A RE-DRAW MUST RE-FOLD THE SEQUENCE THAT CLEARED, AT ITS OWN RANK.** `s119` took `recs[0]`
    (rank 1) for its re-draws, which was correct there because `t151`'s survivors cleared on rank 1.
    ⛔⛔ Here the two designs to re-draw cleared at **ranks 2 and 5**. Folding their rank-1 sequence
    would re-draw **a different molecule** and the result would look like a refutation of a design
    that was never tested. ⇒ the rank is read from `t152`, never defaulted.

⭐ WHICH DESIGNS GET RE-DRAWN IS DERIVED, ⛔ NOT RETYPED. Rule: a `t152` row with `conj == True`
   whose tag is a FRESH rank (`R<n>_`, not `REDRAW_`), for a design that has **no `REDRAW_` rows in
   `t152` at all** -- i.e. it was rescued by a non-rank-1 sequence and has only its original 5
   draws. `[M]` On `t152` that yields exactly 2: `..._19_model_0` @ rank 2 and
   `..._8_model_2` @ rank 5. ⛔ `..._2_model_3` and `..._4_model_0` are EXCLUDED because they
   already carry 20-draw REDRAW rows. The expected count is passed in and asserted, so a change in
   the input cannot slip through as a different cohort.

⚠️ CROSS-ROUND DEPENDENCY, STATED BECAUSE IT IS THE FRAGILE PART: the re-draw sequences live in the
   PREVIOUS round's `s53` outdir and are checked against the PREVIOUS round's `t149`. Both are
   separate flags (`--redraw-mpnn`, `--redraw-t149`) so they can never be silently satisfied by this
   round's files. `[M]` Verified present 2026-11-08 under `$LAB_SCRATCH/c5_2x2_mpnn_v1_20261108`.

⛔ EVERY GUARD FROM `s116`/`s119` IS KEPT, because comparability to `t151`/`t152` depends on the
   inputs being built the same way: the chain remap (design A=binder,B=HLA,C=b2m,D=peptide -> yaml
   A=HLA,B=b2m,C=peptide,D=binder), `msa: empty` on every chain, targets lifted from `s35`'s SOURCE
   by `ast`, no cysteine, and each design's OWN target chains asserted equal to `s35`'s.

⚠️ The manifest is written as `MANIFEST_2X2_COFOLD.tsv` because `s117` reads that exact name. That
   is `s117`'s contract, ⛔ not a claim that this is the 2x2 round.

Usage:
  python3 -I s123_make_mb_depth_cofold_inputs.py --t149 <new t149.tsv> --cohort <s115 cohort.txt>
          --mpnn <new s53 outdir> --s35 <s35 py> --out <fresh dir>
          --expect-backbones <n> --redraw-t152 <t152_per_design.tsv>
          --redraw-t149 <old t149.tsv> --redraw-mpnn <old s53 outdir>
          [--ranks 5] [--expect-redraw 2]
"""
import argparse
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
S71 = os.path.join(HERE, "s71_make_5b_cofold_inputs.py")
S59 = os.path.join(HERE, "s59_make_r4_cofold_inputs.py")
S119 = os.path.join(HERE, "s119_make_seqdepth_cofold_inputs.py")

spec = importlib.util.spec_from_file_location("s71_mod", S71)
s71 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s71)
for need in ("lift_targets", "pdb_chains", "write_yaml"):
    if not hasattr(s71, need):
        sys.exit(f"⛔ s71 lacks `{need}` -- re-audit")

# ⛔ the floor controls must be the ones `s60`/`s118`/`s120` score against, asserted at s59's source
_S59 = open(S59).read()
for _n in ('NONBINDER_crosspaired', 'POSITIVE_crystal_2.05A', 'POSITIVE_crystal_2.27A'):
    if _n not in _S59:
        sys.exit(f"⛔ REFUSING: control `{_n}` absent from s59's source -- the floor would not be "
                 f"the one s60 scores against.")

# ⛔⛔ assert the premise this script is built on: that s119 really does fold ranks AFTER rank 1 and
#    really does take recs[0] for its re-draws. If s119 ever changes, the ⑴/⑵ reasoning above is
#    stale and this script's divergence from it is no longer deliberate.
_S119 = open(S119).read()
for _n in ('how many sequences AFTER rank 1 to fold', 'recs[1:A.ranks + 1]'):
    if _n not in _S119:
        sys.exit(f"⛔ REFUSING: `{_n}` is absent from s119's source -- this script documents itself "
                 f"as diverging from s119 on ranks, and that premise no longer holds. Re-audit.")

ap = argparse.ArgumentParser()
ap.add_argument("--t149", required=True)
ap.add_argument("--cohort", required=True)
ap.add_argument("--mpnn", required=True)
ap.add_argument("--s35", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--expect-backbones", type=int, required=True,
                help="⛔ stated explicitly; a mismatch REFUSES rather than scoring a short cohort")
ap.add_argument("--redraw-t152", required=True)
ap.add_argument("--redraw-t149", required=True)
ap.add_argument("--redraw-mpnn", required=True)
ap.add_argument("--ranks", type=int, default=5, help="fold ranks 1..N INCLUSIVE (⛔ unlike s119)")
ap.add_argument("--expect-redraw", type=int, default=2)
A = ap.parse_args()

if os.path.exists(A.out) and os.listdir(A.out):
    sys.exit(f"⛔ {A.out} exists and is not empty -- use a fresh dir (STOP S5)")

C = s71.lift_targets(A.s35)

ref = {r["design"]: r for r in csv.DictReader(
    [l for l in open(A.t149) if not l.startswith("#")], delimiter="\t")}
ref_prev = {r["design"]: r for r in csv.DictReader(
    [l for l in open(A.redraw_t149) if not l.startswith("#")], delimiter="\t")}
cohort = [l.strip() for l in open(A.cohort) if l.strip()]

# ⛔ this round is ONE cell. The membership test is the ARM NAME from t149, never a filename guess.
bad = sorted({ref[d]["arm"] for d in cohort if ref[d]["arm"] != "MB_nomotif_broad"})
if bad:
    sys.exit(f"⛔ cohort contains non-MB arms {bad} -- this round is MB_nomotif_broad only, and a "
             f"stray arm would be pooled into a single-cell depth round.")
if len(cohort) != A.expect_backbones:
    sys.exit(f"⛔ {len(cohort)} backbones in the cohort, expected {A.expect_backbones} -- refusing. "
             f"A scorer pointed at a short census reports a clean number for an incomplete arm.")

# ── the re-draw set, DERIVED from t152 ────────────────────────────────────────────────────────
t152 = list(csv.DictReader(open(A.redraw_t152), delimiter="\t"))
already = {r["design"] for r in t152 if r["tag"].startswith("REDRAW_")}
cand = {}
for r in t152:
    if r["conj"] != "True" or r["tag"].startswith("REDRAW_"):
        continue
    if r["design"] in already:
        continue                      # ⛔ already validated at 20 draws; not a re-draw target
    k = r["design"]
    m = float(r["margin"])
    if k not in cand or m > cand[k][0]:
        cand[k] = (m, int(r["rank"]))
if len(cand) != A.expect_redraw:
    sys.exit(f"⛔ derived {len(cand)} re-draw targets {sorted(cand)}, expected {A.expect_redraw}. "
             f"⛔ Refusing to guess -- re-read t152 rather than loosening this.")
print(f"✅ re-draw targets DERIVED from t152 (⛔ not retyped): "
      + ", ".join(f"{d} @ rank {rk} (margin {m:+.1f})" for d, (m, rk) in sorted(cand.items())))


def ranked_sequences(fa):
    """All design records (id=), highest overall_confidence first. ⛔ Record 0 of the fasta is the
    INPUT sequence and carries no `id=`; folding it would fold the poly-Ala backbone."""
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


def checked_binder(design, seq, mpnn_dir, reftab):
    """s116's FATAL checks, unchanged -- but parameterised by WHICH round's files to check against,
    because the re-draws come from the previous round."""
    parts = seq.split(":")
    order, lens = s71.pdb_chains(os.path.join(mpnn_dir, "pdb", design + ".pdb"))
    if len(parts) != len(order):
        sys.exit(f"FATAL {design}: {len(parts)} fasta parts vs {len(order)} pdb chains")
    if [len(p) for p in parts] != [lens[c] for c in order]:
        sys.exit(f"FATAL {design}: part lengths differ from pdb chains")
    b = parts[0]
    if len(b) != int(reftab[design]["binder_len"]):
        sys.exit(f"FATAL {design}: binder {len(b)} vs t149 {reftab[design]['binder_len']}")
    if "C" in b:
        sys.exit(f"FATAL {design}: cysteine in a competition design")
    if parts[1] != C["OUR_HLA"] or parts[2] != C["OUR_B2M"][1:] or parts[3] != C["G12D"]:
        sys.exit(f"FATAL {design}: design's own target chains != s35 source")
    return b


os.makedirs(os.path.join(A.out, "yaml"))
man = [["tag", "klass", "cell", "design", "batch", "mpnn_rank", "mpnn_conf", "binder_len",
        "n_residues_total"]]

# ── fresh designs: ranks 1..N INCLUSIVE ───────────────────────────────────────────────────────
n_new = 0
for design in sorted(cohort):
    recs = ranked_sequences(os.path.join(A.mpnn, "mpnn", "seqs", design + ".fa"))
    if len(recs) < A.ranks:
        sys.exit(f"FATAL {design}: {len(recs)} sequences, need {A.ranks} "
                 f"(raise s53's --batches, ⛔ never --seqs, which is a VRAM limit)")
    for rank, (conf, seq, hdr) in enumerate(recs[:A.ranks], start=1):
        b = checked_binder(design, seq, A.mpnn, ref)
        tag = f"R{rank}_" + design.replace("arm_shard", "s")
        tot = s71.write_yaml(A.out, tag, [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]),
                                          ("C", C["G12D"]), ("D", b)])
        man.append([tag, "DESIGN", ref[design]["arm"], design, ref[design]["batch"], str(rank),
                    f"{conf:.4f}", str(len(b)), str(tot)])
        n_new += 1

# ── the re-draws go in as CONTROLS so s117 copies them into EVERY shard -> all 4 seeds ────────
redraw = []
for design, (margin, want_rank) in sorted(cand.items()):
    recs = ranked_sequences(os.path.join(A.redraw_mpnn, "mpnn", "seqs", design + ".fa"))
    if len(recs) < want_rank:
        sys.exit(f"FATAL {design}: {len(recs)} sequences but t152 cleared at rank {want_rank}")
    conf, seq, hdr = recs[want_rank - 1]          # ⛔⛔ THE RANK THAT CLEARED, not rank 1
    b = checked_binder(design, seq, A.redraw_mpnn, ref_prev)
    tag = f"REDRAW_R{want_rank}_" + design.replace("arm_shard", "s")
    tot = s71.write_yaml(A.out, tag, [("A", C["OUR_HLA"]), ("B", C["OUR_B2M"]),
                                      ("C", C["G12D"]), ("D", b)])
    redraw.append([tag, "DESIGN_REDRAW", ref_prev[design]["arm"], design,
                   ref_prev[design]["batch"], str(want_rank), f"{conf:.4f}", str(len(b)), str(tot)])

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
if len(redraw) != A.expect_redraw:
    sys.exit(f"⛔ {len(redraw)} re-draw rows written, expected {A.expect_redraw}")

print(f"✅ {n_yaml} YAMLs -> {A.out}/yaml   ({len(man) - 1} manifest rows)")
print(f"   fresh designs {n_new} = {len(cohort)} MB backbones x ranks 1-{A.ranks} "
      f"(⛔ rank 1 INCLUDED, unlike s119)")
print(f"   re-draws {len(redraw)}: " + ", ".join(r[0] for r in redraw))
print(f"   + 2 POSITIVE crystals + 2 NONBINDER_crosspaired (the in-batch floor)")
print(f"⛔ s117 copies every non-DESIGN row into EVERY shard ⇒ re-draws and controls fold at all 4 "
      f"seeds. Folds = {n_new} + {(len(redraw) + 4) * 4} = {n_new + (len(redraw) + 4) * 4}")
print("⛔ NOTHING IS LAUNCHED OR UPLOADED BY THIS SCRIPT.")
