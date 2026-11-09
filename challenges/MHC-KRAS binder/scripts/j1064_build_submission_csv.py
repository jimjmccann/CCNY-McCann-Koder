#!/usr/bin/env python3
"""j1064 -- build and validate the Challenge 5 submission CSV (KRAS G12D / HLA-A*11:01 pMHC).

⛔⛔ WHAT THIS SCRIPT REFUSES TO DO: retype a design name, a rank, or a sequence. Everything is
DERIVED and then asserted, because a submitted sequence is the one artefact that cannot be corrected
after the fact, and because design BASENAMES ARE NOT UNIQUE ACROSS ROUNDS. `[M]` 2026-11-09:
`MB_nomotif_broad_11_model_1` exists under BOTH `arm_shard03_` (the 2x2 round) and `arm_shard01_`
(the MB-depth round), and `19_model_0`/`4_model_0` each exist under three different arm/shard
prefixes. ⇒ a stem alone cannot identify a design and a stem-matched `.fa` can be the WRONG MOLECULE.

HOW THE COHORT IS DERIVED (never typed):
  Every design ever taken to 20 draws carries `REDRAW_`-prefixed rows in exactly one of the three
  per-design tables. For each such design we read its OWN rank from the table (⛔ not rank 1 -- a
  re-draw re-folds the sequence that CLEARED, which for this cohort is ranks 1,2,3,4,5), count the
  distinct shards it was folded on (= seeds, must be 4), and count how many cleared the
  pre-registered conjunction. 4/4 => VALIDATED. 1-3/4 => INSIDE NOISE, NOT PROMOTED (`J1058` 6).

HOW THE SEQUENCE IS DERIVED: the same way the FOLDED one was, by reusing `s123`'s own logic --
rank by `overall_confidence` descending, take `recs[rank-1]`, split on ":" and take `parts[0]`.
⛔ The sequence is NEVER read from a model output ([[never-source-target-sequence-from-model-output]]).

usage:
  python3 j1064_build_submission_csv.py --out <dir> [--expect 7]
"""
import argparse
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.path.dirname(HERE), "<WORKDIR>")
TAB = os.path.join(WORK, "tables")
SCR = os.path.join(WORK, "scripts")
S71 = os.path.join(SCR, "s71_make_5b_cofold_inputs.py")
S123 = os.path.join(SCR, "s123_make_mb_depth_cofold_inputs.py")

SCRATCH = os.environ.get("LAB_SCRATCH")
if not SCRATCH:
    sys.exit("⛔ LAB_SCRATCH unset -- source "
             "the environment file that defines the scratch volume (⛔ NOT opt/orca_env.sh, "
             "which does not define it).")

# the two MPNN rounds that can hold a design of this cohort; which one holds a given design is
# RESOLVED BY LOOKUP BELOW, never assumed from the shard prefix.
MPNN_DIRS = [os.path.join(SCRATCH, "c5_2x2_mpnn_v1_20261108"),
             os.path.join(SCRATCH, "c5_mb_depth_mpnn_20261108")]

TABLES = [os.path.join(TAB, "t152_seqdepth_20261108_per_design.tsv"),
          os.path.join(TAB, "t154_mb_depth_20261109_per_design.tsv"),
          os.path.join(TAB, "t155_mb_redraw_20261109_per_design.tsv")]

# Target chain lengths, taken from the CHEMISTRY of the complex rather than from any file we are
# reading, so that a mis-assembled input cannot validate itself.
#
# ⛔⛔ THESE ARE THE DESIGN-STAGE LENGTHS, WHICH ARE NOT THE CO-FOLD ONES. In the design-stage
# sequence files the parts are ordered  parts[0]=binder, parts[1]=HLA (276), parts[2]=b2m (99),
# parts[3]=KRAS G12D peptide (10)  -- i.e. the BINDER IS FIRST. In the co-folded complex the
# lettering rotates to A=HLA, B=b2m, C=peptide, D=binder, and b2m there is 100 residues, not 99.
# ⇒ both numbers are correct in their own frame; see methods/README.md for the convention table.
# `s123` compares against OUR_B2M[1:] because OUR_B2M is the 100-residue form and the design-stage
# files carry it without its leading residue, which is exactly the 99 asserted here.
HLA_LEN, B2M_FASTA_LEN, PEP_LEN = 276, 99, 10
MIN_AA, MAX_AA = 10, 250           # the organisers' own per-chain window

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--expect", type=int, default=7)
ap.add_argument("--seeds", type=int, default=4)
A = ap.parse_args()

# ── reuse s71's pdb_chains so the structural check is the SAME code the fold inputs used ──────────
spec = importlib.util.spec_from_file_location("s71_mod", S71)
s71 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s71)
if not hasattr(s71, "pdb_chains"):
    sys.exit("⛔ s71 lacks `pdb_chains` -- re-audit before trusting any length check")

# ⛔ assert s123 still derives a re-draw the way this script assumes, so the rank logic below is
# not quietly diverging from the round that actually produced the data.
_s123 = open(S123).read()
for _n in ("THE RANK THAT CLEARED, not rank 1", "recs[want_rank - 1]", "parts[0]"):
    if _n not in _s123:
        sys.exit(f"⛔ REFUSING: `{_n}` absent from s123 -- the rank/binder convention this script "
                 f"mirrors has changed; re-read s123 before building a submission.")


def ranked_sequences(fa):
    """VERBATIM from s123. Design records (id=), highest overall_confidence first.
    ⛔ Record 0 of the fasta is the INPUT sequence and carries no `id=`."""
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


def truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes")


# ── 1 · DERIVE the cohort from the REDRAW_ rows of all three tables ───────────────────────────────
cohort = {}          # design -> dict
for path in TABLES:
    if not os.path.exists(path):
        sys.exit(f"⛔ missing table {path} -- refusing to build a partial submission")
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row["tag"].startswith("REDRAW"):
                continue
            d = row["design"]
            rec = cohort.setdefault(d, {"design": d, "arm": row["arm"], "ranks": set(),
                                        "shards": set(), "cleared": 0, "table": os.path.basename(path),
                                        "margins": [], "imargins": []})
            if rec["table"] != os.path.basename(path):
                sys.exit(f"⛔ {d} carries REDRAW rows in TWO tables ({rec['table']} and "
                         f"{os.path.basename(path)}) -- ambiguous provenance, refusing.")
            rec["ranks"].add(row["rank"])
            rec["shards"].add(row["shard"])
            rec["margins"].append(float(row["margin"]))
            rec["imargins"].append(float(row["imargin"]))
            if truthy(row["conj"]):
                rec["cleared"] += 1

if len(cohort) != A.expect:
    sys.exit(f"⛔ derived {len(cohort)} re-drawn designs, expected {A.expect}. ⛔ Refusing to guess "
             f"-- if the cohort really changed, pass --expect explicitly and say why in the note.")

for d, r in cohort.items():
    if len(r["ranks"]) != 1:
        sys.exit(f"⛔ {d}: re-draw rows span ranks {sorted(r['ranks'])} -- a re-draw must re-fold "
                 f"ONE rank. Refusing.")
    r["rank"] = int(next(iter(r["ranks"])))
    if len(r["shards"]) != A.seeds:
        sys.exit(f"⛔ {d}: {len(r['shards'])} distinct shards, expected {A.seeds} seeds. "
                 f"⛔ Never rank or tier on partial data.")
    r["tier"] = "VALIDATED" if r["cleared"] == A.seeds else "INSIDE_NOISE_NOT_PROMOTED"

n_val = sum(1 for r in cohort.values() if r["tier"] == "VALIDATED")
print(f"✅ cohort DERIVED from REDRAW_ rows: {len(cohort)} designs "
      f"({n_val} VALIDATED at {A.seeds}/{A.seeds}, {len(cohort)-n_val} at <{A.seeds}/{A.seeds})")

# ── 2 · RESOLVE each design to exactly one MPNN round by lookup ───────────────────────────────────
for d, r in cohort.items():
    hits = [m for m in MPNN_DIRS if os.path.exists(os.path.join(m, "mpnn", "seqs", d + ".fa"))]
    if len(hits) != 1:
        sys.exit(f"⛔ {d}: found in {len(hits)} MPNN rounds {hits} -- a design must resolve to "
                 f"exactly one. ⛔ THIS IS THE BASENAME-COLLISION CASE; refusing to guess.")
    r["mpnn"] = hits[0]
    r["round"] = os.path.basename(hits[0])

# ── 3 · EXTRACT the binder at its OWN rank, with every check ──────────────────────────────────────
for d, r in sorted(cohort.items()):
    recs = ranked_sequences(os.path.join(r["mpnn"], "mpnn", "seqs", d + ".fa"))
    if len(recs) < r["rank"]:
        sys.exit(f"⛔ {d}: {len(recs)} sequences but it cleared at rank {r['rank']}")
    conf, seq, hdr = recs[r["rank"] - 1]
    parts = seq.split(":")
    order, lens = s71.pdb_chains(os.path.join(r["mpnn"], "pdb", d + ".pdb"))
    if len(parts) != len(order):
        sys.exit(f"⛔ {d}: {len(parts)} fasta parts vs {len(order)} pdb chains")
    if [len(p) for p in parts] != [lens[c] for c in order]:
        sys.exit(f"⛔ {d}: part lengths {[len(p) for p in parts]} != pdb chains "
                 f"{[lens[c] for c in order]}")
    binder = parts[0]
    # chemistry-level validation of the TARGET limbs, independent of any file we parsed
    if len(parts) != 4:
        sys.exit(f"⛔ {d}: expected 4 chains (binder + HLA + b2m + peptide), got {len(parts)}")
    for idx, want, what in ((1, HLA_LEN, "HLA"), (2, B2M_FASTA_LEN, "b2m"), (3, PEP_LEN, "peptide")):
        if len(parts[idx]) != want:
            sys.exit(f"⛔ {d}: {what} chain is {len(parts[idx])} aa, chemistry says {want}")
    if "C" in binder:
        sys.exit(f"⛔ {d}: CYSTEINE in a competition design -- s53 should have run --omit_AA C")
    if not (MIN_AA <= len(binder) <= MAX_AA):
        sys.exit(f"⛔ {d}: binder {len(binder)} aa outside the organisers' {MIN_AA}-{MAX_AA} window")
    if not binder.isalpha() or not binder.isupper():
        sys.exit(f"⛔ {d}: binder is not a plain uppercase sequence")
    r.update(binder=binder, conf=conf, binder_len=len(binder),
             target=":".join(parts[1:]), best_margin=max(r["margins"]))

# ⛔ all seven must have been folded against the SAME target, and the binders must be distinct
targets = {r["target"] for r in cohort.values()}
if len(targets) != 1:
    sys.exit(f"⛔ the cohort carries {len(targets)} DIFFERENT target sequences -- these designs were "
             f"not all folded against the same complex. Refusing.")
seqs = {}
for d, r in cohort.items():
    if r["binder"] in seqs:
        sys.exit(f"⛔ {d} and {seqs[r['binder']]} carry the IDENTICAL binder sequence -- a duplicate "
                 f"row would waste a submission slot. Refusing.")
    seqs[r["binder"]] = d
print(f"✅ one shared target across all {len(cohort)} designs; {len(seqs)} distinct binder sequences")

# ── 4 · NAME and RANK the rows ────────────────────────────────────────────────────────────────────
# public-safe, traceable name: target-arm-backbone-rank. ⛔ no internal paths, no person's name.
def public_name(r):
    arm = "mb" if r["arm"].startswith("MB") else ("mn" if r["arm"].startswith("M") else "x")
    stem = r["design"].split("_")
    bb, model = stem[-3], stem[-1]
    return f"kras-mhc-{arm}-{bb}-m{model}-r{r['rank']}"

rows = sorted(cohort.values(),
              key=lambda r: (0 if r["tier"] == "VALIDATED" else 1, -r["best_margin"]))
for r in rows:
    r["name"] = public_name(r)
if len({r["name"] for r in rows}) != len(rows):
    sys.exit("⛔ public names are not unique -- refusing")

os.makedirs(A.out, exist_ok=True)
csv_path = os.path.join(A.out, "CCNY-McCann-Koder_Challenge_5_PROPOSED_20261109.csv")
with open(csv_path, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["name", "sequence", "molecule_class"])
    for r in rows:
        w.writerow([r["name"], r["binder"], "single_chain"])

prov_path = os.path.join(A.out, "j1064_submission_provenance_20261109.tsv")
with open(prov_path, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["row", "name", "design", "arm", "mpnn_round", "mpnn_rank", "mpnn_confidence",
                "binder_len", "seeds_cleared", "seeds_total", "tier", "best_margin"])
    for i, r in enumerate(rows, start=1):
        w.writerow([i, r["name"], r["design"], r["arm"], r["round"], r["rank"],
                    f"{r['conf']:.4f}", r["binder_len"], r["cleared"], A.seeds, r["tier"],
                    f"{r['best_margin']:+.1f}"])

print(f"\n✅ wrote {csv_path}")
print(f"✅ wrote {prov_path}")
print(f"\n{'row':>3}  {'name':<26} {'len':>4} {'tier':<26} {'rank':>4}  design")
for i, r in enumerate(rows, start=1):
    print(f"{i:>3}  {r['name']:<26} {r['binder_len']:>4} {r['tier']:<26} {r['rank']:>4}  "
          f"{r['design']}")
print(f"\n⛔ {len(rows)} rows against the organisers' 20-row cap; every chain inside "
      f"{MIN_AA}-{MAX_AA} aa; molecule_class=single_chain.")
print("⛔ This script does not upload anything. Novelty is scored AFTER upload and the submission "
      "slot is once per 24 hours, so the file is validated here before it is ever sent.")
