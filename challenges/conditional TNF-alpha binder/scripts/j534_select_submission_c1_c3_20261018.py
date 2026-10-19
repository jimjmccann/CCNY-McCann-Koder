#!/usr/bin/env python3
"""j534 -- select the TNFα submission 20 across BOTH arms: C1 + C3, with C3's slots RESERVED.

⛔⛔ WHY NOT `j507`. `j507` hardcodes the **v46 joint11** cohort
(`j496_v46_joint11_fold_20261018/j497_ALL110_...tsv`). The submission is drawn from the **REBIAS**
cohorts (`j513_c1_rebias` / `j514_c3_rebias`) -- that mismatch is `J513`'s highest-risk finding.
Pointing `j507` at a new TSV would also not give it a second arm. Build alongside (`J512`/`J517`
precedent); ⛔ `j507` is untouched and remains the record of the v46 selection.

## THE SELECTION RULES THIS FILE IMPLEMENTS
- The split between the two arms is an argument (`--n-c1` / `--n-c3`); ⛔ nothing is
  hardcoded to a particular split.
- ⭐⭐ **The C3 arm is represented regardless of score** ⇒
  **THE C3 SLOTS ARE RESERVED. SCORE NEVER TAKES ONE.** C1 backfills **only** a slot C3 cannot
  fill on a HARD failure (missing fold · cysteine · length · non-standard residue · duplicate).
  ⛔ Never because a C3 design merely ranks low.

## ⭐ WHAT RANKS THESE -- and what it is NOT
⛔⛔ **This is NOT a binding prediction.** `J497` §9a-FINAL: on this target backbone identity
predicts its own fold and nothing about binding; `J502`/`J505` removed both candidate gates
(`iptm` -- a non-binder outscored our median; epitope enrichment -- the ordering REVERSED).
`J518`: the Baker five-criteria verdict is **0/88**, so ⛔ a binary gate returns all-zeros and
discriminates nothing. We rank on **continuous**
**designed-geometry recovery** -- the one thing we actually specified.

- **C1**: `min(pred_pin_glu107_{B,C,D})` -- the PINNED guidepost ARG to GLU107 at its best-engaged
  protomer. ⛔ The *pinned* Arg, not any Arg: "some ARG reaches GLU107" is satisfiable by accident
  on an MPNN surface. ⛔⛔ **NOT `dbridge_pin_*`** -- that is `pred - design`, a DRIFT not a
  distance, and ranking on it picks designs whose bridge merely *moved least* (`j507`, caught 10-18).
- **C3**: ⛔⛔ **`bridge_copy_best_BEST`, NOT `_worst`** -- see the inversion note below. It is the
  **conjunction**: the WORST copy's BEST bridge. A mean pays partial credit.
  ⚠️ **Watch this one against the data**: a C3 binder may grip at most 2 of 3 protomers, in which case the third copy dominates
  the conjunction and the column stops discriminating. `--c3-rank-col` exists so that is a
  measurement, not a rewrite. ⛔ Look before trusting it.

## ⛔⛔ `j521`'s `_best`/`_worst` LABELS ARE INVERTED FOR EVERY DISTANCE COLUMN (`J535`, 10-19)
`j521.agg()` hardcodes `_best = np.max` and `_worst = np.min`. That is right for pLDDT / TM / ipTM
(higher is better) and ⛔ **BACKWARDS for `bridge_copy_best_*`, `rmsd_binder_copy_*` and
`pae_*`**, where lower is better. verified on `j535`: `_best > _worst` in **96/96** rows of
all four families, so the aggregator is plainly max/min throughout.
⇒ for a DISTANCE, the conjunction (the worst copy) is **`_best`**:
  - `bridge_copy_best_best`  = max over copies of (that copy's min bridge) = ⭐ **the worst copy's
    best bridge = THE CONJUNCTION.** Rank ASCENDING on it.
  - `bridge_copy_best_worst` = min over copies = ⛔ the single MOST-engaged copy = pure partial
    credit, **worse than the mean**, and the exact failure the rule above exists to prevent.
⛔ `j521` is NOT edited: the numbers are correct and only the names mislead, `j533`/`j535` are
already written against those names, and `--mode holo` has yet to run (`J512`/`J517` precedent --
never edit the working scorer mid-deadline). The fix belongs in the CONSUMER, which is this file.
- **Tiebreak, both arms**: self-consistency RMSD. ⛔ Applied as a CONJUNCTION (bridge first, then
  scRMSD), never a mean of normalised fractions.

## HARD CONSTRAINTS -- all asserted BEFORE anything is written
- ⛔ `--n-c1 + --n-c3` ≤ **20**; each chain **10-250 aa**.
- ⛔ **No cysteine** in the DESIGNED chain.
- ⛔ **Backbone diversity cap** (`--per-backbone`): Adaptyv's `batch_duplicate` fails our own rows
  against EACH OTHER at ≥0.95 identity / ≥0.80 coverage, and C1 draws on only **11** backbones.
- ⛔ **C3 homotrimer check**: binder chains B and C must be IDENTICAL to A, or it is not the
  homotrimer we designed and one chain is not a valid submission for it.
- ⛔ The authoritative prior-art verdict is Adaptyv's own prescreen, run separately on the FASTA
  this writes. Our bar is **0.70** (non-antibody, strict arm), ≥8 aa verbatim fragment fails.
"""
import argparse
import csv
import os
import sys

LO, HI = 10, 250
STD = set("ACDEFGHIKLMNPQRSTVWY")


def fnum(row, key):
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return None


def read_fasta_chains(path):
    cur, seqs = None, {}
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split("|")[0].strip()
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line
    return seqs


def c1_key(r):
    vals = [v for v in (fnum(r, f"pred_pin_glu107_{c}") for c in "BCD") if v is not None]
    bridge = min(vals) if vals else 9e9
    sc = fnum(r, "rmsd_binder")
    return (bridge, sc if sc is not None else 9e9)


def make_c3_key(col, tiebreak):
    def key(r):
        b = fnum(r, col)
        sc = fnum(r, tiebreak)
        return (b if b is not None else 9e9, sc if sc is not None else 9e9)
    return key


def vet(name, seq, arm, problems):
    """Return True if this sequence may be submitted; record WHY NOT otherwise."""
    if not seq:
        problems.append(f"{arm} {name}: no sequence"); return False
    if not (LO <= len(seq) <= HI):
        problems.append(f"{arm} {name}: {len(seq)} aa outside {LO}-{HI}"); return False
    bad = sorted(set(seq) - STD)
    if bad:
        problems.append(f"{arm} {name}: non-standard residue(s) {bad}"); return False
    if "C" in seq:
        problems.append(f"{arm} {name}: CYSTEINE -- breaks expression"); return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c1-tsv", required=True)
    ap.add_argument("--c1-inputs", required=True)
    ap.add_argument("--c3-tsv", default="")
    ap.add_argument("--c3-inputs", default="")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--n-c1", type=int, default=10)
    ap.add_argument("--n-c3", type=int, default=10)
    ap.add_argument("--per-backbone", type=int, default=2)
    # ⛔⛔ `_best` IS THE CONJUNCTION FOR A DISTANCE -- j521.agg() is max/min regardless of
    # direction, so for bridge/rmsd/pae the WORST copy lands in `_best`. See the module docstring.
    ap.add_argument("--c3-rank-col", default="bridge_copy_best_best")
    ap.add_argument("--c3-tiebreak-col", default="rmsd_binder_copy_best")
    a = ap.parse_args()

    # ⛔⛔ Tripwire for the inversion. A DISTANCE column named `_worst` is j521's np.min, i.e. the
    # single BEST copy -- ranking on it pays partial credit. Caught once (10-19); refuse, never warn.
    for flag, col in (("--c3-rank-col", a.c3_rank_col), ("--c3-tiebreak-col", a.c3_tiebreak_col)):
        if col.endswith("_worst") and any(k in col for k in ("bridge", "rmsd", "pae", "dist")):
            sys.exit(f"⛔⛔ {flag}={col}: j521's `_worst` is np.min, so for a DISTANCE it is the "
                     f"MOST-engaged copy, not the least. Use {col[:-6]}_best for the conjunction.")

    total = a.n_c1 + a.n_c3
    if total > 20:
        sys.exit(f"⛔ {a.n_c1}+{a.n_c3}={total} exceeds the organiser's 20-protein cap")
    for p in (a.out_tsv, a.out_fasta):
        if os.path.exists(p):
            sys.exit(f"⛔ exists, refusing to overwrite: {p}  (use a fresh --out for each run)")

    problems, picked = [], []

    def take(tsv, inputs, want, arm, keyfn, chains_expected):
        """Rank an arm and take `want` rows that pass every hard filter."""
        if not tsv:
            return []
        rows = list(csv.DictReader(open(tsv), delimiter="\t"))
        ranked = sorted(rows, key=keyfn)
        out, per_bb = [], {}
        for r in ranked:
            if len(out) >= want:
                break
            name = r.get("design", "")
            bb = r.get("backbone", "")
            if per_bb.get(bb, 0) >= a.per_backbone:
                continue
            fp = os.path.join(inputs, name + ".fasta")
            if not os.path.exists(fp):
                problems.append(f"{arm} {name}: missing input FASTA {fp}")
                continue
            seqs = read_fasta_chains(fp)
            if "A" not in seqs:
                problems.append(f"{arm} {name}: no chain A in {fp}")
                continue
            seq = seqs["A"].upper()
            # ⛔ C3 must be a genuine homotrimer, or one chain does not represent it
            if len(chains_expected) > 1:
                bad = [c for c in chains_expected[1:]
                       if seqs.get(c, "").upper() != seq]
                if bad:
                    problems.append(f"{arm} {name}: binder chains {bad} differ from A "
                                    f"-- NOT the homotrimer that was designed")
                    continue
            if not vet(name, seq, arm, problems):
                continue
            out.append(dict(name=name, sequence=seq, arm=arm, backbone=bb,
                            rank_value=keyfn(r)[0], tiebreak=keyfn(r)[1]))
            per_bb[bb] = per_bb.get(bb, 0) + 1
        return out

    # ⭐⭐ C3 FIRST, because its slots are RESERVED and must not be competed for.
    c3 = take(a.c3_tsv, a.c3_inputs, a.n_c3, "C3",
              make_c3_key(a.c3_rank_col, a.c3_tiebreak_col), ["A", "B", "C"])
    shortfall = a.n_c3 - len(c3)
    c1_want = a.n_c1 + shortfall
    c1 = take(a.c1_tsv, a.c1_inputs, c1_want, "C1", c1_key, ["A"])

    picked = c3 + c1
    # dedupe by sequence across arms, keeping the earlier (C3-first) entry
    seen, final = set(), []
    for p in picked:
        if p["sequence"] in seen:
            problems.append(f"{p['arm']} {p['name']}: duplicate sequence, dropped")
            continue
        seen.add(p["sequence"])
        final.append(p)

    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["name", "sequence", "arm", "backbone",
                                       "rank_value", "tiebreak"])
        w.writeheader()
        w.writerows(final)
    with open(a.out_fasta, "w") as fh:
        for p in final:
            fh.write(f">{p['name']}\n{p['sequence']}\n")

    n_c3 = sum(1 for p in final if p["arm"] == "C3")
    n_c1 = sum(1 for p in final if p["arm"] == "C1")
    print(f"wrote {a.out_tsv}  and  {a.out_fasta}")
    print(f"  C3 {n_c3}/{a.n_c3} (RESERVED)   C1 {n_c1}   total {len(final)}/{total}")
    if shortfall > 0 and not a.c3_tsv:
        # ⛔ Distinguish "no C3 data yet" from "C3 designs failed". Reporting a missing
        # input as a hard failure of the designs would misread a pipeline state as a result.
        print(f"  ⛔⛔ NO --c3-tsv GIVEN: all {shortfall} reserved C3 slots went to C1. "
              f"This is NOT a C3 result -- it is a MISSING INPUT, and the output of this "
              f"run is not a valid two-arm submission.")
    elif shortfall > 0:
        print(f"  ⚠️ C3 filled only {a.n_c3 - shortfall}/{a.n_c3}; {shortfall} slot(s) lost "
              f"to HARD failures (never to score) -- C1 backfilled. Reasons listed below.")
    if len(final) < total:
        print(f"  ⛔ {total - len(final)} slots UNFILLED -- an unused slot is discarded option "
              f"value. Investigate before uploading.")
    lens = [len(p["sequence"]) for p in final]
    if lens:
        print(f"  lengths {min(lens)}-{max(lens)} aa (bound {LO}-{HI}), "
              f"{len(set(p['sequence'] for p in final))} distinct sequences, 0 cysteines")
    if problems:
        print(f"\n  {len(problems)} rows rejected or flagged (first 15):")
        for p in problems[:15]:
            print(f"    {p}")
    print("\n⛔ Prior art is NOT checked here -- run Adaptyv's prescreen on the FASTA.")
    print("⛔⛔ The organiser accepts ONE submission per 24 h -- there is no same-day replacement.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
