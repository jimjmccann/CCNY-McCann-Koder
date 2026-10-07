#!/usr/bin/env python3
"""Build a SINGLE-POSITION SATURATION SCAN cycle for egfrA1. NEW FILE 2026-10-05.

⭐⭐ WHY THIS EXISTS, and why it is not the combo builder. `build_combo_cycle_egfrA1.py`
enumerates COMBINATIONS of the cycle-1 accepted set (L16D/A17M/T19Q/L49V/L49K/R55S/S75A).
That walk is **COMPLETE, not stalled** -- 6 positions x 7 accepted mutations = 88 combinations,
cycle 2 took 50 and cycle 3 took the remaining 38, there is no order 7, and the real proposer
returns 0 new at every order. ⛔ So re-running the combo builder
for c4 produces nothing; it is not a bug to fix.

⛔⛔ AND THE ACCEPTED SET CANNOT REACH THE ONE THING LEFT TO OPTIMISE. Under `decisions/0047`
the site-A objective is `min(n_epitope/11, c380)` per draw, and F380 is the binding term:
measured 2026-10-05 on the 98-sequence c2+c3+c90 cohort, `label_A_0047` equals the plain F380
contact rate EXACTLY for 63 of 98 sequences -- and for ALL TEN of the re-measured top designs.
Among designs that fold reliably the objective IS the F380 rate; the other ten footprint
residues are already satisfied whenever F380 is. ⭐ **Binder position 66 is effectively F380's
only handle** (within 5 A of it in 248 of the 310 contacting draws, 80%; position 70 is second
at 7.7%), and NONE of L16/A17/T19/L49/R55/S75 reaches it. 18 of the 19 substitutions at 66 have
never been folded.
⇒ Hence a saturation scan of ONE position, which is exhaustively enumerable (19 candidates) and
therefore needs no surrogate. ⛔ That matters: UCB did NOT beat random on this arm (+0.394 +/-
0.634, t=+0.62; Fisher p=1.000) and OOB RMSE/sd was 1.04-1.14, i.e. at or worse than predicting
the mean. The model is OUT OF THE LOOP by measurement, not by preference.

⭐ THE 0045 REPLICATE RIDE-ALONG IS BUILT IN, because it is standing policy: the top-scoring
designs are re-folded in every cycle. The
scan emits a sibling `<cycle>_replicates.csv` in exactly the format `score_cycle_egfrA1.py`
already reads, so the re-folds are scored, compared against their priors, and ⛔ NEVER merged
into the training table. The `id41` ANCHOR is included by default: it is a known-prior design
low on the scale, which turns the test-retest from a yes/no into a correlation.

⛔ GUARD RAILS, each of which has previously cost either wasted compute or a wrong conclusion:
  - The wild-type residue at the scanned position is ASSERTED against the parent before anything
    is generated. An off-by-one would produce 19 plausible-looking wrong sequences.
  - Every emitted sequence must differ from the parent in EXACTLY one position.
  - The parent's own residue is never emitted as a "substitution" (that would be a silent
    duplicate of the parent, and a free draw spent on nothing).
  - `--exclude` drops substitutions already folded FROM THIS PARENT. ⛔ Do NOT use it to drop
    one folded from a DIFFERENT reference: Q66G was measured from the id41 hub (fef 0.40 -> 0.00
    while folding BETTER than the hub, median_sc 0.587 vs 0.612), which says nothing about Q66G
    from a cycle-2 champion. Re-folding it from a new parent is a control, not a waste.

usage:
  build_sitescan_cycle_egfrA1.py <parent.txt> <out.tsv> <out_map.csv>
      --position 66 [--idbase 4001] [--design <name>] [--exclude G,P]
      [--replicate-from <scored_0047.csv> --replicate-top 5 --replicate-idbase 4901]
      [--anchor id41] [--no-anchor]
"""
import argparse
import csv
import sys

AA = "ACDEFGHIKLMNPQRSTVWY"      # 20; the parent's own residue is dropped, leaving 19


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parent"); ap.add_argument("out_tsv"); ap.add_argument("out_map")
    ap.add_argument("--position", type=int, required=True, help="1-based position in the BINDER")
    ap.add_argument("--idbase", type=int, default=4001)
    ap.add_argument("--design", default="arm_shard11_shard11_24_model_2")
    ap.add_argument("--exclude", default="",
                    help="substitution residues to skip, e.g. G,P -- only if already folded "
                         "FROM THIS PARENT")
    ap.add_argument("--expect-wt", default="",
                    help="optional: the residue you BELIEVE is at --position. If given it is "
                         "asserted too, so a wrong belief fails loudly instead of silently "
                         "scanning the wrong site.")
    ap.add_argument("--replicate-from", help="a scored 0047 CSV to pick re-folds from")
    ap.add_argument("--replicate-top", type=int, default=5)
    ap.add_argument("--replicate-idbase", type=int, default=4901)
    ap.add_argument("--anchor", default="id41",
                    help="a known-prior LOW design, so test-retest yields a correlation")
    ap.add_argument("--no-anchor", action="store_true")
    a = ap.parse_args()

    parent = open(a.parent).read().split()[0]
    pos = a.position
    if not 1 <= pos <= len(parent):
        sys.exit(f"⛔ position {pos} outside 1..{len(parent)} -- refusing")
    wt = parent[pos - 1]
    print(f"parent {len(parent)} aa; position {pos} is {wt!r}")
    if a.expect_wt and a.expect_wt.upper() != wt:
        sys.exit(f"⛔ --expect-wt {a.expect_wt.upper()} but the parent has {wt!r} at position "
                 f"{pos}. REFUSING -- one of the two is wrong and guessing corrupts every "
                 f"sequence in the cycle.")
    drop = {c.strip().upper() for c in a.exclude.split(",") if c.strip()}
    if drop:
        print(f"  --exclude: skipping {sorted(drop)} "
              f"(⛔ only valid if already folded FROM THIS PARENT)")

    subs = [c for c in AA if c != wt and c not in drop]
    print(f"  scanning {len(subs)} substitutions: {''.join(subs)}")

    rows = []
    for i, s in enumerate(subs):
        seq = parent[:pos - 1] + s + parent[pos:]
        rows.append(dict(fold_id=f"id{a.idbase + i}", order=1, n_mut=1,
                         mutations=f"{wt}{pos}{s}", positions=str(pos), sequence=seq))

    # ---- the guards ----
    for r in rows:
        d = [k for k, (x, y) in enumerate(zip(parent, r["sequence"]), 1) if x != y]
        if d != [pos]:
            sys.exit(f"⛔ {r['fold_id']} {r['mutations']}: differs from the parent at {d}, "
                     f"expected exactly [{pos}]. REFUSING.")
    if len({r["sequence"] for r in rows}) != len(rows):
        sys.exit("⛔ duplicate sequences")
    if len({r["fold_id"] for r in rows}) != len(rows):
        sys.exit("⛔ duplicate fold ids")
    if parent in {r["sequence"] for r in rows}:
        sys.exit("⛔ one emitted sequence IS the parent -- refusing to spend a draw on nothing")
    print(f"✅ {len(rows)} sequences, all distinct, each differing from the parent at "
          f"exactly position {pos}, none equal to the parent")

    # ---- decisions/0045: the replicate ride-along ----
    reps = []
    if a.replicate_from:
        scored = list(csv.DictReader(open(a.replicate_from)))
        if "label_A_0047" not in (scored[0].keys() if scored else {}):
            sys.exit(f"⛔ {a.replicate_from} has no label_A_0047 column -- it was not written by "
                     f"score_0047.py. Refusing rather than ranking on a metric this script cannot name.")
        # ⛔ rank on the 0047 objective, NOT on fef. They rank differently (Kendall tau +0.861).
        best, seen = [], set()
        for r in sorted(scored, key=lambda r: -float(r["label_A_0047"])):
            if r["seq_id"] in seen:
                continue
            seen.add(r["seq_id"]); best.append(r)
        pick = best[:a.replicate_top]
        if not a.no_anchor:
            anc = next((r for r in scored if r["seq_id"] == a.anchor), None)
            if anc is None:
                print(f"  ⚠️ anchor {a.anchor} is NOT in {a.replicate_from} -- it has no prior "
                      f"there, so the ride-along goes out WITHOUT an anchor and the test-retest "
                      f"will be a yes/no, not a correlation. ⛔ Not inventing a prior.")
            elif anc["seq_id"] not in {r["seq_id"] for r in pick}:
                pick.append(anc)
                print(f"  ⭐ anchor {a.anchor} added (0047 {anc['label_A_0047']}) -- a LOW known "
                      f"prior is what makes the retest a correlation")
        for i, r in enumerate(pick):
            reps.append(dict(replicate_id=f"id{a.replicate_idbase + i}",
                             source_id=r["seq_id"],
                             prior_fef=r["fef_legacy"],
                             prior_label=r["label_legacy"],
                             prior_median_sc=r["median_sc"],
                             prior_n_draws=r["n_draws"],
                             prior_label_A_0047=r["label_A_0047"],
                             mutations=f"REFOLD_of_{r['seq_id']}",
                             sequence=r["sequence"]))
        missing = [r for r in reps if not r["sequence"]]
        if missing:
            sys.exit(f"⛔ {len(missing)} replicates have NO sequence in the scored table "
                     f"({[r['source_id'] for r in missing]}) -- re-run score_0047.py with the "
                     f"matching --tsv so sequences are carried through. Refusing.")
        print(f"  decisions/0045 ride-along: {len(reps)} re-folds "
              f"({', '.join(r['source_id'] for r in reps)})")

    # ---- write ----
    with open(a.out_tsv, "w") as fh:
        for r in rows + reps:
            fid = r.get("fold_id") or r["replicate_id"]
            fh.write(f"A_{a.design}_{fid}\t{r['sequence']}\n")
    with open(a.out_map, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fold_id", "order", "n_mut", "mutations", "positions"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in w.fieldnames})
    print(f"wrote {a.out_tsv} ({len(rows)+len(reps)} fold entries) and {a.out_map} "
          f"({len(rows)} scan rows)")
    if reps:
        # ⛔ the name is a CONVENTION score_cycle_egfrA1.py derives from the map path. Do not rename.
        if not a.out_map.endswith("_map.csv"):
            sys.exit(f"⛔ {a.out_map} must end in _map.csv -- score_cycle_egfrA1.py finds the "
                     f"replicates file by replacing that suffix, and a re-fold it cannot find "
                     f"gets MERGED into the training table, which is the corruption 0045 exists "
                     f"to prevent.")
        q = a.out_map[:-len("_map.csv")] + "_replicates.csv"
        with open(q, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[k for k in reps[0] if k != "sequence"])
            w.writeheader()
            for r in reps:
                w.writerow({k: v for k, v in r.items() if k != "sequence"})
        print(f"wrote {q} ({len(reps)} re-folds, held OUT of the training merge by 0045)")
    print(f"\n⇒ cycle total: {len(rows)} scan + {len(reps)} replicate = "
          f"{len(rows)+len(reps)} sequences")


if __name__ == "__main__":
    main()
