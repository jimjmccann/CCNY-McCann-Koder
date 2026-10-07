#!/usr/bin/env python3
"""Answer the egfrA1 cycle-3 FORK from the data: do the beneficial singles STACK, or not?

⭐ WHY THIS EXISTS. The campaign plan made cycle 3 conditional -- *"Only launch it if cycle 2 shows
stacking HELPS; if triples are worse than doubles, the walk should go back to singles from the best
double instead. That is a real fork -- read the data."* The proposer enumerates orders 4/5/6
MECHANICALLY and cannot see that triples underperformed, so the fork has to be decided by a human
reading numbers, which is what this script lays out.

⛔ THIS SCRIPT DOES NOT DECIDE. It prints the comparison and the criteria; a human picks the branch.
   That is deliberate -- the decision is a human one, and an automated verdict would quietly convert
   a science fork into a default.

⛔ CONVENTIONS ARE COPIED, NOT REINVENTED, from `score_cycle_egfrA1.py`:
     GATE  = sc_rmsd < 1.5   (a draw that folded)
     FLOOR = gate_frac >= 0.5  -- a CANDIDATE FILTER, never summed into the label
     FULL  = 11 epitope residues = the label's physical CEILING
⛔ READ `full_epitope_frac`, NOT `label`. The label SATURATES: cycle 1's top three sat at identical
   means (dmean +0.00, Fisher p=0.628-1.000), so it can no longer rank. Both are printed, with the
   label marked as the saturated one, because a reader who sees only one metric cannot spot that.

⭐ ORDER IS DERIVED, NOT READ FROM A MAP FILE: order = Hamming distance from the HUB (the row with
   the most Hamming-1 neighbours -- `id41`, 50 vs 4 for the runner-up). Same rule the proposer uses,
   so the two can never disagree about what "a double" means. ⛔ Refuses a thin hub margin.

usage: analyze_stacking_egfrA1.py <training_cN.csv> [--min-hub-margin 3]
"""
import argparse
import collections
import statistics
import sys
from math import comb

GATE, FLOOR, FULL = 1.5, 0.5, 11
MIN_HUB_NEIGHBOURS = 5


def die(m):
    sys.exit(f"⛔ {m}")


def truthy(v):
    return str(v).strip() in ("1", "True", "true")


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def welch(xs, ys):
    """mean difference, SE, t -- no scipy in this env (same helper as the scorer)."""
    if len(xs) < 2 or len(ys) < 2:
        return float("nan"), float("nan"), float("nan")
    mx, my = statistics.mean(xs), statistics.mean(ys)
    se = (statistics.variance(xs) / len(xs) + statistics.variance(ys) / len(ys)) ** 0.5
    return mx - my, se, ((mx - my) / se if se else float("nan"))


def fisher(a, b, c, d):
    n = a + b + c + d
    r1, c1 = a + b, a + c
    if min(n, r1, c1) < 0 or n == 0:
        return float("nan")

    def p(x):
        return comb(r1, x) * comb(n - r1, c1 - x) / comb(n, c1)
    o = p(a)
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    return sum(p(x) for x in range(lo, hi + 1) if p(x) <= o + 1e-12)


def main():
    import csv
    ap = argparse.ArgumentParser()
    ap.add_argument("training")
    ap.add_argument("--min-hub-margin", type=int, default=3)
    a = ap.parse_args()

    rows = [r for r in csv.DictReader(open(a.training)) if r.get("sequence")]
    if not rows:
        die(f"no rows with sequences in {a.training}")
    need = ("seq_id", "sequence", "label", "full_epitope_frac", "passes_floor",
            "n_gated", "n_draws", "median_sc")
    miss = [c for c in need if c not in rows[0]]
    if miss:
        die(f"{a.training} lacks columns {miss}")

    lens = {len(r["sequence"]) for r in rows}
    if len(lens) != 1:
        die(f"mixed sequence lengths {sorted(lens)}")

    # ---- the hub, same rule as the proposer ----
    counts = sorted(((sum(1 for q in rows if q is not r
                          and hamming(r["sequence"], q["sequence"]) == 1), r["seq_id"], r)
                     for r in rows), key=lambda t: (-t[0], t[1]))
    n1, id1, hub = counts[0]
    n2, id2 = (counts[1][0], counts[1][1]) if len(counts) > 1 else (0, "-")
    if n1 < MIN_HUB_NEIGHBOURS:
        die(f"hub {id1} has only {n1} single mutants; no measured walk to read")
    if n1 - n2 < a.min_hub_margin:
        die(f"hub AMBIGUOUS: {id1} ({n1}) vs {id2} ({n2}); refusing to pick a reference")
    H = hub["sequence"]
    print(f"reference (hub): {id1}  [{n1} Hamming-1 neighbours vs {n2} for {id2}]")
    print(f"  {id1}: label {hub['label']}  full_epitope_frac {hub['full_epitope_frac']}  "
          f"gate {hub['n_gated']}/{hub['n_draws']}  median_sc {hub['median_sc']}")

    # ---- the ACCEPTED SET, same rule as the proposer (beats the hub on BOTH metrics) ----
    hl, hf = float(hub["label"]), float(hub["full_epitope_frac"])
    accepted = set()
    for r in rows:
        if r["seq_id"] == id1 or r["label"] == "" or not truthy(r["passes_floor"]):
            continue
        d = [(i + 1, H[i], r["sequence"][i])
             for i in range(len(H)) if H[i] != r["sequence"][i]]
        if len(d) == 1 and float(r["label"]) > hl and float(r["full_epitope_frac"]) > hf:
            accepted.add((d[0][0], d[0][2]))          # (1-based position, substituted residue)
    print(f"accepted set: {len(accepted)} substitutions over "
          f"{len({p for p, _ in accepted})} positions "
          f"{sorted({p for p, _ in accepted})}")
    if not accepted:
        die("accepted set is empty -- nothing beat the reference on both metrics, so there is no "
            "walk to read. That is a RESULT needing a decision, not a cycle.")

    # ---- bucket by ORDER, where order means "a combination of k ACCEPTED substitutions" ----
    # ⛔⛔ RAW HAMMING DISTANCE IS NOT THE ORDER, and using it produced nonsense. Found by running
    #    this on the real `training_c2.csv`: that table also holds **cycle 0's burn-in designs**,
    #    which are NOT mutants of the hub at all -- unrelated sequences sitting 15-31 substitutions
    #    away. Bucketing on distance alone invented "orders 15..31", each with 1-12 members, and the
    #    adjacent-order comparison then walked them and printed `nan`.
    # ⇒ A row counts as order k ONLY IF every one of its k differences from the hub is an ACCEPTED
    #    substitution. Anything else is OFF-WALK and is reported separately, never bucketed as a
    #    high order. This is the proposer's own definition, so the two cannot disagree about what
    #    "a double" means.
    byord = collections.defaultdict(list)
    offwalk = []
    for r in rows:
        if r["seq_id"] == id1:
            continue
        d = [(i + 1, r["sequence"][i])
             for i in range(len(H)) if H[i] != r["sequence"][i]]
        if d and all(x in accepted for x in d) and len({p for p, _ in d}) == len(d):
            byord[len(d)].append(r)
        else:
            offwalk.append(r)
    print(f"on the walk: {sum(len(v) for v in byord.values())} sequences; "
          f"OFF-WALK (not combinations of the accepted set -- e.g. cycle-0 burn-in): "
          f"{len(offwalk)}, excluded from the comparison below")

    def fef(r):
        return float(r["full_epitope_frac"])

    print(f"\n{'order':>5} {'n':>4} {'folds':>7} {'fold-fail':>9} "
          f"{'median fef':>11} {'mean fef':>9} {'best fef':>9} {'mean label*':>12}")
    print("  " + "-" * 76)
    summary = {}
    for o in sorted(byord):
        v = byord[o]
        ok = [r for r in v if truthy(r["passes_floor"])]
        dead = [r for r in v if not truthy(r["passes_floor"])]
        if not ok:
            print(f"{o:>5} {len(v):>4} {0:>7} {len(dead):>9}   -- every one FAILS the fold floor --")
            summary[o] = None
            continue
        f = [fef(r) for r in ok]
        lab = [float(r["label"]) for r in ok if r["label"] != ""]
        summary[o] = dict(n=len(v), ok=ok, f=f,
                          med=statistics.median(f), mean=statistics.mean(f),
                          best=max(ok, key=fef), dead=len(dead),
                          meanlab=statistics.mean(lab) if lab else float("nan"))
        s = summary[o]
        print(f"{o:>5} {len(v):>4} {len(ok):>7} {len(dead):>9} "
              f"{s['med']:>11.3f} {s['mean']:>9.3f} {fef(s['best']):>9.3f} {s['meanlab']:>12.3f}")
    print("  * `label` SATURATES at 11 and cannot rank; it is shown only so a flat column is visible.")
    print("  `folds` = passes the >=50% fold floor;  `fold-fail` = the stacking COST, if any.")

    # ---- the fork ----
    print("\n=== THE FORK: does stacking help? ===")
    orders = [o for o in sorted(summary) if summary[o]]
    if len(orders) < 2:
        die("fewer than two orders have any floor-passing sequence; nothing to compare")
    for lo, hi in zip(orders, orders[1:]):
        if len(summary[hi]["f"]) < 2 or len(summary[lo]["f"]) < 2:
            print(f"  order {hi} vs {lo}:  ⛔ NOT COMPARABLE -- "
                  f"n={len(summary[hi]['f'])} vs {len(summary[lo]['f'])} floor-passing. "
                  f"A one-sequence 'order' has no variance; refusing to print a statistic for it.")
            continue
        d, se, t = welch(summary[hi]["f"], summary[lo]["f"])
        a_hi = sum(1 for x in summary[hi]["f"] if x > 0)
        a_lo = sum(1 for x in summary[lo]["f"] if x > 0)
        p = fisher(a_hi, len(summary[hi]["f"]) - a_hi, a_lo, len(summary[lo]["f"]) - a_lo)
        print(f"  order {hi} vs {lo}:  mean fef {d:+.3f} ± {se:.3f} (t={t:+.2f})   "
              f"any-full {a_hi}/{len(summary[hi]['f'])} vs {a_lo}/{len(summary[lo]['f'])} "
              f"Fisher p={p:.3f}")
        print(f"      fold-floor failures: order {hi} {summary[hi]['dead']}/{summary[hi]['n']}, "
              f"order {lo} {summary[lo]['dead']}/{summary[lo]['n']}")
    top = max(orders, key=lambda o: fef(summary[o]["best"]))
    print(f"\n  best sequence overall is order {top}: {summary[top]['best']['seq_id']} "
          f"fef {fef(summary[top]['best']):.3f} "
          f"(gate {summary[top]['best']['n_gated']}/{summary[top]['best']['n_draws']}, "
          f"median_sc {summary[top]['best']['median_sc']})")
    print("\n  ⇒ READ IT THIS WAY, then choose -- this script does NOT choose:")
    print("    • higher order wins on fef AND fold-failures do not climb  -> stacking HELPS;")
    print("      cycle 3 (orders 4/5/6, 38 seqs) is the right next move.")
    print("    • fef flat or falling with order, or fold-failures climbing -> stacking does NOT help;")
    print("      go back to SINGLES from the best double instead.")
    print("    ⛔ A difference with |t| < 2 or Fisher p > 0.05 is NOT a difference. Cycle 1 already")
    print("       produced three 'leaders' that were statistically identical; do not repeat that")
    print("       mistake one order up.")


if __name__ == "__main__":
    main()
