"""egfrBR50 DEPTH ROUND 2 (`d2`): does id47 or id76 beat id85? With id85 as the IN-BATCH ANCHOR.

NEW FILE 2026-10-06. ⛔ NOT an edit of `bridge_depth_analyze.py` — that file is the committed
provenance of a SCORED round, and its group set is hard-coded to d1's four sequences
(85 / 1104 / 1107 / 1109), so it refuses d2 by its own gate. A sibling, exactly as
`build_bridge_depth_round_d2.py` is a sibling of the d1 builder.

WHAT CHANGES FROM d1, AND WHY

  d1 asked "do the two label leaders hold 1.000, and do they separate?" — a two-leader question
  with two mid-table anchors. d2 asks a DIFFERENT question, so the statistics it needs differ:

    id47  12 replicates x 3 draws = 36 draws   prior 2/3 full_both, never deep-tested
    id76  12 replicates x 3 draws = 36 draws   prior 2/3 full_both, never deep-tested
    id85   8 replicates x 3 draws = 24 draws   ⭐ IN-BATCH ANCHOR, d1 deep 0.936, fb 71/90

  ⛔⛔ THE ANCHOR IS A VETO, NOT A CANDIDATE, and it is a SINGLE anchor rather than d1's two.
  d1 could tell a batch shift from regression to the mean because its two anchors could move in
  OPPOSITE directions. ⇒ d2 cannot do that. All it can do is ask whether id85 reproduces its OWN
  d1 number; if it does not, the round's pre-stated stop condition fires and the comparison is
  VOID. This script reports that first and refuses to rank if it fails.

  ⛔ GROUPS AND PRIORS ARE READ FROM THE MAP, never hard-coded, so the same script serves d3.
  The anchor's reference values are the only numbers typed in, and they are typed in because they
  come from the PREVIOUS round's committed summary rather than from this round's map
  (the map's `prior_label_BR` for id85 is its 3-draw c1 label 1.000, NOT its d1 deep 0.936).

⛔ SEPARATION IS TESTED ON THE REPLICATE SPREAD, never on within-batch draw s.e.: the within-batch
  form UNDERSTATES re-fold noise by 1.56x, measured (within-batch-se-understates-refold-noise).
⛔ Scored SEPARATELY from the c1 table (`0045` item 5) — these are replicates, not new designs.

usage:
  bridge_depth2_analyze.py --labels <label_BR_0046_DEPTH_d2.csv> --map <bridge_depth2_map.csv>
      [--anchor-id 85] [--anchor-deep 0.9360750360750361] [--anchor-fb 71] [--anchor-fb-draws 90]
      [--out <deep_summary_d2.csv>] [--per-rep-out <per_replicate_d2.csv>]
"""
import sys, csv, argparse, statistics, collections


def prop_z(f1, n1, f2, n2):
    p1, p2 = f1 / n1, f2 / n2
    pool = (f1 + f2) / (n1 + n2)
    se = (pool * (1 - pool) * (1 / n1 + 1 / n2)) ** 0.5
    return p1, p2, ((p1 - p2) / se if se else float("nan")), se


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--anchor-id", type=int, default=85)
    ap.add_argument("--anchor-deep", type=float, default=0.9360750360750361,
                    help="the anchor's deep label from the PREVIOUS round's committed summary")
    ap.add_argument("--anchor-se", type=float, default=0.01827314869979729,
                    help="the anchor's REPLICATE s.e. in the previous round. ⛔ Required for a "
                         "correct test: comparing two ESTIMATES combines BOTH s.e.s. Using this "
                         "round's s.e. alone inflates a 1.80-s.e. difference to 3.55 s.e., "
                         "which would declare a valid round void.")
    ap.add_argument("--anchor-fb", type=int, default=71)
    ap.add_argument("--anchor-fb-draws", type=int, default=90)
    ap.add_argument("--anchor-tol-se", type=float, default=1.96,
                    help="how far the anchor may miss its own number before the round is VOID")
    ap.add_argument("--out")
    ap.add_argument("--per-rep-out")
    a = ap.parse_args()

    src, kind, prior = {}, {}, {}
    for r in csv.DictReader(open(a.map)):
        oid = int(r["depth_of_id"])
        src[r["id"]] = oid
        kind[oid] = r["kind"]
        prior[oid] = float(r["prior_label_BR"])

    rows = list(csv.DictReader(open(a.labels)))
    if not rows:
        sys.exit(f"⛔ empty label table {a.labels}")

    groups, unmapped = collections.defaultdict(list), []
    for r in rows:
        if r.get("seq") not in src:
            unmapped.append(r.get("seq")); continue
        groups[src[r["seq"]]].append(r)
    if unmapped:
        sys.exit(f"⛔ {len(unmapped)} label rows not in the map ({unmapped[:5]}). Refusing rather "
                 f"than summarising a subset chosen by whatever went missing.")
    missing = sorted(set(src.values()) - set(groups))
    if missing:
        sys.exit(f"⛔ no label rows for group(s) {missing} that the map declares. Refusing a "
                 f"partial comparison: never rank on partial data.")
    if a.anchor_id not in groups:
        sys.exit(f"⛔ the anchor id{a.anchor_id} is not in this round. Refusing: without the "
                 f"in-batch anchor there is nothing to detect a batch shift with.")

    summ = {}
    per_rep = []
    for oid in sorted(groups):
        g = sorted(groups[oid], key=lambda r: int(r["seq"]))
        lab = [float(r["label_BR"]) for r in g]
        n = len(lab)
        sd = statistics.stdev(lab) if n > 1 else float("nan")
        fb = sum(int(r["n_full_both"]) for r in g)
        draws = sum(int(r["n_draws"]) for r in g)
        folded = sum(int(r["n_folded"]) for r in g)
        summ[oid] = dict(orig_id=oid, kind=kind[oid], reps=n, draws=draws,
                         deep_mean=statistics.mean(lab), rep_sd=sd,
                         rep_se=(sd / n ** 0.5 if n > 1 else float("nan")),
                         rep_min=min(lab), rep_max=max(lab),
                         rep_median=statistics.median(lab),
                         n_at_1000=sum(1 for v in lab if v >= 0.9999),
                         n_zero=sum(1 for v in lab if v == 0.0),
                         fb_draws=fb, fb_rate=fb / draws if draws else float("nan"),
                         fold_rate=folded / draws if draws else float("nan"),
                         prior_3draw=prior[oid])
        for r in g:
            per_rep.append(dict(orig_id=oid, kind=kind[oid], rep_seq=r["seq"],
                                label_BR=r["label_BR"], n_draws=r["n_draws"],
                                n_folded=r["n_folded"], n_full_both=r["n_full_both"],
                                med_sc_rmsd=r["med_sc_rmsd"]))

    print(f"\n=== DEEP ESTIMATES — {a.labels.split('/')[-1]} ===")
    hdr = (f"{'orig':>6} {'kind':<8} {'reps':>5} {'draws':>6} {'deep':>7} {'rep_sd':>7} "
           f"{'s.e.':>6} {'min':>6} {'max':>6} {'@1.000':>7} {'0.000':>6} {'fb_rate':>8} "
           f"{'fold':>6} {'3-draw prior':>13}")
    print(hdr)
    for oid in sorted(summ):
        s = summ[oid]
        print(f"{oid:>6} {s['kind']:<8} {s['reps']:>5} {s['draws']:>6} {s['deep_mean']:>7.3f} "
              f"{s['rep_sd']:>7.3f} {s['rep_se']:>6.3f} {s['rep_min']:>6.3f} {s['rep_max']:>6.3f} "
              f"{s['n_at_1000']:>7} {s['n_zero']:>6} {s['fb_rate']:>8.3f} {s['fold_rate']:>6.3f} "
              f"{s['prior_3draw']:>13.3f}")

    # ---- Q1: THE ANCHOR VETO, FIRST, BEFORE ANY RANKING
    print(f"\n=== Q1. ⛔⛔ THE ANCHOR VETO — does id{a.anchor_id} reproduce its OWN d1 number? ===")
    an = summ[a.anchor_id]
    d_lab = an["deep_mean"] - a.anchor_deep
    se_comb = (an["rep_se"] ** 2 + a.anchor_se ** 2) ** 0.5
    print(f"  LABEL     prev deep {a.anchor_deep:.3f} (s.e. {a.anchor_se:.3f})  ->  this round "
          f"{an['deep_mean']:.3f} (s.e. {an['rep_se']:.3f}, {an['draws']} draws)   "
          f"delta {d_lab:+.3f}")
    print(f"            ⛔ the test COMBINES both s.e.s: sqrt({an['rep_se']:.3f}^2 + "
          f"{a.anchor_se:.3f}^2) = {se_comb:.4f}  ⇒  {abs(d_lab)/se_comb:.2f} s.e. away")
    print(f"            (using this round's s.e. alone would read {abs(d_lab)/an['rep_se']:.2f} "
          f"s.e. and is WRONG -- the previous estimate has uncertainty too)"
          if an["rep_se"] > 0 else "            spread is zero; the s.e. test is undefined")
    p2, p1, z_an, _ = prop_z(an["fb_draws"], an["draws"], a.anchor_fb, a.anchor_fb_draws)
    print(f"  full_both d1 {a.anchor_fb}/{a.anchor_fb_draws} = {p1:.3f}  ->  d2 "
          f"{an['fb_draws']}/{an['draws']} = {p2:.3f}   pooled-proportion z = {z_an:+.2f}")
    void = (se_comb > 0 and abs(d_lab) / se_comb > a.anchor_tol_se) or abs(z_an) > 1.96
    if void:
        print(f"  ⛔⛔ THE ANCHOR MISSED ITS OWN NUMBER beyond {a.anchor_tol_se} s.e. "
              f"⇒ THE BATCH SHIFTED AND THE CROSS-ROUND COMPARISON IS VOID.")
        print(f"     This condition was stated explicitly before the round. Within-round ranking of "
              f"the challengers against THIS batch's id85 is still readable; anything compared to "
              f"d1's 0.936 is not.")
    else:
        print(f"  ⭐ the anchor HOLDS within {a.anchor_tol_se} s.e. on the label and 1.96 s.e. on "
              f"full_both ⇒ no batch shift detected; this round is comparable to the previous one.")
        print(f"  ⚠️ 'holds' is NOT 'identical': the label sits {d_lab:+.3f} "
              f"({abs(d_lab)/se_comb:.2f} s.e.) from its previous value. Quote the direction when "
              f"the margin being judged is of that size. It is not, here -- see Q2.")
    print(f"  ⚠️ ONE anchor cannot distinguish a batch shift from a change in THIS sequence. d1 "
          f"used TWO so they could move in opposite directions. This is a weaker control by "
          f"construction, and it is the control the round was built with.")

    # ---- Q2: do the challengers beat the anchor?
    others = [o for o in sorted(summ) if o != a.anchor_id]
    print(f"\n=== Q2. ⭐⭐ DO THE CHALLENGERS BEAT id{a.anchor_id}? — on the LABEL ===")
    print("  (gap s.e. from the REPLICATE spread of both groups, never within-batch draw s.e.)")
    for oid in others:
        s = summ[oid]
        gap = s["deep_mean"] - an["deep_mean"]
        se = (s["rep_se"] ** 2 + an["rep_se"] ** 2) ** 0.5
        verdict = ("BEATS it" if gap > 1.96 * se else
                   "LOSES to it" if gap < -1.96 * se else "NOT SEPARABLE from it")
        print(f"  id{oid:<5} deep {s['deep_mean']:.3f} vs {an['deep_mean']:.3f}  "
              f"gap {gap:+.4f}  s.e. {se:.4f}  gap/s.e. {gap/se:+.2f}  "
              f"95% CI [{gap-1.96*se:+.4f}, {gap+1.96*se:+.4f}]  ⇒ {verdict}")
    if len(others) >= 2:
        o1, o2 = others[0], others[1]
        g = summ[o1]["deep_mean"] - summ[o2]["deep_mean"]
        se = (summ[o1]["rep_se"] ** 2 + summ[o2]["rep_se"] ** 2) ** 0.5
        print(f"  and against EACH OTHER: id{o1} - id{o2} = {g:+.4f}  s.e. {se:.4f}  "
              f"⇒ {'SEPARATE' if abs(g) > 1.96 * se else 'NOT separable'}")

    # ---- Q3: full_both, the watched criterion
    print(f"\n=== Q3. full_both RATE — the criterion 0046 WATCHES, and it is NOT ceiling-compressed ===")
    for oid in sorted(summ):
        s = summ[oid]
        print(f"    id{oid:<5} {s['kind']:<8} {s['fb_draws']:>3}/{s['draws']:<3} = {s['fb_rate']:.3f}")
    for oid in others:
        s = summ[oid]
        p1, p2, z, _ = prop_z(s["fb_draws"], s["draws"], an["fb_draws"], an["draws"])
        verdict = ("BEATS it" if z > 1.96 else "LOSES to it" if z < -1.96 else "NOT separable")
        print(f"  id{oid} {p1:.3f} vs id{a.anchor_id} {p2:.3f}   diff {p1-p2:+.3f}   "
              f"z = {z:+.2f}  ⇒ {verdict}")
    print("  ⛔ READ Q2 AND Q3 TOGETHER. The label can be ceiling-compressed near the min-floor's")
    print("     top; full_both is a per-draw rate and is not. If they disagree, say which is which.")

    # ---- Q4: the noise this round measures
    print("\n=== Q4. 3-DRAW RE-FOLD NOISE, measured directly in THIS batch ===")
    for oid in sorted(summ):
        s = summ[oid]
        print(f"    id{oid:<5} {s['kind']:<8} n={s['reps']:<3} noise sd = {s['rep_sd']:.3f}  "
              f"(d1 measured 0.069-0.100 for leaders, 0.250-0.271 mid-table)")
    print("  ⚠️ heteroscedastic, as d1 found: a sequence pinned near a ceiling or a floor has")
    print("     near-zero spread by construction, so a small sd is not by itself good news.")

    # ---- Q5: what the 3-draw prior was worth
    print("\n=== Q5. ⛔ WHAT THE 3-DRAW PRIOR WAS WORTH — the round's methodological point ===")
    print(f"  {'id':>6} {'3-draw prior':>13} {'deep':>7} {'delta':>8}")
    for oid in sorted(summ):
        s = summ[oid]
        print(f"  {oid:>6} {s['prior_3draw']:>13.3f} {s['deep_mean']:>7.3f} "
              f"{s['deep_mean']-s['prior_3draw']:>+8.3f}")
    pri = [(summ[o]["prior_3draw"], summ[o]["deep_mean"], o) for o in sorted(summ)]
    rank_prior = [o for _p, _d, o in sorted(pri, key=lambda t: -t[0])]
    rank_deep = [o for _p, _d, o in sorted(pri, key=lambda t: -t[1])]
    print(f"  rank by 3-draw prior: {rank_prior}")
    print(f"  rank at depth       : {rank_deep}")
    print(f"  ⇒ the 3-draw ordering {'HELD' if rank_prior == rank_deep else 'DID NOT HOLD'}.")

    if a.out:
        cols = list(next(iter(summ.values())).keys())
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
            for oid in sorted(summ):
                w.writerow(summ[oid])
        print(f"\nwrote {a.out}: {len(summ)} rows")
    if a.per_rep_out:
        with open(a.per_rep_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(per_rep[0].keys())); w.writeheader()
            w.writerows(per_rep)
        print(f"wrote {a.per_rep_out}: {len(per_rep)} rows")


if __name__ == "__main__":
    main()
