"""BACKBONE-LEVEL ranking for the egfrBR3 BRIDGE arm -- BOTH sites. 2026-10-04.

Same method as the site-A round (`backbone_passrate.py` + `backbone_separation.py`), widened so
each backbone is scored on site A AND site B from the same draws. The question it answers:
WHICH BACKBONES DESERVE A HIGHER-DRAW ROUND.

⛔ METHOD CARRIED OVER UNCHANGED, because these choices already cost this project a reversed winner:
  1. THE INDEPENDENT UNIT IS THE SEQUENCE, not the draw. 3 draws of one sequence are correlated,
     so every CI here is on n = SEQUENCES (10 per backbone), never on 30 draws.
     ⇒ breadth = "how many of this backbone's sequences land >=1 draw at the threshold".
  2. ⛔⛔ NEVER RANK ON PARTIAL DATA -- but MEASURED, not assumed.
     4 of 12 shards FAILED (00/01/04/05: FATAL, or no tarball), so 786 of 1180 sequences landed.
     ⭐ The loss was MEASURED to be clean, and that is what makes this arm rankable at all:
        - EVERY landed sequence has EXACTLY 3 of 3 draws (118/118 backbones) -- so no backbone is
          scored on a truncated draw set, which is the failure mode that reads TOO HIGH.
        - the loss removed WHOLE SEQUENCES, uniformly: every backbone kept 6 or 8 of 10
          (79 backbones at 6, 39 at 8). No backbone lost its draws, none lost all its sequences.
     ⇒ the DRAW-completeness gate is asserted and refuses; sequence coverage of 6/10 or 8/10 is
     recorded per backbone and the ranking is on RATES with Wilson CIs on the OBSERVED n, never on
     raw counts -- a 6-sequence and an 8-sequence backbone cannot be compared by count.
     ⛔ n = 6-8 sequences is SMALL. The CIs below are wide on purpose; read them, not the order.
  3. CONTACT IS THE GATE; `pose_rmsd` is never a cut.
  4. Thresholds are reported at several levels, because the threshold is a choice.

⛔⛔ SITE B IS SCORED ON FP7' (M30 DROPPED). M30 is contacted by 0 of 336 site-B designs and
   r(distance from M30, contact rate) = +0.967, so it was dropped. FP7 as-selected is reported
   alongside for continuity, but a backbone must not be
   penalised for a residue nothing reaches. Neither is hardcoded as "the" gate.

usage: backbone_passrate_bridge.py <contacts_br3_all.csv> <foldset.tsv> [out.csv]
"""
import sys, csv, math, collections, statistics as st
import numpy as np

A_FULL, B_FULL = 11, 7


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * (c - h), 100 * (c + h))


def ztest(k1, n1, k2, n2):
    p1, p2 = k1 / n1, k2 / n2
    se = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    if se == 0:
        return 0.0, 1.0
    z = (p1 - p2) / se
    return z, math.erfc(abs(z) / math.sqrt(2))


def bb_of(design):
    return design.rsplit("_id", 1)[0] if "_id" in design else design


def main(inp, tsv, outp=None):
    # ---- EXPECTED sequences per backbone, from the input the fold actually consumed
    exp = collections.defaultdict(set)
    for line in open(tsv):
        name = line.split("\t")[0].strip()
        if not name:
            continue
        name = name[3:] if name.startswith("BR_") else name
        base, _, sid = name.rpartition("_id")
        exp[base].add(sid)
    print(f"foldset {tsv.split('/')[-1]}: {len(exp)} backbones, "
          f"{sum(len(v) for v in exp.values())} sequences expected\n")

    rows = list(csv.DictReader(open(inp)))
    for r in rows:
        for k in ("n_A", "n_B7", "n_B7p", "n_B8", "n_target_contacts"):
            r[k] = int(r[k])

    # backbone -> seq -> list of draws
    bb = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        bb[bb_of(r["design"])][r["seq"]].append(r)

    recs, partial = [], []
    for k, seqs in bb.items():
        nexp = len(exp.get(k, ()))
        nobs = len(seqs)
        ndraw = sum(len(v) for v in seqs.values())
        # ⛔ THE GATE IS DRAW-COMPLETENESS, measured per sequence. A sequence missing a draw is
        #    the condition that reads too high. Sequence coverage (6/10, 8/10) is NOT a refusal
        #    -- it is recorded, and the ranking uses rates on the observed n instead.
        short_draw = any(len(v) != 3 for v in seqs.values())
        # breadth: fraction of SEQUENCES with at least one draw at the threshold
        fA = sum(1 for v in seqs.values() if max(x["n_A"] for x in v) >= A_FULL)
        f9 = sum(1 for v in seqs.values() if max(x["n_A"] for x in v) >= 9)
        fB = sum(1 for v in seqs.values() if max(x["n_B7p"] for x in v) >= B_FULL)
        fB5 = sum(1 for v in seqs.values() if max(x["n_B7p"] for x in v) >= 5)
        fBoth = sum(1 for v in seqs.values()
                    if any(x["n_A"] >= A_FULL and x["n_B7p"] >= B_FULL for x in v))
        fBothFP7 = sum(1 for v in seqs.values()
                       if any(x["n_A"] >= A_FULL and x["n_B7"] >= B_FULL for x in v))
        allA = [x["n_A"] for v in seqs.values() for x in v]
        allB = [x["n_B7p"] for v in seqs.values() for x in v]
        allT = [x["n_target_contacts"] for v in seqs.values() for x in v]
        rec = dict(backbone=k, n_seqs=nobs, n_seqs_expected=nexp, n_draws=ndraw,
                   complete="YES" if not short_draw else "NO(draws)",
                   cover=f"{nobs}/{nexp}",
                   pct_both=100.0 * fBoth / nobs, pct_full_B=100.0 * fB / nobs,
                   pct_full_A=100.0 * fA / nobs, pct_A_ge9=100.0 * f9 / nobs,
                   seq_full_A=fA, seq_A_ge9=f9, seq_full_B=fB, seq_B_ge5=fB5,
                   seq_both=fBoth, seq_both_fp7=fBothFP7,
                   med_A=float(np.median(allA)), med_B=float(np.median(allB)),
                   med_tgt=float(np.median(allT)),
                   max_A=max(allA), max_B=max(allB),
                   draws_full_A=sum(1 for x in allA if x >= A_FULL),
                   draws_full_B=sum(1 for x in allB if x >= B_FULL))
        (recs if rec["complete"] == "YES" else partial).append(rec)

    print(f"COMPLETE backbones: {len(recs)}   PARTIAL (not ranked): {len(partial)}   "
          f"total seen: {len(bb)}/{len(exp)}\n")

    # ---------------- ranking, complete backbones only
    recs.sort(key=lambda r: (-r["pct_both"], -r["pct_full_B"], -r["pct_full_A"], -r["med_A"]))
    hdr = (f"{'backbone':42s} {'cover':>6} {'BOTH':>9} {'fullB':>9} {'fullA':>9} "
           f"{'A>=9':>9} {'medA':>5} {'medB':>5} {'medT':>5}")
    print("RANKED on RATES (counts are over the OBSERVED sequences, shown as k/n -- a 6-seq and")
    print("an 8-seq backbone are NOT comparable by count). Primary key: spans BOTH sites.")
    print(hdr)
    print("-" * len(hdr))
    for r in recs[:25]:
        print(f"{r['backbone']:42s} {r['cover']:>6} "
              f"{r['seq_both']}/{r['n_seqs']} {r['pct_both']:5.0f}% "
              f"{r['seq_full_B']}/{r['n_seqs']} {r['pct_full_B']:5.0f}% "
              f"{r['seq_full_A']}/{r['n_seqs']} {r['pct_full_A']:5.0f}% "
              f"{r['seq_A_ge9']}/{r['n_seqs']} {r['pct_A_ge9']:5.0f}% "
              f"{r['med_A']:5.1f} {r['med_B']:5.1f} {r['med_tgt']:5.1f}")

    # ---------------- what the ranking can actually support
    ns = sorted({r["n_seqs"] for r in recs})
    nlab = f"{ns[0]}" if len(ns) == 1 else f"{ns[0]}-{ns[-1]}"
    print(f"\n=== CAN THIS RANKING SEPARATE ANYTHING? (n = {nlab} sequences per backbone) ===")
    if not recs:
        print("  ⛔ no draw-complete backbone to test")
    for label, key in (("BOTH sites", "seq_both"),
                       ("full site B (FP7')", "seq_full_B"),
                       ("full site A (11/11)", "seq_full_A"),
                       ("site A >=9/11", "seq_A_ge9")):
        if not recs:
            break
        srt = sorted(recs, key=lambda r: -(r[key] / r["n_seqs"]))
        top = srt[0]
        lo, hi = wilson(top[key], top["n_seqs"])
        print(f"\n  {label}: leader {top['backbone']} {top[key]}/{top['n_seqs']} "
              f"95% CI [{lo:.1f}, {hi:.1f}]%")
        nsep = 0
        for r in srt[1:]:
            z, p = ztest(top[key], top["n_seqs"], r[key], r["n_seqs"])
            if p < 0.05:
                nsep += 1
        print(f"    separable from {nsep} of {len(srt)-1} other backbones at p<0.05")
        if nsep == 0:
            print("    ⛔ NOT SEPARABLE FROM ANY OF THEM -- this ranking orders, it does not prove.")

    # ---------------- variance decomposition: backbone vs sequence
    print("\n=== BACKBONE CHOICE vs SEQUENCE CHOICE (site A, >=9/11 per-draw fraction) ===")
    means, within = [], []
    for k, seqs in bb.items():
        if len(seqs) < 2:
            continue
        fr = [sum(1 for x in v if x["n_A"] >= 9) / len(v) for v in seqs.values()]
        means.append(st.mean(fr))
        if len(fr) > 1:
            within.append(st.variance(fr))
    if len(means) > 1 and within:
        between = st.variance(means)
        wmean = st.mean(within)
        print(f"  BETWEEN-backbone variance : {between:.5f} (sd {math.sqrt(between):.3f})")
        print(f"  WITHIN-backbone variance  : {wmean:.5f} (sd {math.sqrt(wmean):.3f})")
        print(f"  ⇒ within/between = {wmean/between:.1f}x  "
              f"({'SEQUENCE' if wmean > between else 'BACKBONE'} choice dominates)")

    # ---------------- the both-site hits, named
    print("\n=== EVERY SEQUENCE THAT SPANNED BOTH SITES (the whole arm) ===")
    hits = []
    for k, seqs in bb.items():
        for sid, v in seqs.items():
            for x in v:
                if x["n_A"] >= A_FULL and x["n_B7p"] >= B_FULL:
                    hits.append((k, sid, x))
    if not hits:
        print("  none")
    for k, sid, x in sorted(hits):
        tag = "" if bb_of(k) in {r["backbone"] for r in recs} else "  ⛔ backbone PARTIAL"
        print(f"  {k}_id{sid:3s} model {x['model']}  A {x['n_A']}/11  "
              f"B(FP7') {x['n_B7p']}/7  B(FP7) {x['n_B7']}/7  tgt {x['n_target_contacts']}{tag}")

    if partial:
        print(f"\n⛔ PARTIAL BACKBONES -- NOT RANKED ({len(partial)}; the 4 dead shards). "
              f"A short backbone reads TOO HIGH.")
        for r in sorted(partial, key=lambda r: -r["seq_both"])[:12]:
            print(f"   {r['backbone']:42s} {r['n_seqs']}/{r['n_seqs_expected']} seqs  "
                  f"BOTH {r['seq_both']}  fullB {r['seq_full_B']}  fullA {r['seq_full_A']}")
        miss = sorted(set(exp) - set(bb))
        if miss:
            print(f"   plus {len(miss)} backbones with ZERO folds landed")

    if outp:
        cols = list(recs[0].keys())
        with open(outp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, restval="")
            w.writeheader()
            w.writerows(recs + partial)
        print(f"\nwrote {outp} ({len(recs)} complete + {len(partial)} partial)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
