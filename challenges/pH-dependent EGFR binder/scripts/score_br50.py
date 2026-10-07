"""egfrBR50 label per `decisions/0046` — the per-draw FLOOR, averaged over EVERY draw.

    per-draw score = min(n_A/n_A_total, n_B7p/n_B7p_total)   if the draw folds (sc_rmsd < SC_GATE)
                   = 0                                        otherwise
    label_BR       = mean of that over EVERY draw Boltz produced for the sequence

⛔⛔ A MISFOLDED OR MISPOSITIONED DRAW SCORES 0 AND IS NEVER DROPPED. That is the whole decision
   (scoring 0 is different from not reporting a score). The ONLY admissible missing value is
   a draw where Boltz itself failed and wrote no structure -- which in BR50 never happened
   (1200/1200 structures, 10/10 manifests rc=0). ⇒ this script REFUSES on a short join rather than
   averaging over a sample chosen by the thing being averaged, which is the measured bias
   (`corr(n_A, n_B7p) = +0.590`) that `0046` exists to remove.

⛔ FLOOR, NEVER A SUM OR A MEAN OF THE TWO SITES: a min is low if EITHER site is weak, so it
   enforces both sites in the SAME draw by construction. ⛔ Do not "improve" it into a weighted sum.

⭐ `full_both` (n_A == 11 AND n_B7p == 7) is WATCHED and reported, NEVER regressed on (`0046`).

⚠️ The totals are read from the data (`n_A_total`, `n_B7p_total`) and ASSERTED to be the 11 and 7
   the ADR names -- a silent change in site definition would otherwise rescale every label.

⭐ MULTI-SOURCE. `--contacts`/`--sc-dir` may be given more than once, paired in order, to pool
   campaigns that measured the SAME backbone (BR3's ids 1-10 and BR50's ids 11-100 are disjoint
   sequence sets on the same 8 backbones at the same 3 draws). Each source gets a `campaign` tag so
   a batch effect stays visible in the output instead of being averaged away.
   ⛔ POOLING IS NOT FREE: measured 2026-10-05, the 8 backbone means rank consistently between the
   two campaigns (Spearman +0.881) but their absolute levels do not track as well (Pearson +0.487).
   ⇒ Pool for n, and keep `campaign` in any model that uses the pooled table.
   ⛔ A sequence appearing in two sources is REFUSED, not silently deduplicated -- it would be a
   replicate, and 0045 says a replicate must never enter the table as a new design.

usage:
  score_br50.py --contacts <csv> --sc-dir <dir> [--contacts <csv2> --sc-dir <dir2> ...]
                --backbone <design name> [--seqs <tsv> [--seqs <tsv2>]]
                [--out <per_sequence.csv>] [--all-backbones]
"""
import sys, os, csv, glob, argparse, statistics, collections

SC_GATE   = 1.5     # decisions/0046: a draw "folds" iff sc_rmsd < 1.5
EXP_A     = 11      # n_A_total the ADR's formula names
EXP_B7P   = 7       # n_B7p_total the ADR's formula names


def load_sc(sc_dir):
    """-> {(design, seq, model): sc_rmsd}. Refuses on a duplicate key."""
    files = sorted(glob.glob(os.path.join(sc_dir, "sc_shard*_BR.csv")))
    if not files:
        sys.exit(f"⛔ no sc_shard*_BR.csv under {sc_dir}")
    sc = {}
    for p in files:
        for r in csv.DictReader(open(p)):
            k = (r["design"], r["seq"], r["model"])
            if k in sc:
                sys.exit(f"⛔ duplicate sc row {k} in {os.path.basename(p)} -- two shards hold the "
                         f"same draw; refusing rather than picking one")
            sc[k] = float(r["sc_rmsd"])
    return sc, files


def score_backbone(rows, sc, backbone):
    """-> (per_seq list, diagnostics dict). Refuses on any draw with no sc row.

    `rows` may come from several campaigns; each row carries `_campaign`, and a sequence id seen
    under two campaigns is refused rather than merged (0045: a replicate is not a new design).
    """
    per = collections.defaultdict(list)
    seen_campaign = {}
    for r in rows:
        c = r.get("_campaign", "")
        if r["seq"] in seen_campaign and seen_campaign[r["seq"]] != c:
            sys.exit(f"⛔ {backbone} seq {r['seq']} appears in BOTH '{seen_campaign[r['seq']]}' and "
                     f"'{c}'. That is a REPLICATE, not a new design -- 0045 forbids merging it into "
                     f"the table. Score the campaigns separately and compare. Refusing.")
        seen_campaign[r["seq"]] = c
        k = (r["design"], r["seq"], r["model"])
        if k not in sc:
            sys.exit(f"⛔ {backbone} seq {r['seq']} draw {r['model']}: NO self-consistency row. "
                     f"0046 admits a missing value ONLY when Boltz wrote no structure; a contacts "
                     f"row exists, so a structure DOES exist. Refusing -- fix the sc join.")
        if r["n_A_total"] != str(EXP_A) or r["n_B7p_total"] != str(EXP_B7P):
            sys.exit(f"⛔ {backbone} seq {r['seq']}: totals are "
                     f"n_A_total={r['n_A_total']} n_B7p_total={r['n_B7p_total']}, "
                     f"but 0046's formula is min(n_A/{EXP_A}, n_B7p/{EXP_B7P}). The site "
                     f"definition changed; every label would rescale silently. Refusing.")
        folds = sc[k] < SC_GATE
        s = min(int(r["n_A"]) / EXP_A, int(r["n_B7p"]) / EXP_B7P) if folds else 0.0
        per[r["seq"]].append({
            "score": s, "folds": folds, "sc_rmsd": sc[k],
            "n_A": int(r["n_A"]), "n_B7p": int(r["n_B7p"]),
            "full_both": r["full_both"] == "1", "campaign": r.get("_campaign", ""),
        })

    ndraws = {len(v) for v in per.values()}
    if len(ndraws) != 1:
        sys.exit(f"⛔ {backbone}: draws per sequence differ {sorted(ndraws)} -- an uneven draw "
                 f"count makes the mean incomparable across sequences. Refusing.")

    out = []
    for seqid, draws in sorted(per.items(), key=lambda kv: int(kv[0])):
        out.append({
            "design": backbone,
            "campaign": draws[0]["campaign"],
            "seq": seqid,
            "label_BR": sum(d["score"] for d in draws) / len(draws),
            "n_draws": len(draws),
            "n_folded": sum(1 for d in draws if d["folds"]),
            "n_full_both": sum(1 for d in draws if d["full_both"]),
            "best_draw": max(d["score"] for d in draws),
            "max_n_A": max(d["n_A"] for d in draws),
            "max_n_B7p": max(d["n_B7p"] for d in draws),
            "med_sc_rmsd": statistics.median(d["sc_rmsd"] for d in draws),
        })
    labels = [r["label_BR"] for r in out]
    diag = {
        "n_seqs": len(out),
        "n_draws_each": ndraws.pop(),
        "n_zero": sum(1 for x in labels if x == 0.0),
        "n_distinct": len({round(x, 12) for x in labels}),
        "max": max(labels), "mean": statistics.mean(labels),
        "sd": statistics.pstdev(labels) if len(labels) > 1 else 0.0,
        "full_both_draws": sum(r["n_full_both"] for r in out),
        "full_both_seqs": sum(1 for r in out if r["n_full_both"] > 0),
    }
    return out, diag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contacts", required=True, action="append")
    ap.add_argument("--sc-dir", required=True, action="append")
    ap.add_argument("--backbone")
    ap.add_argument("--all-backbones", action="store_true")
    ap.add_argument("--seqs", action="append",
                    help="tsv (BR_<backbone>_id<N>\\t<seq>) to attach sequences; repeatable")
    ap.add_argument("--out")
    a = ap.parse_args()
    if not a.backbone and not a.all_backbones:
        sys.exit("⛔ give --backbone or --all-backbones")

    if len(a.contacts) != len(a.sc_dir):
        sys.exit(f"⛔ {len(a.contacts)} --contacts but {len(a.sc_dir)} --sc-dir; they pair in order")
    rows, sc = [], {}
    for cpath, scdir in zip(a.contacts, a.sc_dir):
        tag = os.path.basename(os.path.dirname(os.path.abspath(cpath))) or os.path.basename(cpath)
        sub = list(csv.DictReader(open(cpath)))
        for r in sub:
            r["_campaign"] = tag
        s_i, scfiles = load_sc(scdir)
        dup = set(s_i) & set(sc)
        if dup:
            sys.exit(f"⛔ {len(dup)} (design,seq,draw) keys appear in two sc sources "
                     f"(e.g. {sorted(dup)[:2]}); sources must be disjoint. Refusing.")
        sc.update(s_i)
        rows += sub
        print(f"source '{tag}': contacts {len(sub)} rows · sc {len(s_i)} rows from "
              f"{len(scfiles)} shard files")
    print(f"pooled: contacts {len(rows)} rows · sc {len(sc)} rows "
          f"· fold gate sc_rmsd < {SC_GATE}")

    seqs = {}
    for sp in (a.seqs or []):
        for line in open(sp):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                seqs.setdefault(p[0], p[1])

    designs = sorted({r["design"] for r in rows}) if a.all_backbones else [a.backbone]
    allout = []
    print(f"\n{'backbone':40s} {'seqs':>5} {'zero':>5} {'distinct':>9} {'max':>7} "
          f"{'mean':>7} {'fb_draws':>9} {'fb_seqs':>8}")
    for d in designs:
        sub = [r for r in rows if r["design"] == d]
        if not sub:
            sys.exit(f"⛔ backbone {d} not in {a.contacts}")
        out, g = score_backbone(sub, sc, d)
        for r in out:
            key = f"BR_{d}_id{r['seq']}"
            r["name"] = key
            r["sequence"] = seqs.get(key, "")
        allout += out
        print(f"{d:40s} {g['n_seqs']:>5} {g['n_zero']:>5} {g['n_distinct']:>9} "
              f"{g['max']:>7.3f} {g['mean']:>7.3f} {g['full_both_draws']:>9} {g['full_both_seqs']:>8}")

    if a.seqs:
        miss = [r["name"] for r in allout if not r["sequence"]]
        if miss:
            print(f"⚠️  {len(miss)} of {len(allout)} rows had no sequence in --seqs "
                  f"(e.g. {miss[:2]})")
    if a.out:
        cols = ["name", "design", "campaign", "seq", "label_BR", "n_draws", "n_folded",
                "n_full_both", "best_draw", "max_n_A", "max_n_B7p", "med_sc_rmsd", "sequence"]
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in allout:
                w.writerow({c: r.get(c, "") for c in cols})
        print(f"\nwrote {a.out}: {len(allout)} rows")


if __name__ == "__main__":
    main()
