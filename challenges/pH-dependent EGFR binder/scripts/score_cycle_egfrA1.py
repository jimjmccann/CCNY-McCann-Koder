#!/usr/bin/env python3
"""Score a landed egfrA1 cycle, answer the UCB-vs-random question, and write the next table.

⭐ WHY THIS EXISTS. Cycle 0's rows were built by an inline heredoc, and
`merge_cycle_into_training.py` only POOLS sequences already present -- it SKIPS new ones
("is new but has NO SEQUENCE in the base table"), which is every sequence in a mutant cycle.
So the first real cycle had no working merge path. This is that path, as a file, because
cycles 2-10 each need it.

⛔ LABEL DEFINITION IS COPIED EXACTLY from the cycle-0 inline script and the burn-in builder.
   Any drift here silently changes what the surrogate is trained on:
     label             = mean n_epitope over GATE-PASSING draws (sc_rmsd < 1.5)
     gate_frac         = n_gated / n_draws ;  passes_floor = gate_frac >= 0.5
     mean_nep_all      = mean n_epitope over ALL draws
     full_epitope_frac = fraction of ALL draws with n_epitope == 11
     median_sc         = median sc_rmsd over ALL draws
     mean_tgt          = mean n_target_contacts over ALL draws
⛔ The floor is a CANDIDATE FILTER, never a term summed into the label (the dftfit lesson).

⭐⭐ THE REPLICATE HOLD-OUT (`decisions/0045`, built 2026-10-05). A cycle may carry RE-FOLDS of
   earlier designs riding in its own slots (`propose_cycle_egfrA1.py` reserves them). Those folds
   are real measurements of sequences ALREADY IN THE TABLE, so:
     ⛔⛔ THEY ARE NEVER MERGED INTO `training_cN.csv`. They are duplicate sequences; merging them
        would put one sequence in the table twice, corrupt the hub detection (which counts
        Hamming-1 neighbours), double-weight one design in anything fitted on the table, and let
        a re-fold become the `parent` for the next cycle.
     ⭐ They are scored like everything else, compared against their PRIOR values, and written to
        `<cycle>_replicates_SCORED.csv` with the test-retest read.
   HOW THEY ARE IDENTIFIED: from `<cycle>_replicates.csv` ONLY -- the file whose name is derived
   from the cycle map path handed in. ⛔ NEVER from the id. "id9" as a prefix already matches ten
   REAL cycle-0 designs (id9, id90, id91..id99, id901, id902); prefix matching would hold out
   genuine data. If the file is absent, this script behaves EXACTLY as it did before the
   ride-along existed, which is what makes cycles 0-3 reproduce byte-for-byte.

usage: score_cycle_egfrA1.py <base_training.csv> <contacts.csv> <sc.csv> <cycle.tsv>
                             <cycle_map.csv> <out_training.csv> <out_parent.txt>
       (the replicate map, if any, is found beside <cycle_map.csv> as <cycle>_replicates.csv)
"""
import csv, os, sys, statistics, collections
from math import comb

GATE, FLOOR, FULL = 1.5, 0.5, 11


def fisher(a, b, c, d):
    n = a + b + c + d; r1 = a + b; c1 = a + c
    def p(x): return comb(r1, x) * comb(n - r1, c1 - x) / comb(n, c1)
    o = p(a); lo = max(0, c1 - (n - r1)); hi = min(r1, c1)
    return sum(p(x) for x in range(lo, hi + 1) if p(x) <= o + 1e-12)


def welch(xs, ys):
    """mean difference, SE, and t -- no scipy in this env."""
    mx, my = statistics.mean(xs), statistics.mean(ys)
    vx = statistics.variance(xs) / len(xs) if len(xs) > 1 else 0.0
    vy = statistics.variance(ys) / len(ys) if len(ys) > 1 else 0.0
    se = (vx + vy) ** 0.5
    return mx - my, se, ((mx - my) / se if se else float("nan"))


def pearson(xs, ys):
    n = len(xs)
    if n < 3: return float("nan")
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sx = sum((x - mx) ** 2 for x in xs) ** 0.5
    sy = sum((y - my) ** 2 for y in ys) ** 0.5
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy) if sx and sy else float("nan")


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v); i = 0
        while i < len(order):                      # average ties
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]: j += 1
            for k in range(i, j + 1): r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))


def load_replicates(mapf):
    """{replicate_id: prior-row} from the sibling replicates file, or {} if there is none.

    ⛔ DERIVED BY CONVENTION from the cycle map path, matching what `propose_cycle_egfrA1.py`
       writes, so the driver needs no new argument and old cycles keep working untouched.
    """
    if not mapf.endswith("_map.csv"):
        return {}, None
    q = mapf[:-len("_map.csv")] + "_replicates.csv"
    if not os.path.exists(q):
        return {}, None
    rep = {r["replicate_id"]: r for r in csv.DictReader(open(q))}
    print(f"replicate map {q}: {len(rep)} re-folds declared "
          f"({', '.join(sorted(rep))})")
    return rep, q


def lsq_slope(xs, ys):
    """Least-squares slope of ys on xs. The regression-to-mean number."""
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = statistics.mean(xs), statistics.mean(ys)
    den = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den if den else float("nan")


def report_replicates(rep, stat, outq):
    """Fresh vs prior for every re-fold: the only noise estimate this campaign has."""
    got = [(rid, m, stat[rid[2:]]) for rid, m in sorted(rep.items()) if rid[2:] in stat]
    missing = sorted(rid for rid in rep if rid[2:] not in stat)
    if missing:
        print(f"  ⚠️ {len(missing)} declared replicates produced NO structures: "
              f"{missing} -- they are simply absent from the read, not zeros")
    if not got:
        print("=== REPLICATES: none landed -- no test-retest this cycle ===")
        return
    print(f"\n=== REPLICATES: fresh vs prior (decisions/0045) ===")
    print(f"  {'replicate':10s} {'source':8s} {'prior_fef':>9s} {'fresh_fef':>9s} "
          f"{'delta':>7s}  {'prior_lab':>9s} {'fresh_lab':>9s}  mutations")
    xs, ys = [], []
    rows_out = []
    for rid, m, st in got:
        try:
            pf = float(m["prior_fef"])
        except (TypeError, ValueError):
            pf = float("nan")
        ff = float(st["full_epitope_frac"])
        if pf == pf:
            xs.append(pf); ys.append(ff)
        print(f"  {rid:10s} {m['source_id']:8s} {pf:9.3f} {ff:9.3f} {ff-pf:+7.3f}  "
              f"{str(m['prior_label']):>9s} {str(st['label']):>9s}  {m['mutations']}")
        rows_out.append({"replicate_id": rid, "source_id": m["source_id"],
                         "mutations": m["mutations"],
                         "prior_fef": m["prior_fef"], "fresh_fef": st["full_epitope_frac"],
                         "delta_fef": round(ff - pf, 3) if pf == pf else "",
                         "prior_label": m["prior_label"], "fresh_label": st["label"],
                         "prior_median_sc": m.get("prior_median_sc", ""),
                         "fresh_median_sc": st["median_sc"],
                         "prior_n_draws": m.get("prior_n_draws", ""),
                         "fresh_n_draws": st["n_draws"],
                         "fresh_gate": f"{st['n_gated']}/{st['n_draws']}"})
    if len(xs) >= 2:
        d = [y - x for x, y in zip(xs, ys)]
        sl = lsq_slope(xs, ys)
        print(f"\n  n={len(xs)}  mean delta {statistics.mean(d):+.3f}  "
              f"mean |delta| {statistics.mean(abs(v) for v in d):.3f}  "
              f"max |delta| {max(abs(v) for v in d):.3f}")
        print(f"  Pearson(prior, fresh) {pearson(xs, ys):+.3f}   slope {sl:+.3f}")
        print("  ⇒ READ IT THIS WAY: slope well BELOW 1.0 is regression to the mean -- the cycle's")
        print("     apparent gain is partly SELECTION, and the champions' priors were optimistic.")
        print("     mean |delta| is the per-design noise floor; a difference smaller than it is")
        print("     ⛔ NOT a difference. A negative mean delta on champions is EXPECTED, not a")
        print("     regression in the designs.")
    else:
        print("  ⚠️ fewer than 2 comparable replicates -- no slope, no noise estimate")
    with open(outq, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0]))
        w.writeheader(); w.writerows(rows_out)
    print(f"  wrote {outq}")


def main(base, contacts, scf, tsv, mapf, outp, parentp):
    rows = list(csv.DictReader(open(base)))
    flds = list(rows[0].keys())
    have = {r["seq_id"] for r in rows}
    seqs = dict(l.rstrip("\n").split("\t") for l in open(tsv) if "\t" in l)
    # fold_id in the tsv is  A_<design>_id<N>;  contacts/sc key on the bare <N>
    by_sid = {k.rsplit("_id", 1)[1]: v for k, v in seqs.items() if "_id" in k}
    meta = {r["fold_id"]: r for r in csv.DictReader(open(mapf))}
    rep, repq = load_replicates(mapf)
    # ⛔ the base table's SEQUENCES, for the corruption guard below. seq_id collisions were already
    #    caught; a SEQUENCE collision under a new id was not, and that is what a stray replicate is.
    base_seq = {}
    for r in rows:
        if r.get("sequence"):
            base_seq.setdefault(r["sequence"], r["seq_id"])

    con = {(r["seq"], r["model"]): r for r in csv.DictReader(open(contacts))}
    sc = {(r["seq"], r["model"]): r for r in csv.DictReader(open(scf))}
    assert set(con) == set(sc), "⛔ contacts/sc key mismatch"

    per = collections.defaultdict(list)
    for k, c in con.items():
        per[k[0]].append(dict(nep=int(c["n_epitope"]), tgt=int(c["n_target_contacts"]),
                              sc=float(sc[k]["sc_rmsd"]), gate=float(sc[k]["sc_rmsd"]) < GATE))

    stat, added, pooled, held = {}, 0, 0, 0
    for sid, v in sorted(per.items()):
        g = [x for x in v if x["gate"]]
        row = {
            "seq_id": f"id{sid}", "sequence": by_sid.get(sid, ""),
            "n_draws": len(v), "n_gated": len(g),
            "gate_frac": round(len(g) / len(v), 3),
            "label": round(statistics.mean(x["nep"] for x in g), 4) if g else "",
            "mean_nep_all": round(statistics.mean(x["nep"] for x in v), 4),
            "full_epitope_frac": round(sum(1 for x in v if x["nep"] == FULL) / len(v), 3),
            "median_sc": round(statistics.median(x["sc"] for x in v), 3),
            "mean_tgt": round(statistics.mean(x["tgt"] for x in v), 1),
            "passes_floor": int(len(g) / len(v) >= FLOOR)}
        stat[sid] = row
        if f"id{sid}" in have:
            print(f"  ⚠️ id{sid} ALREADY in base -- pooling not implemented here; "
                  f"use merge_cycle_into_training.py for repeats"); pooled += 1
            continue
        if not row["sequence"]:
            print(f"  ⛔ id{sid} has NO sequence in {tsv} -- SKIPPED"); continue
        # ⛔⛔ THE HOLD-OUT (decisions/0045 item 5). A declared replicate is scored into `stat`,
        #    where the test-retest reads it, and NEVER appended to the training table.
        if f"id{sid}" in rep:
            held += 1
            continue
        # ⛔⛔ THE CORRUPTION GUARD. A new id carrying a sequence the base table already holds is
        #    either an UNDECLARED replicate or a proposer bug; either way merging it puts one
        #    sequence in the table twice. The table had ZERO duplicate sequences when this was
        #    written (measured on training_c3.csv, 202 rows), so this fires only on real drift.
        if row["sequence"] in base_seq:
            sys.exit(f"⛔ id{sid} carries a sequence ALREADY in the base table as "
                     f"{base_seq[row['sequence']]}, but is NOT declared in the replicate map"
                     + (f" ({repq})" if repq else " (no replicate map found beside the cycle map)")
                     + ". Merging it would hold one sequence twice -- which corrupts the hub "
                       "detection and double-weights one design. Refusing. If this IS a re-fold, "
                       "declare it in the replicates file; if it is not, the proposer's "
                       "never-re-fold exclusion has failed.")
        rows.append(row); added += 1

    # ---- the cycle's own question: did UCB beat random? ----
    # ⛔ ONLY A PROPOSAL-ARM CYCLE HAS THIS QUESTION. Cycle 1's map carried `arm`/`pred_label`
    #    (25 UCB vs 25 random); a COMBINATION cycle's map (cycle 2 onward) carries
    #    fold_id/order/n_mut/mutations/positions and NO arm, because the surrogate is out of the
    #    loop by measurement (OOB RMSE/sd(y) 1.14 on label, 1.04 on full_epitope_frac). Reading
    #    m["arm"] unconditionally raises KeyError on such a map and kills the scoring run
    #    before it writes the training table.
    #    ⇒ Skip the head-to-head when the columns are absent; never invent a comparison.
    arms = collections.defaultdict(list)
    for sid, row in stat.items():
        m = meta.get(f"id{sid}")
        if m and row["label"] != "" and m.get("arm") and m.get("pred_label"):
            arms[m["arm"]].append((sid, float(row["label"]), float(m["pred_label"]), row))
    if not arms:
        print("\n=== no proposal arms in this cycle's map -- UCB-vs-random skipped ===")
        print("  (a combination cycle has no arm split and no surrogate prediction to score)")
    else:
        print("\n=== UCB vs RANDOM -- the head-to-head this cycle was built to answer ===")
    for a in ("ucb", "rand"):
        v = [x[1] for x in arms.get(a, [])]
        if v:
            print(f"  {a:5s} n={len(v):2d}  mean label {statistics.mean(v):6.3f}  "
                  f"median {statistics.median(v):6.3f}  max {max(v):6.3f}  "
                  f"SE {(statistics.variance(v)/len(v))**0.5 if len(v)>1 else 0:.3f}")
    if arms.get("ucb") and arms.get("rand"):
        d, se, t = welch([x[1] for x in arms["ucb"]], [x[1] for x in arms["rand"]])
        print(f"  ⇒ UCB − random = {d:+.3f} ± {se:.3f}  (t = {t:+.2f})")
        fu = sum(1 for x in arms["ucb"] if x[3]["full_epitope_frac"] > 0)
        fr = sum(1 for x in arms["rand"] if x[3]["full_epitope_frac"] > 0)
        print(f"  ⇒ seqs with ANY full-epitope draw: ucb {fu}/{len(arms['ucb'])} vs "
              f"rand {fr}/{len(arms['rand'])}  Fisher p={fisher(fu,len(arms['ucb'])-fu,fr,len(arms['rand'])-fr):.3f}")

    # ---- was the surrogate's prediction worth anything? ----
    allp = [(x[2], x[1]) for a in arms for x in arms[a]]
    if len(allp) > 3:
        px, py = [p for p, _ in allp], [y for _, y in allp]
        print(f"\n=== SURROGATE SKILL on this cycle (pred_label vs measured label, n={len(allp)}) ===")
        print(f"  Pearson {pearson(px,py):+.3f}   Spearman {spearman(px,py):+.3f}")
        print("  ⛔ A near-zero correlation means the ensemble is NOT yet predictive; the walk is")
        print("     then driven by measurement, not by the model. That is a finding, not a failure.")

    # ---- the walk: did anything beat the parent? ----
    # ⛔ THE FLOOR APPLIES TO THE BASELINE TOO. Without it the "previous best" came back as
    #    id93 (label 10.3333 on gate 3/10, median sc_rmsd 1.871) -- a sequence that FAILS the
    #    fold gate 7 draws in 10. Anchoring a cycle's progress on a thin-data outlier flatters
    #    the result; never rank on partial data. The real baseline is the floor-passing best.
    pr = [r for r in rows if r["seq_id"] in have and r["label"] != ""
          and r["passes_floor"] in (1, "1")]
    base_best = max(pr, key=lambda r: float(r["label"])) if pr else None
    rows.sort(key=lambda r: -(float(r["label"]) if r["label"] else -1))
    best = next(r for r in rows if r["label"] != "" and r["passes_floor"] in (1, "1"))
    print(f"\n=== THE WALK ===")
    if base_best:
        print(f"  previous best (in base table): {base_best['seq_id']} label {base_best['label']}")
    print(f"  new best overall:              {best['seq_id']} label {best['label']} "
          f"(gate {best['n_gated']}/{best['n_draws']}, med_sc {best['median_sc']})")
    nb = [r for r in rows if r["label"] != "" and base_best
          and float(r["label"]) > float(base_best["label"]) and r["seq_id"] not in have]
    print(f"  cycle sequences beating the previous best: {len(nb)}")
    for r in nb[:8]:
        m = meta.get(r["seq_id"], {})
        print(f"    {r['seq_id']:8s} {m.get('mutation','?'):8s} {m.get('arm','?'):5s} "
              f"label {r['label']}  full11 {r['full_epitope_frac']}")

    if rep:
        report_replicates(rep, stat, repq[:-len(".csv")] + "_SCORED.csv")

    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=flds); w.writeheader(); w.writerows(rows)
    open(parentp, "w").write(best["sequence"] + "\n")
    ok = [r for r in rows if r["passes_floor"] in (1, "1")]
    print(f"\nwrote {outp}: {len(rows)} sequences ({added} added, {pooled} skipped-as-repeat, "
          f"{held} held out as replicates), {len(ok)} pass the floor")
    print(f"wrote {parentp}: parent for the next cycle = {best['seq_id']}")


if __name__ == "__main__":
    main(*sys.argv[1:8])
