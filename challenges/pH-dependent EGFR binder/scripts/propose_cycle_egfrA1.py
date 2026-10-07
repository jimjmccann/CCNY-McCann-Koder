#!/usr/bin/env python3
"""Propose the next egfrA1 cycle: recompute the accepted set, enumerate the untried combinations.

⭐ WHY THIS EXISTS. It is the proposal step of the egfrA1 cycle loop: recompute the accepted set
and enumerate the untried combinations, so a cycle can be proposed from its training table alone.
Without it the loop has no step 3 and cannot run unattended.

⛔ THERE IS NO MODEL IN THIS LOOP, and that is MEASURED, not a preference. On cycle 1 the
   surrogate's OOB RMSE / sd(y) was 1.14 on `label` and 1.04 on `full_epitope_frac` -- at or worse
   than predicting the mean -- and UCB did not beat random (+0.394 +/- 0.634, t=+0.62, Fisher
   p=1.000 on any-full-epitope). Trained on the saturated label it predicted 11.696 and 12.101,
   above the physical ceiling of 11. Exhaustive coverage needs no ranking, so nothing is lost by
   dropping it. ⇒ Do not re-introduce UCB here without re-measuring OOB first.

⛔ THE LABEL SATURATES. The site-A epitope has 11 residues, so `label` (mean n_epitope over
   gate-passing draws) has a CEILING of 11 and cycle 1's top three sat at identical means
   (dmean +0.00, Fisher p=0.628-1.000). The accepted-set test therefore requires a sequence to
   beat the reference on **BOTH** `label` AND `full_epitope_frac`; `label` alone can no longer
   discriminate and `full_epitope_frac` alone is noisier at 10 draws.

⭐ HOW THE REFERENCE IS CHOSEN -- the one inference in this file, so it is made measurable rather
   than named. The accepted set is only meaningful relative to the sequence whose single mutants
   were actually MEASURED. That sequence is the HUB: the training row with the most Hamming-1
   neighbours in the table. On `training_c2.csv` the hub is `id41` with **50** neighbours against
   **4** for the runner-up -- not a tie, and not a judgement call. ⛔ If the margin is ever thin
   this script REFUSES rather than guessing, because a wrong reference silently redefines every
   mutation in the cycle (the same failure class as the WT2OUT_bridge2/bridge6 swap).

⛔ THE FLOOR APPLIES TO THE REFERENCE AND TO EVERY CANDIDATE. `id93` once topped the table at
   label 10.3333 on gate 3/10 and median sc_rmsd 1.871 -- it fails the fold gate 7 draws in 10.
   Anchoring on a thin-data outlier flatters the result: never rank on partial data.

⛔ ONE SUBSTITUTION PER POSITION. L49V and L49K are both accepted and mutually exclusive; a
   combination takes one or the other, never both.
⛔ NEVER RE-FOLD. Every sequence already in the training table is excluded -- it has been folded
   and scored, and re-folding it costs compute to learn nothing.

⭐ Enumeration order is LOWEST UNTRIED ORDER FIRST, then lexicographic by position tuple, so the
   output is deterministic and a cycle is reproducible from its training table alone. Cycle 2 was
   orders 2+3 (20+30 = 50); cycle 3 is orders 4+5+6 (25+11+2 = 38).

⭐⭐ THE REPLICATE RIDE-ALONG (`decisions/0045`, built 2026-10-05). A few of each cycle's slots are
   RESERVED for RE-FOLDS of the previous cycle's best designs, riding inside the normal cycle
   instead of launching a separate fleet. ⛔ This is the ONE place the "never re-fold" rule above
   is deliberately broken, and it is broken on purpose:
     - WHY: before this, the campaign had **never re-folded anything**. Every design carried
       exactly ONE measurement, so there was NO noise estimate at all, and both metrics have
       SATURATED (`label` at its 11-residue ceiling; `full_epitope_frac` with designs at 1.000).
       A single design's `fef` carries SE ~= 0.16 at 10 draws, so 1.000 vs 0.900 is inside the
       noise. Draws cannot be raised (97.6% VRAM measured, and Boltz OOM exits rc=0), so precision
       can only come from REPEATED 10-draw runs.
     - COST: **zero marginal compute** -- the replicates consume slots the new designs would
       have used, so the draw count is unchanged. It also controls for batch
       effects, since replicate and new designs are folded in the SAME run.
     - THE ANCHOR: one slot re-folds the HUB, whose prior is known on 20 draws. Without an anchor
       a replicate round can only confirm champions; with one it measures REGRESSION TO THE MEAN,
       which is the quantity that says how much of each cycle's apparent gain is selection.
   ⛔⛔ REPLICATES MUST NEVER ENTER `training_cN.csv` AS NEW DESIGNS (`0045` item 5). They are
   duplicate sequences and would corrupt the table the hub detection and the surrogate read. They
   are written to a SEPARATE `cycle<n>_replicates.csv`, which `score_cycle_egfrA1.py` reads to
   hold them out of the merge. ⛔ They are deliberately ABSENT from `cycle<n>_map.csv` -- that
   file defines the cycle's NEW DESIGNS and nothing else.
   ⛔ Identify a replicate ONLY from the replicates file, NEVER from its id. The id prefix "id9"
   already matches REAL cycle-0 designs (id9, id90, id91..id99, id901, id902), so prefix matching
   would silently hold out ten genuine burn-in designs.

usage: propose_cycle_egfrA1.py --training <pooled.csv> --cycle <n> --cap <n>
                               --out-tsv <f> --out-map <f>
                               [--out-replicates <f>] [--replicate-top 5]
                               [--replicate-anchor auto|none|<seq_id>]
                               [--design <name>] [--min-hub-margin 3]
"""
import argparse
import csv
import itertools
import sys

GATE_FLOOR_COL = "passes_floor"
MIN_HUB_NEIGHBOURS = 5          # a reference with fewer measured singles is not a walk
NEEDED = ("seq_id", "sequence", "label", "full_epitope_frac", GATE_FLOOR_COL)
# ⭐ replicate ids live at cycle*1000 + 900 + i, above any plausible new-design count in the same
#    cycle block, so the two never collide and a replicate id reads as one on sight. ⛔ That is a
#    CONVENIENCE, not an identifier -- the replicates FILE is the only authority (see docstring).
REPLICATE_ID_OFFSET = 900
REPLICATE_COLS = ("replicate_id", "source_id", "prior_fef", "prior_label",
                  "prior_median_sc", "prior_n_draws", "mutations")


def die(msg):
    sys.exit(f"⛔ {msg}")


def truthy_floor(v):
    return str(v).strip() in ("1", "True", "true")


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def diffs(ref, seq):
    """[(1-based pos, wt, sub)] where seq differs from ref."""
    return [(i + 1, ref[i], seq[i]) for i in range(len(ref)) if ref[i] != seq[i]]


def load(path):
    try:
        rows = list(csv.DictReader(open(path)))
    except OSError as e:
        die(f"cannot read training table {path}: {e}")
    if not rows:
        die(f"training table {path} is empty")
    missing = [c for c in NEEDED if c not in rows[0]]
    if missing:
        die(f"training table {path} lacks required columns {missing}; "
            f"has {sorted(rows[0])}. Refusing -- a missing label column would silently "
            f"produce an accepted set of nothing.")
    keep = [r for r in rows if r.get("sequence")]
    if not keep:
        die(f"no row in {path} carries a sequence")
    return rows, keep


def pick_hub(scored, min_margin):
    """The reference: the sequence with the most Hamming-1 neighbours. Refuses a thin margin."""
    lens = {len(r["sequence"]) for r in scored}
    if len(lens) != 1:
        die(f"training table mixes sequence lengths {sorted(lens)}; refusing to diff across them")
    counts = []
    for r in scored:
        n = sum(1 for q in scored
                if q is not r and hamming(r["sequence"], q["sequence"]) == 1)
        counts.append((n, r["seq_id"], r))
    counts.sort(key=lambda t: (-t[0], t[1]))
    (n1, id1, hub), = counts[:1]
    n2, id2 = (counts[1][0], counts[1][1]) if len(counts) > 1 else (0, "-")
    print(f"hub detection: {id1} has {n1} Hamming-1 neighbours; runner-up {id2} has {n2}")
    if n1 < MIN_HUB_NEIGHBOURS:
        die(f"hub {id1} has only {n1} measured single mutants (need >= {MIN_HUB_NEIGHBOURS}). "
            f"There is no measured walk to extend. Refusing.")
    if n1 - n2 < min_margin:
        die(f"hub is AMBIGUOUS: {id1} ({n1}) vs {id2} ({n2}), margin {n1-n2} < {min_margin}. "
            f"The reference defines what every mutation in this cycle MEANS, so a near-tie is "
            f"refused rather than broken by sort order.")
    if not truthy_floor(hub[GATE_FLOOR_COL]):
        die(f"hub {id1} FAILS the fold floor ({GATE_FLOOR_COL}={hub[GATE_FLOOR_COL]!r}, "
            f"gate {hub.get('n_gated')}/{hub.get('n_draws')}). Refusing to walk from a sequence "
            f"that does not reliably fold.")
    return hub


def accepted_set(hub, scored):
    """Single mutants of the hub that beat it on BOTH metrics and pass the floor."""
    H = hub["sequence"]
    hl, hf = float(hub["label"]), float(hub["full_epitope_frac"])
    print(f"reference {hub['seq_id']}: label {hl}  full_epitope_frac {hf}  "
          f"gate {hub.get('n_gated')}/{hub.get('n_draws')}")
    acc, rejected = {}, 0
    for r in scored:
        if r["seq_id"] == hub["seq_id"] or r["label"] == "":
            continue
        d = diffs(H, r["sequence"])
        if len(d) != 1:
            continue
        pos, wt, sub = d[0]
        if not truthy_floor(r[GATE_FLOOR_COL]):
            rejected += 1
            continue
        if float(r["label"]) > hl and float(r["full_epitope_frac"]) > hf:
            tok = f"{wt}{pos}{sub}"
            acc.setdefault(pos, []).append((tok, pos, wt, sub, r))
    flat = [t for v in acc.values() for t in v]
    print(f"accepted set: {len(flat)} mutations over {len(acc)} positions "
          f"{sorted(acc)}  ({rejected} single mutants rejected by the fold floor)")
    for pos in sorted(acc):
        for tok, _, _, _, r in sorted(acc[pos]):
            print(f"  ✅ {tok:7s} {r['seq_id']:8s} label {r['label']:>7} "
                  f"full_epitope_frac {r['full_epitope_frac']:>5}")
        if len(acc[pos]) > 1:
            print(f"  ⛔ position {pos} MUTUALLY EXCLUSIVE: "
                  f"{', '.join(sorted(t[0] for t in acc[pos]))}")
    if not flat:
        die("accepted set is EMPTY -- nothing beat the reference on both metrics. "
            "The walk has converged or the metric has saturated; that is a RESULT and needs a "
            "decision, not another cycle. Refusing to propose.")
    return acc


def mutation_string(hub, seq):
    """The hub-relative mutation list, in c90's format: `T19Q+S75A`, or `WT_hub` for the hub."""
    d = diffs(hub["sequence"], seq)
    return "+".join(f"{wt}{pos}{sub}" for pos, wt, sub in d) if d else "WT_hub"


def pick_replicates(scored, cycle, hub, n_top, anchor):
    """Reserve slots to RE-FOLD the previous cycle's best designs (`decisions/0045`).

    ⭐ SELECTION IS MECHANICAL so the reserved slots cannot drift into a judgement call:
    the pool is the PREVIOUS cycle's id block (`id{(cycle-1)*1000+1}`..`+999`), restricted to
    designs that PASS THE FOLD FLOOR and carry a label, ranked by `full_epitope_frac` descending,
    then `median_sc` ascending (fold quality breaks the ties the saturated metric leaves), then
    `seq_id` so the result is deterministic.
    ⛔ THE FLOOR IS APPLIED HERE TOO. Re-measuring a design that fails the fold gate measures the
    noise of a sequence we would not ship. Never rank on partial data.
    ⚠️ RETURNS FEWER THAN ASKED, LOUDLY, rather than dying, when the previous block is thin. A
    missing replicate costs precision; a dead proposer halts an authorized campaign.
    """
    prev = cycle - 1
    lo, hi = prev * 1000, prev * 1000 + 999

    def in_prev_block(r):
        t = r["seq_id"][2:] if r["seq_id"].startswith("id") else ""
        return t.isdigit() and lo <= int(t) <= hi

    def fef(r):
        try:
            return float(r["full_epitope_frac"])
        except (TypeError, ValueError):
            return -1.0

    def msc(r):
        try:
            return float(r["median_sc"])
        except (KeyError, TypeError, ValueError):
            return float("inf")

    out = []
    if n_top > 0:
        pool = [r for r in scored if in_prev_block(r) and r["label"] != ""
                and truthy_floor(r[GATE_FLOOR_COL])]
        pool.sort(key=lambda r: (-fef(r), msc(r), r["seq_id"]))
        if not pool:
            print(f"  ⚠️ NO floor-passing design in the cycle-{prev} id block "
                  f"(id{lo}..id{hi}) -- reserving NO champion replicate slots. The cycle still "
                  f"runs; it just carries no test-retest for this round.")
        elif len(pool) < n_top:
            print(f"  ⚠️ only {len(pool)} floor-passing designs in the cycle-{prev} block, "
                  f"asked for {n_top} -- taking all {len(pool)}")
        out = [dict(src=r, role="top") for r in pool[:n_top]]

    # ⭐ the ANCHOR: a known prior, so the round measures regression to the mean and not only
    #    whether champions repeat (`0045` item 3).
    if anchor != "none":
        want = hub["seq_id"] if anchor == "auto" else anchor
        hit = [r for r in scored if r["seq_id"] == want]
        if not hit:
            die(f"--replicate-anchor {want} is not in the training table. An anchor with no "
                f"prior measures nothing; name one that exists, or pass --replicate-anchor none.")
        if any(r["src"]["seq_id"] == want for r in out):
            print(f"  ⚠️ anchor {want} is ALREADY among the top replicates -- not duplicating it")
        else:
            out.append(dict(src=hit[0], role="anchor"))

    for r in out:
        s = r["src"]
        r["mutations"] = mutation_string(hub, s["sequence"])
    if out:
        print(f"replicate ride-along: {len(out)} slots reserved "
              f"({sum(1 for r in out if r['role']=='top')} top-of-cycle-{prev}, "
              f"{sum(1 for r in out if r['role']=='anchor')} anchor)")
        for r in out:
            s = r["src"]
            print(f"  ♻️ {s['seq_id']:8s} {r['role']:6s} fef {s['full_epitope_frac']:>5} "
                  f"label {s['label']:>7} med_sc {s.get('median_sc','?'):>5} "
                  f"n_draws {s.get('n_draws','?'):>3}  {r['mutations']}")
    return out


def enumerate_untried(hub, acc, known_seqs, cap):
    """Combinations, one substitution per position, lowest untried order first."""
    H = hub["sequence"]
    positions = sorted(acc)
    out, seen = [], set()
    for order in range(2, len(positions) + 1):
        made = skipped = 0
        for ps in itertools.combinations(positions, order):
            for choice in itertools.product(*[sorted(acc[p]) for p in ps]):
                s = list(H)
                for tok, pos, wt, sub, _ in choice:
                    s[pos - 1] = sub
                seq = "".join(s)
                if seq in known_seqs:
                    skipped += 1
                    continue
                if seq in seen:
                    skipped += 1
                    continue
                seen.add(seq)
                out.append(dict(order=order, n_mut=order, sequence=seq,
                                mutations="+".join(t[0] for t in choice),
                                positions="+".join(str(p) for p in ps)))
                made += 1
                if len(out) >= cap:
                    print(f"  order {order}: {made} new, {skipped} already folded "
                          f"-- CAP {cap} reached")
                    return out
        print(f"  order {order}: {made} new, {skipped} already folded")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training", required=True)
    ap.add_argument("--cycle", type=int, required=True)
    ap.add_argument("--cap", type=int, required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--design", default="arm_shard11_shard11_24_model_2")
    ap.add_argument("--min-hub-margin", type=int, default=3)
    ap.add_argument("--replicate-top", type=int, default=5,
                    help="slots reserved to RE-FOLD the previous cycle's best designs "
                         "(decisions/0045). 0 disables the ride-along entirely.")
    ap.add_argument("--replicate-anchor", default="auto",
                    help="'auto' (the hub -- the one sequence with a known prior), 'none', "
                         "or an explicit seq_id. Rides in ADDITION to --replicate-top.")
    ap.add_argument("--out-replicates", default=None,
                    help="replicate map. Default: --out-map with '_map.csv' -> '_replicates.csv', "
                         "which is where score_cycle_egfrA1.py looks for it unaided.")
    a = ap.parse_args()
    if a.cap <= 0:
        die(f"--cap must be positive, got {a.cap}")
    if a.replicate_top < 0:
        die(f"--replicate-top must be >= 0, got {a.replicate_top}")
    # ⛔ DERIVED BY CONVENTION, not passed by the driver. `egfr_backend.propose()` calls this
    #    script with --out-tsv/--out-map only, and `score_cycle_egfrA1.py` re-derives the same
    #    name from the map path it is handed. Keeping the convention in BOTH scripts is what lets
    #    the ride-along work through the autocycle driver with no change to the driver at all.
    repf = a.out_replicates
    if repf is None:
        if not a.out_map.endswith("_map.csv"):
            die(f"--out-map {a.out_map} does not end in '_map.csv', so the replicates filename "
                f"cannot be derived from it. Pass --out-replicates explicitly.")
        repf = a.out_map[:-len("_map.csv")] + "_replicates.csv"

    rows, scored = load(a.training)
    print(f"training table {a.training}: {len(rows)} rows, {len(scored)} with sequences")
    hub = pick_hub(scored, a.min_hub_margin)
    acc = accepted_set(hub, scored)

    reps = pick_replicates(scored, a.cycle, hub, a.replicate_top, a.replicate_anchor)
    cap_new = a.cap - len(reps)
    if cap_new <= 0:
        die(f"--cap {a.cap} leaves no room for new designs after {len(reps)} replicate slots. "
            f"A cycle that is ALL replicates is a re-measurement round, not a cycle -- run it "
            f"deliberately, not by arithmetic accident.")
    if cap_new > REPLICATE_ID_OFFSET:
        die(f"{cap_new} new-design slots would reach id{a.cycle*1000+cap_new}, colliding with "
            f"the replicate id block at id{a.cycle*1000+REPLICATE_ID_OFFSET+1}. Raise "
            f"REPLICATE_ID_OFFSET or lower --cap.")
    print(f"slots: {a.cap} total = {cap_new} new designs + {len(reps)} replicates")

    known = {r["sequence"] for r in scored}
    known.add(hub["sequence"])
    cand = enumerate_untried(hub, acc, known, cap_new)
    if not cand:
        die(f"every combination of the accepted set is ALREADY in the training table "
            f"({len(known)} sequences folded). The exhaustive enumeration is complete -- "
            f"the next move is a NEW accepted set from a new reference, which is a decision. "
            f"Refusing to propose an empty cycle.")

    # ---- ids: cycle-scoped block, asserted not to collide with anything folded ----
    idbase = a.cycle * 1000 + 1
    taken = {r["seq_id"] for r in rows}
    for i, r in enumerate(cand):
        r["fold_id"] = f"id{idbase + i}"
    for i, r in enumerate(reps):
        r["fold_id"] = f"id{a.cycle * 1000 + REPLICATE_ID_OFFSET + 1 + i}"
    clash = sorted(({r["fold_id"] for r in cand} | {r["fold_id"] for r in reps}) & taken)
    if clash:
        die(f"proposed ids collide with the training table: {clash[:8]}. "
            f"An id collision makes the merge pool two different sequences under one key.")

    # ---- assertions: each sequence differs from the reference in exactly n_mut positions ----
    H = hub["sequence"]
    for r in cand:
        d = hamming(H, r["sequence"])
        if d != r["n_mut"]:
            die(f"{r['fold_id']} {r['mutations']}: {d} diffs from the reference, "
                f"expected {r['n_mut']}")
    if len({r["sequence"] for r in cand}) != len(cand):
        die("duplicate sequences in the proposal")
    if len({r["fold_id"] for r in cand}) != len(cand):
        die("duplicate fold_ids in the proposal")
    if known & {r["sequence"] for r in cand}:
        die("a proposed sequence is already in the training table -- the never-re-fold "
            "exclusion failed")
    print(f"✅ {len(cand)} sequences, all distinct, none previously folded, "
          f"each differing from {hub['seq_id']} in exactly n_mut positions")

    # ⛔ THE MIRROR ASSERT. A replicate whose sequence is NOT already folded is not a replicate --
    #    it is a new design that has slipped into the re-measurement slots, and it would be
    #    compared against a "prior" that was never measured.
    bad = [r for r in reps if r["src"]["sequence"] not in known]
    if bad:
        die(f"replicate slots hold sequences never folded before: "
            f"{[r['src']['seq_id'] for r in bad]}. A replicate re-measures an EXISTING "
            f"measurement; there is nothing to compare these against.")
    if len({r["fold_id"] for r in reps}) != len(reps):
        die("duplicate fold_ids among the replicates")
    if len({r["src"]["seq_id"] for r in reps}) != len(reps):
        die("the same source design was reserved twice -- that doubles one design's weight in "
            "the test-retest and measures nothing extra")
    if reps:
        print(f"✅ {len(reps)} replicates, each a sequence ALREADY in the training table, "
              f"each re-folded under a fresh id")

    with open(a.out_tsv, "w") as fh:
        for r in cand:
            fh.write(f"A_{a.design}_{r['fold_id']}\t{r['sequence']}\n")
        # ⭐ replicates ride in the SAME TSV, so build_yamls_from_skeleton.py spreads them
        #    round-robin across the shards with the new designs -- that is the batch-effect
        #    control, and it is why this is a ride-along and not a second fleet.
        for r in reps:
            fh.write(f"A_{a.design}_{r['fold_id']}\t{r['src']['sequence']}\n")
    # ⭐ `mutation` (singular) is emitted ALONGSIDE `mutations` because score_cycle_egfrA1.py's
    #    walk printout reads the singular key; `arm`/`pred_label` are deliberately ABSENT -- there
    #    is no surrogate and no UCB/random split in a combination cycle, and the scorer skips its
    #    head-to-head block when they are missing rather than inventing a comparison.
    cols = ["fold_id", "order", "n_mut", "mutation", "mutations", "positions", "reference"]
    with open(a.out_map, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in cand:
            w.writerow({"fold_id": r["fold_id"], "order": r["order"], "n_mut": r["n_mut"],
                        "mutation": r["mutations"], "mutations": r["mutations"],
                        "positions": r["positions"], "reference": hub["seq_id"]})
    # ⛔ The replicates go in their OWN file, never into cycle<n>_map.csv: that map defines the
    #    cycle's NEW DESIGNS, and score_cycle_egfrA1.py holds out exactly what this file names.
    with open(repf, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(REPLICATE_COLS))
        w.writeheader()
        for r in reps:
            s_ = r["src"]
            w.writerow({"replicate_id": r["fold_id"], "source_id": s_["seq_id"],
                        "prior_fef": s_["full_epitope_frac"], "prior_label": s_["label"],
                        "prior_median_sc": s_.get("median_sc", ""),
                        "prior_n_draws": s_.get("n_draws", ""),
                        "mutations": r["mutations"]})

    byorder = {}
    for r in cand:
        byorder[r["order"]] = byorder.get(r["order"], 0) + 1
    print(f"wrote {a.out_tsv}, {a.out_map} and {repf} ({len(reps)} replicates)")
    print(f"cycle {a.cycle}: {len(cand)} new sequences, orders "
          f"{', '.join(f'{k}-way x{v}' for k, v in sorted(byorder.items()))}, "
          f"reference {hub['seq_id']}"
          + (f"; plus {len(reps)} replicates = {len(cand)+len(reps)} folds"
             if reps else "; NO replicates"))


if __name__ == "__main__":
    main()
