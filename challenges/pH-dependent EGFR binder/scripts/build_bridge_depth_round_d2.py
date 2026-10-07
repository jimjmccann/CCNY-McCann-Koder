"""Build the egfrBR50 DEPTH round 2 (d2) -- order the TOP of the bridge table, not the ceiling pair.

⭐⭐ WHY THIS EXISTS, AND WHY IT IS NOT build_bridge_depth_round.py.
   d1 deep-tested `id85` and `id1104` because they TIED at the label ceiling (both 1.000 on 3
   draws). The d1 round settled that pair -- id85 0.936 vs id1104 0.868,
   separating at 3.1 s.e. on the label and 71/90 vs 26/90 on `full_both` -- and left exactly one
   open item: **"id85 is the campaign's best" is untested**, because only two sequences were ever
   deep-tested and they were chosen for being tied, not for being good.
   ⇒ This round deep-tests the two sequences that look
   most like id85 ON THE MEASURE THAT ACTUALLY SEPARATED (`full_both`, not the label), with id85
   riding along as the IN-BATCH ANCHOR so the three are comparable without crossing rounds.

⭐ THE PICKS, and why they are not the label leaders.
   `id47` and `id76` each hit `full_both` on 2 of 3 draws with clean folds (med sc_rmsd 0.837,
   0.805). The LABEL leaders `id43`/`id59` (0.905) produced ZERO `full_both` on retest -- so the
   campaign's own shortlist, which was ranked on the label, may have deep-tested the wrong pair.
   ⛔ Neither id47 nor id76 has ever been deep-tested.

⛔ THE SEQUENCE SOURCE IS THE POOLED TABLE, NOT bridge_cycle1.tsv. Measured 2026-10-06:
   `bridge_cycle1.tsv` holds 50 rows -- the c1 NEW designs -- and contains `id85` but NOT `id47`
   or `id76`, which are c1 LAND rows. Both live in `label_BR_0046_chosen_POOLED100.csv`, which
   carries a `sequence` column. That is why d1's builder cannot be reused: its `--c1-tsv` input
   is missing two of this round's three sequences. All three are on the SAME backbone
   (`arm_shard11_shard11_41_model_5`), which this script asserts rather than assumes.

⛔ EVERY GUARANTEE OF THE d1 BUILDER IS KEPT, because the scorer depends on them:
   - each replicate gets its OWN id  (score_br50.py keys `per` on the sequence id, so N duplicate
     entries under N distinct ids score as N INDEPENDENT 3-draw labels -- the mean is the deep
     estimate AND the spread is a direct noise measurement; pooling under one id throws the
     spread away)
   - EXACTLY 3 draws per entry (`score_br50.py` refuses a set whose per-sequence draw counts
     differ; depth comes from REPEATING the fold, never from raising draws)
   - refuse on id collision against every id the campaign has used
   - refuse on an existing --out (scripts APPEND, so a reused path merges two campaigns)
   - ⛔ score SEPARATELY from the c1 table (`decisions/0045` item 5; `score_br50.py` enforces it
     by refusing a sequence that appears under two campaign tags)

⭐ SIZING, from d1's own noise measurement rather than from c1's. d1 found the noise is
   HETEROSCEDASTIC: the leaders are quiet (per-replicate sd 0.069-0.100) and the mid-table is loud
   (0.250-0.271). c1's single 0.189 figure, and the ~86-draw guidance built on it, is wrong in
   both directions. At the top, ~12-24 draws orders a pair. 12 replicates = 36 draws carries
   margin over that; the anchor needs only enough to confirm it reproduces its own 0.936.

usage: build_bridge_depth_round_d2.py --pooled <label_BR_0046_chosen_POOLED100.csv>
                                      --out-tsv <d2.tsv> --out-map <d2_map.csv>
                                      [--deep 12] [--anchor 8] [--also-used <prior.tsv> ...]
"""
import argparse, csv, os, re, sys

BACKBONE = "arm_shard11_shard11_41_model_5"
# (source id, kind, why it is here, new id base)
DEEP = [
    ("47", "deep", "2/3 full_both, med sc_rmsd 0.837, NEVER deep-tested; most id85-like on the "
                   "measure that separated", 3000),
    ("76", "deep", "2/3 full_both, med sc_rmsd 0.805, NEVER deep-tested", 3100),
]
ANCHORS = [
    ("85", "anchor", "IN-BATCH ANCHOR: d1 deep value 0.936 over 90 draws, full_both 71/90. A "
                     "batch shift shows up as id85 missing its own d1 number", 3200),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pooled", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--deep", type=int, default=12)
    ap.add_argument("--anchor", type=int, default=8)
    ap.add_argument("--also-used", nargs="*", default=[],
                    help="prior round .tsv files whose ids must also be treated as USED")
    a = ap.parse_args()

    for p in (a.out_tsv, a.out_map):
        if os.path.exists(p):
            sys.exit(f"⛔ {p} exists. Fresh --out per experiment. Refusing.")

    seqs, used, meta = {}, set(), {}
    for r in csv.DictReader(open(a.pooled)):
        if r["design"] != BACKBONE:
            continue
        sid = r["seq"].strip()
        seqs[sid] = r["sequence"].strip()
        meta[sid] = {"label_BR": r["label_BR"], "n_full_both": r["n_full_both"],
                     "med_sc_rmsd": r["med_sc_rmsd"], "campaign": r["campaign"]}
        used.add(int(sid))
    if not seqs:
        sys.exit(f"⛔ no rows for backbone {BACKBONE} in {a.pooled}")

    for p in a.also_used:
        for line in open(p):
            if not line.strip():
                continue
            m = re.match(r"^.*_id(\d+)\t", line)
            if not m:
                sys.exit(f"⛔ cannot parse an id out of {line.split(chr(9))[0]!r} in {p}")
            used.add(int(m.group(1)))

    rows, mapped = [], []
    for sid, kind, why, base in DEEP + ANCHORS:
        n = a.deep if kind == "deep" else a.anchor
        if sid not in seqs:
            sys.exit(f"⛔ id{sid} is not in {a.pooled} on backbone {BACKBONE}. Refusing to guess.")
        seq = seqs[sid]
        for i in range(n):
            nid = base + i
            if nid in used:
                sys.exit(f"⛔⛔ id{nid} IS ALREADY USED. A collision makes the scorer sys.exit "
                         f"AFTER the compute is already consumed. Refusing.")
            used.add(nid)
            nm = f"BR_{BACKBONE}_id{nid}"
            rows.append((nm, seq))
            mapped.append({"name": nm, "id": nid, "replicate_index": i, "depth_of_id": sid,
                           "kind": kind, "prior_label_BR": meta[sid]["label_BR"],
                           "prior_n_full_both": meta[sid]["n_full_both"],
                           "prior_med_sc_rmsd": meta[sid]["med_sc_rmsd"],
                           "prior_campaign": meta[sid]["campaign"], "why": why})

    with open(a.out_tsv, "w") as fh:
        for nm, seq in rows:
            fh.write(f"{nm}\t{seq}\n")
    with open(a.out_map, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(mapped[0].keys()))
        w.writeheader(); w.writerows(mapped)

    print(f"wrote {a.out_tsv}: {len(rows)} fold entries x 3 draws = {len(rows)*3} structures")
    print(f"wrote {a.out_map}")
    for sid, kind, why, base in DEEP + ANCHORS:
        n = a.deep if kind == "deep" else a.anchor
        tag = "DEEP  " if kind == "deep" else "anchor"
        print(f"  {tag} id{sid:<5} -> ids {base}..{base+n-1}  {n} replicates = {n*3} draws "
              f"(prior label {meta[sid]['label_BR'][:5]}, full_both {meta[sid]['n_full_both']}/3)")
    print(f"  all {len(rows)} sequences are {len(rows[0][1])} aa, "
          f"distinct lengths seen: {len({len(s) for _, s in rows})}")
    print(f"  distinct MOLECULES in the round: {len(set(s for _, s in rows))} "
          f"(duplicates are deliberate -- that IS the depth)")


if __name__ == "__main__":
    main()
