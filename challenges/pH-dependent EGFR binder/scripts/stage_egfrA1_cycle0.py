#!/usr/bin/env python3
"""Stage the egfrA1 CYCLE-0 head-to-head TSV: which sequence does the walk start from?

Four candidates, K=10 draws each = 40 structures, ONE shard, one worker:
  consensus_top10   the selected mode -- consensus of the 10 best-scoring designs
  consensus_all100  the all-100 consensus, for contrast (it differs at only 8 of 75 positions)
  id41              the best MEASURED sequence under the loss (label 10.30)
  id15              the best sequence by the BINARY full-epitope rate (6/10 draws at 11/11)

⛔ WHY THIS CYCLE EXISTS. Both consensus sequences are NOVEL -- neither is one of the 100, and
neither has ever been folded. Starting a long (24-cycle) campaign from an unverified sequence 17
mutations away from the only good one we have measured is an avoidable risk for ~35 minutes.
⭐ And it matters MORE than it would otherwise, because the cycle-0 learnability check came back
empty (ridge R2 -0.060 / RF -0.052): with the surrogate uninformative at the start, the starting
point carries more of the weight.

⛔ Consensus sequences are RECOMPUTED here from the burn-in, never copied from a previous printout.

usage: stage_egfrA1_cycle0.py <burnin.csv> <out.tsv>
"""
import csv, collections, sys


def consensus(seqs):
    L = len(seqs[0])
    return "".join(collections.Counter(s[i] for s in seqs).most_common(1)[0][0] for i in range(L))


def main(burnin, outp):
    rows = list(csv.DictReader(open(burnin)))
    ok = [r for r in rows if r["passes_floor"] == "1"]
    ok.sort(key=lambda r: -float(r["label"]))
    by = {r["seq_id"]: r["sequence"] for r in rows}
    L = len(ok[0]["sequence"])

    c10 = consensus([r["sequence"] for r in ok[:10]])
    c100 = consensus([r["sequence"] for r in rows])          # all 100, floor or not
    cands = [("consensus_top10", c10), ("consensus_all100", c100),
             ("id41", by["id41"]), ("id15", by["id15"])]

    seen = {}
    for n, s in cands:
        assert len(s) == L, f"⛔ {n} is {len(s)} aa, expected {L}"
        if s in seen:
            print(f"⚠️ {n} is IDENTICAL to {seen[s]} -- folding it twice would waste a slot")
        seen[s] = n
    known = {r["sequence"]: r["seq_id"] for r in rows}
    for n, s in cands:
        tag = known.get(s)
        d41 = sum(1 for a, b in zip(s, by["id41"]) if a != b)
        print(f"  {n:18s} {'(= ' + tag + ')' if tag else '(NOVEL, never folded)':24s} "
              f"{d41:2d} mutations from id41")

    with open(outp, "w") as fh:
        for n, s in cands:
            fh.write(f"A_egfrA1_c0_{n}\t{s}\n")
    print(f"\nwrote {outp}  ({len(cands)} sequences x 10 draws = {10*len(cands)} structures, 1 shard)")
    print(f"top-10 consensus:\n  {c10}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
