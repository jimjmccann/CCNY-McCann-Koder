#!/usr/bin/env python3
"""Build a COMBINATORIAL egfrA1 cycle from an accepted set of single mutations.

⭐ WHY. Cycle 1 (50 single mutants of id41) ended with SEVEN mutations that beat the parent on
BOTH metrics and were mutually indistinguishable -- pairwise Welch dmean +0.00, Fisher p=0.628-1.000
(measured). The `label` had also SATURATED: the site-A epitope
has 11 residues, every leader draw scored 10 or 11, so the mean can no longer rank candidates.
Picking one winner would have been arbitrary. The chosen approach: use ALL of the
indistinguishable top scorers as the baseline for allowed mutations.

⇒ So this enumerates COMBINATIONS of the accepted set instead of single mutants of one parent.
⛔ It is NOT a surrogate proposal. The ensemble has no skill here -- OOB RMSE / sd(y) was 1.14 on
   `label` and 1.04 on `full_epitope_frac`, i.e. at or worse than predicting the mean, and UCB did
   not beat random (+0.394 +/- 0.634, t=+0.62). Exhaustive coverage needs no ranking, so for this
   cycle the model is out of the loop entirely. That is a measured decision, not a preference.

⛔ ONE SUBSTITUTION PER POSITION. L49V and L49K are both accepted and mutually exclusive; a
   combination may take one or the other, never both.
⛔ EVERY MUTATION'S WILD-TYPE RESIDUE IS ASSERTED against the parent sequence before use. A silent
   off-by-one in the position index would produce 50 plausible-looking wrong sequences.

usage: build_combo_cycle_egfrA1.py <parent.txt> <out.tsv> <out_map.csv>
         [--orders 2,3] [--design <name>] [--idbase 2001]
"""
import argparse, csv, itertools, re, sys

# the accepted set: beat parent id41 on BOTH label and full_epitope_frac, cycle 1
ACCEPTED = ["L16D", "A17M", "T19Q", "L49V", "L49K", "R55S", "S75A"]
MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def parse(m):
    g = MUT.match(m)
    if not g:
        sys.exit(f"⛔ bad mutation token {m!r}")
    return g.group(1), int(g.group(2)), g.group(3)


def apply_muts(seq, muts):
    s = list(seq)
    for wt, pos, sub in muts:
        s[pos - 1] = sub
    return "".join(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parent"); ap.add_argument("out_tsv"); ap.add_argument("out_map")
    ap.add_argument("--orders", default="2,3")
    ap.add_argument("--design", default="arm_shard11_shard11_24_model_2")
    ap.add_argument("--idbase", type=int, default=2001)
    a = ap.parse_args()

    parent = open(a.parent).read().split()[0]
    print(f"parent {len(parent)} aa")

    # ---- assert every accepted mutation against the parent, then group by position ----
    bypos = {}
    for m in ACCEPTED:
        wt, pos, sub = parse(m)
        if not 1 <= pos <= len(parent):
            sys.exit(f"⛔ {m}: position {pos} outside 1..{len(parent)}")
        got = parent[pos - 1]
        if got != wt:
            sys.exit(f"⛔ {m}: parent has {got!r} at position {pos}, mutation claims {wt!r}. "
                     f"REFUSING -- an off-by-one here silently corrupts every sequence.")
        bypos.setdefault(pos, []).append((wt, pos, sub))
        print(f"  ✅ {m}: parent[{pos}] == {wt}")
    positions = sorted(bypos)
    print(f"accepted: {len(ACCEPTED)} mutations over {len(positions)} positions {positions}")
    for p in positions:
        if len(bypos[p]) > 1:
            print(f"  ⛔ position {p} is MUTUALLY EXCLUSIVE: "
                  f"{', '.join(w+str(p)+s for w,p,s in bypos[p])}")

    orders = [int(x) for x in a.orders.split(",")]
    rows, seen = [], {parent: "PARENT"}
    for n in orders:
        made = 0
        for ps in itertools.combinations(positions, n):
            for choice in itertools.product(*[bypos[p] for p in ps]):
                seq = apply_muts(parent, choice)
                tag = "+".join(f"{w}{p}{s}" for w, p, s in choice)
                if seq in seen:
                    print(f"  ⚠️ duplicate sequence {tag} == {seen[seq]} -- skipped"); continue
                seen[seq] = tag
                rows.append(dict(order=n, mutations=tag, n_mut=n, sequence=seq,
                                 positions="+".join(str(p) for p in ps)))
                made += 1
        print(f"  {n}-way: {made}")

    for i, r in enumerate(rows):
        r["fold_id"] = f"id{a.idbase + i}"

    # every sequence must differ from the parent in exactly n_mut positions
    for r in rows:
        d = sum(1 for x, y in zip(parent, r["sequence"]) if x != y)
        assert d == r["n_mut"], f"⛔ {r['fold_id']} {r['mutations']}: {d} diffs, expected {r['n_mut']}"
    assert len({r["sequence"] for r in rows}) == len(rows), "⛔ duplicate sequences"
    assert len({r["fold_id"] for r in rows}) == len(rows), "⛔ duplicate ids"
    print(f"✅ {len(rows)} sequences, all distinct, each differing from the parent in exactly n_mut positions")

    with open(a.out_tsv, "w") as fh:
        for r in rows:
            fh.write(f"A_{a.design}_{r['fold_id']}\t{r['sequence']}\n")
    with open(a.out_map, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fold_id", "order", "n_mut", "mutations", "positions"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in w.fieldnames})
    print(f"wrote {a.out_tsv} and {a.out_map}")


if __name__ == "__main__":
    main()
