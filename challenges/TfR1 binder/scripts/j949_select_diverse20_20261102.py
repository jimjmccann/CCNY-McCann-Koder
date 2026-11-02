#!/usr/bin/env python3
"""j949 -- pick a DIVERSE subset of the N-term arm's MPNN sequences for the final co-fold.

    RUN WITH: python3 (stdlib only, any py3 works)

Local, reads the j744 FASTAs and writes ONE TSV in the schema `j748` consumes
(`design`, `seq_id`, `sequence`, plus the carried confidences). ⛔ Designs nothing, folds nothing.

WHY DIVERSITY AND NOT TOP-N BY CONFIDENCE
-----------------------------------------
`[F]` THE ASK, 2026-11-02: *"Take whatever 20 best sequences we have at this stage (choose a diverse
sequence set among the best backbones?)"* — and the standing reasoning behind it, from `j918`'s
header: *"a diversity of sequences was asked for, on the reasoning that only one has to bind."*
⇒ ⭐ ranking 30 sequences by MPNN confidence and taking the top 20 would concentrate the set on one
basin; **only one has to bind**, so spread is worth more than mean quality.
⛔ AND MPNN CONFIDENCE IS NOT A BINDING PREDICTOR: the sequence-collection analysis measured that it does not even separate
the ARMS, and it was re-measured on the biased cohort with the same result. Using it as the sole
rank would be ranking on a quantity already shown not to carry the signal. It is therefore used only
to SEED the selection, never to filter.

THE RULE, STATED SO IT CAN BE CHECKED
-------------------------------------
Hamming distance over the **27 DESIGNED positions ONLY** (1-25, 27-28). ⛔ The 143 held positions
are identical in every sequence by construction, so including them would divide every distance by
the same constant and make a 2-mutation difference look like 98.8% identity — a diversity number
that is really a measure of how much was held.
Selection = **farthest-point (max-min)**: seed with the single highest `overall_confidence`, then
repeatedly add the candidate whose MINIMUM Hamming distance to the already-chosen set is largest.
⇒ deterministic, no RNG, and it maximises the smallest gap in the chosen set rather than the mean
([[rank-on-conjunction-not-mean-of-fractions]] — a mean pays partial credit).

⛔⛔ NEVER `cat *.fa | grep -c '^>'`. `[M]` MEASURED: LigandMPNN `.fa` files have NO trailing newline,
so `cat` fuses each file's last SEQUENCE line onto the next file's `>` header and the header stops
being at line start — it undercounted by exactly (n_files - 1). ⇒ every file is parsed ON ITS OWN
and its record count is asserted.

⛔ The FASTA holds **ALL CHAINS**, `:`-separated [[mpnn-homo-oligomer-copies-target]]. Only the
FIRST record of each file is the input/native (no `id=`); the designed records carry `id=`.

usage:
  j949_select_diverse20_20261102.py --fasta-dir <dir of *.fa> --fixed-json <fixed_residues_multi.json>
      --n 20 --out <tsv>
"""

import argparse
import csv
import glob
import json
import os
import re
import sys


def parse_fasta(path):
    """Return (native_seq_chainA, [(seq_id, chainA, conf, seq_rec), ...]) for ONE file."""
    txt = open(path).read()
    recs = txt.split(">")
    recs = [r for r in recs if r.strip()]
    if not recs:
        sys.exit("⛔ %s parsed to ZERO records -- an empty parse is not an empty file" % path)
    native, designed = None, []
    for r in recs:
        lines = r.strip().split("\n")
        head, body = lines[0], "".join(lines[1:])
        chain_a = body.split(":")[0].strip()
        if not chain_a:
            sys.exit("⛔ %s: a record has no chain-A segment before ':'" % path)
        if "id=" not in head:
            if native is not None:
                sys.exit("⛔ %s: more than one record without id= -- cannot tell native from design"
                         % path)
            native = chain_a
            continue
        sid = re.search(r"id=(\d+)", head)
        conf = re.search(r"overall_confidence=([0-9.]+)", head)
        srec = re.search(r"seq_rec=([0-9.]+)", head)
        designed.append((sid.group(1) if sid else "?", chain_a,
                         float(conf.group(1)) if conf else float("nan"),
                         float(srec.group(1)) if srec else float("nan")))
    if native is None:
        sys.exit("⛔ %s: no native/input record found (every record carried id=)" % path)
    return native, designed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--fixed-json", required=True,
                    help="j744's fixed_residues_multi.json -- defines the DESIGNED positions")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--out", required=True, help="⛔ fresh; refused if it exists")
    a = ap.parse_args()

    if os.path.exists(a.out):
        sys.exit("⛔ %s exists. Scripts here APPEND; use a fresh --out." % a.out)

    # ---- the designed position set, from j744's own fixed list (never hardcoded)
    fx = json.load(open(a.fixed_json))
    fixed_sets = {frozenset(int(t[1:]) for t in v.split()) for v in fx.values()}
    if len(fixed_sets) != 1:
        sys.exit("⛔ the designs do not share one fixed set (%d distinct) -- a single designed-position "
                 "index would be wrong for some of them. Refusing." % len(fixed_sets))
    fixed = set(next(iter(fixed_sets)))

    files = sorted(glob.glob(os.path.join(a.fasta_dir, "*.fa")))
    if not files:
        sys.exit("⛔ no *.fa under %s -- an empty glob is not an empty cohort" % a.fasta_dir)

    pool, lens, natives = [], set(), {}
    for f in files:
        design = os.path.basename(f)[:-3]
        native, designed = parse_fasta(f)
        natives[design] = native
        lens.add(len(native))
        for sid, seq, conf, srec in designed:
            lens.add(len(seq))
            pool.append({"design": design, "seq_id": sid, "sequence": seq,
                         "overall_confidence": conf, "seq_rec": srec})
        print("  %-44s %2d designed records" % (design, len(designed)))

    if len(lens) != 1:
        sys.exit("⛔ chain-A lengths differ across the cohort (%s) -- Hamming over a fixed index "
                 "would be meaningless. Refusing." % sorted(lens))
    L = next(iter(lens))
    designed_pos = sorted(set(range(1, L + 1)) - fixed)
    print("\n  chain-A length %d | held %d | ⭐ DESIGNED %d: %s"
          % (L, len(fixed), len(designed_pos), designed_pos))
    if not designed_pos:
        sys.exit("⛔ zero designed positions -- there is nothing to be diverse over.")

    # ---- assert the held positions really are constant, then ignore them
    idx = [p - 1 for p in designed_pos]
    held_idx = [p - 1 for p in sorted(fixed)]
    for r in pool:
        for d, nat in natives.items():
            if r["design"] != d:
                continue
            bad = [sorted(fixed)[i] for i, h in enumerate(held_idx)
                   if r["sequence"][h] != nat[h]]
            if bad:
                sys.exit("⛔⛔ %s id=%s changed HELD position(s) %s -- the fixed set was not honoured."
                         % (d, r["seq_id"], bad[:10]))
    print("  ✅ every held position identical to its backbone's native, across %d sequences"
          % len(pool))

    for r in pool:
        r["_v"] = "".join(r["sequence"][i] for i in idx)

    uniq = len({r["_v"] for r in pool})
    print("  distinct designed-region variants: %d of %d" % (uniq, len(pool)))

    n = min(a.n, len(pool))
    ham = lambda x, y: sum(1 for p, q in zip(x, y) if p != q)

    # ---- farthest-point: seed on best confidence, then max-min
    pool.sort(key=lambda r: (-(r["overall_confidence"]
                               if r["overall_confidence"] == r["overall_confidence"] else -1),
                             r["design"], int(r["seq_id"])))
    chosen = [pool[0]]
    rest = pool[1:]
    while len(chosen) < n and rest:
        best, bi = None, None
        for i, c in enumerate(rest):
            d = min(ham(c["_v"], s["_v"]) for s in chosen)
            if best is None or d > best:
                best, bi = d, i
        chosen.append(rest.pop(bi))

    # ---- report what the selection bought, measured
    def pair_stats(rows):
        ds = [ham(rows[i]["_v"], rows[j]["_v"])
              for i in range(len(rows)) for j in range(i + 1, len(rows))]
        if not ds:
            return (0, 0, 0.0)
        return (min(ds), max(ds), sum(ds) / len(ds))

    pmin, pmax, pavg = pair_stats(chosen)
    amin, amax, aavg = pair_stats(pool if len(pool) <= 60 else pool[:60])
    print("\n  SELECTED %d of %d" % (len(chosen), len(pool)))
    print("    pairwise Hamming over the %d designed positions:" % len(designed_pos))
    print("      chosen set : min %d  max %d  mean %.2f" % (pmin, pmax, pavg))
    print("      whole pool : min %d  max %d  mean %.2f" % (amin, amax, aavg))
    print("    ⭐ min pairwise distance is the number farthest-point maximises -- "
          "chosen %d vs pool %d" % (pmin, amin))
    perd = {}
    for r in chosen:
        perd[r["design"]] = perd.get(r["design"], 0) + 1
    print("    per backbone: %s" % dict(sorted(perd.items())))
    cf = [r["overall_confidence"] for r in chosen if r["overall_confidence"] == r["overall_confidence"]]
    if cf:
        print("    overall_confidence of chosen: min %.4f  mean %.4f  max %.4f"
              % (min(cf), sum(cf) / len(cf), max(cf)))

    # ⛔ `length` is REQUIRED by j748, which cross-checks it against len(sequence) rather than
    #    trusting either alone. It is written from the sequence itself, so the two can only
    #    disagree if this file is hand-edited -- which is exactly what that check is for.
    cols = ["design", "seq_id", "sequence", "length", "overall_confidence", "seq_rec"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in chosen:
            r["length"] = len(r["sequence"].strip())
            w.writerow(r)
    print("\n  ⇒ %s" % a.out)
    print("  ⛔ no verdict column; MPNN confidence SEEDED the selection and did not filter it.")


if __name__ == "__main__":
    main()
