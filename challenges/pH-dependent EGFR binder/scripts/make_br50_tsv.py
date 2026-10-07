"""egfrBR50 foldset: 8 shortlisted bridge backbones x 50 NEW sequences, ranked by MPNN confidence.

⭐ WHY ids 1-10 ARE EXCLUDED. Those are exactly the sequences egfrBR3 ALREADY FOLDED, at 3 draws.
   This round is also 3 draws, so re-folding them buys nothing -- it would be a duplicate
   campaign folded twice. ⇒ select from ids 11-100 only; combined with the existing data each backbone
   then has 10 (held) + 50 (new) = 60 sequences at 3 draws.

⭐⭐ WHY RANK ON LigandMPNN `overall_confidence`, MEASURED ON OUR OWN 1,082 SCORED SEQUENCES
   rather than assumed. Spearman vs outcome:
       confidence vs median sc_rmsd  rho = -0.299  (p=8e-24)   lower sc is better => PREDICTIVE
       confidence vs best site A     rho = -0.137  (p=6e-06)   mildly ANTI-correlated
       confidence vs best site B     rho = -0.036  (p=0.23)    nothing
   ⛔ Taken at face value that looks like a reason NOT to use confidence: the LOW-confidence half
   holds MORE full-site-A sequences (172 vs 119). ⛔⛔ THAT IS FOLD-ARTIFACT CONTACT. After the
   sc_rmsd gate the ranking REVERSES on every axis:
       top half by confidence : 410 admissible, 78 gated full-A, 11 gated full-B, 3 BOTH
       bottom half            : 302 admissible, 59 gated full-A,  5 gated full-B, 0 BOTH
   All three known bridges sit in the high-confidence half. ⇒ confidence is the right ranker, but
   ONLY because the gate is applied; ungated it points the wrong way
   (design-vs-prediction-selfconsistency).
   ⚠️ `ligand_confidence` is byte-identical to `overall_confidence` in these files -- it adds nothing.

⛔ TWO FORMAT TRAPS, same guards as the 100-sequence builder:
   record 0 is the INPUT sequence (no `id=`), and the line is `BINDER:TARGET_CONTEXT` on ONE line,
   so a naive read yields 355-375 aa instead of 130-150.
⭐ The extraction CODE PATH is still cross-checked against sequences that actually folded (ids 1-10
   vs bridge118x10.tsv) even though those ids are then excluded -- that is what validates the
   extraction of ids 11-100, which have no ground truth of their own.

usage: make_br50_tsv.py <mpnn_seqs_dir> <reference_tsv> <out.tsv> [n_per_backbone]
"""
import sys, os, re, glob, collections

BACKBONES = [
    "arm_shard07_shard07_10_model_6",
    "arm_shard11_shard11_41_model_5",
    "arm_shard07_shard07_23_model_3",
    "arm_shard12_shard12_30_model_7",
    "arm_shard06_shard06_28_model_7",
    "arm_shard10_shard10_37_model_1",
    "arm_shard01_shard01_2_model_0",
    "arm_shard13_shard13_10_model_0",
]
CTX_LEN = 224
HELD_IDS = set(range(1, 11))      # already folded by egfrBR3 at 3 draws


def read_fa(path):
    """-> {id: (binder, confidence)}; asserts the colon context and skips record 0."""
    out, cur, conf = {}, None, None
    ctx = set()
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            m = re.search(r"\bid=(\d+)", line)
            c = re.search(r"overall_confidence=([0-9.]+)", line)
            cur = int(m.group(1)) if m else None
            conf = float(c.group(1)) if c else None
            continue
        if not line or cur is None:
            continue
        if ":" not in line:
            sys.exit(f"⛔ {os.path.basename(path)} id={cur}: no ':' -- format changed; refusing")
        binder, _, context = line.partition(":")
        ctx.add(context)
        if conf is None:
            sys.exit(f"⛔ {os.path.basename(path)} id={cur}: no overall_confidence to rank on")
        out[cur] = (binder, conf)
    if len(ctx) != 1:
        sys.exit(f"⛔ {os.path.basename(path)}: {len(ctx)} distinct contexts, expected 1")
    c = ctx.pop()
    if len(c) != CTX_LEN:
        sys.exit(f"⛔ {os.path.basename(path)}: context {len(c)} aa, expected {CTX_LEN}")
    return out


def main(seqdir, reftsv, outp, n_per=50):
    ref = {}
    for line in open(reftsv):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            ref[p[0]] = p[1]

    rows, checked, mismatch = [], 0, []
    print(f"{'backbone':40s} {'binder':>7} {'conf range picked':>22} {'held':>5}")
    for bb in BACKBONES:
        f = os.path.join(seqdir, bb + ".fa")
        if not os.path.exists(f):
            sys.exit(f"⛔ missing {f}")
        seqs = read_fa(f)
        if len(seqs) != 100:
            sys.exit(f"⛔ {bb}: {len(seqs)} designs, expected 100")
        L = {len(s) for s, _ in seqs.values()}
        if len(L) != 1:
            sys.exit(f"⛔ {bb}: binder lengths differ: {sorted(L)}")
        blen = L.pop()

        # validate the extraction CODE PATH on the ids that actually folded
        for i in sorted(HELD_IDS):
            name = f"BR_{bb}_id{i}"
            if name in ref:
                checked += 1
                if ref[name] != seqs[i][0]:
                    mismatch.append(name)

        cand = [(i, s, c) for i, (s, c) in seqs.items() if i not in HELD_IDS]
        if len(cand) != 100 - len(HELD_IDS):
            sys.exit(f"⛔ {bb}: {len(cand)} candidates, expected {100-len(HELD_IDS)}")
        cand.sort(key=lambda t: (-t[2], t[0]))        # confidence desc, id asc as tiebreak
        pick = cand[:n_per]
        if len(pick) != n_per:
            sys.exit(f"⛔ {bb}: could only pick {len(pick)} of {n_per}")
        for i, s, c in pick:
            rows.append((f"BR_{bb}_id{i}", s))
        print(f"{bb:40s} {blen:>7} {pick[0][2]:.4f}-{pick[-1][2]:.4f} "
              f"{len(HELD_IDS):>12}")

    print(f"\n⭐ extraction cross-check vs sequences that actually folded: "
          f"{checked} compared, {len(mismatch)} mismatched")
    if mismatch:
        sys.exit(f"⛔⛔ EXTRACTION WRONG -- {mismatch[:3]}; refusing to write")
    if checked < 80:
        sys.exit(f"⛔ only {checked} cross-checked (expected 80); refusing")

    names = [n for n, _ in rows]
    if len(set(names)) != len(names):
        dup = [k for k, v in collections.Counter(names).items() if v > 1]
        sys.exit(f"⛔ duplicate names: {dup[:3]}")
    if set(names) & set(ref):
        overlap = sorted(set(names) & set(ref))
        sys.exit(f"⛔⛔ {len(overlap)} picked sequences were ALREADY FOLDED "
                 f"({overlap[:3]}) -- that is a duplicate campaign; refusing")
    with open(outp, "w") as fh:
        for n, s in rows:
            fh.write(f"{n}\t{s}\n")
    print(f"wrote {outp}: {len(rows)} rows ({len(BACKBONES)} backbones x {n_per}), "
          f"{len({s for _, s in rows})} distinct sequences")
    print(f"⛔ NONE of these overlap the {len(ref)} already-folded sequences -- asserted above.")
    bl = sorted({len(s) for _, s in rows})
    print(f"binder lengths: {bl}  <- VRAM ceiling set by {max(bl)} aa")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3],
         int(sys.argv[4]) if len(sys.argv) > 4 else 50)
