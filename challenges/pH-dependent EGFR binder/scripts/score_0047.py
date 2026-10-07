#!/usr/bin/env python3
"""The `decisions/0047` site-A objective, as a scorer. NEW FILE 2026-10-05.

⭐ WHY THIS EXISTS. `decisions/0047` was accepted 2026-10-05 and NOTHING implemented it:
`score_cycle_egfrA1.py` still ranks on `full_epitope_frac` ("fef") and on `label` (mean
n_epitope over gate-passing draws). Every landed number in every `training_c*.csv` predates
the ADR. ⛔ So label the metric whenever you quote a value from either side.

⭐⭐ THE OBJECTIVE, copied from `decisions/0047` and not reinterpreted:

    c380           = 1 if the draw contacts F380, else 0
    per-draw score = min(n_epitope / 11, c380)   if the draw FOLDS (sc_rmsd < 1.5)
                   = 0                            otherwise
    label_A_0047   = mean of that over EVERY draw Boltz produced for the sequence

⛔⛔ F380 IS INSIDE THE MIN. That is the whole point -- it enforces the same-draw
   requirement. A draw must hold the footprint AND
   F380 simultaneously to score anything. ⛔ NEVER a sum, NEVER a per-residue or per-region
   average first: measured on real draws, `id2039` holds every one of the 11 residues in >=50%
   of its draws yet covers all 11 in 1 draw in 10, and on the bridge `seq 44` scores 0.857 on a
   per-site average with ZERO both-site draws. Averaging marginals rewards a design that never
   once does the thing we want.
⛔ THE MEAN IS OVER **EVERY** DRAW, not over gate-passing draws. A misfolded draw contributes a
   real 0, because scoring 0 is different from not reporting a score: a misfolded or
   mispositioned protein is still a data point. This is the one place where this scorer
   deliberately differs in SHAPE from `label` in `score_cycle_egfrA1.py`, which averages over
   gate-passers only. Do not "fix" it into agreement.

⛔⛔ IT REFUSES RATHER THAN DEGRADES. `c380` must be present in the contacts CSV. A contacts
   file written before 2026-10-05 does NOT have it (`contacts_a21.py` computed the hit set and
   recorded only its size), and the only honest response is to re-run the patched
   `contacts_a21.py` over the structures -- which is free and local. ⛔ There is NO fallback,
   NO inference of c380 from n_epitope, and no silent default: that failure mode silently
   publishes a stale table at rc=0, which is the worst possible outcome.

usage:
  score_0047.py --pair <contacts.csv> <sc.csv> [--pair <contacts.csv> <sc.csv> ...]
                [--tsv <cycle.tsv> ...]          # to carry sequences through
                [--out <per_sequence.csv>]
                [--augment <training.csv> --augment-out <out.csv>]

  Each --pair is tagged by the name of the directory holding its contacts CSV, so pooling
  several cycles keeps them distinguishable and a seq id reused across cycles cannot collide.
"""
import argparse
import collections
import csv
import os
import statistics
import sys

GATE = 1.5          # sc_rmsd < GATE == "the draw folded". Same constant as score_cycle_egfrA1.py.
FULL = 11           # the site-A footprint has 11 residues
REQUIRED = "c380"   # the column this scorer refuses to run without


def load_pair(contacts, scf):
    """Per-draw records for one cycle. Refuses on a missing c380 or a key mismatch."""
    with open(contacts) as fh:
        rd = csv.DictReader(fh)
        if REQUIRED not in (rd.fieldnames or []):
            sys.exit(
                f"⛔ {contacts} has no `{REQUIRED}` column -- columns are "
                f"{rd.fieldnames}.\n"
                f"   This file was written by the PRE-0047 `contacts_a21.py`, which computed the\n"
                f"   contacted-residue set and recorded only its size. decisions/0047 needs the\n"
                f"   F380 bit per draw and it CANNOT be recovered from `n_epitope`.\n"
                f"   ⇒ Re-run the patched scorer over the structures (local, read-only):\n"
                f"       <a python that carries biotite> contacts_a21.py <land_dir> <out.csv>\n"
                f"   ⛔ Refusing rather than guessing.")
        con = {(r["seq"], r["model"]): r for r in rd}
    with open(scf) as fh:
        sc = {(r["seq"], r["model"]): r for r in csv.DictReader(fh)}
    if set(con) != set(sc):
        only_c, only_s = sorted(set(con) - set(sc))[:5], sorted(set(sc) - set(con))[:5]
        sys.exit(f"⛔ contacts/sc key mismatch for {contacts} + {scf}: "
                 f"{len(set(con)-set(sc))} only in contacts (e.g. {only_c}), "
                 f"{len(set(sc)-set(con))} only in sc (e.g. {only_s}). Refusing.")
    out = collections.defaultdict(list)
    for k, c in con.items():
        rmsd = float(sc[k]["sc_rmsd"])
        out[k[0]].append(dict(nep=int(c["n_epitope"]),
                              c380=int(c[REQUIRED]),
                              sc=rmsd,
                              folds=rmsd < GATE))
    return out


def per_draw_0047(d):
    """decisions/0047, one draw. 0 on a misfold; F380 INSIDE the min."""
    if not d["folds"]:
        return 0.0
    return min(d["nep"] / FULL, float(d["c380"]))


def fef(draws):
    """The CURRENT metric, for the side-by-side only: fraction of ALL draws with n_epitope == 11."""
    return sum(1 for d in draws if d["nep"] == FULL) / len(draws)


def legacy_label(draws):
    """`label` as score_cycle_egfrA1.py defines it: mean n_epitope over GATE-PASSING draws."""
    g = [d["nep"] for d in draws if d["folds"]]
    return statistics.mean(g) if g else float("nan")


def kendall_tau(xs, ys):
    """Tau-b. No scipy in this env (same reason score_cycle_egfrA1.py hand-rolls its stats)."""
    n = len(xs)
    con = dis = tx = ty = 0
    for i in range(n):
        for j in range(i + 1, n):
            a, b = xs[i] - xs[j], ys[i] - ys[j]
            if a == 0 and b == 0:
                tx += 1; ty += 1
            elif a == 0:
                tx += 1
            elif b == 0:
                ty += 1
            elif (a > 0) == (b > 0):
                con += 1
            else:
                dis += 1
    d0 = ((con + dis + tx) * (con + dis + ty)) ** 0.5
    return (con - dis) / d0 if d0 else float("nan")


def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("--pair", nargs=2, action="append", metavar=("CONTACTS", "SC"),
                    required=True)
    ap.add_argument("--tsv", action="append", default=[],
                    help="cycle tsv (fold_id<TAB>sequence) to carry sequences through")
    ap.add_argument("--out", help="per-sequence CSV to write")
    ap.add_argument("--alias-map", action="append", default=[], metavar="MAP.CSV",
                    help="a replicates map (replicate_id,source_id,...). Draws folded under a "
                         "replicate id are POOLED ONTO THEIR SOURCE DESIGN, so a design measured "
                         "twice yields one row over all its draws. \u26d4 Without this a re-fold "
                         "appears as a separate design and the top of the table is the top of a "
                         "SINGLE optimistic measurement -- the residual flagged under "
                         "decisions/0045.")
    ap.add_argument("--augment", help="an existing training table to add label_A_0047 to")
    ap.add_argument("--augment-out", help="where to write the augmented table")
    a = ap.parse_args()
    if bool(a.augment) != bool(a.augment_out):
        sys.exit("⛔ --augment and --augment-out go together")

    seqs = {}
    for t in a.tsv:
        for line in open(t):
            if "\t" in line:
                k, v = line.rstrip("\n").split("\t", 1)
                if "_id" in k:
                    seqs[k.rsplit("_id", 1)[1]] = v

    # ⛔ replicate_id -> source_id, so a re-fold pools onto the design it re-measures.
    alias = {}
    for mf in a.alias_map:
        for r in csv.DictReader(open(mf)):
            if "replicate_id" not in r or "source_id" not in r:
                sys.exit(f"\u26d4 {mf} is not a replicates map (needs replicate_id,source_id) "
                         f"-- columns are {list(r)}. Refusing.")
            alias[r["replicate_id"]] = r["source_id"]
    if alias:
        print(f"  alias map: {len(alias)} re-folds pooled onto their sources "
              f"({', '.join(f'{k}->{v}' for k, v in sorted(alias.items()))})")

    rows, per_tag = [], {}
    for contacts, scf in a.pair:
        tag = os.path.basename(os.path.dirname(os.path.abspath(contacts))) or "pair"
        got = load_pair(contacts, scf)
        per_tag[tag] = got
        ndraw = sum(len(v) for v in got.values())
        print(f"  {tag:14s} {len(got):4d} sequences  {ndraw:5d} draws   ({contacts})")

    print(f"\n=== decisions/0047 site-A objective: mean over EVERY draw of "
          f"min(n_epitope/11, c380), 0 on a misfold ===")
    # pool by DESIGN. Without an alias map this is one row per fold id, exactly as before.
    pool, tags, poolseq = collections.defaultdict(list), collections.defaultdict(set), {}
    for tag, got in per_tag.items():
        for sid, draws in got.items():
            design = alias.get(f"id{sid}", f"id{sid}")
            pool[design] += draws
            tags[design].add(tag)
            # ⛔ a re-fold carries the SAME sequence as its source; keep whichever we have.
            if seqs.get(sid):
                poolseq.setdefault(design, seqs[sid])
    for design, draws in sorted(pool.items()):
        sid = design[2:]
        tag = ",".join(sorted(tags[design]))
        d47 = [per_draw_0047(d) for d in draws]
        rows.append(dict(
            tag=tag, seq_id=design, sequence=poolseq.get(design, seqs.get(sid, "")),
            n_draws=len(draws),
            n_folded=sum(1 for d in draws if d["folds"]),
            n_c380=sum(d["c380"] for d in draws),
            c380_frac=round(sum(d["c380"] for d in draws) / len(draws), 3),
            # ⭐ THE 0047 METRIC. Name it in full in every table so it can never be
            #    confused with `label` or `fef`, which rank differently.
            label_A_0047=round(statistics.mean(d47), 4),
            n_draws_scoring_zero=sum(1 for v in d47 if v == 0.0),
            # the two legacy metrics, for the side-by-side ONLY
            fef_legacy=round(fef(draws), 3),
            label_legacy=(round(legacy_label(draws), 4)
                          if legacy_label(draws) == legacy_label(draws) else ""),
            median_sc=round(statistics.median(d["sc"] for d in draws), 3),
            n_measurements=len(tags[design])))

    rows.sort(key=lambda r: -r["label_A_0047"])

    # ---- resolution: 0047's own table, recomputed on whatever was handed in ----
    n = len(rows)
    v47 = [r["label_A_0047"] for r in rows]
    vfef = [r["fef_legacy"] for r in rows]
    print(f"\n  pooled sequences: {n}   total draws: {sum(r['n_draws'] for r in rows)}")
    print(f"  {'metric':24s} {'distinct':>9s} {'zeros':>7s} {'median':>8s} {'max':>7s}")
    for nm, v in (("label_A_0047 (0047)", v47), ("fef (the CURRENT metric)", vfef)):
        print(f"  {nm:24s} {len(set(v)):9d} {sum(1 for x in v if x == 0):7d} "
              f"{statistics.median(v):8.3f} {max(v):7.3f}")
    print(f"  Kendall tau(label_A_0047, fef) = {kendall_tau(v47, vfef):+.3f}"
          "   ⇒ below 1.0 means the two metrics RANK DIFFERENTLY; re-scoring renumbers the table.")

    print(f"\n  F380 contact rate over all pooled draws: "
          f"{sum(r['n_c380'] for r in rows)}/{sum(r['n_draws'] for r in rows)} = "
          f"{100*sum(r['n_c380'] for r in rows)/sum(r['n_draws'] for r in rows):.1f}%")

    print(f"\n=== TOP 12 under decisions/0047 (⛔ label_A_0047, NOT fef) ===")
    print(f"  {'rank':>4s} {'seq_id':9s} {'tag':26s} {'0047':>7s} {'fef':>6s} "
          f"{'c380':>6s} {'fold':>7s} {'med_sc':>7s}  fef_rank")
    fef_order = {r["seq_id"]: i + 1 for i, r in
                 enumerate(sorted(rows, key=lambda r: -r["fef_legacy"]))}
    for i, r in enumerate(rows[:12], 1):
        print(f"  {i:4d} {r['seq_id']:9s} {r['tag']:26s} {r['label_A_0047']:7.3f} "
              f"{r['fef_legacy']:6.2f} {r['c380_frac']:6.2f} "
              f"{r['n_folded']:3d}/{r['n_draws']:<3d} {r['median_sc']:7.3f}  "
              f"{fef_order[r['seq_id']]:>4d}")

    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\nwrote {a.out}: {len(rows)} sequences")

    if a.augment:
        # ⛔ keyed on seq_id. A sequence with no fresh draws gets an EMPTY cell, never a 0 --
        #    "not measured" and "measured as zero" are different facts (0046/0047).
        by_id = {}
        for r in rows:
            by_id.setdefault(r["seq_id"], r)
        base = list(csv.DictReader(open(a.augment)))
        flds = list(base[0].keys())
        for col in ("label_A_0047", "c380_frac"):
            if col not in flds:
                flds.append(col)
        hit = 0
        for b in base:
            r = by_id.get(b["seq_id"])
            if r:
                b["label_A_0047"] = r["label_A_0047"]; b["c380_frac"] = r["c380_frac"]; hit += 1
            else:
                b.setdefault("label_A_0047", ""); b.setdefault("c380_frac", "")
        with open(a.augment_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=flds)
            w.writeheader(); w.writerows(base)
        print(f"wrote {a.augment_out}: {len(base)} rows, {hit} carry label_A_0047, "
              f"{len(base)-hit} left EMPTY (no per-draw data handed in -- ⛔ not zero)")


if __name__ == "__main__":
    main()
