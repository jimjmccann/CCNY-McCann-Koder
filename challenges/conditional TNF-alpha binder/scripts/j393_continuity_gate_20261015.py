#!/usr/bin/env python3
"""j393_continuity_gate_20261015.py -- is the designed chain ONE CONNECTED CHAIN?

WHY THIS EXISTS: `J392` §1 measured that **265/384 (69.0%)** of TNFa v36 carry a CA-CA step over
10 A, median gap 56 A, and that the `bb_min >= 2.0` clash floor **REWARDS** the break -- 4 of the
top 5 by `bb_min` are disconnected, because a chain with half of itself 56 A away has fewer close
contacts. So `J390`'s "20/384 clear the floor" was right as a CLASH statement and wrong as a
SURVIVOR count (real figure 9/384). ⛔⛔ Continuity is therefore gated **BEFORE** any clash or
pin metric, never after.

`J219` had already measured this once (zero steps > 4.2 A: 9/24) and it was not carried forward.
`j392_continuity_20261015.py` measured it again as a ONE-OFF with the chain hardcoded to "A", the
path hardcoded to the v36 extract, and its output written to a temporary directory that is now
gone. This is that measurement made REUSABLE, so the GLP1R arms and every later arm get it for
free.

⛔ A SCORER MUST NEVER NAME A CHAIN (fourth incident -- a hardcoded
chain pair compared a designed chain against the TARGET and reported it as a self-interface, with
no error). Roles here are DETECTED by `j261`'s discriminator -- the fraction of residues carrying
`is_motif_atom_with_fixed_seq`, which is ~100% for a target (fixed by definition) and a few % for
a designed binder (just its pinned guideposts). The chains actually used are PRINTED and written
into the TSV, so a wrong assignment is visible in the output rather than silent.

⛔ ONE ROW PER (design, binder chain) -- never collapsed. A C3 arm has three binder copies and
partial data must not be ranked: a trimer whose third copy is broken is not 2/3 good.

⛔ THE PRIMARY COLUMN IS CONTINUOUS (`max_step`), not the pass/fail flag
(a binary gate at small n returns all-zeros and ranks
nothing). `continuous` is DERIVED for counting, never for ranking.

⛔ NUMBERING GAPS ARE NOT BREAKS. Steps are measured only between residues whose numbers are
CONSECUTIVE; a jump in residue numbering is reported in its own column (`n_numbering_gaps`)
because a missing residue is a different defect from a chain flying apart, and the one-off script
conflated them by walking sorted order regardless.

⛔ Local, no spend, no model, nothing launched.

usage:
  j393_continuity_gate_20261015.py --out <fresh.tsv> '<glob>' ['<glob>' ...]
    [--step-warn 4.2] [--step-break 10.0] [--force]
"""

import argparse
import glob as globmod
import gzip
import math
import os
import statistics
import sys
from collections import defaultdict

# The column holds the STRINGS "True"/"False", not 1/0. ⛔ ALLOWLIST, never a denylist of falsy
# spellings -- a denylist read "False" as TRUE and called all six chains TARGET at 100% (j261).
TRUTHY = ("True", "true", "TRUE", "1", "Y", "y", "yes")

# Ideal trans CA-CA is 3.80-3.81 A. These are the J219/J392 bars, named here so they are
# overridable and RECORDED in the output rather than inherited silently
# (never inherit a default: record the resolved config).
DEFAULT_STEP_WARN = 4.2
DEFAULT_STEP_BREAK = 10.0


def read_cif_atoms(path):
    """Minimal CIF atom_site reader -- same parser as j261, no biotite dependency."""
    opener = gzip.open if path.endswith(".gz") else open
    cols, rows, in_loop = [], [], False
    with opener(path, "rt") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1])
                in_loop = True
                continue
            if in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
                f = s.split()
                if len(f) == len(cols):
                    rows.append(dict(zip(cols, f)))
            elif in_loop and rows and (s.startswith("#") or s.startswith("loop_")):
                break
    return cols, rows


def classify_chains(cols, rows, ch_key, seq_key, fix_key):
    """j261's discriminator: fixed-seq fraction per chain. Returns (binder, target, ambiguous, stats)."""
    res, fixed = defaultdict(set), defaultdict(set)
    for r in rows:
        c = r[ch_key]
        res[c].add(r[seq_key])
        if fix_key and r.get(fix_key, "") in TRUTHY:
            fixed[c].add(r[seq_key])
    binder, target, ambiguous, stats = [], [], [], {}
    for c in sorted(res):
        n, nf = len(res[c]), len(fixed.get(c, ()))
        pct = 100.0 * nf / n if n else 0.0
        stats[c] = (n, nf, pct)
        if pct > 90:
            target.append(c)
        elif pct < 20:
            binder.append(c)
        else:
            ambiguous.append(c)
    return binder, target, ambiguous, stats


def chain_continuity(rows, chain, ch_key, seq_key, fix_key, step_warn, step_break):
    """CA ladder for ONE chain. Steps only between CONSECUTIVELY NUMBERED residues."""
    ca, pins = {}, set()
    for r in rows:
        if r[ch_key] != chain:
            continue
        try:
            n = int(r[seq_key])
        except (ValueError, TypeError):
            continue  # label_seq_id is "." for HETATM
        if fix_key and r.get(fix_key, "") in TRUTHY:
            pins.add(n)
        if r["label_atom_id"] == "CA":
            ca[n] = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))

    ks = sorted(ca)
    if len(ks) < 2:
        return None
    pins = sorted(pins)

    steps, gaps = [], 0
    for a, b in zip(ks, ks[1:]):
        if b - a == 1:
            steps.append((a, math.dist(ca[a], ca[b])))
        else:
            gaps += 1
    if not steps:
        return None

    at_resi, max_step = max(steps, key=lambda t: t[1])
    n_warn = sum(1 for _, d in steps if d > step_warn)
    n_break = sum(1 for _, d in steps if d > step_break)

    # Where did the worst step fall relative to the PINNED span? J392: the breaks are ONLY in the
    # unanchored tails (N-term 132 / C-term 133 of 265), just 4/384 between pins -- so the pinned
    # scaffold is intact and the fix is a CONTIG change, not relaxation.
    if not pins:
        zone = "no_pins"
    elif at_resi < pins[0]:
        zone = "N_TERM_FREE"
    elif at_resi >= pins[-1]:
        zone = "C_TERM_FREE"
    else:
        zone = "BETWEEN_PINS"

    return {
        "n_res": len(ks),
        "first_res": ks[0],
        "last_res": ks[-1],
        "n_pins": len(pins),
        "first_pin": pins[0] if pins else 0,
        "last_pin": pins[-1] if pins else 0,
        "n_free_nterm": (pins[0] - ks[0]) if pins else 0,
        "n_free_cterm": (ks[-1] - pins[-1]) if pins else 0,
        "max_step": max_step,
        "at_resi": at_resi,
        "median_step": statistics.median(d for _, d in steps),
        "n_gt_warn": n_warn,
        "n_gt_break": n_break,
        "n_numbering_gaps": gaps,
        "break_zone": zone,
        "continuous": int(n_warn == 0 and gaps == 0),
    }


COLUMNS = [
    "design", "binder_chain", "chains_binder", "chains_target", "n_res", "first_res", "last_res",
    "n_pins", "first_pin", "last_pin", "n_free_nterm", "n_free_cterm", "max_step", "at_resi",
    "median_step", "n_gt_warn", "n_gt_break", "n_numbering_gaps", "break_zone", "continuous",
]


def main(argv):
    ap = argparse.ArgumentParser(description="backbone continuity gate for RFd3 designed chains")
    ap.add_argument("--out", required=True, help="FRESH tsv path (refuses to overwrite)")
    ap.add_argument("--step-warn", type=float, default=DEFAULT_STEP_WARN)
    ap.add_argument("--step-break", type=float, default=DEFAULT_STEP_BREAK)
    ap.add_argument("--force", action="store_true", help="allow overwriting --out")
    ap.add_argument("globs", nargs="+")
    a = ap.parse_args(argv)

    # ⛔ Scripts APPEND to --out, so pointing a run at an existing results file silently
    #    merges two campaigns. Fresh --out per experiment; refuse rather than merge.
    if os.path.exists(a.out) and not a.force:
        print(f"⛔ --out already exists: {a.out}\n   Use a FRESH --out per experiment, "
              f"or pass --force deliberately.")
        return 2

    files = []
    for g in a.globs:
        files.extend(sorted(globmod.glob(g)))
    if not files:
        print(f"⛔ no files matched {a.globs} -- ⛔ an empty glob is NOT a result. "
              f"Re-check the path WITHOUT 2>/dev/null and with an ABSOLUTE path.")
        return 2

    print(f"=== j393 continuity gate ===")
    print(f"  files matched      : {len(files)}")
    print(f"  step-warn (break?) : {a.step_warn} A      step-break (catastrophic): {a.step_break} A")
    print(f"  out                : {a.out}")

    rows_out, skipped, layouts = [], [], defaultdict(int)
    for path in files:
        cols, atoms = read_cif_atoms(path)
        if not atoms:
            skipped.append((path, f"no atom_site rows (columns seen: {len(cols)})"))
            continue
        ch_key = "label_asym_id" if "label_asym_id" in cols else "auth_asym_id"
        seq_key = "label_seq_id" if "label_seq_id" in cols else "auth_seq_id"
        fix_key = next((c for c in cols if "fixed_seq" in c), None)
        if not fix_key:
            skipped.append((path, "⛔ no *fixed_seq* column -- roles UNDETECTABLE, and J248 §4 says "
                                  "length alone is not a sound discriminator. NOT scored."))
            continue

        models = {r.get("pdbx_PDB_model_num", "1") for r in atoms}
        if len(models) > 1:
            skipped.append((path, f"⛔ {len(models)} models in one file -- ambiguous. NOT scored."))
            continue

        binder, target, ambiguous, stats = classify_chains(cols, atoms, ch_key, seq_key, fix_key)
        if ambiguous:
            skipped.append((path, f"⛔⛔ AMBIGUOUS chains {ambiguous} "
                                  f"(fixed% {[round(stats[c][2], 1) for c in ambiguous]}) -- NOT scored"))
            continue
        if not binder:
            skipped.append((path, f"⛔ no binder chain detected (all target: {target}) -- NOT scored"))
            continue

        layouts[f"binder={','.join(binder)} target={','.join(target) or 'NONE'}"] += 1
        for c in binder:
            m = chain_continuity(atoms, c, ch_key, seq_key, fix_key, a.step_warn, a.step_break)
            if m is None:
                skipped.append((path, f"chain {c}: fewer than 2 consecutively-numbered CA atoms"))
                continue
            m["design"] = os.path.basename(path)
            m["binder_chain"] = c
            m["chains_binder"] = ",".join(binder)
            m["chains_target"] = ",".join(target) or "NONE"
            rows_out.append(m)

    with open(a.out, "w") as fo:
        # resolved config recorded IN THE ARTEFACT, never inherited
        fo.write(f"# j393_continuity_gate_20261015.py  step_warn={a.step_warn} "
                 f"step_break={a.step_break} n_files={len(files)} n_rows={len(rows_out)}\n")
        fo.write("\t".join(COLUMNS) + "\n")
        for m in rows_out:
            # ⛔⛔ THREE DECIMALS, NOT TWO, AND IT IS NOT COSMETIC. `J392` reported 93/384
            #    continuous where its own script's flag said 91: two designs have a true max_step
            #    in (4.200, 4.205] and PRINTED as "4.20" at 2 dp, so a filter on the printed
            #    column let them through. ⛔ Filter on the `n_gt_warn` FLAG, never on a rounded
            #    distance column: a count is not a measurement.
            fo.write("\t".join(
                f"{m[k]:.3f}" if isinstance(m[k], float) else str(m[k]) for k in COLUMNS) + "\n")

    # ---- summary. ⛔ Report the SERIES, not just the headline.
    print(f"\n  LAYOUTS OBSERVED (⛔ detected, never assumed):")
    for k, v in sorted(layouts.items(), key=lambda t: -t[1]):
        print(f"    {v:>5}  {k}")

    if not rows_out:
        print("\n⛔ ZERO rows scored -- that is a FAILURE, not a result of 0% continuity.")
        for p, why in skipped[:10]:
            print(f"    {os.path.basename(p)}: {why}")
        return 1

    n = len(rows_out)
    cont = sum(m["continuous"] for m in rows_out)
    brk = sum(1 for m in rows_out if m["n_gt_break"] > 0)
    steps = sorted(m["max_step"] for m in rows_out)
    print(f"\n  rows (design x binder chain) : {n}")
    print(f"  CONTINUOUS (no step > {a.step_warn} A, no numbering gap) : {cont}/{n} "
          f"({100.0 * cont / n:.1f}%)")
    print(f"  any step > {a.step_break} A (catastrophic)              : {brk}/{n} "
          f"({100.0 * brk / n:.1f}%)")
    print(f"  max_step  min / median / max : {steps[0]:.3f} / {statistics.median(steps):.3f} "
          f"/ {steps[-1]:.3f} A")

    # ⛔⛔ THE BAR IS FRAGILE AND THE LADDER MAKES THAT VISIBLE. Measured on v36: 9 designs clear
    #    4.0 A and 116 clear 4.3 A -- a 13x swing over 0.3 A, because the bar cuts straight
    #    through the bulk of the distribution. ⇒ never quote a continuity COUNT without the bar,
    #    and rank on `max_step`.
    print(f"  ⚠️ BAR SENSITIVITY -- rows with NO step above each bar (the count IS the bar):")
    for bar in (4.0, 4.1, a.step_warn, 4.3, 4.5, 5.0, a.step_break):
        k = sum(1 for m in rows_out if m["max_step"] <= bar)
        tag = "  <-- --step-warn" if bar == a.step_warn else (
            "  <-- --step-break" if bar == a.step_break else "")
        print(f"    <= {bar:>5.1f} A : {k:>5}/{n} ({100.0 * k / n:>5.1f}%){tag}")

    # ⛔ TWO DENOMINATORS, BOTH NAMED. J392 quoted 132 N-term / 133 C-term, which is the
    #    >step_break set; the >step_warn set is a different (larger) cut and the two must never be
    #    swapped.
    for label, pred in ((f"a step > {a.step_break} A (CATASTROPHIC)", lambda m: m["n_gt_break"] > 0),
                        (f"a step > {a.step_warn} A (warn bar)", lambda m: m["n_gt_warn"] > 0)):
        sel = [m for m in rows_out if pred(m)]
        z = defaultdict(int)
        for m in sel:
            z[m["break_zone"]] += 1
        print(f"  WHERE the worst step falls, among the {len(sel)} with {label}:")
        for zz in ("N_TERM_FREE", "BETWEEN_PINS", "C_TERM_FREE", "no_pins"):
            if z.get(zz):
                print(f"    {z[zz]:>5}  {zz}")
    frees = [m["n_free_nterm"] for m in rows_out]
    frec = [m["n_free_cterm"] for m in rows_out]
    print(f"  UNANCHORED residues before first pin : {min(frees)}-{max(frees)} "
          f"(median {statistics.median(frees):.0f})")
    print(f"  UNANCHORED residues after  last pin  : {min(frec)}-{max(frec)} "
          f"(median {statistics.median(frec):.0f})")
    if skipped:
        print(f"\n  ⚠️ {len(skipped)} file(s) NOT scored -- ⛔ absent from the TSV is not a pass:")
        for p, why in skipped[:10]:
            print(f"    {os.path.basename(p)}: {why}")
        if len(skipped) > 10:
            print(f"    ... and {len(skipped) - 10} more")
    print(f"\n  ⇒ wrote {a.out}")
    print(f"  ⛔ RANK ON max_step (continuous); 'continuous' is for COUNTING only.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
