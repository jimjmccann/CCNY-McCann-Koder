#!/usr/bin/env python3
"""s75 -- score a SCAFFOLDED RFd3 arm. Reuses `s51`'s GATES VERBATIM; replaces ONLY role resolution.

⛔⛔ WHY THIS EXISTS AND WHY `s51` IS NOT EDITED (`J1001`).
   `s51` refused **320 of 320** `hit02` designs with
       `expected exactly 1 binder chain, got [] (present=[A,B,C,D], target=[A,B,C,D])`
   and ⭐ the refusal was CORRECT. `s51` resolves the binder as *the output chain with NO named
   contig tokens*. `[M]` In a SCAFFOLDED arm the binder element is a MIX:
       26-36,E173-181,8-14,E198-213,26-35,/0,A1-276,/0,B1-99,/0,C1-10
   the generated bands PLUS named scaffold segments from chain `E` (AF-Q5F8R2 = `hit02`, carrying
   the four amide guideposts). `[M]` `E173 -> A35`, so chain `E` maps into output `A` -- the BINDER
   -- and `s51` reads every *named* segment as target ⇒ all four output chains classify as target
   and the binder set is empty.
   ⛔ `s51` is left UNMODIFIED because `J997`/`ADR 0063` were adjudicated on it; and ⛔ its own
   docstring already refused this widening for its predecessor `j948`: *"a rule loose enough to
   admit this construct would admit a genuine target too."* The same refusal applies to `s51`.

⭐⭐ THE CONSTRUCTION THAT MAKES THIS SAFE: this file does NOT re-implement a single gate.
   It imports `s51` as a module and rebinds ONE name, `s51.roles`. `score_one` and every helper
   (`read_cif`, `ca_of`, `heavy_of`, `mindist`, `contacts`, `res_atoms`, `seq_of`) are the SAME
   OBJECTS from the same file on disk. ⇒ gate drift between the matrix arm and this arm is
   IMPOSSIBLE BY CONSTRUCTION, not merely tested for. A copy-paste scorer would have made the two
   arms incomparable, which is the whole reason the matrix numbers exist.
   ⚠️ `s51` runs its scoring loop at import, so it is imported with `sys.argv` stubbed to an empty
   glob list and its `OUT` pointed at a throwaway temp file: it scores 0, writes a header, returns.
   ⛔ That temp file is NOT a result and is deleted.

THE NEW ROLE RULE -- from the contig's own `/0` SEGMENTATION, never a chain letter, never "is it
named" (`J1001` §5):
  element 0            = the BINDER   (generated bands and/or named scaffold segments)
  elements 1..n        = the TARGET chains (`A1-276`, `B1-99`, `C1-10`)
  each element's named tokens map through `diffused_index_map` to exactly one OUTPUT chain.
  ⭐ STRICT SUPERSET OF `s51`: if element 0 has NO named tokens (the `arm3_hotspot_only` regime)
     this falls back to `s51`'s own rule -- the structure chain that is not a target chain -- so
     the matrix arm resolves IDENTICALLY.
  ⛔⛔ CORRECTED 2026-11-04 -- THIS DOCSTRING HAD THE LOGIC INVERTED. It used to end *"That is
     what makes `--compare` a real check."* The superset property is what makes `--compare` a
     **TAUTOLOGY** on the matrix arm, NOT what makes it a real check. `[M]` All 384 matrix
     reference rows have ZERO named tokens in contig element 0 (only 2 distinct matrix contigs
     exist, `65-85,/0,...` and `85-110,/0,...`), so `binder_named` is empty on every one of them
     and `roles_scaffolded` takes its `else:` fallback -- which is `s51`'s rule verbatim.
     ⇒ 0 diffs is ENTAILED BY THE CONSTRUCTION, and the new SCAFFOLDED branch -- the branch the
     320 `hit02` numbers actually rest on -- is exercised on **0 of 384** reference rows.
     ⭐ WHAT `--compare` DOES PROVE: a REGRESSION test that `s75` does not break `s51`'s own
        path (and it proves it byte-identically, which is stronger than it was claimed).
     ⛔ WHAT IT DOES NOT PROVE: the scaffolded role resolution itself. `J1002`'s *"proven before
        it was used"* was overstated and is bannered there.
     ⭐ THE GENUINE POSITIVE CONTROL FOR THE SCAFFOLDED BRANCH: resolved `binder_len` spans
        **68-104 over 18 distinct values** across the two designed length bands, where any
        mis-resolved TARGET chain would be a CONSTANT (HLA 276, b2m 99, peptide 10). `[M]` on
        `t90`: **0 designs at 276, 0 at 10.** A mis-resolved binder cannot produce the designed
        length distribution.
     ⭐ AND THE CONTROL IS STRONGER THAN `J1002`'s BANNER CLAIMED -- measured 2026-11-04, the
        two bands SEPARATE and do not overlap: `L65_85` gives **68-80** (8 distinct, n=160) and
        `L85_110` gives **92-104** (10 distinct, n=160). b2m's 99 is therefore **UNREACHABLE in
        the `L65_85` band**, so on those 160 designs the b2m ambiguity is excluded OUTRIGHT, not
        merely argued away. ⇒ `J1002`'s *"not airtight against b2m alone"* caveat is correct but
        understated: it applies only to the `L85_110` half, where `[M]` 24 designs legitimately
        resolve to 99 -- and b2m would be 99 for ALL 320, not 24 of 160 in one band.
  ⚠️ `[M]` the `/0` here are COMMA-delimited (`,/0,`) so elements split cleanly; a `/0` without a
     comma is swallowed and silently changes the output layout
     ([[rfd3-contig-chain-break-is-swallowed-without-a-comma]]) ⇒ asserted below, per design.

GATES, ORDER AND UNITS ARE `s51`'S AND UNCHANGED: continuity FIRST
([[gate-backbone-continuity-before-any-clash-metric]]), raw floor `bb_min >= 2.0 A`
([[raw-rfd3-output-clashes-normally]]), contacts STRATIFIED peptide-vs-HLA
([[footprint-engagement-not-epitope-validation]]), rank on the CONJUNCTION
([[rank-on-conjunction-not-mean-of-fractions]]), and ⛔⛔ **n is BATCHES, not designs**
([[rfd3-designs-are-clustered-by-batch]]) with key `s<shard>b<batch>`.

Usage:
  python3 -I s75_score_scaffolded_arm.py <out.tsv> <glob> [<glob> ...]
  python3 -I s75_score_scaffolded_arm.py <out.tsv> <glob> --compare <reference.tsv>
      ⭐ `--compare` is the PRE-REGISTERED EQUIVALENCE GATE (`J1001` §6): score the MATRIX tree and
        require EXACT per-row, per-column agreement with `s51`'s `t80`. Any disagreement ⇒ this
        instrument is wrong and the scaffolded numbers DO NOT SHIP.
"""
import importlib.util
import json
import os
import re
import sys
import tempfile
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S51 = os.path.join(HERE, "s51_score_matrix.py")

argv = sys.argv[1:]
if "--compare" in argv:
    i = argv.index("--compare")
    REF = argv[i + 1]
    argv = argv[:i] + argv[i + 2:]
else:
    REF = None
if len(argv) < 2:
    sys.exit(__doc__)
OUT, GLOBS = argv[0], argv[1:]


def load_s51():
    """Import s51 with argv stubbed so its module-level loop scores nothing.

    ⛔ Do NOT 'fix' this by refactoring s51 -- it must stay byte-identical on disk."""
    if not os.path.exists(S51):
        raise SystemExit(f"⛔ s51 not found at {S51}")
    fd, tmp = tempfile.mkstemp(prefix="s75_s51_import_", suffix=".tsv")
    os.close(fd)
    saved = sys.argv
    sys.argv = ["s51_score_matrix.py", tmp]        # no globs => 0 files => 0 rows
    try:
        spec = importlib.util.spec_from_file_location("s51_score_matrix", S51)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
        try:
            os.unlink(tmp)
        except OSError:
            pass
    for need in ("roles", "score_one", "read_cif", "seq_of", "PEPTIDE_SEQ", "COLS"):
        if not hasattr(mod, need):
            raise SystemExit(f"⛔ s51 lacks `{need}` -- it has changed shape; re-audit before use")
    return mod


s51 = load_s51()
SEG = re.compile(r"^([A-Za-z]+)(\d+)-(\d+)$")      # a NAMED segment, e.g. E173-181 / A1-276
TOK = re.compile(r"^([A-Za-z]+)(\d+)$")            # a single token, e.g. E173 / A35

# `J1001` §5 ⑴ -- the contig must split into exactly ONE binder element plus the THREE target
# chains of this pMHC target. ⛔ TARGET-SPECIFIC, exactly like `s51.PEPTIDE_SEQ`: it encodes
# "HLA + b2m + peptide", so a different target needs this changed DELIBERATELY, not silently.
# ⛔ Do NOT relax it to a `>= 2` test -- that was the defect this constant replaces.
N_CONTIG_ELEMENTS = 4


def contig_elements(contig):
    """Split a contig into its `/0`-delimited chain elements.

    ⛔ Asserted rather than assumed: a `/0` that is not its own comma-delimited token has been
    SWALLOWED, which silently changes the output layout. Refuse instead of guessing."""
    raw = [t.strip() for t in contig.split(",")]
    for t in raw:
        if "/0" in t and t != "/0":
            raise ValueError(f"⛔ swallowed chain break {t!r} in contig {contig!r} -- "
                             f"the `/0` is not its own comma-delimited token")
    els, cur = [], []
    for t in raw:
        if t == "/0":
            els.append(cur)
            cur = []
        elif t:
            cur.append(t)
    els.append(cur)
    return [e for e in els if e]


def roles_scaffolded(jpath, atoms):
    """Role from the contig's `/0` SEGMENTATION. Signature matches s51.roles exactly."""
    d = json.load(open(jpath))
    spec = d.get("specification") or {}
    imap = d.get("diffused_index_map") or {}
    contig = spec.get("contig") or ""
    if not imap or not contig:
        raise ValueError("spec lacks diffused_index_map or contig")

    els = contig_elements(contig)
    # `J1001` §5 ⑴ -- pre-registered, BUILT 2026-11-04. ⛔ The old test was `len(els) < 2`
    # ("a binder plus >=1 target"), which is NOT the assertion that was pre-registered and would
    # have scored an arm whose target chains had silently collapsed or multiplied.
    if len(els) != N_CONTIG_ELEMENTS:
        raise ValueError(f"contig splits into {len(els)} element(s), expected exactly "
                         f"{N_CONTIG_ELEMENTS} (binder + HLA + b2m + peptide): {contig!r}")

    # input chain -> output chain(s), exactly as s51 does it
    in2out = defaultdict(set)
    for tok, o in imap.items():
        mt, mo = TOK.match(str(tok)), TOK.match(str(o))
        if mt and mo:
            in2out[mt.group(1)].add(mo.group(1))

    def named_out(el):
        o = set()
        for tok in el:
            m = SEG.match(tok)
            if m:
                o |= in2out.get(m.group(1), set())
        return o

    binder_named = named_out(els[0])

    # `J1001` §5 ⑶ -- pre-registered, BUILT 2026-11-04. ⛔ The old code UNIONED elements 1..
    # straight into `tgt_out` and asserted only that the union was non-empty. A union hides
    # exactly the two defects worth catching: an element resolving to ZERO output chains (its
    # named tokens absent from `diffused_index_map`) and TWO elements collapsing onto the SAME
    # output chain. Both leave a non-empty union and a plausible-looking table.
    tgt_out, seen = set(), {}
    for i, el in enumerate(els[1:], start=1):
        o = named_out(el)
        if len(o) != 1:
            raise ValueError(f"target element {i} ({el!r}) maps to {sorted(o)}, expected exactly "
                             f"1 output chain (contig {contig!r})")
        c = next(iter(o))
        if c in seen:
            raise ValueError(f"target elements {seen[c]} and {i} both map to output chain {c!r} "
                             f"-- the target chains are not separable (contig {contig!r})")
        seen[c] = i
        tgt_out |= o
    if not tgt_out:
        raise ValueError(f"no target output chain resolved from elements 1.. of {contig!r}")

    present = []
    for a in atoms:
        if a["ch"] not in present:
            present.append(a["ch"])

    if binder_named:
        # SCAFFOLDED: the binder element names real input segments
        binder = sorted(binder_named)
        if len(binder) != 1:
            raise ValueError(f"binder element maps to {binder}, expected exactly 1 "
                             f"(contig {contig!r})")
        if binder[0] in tgt_out:
            raise ValueError(f"binder chain {binder[0]} is also a target chain -- "
                             f"roles are not separable from {contig!r}")
        if binder[0] not in present:
            raise ValueError(f"binder chain {binder[0]} resolved from spec is absent from the "
                             f"structure (present={present})")
    else:
        # GENERATED binder (the arm3 regime): fall back to s51's own rule, verbatim
        binder = [c for c in present if c not in tgt_out]
        if len(binder) != 1:
            raise ValueError(f"expected exactly 1 binder chain, got {binder} "
                             f"(present={present}, target={sorted(tgt_out)})")

    # peptide identified by the HARD-CODED sequence and ASSERTED; never learned from the model
    pep = None
    for c in sorted(tgt_out):
        if s51.seq_of(atoms, c) == s51.PEPTIDE_SEQ:
            pep = c
    if pep is None:
        raise ValueError(f"no target chain has the peptide sequence {s51.PEPTIDE_SEQ}")
    hla = [c for c in sorted(tgt_out) if c != pep]
    hla = max(hla, key=lambda c: len(s51.seq_of(atoms, c))) if hla else None

    _sfa = spec.get("select_fixed_atoms")           # ⛔ type is NOT stable across arms (s51:156)
    pins = list(_sfa.keys()) if isinstance(_sfa, dict) else []
    pin_out = sorted({imap[p] for p in pins if p in imap and not SEG.match(p)})
    return binder[0], sorted(tgt_out), pep, hla, pin_out


# ⭐ THE ONLY REBINDING. Every gate below is s51's own code, unchanged.
s51.roles = roles_scaffolded

files = []
for g in GLOBS:
    files += sorted(__import__("glob").glob(g))
files = [f for f in files if f.endswith(".cif.gz")]
if not files:
    raise SystemExit(f"⛔ no .cif.gz matched {GLOBS!r} -- refusing to write an empty table")

rows, bad = [], []
for f in files:
    j = f.replace(".cif.gz", ".json")
    try:
        rows.append(s51.score_one(f, j))
    except Exception as e:
        bad.append((f, f"{type(e).__name__}: {e}"))

for r in rows:
    m = re.match(r"^arm_shard(\d+)_(.+?)_(\d+)_model_(\d+)$", r["design"])
    r["arm"] = m.group(2) if m else "?"
    r["batch"] = f"s{m.group(1)}b{m.group(3)}" if m else "?"

COLS = s51.COLS
with open(OUT, "w") as fh:
    fh.write("# s75 -- scaffolded RFd3 arm. Roles from the contig's /0 segmentation; "
             "GATES ARE s51'S OWN CODE, rebound role fn only.\n")
    fh.write("# ⛔ CONTINUITY IS GATED FIRST -- a clash floor rewards a broken chain.\n")
    fh.write(f"# raw clash floor bb_min >= {s51.RAW_BB_FLOOR} A (raw RFd3 output clashes normally)\n")
    fh.write(f"# n is BATCHES, not designs. scored {len(rows)}, refused {len(bad)}\n")
    fh.write("\t".join(COLS) + "\n")
    for r in rows:
        fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")

print(f"=== s75: {len(rows)} scored, {len(bad)} refused -> {OUT}")
if bad:
    print("⚠️ ⛔ ABSENT FROM THE TSV IS NOT A PASS:")
    for f, e in bad[:10]:
        print(f"   {os.path.basename(f)}: {e}")
if rows:
    first = rows[0]
    print(f"  resolved layout: binder={first['binder']} target={first['target']} "
          f"peptide={first['pep_ch']} hla={first['hla_ch']} (⭐ asserted by sequence)")
    nb = len({r["batch"] for r in rows})
    ns = len({r["batch"][:3] for r in rows})
    print(f"  ⛔ n = {nb} BATCHES over {len(rows)} designs, {ns} shards")
    cont = sum(r["continuous"] for r in rows)
    clash = sum(r["pass_raw_clash"] for r in rows)
    fold = sum(r["one_fold"] for r in rows)
    conj = sum(1 for r in rows if r["continuous"] and r["pass_raw_clash"] and r["one_fold"])
    print(f"  G-a continuous           {cont}/{len(rows)}")
    print(f"  G-c raw clash >= {s51.RAW_BB_FLOOR}      {clash}/{len(rows)}")
    print(f"  G-e one fold             {fold}/{len(rows)}")
    print(f"  ⇒ CONJUNCTION survivors  {conj}/{len(rows)}  "
          f"(⛔ conjunction, never a mean of fractions)")
    print(f"  {'cell':<36} {'surv':>8}  {'cont_pep':>9}  {'cont_hla':>9}  {'n_pins':>6}")
    for arm in sorted({r["arm"] for r in rows}):
        sub = [r for r in rows if r["arm"] == arm]
        s = sum(1 for r in sub if r["continuous"] and r["pass_raw_clash"] and r["one_fold"])
        print(f"  {arm:<36} {s:>3}/{len(sub):<4}  "
              f"{np.mean([r['cont_pep'] for r in sub]):>9.2f}  "
              f"{np.mean([r['cont_hla'] for r in sub]):>9.2f}  "
              f"{np.mean([r['n_pins'] for r in sub]):>6.2f}")

# ---- the pre-registered equivalence gate ----
if REF:
    print(f"\n=== --compare: EXACT equivalence against {REF}")
    with open(REF) as fh:
        lines = [ln.rstrip("\n") for ln in fh if not ln.startswith("#")]
    hdr = lines[0].split("\t")
    ref = {}
    for ln in lines[1:]:
        if not ln.strip():
            continue
        v = ln.split("\t")
        ref[v[hdr.index("design")]] = dict(zip(hdr, v))
    mine = {r["design"]: r for r in rows}
    only_ref = sorted(set(ref) - set(mine))
    only_mine = sorted(set(mine) - set(ref))
    shared = sorted(set(ref) & set(mine))
    diffs = []
    for dn in shared:
        for c in hdr:
            if c in ("arm", "batch", "design"):
                continue
            a, b = str(ref[dn].get(c, "")), str(mine[dn].get(c, ""))
            if a != b:
                diffs.append((dn, c, a, b))
    print(f"  rows: reference {len(ref)}, mine {len(rows)}, shared {len(shared)}")
    print(f"  only in reference: {len(only_ref)}   only in mine: {len(only_mine)}")
    print(f"  per-cell disagreements: {len(diffs)}")
    for dn, c, a, b in diffs[:12]:
        print(f"    {dn}  {c}: ref={a!r} mine={b!r}")
    ok = (not diffs) and (not only_ref) and (not only_mine)
    print(f"  ⇒ {'✅ EQUIVALENT -- the gate PASSES' if ok else '⛔⛔ NOT EQUIVALENT -- DO NOT SHIP'}")
    if not ok:
        sys.exit(1)
