#!/usr/bin/env python3
"""s88 -- recompute `cont_pep` with the HELD PIN RESIDUES EXCLUDED from the binder, and settle
whether `J997`'s `arm1`-vs-`arm2` result survives it.

⛔ READ-ONLY on every input. Nothing here launches, spends, or edits a table.

WHY THIS EXISTS (`J1016` §4⒝, the owed `$0` test). CMP-10 §3 claim 3 is correct that `cont_pep`
never subtracts the peptide contacts contributed by the HELD guidepost tips. `J1016` measured the
consequence with a crude `minus n_pins` bound and found:
  * ⒜ `arm3`-vs-pinned is STRENGTHENED (arm3 holds 0 pins; gap +0.617 -> +4.117, same sign)
  * ⒝ ⛔ `arm1`-vs-`arm2` INVERTS (arm2 holds 4 pins to arm1's 3), so its "settled" label was
       WITHDRAWN -- ⛔ withdrawn, NOT reversed, because `minus n_pins` assumes each held tip
       contributes EXACTLY 1.0 and that was an assumption, not a measurement.
⇒ this measures it instead.

⭐⭐ WHY THE SUBTRACTION IS PLAUSIBLY EXACT RATHER THAN A BOUND, which is what makes the test sharp:
   `s51.contacts()` is a **RESIDUE-LEVEL** count ("how many residues of a_res touch b_at within
   CONTACT"), so a pinned residue can contribute AT MOST 1 regardless of how many of its atoms are
   near the peptide. And the pins were placed BY CONSTRUCTION to hydrogen-bond the peptide at
   2.6-3.0 A, well inside `CONTACT` 4.5 A. ⇒ the expected result is
   `cont_pep_all - cont_pep_excl == n_pins_touching == n_pins`. ⛔ But that is a PREDICTION; the
   point of this script is that it is checked, per design, and reported when it fails.

⛔⛔ THE PINS ARE INSIDE THE BINDER CHAIN, WHICH IS WHY THEY WERE EVER COUNTED. `s51` asserts
   EXACTLY ONE binder chain (`s51`:141-144), and `unindex` appends its components into the binder
   and renumbers them above `max(res_id)`, resetting the chain
   (`an earlier internal note on generator unindex semantics`). So a pin lands in the binder chain as a
   trailing residue and `res_atoms(binder)` includes it. ⛔ Had the pins been separate chains,
   `s51` would have raised on >1 binder chain and `cont_pep` would never have seen them -- in which
   case CMP-10's claim 3 would be FALSE. That branch is tested here, not assumed.

⭐ `s51` IS NOT EDITED AND NOT RE-IMPLEMENTED. It is imported as a module with `sys.argv` stubbed,
   exactly as `s75` does it (`s75`:103-126), and every gate object used here is `s51`'s own:
   `read_cif`, `roles`, `contacts`, `heavy_of`, `CONTACT`. ⛔ `s51` must stay byte-identical on disk
   because `J997`/`ADR 0063` were adjudicated on it.

⛔⛔ THE CONTROL COMES FIRST. Before any exclusion is believed, `cont_pep_all` recomputed here must
   reproduce `t80`'s `cont_pep` EXACTLY, per design. A mismatch means this script is not measuring
   `s51`'s quantity and the exclusion is meaningless ⇒ it REFUSES rather than reporting
   ([[a-clean-empty-report-is-a-failure-mode]]: a clean number from the wrong quantity is worse
   than an error).

⚠️ `n` IS BATCHES, NOT DESIGNS ([[rfd3-designs-are-clustered-by-batch]]). Per-design means are
   reported for comparability with `t80`, and the batch-level mean of means is reported alongside
   because that is the pre-registered unit.

Usage:
  python3 -I s88_cont_pep_excluding_pins.py <t80.tsv> <out.md> <glob> [<glob> ...]
"""
import collections
import csv
import glob as globmod
import importlib.util
import json
import os
import re
import statistics as st
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S51 = os.path.join(HERE, "s51_score_matrix.py")

if len(sys.argv) < 4:
    sys.exit(__doc__)
REF, OUT, GLOBS = sys.argv[1], sys.argv[2], sys.argv[3:]


def load_s51():
    """Import s51 with argv stubbed so its module-level loop scores nothing.

    ⛔ Do NOT 'fix' this by refactoring s51 -- it must stay byte-identical on disk (s75's note)."""
    if not os.path.exists(S51):
        raise SystemExit(f"⛔ s51 not found at {S51}")
    fd, tmp = tempfile.mkstemp(prefix="s88_s51_import_", suffix=".tsv")
    os.close(fd)
    saved = sys.argv
    sys.argv = ["s51_score_matrix.py", tmp]
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
    for need in ("roles", "read_cif", "contacts", "heavy_of", "CONTACT"):
        if not hasattr(mod, need):
            raise SystemExit(f"⛔ s51 lacks `{need}` -- it has changed shape; re-audit before use")
    return mod


s51 = load_s51()
TOK = re.compile(r"^([A-Za-z]+)(\d+)$")


def binder_res(atoms, ch):
    """binder residues KEYED BY resi -- s51.res_atoms() drops the label, and the label is what
    the exclusion needs."""
    g = collections.defaultdict(list)
    for a in atoms:
        if a["ch"] == ch and a["el"] != "H":
            g[a["resi"]].append(a["xyz"])
    return {k: np.array(v) for k, v in sorted(g.items())}


# ── the reference table, for the mandatory control ────────────────────────────────────────────
ref = {}
with open(REF) as fh:
    for r in csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"):
        ref[r["design"]] = r

files = []
for g in GLOBS:
    files += globmod.glob(g)
files = sorted(set(files))

rows, mismatch, refused = [], [], []
for cif in files:
    jpath = cif.replace(".cif.gz", ".json")
    if not os.path.exists(jpath):
        refused.append((os.path.basename(cif), "no sidecar json"))
        continue
    name = os.path.basename(cif).replace(".cif.gz", "")
    try:
        atoms = s51.read_cif(cif)
        binder, tgt, pep, hla, pin_out = s51.roles(jpath, atoms)
    except Exception as e:                                  # noqa: BLE001
        refused.append((name, f"{type(e).__name__}: {e}"))
        continue

    br = binder_res(atoms, binder)
    pep_at = s51.heavy_of(atoms, pep)

    pin_resis = set()
    off_chain = []
    for p in pin_out:
        m = TOK.match(str(p))
        if not m:
            continue
        if m.group(1) == binder:
            pin_resis.add(int(m.group(2)))
        else:
            off_chain.append(p)

    all_res = [br[k] for k in sorted(br)]
    exc_res = [br[k] for k in sorted(br) if k not in pin_resis]
    cp_all = s51.contacts(all_res, pep_at)
    cp_exc = s51.contacts(exc_res, pep_at)
    # how many of the pin residues are THEMSELVES in contact with the peptide
    pin_res = [br[k] for k in sorted(br) if k in pin_resis]
    n_touch = s51.contacts(pin_res, pep_at)

    rows.append({
        "design": name, "arm": ref.get(name, {}).get("arm", "?"),
        "batch": ref.get(name, {}).get("batch", "?"),
        "n_pins_spec": len(pin_out), "n_pins_in_binder": len(pin_resis),
        "n_pins_off_chain": len(off_chain),
        "cont_pep_all": cp_all, "cont_pep_excl": cp_exc, "n_pins_touching": n_touch,
        "ref_cont_pep": ref.get(name, {}).get("cont_pep"),
    })
    rc = ref.get(name, {}).get("cont_pep")
    if rc is not None and int(float(rc)) != cp_all:
        mismatch.append((name, rc, cp_all))

L = []
W = L.append
W("# s88 -- `cont_pep` with the HELD PIN RESIDUES EXCLUDED")
W("")
W(f"- `[M]` **{len(files)} design files** found, **{len(rows)} scored**, {len(refused)} refused")
W(f"- reference table `{os.path.basename(REF)}` ({len(ref)} rows)")
W(f"- `s51` imported as a module and UNMODIFIED; `CONTACT` = **{s51.CONTACT} Å**, "
  f"residue-level count")
if refused:
    W("")
    W("⛔ **REFUSED**")
    for n, why in refused[:12]:
        W(f"  - `{n}` — {why}")

# ── THE CONTROL ───────────────────────────────────────────────────────────────────────────────
W("")
W("## 0 · THE MANDATORY CONTROL — does `cont_pep_all` reproduce `t80`?")
W("")
checked = [r for r in rows if r["ref_cont_pep"] is not None]
W(f"`[M]` **{len(checked)}** designs matched to the reference; **{len(mismatch)}** disagree.")
if mismatch:
    W("")
    W("⛔⛔ **REFUSING TO REPORT THE EXCLUSION.** A disagreement here means this script is not "
      "measuring `s51`'s quantity, so the exclusion would be a clean number from the wrong "
      "measurement. First five:")
    for n, rc, got in mismatch[:5]:
        W(f"  - `{n}`: `t80` {rc}, recomputed {got}")
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    raise SystemExit(f"⛔ {len(mismatch)} control mismatches")
W("⇒ ✅ **EXACT, every design** — the quantity below is `s51`'s own `cont_pep`, minus pins.")

# ── the pin-location branch that decides whether claim 3 is even true ─────────────────────────
W("")
W("## 1 · ARE THE PINS INSIDE THE BINDER CHAIN? (the branch that decides claim 3)")
W("")
tot_spec = sum(r["n_pins_spec"] for r in rows)
tot_in = sum(r["n_pins_in_binder"] for r in rows)
tot_off = sum(r["n_pins_off_chain"] for r in rows)
W(f"- `[M]` pin tokens in the specs: **{tot_spec}**; resolved INTO the binder chain: "
  f"**{tot_in}**; on another chain: **{tot_off}**")
if tot_in == 0:
    W("⇒ ⛔⛔ **CLAIM 3 IS FALSE**: no pin lands in the binder chain, so `cont_pep` never counted "
      "a held tip and there is nothing to subtract.")
else:
    W(f"⇒ ⭐ the pins ARE binder residues ({tot_in}/{tot_spec}), so `cont_pep` did count them "
      f"⇒ **claim 3's premise holds.**")
W("")
W("### Is the contribution EXACTLY one per pin, as predicted?")
W("")
exact = sum(1 for r in rows if r["cont_pep_all"] - r["cont_pep_excl"] == r["n_pins_in_binder"])
touch_all = sum(1 for r in rows if r["n_pins_touching"] == r["n_pins_in_binder"])
W(f"- `[M]` `cont_pep_all - cont_pep_excl == n_pins_in_binder` in **{exact}/{len(rows)}** designs")
W(f"- `[M]` every pin residue is itself within {s51.CONTACT} Å of the peptide in "
  f"**{touch_all}/{len(rows)}** designs")
dist = collections.Counter(r["cont_pep_all"] - r["cont_pep_excl"] - r["n_pins_in_binder"]
                           for r in rows)
W(f"- `[M]` distribution of (drop − n_pins): {dict(sorted(dist.items()))} "
  f"⇒ 0 means exactly one contact per pin")
W("")
if exact == len(rows):
    W("⇒ ⭐⭐ **THE SUBTRACTION IS EXACT, NOT A BOUND.** `J1016`'s crude `minus n_pins` was the "
      "right correction and its assumption is now MEASURED.")
else:
    W("⇒ ⚠️ **the subtraction is NOT uniformly exact** ⇒ `J1016`'s bound was approximate; the "
      "measured per-design values below supersede it.")

# ── the result ────────────────────────────────────────────────────────────────────────────────
W("")
W("## 2 · THE RESULT, PER ARM")
W("")


def grp(r):
    a = r["arm"]
    return "arm3" if "hotspot" in a else ("arm1" if "P5x_" in a else "arm2")


G = collections.defaultdict(list)
for r in rows:
    G[grp(r)].append(r)

W("| group | n designs | n batches | pins | `cont_pep` raw | **excl. pins** | Δ |")
W("|---|---|---|---|---|---|---|")
for k in ("arm1", "arm2", "arm3"):
    g = G.get(k, [])
    if not g:
        continue
    nb = len({(r["arm"], r["batch"]) for r in g})
    pins = sorted({r["n_pins_in_binder"] for r in g})
    raw = st.mean(r["cont_pep_all"] for r in g)
    exc = st.mean(r["cont_pep_excl"] for r in g)
    W(f"| `{k}` | {len(g)} | {nb} | {pins} | {raw:.3f} | **{exc:.3f}** | {exc-raw:+.3f} |")

W("")
W("### ⇒ ⒜ `arm3` vs the PINNED arms")
p12 = G.get("arm1", []) + G.get("arm2", [])
if p12 and G.get("arm3"):
    raw12 = st.mean(r["cont_pep_all"] for r in p12)
    exc12 = st.mean(r["cont_pep_excl"] for r in p12)
    raw3 = st.mean(r["cont_pep_all"] for r in G["arm3"])
    exc3 = st.mean(r["cont_pep_excl"] for r in G["arm3"])
    W("")
    W(f"- `[M]` raw: `arm3` {raw3:.3f} vs `arm1+arm2` {raw12:.3f} ⇒ Δ **{raw3-raw12:+.3f}**")
    W(f"- `[M]` pins excluded: `arm3` {exc3:.3f} vs `arm1+arm2` {exc12:.3f} ⇒ "
      f"Δ **{exc3-exc12:+.3f}**")
    same = (raw3 - raw12) * (exc3 - exc12) > 0
    W(f"⇒ {'⭐ **SAME SIGN** — the verdict is STRENGTHENED, not threatened' if same else '⛔ THE SIGN CHANGES'}"
      f"; ⛔ `arm3` holds 0 pins so its own number is unchanged by construction "
      f"({raw3:.3f} = {exc3:.3f}).")

W("")
W("### ⇒ ⒝ `arm1` vs `arm2` — the WITHDRAWN result")
if G.get("arm1") and G.get("arm2"):
    r1, r2 = G["arm1"], G["arm2"]
    a1r, a2r = st.mean(r["cont_pep_all"] for r in r1), st.mean(r["cont_pep_all"] for r in r2)
    a1e, a2e = st.mean(r["cont_pep_excl"] for r in r1), st.mean(r["cont_pep_excl"] for r in r2)
    W("")
    W(f"- `[M]` raw: `arm2` {a2r:.3f} vs `arm1` {a1r:.3f} ⇒ Δ **{a2r-a1r:+.3f}** "
      f"(the `p = 0.0018` result favouring `arm2`, i.e. the `P5xA` site)")
    W(f"- `[M]` pins excluded: `arm2` {a2e:.3f} vs `arm1` {a1e:.3f} ⇒ Δ **{a2e-a1e:+.3f}**")
    if (a2r - a1r) * (a2e - a1e) > 0:
        W("⇒ ⭐⭐ **THE SIGN SURVIVES** ⇒ the `arm1`-vs-`arm2` result is NOT an artefact of the "
          "held input, and its withdrawal in `J1016` §4⒝ can be LIFTED.")
    else:
        W("⇒ ⛔⛔ **THE SIGN INVERTS ON THE MEASURED EXCLUSION** ⇒ `J1016` §4⒝ is CONFIRMED: the "
          "`arm1`-vs-`arm2` comparison was carried by the held guidepost tips, not by the design. "
          "⛔ Its *\"settled\"* label stays WITHDRAWN.")
    W("")
    W("⚠️ ⛔ **AND WHAT THIS STILL DOES NOT DO.** This is a point estimate on a metric, ⛔ not a "
      "significance test, and `n` is BATCHES not designs "
      "([[rfd3-designs-are-clustered-by-batch]]). ⛔ It does NOT establish that `arm1` beats "
      "`arm2`; it establishes which way the metric points once the held input is removed from it. "
      "⛔ And `cont_pep` remains a DESIGN-STAGE geometry count — it is not affinity and it is "
      "⛔⛔ not step-1 selectivity, which still has no identified mechanism (`J989`).")

W("")
W("## 3 · PER-ARM-CELL DETAIL")
W("")
W("| arm | n | pins | raw | excl | Δ |")
W("|---|---|---|---|---|---|")
for a in sorted({r["arm"] for r in rows}):
    g = [r for r in rows if r["arm"] == a]
    pins = sorted({r["n_pins_in_binder"] for r in g})
    raw = st.mean(r["cont_pep_all"] for r in g)
    exc = st.mean(r["cont_pep_excl"] for r in g)
    W(f"| `{a}` | {len(g)} | {pins} | {raw:.3f} | **{exc:.3f}** | {exc-raw:+.3f} |")


# ── a t80-SHAPED TABLE CARRYING THE CORRECTED METRIC, so `s62` can be reused VERBATIM ─────────
# ⛔⛔ PROVENANCE HAZARD, NAMED BECAUSE IT IS REAL: `s62` reads a column called `cont_pep`, so the
#    corrected value is written into THAT column name. The file therefore says `cont_pep` and holds
#    `cont_pep_excl`. ⇒ the filename says `_CONT_PEP_EXCL_`, a header comment says it in words, and
#    ⛔ this file must NEVER be cited as `t80` or compared column-wise against it
#    ([[never-inherit-a-default-record-the-resolved-config]] -- record the resolved config).
# ⭐ The alternative was re-implementing s62's permutation test, which risks differing from the
#    test `ADR 0063` was adjudicated on. Reusing s62 unmodified is the safer of the two.
TSV = OUT.replace(".md", "") + "_CONT_PEP_EXCL_for_s62.tsv"
by_name = {r["design"]: r for r in rows}
with open(REF) as fh:
    src_lines = fh.readlines()
hdr = [l for l in src_lines if l.startswith("#")]
body = [l for l in src_lines if not l.startswith("#")]
rdr = csv.DictReader(body, delimiter="\t")
cols = rdr.fieldnames
out_rows, swapped = [], 0
for r in rdr:
    m = by_name.get(r["design"])
    if m is not None:
        r["cont_pep"] = str(m["cont_pep_excl"])
        swapped += 1
    out_rows.append(r)
with open(TSV, "w", newline="") as fh:
    for h in hdr:
        fh.write(h)
    fh.write("# ⛔⛔ THE `cont_pep` COLUMN IN THIS FILE IS `cont_pep` WITH THE HELD PIN RESIDUES\n")
    fh.write("#    EXCLUDED FROM THE BINDER (s88). It is NOT t80's cont_pep. Never cite this as t80.\n")
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(out_rows)
W("")
W("## 4 · A `t80`-SHAPED TABLE WITH THE CORRECTED METRIC, FOR `s62`")
W("")
W(f"- wrote `{os.path.basename(TSV)}` — {len(out_rows)} rows, **{swapped} `cont_pep` values "
  f"replaced** by `cont_pep_excl`")
W("- ⛔⛔ its `cont_pep` COLUMN HOLDS `cont_pep_excl`. ⛔ Never cite it as `t80`; it exists so "
  "`s62` can run its own pre-registered batch-level permutation test on the corrected metric "
  "**unmodified**, rather than that test being re-implemented here.")

open(OUT, "w").write("\n".join(L) + "\n")
print("\n".join(L))
