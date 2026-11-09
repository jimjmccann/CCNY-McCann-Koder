#!/usr/bin/env python3
"""s114 -- the option-③ designability 2x2, scored at the BATCH level with the held pin EXCLUDED.

⛔ READ-ONLY on every input. Nothing here launches, spends, or edits a table.

⛔⛔ WHY THIS EXISTS AND WHY `s88` COULD NOT BE USED (measured 2026-11-08, not supposed).
   `s88` is the committed instrument for "recompute `cont_pep` with the held pins excluded", and it
   is the right *logic*. But it resolves roles with **`s51.roles`**, and `s51`'s rule is exactly the
   one that CANNOT read a scaffolded contig (`J1001`; that is the whole reason `s75` exists). Run on
   this round it raised

       ValueError: expected exactly 1 binder chain, got [] (present=['A','B','C','D'], ...)

   on **all 160 motif-arm designs** and scored only the **160 no-motif** ones -- the half with
   **ZERO pins**. Its report then read perfectly clean and concluded *"CLAIM 3 IS FALSE: no pin lands
   in the binder chain"* and *"the subtraction is exact, 160/160"*, both of which are artefacts of
   having silently dropped every design that HAS a pin
   ([[a-clean-empty-report-is-a-failure-mode]], [[never-rank-on-partial-data]]).
   ⇒ ⛔ that output must not be cited. `s88` is left BYTE-IDENTICAL on disk because `J1016` was
   adjudicated on it; this file is the scaffolded-capable equivalent.

⭐⭐ THE CONSTRUCTION: this file does NOT re-implement a gate or a metric. It imports **`s75`**,
   whose module body rebinds `s51.roles = roles_scaffolded`, so `s75.s51` IS `s51` with the
   scaffolded resolver installed and every metric object (`read_cif`, `contacts`, `heavy_of`,
   `CONTACT`) is `s51`'s own. `binder_res` is copied from `s88` and ASSERTED against `s88`'s source
   at run time, so a silent divergence from the instrument this correction comes from is impossible.
   ⚠️ `s75` scores at import, so it is imported with `sys.argv` stubbed to a throwaway `--out` and a
   glob matching ONE file, with stdout captured. ⛔ That temp file is NOT a result and is deleted.

⛔⛔ WHY THE EXCLUSION IS LOAD-BEARING HERE, i.e. why this is not bookkeeping.
   `[M]` The motif arms hold **n_pins = 1** -- `E101`, whose `CG,OD1,ND2` are held FROM THE INPUT and
   were placed by construction to hydrogen-bond the peptide. The no-motif arms hold **0**. So
   `cont_pep` gives the motif arms up to **+1.0** contact from an ATOM SET THE OTHER ARMS DO NOT HAVE
   ([[a-metric-that-counts-held-input-cannot-compare-arms]]). Raw `cont_pep` is 5.88/5.86 for the
   motif arms against 5.28/6.24 -- so a ~1.0 correction does not shrink the motif factor, it can
   INVERT it. ⇒ the uncorrected comparison is not a weaker version of the answer, it is a different
   one, and `s51.contacts` being RESIDUE-level bounds one pin's contribution at exactly 1.

⛔⛔ `n` IS BATCHES, NOT DESIGNS ([[rfd3-designs-are-clustered-by-batch]]): all 4 designs in a
   diffusion batch share one sampled contig. 20 batches per arm, 40 per factor level. Per-design
   means are reported for comparability only; ⛔ every TEST is on the 20 batch-level values, because
   80 designs on 20 contigs are not 80 observations (the pseudo-replication that nearly shipped a
   3.4x enrichment off one backbone).

⛔ THIS IS THE DESIGN-STAGE GATE, NOT THE PRE-REGISTERED ENDPOINT. `J1055` §3 registers the
   falsifier on the REFOLD conjunction (peptide contacts >= an in-batch cross-paired non-binder
   floor AND `pair(3,2)` above that floor). Nothing here can satisfy or refute it; survival and
   contact counts are design-stage quantities only.

Usage:
  python3 -I s114_designability_2x2.py <t149.tsv> <out.md> <glob> [<glob> ...]
"""
import collections
import contextlib
import csv
import glob as globmod
import importlib.util
import io
import os
import random
import re
import statistics as st
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S75 = os.path.join(HERE, "s75_score_scaffolded_arm.py")
S88 = os.path.join(HERE, "s88_cont_pep_excluding_pins.py")

if len(sys.argv) < 4:
    sys.exit(__doc__)
REF, OUT, GLOBS = sys.argv[1], sys.argv[2], sys.argv[3:]

files = []
for g in GLOBS:
    files += globmod.glob(g)
files = sorted({f for f in files if f.endswith(".cif.gz")})
if not files:
    sys.exit(f"⛔ no .cif.gz matched {GLOBS!r} -- refusing to write an empty report")

TOK = re.compile(r"^([A-Za-z]+)(\d+)$")


def load_s75(one_file):
    """Import s75 with argv stubbed so its module-level loop scores ONE throwaway design.

    ⛔ s75 raises SystemExit on an empty glob (it refuses to write an empty table), so unlike
    `s88`'s treatment of `s51` it cannot be stubbed with no files at all."""
    fd, tmp = tempfile.mkstemp(prefix="s114_s75_import_", suffix=".tsv")
    os.close(fd)
    spec = importlib.util.spec_from_file_location("s75_mod", S75)
    mod = importlib.util.module_from_spec(spec)
    saved = sys.argv
    sys.argv = ["s75_score_scaffolded_arm.py", tmp, one_file]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
        try:
            os.unlink(tmp)
        except OSError:
            pass
    if not hasattr(mod, "s51") or not hasattr(mod, "roles_scaffolded"):
        raise SystemExit("⛔ s75 has changed shape -- re-audit before use")
    if mod.s51.roles is not mod.roles_scaffolded:
        raise SystemExit("⛔ s75 did not rebind s51.roles -- the scaffolded arms would refuse")
    for need in ("read_cif", "contacts", "heavy_of", "CONTACT"):
        if not hasattr(mod.s51, need):
            raise SystemExit(f"⛔ s51 lacks `{need}` -- re-audit before use")
    return mod.s51


s51 = load_s75(files[0])


# ⛔ copied from `s88`:101-108 and ASSERTED against its source, so the correction cannot silently
#   diverge from the instrument it is taken from.
def binder_res(atoms, ch):
    g = collections.defaultdict(list)
    for a in atoms:
        if a["ch"] == ch and a["el"] != "H":
            g[a["resi"]].append(a["xyz"])
    return {k: np.array(v) for k, v in sorted(g.items())}


_needle = 'if a["ch"] == ch and a["el"] != "H":'
if _needle not in open(S88).read():
    sys.exit("⛔ REFUSING: `binder_res` no longer matches s88's source -- the pin exclusion would "
             "not be the one J1016 was corrected with, and the drift would be silent.")

ref = {}
with open(REF) as fh:
    for r in csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"):
        ref[r["design"]] = r

rows, refused, mismatch = [], [], []
for cif in files:
    jpath = cif.replace(".cif.gz", ".json")
    name = os.path.basename(cif).replace(".cif.gz", "")
    if not os.path.exists(jpath):
        refused.append((name, "no sidecar json"))
        continue
    try:
        atoms = s51.read_cif(cif)
        binder, tgt, pep, hla, pin_out = s51.roles(jpath, atoms)
    except Exception as e:                                       # noqa: BLE001
        refused.append((name, f"{type(e).__name__}: {e}"))
        continue

    br = binder_res(atoms, binder)
    pep_at = s51.heavy_of(atoms, pep)

    pin_resis, off_chain = set(), []
    for p in pin_out:
        m = TOK.match(str(p))
        if not m:
            continue
        if m.group(1) == binder:
            pin_resis.add(int(m.group(2)))
        else:
            off_chain.append(p)

    cp_all = s51.contacts([br[k] for k in sorted(br)], pep_at)
    cp_exc = s51.contacts([br[k] for k in sorted(br) if k not in pin_resis], pep_at)
    n_touch = s51.contacts([br[k] for k in sorted(br) if k in pin_resis], pep_at)

    r = ref.get(name)
    if r is None:
        refused.append((name, "absent from the reference table"))
        continue
    if int(float(r["cont_pep"])) != cp_all:
        mismatch.append((name, r["cont_pep"], cp_all))

    rows.append({
        "design": name, "arm": r["arm"], "batch": r["batch"],
        "motif": "held" if r["arm"].endswith(("_motif_narrow", "_motif_broad")) else "none",
        "breadth": "broad" if r["arm"].endswith("broad") else "narrow",
        "n_pins_binder": len(pin_resis), "n_pins_offchain": len(off_chain),
        "cont_pep_all": cp_all, "cont_pep_excl": cp_exc, "n_pins_touching": n_touch,
        "cont_hla": float(r["cont_hla"]),
        "surv": int(r["continuous"]) and int(r["pass_raw_clash"]) and int(r["one_fold"]),
        "binder_len": int(r["binder_len"]), "helix_frac": float(r["helix_frac"]),
        "rg_ratio": float(r["rg_ratio"]),
    })


def perm_p(a, b, n=20000, seed=20261108):
    """Two-sided permutation test on the difference of means. ⛔ Assumption-free: at n=40 per level
    a normal approximation is not worth defending, and the labels are exchangeable under the null."""
    a, b = list(a), list(b)
    obs = abs(st.mean(a) - st.mean(b))
    pool = a + b
    na = len(a)
    rnd = random.Random(seed)
    hits = 0
    for _ in range(n):
        rnd.shuffle(pool)
        if abs(st.mean(pool[:na]) - st.mean(pool[na:])) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (n + 1)


L = []
W = L.append
W("# s114 -- the designability 2x2, at the BATCH level, held pin EXCLUDED")
W("")
W(f"- `[M]` **{len(files)} design files**, **{len(rows)} scored**, {len(refused)} refused")
W(f"- roles from `s75.roles_scaffolded`; metrics are `s51`'s own, `CONTACT` = **{s51.CONTACT} Å**, "
  f"residue-level")
W(f"- reference `{os.path.basename(REF)}` ({len(ref)} rows)")
W("- ⛔ DESIGN-STAGE ONLY. The `J1055` §3 falsifier lives on the REFOLD conjunction and is "
  "untouched by anything here.")
if refused:
    W("")
    W("⛔ **REFUSED**")
    for n, why in refused[:12]:
        W(f"  - `{n}` — {why}")

# ── 0 · the mandatory control ─────────────────────────────────────────────────────────────────
W("")
W("## 0 · ⛔⛔ MANDATORY CONTROL — does `cont_pep_all` reproduce `t149` EXACTLY?")
W("")
W(f"`[M]` **{len(rows)}** designs checked; **{len(mismatch)}** disagree.")
if mismatch or len(rows) != len(files):
    W("")
    W("⛔⛔ **REFUSING TO REPORT.** Either a design was dropped or the recomputed metric is not "
      "`s51`'s own quantity, and a clean number from the wrong quantity is worse than an error.")
    for n, a, b in mismatch[:10]:
        W(f"  - `{n}`: reference {a} vs recomputed {b}")
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    sys.exit(1)
W("⇒ ✅ **EXACT on every design, all four arms** — the quantity below is `s51`'s `cont_pep`, "
  "minus the held pin.")

# ── 1 · the confound, measured ────────────────────────────────────────────────────────────────
W("")
W("## 1 · THE CONFOUND, MEASURED PER ARM (⛔ not assumed, and it is not zero)")
W("")
W("| arm | n | pins in binder | `cont_pep` raw | **excl. pin** | Δ | pin itself touching |")
W("|---|---|---|---|---|---|---|")
by_arm = collections.defaultdict(list)
for r in rows:
    by_arm[r["arm"]].append(r)
for arm in sorted(by_arm):
    g = by_arm[arm]
    pins = sorted({x["n_pins_binder"] for x in g})
    ra, ex = st.mean(x["cont_pep_all"] for x in g), st.mean(x["cont_pep_excl"] for x in g)
    tch = st.mean(x["n_pins_touching"] for x in g)
    W(f"| `{arm}` | {len(g)} | {pins} | {ra:.3f} | **{ex:.3f}** | {ex - ra:+.3f} | {tch:.3f} |")
W("")
drop_ok = all(x["cont_pep_all"] - x["cont_pep_excl"] == x["n_pins_touching"] for x in rows)
W(f"- `[M]` `cont_pep_all - cont_pep_excl == n_pins_touching` in "
  f"**{sum(1 for x in rows if x['cont_pep_all'] - x['cont_pep_excl'] == x['n_pins_touching'])}"
  f"/{len(rows)}** designs ⇒ {'✅ the subtraction is EXACT, not a bound' if drop_ok else '⛔ NOT exact'}")
_m = [x for x in rows if x["motif"] == "held"]
W(f"- `[M]` on the motif arms the held pin is itself within {s51.CONTACT} Å of the peptide in "
  f"**{sum(x['n_pins_touching'] for x in _m)}/{len(_m)}** designs "
  f"({100.0 * sum(x['n_pins_touching'] for x in _m) / len(_m):.1f}%) ⇒ that is the free contact.")

# ── 2 · batch-level aggregation ───────────────────────────────────────────────────────────────
batches = collections.defaultdict(list)
for r in rows:
    batches[(r["arm"], r["batch"])].append(r)

B = []
for (arm, b), g in sorted(batches.items()):
    surv = [x for x in g if x["surv"]]
    B.append({
        "arm": arm, "batch": b, "motif": g[0]["motif"], "breadth": g[0]["breadth"],
        "n": len(g), "surv_frac": len(surv) / len(g),
        "cp_excl_surv": st.mean(x["cont_pep_excl"] for x in surv) if surv else None,
        "chla_surv": st.mean(x["cont_hla"] for x in surv) if surv else None,
        "cp_excl_all": st.mean(x["cont_pep_excl"] for x in g),
    })

W("")
W("## 2 · THE 2x2, CELL BY CELL (n = 20 BATCHES per cell, ⛔ never 80 designs)")
W("")
W("| cell | motif | breadth | batches | survival | `cont_pep` excl (survivors) | "
  "`cont_hla` (survivors) |")
W("|---|---|---|---|---|---|---|")
for arm in sorted(by_arm):
    g = [x for x in B if x["arm"] == arm]
    sv = st.mean(x["surv_frac"] for x in g)
    cp = [x["cp_excl_surv"] for x in g if x["cp_excl_surv"] is not None]
    ch = [x["chla_surv"] for x in g if x["chla_surv"] is not None]
    W(f"| `{arm}` | {g[0]['motif']} | {g[0]['breadth']} | {len(g)} | "
      f"**{100 * sv:.1f}%** | {st.mean(cp):.2f} | {st.mean(ch):.2f} |")

# ── 3 · main effects ──────────────────────────────────────────────────────────────────────────
W("")
W("## 3 · MAIN EFFECTS — permutation test on the 40 vs 40 batch-level values")
W("")
W("| outcome | factor | level A | level B | mean A | mean B | Δ | p |")
W("|---|---|---|---|---|---|---|---|")
for oname, key in (("survival", "surv_frac"),
                   ("cont_pep excl (survivors)", "cp_excl_surv"),
                   ("cont_pep excl (all)", "cp_excl_all"),
                   ("cont_hla (survivors)", "chla_surv")):
    for fname, la, lb in (("motif", "held", "none"), ("breadth", "broad", "narrow")):
        a = [x[key] for x in B if x[fname] == la and x[key] is not None]
        b = [x[key] for x in B if x[fname] == lb and x[key] is not None]
        p = perm_p(a, b)
        W(f"| {oname} | {fname} | {la} (n={len(a)}) | {lb} (n={len(b)}) | {st.mean(a):.3f} | "
          f"{st.mean(b):.3f} | {st.mean(a) - st.mean(b):+.3f} | **{p:.4f}** |")

# ── 4 · the interaction, described and NOT over-tested ────────────────────────────────────────
W("")
W("## 4 · THE INTERACTION — described, ⛔ not tested")
W("")
W("⛔ At 20 batches per cell an interaction test is underpowered and would be the third test on "
  "the same 80 batches. The cell means are given so the SHAPE is visible; ⛔ no p-value is quoted.")
W("")
W("| | narrow | broad |")
W("|---|---|---|")
for mo in ("held", "none"):
    cells = []
    for br in ("narrow", "broad"):
        g = [x["surv_frac"] for x in B if x["motif"] == mo and x["breadth"] == br]
        cells.append(f"{100 * st.mean(g):.1f}%")
    W(f"| **motif {mo}** | {cells[0]} | {cells[1]} |")

# ── 5 · what else moved, so a single gate is not read as the whole story ──────────────────────
W("")
W("## 5 · THE OTHER DESIGN-STAGE DESCRIPTORS (per design, context only)")
W("")
W("| arm | binder_len median | helix_frac | rg_ratio | continuity | raw-clash pass |")
W("|---|---|---|---|---|---|")
for arm in sorted(by_arm):
    g = by_arm[arm]
    rr = [ref[x["design"]] for x in g]
    W(f"| `{arm}` | {st.median(x['binder_len'] for x in g):.0f} | "
      f"{st.mean(x['helix_frac'] for x in g):.3f} | {st.mean(x['rg_ratio'] for x in g):.3f} | "
      f"{sum(int(x['continuous']) for x in rr)}/{len(g)} | "
      f"{sum(int(x['pass_raw_clash']) for x in rr)}/{len(g)} |")

TSV = OUT.replace(".md", "") + "_per_batch.tsv"
with open(TSV, "w") as fh:
    cols = ["arm", "batch", "motif", "breadth", "n", "surv_frac", "cp_excl_surv", "chla_surv",
            "cp_excl_all"]
    fh.write("\t".join(cols) + "\n")
    for x in B:
        fh.write("\t".join("" if x[c] is None else str(x[c]) for c in cols) + "\n")
W("")
W(f"- wrote `{os.path.basename(TSV)}` — {len(B)} BATCH-level rows, the unit every test above uses")

open(OUT, "w").write("\n".join(L) + "\n")
print("\n".join(L))
