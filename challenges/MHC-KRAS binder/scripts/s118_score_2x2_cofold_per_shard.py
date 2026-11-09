#!/usr/bin/env python3
"""s118 -- score the J1056 2x2 refold by running `s60` ONCE PER SHARD, then aggregate the 2x2.

⛔ READ-ONLY on every input. Nothing launches, spends, or edits a table.

⛔⛔ WHY PER SHARD, AND WHY POOLING WOULD HAVE BEEN WRONG (verified at `s60`'s source, not supposed).
   `s60`:161 collects EVERY `NONBINDER*` row under whatever root it is handed, and `s60`:177-178
   take ONE global floor as the max over them. ⇒ handed the pooled tree, it would build a single
   floor from **8 non-binders spanning 4 SEEDS** (the seed IS the shard, `J1039`) and compare every
   design against it. That is exactly the cross-seed comparison `s117` paid 12 extra folds to avoid,
   and `J986` §3 rule 5 forbids: the floor MOVES with the batch.
   ⇒ `s60` is invoked **once per shard root with that shard's own manifest**, so each shard's 20
   designs are floored by the 2 non-binders folded beside them at their own seed.

⭐⭐ `s60` IS NOT EDITED AND NOT RE-IMPLEMENTED. It is loaded FRESH per shard via `importlib` with
   `sys.argv` set for that shard, so each load executes its module body with its own ROOT/MAN/OUT.
   Every metric, the floor rule and the conjunction are `s60`'s own objects. ⛔ `s60` must stay
   byte-identical on disk because `J997`/`ADR 0063` were adjudicated on it. This file adds ONLY
   aggregation across the four per-shard results.

⛔⛔ A PROVENANCE DEFECT IN `s60`'s OUTPUT THAT THIS FILE MUST NOT PROPAGATE. `s60`:107-108 HARD-CODES
   into its report header *"the ADOPTED r4 bias cohort (A:-0.495,E:-0.354,L:-0.092)"* and
   *"seed 20261104"*. ⛔ BOTH ARE FALSE FOR THIS ROUND: `[M]` `s53` ran this cohort with **NO
   `--bias_AA`** (the deliberate unbiased measurement), and the seeds here are **01-04**, one per
   shard. ⇒ the four per-shard `s60` reports are SCRATCH INTERMEDIATES and ⛔ must never be committed
   or cited as this round's record; this file's own output carries the true resolved config.
   (Same class of hazard as `s83`'s hard-coded banner in `t108`/`t110`, still open.)

⭐ `n` IS BATCHES ([[rfd3-designs-are-clustered-by-batch]]) -- and here that costs nothing, because
   `s115` selected **one design per diffusion batch**. ⇒ the 20 designs per arm ARE 20 independent
   contigs and design-level aggregation IS batch-level. ⛔ That is a property of this cohort's
   SELECTION, not a general licence.

⛔ THE ENDPOINT IS `J1055` §3's CONJUNCTION, PRE-REGISTERED BEFORE THE DATA EXISTED:
   peptide contacts > the in-batch non-binder floor **AND** binder<->peptide `pair(3,2)` above the
   same floor. ⛔ Never either alone, never a ratio. `pair(3,2)` is confounded by the shared anchor
   chain (`J1054` §2: Spearman +0.681 to HLA contacts vs +0.351 to peptide) and the confound only
   INFLATES ⇒ ipTM can REJECT but can never PROMOTE, so it is used here as a TIE-BREAK only and the
   ranking quantity is the CONTACT MARGIN.
⛔ The STRICT floor (the HIGHER non-binder) is the ONLY headline. The lenient floor is reported as a
   sensitivity line. ⛔ Switching to the lenient one after seeing the count is forbidden.
⭐ ZERO passing designs is a legitimate, pre-registered outcome (`J1055` §3/§7) and is reported as
   zero -- ⛔ not backfilled, not relaxed, not re-folded.

Usage:
  python3 -I s118_score_2x2_cofold_per_shard.py <extract_root> <master_manifest.tsv> <out.md>
"""
import collections
import contextlib
import csv
import importlib.util
import io
import os
import statistics as st
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S60 = os.path.join(HERE, "s60_score_r4_cofold.py")

if len(sys.argv) != 4:
    sys.exit(__doc__)
ROOT, MASTER, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

# ⛔ assert the defects this file is built around are still the ones in s60's source
_SRC = open(S60).read()
for _needle, _what in (("FLOOR_C = max(mean(r[\"cpep\"]) for r in nb)", "the global-max floor rule"),
                       ("seed 20261104", "the stale hard-coded seed in the header"),
                       ("ADOPTED r4 bias cohort", "the stale hard-coded bias cohort")):
    if _needle not in _SRC:
        sys.exit(f"⛔ REFUSING: s60's source no longer contains {_what} ({_needle!r}). This file's "
                 f"whole justification is built on it -- re-audit before use.")

master = {r["tag"]: r for r in csv.DictReader(open(MASTER), delimiter="\t")}
shards = sorted(d for d in os.listdir(ROOT) if d.startswith("boltz_results_shard"))
if not shards:
    sys.exit(f"⛔ no boltz_results_shard* under {ROOT}")

mean = lambda v: float(np.nanmean(v))


def load_s60(root, man, out):
    """Load s60 FRESH so its module body runs with this shard's argv."""
    spec = importlib.util.spec_from_file_location(f"s60_{os.path.basename(root)}", S60)
    mod = importlib.util.module_from_spec(spec)
    saved = sys.argv
    sys.argv = ["s60_score_r4_cofold.py", root, man, out]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            try:
                spec.loader.exec_module(mod)
            except SystemExit:
                pass
    finally:
        sys.argv = saved
    return mod


L = []
W = L.append
W("# s118 -- J1056 2x2 REFOLD, scored PER SHARD with each shard's own in-batch floor")
W("")
W("## ⭐ THE RESOLVED CONFIG, recorded here because `s60`'s own header is STALE")
W("- `[M]` sequence design: SolubleMPNN via `s53`, **NO `--bias_AA`** (the unbiased measurement), "
  "`--omit_AA C`, T=0.1, 8 seqs/backbone, best by `overall_confidence`")
W("- `[M]` fold: Boltz, `msa: empty` on every chain, `diffusion_samples 5`, recycling 3, "
  "**seeds 01-04, one per shard**")
W("- `[M]` chains: A=**HLA-A\\*11:01**(276) B=b2m(100) C=peptide G12D `VVVGADGVGK`(10) D=binder "
  "⇒ binder<->peptide is `pair_chains_iptm(3,2)`")
W("- ⛔⛔ **THE ALLELE LABEL IN THIS ROUND'S OWN `manifest.json` IS WRONG.** `[M]` "
  "`the pinned fold runner`:381 hard-codes "
  "*\"A=HLA-A\\*03:01\"*, and this round reused that runner **byte-for-byte** (deliberately — its "
  "sha1 is pinned and pods refuse on `FATAL_SHA1_MISMATCH`), so all four shard manifests carry the "
  "wrong allele. ⭐ The true allele is **A\\*11:01** (`J960`:1; `s35`:15-16). ⛔ The SEQUENCES are "
  "correct and unaffected — they are lifted from `s35`'s source by `ast` and `s116` asserted "
  "per design that the folded HLA/b2m/peptide chains match it exactly — so this is a FALSE LABEL in "
  "an artefact, ⛔ not a scientific error. Same class as `s83`'s hard-coded banner; ⛔ do not "
  "\"fix\" the runner, its bytes are load-bearing.")
W("- ⛔ `s60`:107-108 prints *\"ADOPTED r4 bias cohort\"* and *\"seed 20261104\"*. ⛔⛔ BOTH FALSE "
  "here. The per-shard `s60` reports are scratch intermediates and must not be cited.")
W("- ⭐ `s115` picked **one design per diffusion batch**, so the 20 designs per arm are 20 "
  "INDEPENDENT contigs ⇒ design-level == batch-level.")
W("")

per_shard, allrows = {}, []
W("## 0 · PER-SHARD FLOORS — ⛔ never pooled (the seed IS the shard)")
W("")
W(f"| shard | seed | designs | non-binders | STRICT floor cont_pep | STRICT floor D-C ipTM | "
  f"lenient cont_pep | lenient ipTM |")
W("|---|---|---|---|---|---|---|---|")
for sd in shards:
    sn = sd.replace("boltz_results_shard", "")
    root = os.path.join(ROOT, sd)
    preds = os.path.join(root, "predictions")
    tags = sorted(os.listdir(preds)) if os.path.isdir(preds) else []
    man_path = os.path.join(ROOT, f"manifest_shard{sn}.tsv")
    cols = list(next(iter(master.values())).keys())
    with open(man_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for t in tags:
            if t in master:
                fh.write("\t".join(master[t][c] for c in cols) + "\n")
    mod = load_s60(root, man_path, os.path.join(ROOT, f"s60_shard{sn}_SCRATCH.txt"))
    if not hasattr(mod, "rows") or not getattr(mod, "HAVE_FLOOR", False):
        sys.exit(f"⛔ shard{sn}: s60 produced no rows or NO FLOOR -- designs would be "
                 f"uninterpretable; refusing to aggregate")
    nb = [r for r in mod.rows if r["klass"].startswith("NONBINDER")]
    if len(nb) != 2:
        sys.exit(f"⛔ shard{sn} has {len(nb)} non-binders, expected 2 -- ⛔ a shard missing a floor "
                 f"must NEVER borrow another shard's ([[never-rank-on-partial-data]])")
    fC, fI = mod.FLOOR_C, mod.FLOOR_I
    lC, lI = mod.LEN_C, mod.LEN_I
    per_shard[sn] = dict(fC=fC, fI=fI, lC=lC, lI=lI)
    for r in mod.rows:
        if r["klass"] != "DESIGN":
            continue
        cp, dp = mean(r["cpep"]), mean(r["dpep"])
        allrows.append(dict(tag=r["tag"], arm=r["cell"], shard=sn, cpep=cp, dpep=dp,
                            chla=mean(r["chla"]), bpl=mean(r["bpl"]),
                            margin=cp - fC, imargin=dp - fI,
                            pass_c=cp > fC, pass_i=dp > fI,
                            conj=(cp > fC and dp > fI),
                            conj_len=(cp > lC and dp > lI)))
    nd = sum(1 for r in mod.rows if r["klass"] == "DESIGN")
    W(f"| {sn} | {int(sn)} | {nd} | {len(nb)} | {fC:.1f} | {fI:.3f} | {lC:.1f} | {lI:.3f} |")

if len(allrows) != 80:
    W("")
    W(f"⛔⛔ {len(allrows)} designs aggregated, expected 80 -- REFUSING to read the 2x2.")
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    sys.exit(1)

W("")
W(f"- `[M]` **{len(allrows)} designs** aggregated over {len(per_shard)} shards, each floored by its "
  f"OWN 2 non-binders.")
W(f"- ⚠️ At n=2 per shard a floor is POORLY DETERMINED. The STRICT floor (higher non-binder) is the "
  f"only headline; the lenient column is a sensitivity line and ⛔ was fixed as such before the "
  f"counts below were seen.")

# ── the endpoint ──────────────────────────────────────────────────────────────────────────────
W("")
W("## 1 · ⭐⭐ THE PRE-REGISTERED CONJUNCTION (`J1055` §3)")
W("")
W("| arm | motif | breadth | n | **conj (STRICT)** | contacts alone | ipTM alone | conj (lenient) |")
W("|---|---|---|---|---|---|---|---|")
by_arm = collections.defaultdict(list)
for r in allrows:
    by_arm[r["arm"]].append(r)
for arm in sorted(by_arm):
    g = by_arm[arm]
    mo = "held" if arm.endswith(("_motif_narrow", "_motif_broad")) else "none"
    br = "broad" if arm.endswith("broad") else "narrow"
    W(f"| `{arm}` | {mo} | {br} | {len(g)} | **{sum(r['conj'] for r in g)}** | "
      f"{sum(r['pass_c'] for r in g)} | {sum(r['pass_i'] for r in g)} | "
      f"{sum(r['conj_len'] for r in g)} |")
tot = sum(r["conj"] for r in allrows)
W("")
W(f"### ⇒ **{tot} of {len(allrows)} designs clear the pre-registered conjunction** against their "
  f"own in-batch floor.")
if tot == 0:
    W("")
    W("## ⛔⛔ THE PRE-REGISTERED FALSIFIER HAS FIRED")
    W("`J1055` §3, written before the data existed: *\"If no arm clears the conjunction above its "
      "own in-batch non-binder floor, then backbone identity is **not** the binding constraint "
      "either, and the failure is upstream of backbone generation — in the **ask itself** (the "
      "anchor set, or the premise that this peptide presents a designable epitope at all).\"*")
    W("⛔ That is a RESULT, not a reason to re-score, relax the floor, or re-fold. ⛔ It must not be "
      "re-described as anything else.")

# ── the 2x2 on the endpoint ───────────────────────────────────────────────────────────────────
W("")
W("## 2 · THE 2x2 ON THE ENDPOINT, and on the CONTACT MARGIN")
W("")
W("| | narrow | broad |")
W("|---|---|---|")
for mo in ("held", "none"):
    cells = []
    for br in ("narrow", "broad"):
        g = [r for r in allrows
             if (("held" if r["arm"].endswith(("_motif_narrow", "_motif_broad")) else "none") == mo)
             and (("broad" if r["arm"].endswith("broad") else "narrow") == br)]
        cells.append(f"{sum(r['conj'] for r in g)}/{len(g)}")
    W(f"| **motif {mo}** | {cells[0]} | {cells[1]} |")
W("")
W("| factor | level A | level B | mean contact margin A | B | Δ | mean ipTM margin A | B |")
W("|---|---|---|---|---|---|---|---|")
for fname, la, lb, key in (("motif", "held", "none", "mo"), ("breadth", "broad", "narrow", "br")):
    def lev(r):
        if key == "mo":
            return "held" if r["arm"].endswith(("_motif_narrow", "_motif_broad")) else "none"
        return "broad" if r["arm"].endswith("broad") else "narrow"
    a = [r for r in allrows if lev(r) == la]
    b = [r for r in allrows if lev(r) == lb]
    W(f"| {fname} | {la} (n={len(a)}) | {lb} (n={len(b)}) | "
      f"{st.mean(r['margin'] for r in a):+.2f} | {st.mean(r['margin'] for r in b):+.2f} | "
      f"{st.mean(r['margin'] for r in a) - st.mean(r['margin'] for r in b):+.2f} | "
      f"{st.mean(r['imargin'] for r in a):+.3f} | {st.mean(r['imargin'] for r in b):+.3f} |")

# ── the promotion list ────────────────────────────────────────────────────────────────────────
W("")
W("## 3 · ⭐ THE PROMOTION LIST — ranked on CONTACT MARGIN, ipTM as TIE-BREAK ONLY")
W("")
W("⛔ Ranking on `pair(3,2)` would rank on the confounded quantity. The confound only INFLATES, so "
  "ipTM can REJECT and never PROMOTE ⇒ it breaks ties and nothing more.")
W("⭐ Clearing the floor makes a design a **CANDIDATE, never a result** (`J1055` §3).")
W("")
ok = sorted([r for r in allrows if r["conj"]],
            key=lambda r: (-r["margin"], -r["imargin"], r["tag"]))
if not ok:
    W("⛔ **EMPTY — no design clears the conjunction.** Reported as zero; ⛔ not backfilled.")
else:
    W(f"| # | tag | arm | shard | cont_pep | floor | **margin** | D-C ipTM | floor | binder pLDDT |")
    W("|---|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(ok, 1):
        ps = per_shard[r["shard"]]
        W(f"| {i} | `{r['tag']}` | {r['arm']} | {r['shard']} | {r['cpep']:.1f} | {ps['fC']:.1f} | "
          f"**{r['margin']:+.1f}** | {r['dpep']:.3f} | {ps['fI']:.3f} | {r['bpl']:.1f} |")

TSV = OUT.replace(".md", "") + "_per_design.tsv"
with open(TSV, "w") as fh:
    cols = ["tag", "arm", "shard", "cpep", "dpep", "chla", "bpl", "margin", "imargin",
            "pass_c", "pass_i", "conj", "conj_len"]
    fh.write("\t".join(cols) + "\n")
    for r in sorted(allrows, key=lambda x: x["tag"]):
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
W("")
W(f"- wrote `{os.path.basename(TSV)}` — {len(allrows)} per-design rows")
W("- ⛔ the four `s60_shard*_SCRATCH.txt` files carry `s60`'s STALE header and are NOT the record.")

open(OUT, "w").write("\n".join(L) + "\n")
print("\n".join(L))
