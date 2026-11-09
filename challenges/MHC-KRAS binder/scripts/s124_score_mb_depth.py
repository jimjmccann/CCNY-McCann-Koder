#!/usr/bin/env python3
"""s124 -- score the J1058 MB-DEPTH round, and REPLICATE the rank-1-under-count result within one cohort.

⛔ READ-ONLY on every input. Nothing launches, spends, or edits a table.

⛔⛔ WHY NOT `s120` (which is NOT edited here). `s120` produced `t152` and is left byte-identical so
   `t152` stays reproducible. Two of its assumptions are false for this round and would make it
   REFUSE rather than mislead -- but refusing is still not scoring:
     ⑴ `s120`:148 asserts exactly **160 fresh + 8 re-draw** folds; this round's fresh count is
        `backbones x 5` and is passed in explicitly instead of hard-coded.
     ⑵ `s120`:174-178 splits a **breadth** factor on `arm.endswith("broad")`; here the arm is
        **constant `MB_nomotif_broad`**, so the `narrow` group is n=0 and `st.mean` of an empty
        sequence raises.
   ⇒ `s60` is still the metric and the floor rule, loaded FRESH PER SHARD exactly as `s118`/`s120`
   do it (⛔ `s60`:161,177-178 take ONE GLOBAL floor over whatever root they are handed, and
   **the seed IS the shard**). Both shared loaders are asserted against their sources at run time.

⭐⭐ WHAT THIS ROUND CAN DO THAT `t152` COULD NOT, AND IT IS THE REASON TO SCORE IT THIS WAY.
   `t152` compared **ranks 2-5 here** against **rank 1 in a DIFFERENT round** (`t151`) -- legitimate,
   because each design was floored by non-binders folded beside it, but still a cross-round union.
   This round folds **ranks 1-5 together on the same backbones in the same batch**, so the
   under-count can be measured **WITHIN ONE COHORT**:
       backbones clearing on **rank 1 alone**   vs   backbones clearing on **any of ranks 1-5**
   ⇒ an independent replication of `t152`'s headline (`[M]` 2/40 -> 4/40, i.e. 2x) with the
   cross-round step removed. ⛔ It is a REPLICATION, not the same measurement: a different cohort
   can legitimately give a different ratio, and ⛔ one cohort agreeing does not make either a law.

⛔ ENDPOINT UNCHANGED, `J1055` §3: peptide contacts > the in-batch non-binder floor AND
   binder<->peptide `pair(3,2)` above the same floor. ⛔ Never either alone, never a ratio. ipTM is
   confounded by the shared anchor chain and the confound only INFLATES ⇒ it can REJECT, never
   PROMOTE ⇒ ranking is on CONTACT MARGIN, ipTM as tie-break only. STRICT floor (the higher
   non-binder) is the only headline; lenient is a sensitivity line, fixed as such in advance.
⭐ ZERO is a legitimate pre-registered outcome and is reported as zero.
⛔ Clearing a floor makes a design a CANDIDATE, never a result.

Usage:
  python3 -I s124_score_mb_depth.py <extract_root> <manifest.tsv> <out.md>
          <expect_fresh> <expect_redraw>
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
S118 = os.path.join(HERE, "s118_score_2x2_cofold_per_shard.py")
S120 = os.path.join(HERE, "s120_score_seqdepth.py")

if len(sys.argv) != 6:
    sys.exit(__doc__)
ROOT, MASTER, OUT = sys.argv[1:4]
EXP_FRESH, EXP_REDRAW = int(sys.argv[4]), int(sys.argv[5])

# ⛔ the metric and the floor rule must be s60's own, unmodified (ADR 0063)
_S60 = open(S60).read()
if 'FLOOR_C = max(mean(r["cpep"]) for r in nb)' not in _S60:
    sys.exit("⛔ REFUSING: s60's global-max floor rule has changed -- re-audit before use.")
if "s51.roles = roles_scaffolded" in _S60:
    sys.exit("⛔ REFUSING: s60 appears modified -- it must stay byte-identical (ADR 0063).")
if 'sys.argv = ["s60_score_r4_cofold.py", root, man, out]' not in open(S118).read():
    sys.exit("⛔ REFUSING: s118's per-shard loader has changed -- this file copies it and the drift "
             "would be silent.")
# ⛔⛔ assert the premises this file documents itself as diverging from. If s120 changes, the ⑴/⑵
#    reasoning above is stale and the divergence is no longer deliberate.
_S120 = open(S120).read()
for _n in ('if len(fresh) != 160 or len(redraw) != 8:', 'r["arm"].endswith("broad")'):
    if _n not in _S120:
        sys.exit(f"⛔ REFUSING: `{_n}` is absent from s120's source -- this file states its "
                 f"divergence from s120 in terms of it. Re-audit.")

master = {r["tag"]: r for r in csv.DictReader(open(MASTER), delimiter="\t")}
shards = sorted(d for d in os.listdir(ROOT) if d.startswith("boltz_results_shard"))
if not shards:
    sys.exit(f"⛔ no boltz_results_shard* under {ROOT}")

mean = lambda v: float(np.nanmean(v))


def load_s60(root, man, out):
    """s118's per-shard loader, copied and asserted above. ⛔ s60 must see ONE shard's root only."""
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
W("# s124 -- J1058 MB DEPTH: new MB backbones, ranks 1-5, and a within-cohort replication")
W("")
W("- `[M]` fold: Boltz, `msa: empty` every chain, `diffusion_samples 5`, recycling 3, seeds 01-04")
W("- `[M]` MPNN: SolubleMPNN, **NO `--bias_AA`**, `--omit_AA C`, T=0.1 — ranks **1-5** per backbone")
W("- `[M]` design: RFd3, 2 shards of `MB_nomotif_broad`, spec byte-identical to the 2×2's arm 03 "
  "(`spec_sha1 23d9c25c…`), `RFD3_NBATCH=20 × diffusion_batch_size=4`")
W("- ⛔ `s60`'s report header hard-codes *\"r4 bias cohort\"* / *\"seed 20261104\"*, both FALSE here "
  "⇒ the per-shard `s60` files are SCRATCH, never cite them")
W("- ⛔ the shard `manifest.json` carries `A=HLA-A*03:01` from the pinned runner; the true allele is "
  "**A\\*11:01** and the SEQUENCES are correct (asserted per design by `s123`)")
W("- ⛔ SCOPE: **ONE cell (`MB_nomotif_broad`).** There is no breadth or motif contrast in this "
  "round by construction ⇒ ⛔ nothing here speaks to either.")
W("")

per_shard, fresh, redraw = {}, [], []
W("## 0 · PER-SHARD FLOORS — ⛔ never pooled (the seed IS the shard)")
W("")
W("| shard | folds | non-binders | STRICT cont_pep | STRICT D-C ipTM | lenient cont_pep | lenient ipTM |")
W("|---|---|---|---|---|---|---|")
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
        sys.exit(f"⛔ shard{sn}: no rows or NO FLOOR -- designs uninterpretable, refusing")
    nb = [r for r in mod.rows if r["klass"].startswith("NONBINDER")]
    if len(nb) != 2:
        sys.exit(f"⛔ shard{sn} has {len(nb)} non-binders, expected 2 -- ⛔ a shard must never "
                 f"borrow another's floor")
    fC, fI, lC, lI = mod.FLOOR_C, mod.FLOOR_I, mod.LEN_C, mod.LEN_I
    per_shard[sn] = dict(fC=fC, fI=fI, lC=lC, lI=lI)
    for r in mod.rows:
        if r["klass"] not in ("DESIGN", "DESIGN_REDRAW"):
            continue
        cp, dp = mean(r["cpep"]), mean(r["dpep"])
        m = master[r["tag"]]
        rec = dict(tag=r["tag"], arm=r["cell"], shard=sn, design=m["design"],
                   rank=int(m["mpnn_rank"]), cpep=cp, dpep=dp, chla=mean(r["chla"]),
                   bpl=mean(r["bpl"]), margin=cp - fC, imargin=dp - fI,
                   pass_c=cp > fC, pass_i=dp > fI, conj=(cp > fC and dp > fI),
                   conj_len=(cp > lC and dp > lI))
        (redraw if r["klass"] == "DESIGN_REDRAW" else fresh).append(rec)
    W(f"| {sn} | {sum(1 for r in mod.rows)} | {len(nb)} | {fC:.1f} | {fI:.3f} | {lC:.1f} | {lI:.3f} |")

W("")
W(f"- `[M]` **{len(fresh)} fresh designs** (ranks 1-5) and **{len(redraw)} re-draw folds** "
  f"aggregated over {len(per_shard)} shards.")
if len(fresh) != EXP_FRESH or len(redraw) != EXP_REDRAW:
    W(f"⛔⛔ expected {EXP_FRESH} fresh + {EXP_REDRAW} re-draw folds; got {len(fresh)} + "
      f"{len(redraw)}. REFUSING — a scorer on a short census reports a clean number for an "
      f"incomplete round.")
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    sys.exit(1)

# ── 1 · the cohort on the pre-registered conjunction ──────────────────────────────────────────
W("")
W("## 1 · RANKS 1-5 ON THE PRE-REGISTERED CONJUNCTION")
W("")
W("| rank | n | **conj (STRICT)** | rate | contacts alone | ipTM alone | mean margin |")
W("|---|---|---|---|---|---|---|")
by_rank = collections.defaultdict(list)
for r in fresh:
    by_rank[r["rank"]].append(r)
for rk in sorted(by_rank):
    g = by_rank[rk]
    W(f"| {rk} | {len(g)} | **{sum(r['conj'] for r in g)}** | "
      f"{100.0 * sum(r['conj'] for r in g) / len(g):.1f}% | {sum(r['pass_c'] for r in g)} | "
      f"{sum(r['pass_i'] for r in g)} | {st.mean(r['margin'] for r in g):+.2f} |")
W("")
W(f"### ⇒ **{sum(r['conj'] for r in fresh)} of {len(fresh)}** designs clear the conjunction against "
  f"their own in-batch floor.")
W("⛔ Per-rank rates are DESCRIPTIVE. The ranks are not independent draws — they are 5 sequences on "
  "the same backbone — so ⛔ do not read a rank-vs-rank difference as a result.")

# ── 2 · THE WITHIN-COHORT REPLICATION ─────────────────────────────────────────────────────────
W("")
W("## 2 · ⭐⭐ WITHIN-COHORT: DOES RANK-1-ONLY FOLDING UNDER-COUNT?")
W("")
by_design = collections.defaultdict(list)
for r in fresh:
    by_design[r["design"]].append(r)
n_bb = len(by_design)
bb_r1, bb_any, rescued = 0, 0, []
for design, rows in sorted(by_design.items()):
    r1 = [r for r in rows if r["rank"] == 1]
    was = any(r["conj"] for r in r1)
    now = any(r["conj"] for r in rows)
    if was:
        bb_r1 += 1
    if now:
        bb_any += 1
    if now and not was:
        best = max((r for r in rows if r["conj"]), key=lambda r: r["margin"])
        rescued.append((design, best))
W(f"- `[M]` new MB backbones folded: **{n_bb}**, at 5 sequences each")
W(f"- `[M]` clear on **rank 1 alone**: **{bb_r1}/{n_bb}**")
W(f"- `[M]` clear on **any of ranks 1-5**: **{bb_any}/{n_bb}**")
W(f"- ⭐⭐ `[M]` **backbones found ONLY by a non-rank-1 sequence: {len(rescued)}**")
if bb_r1:
    W(f"- `[R]` under-count factor **{bb_any / bb_r1:.2f}×** "
      f"(`t152` measured **2.00×** across rounds)")
else:
    W("- ⛔ rank 1 cleared **ZERO** backbones, so no ratio is defined — ⛔ do not report a factor; "
      "the union count stands on its own.")
W("")
if not rescued:
    W("### ⇒ ⛔⛔ **NO WITHIN-COHORT UNDER-COUNT.** Every backbone that cleared, cleared on rank 1 "
      "too ⇒ `t152`'s 2× does **NOT** replicate here.")
    W("⛔ That does not retract `t152` — it was a real measurement on its own cohort — but it does "
      "mean the 2× is ⛔ **NOT a planning constant**, and folding ranks 2-5 cannot be justified on "
      "this round's evidence.")
else:
    W("| design | shard | rank | cont_pep | floor | **margin** | D-C ipTM | floor |")
    W("|---|---|---|---|---|---|---|---|")
    for design, r in sorted(rescued, key=lambda x: -x[1]["margin"]):
        ps = per_shard[r["shard"]]
        W(f"| `{design}` | {r['shard']} | {r['rank']} | {r['cpep']:.1f} | {ps['fC']:.1f} | "
          f"**{r['margin']:+.1f}** | {r['dpep']:.3f} | {ps['fI']:.3f} |")
    W("")
    W(f"### ⇒ ⭐⭐ **{len(rescued)} backbone(s) are found ONLY by a non-rank-1 sequence** ⇒ "
      f"`t152`'s result REPLICATES within a single cohort, with the cross-round step removed.")

# ── 3 · the re-draw of t152's two 5-draw designs ──────────────────────────────────────────────
W("")
W("## 3 · ⭐⭐ THE RE-DRAW — are `t152`'s two RESCUED designs robust at 20 draws?")
W("")
W("⚠️ `t152` scored them on **5 draws each**, i.e. exactly the position "
  "`..._4_model_0` was in before its own re-draw demoted it to 3/4 seeds. `[M]` re-fold noise runs "
  "~1.56× the within-batch s.e. ⭐ Each is re-folded at **its OWN clearing rank**, ⛔ not rank 1 — "
  "folding rank 1 would re-draw a different molecule.")
W("")
if not redraw:
    W("⛔ no re-draw rows present -- refusing to call the candidates validated.")
else:
    W("| candidate | rank | folds | draws | mean cont_pep | sd | mean floor | **mean margin** | "
      "mean ipTM | clears |")
    W("|---|---|---|---|---|---|---|---|---|---|")
    by_cand = collections.defaultdict(list)
    for r in redraw:
        by_cand[r["design"]].append(r)
    for design, rows in sorted(by_cand.items()):
        cps = [r["cpep"] for r in rows]
        fCs = [per_shard[r["shard"]]["fC"] for r in rows]
        mg = [r["margin"] for r in rows]
        W(f"| `{design}` | {rows[0]['rank']} | {len(rows)} | {len(rows) * 5} | {st.mean(cps):.1f} | "
          f"{(st.stdev(cps) if len(cps) > 1 else 0.0):.2f} | {st.mean(fCs):.1f} | "
          f"**{st.mean(mg):+.1f}** | {st.mean(r['dpep'] for r in rows):.3f} | "
          f"{sum(r['conj'] for r in rows)}/{len(rows)} |")
    W("")
    for design, rows in sorted(by_cand.items()):
        k = sum(r["conj"] for r in rows)
        verdict = ("⭐⭐ **HOLDS** at every seed ⇒ validated at 20 draws" if k == len(rows)
                   else ("⚠️ **PARTIAL** — clears at some seeds only ⇒ ⛔ inside noise, do NOT "
                         "promote it" if k else "⛔⛔ **DOES NOT SURVIVE** the re-draw"))
        W(f"- `{design}`: clears **{k}/{len(rows)}** seeds ⇒ {verdict}")
W("⛔ Clearing a floor makes a design a CANDIDATE, never a result (`J1055` §3).")

TSV = OUT.replace(".md", "") + "_per_design.tsv"
with open(TSV, "w") as fh:
    cols = ["tag", "design", "arm", "shard", "rank", "cpep", "dpep", "chla", "bpl", "margin",
            "imargin", "pass_c", "pass_i", "conj", "conj_len"]
    fh.write("\t".join(cols) + "\n")
    for r in sorted(fresh + redraw, key=lambda x: x["tag"]):
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
W("")
W(f"- wrote `{os.path.basename(TSV)}` — {len(fresh) + len(redraw)} rows")

open(OUT, "w").write("\n".join(L) + "\n")
print("\n".join(L))
