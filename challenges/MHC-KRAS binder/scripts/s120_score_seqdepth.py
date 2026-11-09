#!/usr/bin/env python3
"""s120 -- score the J1057 SEQUENCE-DEPTH round, and answer the backbone-vs-sequence question.

⛔ READ-ONLY on every input. Nothing launches, spends, or edits a table.

⛔⛔ WHY NOT `s118` (which I wrote hours earlier, and which is NOT edited here).
   `s118` produced `t151` and is left byte-identical so `t151` stays reproducible. Three of its
   assumptions are false for this round, each of which would make it REFUSE or CRASH rather than
   mislead -- but refusing is still not scoring:
     ⑴ `s118` asserts `len(allrows) != 80`; this round has **160** fresh designs.
     ⑵ `s118` tests a `motif` factor with two levels; here motif is **constant (none)** -- the
        held-motif arms were deliberately not deepened -- so a `held` group of n=0 would crash on
        `st.mean` of an empty sequence.
     ⑶ `s118` ignores `klass == "DESIGN_REDRAW"`, which is exactly the 2 survivors' 20-draw
        re-fold ⇒ the headline robustness check would be silently dropped.
   ⇒ `s60` is still the metric and the floor rule, loaded FRESH PER SHARD exactly as `s118` does it
   (⛔ `s60`:161,177-178 take ONE GLOBAL floor over whatever root they are handed, and the seed IS
   the shard). The shared loader is asserted against `s118`'s source at run time.

⭐⭐ THE QUESTION, and the unit it must be answered on. `t151` folded **MPNN rank 1 only**, so
   `[M]` 2 of 40 no-motif BACKBONES produced a clearing sequence. This round adds **ranks 2-5** on
   those same 40 backbones. ⇒ the test is at the **BACKBONE** level:
       does a backbone that failed on rank 1 clear on some other sequence?
   ⛔ Comparing design COUNTS would be meaningless (20/arm at rank 1 against 80/arm here); the
   comparable quantities are the per-arm RATE and the per-BACKBONE union.
   ⭐ Each design is scored against ITS OWN round's in-batch floor, so the union across rounds is
   legitimate: every member was floored by non-binders folded beside it.

⛔ ENDPOINT UNCHANGED, `J1055` §3: peptide contacts > the in-batch non-binder floor AND
   binder<->peptide `pair(3,2)` above the same floor. ⛔ Never either alone, never a ratio. ipTM is
   confounded by the shared anchor chain and the confound only INFLATES ⇒ it can REJECT, never
   PROMOTE ⇒ ranking is on CONTACT MARGIN, ipTM as tie-break only. STRICT floor (the higher
   non-binder) is the only headline; lenient is a sensitivity line, fixed as such in advance.
⭐ ZERO is a legitimate pre-registered outcome and is reported as zero.

Usage:
  python3 -I s120_score_seqdepth.py <extract_root> <manifest.tsv> <t151_per_design.tsv> <out.md>
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

if len(sys.argv) != 5:
    sys.exit(__doc__)
ROOT, MASTER, T151, OUT = sys.argv[1:5]

_S60 = open(S60).read()
if 'FLOOR_C = max(mean(r["cpep"]) for r in nb)' not in _S60:
    sys.exit("⛔ REFUSING: s60's global-max floor rule has changed -- re-audit before use.")
if "s51.roles = roles_scaffolded" in _S60:
    sys.exit("⛔ REFUSING: s60 appears modified -- it must stay byte-identical (ADR 0063).")
if 'sys.argv = ["s60_score_r4_cofold.py", root, man, out]' not in open(S118).read():
    sys.exit("⛔ REFUSING: s118's per-shard loader has changed -- this file copies it and the drift "
             "would be silent.")

master = {r["tag"]: r for r in csv.DictReader(open(MASTER), delimiter="\t")}
prev = {r["tag"]: r for r in csv.DictReader(open(T151), delimiter="\t")}
shards = sorted(d for d in os.listdir(ROOT) if d.startswith("boltz_results_shard"))
if not shards:
    sys.exit(f"⛔ no boltz_results_shard* under {ROOT}")

mean = lambda v: float(np.nanmean(v))


def load_s60(root, man, out):
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
W("# s120 -- J1057 SEQUENCE DEPTH: does another sequence rescue a failed backbone?")
W("")
W("- `[M]` fold: Boltz, `msa: empty` every chain, `diffusion_samples 5`, recycling 3, seeds 01-04")
W("- `[M]` MPNN: SolubleMPNN, **NO `--bias_AA`**, `--omit_AA C`, T=0.1 — ranks **2-5** per backbone")
W("- ⛔ `s60`'s report header hard-codes *\"r4 bias cohort\"* / *\"seed 20261104\"*, both FALSE here "
  "⇒ the per-shard `s60` files are SCRATCH, never cite them")
W("- ⛔ the shard `manifest.json` carries `A=HLA-A*03:01` from the pinned runner; the true allele is "
  "**A\\*11:01** and the SEQUENCES are correct (asserted per design by `s119`)")
W("- ⛔ SCOPE: **no-motif arms only.** The held-motif arms were refuted on two axes and deliberately "
  "not deepened ⇒ ⛔ nothing here speaks to them.")
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
                   rank=m["mpnn_rank"], cpep=cp, dpep=dp, chla=mean(r["chla"]),
                   bpl=mean(r["bpl"]), margin=cp - fC, imargin=dp - fI,
                   pass_c=cp > fC, pass_i=dp > fI, conj=(cp > fC and dp > fI),
                   conj_len=(cp > lC and dp > lI))
        (redraw if r["klass"] == "DESIGN_REDRAW" else fresh).append(rec)
    W(f"| {sn} | {sum(1 for r in mod.rows)} | {len(nb)} | {fC:.1f} | {fI:.3f} | {lC:.1f} | {lI:.3f} |")

W("")
W(f"- `[M]` **{len(fresh)} fresh designs** (ranks 2-5) and **{len(redraw)} re-draw folds** "
  f"aggregated over {len(per_shard)} shards.")
if len(fresh) != 160 or len(redraw) != 8:
    W(f"⛔⛔ expected 160 fresh + 8 re-draw folds; got {len(fresh)} + {len(redraw)}. REFUSING.")
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    sys.exit(1)

# ── 1 · the fresh designs ─────────────────────────────────────────────────────────────────────
W("")
W("## 1 · RANKS 2-5 ON THE PRE-REGISTERED CONJUNCTION")
W("")
W("| arm | n | **conj (STRICT)** | rate | contacts alone | ipTM alone | conj (lenient) |")
W("|---|---|---|---|---|---|---|")
by_arm = collections.defaultdict(list)
for r in fresh:
    by_arm[r["arm"]].append(r)
for arm in sorted(by_arm):
    g = by_arm[arm]
    W(f"| `{arm}` | {len(g)} | **{sum(r['conj'] for r in g)}** | "
      f"{100.0 * sum(r['conj'] for r in g) / len(g):.1f}% | {sum(r['pass_c'] for r in g)} | "
      f"{sum(r['pass_i'] for r in g)} | {sum(r['conj_len'] for r in g)} |")
W("")
W(f"### ⇒ **{sum(r['conj'] for r in fresh)} of {len(fresh)}** ranks-2-5 designs clear the "
  f"conjunction against their own in-batch floor.")
W("")
W("| factor | level A | level B | mean contact margin A | B | Δ |")
W("|---|---|---|---|---|---|")
a = [r for r in fresh if r["arm"].endswith("broad")]
b = [r for r in fresh if not r["arm"].endswith("broad")]
W(f"| breadth | broad (n={len(a)}) | narrow (n={len(b)}) | {st.mean(r['margin'] for r in a):+.2f} | "
  f"{st.mean(r['margin'] for r in b):+.2f} | "
  f"{st.mean(r['margin'] for r in a) - st.mean(r['margin'] for r in b):+.2f} |")
W("⛔ `motif` is CONSTANT (none) in this round and is therefore not a factor here.")

# ── 2 · THE BACKBONE-LEVEL ANSWER ─────────────────────────────────────────────────────────────
W("")
W("## 2 · ⭐⭐ THE ANSWER: BACKBONE vs SEQUENCE, at the BACKBONE level")
W("")
r1 = {t: r for t, r in prev.items() if r["conj"] in ("True", "False")}
# map t151 tags back to their design name the way s116 built them
r1_by_design = {}
for t, r in r1.items():
    r1_by_design.setdefault(r["arm"], None)
by_design = collections.defaultdict(list)
for r in fresh:
    by_design[r["design"]].append(r)
# rank-1 outcome per design, from t151's own tag convention D_<design with arm_shard->s>
rank1 = {}
for t, r in r1.items():
    rank1[t] = r["conj"] == "True"
n_bb = len(by_design)
bb_r1_pass, bb_any_pass, rescued = 0, 0, []
for design, rows in sorted(by_design.items()):
    t = "D_" + design.replace("arm_shard", "s")
    was = rank1.get(t)
    now = any(r["conj"] for r in rows)
    if was:
        bb_r1_pass += 1
    if was or now:
        bb_any_pass += 1
    if now and not was:
        best = max((r for r in rows if r["conj"]), key=lambda r: r["margin"])
        rescued.append((design, best))
W(f"- `[M]` backbones deepened: **{n_bb}** (the 40 no-motif backbones, 4 sequences each)")
W(f"- `[M]` cleared on **rank 1 alone** (`t151`): **{bb_r1_pass}/{n_bb}**")
W(f"- `[M]` cleared on **any of ranks 1-5**: **{bb_any_pass}/{n_bb}**")
W(f"- ⭐⭐ `[M]` **backbones RESCUED by a non-rank-1 sequence: {len(rescued)}**")
W("")
if not rescued:
    W("### ⇒ ⛔⛔ **SEQUENCE DEPTH BUYS NOTHING HERE.** Four extra sequences per backbone rescued "
      "**ZERO** backbones that rank 1 had failed.")
    W("⇒ ⭐ the constraint is the **BACKBONE**, not the sequence — which is the same direction the "
      "design-stage 2×2 and the refold already pointed, and it is now measured on its own axis.")
    W("⛔ It does NOT say ranks 2-5 are worse than rank 1; it says they are not DIFFERENT enough to "
      "cross an in-batch floor this high.")
else:
    W("| design | arm | shard | rank | cont_pep | floor | **margin** | D-C ipTM | floor |")
    W("|---|---|---|---|---|---|---|---|---|")
    for design, r in sorted(rescued, key=lambda x: -x[1]["margin"]):
        ps = per_shard[r["shard"]]
        W(f"| `{design}` | {r['arm']} | {r['shard']} | {r['rank']} | {r['cpep']:.1f} | "
          f"{ps['fC']:.1f} | **{r['margin']:+.1f}** | {r['dpep']:.3f} | {ps['fI']:.3f} |")
    W("")
    W(f"### ⇒ ⭐⭐ **{len(rescued)} backbone(s) that failed on rank 1 clear on another sequence** ⇒ "
      f"the sequence matters and rank-1-only folding UNDER-COUNTS the cohort.")

# ── 3 · the re-draw, i.e. are the 2 survivors real ────────────────────────────────────────────
W("")
W("## 3 · ⭐⭐ THE RE-DRAW — are `t151`'s 2 survivors robust at 20 draws?")
W("")
W("⚠️ `t151` scored them on **5 draws**. `[M]` re-fold noise runs ~1.56× the within-batch s.e. ⇒ "
  "this is the check that decides whether they survive at all.")
W("")
W("| candidate | folds | draws | mean cont_pep | sd | mean floor | **mean margin** | mean ipTM | "
  "clears per shard |")
W("|---|---|---|---|---|---|---|---|---|")
by_cand = collections.defaultdict(list)
for r in redraw:
    by_cand[r["design"]].append(r)
for design, rows in sorted(by_cand.items()):
    cps = [r["cpep"] for r in rows]
    fCs = [per_shard[r["shard"]]["fC"] for r in rows]
    mg = [r["margin"] for r in rows]
    W(f"| `{design}` | {len(rows)} | {len(rows) * 5} | {st.mean(cps):.1f} | "
      f"{(st.stdev(cps) if len(cps) > 1 else 0.0):.2f} | {st.mean(fCs):.1f} | "
      f"**{st.mean(mg):+.1f}** | {st.mean(r['dpep'] for r in rows):.3f} | "
      f"{sum(r['conj'] for r in rows)}/{len(rows)} |")
W("")
for design, rows in sorted(by_cand.items()):
    k = sum(r["conj"] for r in rows)
    verdict = ("⭐⭐ **HOLDS** at every seed" if k == len(rows)
               else ("⚠️ **PARTIAL** — clears at some seeds only ⇒ ⛔ inside noise, do not promote "
                     "it as a result" if k else "⛔⛔ **DOES NOT SURVIVE** the re-draw"))
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
