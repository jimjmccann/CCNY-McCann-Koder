#!/usr/bin/env python3
"""j521 -- score Boltz folds of the C3 (homotrimer binder) designs. Local, no spend.

WHY THIS EXISTS, AND WHY IT IS NOT `j465`. `j465_boltz_selfconsistency_20261017.py` reads
**only the GLOBAL `iptm`** (`:158`, `cj["iptm"]`) and its docstring's `chains_ptm` promise
is NOT implemented. On C3 that confound is WORSE than on C1: 3 binder copies + 3 target
protomers = **6 chains, 15 pairs, and only 9 of them are binder<->target**. The other 6 are
3 binder<->binder (our own designed trimer interface) and 3 target<->target (TNFa's own
crystal trimer, which Boltz rebuilds near-perfectly off a deep MSA). A global `iptm` here is
mostly a measurement of TNFa folding.
-- that confound already forced the retraction of a PUSHED claim once.

WHY NOT JUST FIX `j465`. ⛔⛔ `j497` imports `j465`'s math **under a SHA1 pin** and exits
FATAL on drift. Editing `j465` in place would break the C1 scorer in the middle of a
deadline. ⇒ this is a NEW script that imports the SAME pinned math, so C1 and C3 numbers
stay comparable by construction.

⛔⛔ THE C3 AGGREGATION RULE. A C3 binder grips **at most 2 of 3** protomers, so the mean over copies pays partial
credit for a design where one copy does all the work. Every per-copy quantity is therefore
reported THREE ways -- `_mean`, `_best`, and `_worst` -- and ⭐ **`_worst` IS THE
CONJUNCTION**: it is the number to rank on when the claim is "all three copies engage".

⛔ LAYOUT IS AN ARGUMENT, NEVER A DEFAULT.
and. The resolved layout is ECHOED before any row is written.

2026-10-18, read from the files, never assumed:
  design .cif.gz : binder A,B,C (160 CA each) + target a,b,c (152 each)
  c3_apo    pred : A           (one protomer alone)
  c3_trimer pred : A,B,C       (the designed homotrimer, NO target -- never tested before)
  c3_holo   pred : A,B,C + D,E,F   => design map D->a, E->b, F->c
  chains_ptm / pair_chains_iptm are keyed "0".."n-1" by INPUT ORDER.

⭐⭐ VALIDATION BEFORE TRUST. This script is general over chain lists, so running it as
*1 binder + 3 target* on `c1_holo` must reproduce `j497` EXACTLY. `--selftest-vs-j497 <tsv>`
does that comparison and is the only evidence that the per-pair indexing is right.

⛔ `--pred-dir` is the `predictions/` directory, i.e. `out_<arm>/boltz_results_<arm>/predictions`.
   Pointed one level up it returns n=1 / OK=0 MISSING_OUTPUT, which is NOT an empty result.
"""
import argparse
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
J497 = HERE / "j497_c1_boltz_selfconsistency_20261018.py"
J497_SHA1 = "760d8f63919d8e1aeace12492f71c94d1f00f58b"


def load_j497(allow_drift=False):
    """Import j497's readers + geometry under a SHA1 pin. No copy, no drift."""
    if not J497.exists():
        sys.exit(f"FATAL: j497 not found at {J497} -- the readers are imported, not copied.")
    got = hashlib.sha1(J497.read_bytes()).hexdigest()
    if got != J497_SHA1:
        msg = (f"FATAL: j497 SHA1 drifted.\n  expected {J497_SHA1}\n  got      {got}\n"
               "  C1 and C3 scores are compared to each other, so the shared readers are\n"
               "  pinned. Re-verify j497, then update J497_SHA1 here deliberately.")
        if not allow_drift:
            sys.exit(msg)
        print("⚠️ " + msg.replace("FATAL", "WARNING") + "\n  --allow-sha-drift given: continuing.")
    spec = importlib.util.spec_from_file_location("j497_mod", J497)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)      # safe: j497 guards main() behind __main__
    return mod, got


def agg(vals, prefix, out, nd=4):
    """⭐ mean / best / WORST. `_worst` is the conjunction -- rank on it, not the mean."""
    v = [x for x in vals if x is not None]
    if not v:
        return
    out[f"{prefix}_mean"] = round(float(np.mean(v)), nd)
    out[f"{prefix}_best"] = round(float(np.max(v)), nd)
    out[f"{prefix}_worst"] = round(float(np.min(v)), nd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", required=True,
                    help="the predictions/ dir, NOT out_<arm>")
    ap.add_argument("--design-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", required=True, choices=["apo", "trimer", "holo"])
    ap.add_argument("--binder-chains", default="A,B,C")
    ap.add_argument("--target-chains", default="D,E,F")
    ap.add_argument("--design-binder-chains", default=None)
    ap.add_argument("--design-target-chains", default="a,b,c")
    ap.add_argument("--pin-resnums", default="",
                    help="comma list; restricts the bridge to the PINNED Arg (⛔ a free Arg "
                         "elsewhere is not the designed mechanism). ⛔ ONE list for the whole "
                         "cohort is almost always WRONG -- prefer --plan-json.")
    ap.add_argument("--plan-json", default="",
                    help="⭐ the LigandMPNN plan.json. Pins are read PER BACKBONE from its `gp` "
                         "field. ⛔⛔ REQUIRED in practice for C3: measured 2026-10-18, the 12 "
                         "v50 backbones have 12 DISTINCT pin sets and lengths 159-167, so a "
                         "single hardcoded list is wrong for 11 of them. "
                         "")
    ap.add_argument("--allow-sha-drift", action="store_true")
    ap.add_argument("--selftest-vs-j497", default="",
                    help="a j497 TSV; compare shared columns row-for-row and exit")
    a = ap.parse_args()

    j, sha = load_j497(a.allow_sha_drift)
    # ⛔ kabsch_rmsd/tm_score live in j465 and j497 loads them INSIDE main(), so they are
    # not attributes of the imported module. Load them the same pinned way j497 does.
    kabsch_rmsd, tm_score, j465_sha = j.load_j465_math(a.allow_sha_drift)

    split = lambda s: [c for c in s.split(",") if c]
    bchains = split(a.binder_chains)
    tchains = split(a.target_chains) if a.mode == "holo" else []
    dbchains = split(a.design_binder_chains) if a.design_binder_chains else list(bchains)
    dtchains = split(a.design_target_chains) if a.mode == "holo" else []
    if a.mode == "apo":
        bchains = bchains[:1]
        dbchains = dbchains[:1]
    assert len(dbchains) == len(bchains), (dbchains, bchains)
    assert len(dtchains) == len(tchains), (dtchains, tchains)
    pins = {int(x) for x in split(a.pin_resnums)} or None

    # ⭐ per-backbone pins. A cohort-wide list cannot be right when every backbone differs.
    plan_pins = {}
    if a.plan_json:
        for e in json.load(open(a.plan_json)):
            if e.get("gp"):
                plan_pins[e["base"]] = {int(x) for x in e["gp"]}
        n_distinct = len({frozenset(v) for v in plan_pins.values()})
        print(f"  plan.json: pins for {len(plan_pins)} backbones, {n_distinct} DISTINCT sets")
        if pins:
            sys.exit("⛔ pass --plan-json OR --pin-resnums, not both -- they would disagree.")

    pred_chains = bchains + tchains
    design_chains = dbchains + dtchains

    # ⭐ ECHO THE RESOLVED LAYOUT. Never inherit a default silently.
    print(f"RESOLVED LAYOUT  mode={a.mode}")
    print(f"  prediction: binder={','.join(bchains)}"
          f"  target={','.join(tchains) if tchains else '(none)'}")
    print(f"  design:     binder={','.join(dbchains)}"
          f"  target={','.join(dtchains) if dtchains else '(none)'}")
    print(f"  copies={len(bchains)}  pairs={len(pred_chains)*(len(pred_chains)-1)//2}"
          f"  binder<->target pairs={len(bchains)*len(tchains)}")
    print(f"  j497 readers imported, sha1={sha}")
    print(f"  pinned Arg resnums: {sorted(pins) if pins else '(ALL Arg -- not the pins)'}")
    print("  ⭐ _worst IS THE CONJUNCTION. Rank on it, not on _mean.")

    design_dir = Path(a.design_dir)
    dcache = {}

    def find_design(tag):
        if tag in dcache:
            return dcache[tag]
        hits = []
        for ext in ("cif.gz", "cif", "pdb.gz", "pdb"):
            hits += sorted(design_dir.rglob(f"{tag}.{ext}"))
        if not hits:
            hits = sorted(p for p in design_dir.rglob(f"*{tag}*")
                          if p.is_file() and p.suffix in (".pdb", ".cif", ".gz"))
        dcache[tag] = hits[0] if hits else None
        return dcache[tag]

    rows, notes = [], []
    for sub in sorted(Path(a.pred_dir).iterdir()):
        if not sub.is_dir():
            continue
        name = sub.name
        pdbs = sorted(sub.glob("*_model_0.pdb"))
        conf = sorted(sub.glob("confidence_*_model_0.json"))
        plddt_f = sorted(sub.glob("plddt_*_model_0.npz"))
        pae_f = sorted(sub.glob("pae_*_model_0.npz"))
        if not (pdbs and conf and plddt_f):
            rows.append({"design": name, "status": "MISSING_OUTPUT"})
            continue
        if pdbs[0].stat().st_size == 0:
            rows.append({"design": name, "status": "EMPTY_PDB"})
            continue

        tag = re.sub(r"_seq\d+(_bound|_apo|_trimer|_holo)?$", "", name)
        dpath = find_design(tag)

        pred = j.read_atoms(pdbs[0], pred_chains)
        missing = [c for c in pred_chains if len(pred[c]["ca"]) == 0]
        if missing:
            rows.append({"design": name, "status": f"PRED_CHAIN_EMPTY:{','.join(missing)}"})
            continue
        nres = {c: len(pred[c]["ca"]) for c in pred_chains}
        cj = json.load(open(conf[0]))
        plddt = np.load(plddt_f[0])["plddt"]
        pae = np.load(pae_f[0])["pae"] if pae_f else None

        off, per_chain, bounds = 0, {}, {}
        for c in pred_chains:
            n = nres[c]
            per_chain[c] = float(plddt[off:off + n].mean())
            bounds[c] = (off, off + n)
            off += n
        assert off == len(plddt), (off, len(plddt), name)

        row = {
            "design": name,
            "backbone": tag,
            "status": "OK",
            "design_file": dpath.name if dpath else "MISSING",
            "n_copies": len(bchains),
            "n_res_total": int(off),
            "confidence_score": round(float(cj["confidence_score"]), 4),
            "ptm": round(float(cj["ptm"]), 4),
            "iptm_GLOBAL_CONFOUNDED": round(float(cj["iptm"]), 4),
            "complex_plddt": round(float(cj["complex_plddt"]), 4),
            "complex_iplddt": round(float(cj["complex_iplddt"]), 4),
        }
        for c in pred_chains:
            row[f"plddt_{c}"] = round(per_chain[c], 4)
        agg([per_chain[c] for c in bchains], "plddt_binder", row)

        # ---- ⭐⭐ DECONFOUND iptm: 15 pairs, only 9 are ours ----
        cp = cj.get("chains_ptm") or {}
        pci = cj.get("pair_chains_iptm") or {}
        if cp and len(cp) == len(pred_chains):
            idx = {c: str(i) for i, c in enumerate(pred_chains)}
            pair = lambda x, y: float(np.mean([float(pci[idx[x]][idx[y]]),
                                               float(pci[idx[y]][idx[x]])]))
            agg([float(cp[idx[c]]) for c in bchains], "ptm_binder", row)

            if tchains:
                per_copy_best, per_copy_mean, bt_all = [], [], []
                for bc in bchains:
                    vals = [pair(bc, tc) for tc in tchains]
                    for tc, v in zip(tchains, vals):
                        row[f"iptm_{bc}_{tc}"] = round(v, 4)
                    per_copy_best.append(max(vals))
                    per_copy_mean.append(float(np.mean(vals)))
                    bt_all += vals
                # ⭐ THE CONJUNCTION: _worst = the least-engaged copy's best contact
                agg(per_copy_best, "iptm_copy_bestcontact", row)
                agg(per_copy_mean, "iptm_copy_meancontact", row)
                row["iptm_binder_target_mean_all9"] = round(float(np.mean(bt_all)), 4)
                row["n_copies_engaged_iptm_ge_0.3"] = sum(1 for v in per_copy_best if v >= 0.3)

                # ⛔ ONE DIRECTION, NOT SYMMETRISED -- this mirrors j497 (:376) exactly.
                # `pair_chains_iptm` is NOT perfectly symmetric (max |Δ| 1e-3 measured over
                # 30 c1_holo rows), and j497 symmetrises binder<->target but NOT
                # target<->target. Following its convention is deliberate: C1 and C3 numbers
                # are compared to each other, so the asymmetry must be inherited, not fixed
                # here.
                tt = [float(pci[idx[tchains[i]]][idx[tchains[k]]])
                      for i in range(len(tchains)) for k in range(i + 1, len(tchains))]
                if tt:
                    row["iptm_target_target_mean"] = round(float(np.mean(tt)), 4)
                    row["iptm_confound_gap"] = round(
                        float(np.mean(tt)) - float(np.mean(bt_all)), 4)

            if len(bchains) > 1:
                # one direction, same inherited convention as the target<->target block
                bb = [float(pci[idx[bchains[i]]][idx[bchains[k]]])
                      for i in range(len(bchains)) for k in range(i + 1, len(bchains))]
                # ⭐ the designed homotrimer interface itself -- never measured before
                row["iptm_binder_binder_mean"] = round(float(np.mean(bb)), 4)
                row["iptm_binder_binder_worst"] = round(float(np.min(bb)), 4)
        elif cp:
            notes.append(f"{name}: chains_ptm has {len(cp)} entries, {len(pred_chains)} "
                         "chains scored -- per-pair iptm SKIPPED rather than mis-indexed")

        # ---- PAE, per chain, never pooled ----
        if pae is not None:
            intra = []
            for bc in bchains:
                s, e = bounds[bc]
                intra.append(float(pae[s:e, s:e].mean()))
            agg(intra, "pae_binder_intra", row, nd=3)
            if tchains:
                per_copy = []
                for bc in bchains:
                    sb, eb = bounds[bc]
                    v = []
                    for tc in tchains:
                        st, et = bounds[tc]
                        v.append(float(np.mean([pae[sb:eb, st:et].mean(),
                                                pae[st:et, sb:eb].mean()])))
                    per_copy.append(min(v))
                agg(per_copy, "pae_binder_target_bestcontact", row, nd=3)

        # ---- self-consistency: each binder copy vs its OWN design chain ----
        if dpath is not None:
            try:
                des = j.read_atoms(dpath, design_chains)
                rms, tms = [], []
                for pc, dc in zip(bchains, dbchains):
                    P, D = pred[pc]["ca"], des[dc]["ca"]
                    n = min(len(P), len(D))
                    if n < 10:
                        continue
                    rms.append(float(kabsch_rmsd(P[:n], D[:n])))
                    tms.append(float(tm_score(P[:n], D[:n])))
                agg([x for x in rms if x is not None], "rmsd_binder_copy", row, nd=3)
                agg([x for x in tms if x is not None], "tm_binder_copy", row)
                row["sc_status"] = "OK"
            except Exception as exc:                       # noqa: BLE001
                row["sc_status"] = f"SC_FAIL:{type(exc).__name__}"
                notes.append(f"{name}: self-consistency failed -- {exc}")
        else:
            row["sc_status"] = "NO_DESIGN"

        # ---- the designed bridge, PER COPY ----
        if tchains:
            try:
                # ⛔ pins are PER BACKBONE; a cohort-wide list is wrong for 11 of 12 here.
                rpins = plan_pins.get(tag, pins) if plan_pins else pins
                if plan_pins and tag not in plan_pins:
                    notes.append(f"{name}: backbone {tag} NOT in plan.json -- bridge measured "
                                 "over ALL Arg, which is NOT the designed mechanism")
                row["pins_used"] = ",".join(str(x) for x in sorted(rpins)) if rpins else "ALL_ARG"
                for bc in bchains:
                    blk = j.bridge_block(pred, bc, tchains, name, notes,
                                         only_resnums=rpins, prefix=f"pred_{bc}_")
                    row.update(blk)
                percopy = []
                for bc in bchains:
                    v = [row.get(f"pred_{bc}_glu107_{tc}") for tc in tchains]
                    v = [x for x in v if x is not None]
                    if v:
                        percopy.append(min(v))
                # ⛔ _worst here is the WORST copy's BEST bridge = the conjunction
                agg(percopy, "bridge_copy_best", row, nd=3)
                if percopy:
                    row["n_copies_bridged_under_9A"] = sum(1 for v in percopy if v < 9.0)
                    row["n_copies_bridged_under_4A"] = sum(1 for v in percopy if v < 4.0)
            except Exception as exc:                       # noqa: BLE001
                notes.append(f"{name}: bridge block failed -- {exc}")

        rows.append(row)

    keys, seen = [], set()
    for r in rows:
        for k in r:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with open(a.out, "w") as fh:                           # ⛔ overwrite, never append
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")

    ok = [r for r in rows if r.get("status") == "OK"]
    print(f"\nwrote {a.out}  n={len(rows)}  OK={len(ok)}")
    for n in notes[:12]:
        print(f"   ⚠️ {n}")
    if len(notes) > 12:
        print(f"   ... and {len(notes)-12} more notes")
    if not ok:
        print("⛔ ZERO OK rows -- check --pred-dir points at predictions/, not out_<arm>")

    if a.selftest_vs_j497:
        selftest(a.out, a.selftest_vs_j497)


def selftest(mine, theirs):
    """⭐ The only evidence the per-pair indexing is right: reproduce j497 on C1."""
    def load(p):
        h = open(p).readline().rstrip("\n").split("\t")
        return {r["design"]: r for r in
                (dict(zip(h, l.rstrip("\n").split("\t")))
                 for l in open(p).readlines()[1:] if l.strip())}
    A, B = load(mine), load(theirs)
    common = sorted(set(A) & set(B))
    print(f"\n=== SELFTEST vs j497 ===\n  rows mine={len(A)} theirs={len(B)} common={len(common)}")
    if not common:
        print("  ⛔ NO COMMON ROWS -- cannot validate")
        return
    # j497's single-binder names vs this script's generalised names
    pairs = [("iptm_copy_bestcontact_best", "iptm_binder_target_best"),
             ("iptm_copy_meancontact_best", "iptm_binder_target_mean"),
             ("iptm_target_target_mean", "iptm_target_target_mean"),
             ("iptm_confound_gap", "iptm_confound_gap"),
             ("complex_plddt", "complex_plddt"),
             ("ptm", "ptm")]
    bad = 0
    for km, kt in pairs:
        diffs = []
        for d in common:
            try:
                x, y = float(A[d][km]), float(B[d][kt])
            except (KeyError, ValueError):
                continue
            diffs.append(abs(x - y))
        if not diffs:
            print(f"  ⚠️ {km:<34} no comparable values")
            continue
        m = max(diffs)
        flag = "✅" if m < 1e-6 else "⛔"
        if m >= 1e-6:
            bad += 1
        print(f"  {flag} {km:<34} vs {kt:<28} max|Δ| = {m:.2e}  (n={len(diffs)})")
    print("  ⇒ " + ("✅ IDENTICAL -- per-pair indexing validated"
                    if not bad else f"⛔ {bad} column(s) DIFFER -- do NOT trust this on C3"))


if __name__ == "__main__":
    main()
