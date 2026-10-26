#!/usr/bin/env python3
"""j637 -- apply the Baker cysteine-protease paper's ROUND-3 filter to the C3 (GLP-1R) folds.

WHY NOW. The round-3 in-silico filter from a recent de novo enzyme-design paper was run once
before, on the challenge-2 cohort (`J509`, record `J507`) -- 0/110, with a published de novo binder
of measured affinity passing as a positive control.
⛔⛔ That record says **"C3 CANNOT BE FILTERED THIS WAY AT ALL"** because `v50` had no
Boltz fold of any kind. **That is no longer true.** `J628` folded 96 designs x 10 seeds against the
continuous 409-aa receptor, so every channel the filter needs now exists on disk. ⇒ this script is
the GLP-1R arm of J507, unblocked by J628. It reads folds already on disk; nothing is run.

## THE PUBLISHED THRESHOLDS [Methods p10 and Extended Data Fig. 1b p15 of that paper,
##   recorded in this group's own reading note on it]
  complex pLDDT > 90 . ipTM > 0.90 . min PAE < 1.2 . apo pLDDT > 90 . apo-to-complex CA rmsd < 0.8 A
⛔ Only the PUBLISHED THRESHOLDS are used. Their repository is under a standing **NEVER RUN**
prohibition in this group's own working rules and is not touched.

## ⛔⛔ THE CONFOUND THAT MUST BE PRINTED WITH EVERY NUMBER -- ⛔ do NOT quote the control gap alone
The only known-binder controls we hold were folded with the **receptor MSA ON (7,938 sequences)**;
all 960 of our designs are **single-sequence, 0 MSA on both chains**, verified in the
input bytes. pLDDT, ipTM and PAE are all MSA-sensitive, so a control-vs-design gap
here confounds design quality with the MSA condition. ⛔ This is NOT the same question an
earlier record settled: that ruled the receptor MSA out as the cause of *epitope instability*, which
says nothing about whether an MSA inflates *confidence* channels. ⇒ the control establishes that
the filter DISCRIMINATES under Boltz-2 on this target; it does ⛔ NOT establish that our designs
would fail a matched-MSA comparison. A matched single-sequence GLP-1 control is the measurement
that would close this, and it is not run here.

## ⛔ "min PAE" IS AMBIGUOUS -- both readings are computed, neither hidden (as in j509)
The paper prints *"min PAE < 1.2"* without defining the reduction. Reported here: the literal **MIN**
element over the binder x target block, and the **MEAN** (what the Baker binder literature usually
calls `pae_interaction`). j509 found the literal min is the reading that reproduces the paper's own
behaviour, because a real binder (LCB1) passed on min and failed on mean.

## ⛔ n IS 96 DESIGNS, NOT 960 FOLDS
The 10 shards are 10 SEEDS of the same 96 designs. A per-fold pass rate would count each design ten
times. Both are printed: per-fold (960) and per-design best-of-10 (96), which is the generous
reading and the one to quote against the designs.

usage: j637_baker_filter_c3_20261024.py --out <fresh.tsv>
"""
import argparse
import csv
import glob
import json
import os
import re

import numpy as np

SCRATCH = "/path/to/your/working/tree"  # PLACEHOLDER PATH -- not a working default; set it for your own tree.
J628 = f"{SCRATCH}/j628_unpack"
J630 = f"{SCRATCH}/j630_glp1_cofold_21"
J318 = f"{SCRATCH}/glp1r_cofold_20261013"
TGT_LEN = 409                      # the continuous GLP-1R construct, UniProt 24-432
CUTS = dict(plddt=0.90, iptm=0.90, pae=1.2)


def iface_pae(pae, bs, nb, ts, tgt):
    """-> (MIN element, MEAN) over the binder x target block, symmetrised.

    ⛔⛔ THE CHAIN ORDER IS NOT THE SAME IN EVERY COHORT AND ASSUMING IT IS WAS A REAL BUG.
    Our design fastas are BINDER-first (`>A` binder, `>B` receptor); both control yamls are
    RECEPTOR-first (`id: A` = the 409-mer, `id: P` = the peptide). The first version of this
    script took rows [0:nb) as the binder for everything, so for the controls it measured an
    INTRA-RECEPTOR block -- adjacent residues, PAE ~0.3 A -- and reported it as a spectacular
    interface. ⇒ the blocks are passed in explicitly by the caller, never inferred from position
    (chain order has silently changed a conclusion in this group's work before).
    """
    A, B = pae[bs:bs + nb, ts:ts + tgt], pae[ts:ts + tgt, bs:bs + nb]
    return float(min(A.min(), B.min())), float(np.mean([A.mean(), B.mean()]))


def score(conf_path, pae_path, cohort, name, seed, binder_first):
    c = json.load(open(conf_path))
    pae = np.load(pae_path)["pae"]
    n = pae.shape[0]
    nb = n - TGT_LEN
    assert 20 < nb < 120, f"{name}: implausible binder length {nb} from PAE edge {n}"
    # explicit block offsets, from the cohort's OWN input layout
    bs, ts = (0, nb) if binder_first else (TGT_LEN, 0)
    # ⭐ pair_chains_iptm, not the global iptm, which averages EVERY chain pair.
    # Here the target is a SINGLE chain so the two nearly coincide, but the pair value is the
    # one that means "binder vs target" and it is what j509 used.
    pc = c.get("pair_chains_iptm", {})
    pair = max(float(pc["0"]["1"]), float(pc["1"]["0"])) if pc else float(c["iptm"])
    pmin, pmean = iface_pae(pae, bs, nb, ts, TGT_LEN)
    return dict(cohort=cohort, name=name, seed=seed, binder_len=nb,
                binder_first=int(binder_first),
                complex_plddt=round(float(c["complex_plddt"]), 4),
                iptm_global=round(float(c["iptm"]), 4),
                iptm_pair=round(pair, 4),
                pae_iface_min=round(pmin, 3), pae_iface_mean=round(pmean, 3))


def assert_same_receptor():
    """⛔ The whole control comparison rests on both cohorts holding the SAME receptor sequence.
    If they differ, a confidence gap is not attributable to the MSA or to the designs. Asserted
    from the input bytes, never eyeballed."""
    design_fa = sorted(glob.glob(f"{SCRATCH}/j628_continuous/inputs/*.fasta"))[0]
    seqs, cur = {}, None
    for ln in open(design_fa):
        ln = ln.strip()
        if ln.startswith(">"):
            cur = ln[1:].split("|")[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += ln
    ours = seqs["B"]
    ctrl = None
    for ln in open(f"{SCRATCH}/j630_glp1_cofold_21/glp1_glp1r_F260.yaml"):
        s = ln.strip()
        if s.startswith("sequence:") and len(s.split(None, 1)[1]) == TGT_LEN:
            ctrl = s.split(None, 1)[1]
            break
    assert ours is not None and ctrl is not None, "could not read both receptor sequences"
    assert len(ours) == TGT_LEN, f"design receptor is {len(ours)} aa, expected {TGT_LEN}"
    assert ours == ctrl, ("⛔⛔ RECEPTOR SEQUENCES DIFFER between designs and control -- the "
                          "control comparison is invalid. Resolve before reading any number.")
    print(f"✅ receptor sequence identical across cohorts ({TGT_LEN} aa, F260 at "
          f"{ours.find('F', 230)}): the control differs from our folds in MSA and binder only.\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    assert not os.path.exists(a.out), f"⛔ exists, pick a fresh --out: {a.out}"

    assert_same_receptor()
    rows = []

    # ---- (1) our 960 C3 folds: 96 designs x 10 seed-shards, SINGLE SEQUENCE ----------------
    for conf in sorted(glob.glob(f"{J628}/shard*/*/predictions/*/confidence_*.json")):
        d = os.path.dirname(conf)
        pae = glob.glob(f"{d}/pae_*.npz")
        assert len(pae) == 1, f"{d}: expected 1 pae npz, found {len(pae)}"
        stem = os.path.basename(conf)[len("confidence_"):-len("_model_0.json")]
        design = re.sub(r"__id1$", "", stem)
        seed = re.search(r"/(shard\d+)/", conf).group(1)
        rows.append(score(conf, pae[0], "c3_design", design, seed, binder_first=True))
    ours = [r for r in rows if r["cohort"] == "c3_design"]
    assert len(ours) == 960, f"expected 960 design folds, found {len(ours)}"
    assert len({r["name"] for r in ours}) == 96, "expected 96 distinct designs"

    # ---- (2) controls. ⛔⛔ ALL MSA-ON -- see the docstring confound ------------------------
    for conf in sorted(glob.glob(f"{J630}/out/seed*/*/predictions/*/confidence_*.json")):
        d = os.path.dirname(conf)
        pae = glob.glob(f"{d}/pae_*.npz")
        assert len(pae) == 1, f"{d}: expected 1 pae npz, found {len(pae)}"
        seed = re.search(r"/(seed\d+)/", conf).group(1)
        # ⛔ RECEPTOR-FIRST: `id: A` is the 409-mer, `id: P` the peptide. See iface_pae's note.
        rows.append(score(conf, pae[0], "ctrl_glp1_msa", "glp1_natural", seed,
                          binder_first=False))

    for conf in sorted(glob.glob(f"{J318}/**/confidence_*.json", recursive=True)):
        d = os.path.dirname(conf)
        stem = os.path.basename(conf)[len("confidence_"):-len(".json")]
        pae = glob.glob(f"{d}/pae_{stem}.npz")
        assert len(pae) == 1, f"{d}: expected pae for {stem}, found {len(pae)}"
        arm = "glp1" if stem.startswith("glp1") else "lira"
        rows.append(score(conf, pae[0], f"ctrl_{arm}_msa_j318", f"j318_{arm}", stem,
                          binder_first=False))   # ⛔ receptor-first, same as j630

    for r in rows:
        r["pass_plddt"] = int(r["complex_plddt"] > CUTS["plddt"])
        r["pass_iptm"] = int(r["iptm_pair"] > CUTS["iptm"])
        r["pass_pae_min"] = int(r["pae_iface_min"] < CUTS["pae"])
        r["pass_pae_mean"] = int(r["pae_iface_mean"] < CUTS["pae"])
        r["pass_all_minpae"] = int(r["pass_plddt"] and r["pass_iptm"] and r["pass_pae_min"])

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ------------------------------- report -------------------------------------------------
    print(f"Baker round-3 thresholds: complex pLDDT > {CUTS['plddt']} . "
          f"ipTM > {CUTS['iptm']} . min PAE < {CUTS['pae']} A\n")

    print("== OUR C3 DESIGNS, PER FOLD (n=960 folds) ==")
    for lbl, key, best in [
            ("complex pLDDT > 0.90", "pass_plddt", max(r["complex_plddt"] for r in ours)),
            ("ipTM binder-target > 0.90", "pass_iptm", max(r["iptm_pair"] for r in ours)),
            ("interface PAE MIN  < 1.2 A", "pass_pae_min", min(r["pae_iface_min"] for r in ours)),
            ("interface PAE MEAN < 1.2 A", "pass_pae_mean", min(r["pae_iface_mean"] for r in ours))]:
        print(f"  {lbl:<30} PASS {sum(r[key] for r in ours):>4}/960   best {best:8.3f}")
    print(f"  {'CONJUNCTION (with MIN PAE)':<30} PASS "
          f"{sum(r['pass_all_minpae'] for r in ours):>4}/960")

    # ⭐ per-DESIGN best-of-10: the generous reading, and the one where n is honest
    by_design = {}
    for r in ours:
        d = by_design.setdefault(r["name"], [])
        d.append(r)
    print(f"\n== OUR C3 DESIGNS, BEST-OF-10-SEEDS PER DESIGN (n={len(by_design)} designs) ==")
    dbest = []
    for name, rs in by_design.items():
        dbest.append(dict(name=name,
                          plddt=max(r["complex_plddt"] for r in rs),
                          iptm=max(r["iptm_pair"] for r in rs),
                          pmin=min(r["pae_iface_min"] for r in rs),
                          anypass=max(r["pass_all_minpae"] for r in rs)))
    print(f"  {'complex pLDDT > 0.90':<30} PASS "
          f"{sum(d['plddt'] > CUTS['plddt'] for d in dbest):>4}/96   best "
          f"{max(d['plddt'] for d in dbest):8.3f}")
    print(f"  {'ipTM binder-target > 0.90':<30} PASS "
          f"{sum(d['iptm'] > CUTS['iptm'] for d in dbest):>4}/96   best "
          f"{max(d['iptm'] for d in dbest):8.3f}")
    print(f"  {'interface PAE MIN  < 1.2 A':<30} PASS "
          f"{sum(d['pmin'] < CUTS['pae'] for d in dbest):>4}/96   best "
          f"{min(d['pmin'] for d in dbest):8.3f}")
    print(f"  {'CONJUNCTION, any seed':<30} PASS "
          f"{sum(d['anypass'] for d in dbest):>4}/96")

    print("\n== CONTROLS -- ⛔⛔ MSA-ON (7,938 receptor seqs) vs our 0. NOT a matched comparison ==")
    print(f"  {'arm':<22} {'n':>3} {'plddt':>7} {'ipTM':>7} {'PAEmin':>8} {'PAEmean':>8}  pass")
    for coh in ("ctrl_glp1_msa", "ctrl_glp1_msa_j318", "ctrl_lira_msa_j318"):
        g = [r for r in rows if r["cohort"] == coh]
        if not g:
            continue
        print(f"  {coh:<22} {len(g):>3} "
              f"{np.median([r['complex_plddt'] for r in g]):>7.3f} "
              f"{np.median([r['iptm_pair'] for r in g]):>7.3f} "
              f"{np.median([r['pae_iface_min'] for r in g]):>8.3f} "
              f"{np.median([r['pae_iface_mean'] for r in g]):>8.3f}  "
              f"{sum(r['pass_all_minpae'] for r in g)}/{len(g)}")
    print("  (medians; per-fold rows are all in the TSV)")

    print("\n⇒ the 10 best C3 designs by worst-channel shortfall:")
    def gap(d):
        return max((CUTS["plddt"] - d["plddt"]) / CUTS["plddt"],
                   (CUTS["iptm"] - d["iptm"]) / CUTS["iptm"],
                   (d["pmin"] - CUTS["pae"]) / CUTS["pae"])
    for i, d in enumerate(sorted(dbest, key=gap)[:10], 1):
        worst = max([("pLDDT", (CUTS["plddt"] - d["plddt"]) / CUTS["plddt"]),
                     ("ipTM", (CUTS["iptm"] - d["iptm"]) / CUTS["iptm"]),
                     ("PAE", (d["pmin"] - CUTS["pae"]) / CUTS["pae"])], key=lambda x: x[1])
        print(f"  {i:>2}. {d['name']:<42} pLDDT {d['plddt']:.3f}  ipTM {d['iptm']:.3f}  "
              f"PAEmin {d['pmin']:7.3f}   worst {worst[0]} +{worst[1]*100:.0f}%")

    print(f"\n-> {a.out}")
    print("\n⛔ apo pLDDT and apo-to-complex CA rmsd are NOT evaluated: they need a second fold of")
    print("   each binder alone. Whether that can change the verdict depends on the result above.")


if __name__ == "__main__":
    main()
