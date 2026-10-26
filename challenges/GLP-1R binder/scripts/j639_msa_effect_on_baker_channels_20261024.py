#!/usr/bin/env python3
"""j639 -- does the MSA, not the design, drive the Baker channels? The matched control, found on disk.

A matched single-sequence control for this comparison already existed on disk. ⭐ It is NOT the
`j630` arm (that one is MSA-ON: 8,276 receptor rows, 2,548 peptide rows, and its builder asserts
that no chain carries an empty alignment). It is **`J523`**,
2026-10-18, `boltz/j523_ratio_20261018/` -- which holds the SAME GLP-1 + GLP-1R co-fold run under
several conditions, 21 draws each.

## ⭐⭐ WHY THIS ANSWERS `J638` WITHOUT A SINGLE NEW FOLD
`J638` measured a near-total separation between our designs and GLP-1 on ipTM and interface PAE,
but could not attribute it: every control we had was MSA-ON (receptor `ptm` 0.977) while our 960
folds are single-sequence (0.860), on a byte-identical receptor. `J523` breaks that tie because two
of its arms differ in **exactly one factor**:

    active_f260    receptor MSA ON  (apo MSA) + forced 6X18 template
    active_nomsa   receptor MSA OFF (`empty`) + forced 6X18 template   <- the matched arm

⇒ the `active_f260` - `active_nomsa` delta IS the MSA effect on each Baker channel, measured on the
natural ligand, with everything else held.

## ⛔ WHAT THIS STILL DOES NOT MATCH, AND IT MATTERS
⛔ **Both `J523` arms carry a FORCED TEMPLATE** (`glp1r_active_6X18_matched.cif`, `force: true`,
threshold 2.0); our 960 folds have **no template**. So `active_nomsa` is *not* a like-for-like
control for our designs -- it is a like-for-like control for the MSA QUESTION only. ⛔ Do not quote
its absolute values against our designs as though the template were absent.
⚠️ The peptide is 30 aa here (`...WLVKGR`) vs 31 in `j630` (`...WLVKGRG`), and the receptor MSA is
the **apo** MSA, not `j630`'s sliced F260 one. Neither affects the within-J523 delta.
⚠️ **both J523 arms carry L260, while our 960 folds carry F260.** Stated because it is a
construct difference, ⛔ NOT a proposal to change anything;
it is identical across the pair, so the delta is clean.
**MSA use verified in the runs' own bytes**, not from the yaml: `out_active/processed/msa`
holds **11** files, `out_active_nomsa/processed/msa` holds **0**.

## ⛔ CHAIN ORDER -- the bug `J638` records
`J523` yamls are **RECEPTOR-FIRST** (`id: A` = 409-mer, `id: B` = peptide), like the other controls
and unlike our design fastas. Block offsets are passed explicitly; nothing is inferred from position.

⛔ Reads folds already on disk. No folding, nothing launched.

usage: j639_msa_effect_on_baker_channels_20261024.py --out <fresh.tsv>
"""
import argparse
import csv
import glob
import importlib.util as iu
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
C3 = os.path.dirname(HERE)
J523 = f"{C3}/boltz/j523_ratio_20261018"
TGT_LEN = 409
PEP_LEN = 30
CUTS = dict(plddt=0.90, iptm=0.90, pae=1.2)

_spec = iu.spec_from_file_location(
    "j637", os.path.join(HERE, "j637_baker_filter_c3_20261024.py"))
j637 = iu.module_from_spec(_spec)
_spec.loader.exec_module(j637)          # safe: it has a __main__ guard

# the J523 arms worth scoring, and what each one IS
# ⛔ `active_f260`/`inactive_f260` were WRITTEN but never RUN -- only `out_active`, `out_active_rep`,
# `out_inactive`, `out_probe` and `out_active_nomsa` exist. The MSA-ON partner is therefore
# `out_active`, and it pairs cleanly: its receptor and peptide are BYTE-IDENTICAL to
# `active_nomsa`'s (409 + 30 aa, 0 mismatches), so the arms differ in the MSA alone.
ARMS = [
    ("j523_active_msa_tmpl",   "out_active",       "MSA ON  + template"),
    ("j523_active_NOMSA_tmpl", "out_active_nomsa", "MSA OFF + template  <- MATCHED ARM"),
]


def score_arm(outdir, label):
    rows = []
    pat = f"{J523}/{outdir}/*/predictions/L01_GLP-17-36NH2/confidence_*.json"
    for conf in sorted(glob.glob(pat)):
        d = os.path.dirname(conf)
        stem = os.path.basename(conf)[len("confidence_"):-len(".json")]
        pae = f"{d}/pae_{stem}.npz"
        assert os.path.exists(pae), f"missing pae for {stem}"
        m = re.search(r"_model_(\d+)$", stem)
        # ⛔ receptor-first: binder block starts AFTER the 409-mer
        r = j637.score(conf, pae, label, "glp1_7-36", f"model_{m.group(1)}",
                       binder_first=False)
        assert r["binder_len"] == PEP_LEN, \
            f"{stem}: binder_len {r['binder_len']} != {PEP_LEN} -- check the chain layout"
        rows.append(r)
    assert rows, f"⛔ no folds found under {pat}"
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    assert not os.path.exists(a.out), f"⛔ exists, pick a fresh --out: {a.out}"

    rows = []
    for label, outdir, _desc in ARMS:
        rows += score_arm(outdir, label)

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    def med(lab, k):
        v = [r[k] for r in rows if r["cohort"] == lab]
        return float(np.median(v)), len(v)

    print("GLP-1 (7-36) + GLP-1R, 409-mer, SAME peptide, SAME template -- only the MSA differs.\n")
    print(f"  {'arm':<24} {'n':>3} {'cplDDT':>8} {'ipTM':>8} {'PAEmin':>9} {'PAEmean':>9}")
    for label, _o, desc in ARMS:
        n = med(label, "iptm_pair")[1]
        print(f"  {label:<24} {n:>3} "
              f"{med(label,'complex_plddt')[0]:>8.3f} {med(label,'iptm_pair')[0]:>8.3f} "
              f"{med(label,'pae_iface_min')[0]:>9.3f} {med(label,'pae_iface_mean')[0]:>9.3f}   {desc}")

    on, off = ARMS[0][0], ARMS[1][0]
    print(f"\n  {'Δ (MSA ON - MSA OFF)':<24} {'':>3} "
          f"{med(on,'complex_plddt')[0]-med(off,'complex_plddt')[0]:>8.3f} "
          f"{med(on,'iptm_pair')[0]-med(off,'iptm_pair')[0]:>8.3f} "
          f"{med(on,'pae_iface_min')[0]-med(off,'pae_iface_min')[0]:>9.3f} "
          f"{med(on,'pae_iface_mean')[0]-med(off,'pae_iface_mean')[0]:>9.3f}")

    # ---- the question J638 could not answer ------------------------------------------------
    print("\n⇒ THE TEST: does single-sequence GLP-1 still clear the two channels that survived")
    print("  their controls in J638 (ipTM > 0.90, interface PAE min < 1.2 A)?")
    offr = [r for r in rows if r["cohort"] == off]
    pi = sum(r["iptm_pair"] > CUTS["iptm"] for r in offr)
    pp = sum(r["pae_iface_min"] < CUTS["pae"] for r in offr)
    pb = sum(r["iptm_pair"] > CUTS["iptm"] and r["pae_iface_min"] < CUTS["pae"] for r in offr)
    print(f"    ipTM > 0.90          {pi}/{len(offr)}   (range "
          f"{min(r['iptm_pair'] for r in offr):.3f}-{max(r['iptm_pair'] for r in offr):.3f})")
    print(f"    PAE min < 1.2 A      {pp}/{len(offr)}   (range "
          f"{min(r['pae_iface_min'] for r in offr):.3f}-{max(r['pae_iface_min'] for r in offr):.3f})")
    print(f"    BOTH                 {pb}/{len(offr)}")

    # ---- and against OUR designs' best, from j637's own output -----------------------------
    tsv = f"{C3}/data/j637_baker_filter_c3_20261024_v2.tsv"
    if os.path.exists(tsv):
        ours = [r for r in csv.DictReader(open(tsv), delimiter="\t")
                if r["cohort"] == "c3_design"]
        bi = max(float(r["iptm_pair"]) for r in ours)
        bp = min(float(r["pae_iface_min"]) for r in ours)
        print(f"\n  our 960 single-sequence design folds, BEST: ipTM {bi:.3f} · PAE min {bp:.3f}")
        wi = min(r["iptm_pair"] for r in offr)
        wp = max(r["pae_iface_min"] for r in offr)
        print(f"  single-sequence GLP-1,      WORST: ipTM {wi:.3f} · PAE min {wp:.3f}")
        print(f"  ⇒ our best still below GLP-1's worst?  "
              f"ipTM {'YES' if bi < wi else '⛔ NO'} · PAE min {'YES' if bp > wp else '⛔ NO'}")
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
