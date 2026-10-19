#!/usr/bin/env python3
"""j506 -- score the DE NOVO control (`j505`): epitope recovery + GROUND-TRUTH pose RMSD.

WHY A NEW SCORER. `J502` §0's epitope-recovery table was computed by a script that was **never
saved** -- `scripts/` holds `j501` (inputs) and `j497` (confidence metrics), and no epitope scorer.
⇒ those numbers are not reproducible from the repo as it stands. This script is the durable
instrument, and it is **self-validated**: scoring the CRYSTAL complex against its own reference
epitope must return recall **1.000**, which is asserted. A scorer that cannot recover a perfect
input cannot be trusted on a prediction.

## THE METRIC, stated so it can be checked -- the same FORM as `J502` §0
reference epitope  = target residues within `CUT` (5.0 Å, heavy atom) of the binder IN THE CRYSTAL
predicted epitope  = target residues within `CUT` of the binder IN THE PREDICTION
recall             = |predicted ∩ reference| / |reference|
chance             = |predicted| / |target residues|      ⇒ random placement of the same-size patch
enrichment         = recall / chance
⭐ `chance` is normalised by the PREDICTED patch size on purpose: a prediction that smears contacts
over the whole target would score high recall for free, and dividing by its own footprint removes
that. ⛔ Verified against `J502` §0's published arithmetic: positive 57/456 = 0.125 chance and
0.781/0.125 = 6.25×; negative 23/456 = 0.050 and 0.250/0.050 = 4.96×. Same definition.

## ⭐⭐ AND THE METRIC `J502` COULD NOT HAVE -- ground truth
`7JZU` is a solved LCB1–RBD complex, so after superposing the PREDICTED target onto the CRYSTAL
target (CA, Kabsch) the predicted binder can be compared to the crystal binder directly:
  `rmsd_binder_truth` = CA RMSD of the predicted binder to the crystal binder, in the target frame
This is a **direct** test of whether the binding mode was reproduced, not a proxy. ⛔ It exists only
because a crystal exists; for our own designs it does not, which is exactly why `J497` had to
measure displacement from the *designed* site instead. ⚠️ Defined only for the positive arm.

⛔ ROLE DETECTION BY SEQUENCE IDENTITY to the known target, never by chain letter. Residue mapping prediction→crystal is by INDEX into the parsed
chain, and the chain's sequence is asserted identical to the crystal's first.

⛔ WHAT THIS IS NOT. A prediction is not an affinity measurement. a global ipTM on a multimeric target and `J502` §4 already ruled `iptm` out as a gate; it is reported for
completeness only. Nothing here is a submission gate.

usage: j506_score_denovo_control_20261018.py --out <fresh.tsv>
"""
import argparse
import csv
import json
import os

import numpy as np

REPO = "<placeholder>"
SEVENJZU = f"{REPO}/data/J130_inputs_20261010_7jzu.pdb"
PRED = f"{REPO}/boltz/j505_denovo_control_20261018/control_out/boltz_results_inputs/predictions"
CUT = 5.0
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
ARMS = [("pos_lcb1_denovo", "POSITIVE de novo -- LCB1, designed, measured affinity"),
        ("neg_tnfa_vs_rbd", "NEGATIVE -- TNFa protomer, folded, no RBD affinity")]


def read_pdb(path):
    """-> {chain: {'seq': str, 'resnums': [...], 'heavy': {resnum: Nx3}, 'ca': {resnum: 3}}}"""
    ch = {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        el = line[76:78].strip() or name[:1]
        if el == "H":
            continue
        c, rn, res = line[21], int(line[22:26]), line[17:20].strip()
        if res not in AA3:
            continue
        d = ch.setdefault(c, dict(seq={}, heavy={}, ca={}))
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        d["seq"][rn] = AA3[res]
        d["heavy"].setdefault(rn, []).append(xyz)
        if name == "CA":
            d["ca"][rn] = xyz
    out = {}
    for c, d in ch.items():
        rns = sorted(d["seq"])
        out[c] = dict(seq="".join(d["seq"][r] for r in rns), resnums=rns,
                      heavy={r: np.array(v) for r, v in d["heavy"].items()}, ca=d["ca"])
    return out


def contacts(binder, target):
    """Target resnums with any heavy atom within CUT of any binder heavy atom."""
    B = np.vstack([binder["heavy"][r] for r in binder["resnums"]])
    hit = set()
    for r in target["resnums"]:
        T = target["heavy"][r]
        if np.min(np.linalg.norm(T[:, None, :] - B[None, :, :], axis=-1)) <= CUT:
            hit.add(r)
    return hit


def kabsch(P, Q):
    """Rotation+translation taking P onto Q (both Nx3, paired)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def detect(pred, target_seq):
    """⛔ Role by SEQUENCE IDENTITY to the known target, not by letter."""
    tgt = [c for c, d in pred.items() if d["seq"] == target_seq]
    bnd = [c for c, d in pred.items() if d["seq"] != target_seq]
    assert len(tgt) == 1 and len(bnd) == 1, (
        f"target {tgt} binder {bnd}; lengths { {c: len(d['seq']) for c, d in pred.items()} }")
    return bnd[0], tgt[0]


def conf(tag):
    p = f"{PRED}/{tag}_bound/confidence_{tag}_bound_model_0.json"
    return json.load(open(p)) if os.path.exists(p) else {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="⛔ FRESH path -- scripts APPEND (STOP S5)")
    a = ap.parse_args()
    assert not os.path.exists(a.out), f"⛔ --out exists, refusing to merge campaigns: {a.out}"

    xtal = read_pdb(SEVENJZU)
    lcb1_x, rbd_x = xtal["A"], xtal["B"]
    assert len(lcb1_x["seq"]) == 55 and len(rbd_x["seq"]) == 193, "7JZU chain lengths unexpected"
    ref = contacts(lcb1_x, rbd_x)
    n_tgt = len(rbd_x["resnums"])
    print(f"CRYSTAL 7JZU: LCB1 {len(lcb1_x['seq'])} aa, RBD {n_tgt} aa")
    print(f"reference epitope (RBD within {CUT} A heavy-atom of LCB1): {len(ref)} residues")
    print(f"  resnums: {sorted(ref)}\n")

    # ⭐ SELF-VALIDATION: the crystal scored against its own reference MUST be perfect.
    pred_x = contacts(lcb1_x, rbd_x)
    rec_x = len(pred_x & ref) / len(ref)
    print(f"== scorer self-validation ==\n  crystal vs its own reference: recall {rec_x:.3f}")
    assert abs(rec_x - 1.0) < 1e-9, f"⛔⛔ scorer cannot recover a perfect input ({rec_x})"
    print("  ⇒ OK, recall 1.000 on a perfect input.\n")

    rows = []
    for tag, label in ARMS:
        pdb = f"{PRED}/{tag}_bound/{tag}_bound_model_0.pdb"
        assert os.path.exists(pdb), f"⛔ missing prediction: {pdb}"
        pr = read_pdb(pdb)
        bc, tc = detect(pr, rbd_x["seq"])
        pb, pt = pr[bc], pr[tc]
        assert pt["seq"] == rbd_x["seq"], "predicted target sequence != crystal RBD"

        # predicted contacts, mapped to CRYSTAL resnums by index (sequences are identical)
        hit = contacts(pb, pt)
        idx = {r: i for i, r in enumerate(pt["resnums"])}
        hit_x = {rbd_x["resnums"][idx[r]] for r in hit}
        ov = hit_x & ref
        recall = len(ov) / len(ref)
        chance = len(hit_x) / n_tgt
        enr = recall / chance if chance else float("nan")

        # ground truth: superpose predicted RBD -> crystal RBD, then binder RMSD (positive only)
        rmsd_truth = ""
        if tag.startswith("pos"):
            P = np.array([pt["ca"][r] for r in pt["resnums"]])
            Q = np.array([rbd_x["ca"][r] for r in rbd_x["resnums"]])
            R, t = kabsch(P, Q)
            assert len(pb["resnums"]) == len(lcb1_x["resnums"]), "binder length mismatch vs crystal"
            Bp = np.array([pb["ca"][r] for r in pb["resnums"]]) @ R.T + t
            Bx = np.array([lcb1_x["ca"][r] for r in lcb1_x["resnums"]])
            rmsd_truth = round(float(np.sqrt(((Bp - Bx) ** 2).sum(1).mean())), 2)
            tgt_fit = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
            print(f"  [{tag}] target superposition RMSD {tgt_fit:.2f} A "
                  f"(the frame the binder RMSD is measured in)")

        c = conf(tag)
        rows.append(dict(arm=tag, label=label, binder_chain=bc, target_chain=tc,
                         binder_len=len(pb["seq"]), n_ref=len(ref), n_pred=len(hit_x),
                         n_overlap=len(ov), recall=round(recall, 3), chance=round(chance, 3),
                         enrichment=round(enr, 2), rmsd_binder_truth=rmsd_truth,
                         iptm=round(c.get("iptm", float("nan")), 3),
                         ptm=round(c.get("ptm", float("nan")), 3),
                         complex_plddt=round(c.get("complex_plddt", float("nan")), 3)))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\n== DE NOVO CONTROL -- epitope recovery (ref {len(ref)} res, target {n_tgt} res) ==")
    print(f"{'arm':<18} {'pred':>5} {'ov':>4} {'recall':>7} {'chance':>7} {'enrich':>7} "
          f"{'RMSDtruth':>10} {'iptm':>6}")
    for r in rows:
        print(f"{r['arm']:<18} {r['n_pred']:>5} {r['n_overlap']:>4} {r['recall']:>7.3f} "
              f"{r['chance']:>7.3f} {r['enrichment']:>7.2f} {str(r['rmsd_binder_truth']):>10} "
              f"{r['iptm']:>6.3f}")
    print(f"\n-> {a.out}")
    print("\n⚠️ For comparison, J502 §0 (computed by a script that was NOT saved -- quoted, not")
    print("   reproduced here): TNFR2 6.25x, RBD-vs-TNFa 4.96x, our 110 v46 folds 2.78x median.")
    print("⛔ A prediction is not an affinity measurement, and iptm is not a gate (J502 §4).")


if __name__ == "__main__":
    main()
