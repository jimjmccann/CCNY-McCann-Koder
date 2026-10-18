#!/usr/bin/env python3
"""J465 -- score Boltz apo/bound folds of the v41c C3 SolubleMPNN designs.

Measures, per design:
  - complex and per-chain pLDDT (from the Boltz plddt npz, sliced by chain length)
  - ptm / iptm / chains_ptm (from the confidence json)
  - mean intra- and inter-chain PAE
  - design-vs-prediction self-consistency CA-RMSD against the RFd3 v41c design
    backbone, chains named EXPLICITLY.

Chain convention (verified from the files, never assumed):
  design PDB  input_best_c3_s00_11.pdb : binder A,B,C (147 aa) + target a,b,c (152 aa)
  Boltz apo   : A,B,C  (binder only)
  Boltz bound : A,B,C binder + D,E,F target

Usage: j465_score.py <predictions_dir> <mode: apo|bound> <out_tsv>
"""
import json
import sys
from pathlib import Path

import numpy as np

# PLACEHOLDER PATH -- not a working default; set it for your own tree.
DESIGN = "/path/to/the/C3/design/input_best_c3_s00_11.pdb"
BINDER_CHAINS = ["A", "B", "C"]
# Boltz renames the target chains; map prediction chain -> design chain.
TARGET_MAP = {"D": "a", "E": "b", "F": "c"}


def read_ca(path, chains):
    """CA coords per chain, in file order. Returns {chain: (N,3) array}."""
    out = {c: [] for c in chains}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16] != " CA ":
                continue
            ch = line[21]
            if ch not in out:
                continue
            out[ch].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return {c: np.asarray(v, dtype=float) for c, v in out.items()}


def kabsch_rmsd(P, Q):
    """RMSD of P onto Q after optimal superposition. P,Q are (N,3), paired."""
    assert P.shape == Q.shape, (P.shape, Q.shape)
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    Pr = Pc @ R
    return float(np.sqrt(((Pr - Qc) ** 2).sum() / len(P)))


def _superpose(P, Q):
    """Return rotation R and translations so that (P - pc) @ R + qc fits Q."""
    pc, qc = P.mean(0), Q.mean(0)
    V, S, Wt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, pc, qc


def tm_score(P, Q):
    """TM-score of a FIXED 1:1 correspondence, normalised by L.

    Zhang-Skolnick iterative-extension heuristic: seed on fragments, then
    re-superpose on the subset within a growing cutoff until it converges.
    A plain RMSD superposition underestimates TM, so the search matters.
    """
    L = len(P)
    d0 = max(1.24 * (L - 15) ** (1.0 / 3.0) - 1.8, 0.5)

    def score(R, pc, qc):
        d = np.linalg.norm((P - pc) @ R + qc - Q, axis=1)
        return float((1.0 / (1.0 + (d / d0) ** 2)).sum() / L), d

    best = 0.0
    seeds = []
    Lf = L
    while Lf >= 4:
        for s in range(0, L - Lf + 1, max(1, Lf // 2)):
            seeds.append((s, s + Lf))
        Lf //= 2
    for s, e in seeds:
        idx = np.arange(s, e)
        for _ in range(20):
            R, pc, qc = _superpose(P[idx], Q[idx])
            sc, d = score(R, pc, qc)
            best = max(best, sc)
            cut = d0 + 1.0
            new = np.where(d < cut)[0]
            while len(new) < 4:
                cut += 0.5
                new = np.where(d < cut)[0]
            if len(new) == len(idx) and np.array_equal(new, idx):
                break
            idx = new
    return best


def main():
    pred_dir = Path(sys.argv[1])
    mode = sys.argv[2]
    out_tsv = Path(sys.argv[3])

    des = read_ca(DESIGN, BINDER_CHAINS + list(TARGET_MAP.values()))
    for c in BINDER_CHAINS:
        assert len(des[c]) == 147, (c, len(des[c]))

    rows = []
    for sub in sorted(pred_dir.iterdir()):
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
        pdb = pdbs[0]
        if pdb.stat().st_size == 0:
            rows.append({"design": name, "status": "EMPTY_PDB"})
            continue

        chains = BINDER_CHAINS + (list(TARGET_MAP) if mode == "bound" else [])
        pred = read_ca(pdb, chains)
        nres = {c: len(pred[c]) for c in chains}

        cj = json.load(open(conf[0]))
        plddt = np.load(plddt_f[0])["plddt"]
        pae = np.load(pae_f[0])["pae"] if pae_f else None

        # slice plddt by chain in file order
        off = 0
        per_chain_plddt = {}
        bounds = {}
        for c in chains:
            n = nres[c]
            per_chain_plddt[c] = float(plddt[off:off + n].mean())
            bounds[c] = (off, off + n)
            off += n
        assert off == len(plddt), (off, len(plddt), name)

        row = {
            "design": name,
            "status": "OK",
            "n_res_total": int(off),
            "complex_plddt": round(float(cj["complex_plddt"]), 4),
            "complex_iplddt": round(float(cj["complex_iplddt"]), 4),
            "ptm": round(float(cj["ptm"]), 4),
            "iptm": round(float(cj["iptm"]), 4),
            "confidence_score": round(float(cj["confidence_score"]), 4),
        }
        for c in chains:
            row[f"plddt_{c}"] = round(per_chain_plddt[c], 4)
        row["plddt_binder_mean"] = round(
            float(np.mean([per_chain_plddt[c] for c in BINDER_CHAINS])), 4)

        if pae is not None:
            # binder intra (block-diagonal A,B,C) vs binder-binder inter
            intra, inter = [], []
            for i, ci in enumerate(BINDER_CHAINS):
                si, ei = bounds[ci]
                for j, cj2 in enumerate(BINDER_CHAINS):
                    sj, ej = bounds[cj2]
                    blk = pae[si:ei, sj:ej]
                    (intra if i == j else inter).append(blk.mean())
            row["pae_binder_intra"] = round(float(np.mean(intra)), 3)
            row["pae_binder_inter"] = round(float(np.mean(inter)), 3)
            if mode == "bound":
                bt = []
                for ci in BINDER_CHAINS:
                    si, ei = bounds[ci]
                    for ct in TARGET_MAP:
                        sj, ej = bounds[ct]
                        bt.append(pae[si:ei, sj:ej].mean())
                        bt.append(pae[sj:ej, si:ei].mean())
                row["pae_binder_target"] = round(float(np.mean(bt)), 3)

        # ---- assembly sanity: does the C3 trimer stay together? ----
        ABC = np.vstack([pred[c] for c in BINDER_CHAINS])
        row["rg_binder_pred"] = round(
            float(np.sqrt(((ABC - ABC.mean(0)) ** 2).sum(1).mean())), 2)
        dABC = np.vstack([des[c] for c in BINDER_CHAINS])
        row["rg_binder_design"] = round(
            float(np.sqrt(((dABC - dABC.mean(0)) ** 2).sum(1).mean())), 2)
        dmin = []
        for i, ci in enumerate(BINDER_CHAINS):
            for cj2 in BINDER_CHAINS[i + 1:]:
                D = np.linalg.norm(pred[ci][:, None, :] - pred[cj2][None, :, :], axis=-1)
                dmin.append(D.min())
        row["ca_min_interchain_pred"] = round(float(np.min(dmin)), 2)

        # ---- self-consistency RMSD, chains named explicitly ----
        # (a) per-protomer: pred chain X CA vs design chain X CA (147 pairs)
        for c in BINDER_CHAINS:
            if nres[c] == 147:
                row[f"rmsd_mono_{c}"] = round(kabsch_rmsd(pred[c], des[c]), 3)
                row[f"tm_mono_{c}"] = round(tm_score(pred[c], des[c]), 3)
            else:
                row[f"rmsd_mono_{c}"] = None
        mono = [row[f"rmsd_mono_{c}"] for c in BINDER_CHAINS
                if row.get(f"rmsd_mono_{c}") is not None]
        row["rmsd_mono_mean"] = round(float(np.mean(mono)), 3) if mono else None
        row["rmsd_mono_min"] = round(float(np.min(mono)), 3) if mono else None
        tms = [row[f"tm_mono_{c}"] for c in BINDER_CHAINS if row.get(f"tm_mono_{c}") is not None]
        row["tm_mono_mean"] = round(float(np.mean(tms)), 3) if tms else None
        row["tm_mono_max"] = round(float(np.max(tms)), 3) if tms else None

        # (b) C3 trimer assembly: all of A+B+C at once, one superposition
        if all(nres[c] == 147 for c in BINDER_CHAINS):
            P = np.vstack([pred[c] for c in BINDER_CHAINS])
            Q = np.vstack([des[c] for c in BINDER_CHAINS])
            row["rmsd_trimer_ABC"] = round(kabsch_rmsd(P, Q), 3)
        else:
            row["rmsd_trimer_ABC"] = None

        # (c) bound only: full 6-chain complex, binder A,B,C + target D,E,F->a,b,c
        if mode == "bound" and all(nres[c] == 147 for c in BINDER_CHAINS):
            try:
                P = np.vstack([pred[c] for c in BINDER_CHAINS]
                              + [pred[c] for c in TARGET_MAP])
                Q = np.vstack([des[c] for c in BINDER_CHAINS]
                              + [des[TARGET_MAP[c]] for c in TARGET_MAP])
                row["rmsd_complex_6ch"] = round(kabsch_rmsd(P, Q), 3)
                # target-aligned binder RMSD: superpose on target a,b,c, then
                # measure binder A,B,C without re-fitting -- the pose metric.
                Pt = np.vstack([pred[c] for c in TARGET_MAP])
                Qt = np.vstack([des[TARGET_MAP[c]] for c in TARGET_MAP])
                Ptc, Qtc = Pt - Pt.mean(0), Qt - Qt.mean(0)
                V, S, Wt = np.linalg.svd(Ptc.T @ Qtc)
                d = np.sign(np.linalg.det(V @ Wt))
                R = V @ np.diag([1.0, 1.0, d]) @ Wt
                Pb = (np.vstack([pred[c] for c in BINDER_CHAINS]) - Pt.mean(0)) @ R
                Qb = np.vstack([des[c] for c in BINDER_CHAINS]) - Qt.mean(0)
                row["rmsd_binder_on_target_frame"] = round(
                    float(np.sqrt(((Pb - Qb) ** 2).sum() / len(Pb))), 3)
            except Exception as e:  # noqa: BLE001
                row["rmsd_complex_6ch"] = f"ERR:{e}"
        rows.append(row)

    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(k) is None else str(r.get(k, "")) for k in keys) + "\n")
    print(f"wrote {out_tsv}  n={len(rows)}")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
