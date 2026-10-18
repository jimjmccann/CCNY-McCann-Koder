#!/usr/bin/env python3
"""J497 -- score Boltz-2 BOUND folds of the C1 `c1core_v46` SolubleMPNN designs.

WHY THIS IS NOT j465. `j465_boltz_selfconsistency_20261017.py` carries its layout as
MODULE CONSTANTS, not flags: `BINDER_CHAINS = ["A","B","C"]` (:25) and
`TARGET_MAP = {"D":"a","E":"b","F":"c"}` (:26). That is the **C3** layout, where the
binder is a homotrimer. C1 `v46` is binder **A** alone + target **B,C,D**, so j465 run
unchanged here would fold two TARGET protomers into its binder average and look for a
design chain `a` that does not exist. It also pins the binder at 147 aa, while the v46
cohort spans 187-204 aa in 9 distinct lengths. ⇒ layout and length are ARGUMENTS here,
echoed into the output, and asserted.

THE MATH IS NOT COPIED. `kabsch_rmsd` and `tm_score` are IMPORTED from j465 under a
SHA1 pin, so the C1 numbers stay directly comparable to the C3 ones and cannot drift
from them silently. If j465 changes, this script REFUSES rather than quietly scoring
with different math (--allow-sha-drift to override deliberately).

MEASURES, per prediction:
  - per-chain and complex pLDDT (npz sliced by chain in file order, asserted exact)
  - ptm / iptm / confidence_score / complex_iplddt (confidence json)
  - mean PAE: binder-intra, and binder-to-target PER TARGET CHAIN (stratified, never
    pooled -- a pooled interface number is confounded, J476 / the footprint memory)
  - self-consistency vs the RFd3 design backbone:
      rmsd_binder / tm_binder          -- is it the SAME FOLD (binder refit alone)
      rmsd_binder_on_target_frame      -- is it in the SAME PLACE (superpose on the
                                          target, then measure the binder WITHOUT
                                          refitting). ⭐ THIS is the pose metric.
      rmsd_complex_4ch                 -- whole assembly, one superposition
  - interface geometry, PER TARGET CHAIN: min CA-CA distance, and CA contacts <8 A
  - ⭐⭐ ARG-GLU107 BRIDGE SURVIVAL, in the prediction AND in the design, per target
    chain: min distance from any binder ARG guanidinium N (NE/NH1/NH2) to that chain's
    GLU107 carboxylate O (OE1/OE2). Affinity lives in these bridges and nowhere else
    (the His is a destabilising push, not a binder), so a fold that keeps iptm but
    loses the bridge has lost the thing we are designing for.
    ⛔ HEAVY ATOMS, never CA: the pins are whole residues and every atom is held.

⛔ WHAT THIS DOES NOT MEASURE: affinity. iptm and pLDDT are Boltz's confidence in its
own prediction; the bridge columns are geometry, not energy. No column here is a Kd.

TARGET RESIDUE ADDRESSING IS BY SEQUENTIAL INDEX, NOT RESNUM, because Boltz renumbers
the target 1..152 while the design may carry TNFa numbering. MEASURED 2026-10-18 on
stage1 output: pred resnum = TNFa resnum - 5, verified at THREE points (K11->6,
E107->102, E110->105). The sequential index is renumbering-proof, and the script
ASSERTS the residue it lands on is actually GLU in both files before using it.

Usage:
  j497_c1_boltz_selfconsistency_20261018.py \
      --pred-dir <boltz_results_*/predictions> \
      --design-dir <dir holding the RFd3 design .pdb per backbone> \
      --out <tsv> \
      [--binder-chain A] [--target-chains B,C,D] [--mode bound|apo]
      [--design-binder-chain A] [--design-target-chains B,C,D]
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

J465 = Path(__file__).resolve().parent / "j465_boltz_selfconsistency_20261017.py"
# ⚠️ THIS PIN WAS RE-POINTED FOR THE PUBLIC DEPOSIT, and the reason matters.
#    The pin exists so the C1 and C3 numbers cannot drift apart silently. One line of j465
#    changed on deposit -- an absolute input path became a placeholder -- which changed its
#    SHA1 while leaving `kabsch_rmsd` and `tm_score` byte-identical. The pin is therefore
#    re-pointed at the DEPOSITED j465, not loosened: it still refuses if the math moves.
#    The hash the published numbers were produced under was
#    e2c7c22c01b73eb293f87088990a43a6a3a45d4c (the pre-deposit j465).
J465_SHA1 = "6e4673f2ca6e32305ec77c8ca3c7593d5aac2536"

# TNFa GLU107 / GLU110 as a 1-BASED SEQUENTIAL INDEX inside a 152-residue target chain
# (contig A6-157 => index = resnum - 5). Asserted to be GLU at load, never trusted.
GLU107_IDX = 102
GLU110_IDX = 105
ARG_TIP = (" NE ", " NH1", " NH2")
GLU_TIP = (" OE1", " OE2")


def load_j465_math(allow_drift=False):
    """Import kabsch_rmsd/tm_score from j465 under a SHA1 pin. No copy, no drift."""
    if not J465.exists():
        sys.exit(f"FATAL: j465 not found at {J465} -- the math is imported, not copied.")
    got = hashlib.sha1(J465.read_bytes()).hexdigest()
    if got != J465_SHA1:
        msg = (f"FATAL: j465 SHA1 drifted.\n  expected {J465_SHA1}\n  got      {got}\n"
               "  The C1 and C3 scores are compared to each other, so the shared math is\n"
               "  pinned. Re-verify j465, then update J465_SHA1 in this file deliberately.")
        if not allow_drift:
            sys.exit(msg)
        print("⚠️ " + msg.replace("FATAL", "WARNING") + "\n  --allow-sha-drift given: continuing.")
    spec = importlib.util.spec_from_file_location("j465_math", J465)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # safe: j465 guards its main() behind __main__
    return mod.kabsch_rmsd, mod.tm_score, got


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


MOTIF_COL = "_atom_site.is_motif_atom_with_fixed_seq"


def read_atoms_cif(path, chains):
    """mmCIF reader for RFd3 output. Columns are read BY DECLARED NAME, never by
    position -- the loop header is explicit, so hardcoded field indices are an
    avoidable way to be silently wrong.

    ⭐ Also captures RFd3's non-standard `is_motif_atom_with_fixed_seq` column, which
    is where the GUIDEPOSTS live: there is no REMARK/HETATM/segID record of them. On a
    v46 design, MEASURED: the binder has 43 True atoms = HIS(10) + 3xARG(11) = the 4 pins,
    while all three target chains are entirely True.
    """
    out = {c: {"ca": [], "res": [], "atoms": [], "motif": set()} for c in chains}
    cols, in_loop = [], False
    with _open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split()[0])
                in_loop = True
                continue
            if not (s.startswith("ATOM") or s.startswith("HETATM")):
                continue
            if not in_loop or not cols:
                raise ValueError(f"{path}: ATOM records before the _atom_site loop header")
            f = s.split()
            if len(f) != len(cols):
                raise ValueError(f"{path}: {len(f)} fields vs {len(cols)} declared columns")
            g = dict(zip(cols, f))
            ch = g.get("_atom_site.auth_asym_id", g.get("_atom_site.label_asym_id"))
            if ch not in out:
                continue
            aname = g["_atom_site.label_atom_id"]
            resname = g.get("_atom_site.auth_comp_id", g["_atom_site.label_comp_id"])
            try:
                resnum = int(g.get("_atom_site.auth_seq_id", g["_atom_site.label_seq_id"]))
                xyz = [float(g["_atom_site.Cartn_x"]), float(g["_atom_site.Cartn_y"]),
                       float(g["_atom_site.Cartn_z"])]
            except (ValueError, KeyError):
                continue
            # normalise the atom name to the 4-col PDB spelling the rest of this file uses
            pad = f" {aname:<3}" if len(aname) < 4 else aname
            out[ch]["atoms"].append((resname, resnum, pad, xyz))
            if g.get(MOTIF_COL) == "True":
                out[ch]["motif"].add(resnum)
            if aname == "CA":
                out[ch]["ca"].append(xyz)
                out[ch]["res"].append((resname, resnum))
    for c in chains:
        out[c]["ca"] = np.asarray(out[c]["ca"], dtype=float)
    return out


def read_atoms(path, chains):
    """{chain: {"ca": (N,3), "res": [(resname, resnum)], "atoms": [...], "motif": set}}

    CA list is in FILE ORDER, which is what the pLDDT slicing depends on.
    Dispatches to the mmCIF reader for .cif/.cif.gz -- the RFd3 designs are CIF, and a
    PDB column parser reads a CIF as ZERO atoms, which would look like an empty file.
    """
    low = str(path).lower()
    if low.endswith(".cif") or low.endswith(".cif.gz"):
        return read_atoms_cif(path, chains)
    out = {c: {"ca": [], "res": [], "atoms": [], "motif": set()} for c in chains}
    with _open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            ch = line[21]
            if ch not in out:
                continue
            aname = line[12:16]
            resname = line[17:20].strip()
            try:
                resnum = int(line[22:26])
                xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError:
                continue
            out[ch]["atoms"].append((resname, resnum, aname, xyz))
            if aname == " CA ":
                out[ch]["ca"].append(xyz)
                out[ch]["res"].append((resname, resnum))
    for c in chains:
        out[c]["ca"] = np.asarray(out[c]["ca"], dtype=float)
    return out


def glu_tip_coords(chain_rec, idx, label, src):
    """Carboxylate O coords of the residue at 1-BASED SEQUENTIAL index `idx`.

    ⛔ Asserts the residue is GLU. A silent miss here would turn the bridge column
    into a distance to an arbitrary residue, which would read exactly like a result.
    """
    res = chain_rec["res"]
    if len(res) < idx:
        return None, f"chain has {len(res)} residues, need index {idx}"
    resname, resnum = res[idx - 1]
    if resname != "GLU":
        return None, f"{label} index {idx} is {resname}{resnum}, not GLU ({src})"
    want = resnum
    co = [a[3] for a in chain_rec["atoms"] if a[1] == want and a[2] in GLU_TIP]
    if not co:
        return None, f"{label} GLU{resnum} has no OE1/OE2 ({src})"
    return np.asarray(co, dtype=float), None


def arg_tip_coords(chain_rec, only_resnums=None):
    """Guanidinium N coords of every ARG in a chain, with their resnums.

    `only_resnums` restricts to specific residues -- used to ask the sharper question
    "does the PINNED guidepost arginine reach GLU107", as opposed to "does ANY
    arginine happen to", which an MPNN-designed surface can satisfy by accident.
    """
    pts, who = [], []
    for resname, resnum, aname, xyz in chain_rec["atoms"]:
        if resname != "ARG" or aname not in ARG_TIP:
            continue
        if only_resnums is not None and resnum not in only_resnums:
            continue
        pts.append(xyz)
        who.append(resnum)
    return (np.asarray(pts, dtype=float) if pts else None), who


def min_dist(A, B):
    if A is None or B is None or len(A) == 0 or len(B) == 0:
        return None
    return float(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1).min())


def bridge_block(rec, binder_chain, target_chains, src, note, only_resnums=None,
                 prefix=""):
    """min ARG-NH -> GLU107-OE per target chain, plus which ARG and GLU110 as a control."""
    out = {}
    argp, argwho = arg_tip_coords(rec[binder_chain], only_resnums)
    if argp is None:
        which = "PINNED ARG" if only_resnums is not None else "ARG"
        note.append(f"{src}: NO {which} in binder chain {binder_chain}"
                    + (f" at {sorted(only_resnums)}" if only_resnums else ""))
        return out
    for tc in target_chains:
        for idx, tag in ((GLU107_IDX, "glu107"), (GLU110_IDX, "glu110")):
            co, err = glu_tip_coords(rec[tc], idx, tag, src)
            if err:
                note.append(err)
                continue
            d = np.linalg.norm(argp[:, None, :] - co[None, :, :], axis=-1)
            out[f"{prefix}{tag}_{tc}"] = round(float(d.min()), 3)
            if tag == "glu107":
                out[f"{prefix}glu107_{tc}_arg"] = argwho[
                    int(np.unravel_index(d.argmin(), d.shape)[0])]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", required=True)
    ap.add_argument("--design-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", default="bound", choices=["bound", "apo"])
    ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--target-chains", default="B,C,D")
    ap.add_argument("--design-binder-chain", default=None,
                    help="defaults to --binder-chain")
    ap.add_argument("--design-target-chains", default=None,
                    help="defaults to --target-chains")
    ap.add_argument("--contact-cut", type=float, default=8.0)
    ap.add_argument("--allow-sha-drift", action="store_true")
    a = ap.parse_args()

    kabsch_rmsd, tm_score, sha = load_j465_math(a.allow_sha_drift)

    bchain = a.binder_chain
    tchains = [c for c in a.target_chains.split(",") if c]
    dbchain = a.design_binder_chain or bchain
    dtchains = ([c for c in a.design_target_chains.split(",") if c]
                if a.design_target_chains else list(tchains))
    assert len(dtchains) == len(tchains), (dtchains, tchains)
    pred_chains = [bchain] + (tchains if a.mode == "bound" else [])

    # ⭐ ECHO THE RESOLVED LAYOUT. Never inherit a default silently.
    print(f"RESOLVED LAYOUT  mode={a.mode}")
    print(f"  prediction: binder={bchain}  target={','.join(tchains) if a.mode=='bound' else '(apo)'}")
    print(f"  design:     binder={dbchain}  target={','.join(dtchains)}")
    print(f"  j465 math imported, sha1={sha}")
    print(f"  GLU107 at sequential index {GLU107_IDX}, GLU110 at {GLU110_IDX} (asserted GLU)")

    design_dir = Path(a.design_dir)
    cache = {}

    def find_design(tag):
        if tag in cache:
            return cache[tag]
        hits = []
        for ext in ("pdb", "cif", "pdb.gz", "cif.gz"):
            hits += sorted(design_dir.rglob(f"{tag}.{ext}"))
        if not hits:
            hits = sorted(p for p in design_dir.rglob(f"*{tag}*")
                          if p.is_file() and p.suffix in (".pdb", ".cif", ".gz"))
        cache[tag] = hits[0] if hits else None
        if len(hits) > 1:
            print(f"⚠️ {tag}: {len(hits)} design files match, using {hits[0].name}")
        return cache[tag]

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

        # tag = the design backbone this sequence came from: strip _seqNN[_bound]
        tag = re.sub(r"_seq\d+(_bound|_apo)?$", "", name)
        dpath = find_design(tag)

        pred = read_atoms(pdbs[0], pred_chains)
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
            "binder_len_pred": nres[bchain],
            "n_res_total": int(off),
            "confidence_score": round(float(cj["confidence_score"]), 4),
            "ptm": round(float(cj["ptm"]), 4),
            "iptm": round(float(cj["iptm"]), 4),
            "complex_plddt": round(float(cj["complex_plddt"]), 4),
            "complex_iplddt": round(float(cj["complex_iplddt"]), 4),
            "plddt_binder": round(per_chain[bchain], 4),
        }
        for c in pred_chains:
            row[f"plddt_{c}"] = round(per_chain[c], 4)

        # ---- ⭐⭐ DECONFOUND iptm. THIS IS THE POINT OF THIS SCRIPT. ----
        # The GLOBAL iptm averages over ALL chain pairs. With 1 binder + 3 target
        # protomers, 3 of the 6 pairs are TARGET-TARGET -- the real TNFa trimer, which
        # Boltz predicts near-perfectly off a 2198-row MSA. MEASURED 2026-10-18: global
        # iptm 0.72 while binder<->target pairs were 0.30-0.45 and target<->target were
        # 0.96-0.98. ⇒ ⛔ THE GLOBAL iptm ON A TRIMERIC TARGET IS NOT AN INTERFACE
        # SCORE FOR OUR BINDER. Chain order in these dicts is INPUT order = pred_chains.
        cp = cj.get("chains_ptm") or {}
        pci = cj.get("pair_chains_iptm") or {}
        if cp and len(cp) == len(pred_chains):
            idx = {c: str(i) for i, c in enumerate(pred_chains)}
            row["ptm_binder_chain"] = round(float(cp[idx[bchain]]), 4)
            if a.mode == "bound":
                bt, tt = [], []
                for tc in tchains:
                    v = float(np.mean([float(pci[idx[bchain]][idx[tc]]),
                                       float(pci[idx[tc]][idx[bchain]])]))
                    row[f"iptm_binder_{tc}"] = round(v, 4)
                    bt.append(v)
                    row[f"ptm_{tc}"] = round(float(cp[idx[tc]]), 4)
                for i, ti in enumerate(tchains):
                    for tj in tchains[i + 1:]:
                        tt.append(float(pci[idx[ti]][idx[tj]]))
                row["iptm_binder_target_mean"] = round(float(np.mean(bt)), 4)
                row["iptm_binder_target_best"] = round(float(np.max(bt)), 4)
                row["iptm_target_target_mean"] = round(float(np.mean(tt)), 4) if tt else None
                # how much of the global iptm is the target's own trimer, not us
                row["iptm_confound_gap"] = round(
                    float(np.mean(tt)) - float(np.mean(bt)), 4) if tt else None
        elif cp:
            notes.append(f"{name}: chains_ptm has {len(cp)} entries, {len(pred_chains)} "
                         "chains scored -- per-pair iptm SKIPPED rather than mis-indexed")

        # ---- PAE: binder-intra, and binder->target PER CHAIN (never pooled) ----
        if pae is not None:
            sb, eb = bounds[bchain]
            row["pae_binder_intra"] = round(float(pae[sb:eb, sb:eb].mean()), 3)
            if a.mode == "bound":
                per_t = []
                for tc in tchains:
                    st, et = bounds[tc]
                    v = float(np.mean([pae[sb:eb, st:et].mean(), pae[st:et, sb:eb].mean()]))
                    row[f"pae_binder_{tc}"] = round(v, 3)
                    per_t.append(v)
                row["pae_binder_target_best"] = round(float(np.min(per_t)), 3)
                row["pae_binder_target_mean"] = round(float(np.mean(per_t)), 3)

        # ---- binder compactness (chain A alone; C1 has no trimer to hold together) ----
        B = pred[bchain]["ca"]
        row["rg_binder_pred"] = round(
            float(np.sqrt(((B - B.mean(0)) ** 2).sum(1).mean())), 2)

        # ---- interface geometry, PER TARGET CHAIN ----
        if a.mode == "bound":
            mins, ncon = [], []
            for tc in tchains:
                T = pred[tc]["ca"]
                D = np.linalg.norm(B[:, None, :] - T[None, :, :], axis=-1)
                row[f"ca_min_{tc}"] = round(float(D.min()), 2)
                row[f"ncontact_{tc}"] = int((D < a.contact_cut).sum())
                mins.append(float(D.min()))
                ncon.append(int((D < a.contact_cut).sum()))
            row["ca_min_target"] = round(float(np.min(mins)), 2)
            row["ncontact_total"] = int(np.sum(ncon))
            row["n_target_chains_engaged"] = int(sum(1 for m in mins if m < a.contact_cut))

        # ---- ARG -> GLU107 bridges in the PREDICTION ----
        if a.mode == "bound":
            for k, v in bridge_block(pred, bchain, tchains, f"pred:{name}", notes).items():
                row[f"pred_{k}"] = v

        # ---- self-consistency vs the design backbone ----
        if dpath is None:
            row["sc_status"] = "NO_DESIGN_FILE"
            rows.append(row)
            continue
        des = read_atoms(dpath, [dbchain] + dtchains)
        dn = len(des[dbchain]["ca"])
        row["binder_len_design"] = dn
        if dn != nres[bchain]:
            # ⛔ a length mismatch makes every paired metric meaningless -- refuse it
            row["sc_status"] = f"LEN_MISMATCH pred={nres[bchain]} design={dn}"
            rows.append(row)
            continue
        row["sc_status"] = "OK"

        row["rmsd_binder"] = round(kabsch_rmsd(B, des[dbchain]["ca"]), 3)
        row["tm_binder"] = round(tm_score(B, des[dbchain]["ca"]), 3)
        dB = des[dbchain]["ca"]
        row["rg_binder_design"] = round(
            float(np.sqrt(((dB - dB.mean(0)) ** 2).sum(1).mean())), 2)

        if a.mode == "bound":
            dt_ok = all(len(des[dc]["ca"]) == nres[tc]
                        for dc, tc in zip(dtchains, tchains))
            if not dt_ok:
                row["pose_status"] = "TARGET_LEN_MISMATCH " + ",".join(
                    f"{dc}:{len(des[dc]['ca'])}vs{tc}:{nres[tc]}"
                    for dc, tc in zip(dtchains, tchains))
            else:
                row["pose_status"] = "OK"
                P = np.vstack([B] + [pred[tc]["ca"] for tc in tchains])
                Q = np.vstack([dB] + [des[dc]["ca"] for dc in dtchains])
                row["rmsd_complex_4ch"] = round(kabsch_rmsd(P, Q), 3)

                # ⭐⭐ THE POSE METRIC: superpose on the TARGET only, then measure the
                # binder with NO refit. Large here + small rmsd_binder = right fold,
                # WRONG PLACE.
                # ⛔⛔ THE TARGET IS A C3 TRIMER, SO THERE ARE THREE EQUIVALENT SITES.
                # Under the identity chain mapping a binder docked at a symmetry-
                # equivalent site reads as a ~40 A displacement while being the SAME
                # binding mode. So try all THREE cyclic mappings of the target chains
                # and keep the BEST; `pose_frame` records which one won. Scoring only
                # the identity mapping would manufacture a failure.
                # ⛔ THE POSE NUMBER IS ONLY MEANINGFUL IF THE FRAME IS GOOD. A large
                # binder displacement measured in a badly-superposed target frame is an
                # artefact of the frame, not a finding, so report the target-only RMSD
                # alongside it: this must be small (the target is the same trimer in
                # both files) or the pose column must not be read.
                best, bestk = None, None
                for k in range(len(tchains)):
                    perm = tchains[k:] + tchains[:k]
                    Pt = np.vstack([pred[tc]["ca"] for tc in perm])
                    Qt = np.vstack([des[dc]["ca"] for dc in dtchains])
                    if Pt.shape != Qt.shape:
                        continue
                    row[f"rmsd_target_map{k}"] = round(kabsch_rmsd(Pt, Qt), 3)
                    V, S, Wt = np.linalg.svd((Pt - Pt.mean(0)).T @ (Qt - Qt.mean(0)))
                    d = np.sign(np.linalg.det(V @ Wt))
                    R = V @ np.diag([1.0, 1.0, d]) @ Wt
                    Pb = (B - Pt.mean(0)) @ R
                    Qb = dB - Qt.mean(0)
                    v = float(np.sqrt(((Pb - Qb) ** 2).sum() / len(Pb)))
                    row[f"pose_rmsd_map{k}"] = round(v, 3)
                    if best is None or v < best:
                        best, bestk = v, k
                row["rmsd_binder_on_target_frame"] = round(best, 3) if best is not None else None
                row["pose_frame"] = ("".join(tchains[bestk:] + tchains[:bestk])
                                     if bestk is not None else None)
                if bestk is not None:
                    row["rmsd_target_frame"] = row.get(f"rmsd_target_map{bestk}")
                    if row["rmsd_target_frame"] is not None and row["rmsd_target_frame"] > 5.0:
                        row["pose_status"] = (
                            f"FRAME_SUSPECT target_rmsd={row['rmsd_target_frame']}")
                        notes.append(
                            f"pred:{name}: target superposes at "
                            f"{row['rmsd_target_frame']} A -- the pose column is NOT "
                            "trustworthy for this row")

            # ---- the SAME bridges in the DESIGN, so survival is a comparison ----
            for k, v in bridge_block(des, dbchain, dtchains,
                                     f"design:{dpath.name}", notes).items():
                row[f"des_{k}"] = v
            for tc, dc in zip(tchains, dtchains):
                p, q = row.get(f"pred_glu107_{tc}"), row.get(f"des_glu107_{dc}")
                if p is not None and q is not None:
                    row[f"dbridge_{tc}"] = round(p - q, 3)

            # ---- ⭐⭐ THE SHARP QUESTION: the PINNED guidepost ARG specifically ----
            # "any ARG reaches GLU107" can be satisfied by accident on an MPNN-designed
            # surface. The guideposts are the residues RFd3 was told to place, recorded
            # ONLY in the cif's is_motif_atom_with_fixed_seq column, so read them from
            # there and never infer them.
            pins = des[dbchain].get("motif") or set()
            pin_res = {rn: nm for nm, rn in des[dbchain]["res"] if rn in pins}
            pin_args = {rn for rn, nm in pin_res.items() if nm == "ARG"}
            row["n_pins"] = len(pin_res)
            row["pins"] = ",".join(f"{pin_res[r]}{r}" for r in sorted(pin_res))
            if not pin_args:
                notes.append(f"design:{dpath.name}: no PINNED ARG among pins "
                             f"{row['pins']} -- pin-restricted bridge SKIPPED")
            else:
                for k, v in bridge_block(des, dbchain, dtchains,
                                         f"design-pin:{dpath.name}", notes,
                                         only_resnums=pin_args, prefix="pin_").items():
                    row[f"des_{k}"] = v
                # the SAME residue numbers in the PREDICTION. ⛔ assert the identity
                # survived MPNN (pinned positions are fixed, so it must be ARG there).
                pred_names = {rn: nm for nm, rn in pred[bchain]["res"]}
                bad = {r: pred_names.get(r) for r in sorted(pin_args)
                       if pred_names.get(r) != "ARG"}
                if bad:
                    notes.append(f"pred:{name}: pinned ARG resnum(s) are not ARG in the "
                                 f"prediction: {bad} -- pin bridge SKIPPED, numbering "
                                 "may not correspond")
                else:
                    for k, v in bridge_block(pred, bchain, tchains, f"pred-pin:{name}",
                                             notes, only_resnums=pin_args,
                                             prefix="pin_").items():
                        row[f"pred_{k}"] = v
                    for tc, dc in zip(tchains, dtchains):
                        p = row.get(f"pred_pin_glu107_{tc}")
                        q = row.get(f"des_pin_glu107_{dc}")
                        if p is not None and q is not None:
                            row[f"dbridge_pin_{tc}"] = round(p - q, 3)
        rows.append(row)

    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(a.out, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(k) is None else str(r.get(k, ""))
                               for k in keys) + "\n")

    ok = [r for r in rows if r.get("status") == "OK"]
    print(f"\nwrote {a.out}  n={len(rows)}  OK={len(ok)}")
    if notes:
        print(f"\n⚠️ {len(notes)} assertion note(s) -- these are NOT silently dropped:")
        for n in sorted(set(notes)):
            print("   " + n)
    for r in rows:
        if r.get("status") != "OK":
            print(f"   ⛔ {r['design']}: {r.get('status')}")
        elif r.get("sc_status") not in (None, "OK"):
            print(f"   ⚠️ {r['design']}: sc_status={r.get('sc_status')}")


if __name__ == "__main__":
    main()
