#!/usr/bin/env python3
"""POST-FOLD DESIGN QC on the submission candidates: backbone geometry, buried unsatisfied charges,
and solvent-exposed hydrophobics — the axis that bears on EXPRESSION, which is assayed before
affinity.

NEW FILE 2026-10-06. Local, CPU only, no network, no GPU. Writes only its own --out.
⛔ Interpreter: it must be a python that carries **freesasa**; biotite alone is not enough.
   ⛔ An environment whose site-packages were built against a different python minor version will
   import and then mis-resolve freesasa — check `freesasa.__version__` before trusting a run.

WHY: the intended post-fold QC axis is backbone phi/psi/tau + MolProbity, buried unsatisfied
charges, and solvent-exposed hydrophobics (freesasa RSA), post-fold only. It had never been run
on any candidate. Every ranking in this programme scores CONTACT and FOLD
SELF-CONSISTENCY. Neither notices a design that will not express.

⛔⛔ WHAT THIS IS **NOT**, STATED BECAUSE THE INTENT SAID "MolProbity"
  * **NOT a MolProbity clashscore.** That is defined on `probe` with hydrogens; `probe` is NOT
    available here (`reduce` 4.10 is, but reduce alone does not give the statistic). What is reported is
    a **heavy-atom clash count** with its own definition, below. ⛔ Do not quote it as a clashscore.
  * **NOT a MolProbity Ramachandran test.** That needs the MolProbity contour grids, which are not
    available here. What is reported is a **coarse favoured-region test** whose regions are written
    out explicitly below, so a reader can see exactly how generous it is.
  ⇒ Both are **relative** measures. They are only usable because of the control in the next block.

⭐⭐ THE INTERNAL CONTROL, AND IT IS WHAT MAKES THE NUMBERS READABLE
  Every model carries chain B = the **ecto621 target**, which is crystal-derived. So each metric is
  computed on the DESIGN (chain A) and on the TARGET (chain B) from the SAME file, by the SAME
  code. The target's rate is what "normal" looks like for this code at this resolution.
  ⇒ a design metric is reported as a ratio to its own file's target metric, never as an absolute
  against a literature threshold ⛔ (no calibrated threshold for these metrics
  exists for this pipeline, and inventing one is how a screen becomes a fake gate).

METRICS, per model, chain A (design) and chain B (target control)
  rama_out_frac   fraction of residues whose (phi,psi) falls outside EVERY coarse favoured region
  tau_out_frac    fraction with |tau - 111.0 deg| > 7.5 deg   (Berkholz 2009 ideal 111.0, sd ~2.8)
  clash_per_1k    heavy-atom pairs < 2.2 A, excluding bonded/1-2/1-3 neighbours and S-S, per 1000
                  heavy atoms  ⛔ a count under THIS SCRIPT'S definition, not a clashscore
  n_bur_unsat_chg charged side-chain groups with side-chain RSA < 10% AND no counter-charge heavy
                  atom within 4.0 A AND no polar N/O within 3.5 A
  n_exp_phobic    apolar side chains (A V L I M F W) with side-chain RSA > 50%
  exp_apolar_A2   total SASA of those apolar side chains
  ⚠️ RSA is computed IN THE COMPLEX (freesasa on the whole model), so a residue buried by the
  TARGET counts as buried. That is the right frame for expression of the complex but the WRONG
  frame for expression of the binder alone — ⛔ so `--apo` recomputes chain A by itself and both
  are reported. The difference is the interface.

Run: <a python that carries freesasa> scripts/post_fold_qc.py --label id85 \
        --dir <land dir> --id-filter 3200,3201,... --out <fresh>.csv
"""
import argparse, collections, csv, glob, math, os, re, sys
import numpy as np
import freesasa

BB = ("N", "CA", "C", "O", "OXT")
# coarse favoured (phi,psi) boxes in degrees -- WRITTEN OUT so the reader can see the generosity.
# alpha-R, beta/extended, alpha-L, and a wide bridge. ⛔ NOT MolProbity contours.
FAVOURED = [(-180, -20, -120, 50),      # alpha-R + bridge (phi -180..-20, psi -120..50)
            (-180, -20, 60, 180),       # beta / extended
            (20, 100, -20, 90),         # alpha-L
            (-180, -20, -180, -150)]    # psi wrap-around tail of alpha-R
CHARGED = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"), "LYS": ("NZ",),
           "ARG": ("NE", "NH1", "NH2"), "HIS": ("ND1", "NE2")}
POS = {"LYS", "ARG", "HIS"}
NEG = {"ASP", "GLU"}
APOLAR = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP"}
# Tien 2013 theoretical max side-chain+backbone ASA, used for RSA
MAXASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223,
          "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240,
          "PRO": 159, "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}
TAU_IDEAL, TAU_TOL = 111.0, 7.5


def read(path):
    at = []
    for ln in open(path):
        if not ln.startswith("ATOM"):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[0]).upper()
        if el == "H":
            continue
        at.append(dict(name=ln[12:16].strip(), resn=ln[17:20].strip(), ch=ln[21],
                       resi=int(ln[22:26]), el=el,
                       xyz=np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return at


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def angle(a, b, c):
    v1, v2 = a - b, c - b
    cs = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, cs))))


def in_favoured(phi, psi):
    for lo1, hi1, lo2, hi2 in FAVOURED:
        if lo1 <= phi <= hi1 and lo2 <= psi <= hi2:
            return True
    return False


def chain_metrics(at, ch, rsa):
    sel = [a for a in at if a["ch"] == ch]
    byres = collections.defaultdict(dict)
    resn = {}
    for a in sel:
        byres[a["resi"]][a["name"]] = a["xyz"]
        resn[a["resi"]] = a["resn"]
    keys = sorted(byres)
    rama_tot = rama_out = tau_tot = tau_out = 0
    for i, r in enumerate(keys):
        cur = byres[r]
        if not {"N", "CA", "C"} <= set(cur):
            continue
        tau_tot += 1
        if abs(angle(cur["N"], cur["CA"], cur["C"]) - TAU_IDEAL) > TAU_TOL:
            tau_out += 1
        if i == 0 or i == len(keys) - 1:
            continue
        pv, nx = byres[keys[i - 1]], byres[keys[i + 1]]
        if "C" not in pv or "N" not in nx:
            continue
        phi = dihedral(pv["C"], cur["N"], cur["CA"], cur["C"])
        psi = dihedral(cur["N"], cur["CA"], cur["C"], nx["N"])
        rama_tot += 1
        if not in_favoured(phi, psi):
            rama_out += 1
    # heavy-atom clashes within this chain
    from scipy.spatial import cKDTree
    P = np.array([a["xyz"] for a in sel])
    tr = cKDTree(P)
    pairs = tr.query_pairs(2.2, output_type="ndarray")
    clash = 0
    for u, v in pairs:
        au, av = sel[u], sel[v]
        if abs(au["resi"] - av["resi"]) <= 1:
            continue                      # bonded / 1-2 / 1-3 across the peptide bond
        if au["el"] == "S" and av["el"] == "S":
            continue                      # disulfide
        clash += 1
    # charges and apolars, using the RSA passed in
    bur_unsat, exp_phob, exp_area = 0, 0, 0.0
    allat = at
    AP = np.array([a["xyz"] for a in allat])
    atr = cKDTree(AP)
    for r in keys:
        rn = resn[r]
        key = (ch, r)
        if key not in rsa:
            continue
        if rn in CHARGED and rsa[key] < 10.0:
            grp = [byres[r][n] for n in CHARGED[rn] if n in byres[r]]
            if not grp:
                continue
            sat = False
            for g in grp:
                for j in atr.query_ball_point(g, 4.0):
                    o = allat[j]
                    if o["resi"] == r and o["ch"] == ch:
                        continue
                    opp = (rn in POS and o["resn"] in NEG and o["name"] in CHARGED.get(o["resn"], ())) \
                        or (rn in NEG and o["resn"] in POS and o["name"] in CHARGED.get(o["resn"], ()))
                    if opp:
                        sat = True; break
                    if o["el"] in ("N", "O") and np.linalg.norm(g - o["xyz"]) <= 3.5:
                        sat = True; break
                if sat:
                    break
            if not sat:
                bur_unsat += 1
        if rn in APOLAR and rsa[key] > 50.0:
            exp_phob += 1
            exp_area += rsa[key] * MAXASA.get(rn, 200) / 100.0
    return dict(n_res=len(keys),
                rama_out_frac=(rama_out / rama_tot if rama_tot else float("nan")),
                tau_out_frac=(tau_out / tau_tot if tau_tot else float("nan")),
                clash_per_1k=1000.0 * clash / max(1, len(sel)),
                n_bur_unsat_chg=bur_unsat, n_exp_phobic=exp_phob,
                exp_apolar_A2=exp_area)


def rsa_of(path, chains=None):
    """per-(chain,resi) RSA % from freesasa on the file as written (or on a chain subset)."""
    if chains is None:
        st = freesasa.Structure(path)
    else:
        tmp = path + f".only{''.join(chains)}.pdb"
        with open(tmp, "w") as fh:
            for ln in open(path):
                if ln.startswith("ATOM") and ln[21] in chains:
                    fh.write(ln)
            fh.write("END\n")
        st = freesasa.Structure(tmp)
        os.unlink(tmp)
    res = freesasa.calc(st).residueAreas()
    out = {}
    for ch, d in res.items():
        for ri, ra in d.items():
            try:
                out[(ch, int(ri))] = float(ra.relativeTotal) * 100.0
            except (TypeError, ValueError):
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--dir", required=True, action="append")
    ap.add_argument("--id-filter", default="", help="comma list of design ids to keep")
    ap.add_argument("--binder", default="A")
    ap.add_argument("--target", default="B")
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit(f"⛔ --out exists: {a.out}. Fresh --out per experiment.")
    keep = set(x.strip() for x in a.id_filter.split(",") if x.strip())
    files = []
    for d in a.dir:
        files += sorted(glob.glob(os.path.join(d, "**", "*_model_*.pdb"), recursive=True))
    if keep:
        files = [f for f in files
                 if (m := re.search(r'_id(\d+)_model_', os.path.basename(f))) and m.group(1) in keep]
    if a.limit:
        files = files[:a.limit]
    if not files:
        sys.exit(f"⛔ no model PDBs matched (dirs={a.dir} filter={sorted(keep)})")
    print(f"{a.label}: {len(files)} models")

    rows = []
    for n, f in enumerate(files):
        at = read(f)
        chs = {x["ch"] for x in at}
        if a.binder not in chs:
            print(f"  ⛔ {os.path.basename(f)}: no chain {a.binder}, skipped"); continue
        rsa_cplx = rsa_of(f)
        rsa_apo = rsa_of(f, chains=(a.binder,))
        d = chain_metrics(at, a.binder, rsa_cplx)
        dapo = chain_metrics(at, a.binder, rsa_apo)
        row = {"label": a.label, "model": os.path.basename(f)}
        row.update({f"d_{k}": v for k, v in d.items()})
        row["d_apo_n_bur_unsat_chg"] = dapo["n_bur_unsat_chg"]
        row["d_apo_n_exp_phobic"] = dapo["n_exp_phobic"]
        row["d_apo_exp_apolar_A2"] = dapo["exp_apolar_A2"]
        if a.target in chs:
            t = chain_metrics(at, a.target, rsa_cplx)
            row.update({f"t_{k}": v for k, v in t.items()})
        rows.append(row)
        if (n + 1) % 20 == 0:
            print(f"   {n+1}/{len(files)}")
    cols = sorted({k for r in rows for k in r}, key=lambda k: (k != "label", k != "model", k))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {a.out}  ({len(rows)} rows)")
    # summary, with the target control beside every design number
    import statistics as st
    def m(k):
        v = [r[k] for r in rows if k in r and r[k] == r[k]]
        return st.mean(v) if v else float("nan")
    print(f"\n  {'metric':24s} {'DESIGN (chain A)':>18s} {'TARGET control (B)':>20s} {'ratio':>8s}")
    for k in ("rama_out_frac", "tau_out_frac", "clash_per_1k"):
        dv, tv = m("d_" + k), m("t_" + k)
        print(f"  {k:24s} {dv:18.4f} {tv:20.4f} {(dv/tv if tv else float('nan')):8.2f}")
    for k in ("n_bur_unsat_chg", "n_exp_phobic", "exp_apolar_A2"):
        print(f"  {k:24s} {m('d_'+k):18.2f} {m('t_'+k):20.2f}")
    print(f"  {'apo n_bur_unsat_chg':24s} {m('d_apo_n_bur_unsat_chg'):18.2f}")
    print(f"  {'apo n_exp_phobic':24s} {m('d_apo_n_exp_phobic'):18.2f}")
    print(f"  {'apo exp_apolar_A2':24s} {m('d_apo_exp_apolar_A2'):18.2f}")
    print("\n  ⛔ ratios are to THIS FILE'S OWN crystal-derived target, not to a literature")
    print("     threshold. No calibrated threshold for these metrics exists in this repo.")


if __name__ == "__main__":
    main()
