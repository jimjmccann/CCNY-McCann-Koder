#!/usr/bin/env python3
"""THE FULL-RECEPTOR GATE — superpose each design onto the whole receptor and measure what the
generation target did not contain.

⛔⛔ WHY THIS EXISTS. A design campaign built against a CROP of the receptor can place designs
inside receptor that was deleted from the generation target — backbone overlap, not contact. A crop
validated by INCLUSION ("all 11 footprint residues are present") and never by EXCLUSION leaves that
failure unfalsifiable from inside the pipeline, because every downstream stage treats the crop as
the target. This script is the missing EXCLUSION check. Replacing a crop with a sphere cut is an
improvement, not a fix: the radius is a GUESS, and this is the only thing that proves the guess was
generous enough.

⛔⛔ IT MUST SUPERPOSE, AND THAT IS NOT A STYLE CHOICE. RFd3 recenters its output by **66.93 A**
(measured; the same offset on every model). A gate that diffs raw coordinates calls every design a
catastrophic clash; one that "passes" everything is equally wrong. So superpose on the residues the
design and the receptor SHARE, then measure the binder against everything else.

⭐ The superposition is well-posed because RFd3 holds the target RIGID: the sphere comes back at
0.0402 A over 978/978 atoms. A fit RMSD materially above that means the match is wrong — so the
RMSD is itself a tripwire, reported per design and asserted against --max-fit-rmsd.

WHAT IT MEASURES, per design:
  fit_rmsd      superposition RMSD on the shared CA atoms      (tripwire: should be ~0.0-0.1 A)
  d_omitted     min binder-heavy-atom distance to receptor protein NOT in the generation target
  d_glycan      min binder-heavy-atom distance to any glycan atom
  glyc_sequon   WHICH sequon that nearest sugar belongs to (attributed by nearest ASN ND2)
  d_N<nnn>      per-sequon minimum, one column per sequon within --sequon-report of either site
  d_sphere      min binder-heavy-atom distance to the residues that WERE in the target (context)
  n_omit_2A     count of omitted-receptor atoms within 2.0 A of the binder   <- the voiding metric
  n_glyc_2A     count of glycan atoms within 2.0 A of the binder
  verdict       PASS / WARN_GLYCAN / CLASH_OMITTED / CLASH_GLYCAN

⚠️⚠️ WHY THERE IS A **WARN** TIER FOR GLYCANS AND NOT JUST A CLASH. The crystal resolves only
1-7 sugars per tree; native HEK293 trees are ~15-20 residues and sialylated, and the antennae are
absent from every structure. So a design sitting 3.5 A from a TRUNCATED tree may well be inside
the real one. ⇒ a glycan distance is a FLOOR on the clash, never the clash itself, and anything
inside --warn-glycan is flagged rather than passed silently. MEASURED on this model:
    sequon  sugars modelled   -> site A   -> site B
    N420          2            7.47 A      35.91
    N328          6 (biggest)  11.70       32.54
    N32           1            26.14        3.91
⛔ TWO WAYS A GLYCAN CHECK MISSES. (1) N328 lies outside a 106 aa crop, so a check that looks only
inside the crop never sees the biggest modelled tree. (2) A check that transfers ONLY THE SUGAR onto
each design, rather than the whole receptor, cannot detect a protein-level crop failure at all.
This script carries the entire receptor, so both are closed. ⛔ And a sequon count is a property of
the region examined, not of the receptor: a crop with ZERO sequons becomes a sphere with N32 sitting
3.91 A from site B, where 3 of 23 site-B designs clash with it. That risk is invisible from inside
the crop.

Usage:
  full_receptor_gate.py <design-dir> --receptor <full.pdb> --wt2out <WT2OUT.tsv> \
      [--csv out.csv] [--write-superposed DIR] [--clash 2.0] [--max-fit-rmsd 0.5]
"""
import argparse, glob, gzip, os, sys
import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure.io.pdb import PDBFile
from biotite.structure import superimpose


def load_design(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        st = get_structure(CIFFile.read(fh), model=1)
    return st[st.element != "H"]


def load_receptor(path):
    st = PDBFile.read(path).get_structure(model=1)
    return st[st.element != "H"]


def read_map(tsv):
    """-> (wt_list, out_list, resname_list) in output order."""
    wt, out, rn = [], [], []
    for l in open(tsv):
        if l.startswith("#") or l.startswith("wt_resnum"):
            continue
        a, b, c, _ = l.split()
        wt.append(int(a)); out.append(int(b)); rn.append(c)
    return wt, out, rn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--wt2out", required=True)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--write-superposed", default=None,
                    help="dir for designs transformed into the RECEPTOR's frame")
    ap.add_argument("--clash", type=float, default=2.0,
                    help="an atom this close counts as overlap, not contact")
    ap.add_argument("--max-fit-rmsd", type=float, default=0.5,
                    help="refuse a design whose superposition is worse than this -- the match is wrong")
    ap.add_argument("--glycan-chain", default="G")
    ap.add_argument("--warn-glycan", type=float, default=4.0,
                    help="flag WARN_GLYCAN inside this distance. The modelled trees are TRUNCATED "
                         "(1-7 sugars vs ~15-20 native), so a near miss on this model is not a "
                         "near miss on the real receptor. Set 0 to disable.")
    ap.add_argument("--sequon-report", type=float, default=15.0,
                    help="emit a per-sequon distance column for every sequon whose tree lies "
                         "within this distance of ANY mapped target residue")
    a = ap.parse_args()

    wt, out, rn3 = read_map(a.wt2out)
    rec = load_receptor(a.receptor)
    prot = rec[rec.chain_id == "A"]
    glyc = rec[rec.chain_id == a.glycan_chain]
    print(f"receptor {os.path.basename(a.receptor)}: {len(np.unique(prot.res_id))} protein residues, "
          f"{len(glyc)} glycan atoms on chain {a.glycan_chain}")
    print(f"map {os.path.basename(a.wt2out)}: {len(wt)} shared residues "
          f"⇒ {len(np.unique(prot.res_id)) - len(wt)} receptor residues were OMITTED from generation")

    # receptor CA for the shared residues, in map order; assert type identity
    rec_ca, bad = [], []
    for w, r in zip(wt, rn3):
        s = prot[(prot.res_id == w) & (prot.atom_name == "CA")]
        if not len(s):
            bad.append((w, "ABSENT")); continue
        if prot[prot.res_id == w].res_name[0] != r:
            bad.append((w, prot[prot.res_id == w].res_name[0]))
        rec_ca.append(s.coord[0])
    if bad:
        sys.exit(f"⛔ {len(bad)} shared residues missing/mismatched in the receptor: {bad[:5]}")
    rec_ca = np.array(rec_ca)
    print(f"✅ all {len(rec_ca)} shared residues present in the receptor, types identical")

    # ---- attribute every glycan residue to its parent ASN, so a hit is NAMED, not anonymous ----
    asn = prot[(prot.res_name == "ASN") & (prot.atom_name == "ND2")]
    seq_of = {}
    for g in np.unique(glyc.res_id):
        gs = glyc[glyc.res_id == g]
        if not len(asn):
            continue
        d = np.linalg.norm(gs.coord[:, None] - asn.coord[None], axis=-1)
        seq_of[int(g)] = int(asn.res_id[int(np.unravel_index(d.argmin(), d.shape)[1])])
    seqs = sorted(set(seq_of.values()))
    shared_atoms = prot[np.isin(prot.res_id, wt)]
    near_seqs = []
    for n in seqs:
        gm = glyc[np.isin(glyc.res_id, [g for g, p in seq_of.items() if p == n])]
        dd = float(np.linalg.norm(shared_atoms.coord[:, None] - gm.coord[None], axis=-1).min())
        if dd <= a.sequon_report:
            near_seqs.append((n, len(np.unique(gm.res_id)), dd))
    print(f"   glycan trees: {len(seqs)} sequons, "
          f"{[f'N{n}({len(np.unique(glyc[np.isin(glyc.res_id,[g for g,p in seq_of.items() if p==n])].res_id))})' for n in seqs]}")
    if near_seqs:
        print(f"   ⚠️ within {a.sequon_report} A of this target: "
              + ", ".join(f"N{n} ({k} sugars, {d:.2f} A)" for n, k, d in near_seqs))
    else:
        print(f"   ✅ no sequon tree within {a.sequon_report} A of this target")

    omit_mask = ~np.isin(prot.res_id, wt)
    omitted = prot[omit_mask]
    print(f"   omitted receptor atoms: {len(omitted)} over "
          f"{len(np.unique(omitted.res_id))} residues\n")

    files = sorted(glob.glob(os.path.join(a.indir, "*.cif.gz")) +
                   glob.glob(os.path.join(a.indir, "*.cif")) +
                   glob.glob(os.path.join(a.indir, "*.pdb")))
    if not files:
        sys.exit(f"⛔ no designs in {a.indir}")
    if a.write_superposed:
        os.makedirs(a.write_superposed, exist_ok=True)

    rows = []
    for f in files:
        st = load_design(f) if not f.endswith(".pdb") else \
             PDBFile.read(f).get_structure(model=1)
        st = st[st.element != "H"]
        B = st[st.chain_id == "B"]
        A = st[st.chain_id == "A"]
        if not len(A) or not len(B):
            print(f"  SKIP {os.path.basename(f)}: missing chain A or B"); continue
        # design CA for the shared residues, in the SAME order as rec_ca
        dca = []
        okay = True
        for o in out:
            s = B[(B.res_id == o) & (B.atom_name == "CA")]
            if not len(s):
                okay = False; break
            dca.append(s.coord[0])
        if not okay:
            print(f"  SKIP {os.path.basename(f)}: target chain B does not carry all mapped residues")
            continue
        dca = np.array(dca)

        # Kabsch: design -> receptor frame, fitted ONLY on the shared target residues
        from biotite.structure import AtomArray
        tmp = AtomArray(len(dca)); tmp.coord = dca
        ref = AtomArray(len(rec_ca)); ref.coord = rec_ca
        _, transform = superimpose(ref, tmp)
        st_sup = st.copy()
        st_sup.coord = transform.apply(st.coord)
        fit = float(np.sqrt(np.mean(np.sum(
            (transform.apply(dca) - rec_ca) ** 2, axis=-1))))

        Asup = st_sup[st_sup.chain_id == "A"]
        def mind(target):
            if not len(target):
                return float("nan"), 0
            D = np.linalg.norm(Asup.coord[:, None, :] - target.coord[None, :, :], axis=-1)
            return float(D.min()), int((D < a.clash).sum())
        d_om, n_om = mind(omitted)
        d_gl, n_gl = mind(glyc)
        d_sp, _ = mind(prot[~omit_mask])

        # which sequon owns the nearest sugar, and the per-sequon minima
        per_seq, worst_seq = {}, None
        for n, _k, _d in near_seqs:
            gm = glyc[np.isin(glyc.res_id, [g for g, p in seq_of.items() if p == n])]
            dn, _ = mind(gm)
            per_seq[f"d_N{n}"] = round(dn, 2)
            if worst_seq is None or dn < per_seq[f"d_N{worst_seq}"]:
                worst_seq = n
        if not len(glyc):
            worst_seq = None

        verdict = "PASS"
        if fit > a.max_fit_rmsd:
            verdict = "BAD_FIT"
        elif n_om:
            verdict = "CLASH_OMITTED"
        elif n_gl:
            verdict = "CLASH_GLYCAN"
        elif a.warn_glycan and d_gl == d_gl and d_gl < a.warn_glycan:
            verdict = "WARN_GLYCAN"
        rows.append(dict(design=os.path.basename(f), fit_rmsd=round(fit, 4),
                         d_omitted=round(d_om, 2), n_omit_2A=n_om,
                         d_glycan=round(d_gl, 2), n_glyc_2A=n_gl,
                         glyc_sequon=(f"N{worst_seq}" if worst_seq is not None else ""),
                         d_sphere=round(d_sp, 2), verdict=verdict, **per_seq))
        if a.write_superposed:
            o = os.path.join(a.write_superposed,
                             os.path.basename(f).replace(".cif.gz", "").replace(".pdb", "")
                             + "_inreceptor.pdb")
            pf = PDBFile(); pf.set_structure(st_sup); pf.write(o)

    if not rows:
        sys.exit("⛔ nothing measured")
    import csv as _csv
    hdr = list(rows[0])
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=hdr); w.writeheader(); w.writerows(rows)
    print(f"  {'design':40s} {'fit':>6} {'d_omit':>7} {'n<2A':>5} {'d_glyc':>7} {'n<2A':>5} "
          f"{'sequon':>7} {'verdict':>14}")
    for r in rows:
        print(f"  {r['design'][:40]:40s} {r['fit_rmsd']:>6} {r['d_omitted']:>7} "
              f"{r['n_omit_2A']:>5} {r['d_glycan']:>7} {r['n_glyc_2A']:>5} "
              f"{r['glyc_sequon']:>7} {r['verdict']:>14}")
    fits = np.array([r["fit_rmsd"] for r in rows])
    print(f"\n  fit_rmsd: min {fits.min():.4f} median {np.median(fits):.4f} max {fits.max():.4f} "
          f"(⭐ RFd3 holds the target rigid to ~0.04 A; a big value means the match is wrong)")
    from collections import Counter
    for k, v in Counter(r["verdict"] for r in rows).items():
        print(f"  {k:16s} {v}/{len(rows)}")
    if a.csv:
        print(f"  wrote {a.csv}")


if __name__ == "__main__":
    main()
