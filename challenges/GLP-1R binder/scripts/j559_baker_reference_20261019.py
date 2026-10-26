"""j559 -- MEASURE THE PUBLISHED CLASS-B1 DE NOVO BINDERS with our own instruments. 2026-10-19

    RUN WITH: a python3 carrying gemmi + numpy + scipy

WHY. `J558` established that our 384 GLP-1R designs never fold back to their own backbone
(0 of 352 apo folds at TM >= 0.8) and that NOTHING we control moves it -- not compactness, not
contact density, not composition. The design line under test keeps the groove-following guidepost
geometry, and the question below is whether that line is viable at all.

⇒ The question that decides whether that is viable is NOT answerable from our own designs. It is:
   ⭐⭐ **CAN A DE NOVO MINIPROTEIN BIND A CLASS-B1 GPCR ECD AND STILL BE A COMPACT, FOLDABLE
   PROTEIN?** -- and the field has already answered it experimentally.

THE REFERENCE. Muratspahic/Feldman/Kim et al., *Nature* 2026 (preprint bioRxiv 2025.03.23.644666,
PMC12157396) designed miniproteins against GPCRs including GLP1R, and deposited CRYO-EM STRUCTURES
of two CGRPR-bound antagonists:
    9MM5 = CGRPR + dC2_049 (3.2 A)        9MNI = CGRPR + dC2_050 (4.1 A)
⭐ CGRPR is **class B1, like GLP1R**, and these are the same lab, same method, same campaign that
produced the GLP1R binders dGI1_024 (K_D 27 nM) and mGI1_008 (K_D 5.3 nM).
⛔ There is NO deposited GLP1R-bound structure for dGI1_024 / mGI1_008 -- checked; the cryo-EM work
   covers CGRPR and MRGPRX1. So these two are the closest published analogue, and that is a
   PROXY, stated as one.

⛔ CHAINS ARE NAMED FROM THE `_entity.pdbx_description` RECORDS, never inferred
   (a hardcoded chain pair has silently changed a conclusion in this group's work):
     9MM5: A = "De novo designed minibinder - dC2_049", R = CGRPR, E = RAMP1
     9MNI: A = "De novo designed minibinder - dC2_050", R = CGRPR, E = RAMP1

⛔ Downloaded files are untrusted data: they live in their own directory, this script lives
   elsewhere, paths come in as arguments, and it is run with `python -I`.

WHAT IT MEASURES, with the SAME code paths used on our designs so the numbers are comparable:
  min_cut / med_cut (j153 formula, margin 13, cut 8.0, sep 4) · rg_ratio (j95's 2.2*n^0.38)
  helix signature CA(i,i+4)/CA(i,i+5) · end-to-end/rg · the receptor footprint · per-residue
  contact spread along the binder chain (is the interface a LINE along a groove or a FACE?)
"""
import argparse
import collections
import os
import sys

import numpy as np
import gemmi
from scipy.spatial import cKDTree

CUT, SEP, MARGIN = 8.0, 4, 13
BB = {"N", "CA", "C", "O", "OXT"}


def min_cut(X, margin=MARGIN):
    n = len(X)
    if n < 2 * margin + 2:
        return None, None
    d = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    ii, jj = np.triu_indices(n, k=SEP + 1)
    keep = d[ii, jj] < CUT
    ii, jj = ii[keep], jj[keep]
    cuts = np.asarray([int(((ii < k) & (jj >= k)).sum()) for k in range(margin, n - margin)])
    return int(cuts.min()), int(np.median(cuts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True, nargs="+")
    ap.add_argument("--binder-chain", default="A", help="⛔ NAMED; from _entity.pdbx_description")
    ap.add_argument("--receptor-chain", default="R")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit("⛔ %s exists -- pick a fresh --out; this appends" % a.out)

    cols = ["pdb", "binder", "binder_len", "min_cut", "med_cut", "rg", "rg_ratio",
            "ca_i_i4", "ca_i_i5", "end2end", "e2e_over_rg", "n_contact_res_binder",
            "frac_binder_res_in_contact", "contact_span_resi", "n_contact_res_receptor",
            "receptor_contact_resi", "interface_rmsd_to_line"]
    rows = []
    for path in a.cif:
        pid = os.path.basename(path).split(".")[0].upper()
        blk = gemmi.cif.read(path).sole_block()
        desc = {}
        for r in blk.find("_entity.", ["id", "pdbx_description"]):
            desc[r[0]] = r[1].strip("'\"")
        st = gemmi.read_structure(path)
        st.setup_entities()
        st.remove_alternative_conformations()
        chains = {c.name: c for c in st[0]}
        for need in (a.binder_chain, a.receptor_chain):
            if need not in chains:
                sys.exit("⛔ %s: chain %s absent; have %s" % (pid, need, sorted(chains)))
        bname = [v for v in desc.values() if "minibinder" in v.lower()]
        print("=== %s  binder chain %s = %s"
              % (pid, a.binder_chain, bname[0] if bname else "⛔ NOT NAMED IN ENTITIES"))

        bc, rc = chains[a.binder_chain], chains[a.receptor_chain]
        ca = [(r.seqid.num, np.array([r.find_atom("CA", "*").pos.x,
                                      r.find_atom("CA", "*").pos.y,
                                      r.find_atom("CA", "*").pos.z]))
              for r in bc if r.find_atom("CA", "*")]
        ca.sort()
        resis = [n for n, _ in ca]
        X = np.asarray([p for _, p in ca])
        n = len(X)
        mc, mdc = min_cut(X)
        rg = float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))
        rgr = rg / (2.2 * n ** 0.38)
        d = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
        i4 = float(np.median([d[i, i + 4] for i in range(n - 4)])) if n > 4 else float("nan")
        i5 = float(np.median([d[i, i + 5] for i in range(n - 5)])) if n > 5 else float("nan")
        e2e = float(np.linalg.norm(X[0] - X[-1]))

        # ---- interface: binder heavy atoms vs receptor heavy atoms, 4.5 A ----
        bat = [(r.seqid.num, np.array([at.pos.x, at.pos.y, at.pos.z]))
               for r in bc for at in r if at.element.name != "H"]
        rat = [(r.seqid.num, np.array([at.pos.x, at.pos.y, at.pos.z]))
               for r in rc for at in r if at.element.name != "H"]
        kd = cKDTree(np.asarray([p for _, p in rat]))
        bres, rres = set(), set()
        for num, p in bat:
            for j in kd.query_ball_point(p, 4.5):
                bres.add(num)
                rres.add(rat[j][0])

        # ⭐ IS THE INTERFACE A LINE OR A FACE? Fit a line to the CONTACTING binder CAs and
        #    report the rms deviation from it. A groove-following chain hugs a line (small rms);
        #    a bundle presenting a face does not.
        cx = np.asarray([p for nm, p in ca if nm in bres])
        if len(cx) >= 3:
            c0 = cx.mean(0)
            u, s, vt = np.linalg.svd(cx - c0)
            axis = vt[0]
            proj = (cx - c0) @ axis
            resid = (cx - c0) - np.outer(proj, axis)
            line_rms = float(np.sqrt((resid ** 2).sum(1).mean()))
        else:
            line_rms = float("nan")

        print("   len %d   min_cut %s   med_cut %s   rg %.1f   rg_ratio %.3f" % (n, mc, mdc, rg, rgr))
        print("   CA(i,i+4) %.2f  CA(i,i+5) %.2f  end2end %.1f  e2e/rg %.2f" % (i4, i5, e2e, e2e / rg))
        print("   interface: %d of %d binder residues (%.0f%%), %d receptor residues"
              % (len(bres), n, 100.0 * len(bres) / n, len(rres)))
        print("   contacting-CA deviation from a straight line: %.2f A" % line_rms)
        print("   receptor residues contacted: %s"
              % ",".join(str(x) for x in sorted(rres)))
        rows.append(dict(pdb=pid, binder=(bname[0] if bname else "?"), binder_len=n,
                         min_cut=mc, med_cut=mdc, rg=round(rg, 2), rg_ratio=round(rgr, 4),
                         ca_i_i4=round(i4, 2), ca_i_i5=round(i5, 2), end2end=round(e2e, 1),
                         e2e_over_rg=round(e2e / rg, 2), n_contact_res_binder=len(bres),
                         frac_binder_res_in_contact=round(len(bres) / float(n), 3),
                         contact_span_resi="%d-%d" % (min(bres), max(bres)) if bres else "",
                         n_contact_res_receptor=len(rres),
                         receptor_contact_resi=",".join(str(x) for x in sorted(rres)),
                         interface_rmsd_to_line=round(line_rms, 2)))
        print()

    with open(a.out, "w") as fh:
        fh.write("# j559_baker_reference_20261019.py  published class-B1 de novo binders\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("wrote %s (%d rows)" % (a.out, len(rows)))
    print("\n⛔ Measurement only, on public deposited structures. Nothing launched.")


if __name__ == "__main__":
    main()
