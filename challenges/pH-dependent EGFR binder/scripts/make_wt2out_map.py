#!/usr/bin/env python3
"""Persist the WT -> RFd3-output residue map for a sphere-cut target.

⛔⛔ WHY THIS EXISTS.
RFd3 renumbers the target chain 1..N CONTIGUOUS and DISCARDS both the WT numbers and every
segment gap.  A sphere cut is multi-segment (site A 15 segs, site B 26), so **no constant
offset exists** -- the map is piecewise.  Every scorer that hardcoded one (`select_all31.py`
W2O=-369, `select_siteB.py` A2O=-3, `mpnn_egfr.py` anchor=(409-369)) reads the WRONG residue,
and site B's happened to be right only because all 8 of its footprint residues fall in the
sphere's FIRST segment.  Nothing persisted the real map: the builder computed `keep` and threw
it away, so it was recoverable only by re-running the cut at the identical radius.

⭐ VERIFIED against real RFd3 output, not assumed: output chain B residue i (1-based, file
order) IS input residue i (1-based, file order).  Checked by residue-type identity on the
site-A 22 A sphere, 231/231 -- see `--verify-cif`.
⛔ Chain roles in RFd3 output: BINDER = chain A, TARGET = chain B.  (Opposite to RFd2.)

Usage:
  make_wt2out_map.py <sphere.pdb> --out WT2OUT.tsv [--verify-cif <rfd3 out .cif.gz>]
"""
import argparse, gzip, sys


def read_sphere(pdb):
    """Residues in FILE ORDER -- that order is what RFd3 renumbers. -> [(wt, resname, chain)]."""
    out, seen = [], None
    for l in open(pdb):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        key = (l[21], int(l[22:26]), l[26])
        if key != seen:
            seen = key
            out.append((int(l[22:26]), l[17:20].strip(), l[21]))
    return out


def read_cif_chain(cif, chain="B"):
    """(auth_seq_id, auth_comp_id) of one chain of an RFd3 output cif, in file order."""
    rows, inloop, hdr = [], False, []
    op = gzip.open if cif.endswith(".gz") else open
    for l in op(cif, "rt"):
        if l.startswith("_atom_site."):
            inloop = True; hdr.append(l.strip())
        elif inloop and l.startswith(("ATOM", "HETATM")):
            rows.append(l.split())
        elif inloop and rows and (l.startswith("#") or not l.strip()):
            break
    i_seq = hdr.index("_atom_site.auth_seq_id")
    i_comp = hdr.index("_atom_site.auth_comp_id")
    i_ch = hdr.index("_atom_site.auth_asym_id")
    res, seen = [], None
    for r in rows:
        if r[i_ch] != chain:
            continue
        if r[i_seq] != seen:
            seen = r[i_seq]
            res.append((int(r[i_seq]), r[i_comp]))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--out", required=True)
    ap.add_argument("--verify-cif", help="an RFd3 output cif(.gz) built from THIS target")
    a = ap.parse_args()

    res = read_sphere(a.pdb)
    # segment id: increments wherever the WT numbering is not consecutive
    segs, s = [], 0
    for i, (wt, _, _) in enumerate(res):
        if i and wt != res[i - 1][0] + 1:
            s += 1
        segs.append(s)

    if a.verify_cif:
        cif = read_cif_chain(a.verify_cif, "B")
        if len(cif) != len(res):
            sys.exit(f"⛔ VERIFY FAILED: cif chain B has {len(cif)} res, target pdb has {len(res)}")
        bad = [(i + 1, res[i][1], c[1]) for i, c in enumerate(cif) if c[1] != res[i][1]]
        if bad:
            sys.exit(f"⛔ VERIFY FAILED: {len(bad)} residue-type mismatches, first {bad[:5]}")
        if [c[0] for c in cif] != list(range(1, len(cif) + 1)):
            sys.exit("⛔ VERIFY FAILED: cif chain B is not numbered 1..N contiguous")
        print(f"✅ VERIFIED against {a.verify_cif}: {len(cif)}/{len(cif)} residue types identical, "
              f"chain B numbered 1..{len(cif)}")

    with open(a.out, "w") as fh:
        fh.write("# WT -> RFd3 output residue map for %s\n" % a.pdb.split("/")[-1])
        fh.write("# out_resnum = 1-based rank in file order; TARGET is chain B in RFd3 output.\n")
        fh.write("# BINDER is chain A. Piecewise over %d segments: NO constant offset exists.\n"
                 % (segs[-1] + 1))
        fh.write("wt_resnum\tout_resnum\tresname\tseg\n")
        for i, (wt, rn, _) in enumerate(res):
            fh.write(f"{wt}\t{i+1}\t{rn}\t{segs[i]}\n")
    print(f"wrote {a.out}: {len(res)} residues, {segs[-1]+1} segments, "
          f"WT {res[0][0]}-{res[-1][0]} -> out 1-{len(res)}")


if __name__ == "__main__":
    main()
