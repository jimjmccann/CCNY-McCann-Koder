#!/usr/bin/env python3
"""j894 -- enumerate CORE OPTIONS for a 5XOX-scaffolded RFd3 run against TfR1.

THE CONSTRUCT THIS IMPLEMENTS
-----------------------------
Hold whole SECONDARY-STRUCTURE ELEMENTS of 5XOX chain B as a rigid motif, hold the six
functional pins, and hold NOTHING in between -- RFd3 rebuilds every linker with TfR1 present.
Produce ~10 options that differ in WHICH elements are kept. Target total length 120-160 aa,
benchmarked against 5XOX's own size. One option keeps element 10 deliberately as a control.

WHY ELEMENT 10 IS A CONTROL AND NOT A CANDIDATE
-----------------------------------------------
Under the six-pin placement, element 10 (140-169, the 30-residue helix carrying pin 155) is
337 of the 679 heavy-atom collisions with the TfR1 crop, closest approach 0.40 A. Element 12
(187-195) is another 105 with no pin. Both are excluded from the candidate pool; element 10
appears in exactly one option so the run measures what RFd3 does with a motif that starts
inside the receptor.

PLACEMENT, AND WHY IT IS RE-SOLVED PER OPTION
---------------------------------------------
The held elements keep 5XOX's NATIVE relative geometry -- that rigid arrangement is the whole
reason to use a real scaffold, because it supplies a span a short de novo chain measurably
cannot reach. But the fragment's position is set by fitting its CONTAINED pins onto the design
spec's guidepost positions, so it changes with the element set:

  >= 3 contained pins -> Kabsch, fully determined.
  == 2 contained pins -> one rotational degree of freedom survives (two points do not fix a
     rigid body). That slack is REAL and is used: sweep the rotation about the axis through the
     two pin centroids and keep the angle that minimises collision with the target.
  <= 1 contained pin  -> rejected. The placement would be arbitrary.

A pin that is NOT inside a held element stays a bare single-residue guidepost at its spec
position -- which is sound here because every pin is near-collision-free on its own (measured:
0 atoms for 32/34/38/211, 3 for 155, 1 for 180; it is the surrounding helix that collides).
A pin that IS inside a held element is defined by that element's own coordinates: the same
residue cannot be held in two places, and the resulting deviation from spec is just the
scaffold's 1.54 A constellation match quality. Both numbers are reported per option.

LINKERS
-------
Reported as a MINIMUM from the placed geometry: a helix advances ~1.5 A per residue, so the
requirement is round(CA-CA distance / 1.5) with a floor, NOT an extended-chain 3.5 A/residue.
5XOX's own path from pin 32 to pin 211 takes 178 residues to cover 24.2 A, so native loop
lengths are a gross overestimate and are NOT used as the budget. Each linker range is
[min_required, min_required + slack], and the WHOLE band must fit the length window: the option
is kept only if (held + sum of minima) >= lmin and (held + sum of maxima) <= lmax. ⛔ Not a
midpoint test -- the reported total is the low end of the band.

⛔ This script ENUMERATES and SIZES options. It does NOT write design inputs and runs no
generative model. Local CPU only; no GPU and no accelerator.
"""
import argparse
import itertools
import json
import math
import os

import gemmi
import numpy as np

# 5XOX chain B secondary-structure elements, >= 3 residues, from DSSP on the chain itself.
ELEMENTS = {
    1: ("helix", 7, 12), 2: ("strand", 23, 30), 3: ("helix", 33, 40),
    4: ("helix", 48, 64), 5: ("strand", 69, 74), 6: ("strand", 78, 83),
    7: ("helix", 89, 92), 8: ("helix", 94, 116), 9: ("strand", 130, 137),
    10: ("helix", 140, 169), 11: ("helix", 175, 182), 12: ("helix", 187, 195),
    13: ("helix", 202, 204), 14: ("helix", 207, 210), 15: ("strand", 213, 216),
    16: ("helix", 229, 232),
}
# ⛔⛔ THE PIN ORDER IS READ FROM THE HIT LIST, NEVER TYPED IN.
# The first version of this script hardcoded [32,34,38,155,180,211] -- the pins SORTED
# NUMERICALLY -- with a comment claiming it was query order. It is not. The true query order,
# which is the order of the resid_1..6 columns and of sorted(select_fixed_atoms), is
# [180,34,38,32,155,211]. The consequence: pins 34/38/211 happened to map to the right
# guidepost and 32/155/180 were silently scrambled, so every option containing one of those
# three was scored against the wrong target position. Two of the three mismatches were also
# RESIDUE-TYPE mismatches (155 is HIS onto a CATION guidepost, 180 is LYS onto an IMID one),
# which the type assertion used in j892 and j893 -- and omitted here -- would have caught on
# the first subset. assert_pin_mapping() below is that guard, restored.
PIN_LABEL = ["B1", "C1", "D1", "E1", "F1", "G1"]  # = sorted(select_fixed_atoms)
ARG_AT = ("NE", "NH1", "NH2")
HIS_AT = ("ND1", "NE2", "CG")
# ⛔⛔ GAP SIZING IS HELICAL, and that is the lab's MEASURED convention, not a choice here.
# The first version of this script used 3.5 A/residue -- essentially the EXTENDED floor -- which
# the single-residue guidepost builder names as the known-bad direction: it asks an all-helical chain to cross
# each hop in half the residues a helix needs, and a length-band assertion caught exactly that
# error once already. Our designs come out all-helical, so 1.5 A/residue is the REALISTIC size of
# a hop and 3.3 is only a floor below which the hop cannot be spanned at all.
RISE_HELIX = 1.5
RISE_EXT = 3.3
CLASH_CUT = 2.5
HEAVY_EXCLUDE = {10, 12}


def centroid(resname, at):
    if resname == "ARG":
        return np.mean([at[a] for a in ARG_AT], axis=0)
    if resname == "LYS":
        return np.asarray(at["NZ"], dtype=float)
    if resname == "HIS":
        return np.mean([at[a] for a in HIS_AT], axis=0)
    raise ValueError(resname)


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Pc, Qc


def rot_about_axis(axis, ang):
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * (K @ K)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scaffold", required=True, help="5XOX chain B, any frame")
    ap.add_argument("--chain", default="B")
    ap.add_argument("--design-pdb", required=True, help="the RFd3 design input (target + guideposts)")
    ap.add_argument("--spec", required=True)
    ap.add_argument("--hits", required=True, action="append",
                    help="j891 hit TSV(s); the pin order is READ from the scaffold's row")
    ap.add_argument("--hit-chain", required=True, help="e.g. 5xox_B")
    ap.add_argument("--lmin", type=int, default=120)
    ap.add_argument("--lmax", type=int, default=160)
    ap.add_argument("--slack", type=int, default=3)
    ap.add_argument("--rise-helix", type=float, default=RISE_HELIX)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # ---- spec guidepost positions and the target
    spec = json.load(open(a.spec))
    spec = spec[list(spec)[0]]
    sel = spec["select_fixed_atoms"]
    dst = gemmi.read_structure(a.design_pdb)
    dst.remove_alternative_conformations()
    dst.remove_hydrogens()
    qcen = []
    for key in sorted(sel):
        want = [x.strip() for x in sel[key].split(",")]
        got = None
        for ch in dst[0]:
            if ch.name != key[0]:
                continue
            for r in ch:
                at = {x.name: (x.pos.x, x.pos.y, x.pos.z) for x in r}
                if all(w in at for w in want):
                    got = np.mean([at[w] for w in want], axis=0)
        if got is None:
            raise SystemExit("guidepost %s not found" % key)
        qcen.append(got)
    qcen = np.asarray(qcen, dtype=np.float64)
    qtype = ["IMID" if "ND1" in sel[k] else "CATION" for k in sorted(sel)]
    # guidepost CA positions: where a LONE pin's single-residue guidepost actually sits in the
    # RFd3 input. ⛔ NOT the scaffold's transformed copy of that residue -- a pin outside every
    # held element is placed at the SPEC position, so a linker reaching it must be measured to
    # the spec point. The first version measured to the transformed scaffold copy, which is a
    # different place.
    qca = []
    for key in sorted(sel):
        hit = None
        for ch in dst[0]:
            if ch.name != key[0]:
                continue
            for r in ch:
                for x in r:
                    if x.name == "CA":
                        hit = np.array([x.pos.x, x.pos.y, x.pos.z])
        if hit is None:
            raise SystemExit("guidepost %s has no CA" % key)
        qca.append(hit)
    qca = np.asarray(qca, dtype=np.float64)
    targ = np.asarray([[x.pos.x, x.pos.y, x.pos.z]
                       for ch in dst[0] if ch.name == "A" for r in ch for x in r])

    # ---- pin order, READ from the hit list
    row = None
    for t in a.hits:
        for ln in list(open(t))[1:]:
            f = ln.rstrip("\n").split("\t")
            if f[1] == a.hit_chain:
                row = f
    if row is None:
        raise SystemExit("chain %r not in any --hits file" % a.hit_chain)
    PIN_RESI = [int(v) for v in row[4:10]]
    PIN_TYPE = row[10:16]
    print("[pins] query order from the hit list: %s" % list(zip(PIN_LABEL, PIN_RESI, PIN_TYPE)))

    # ---- scaffold chain
    st = gemmi.read_structure(a.scaffold)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    res = {}
    for ch in st[0]:
        if ch.name != a.chain:
            continue
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info is None or not info.is_amino_acid():
                continue
            at = {x.name: np.array([x.pos.x, x.pos.y, x.pos.z]) for x in r}
            if "CA" in at:
                res[r.seqid.num] = (r.name, at)
    if not res:
        raise SystemExit("chain %r not found in %s" % (a.chain, a.scaffold))
    pin_cen = {p: centroid(res[p][0], res[p][1]) for p in PIN_RESI}
    print("[scaffold] %s chain %s: %d residues %d-%d"
          % (os.path.basename(a.scaffold), a.chain, len(res), min(res), max(res)))

    # ---- THE GUARD THAT WAS MISSING. A pin must agree in residue TYPE with the guidepost it is
    # being fitted onto, at the SAME index. This is what catches a scrambled pin order.
    bad = []
    for i, (p, ty) in enumerate(zip(PIN_RESI, PIN_TYPE)):
        rn = res[p][0]
        want = {"CATION": {"ARG", "LYS"}, "IMID": {"HIS"}}[qtype[i]]
        if ty != qtype[i] or rn not in want:
            bad.append("%s(idx %d, %s) <- resi %d %s recorded %s"
                       % (PIN_LABEL[i], i, qtype[i], p, rn, ty))
    if bad:
        raise SystemExit("REFUSING: pin/guidepost type mismatch -- the pin order is wrong:\n  "
                         + "\n  ".join(bad))
    print("[guard] pin/guidepost type agreement: %d/%d PASS" % (len(PIN_RESI), len(PIN_RESI)))

    def elem_resi(e):
        _, lo, hi = ELEMENTS[e]
        return [n for n in range(lo, hi + 1) if n in res]

    def place(subset):
        """Return (R, Pc, Qc, pin_rmsd, nclash, dmin, contained, mode)."""
        held = sorted(set(itertools.chain.from_iterable(elem_resi(e) for e in subset)))
        contained = [p for p in PIN_RESI if p in held]
        if len(contained) < 2:
            return None
        ci = [PIN_RESI.index(p) for p in contained]
        P = np.asarray([pin_cen[p] for p in contained])
        Q = qcen[ci]
        R, Pc, Qc = kabsch(P, Q)
        heldxyz = np.vstack([np.array(list(res[n][1].values())) for n in held])

        def score(Rm):
            X = (heldxyz - Pc) @ Rm.T + Qc
            d = np.linalg.norm(X[:, None, :] - targ[None, :, :], axis=-1)
            return int((d < CLASH_CUT).sum()), float(d.min())

        if len(contained) >= 3:
            nc, dmin = score(R)
            rms = float(np.sqrt((((P - Pc) @ R.T + Qc - Q) ** 2).sum(1).mean()))
            return R, Pc, Qc, rms, nc, dmin, contained, "kabsch"
        # exactly 2 contained pins: one rotational DOF about the pin-pin axis survives
        axis = Q[1] - Q[0]
        best = None
        for deg in range(0, 360, 5):
            Rm = rot_about_axis(axis, math.radians(deg)) @ R
            nc, dmin = score(Rm)
            rms = float(np.sqrt((((P - Pc) @ Rm.T + Qc - Q) ** 2).sum(1).mean()))
            if best is None or nc < best[4]:
                best = (Rm, Pc, Qc, rms, nc, dmin, contained, "sweep%d" % deg)
        return best

    rise_helix = a.rise_helix

    def segments(subset, R, Pc, Qc):
        """Held segments in sequence order: elements, plus lone pins as 1-residue guideposts."""
        segs = [(ELEMENTS[e][1], ELEMENTS[e][2], "E%d" % e) for e in subset]
        held = set(itertools.chain.from_iterable(elem_resi(e) for e in subset))
        for p in PIN_RESI:
            if p not in held:
                segs.append((p, p, "pin%d" % p))
        segs.sort()
        def xf(v):
            return (v - Pc) @ R.T + Qc

        def anchor(resi, tag):
            """Where this segment endpoint actually sits in the RFd3 input.
            A held element is the scaffold, transformed. A LONE pin is its guidepost, at the
            SPEC position -- not the scaffold's transformed copy of that residue."""
            if tag.startswith("pin"):
                return qca[PIN_RESI.index(resi)]
            return xf(res[resi][1]["CA"])

        out, linkers = [], []
        for i, (lo, hi, tag) in enumerate(segs):
            out.append(dict(lo=lo, hi=hi, tag=tag, n=hi - lo + 1))
            if i + 1 < len(segs):
                nlo, nhi, ntag = segs[i + 1]
                d = float(np.linalg.norm(anchor(hi, tag) - anchor(nlo, ntag)))
                gmin = max(2, int(math.ceil(d / rise_helix)))
                floor = max(2, int(math.ceil(d / RISE_EXT)))
                if gmin < floor:
                    raise SystemExit("gap below the extended floor -- unspannable")
                linkers.append((d, gmin))
        return out, linkers

    # ---- candidate pool and enumeration
    pool = [e for e in ELEMENTS if e not in HEAVY_EXCLUDE]
    anchor_A = {2, 3}            # the 32/34/38 station
    anchor_B = {13, 14, 15}      # the 155/211 station
    rows = []
    for k in range(2, 8):
        for subset in itertools.combinations(pool, k):
            s = set(subset)
            if 3 not in s:                       # element 3 carries pins 34 and 38
                continue
            if not (s & anchor_B):               # must hold structure near pin 211
                continue
            pl = place(sorted(s))
            if pl is None:
                continue
            R, Pc, Qc, rms, nc, dmin, contained, mode = pl
            segs, linkers = segments(sorted(s), R, Pc, Qc)
            nheld = sum(x["n"] for x in segs)
            lo_len = nheld + sum(mn for _, mn in linkers)
            hi_len = nheld + sum(mn + a.slack for _, mn in linkers)
            # ⛔ the WHOLE band must fit the window, not a midpoint
            if not (a.lmin <= lo_len and hi_len <= a.lmax):
                continue
            nlink = lo_len - nheld
            total = lo_len
            rows.append(dict(subset=sorted(s), nheld=nheld, nlink=nlink, total=total,
                             len_min=lo_len, len_max=hi_len,
                             clash=nc, dmin=round(dmin, 2), pin_rmsd=round(rms, 3),
                             contained=contained, mode=mode, nseg=len(segs),
                             linkers=[(round(d, 1), mn) for d, mn in linkers]))
    rows.sort(key=lambda r: (r["clash"], r["pin_rmsd"], -r["nheld"]))
    print("[enum] %d element subsets land in %d-%d aa" % (len(rows), a.lmin, a.lmax))

    # ---- the element-10 control, priced the same way
    ctrl = None
    for k in range(2, 7):
        for subset in itertools.combinations([e for e in ELEMENTS if e != 12], k):
            s = set(subset)
            if 10 not in s or 3 not in s or not (s & anchor_B):
                continue
            pl = place(sorted(s))
            if pl is None:
                continue
            R, Pc, Qc, rms, nc, dmin, contained, mode = pl
            segs, linkers = segments(sorted(s), R, Pc, Qc)
            nheld = sum(x["n"] for x in segs)
            lo_len = nheld + sum(mn for _, mn in linkers)
            hi_len = nheld + sum(mn + a.slack for _, mn in linkers)
            if not (a.lmin <= lo_len and hi_len <= a.lmax):
                continue
            nlink = lo_len - nheld
            total = lo_len
            cand = dict(subset=sorted(s), nheld=nheld, nlink=nlink, total=total,
                        len_min=lo_len, len_max=hi_len, clash=nc,
                        dmin=round(dmin, 2), pin_rmsd=round(rms, 3), contained=contained,
                        mode=mode, nseg=len(segs),
                        linkers=[(round(d, 1), mn) for d, mn in linkers], control=True)
            if ctrl is None or cand["clash"] < ctrl["clash"]:
                ctrl = cand

    # ---- pick 10: spread over element composition, not just the top of the ranking
    picked, seen_sig = [], set()
    for r in rows:
        sig = (len(r["subset"]), tuple(sorted(set(r["subset"]) & {4, 8})), r["nseg"])
        if sig in seen_sig:
            continue
        seen_sig.add(sig)
        picked.append(r)
        if len(picked) >= 9:
            break
    for r in rows:
        if len(picked) >= 9:
            break
        if r not in picked:
            picked.append(r)
    if ctrl:
        ctrl["note"] = "CONTROL: holds element 10, which is inside the receptor"
        picked.append(ctrl)

    out = dict(scaffold=os.path.basename(a.scaffold), chain=a.chain,
               length_window=[a.lmin, a.lmax], elements=ELEMENTS, pins=PIN_RESI,
               pin_labels=PIN_LABEL, excluded_from_pool=sorted(HEAVY_EXCLUDE),
               n_subsets_in_window=len(rows), options=picked)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("[out] %s" % a.out)

    print("\n%-3s %-26s %5s %5s %5s %6s %6s %7s %5s %s"
          % ("#", "elements held", "held", "link", "tot", "clash", "minD", "pinRMS", "nseg", "mode"))
    for i, r in enumerate(picked, 1):
        tag = "  <- CONTROL (element 10)" if r.get("control") else ""
        print("%-3d %-26s %5d %5d %5d %6d %6.2f %7.3f %5d %s%s"
              % (i, "+".join("E%d" % e for e in r["subset"]), r["nheld"], r["nlink"],
                 r["total"], r["clash"], r["dmin"], r["pin_rmsd"], r["nseg"], r["mode"], tag))
    print("\npins held INSIDE an element, per option (the rest stay bare guideposts at spec):")
    for i, r in enumerate(picked, 1):
        print("  %-3d contained %s" % (i, r["contained"]))


main()
