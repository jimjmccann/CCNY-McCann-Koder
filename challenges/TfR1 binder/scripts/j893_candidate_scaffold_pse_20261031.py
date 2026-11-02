"""j893 -- one .pse: the four face-gate candidate scaffolds SUPERIMPOSED on the fixed-residue
constellation, with the TfR1 crop drawn so it is visible whether they could sit there at all.

RUN:  pymol -cq j893_candidate_scaffold_pse_20261031.py -- \
        --design-pdb <guidepost constellation .pdb> \
        --spec <design spec .json> \
        --hits ../data/j891_sixpin_hits_fx_20261031.tsv \
        --hits ../data/j891_sixpin_hits_afdb_20261031.tsv \
        --mirror <local mmCIF mirror> --afdb-dir <predicted-model extract> \
        --cand 5xox_B --cand 9gyt_N --cand 9t6h_b --cand AF-P56144-F1-model_v6_A \
        --workdir <scratch dir> --out <file.pse>

⛔ NO `if __name__ == "__main__"` GUARD. Under `pymol -cq` this module is NOT __main__, so a
  guarded script exits 0 having written NOTHING.

⛔⛔ EVERY `scene(message=)` AND EVERY LABEL IS PURE ASCII. A non-ASCII mark corrupts PyMOL's
  heap and the .pse SEGFAULTS ON LOAD.

HOW THE SUPERPOSITION IS DEFINED, and the gate on it
----------------------------------------------------
The transform is Kabsch on SIX POINTS ONLY: each candidate's matched functional-group centroid
onto the corresponding guidepost centroid, in the query's own pin order (B1,C1,D1,E1,F1,G1 =
sorted(select_fixed_atoms) = the order the hit list's resid_1..6 columns were written in).
Centroids use
j851's recipe exactly: ARG = mean(NE,NH1,NH2), LYS = NZ, HIS = mean(ND1,NE2,CG).

⛔⛔ GATE: the RMSD this script achieves MUST equal the `best_rmsd` the hit list recorded for
  that chain, to 1e-3 A. If any candidate misses, the .pse is NOT SAVED. A superposition that
  does not reproduce the published number is a different alignment, and a picture built on it
  would show a measurement nobody made.

⛔ The residue NAME at each guidepost and at each matched residue is ASSERTED against the recorded
  pin type before anything is drawn. A silent off-by-one would
  superimpose on the wrong atoms and look perfect.

⛔ THE .pse IS A VIEW, NOT EVIDENCE. The clash numbers printed and shown are measured here and are
  descriptive only; nothing in this script ranks, scores or re-measures the face gate. ⛔ A
  maintained geometry is never evidence of affinity.
"""
import gzip
import os
import sys

import numpy as np

from pymol import cmd

argv = sys.argv[1:]
if "--" in argv:
    argv = argv[argv.index("--") + 1:]


def arg(name, default=None):
    if name in argv:
        return argv[argv.index(name) + 1]
    if default is None:
        sys.exit("missing %s" % name)
    return default


def args_all(name):
    return [argv[i + 1] for i, v in enumerate(argv) if v == name]


design_pdb = arg("--design-pdb")
spec_json = arg("--spec")
hit_tsvs = args_all("--hits")
mirror = arg("--mirror")
afdb_dir = arg("--afdb-dir")
cands = args_all("--cand")
workdir = arg("--workdir")
out = arg("--out")

if os.path.exists(out):
    sys.exit("REFUSING: %s already exists" % out)
if not cands:
    sys.exit("REFUSING: no --cand given")
os.makedirs(workdir, exist_ok=True)

CATION_RES = {"ARG", "LYS"}
IMID_RES = {"HIS"}
ARG_AT = ("NE", "NH1", "NH2")
HIS_AT = ("ND1", "NE2", "CG")
COLOURS = ["ccny_purple", "teal", "orange", "forest", "magenta", "slate"]

import json

import gemmi

spec = json.load(open(spec_json))
spec = spec[list(spec)[0]]
sel = spec["select_fixed_atoms"]
PINS = sorted(sel)                      # B1,C1,D1,E1,F1,G1 -- the resid_1..6 order


def centroid(resname, at):
    if resname == "ARG" and all(a in at for a in ARG_AT):
        return np.mean([at[a] for a in ARG_AT], axis=0)
    if resname == "LYS" and "NZ" in at:
        return np.asarray(at["NZ"], dtype=float)
    if resname == "HIS" and all(a in at for a in HIS_AT):
        return np.mean([at[a] for a in HIS_AT], axis=0)
    return None


def read_chain(path, chain_name):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            doc = gemmi.cif.read_string(fh.read())
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        st = gemmi.read_structure(path)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    for ch in st[0]:
        if ch.name != chain_name:
            continue
        res = []
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info is None or not info.is_amino_acid():
                continue
            at = {a.name: (a.pos.x, a.pos.y, a.pos.z) for a in r}
            if "CA" in at:
                res.append((r.name, r.seqid.num, at))
        return st, ch.name, res
    return None, None, None


# ---------------------------------------------------------------- query constellation
dst = gemmi.read_structure(design_pdb)
dst.remove_alternative_conformations()
dst.remove_hydrogens()
qcen, qtype = [], []
for key in PINS:
    chname, want = key[0], [a.strip() for a in sel[key].split(",")]
    hit = None
    for ch in dst[0]:
        if ch.name != chname:
            continue
        for r in ch:
            at = {a.name: (a.pos.x, a.pos.y, a.pos.z) for a in r}
            if all(a in at for a in want):
                hit = (r.name, at)
                break
    if hit is None:
        sys.exit("REFUSING: guidepost %s not found in %s" % (key, design_pdb))
    resname, at = hit
    ty = "IMID" if "ND1" in want else "CATION"
    exp = IMID_RES if ty == "IMID" else CATION_RES
    if resname not in exp:
        sys.exit("REFUSING: guidepost %s is %s, expected %s" % (key, resname, sorted(exp)))
    qcen.append(np.mean([at[a] for a in want], axis=0))
    qtype.append(ty)
    print("[query] %s %s %s" % (key, resname, ty))
qcen = np.asarray(qcen, dtype=np.float64)
span = float(np.linalg.norm(qcen[:, None, :] - qcen[None, :, :], axis=-1).max())
print("[query] %d pins, span %.2f A" % (len(qcen), span))

# ---------------------------------------------------------------- hit rows
rows = {}
for t in hit_tsvs:
    for ln in list(open(t))[1:]:
        f = ln.rstrip("\n").split("\t")
        rows[f[1]] = f
missing = [c for c in cands if c not in rows]
if missing:
    sys.exit("REFUSING: not in any --hits file: %s" % missing)


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    Pa = (P - Pc) @ R.T + Qc
    rms = float(np.sqrt(((Pa - Q) ** 2).sum(1).mean()))
    return R, Pc, Qc, rms


# ---------------------------------------------------------------- target atoms, for clash
targ = []
for ch in dst[0]:
    if ch.name != "A":
        continue
    for r in ch:
        for a in r:
            targ.append((a.pos.x, a.pos.y, a.pos.z))
targ = np.asarray(targ, dtype=np.float64)
print("[target] chain A heavy atoms %d" % len(targ))

gate_fail, built = [], []
for ci, cname in enumerate(cands):
    f = rows[cname]
    best_rmsd = float(f[2])
    resids = [int(v) for v in f[4:10]]
    types = f[10:16]
    stem, chain_name = cname.rsplit("_", 1)
    if stem.startswith("AF-"):
        path = os.path.join(afdb_dir, "%s.cif.gz" % stem)
    else:
        path = os.path.join(mirror, stem.lower()[1:3], "%s.cif.gz" % stem.lower())
    if not os.path.exists(path):
        sys.exit("REFUSING: no source structure for %s at %s" % (cname, path))
    st, chn, res = read_chain(path, chain_name)
    if res is None:
        sys.exit("REFUSING: chain %r absent from %s" % (chain_name, path))
    by = {}
    for rn, sq, at in res:
        by.setdefault(sq, (rn, at))
    P = []
    for sq, ty in zip(resids, types):
        if sq not in by:
            sys.exit("REFUSING: %s residue %d absent" % (cname, sq))
        rn, at = by[sq]
        exp = IMID_RES if ty == "IMID" else CATION_RES
        if rn not in exp:
            sys.exit("REFUSING: %s residue %d is %s, recorded %s" % (cname, sq, rn, ty))
        g = centroid(rn, at)
        if g is None:
            sys.exit("REFUSING: %s residue %d %s has no complete group" % (cname, sq, rn))
        P.append(g)
    P = np.asarray(P, dtype=np.float64)
    R, Pc, Qc, rms = kabsch(P, qcen)
    ok = abs(rms - best_rmsd) < 1e-3
    print("[gate] %-26s kabsch %.4f  published %.4f  %s"
          % (cname, rms, best_rmsd, "PASS" if ok else "FAIL"))
    if not ok:
        gate_fail.append(cname)
        continue

    # transform the whole chain and write it out
    allxyz, idx = [], []
    for rn, sq, at in res:
        for an, xyz in at.items():
            allxyz.append(xyz)
            idx.append((rn, sq, an))
    A = (np.asarray(allxyz, dtype=np.float64) - Pc) @ R.T + Qc
    d = np.linalg.norm(A[:, None, :] - targ[None, :, :], axis=-1)
    nclash = int((d < 2.5).sum())
    dmin = float(d.min())
    pdbp = os.path.join(workdir, "j893_%s_aligned.pdb" % cname.replace("-", "_"))
    with open(pdbp, "w") as fh:
        ser = 1
        for (rn, sq, an), (x, y, z) in zip(idx, A):
            el = an[0] if not an[0].isdigit() else (an[1] if len(an) > 1 else "C")
            fh.write("ATOM  %5d %-4s%4s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
                     % (ser, an, rn, "X", sq, x, y, z, el))
            ser += 1
        fh.write("END\n")
    built.append(dict(name=cname, pdb=pdbp, rmsd=rms, resids=resids, nres=len(res),
                      nclash=nclash, dmin=dmin, colour=COLOURS[ci % len(COLOURS)]))
    print("       %s  %d res  clash<2.5A %d atoms  min target dist %.2f A"
          % (cname, len(res), nclash, dmin))

if gate_fail:
    sys.exit("REFUSING TO SAVE: superposition gate FAILED for %s -- the alignment does not "
             "reproduce the published best_rmsd, so the picture would not be the measurement"
             % gate_fail)
if not built:
    sys.exit("REFUSING: nothing built")

# ---------------------------------------------------------------- the scene
cmd.set("retain_order", 1)
cmd.set("pdb_retain_ids", 0)
cmd.bg_color("white")
cmd.set_color("ccny_purple", [0x7D / 255.0, 0x55 / 255.0, 0xC7 / 255.0])

cmd.load(design_pdb, "design_input")
cmd.create("target", "design_input and chain A")
cmd.create("guideposts", "design_input and not chain A")
cmd.delete("design_input")
cmd.hide("everything")

cmd.show("cartoon", "target")
cmd.color("gray70", "target")
cmd.set("cartoon_transparency", 0.55, "target")

cmd.show("sticks", "guideposts")
cmd.show("spheres", "guideposts")
cmd.color("red", "guideposts")
cmd.set("sphere_scale", 0.30, "guideposts")
nq = cmd.count_atoms("guideposts")
if nq == 0:
    sys.exit("REFUSING: guideposts selection is empty")
print("[pymol] guidepost atoms %d" % nq)

for b in built:
    obj = "cand_%s" % b["name"].replace("-", "_")
    cmd.load(b["pdb"], obj)
    cmd.show("cartoon", obj)
    cmd.color(b["colour"], obj)
    pin_sel = "%s and resi %s" % (obj, "+".join(str(r) for r in b["resids"]))
    cmd.show("sticks", pin_sel)
    cmd.create("pins_%s" % obj, pin_sel)
    cmd.show("sticks", "pins_%s" % obj)
    cmd.color("yellow", "pins_%s" % obj)
    n = cmd.count_atoms(obj)
    npin = cmd.count_atoms("pins_%s" % obj)
    if n == 0 or npin == 0:
        sys.exit("REFUSING: empty selection for %s (obj %d, pins %d)" % (obj, n, npin))
    print("[pymol] %s atoms %d, pin atoms %d" % (obj, n, npin))
    b["obj"] = obj

cmd.set("label_size", 14)
cmd.set("label_color", "black")

allc = " or ".join(b["obj"] for b in built)
cmd.orient("guideposts or target")
cmd.scene("01_all_four_on_the_constellation", "store",
          message=("All four candidate scaffolds superimposed on the six fixed residues "
                   "(red spheres, span %.1f A). TfR1 crop in grey. Yellow = each candidate's "
                   "own six matched residues. Superposition is Kabsch on those six points "
                   "only." % span))

cmd.disable(allc)
cmd.orient("guideposts")
cmd.scene("02_constellation_alone", "store",
          message=("The six fixed residues alone, against the TfR1 crop. Four ARG and two HIS; "
                   "max pairwise distance %.2f A. This is the geometry every candidate had to "
                   "present." % span))

for b in built:
    cmd.disable(allc)
    cmd.enable(b["obj"])
    cmd.orient("%s or guideposts" % b["obj"])
    cmd.scene("03_%s" % b["obj"], "store",
              message=("%s: %d residues, constellation rmsd %.4f A against the 2.0 A bar. "
                       "Placed on the constellation it puts %d close atom PAIRS (not distinct atoms) within 2.5 A of the "
                       "TfR1 crop; closest approach %.2f A."
                       % (b["name"], b["nres"], b["rmsd"], b["nclash"], b["dmin"])))

cmd.enable(allc)
cmd.orient("guideposts or target")
cmd.scene("01_all_four_on_the_constellation", "recall")
cmd.save(out)
print("[out] %s" % out)
print("[scenes] %s" % ", ".join(cmd.get_scene_list()))
