#!/usr/bin/env python3
"""j897 -- validate the 5XOX-scaffolded arms on the properties that matter for THIS construct.

WHY THIS EXISTS AND THE LIBRARY'S OWN BUILD CHECK IS NOT ENOUGH
--------------------------------------
That check PASSES on these arms -- every spec parses and builds on the library's
own path, and the target stays two chains -- but two of its checks FAIL because their precondition
does not hold here, not because anything is wrong:

  "every guidepost holds exactly 3 atoms"  -> it treats EVERY fixed residue as a guidepost. These
      arms also hold whole motif residues, which are 4-14 atoms each, so the size set is
      [3,4,...,14] by design.
  "ori_token check could run"              -> it expects the fixed set to be ONLY pin atoms
      (12 = 4 x 3). Here it is 313-436 atoms, so its arithmetic test was SKIPPED, which is not
      the same as passing.

⛔⛔ An inapplicable check is not a passing check. So the two underlying properties are verified
HERE, directly, instead of being waved away:
  G1 each BARE guidepost still holds exactly its 3 functional atoms -- a SUBSET, not the whole
     residue. Guideposts are identified by ATOM-NAME SET ({NE,NH1,NH2} or {ND1,NE2,CG}), ⛔ never
     by "has 3 fixed atoms", which would be a coincidence test.
  G2 ori_token was consumed as a PURE TRANSLATION, computed over the GUIDEPOST atoms only:
     built_centroid == input_centroid - ori_token (that same test, restricted correctly).
  G3 the held elements are held as WHOLE residues, and the count equals the element lengths.
  G4 select_unfixed_sequence took effect: element residues are NOT sequence-fixed while the
     guideposts ARE. ⛔ This is the D2 construct decision; if it silently did not apply, the arm
     would import 5XOX's sequence wholesale and nobody would notice.
  G5 the realised binder length sits inside the band the provenance sidecar declares.

Local only. Needs the generator's own environment, because validation runs through the
library's own code path. ⛔ Builds specs; runs no generative model.
"""
import json
import os
import sys

import numpy as np

from rfd3.inference.input_parsing import DesignInputSpecification

FUNC_SETS = {frozenset(("NE", "NH1", "NH2")): "CATION/ARG",
             frozenset(("ND1", "NE2", "CG")): "IMID/HIS",
             frozenset(("NZ",)): "CATION/LYS"}
BINDER = "A"


class Check:
    def __init__(self):
        self.n, self.fails = 0, []

    def __call__(self, ok, label, detail=""):
        self.n += 1
        print("    %s  %-52s %s" % ("PASS" if ok else "FAIL", label, detail))
        if not ok:
            self.fails.append(label)


def pdb_centroid(path, chain, resi, names):
    got = []
    for ln in open(path):
        if not ln.startswith("ATOM"):
            continue
        if ln[21] != chain or int(ln[22:26]) != resi:
            continue
        if ln[12:16].strip() in names:
            got.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.mean(got, axis=0) if got else None


def one(spec_path, ck):
    name = os.path.basename(spec_path)
    print("\n" + "=" * 74)
    print("ARM %s" % name)
    doc = json.load(open(spec_path))
    key, arm = next(iter(doc.items()))
    here = os.path.dirname(os.path.abspath(spec_path))
    prov_p = os.path.join(here, key + "_provenance.json")
    prov = json.load(open(prov_p))
    kwargs = dict(arm)
    pdb_abs = os.path.join(here, arm["input"])
    kwargs["input"] = pdb_abs

    built = DesignInputSpecification(**kwargs).build()
    aa = built[0] if isinstance(built, tuple) else built
    fixed = aa.get_annotation("is_motif_atom_with_fixed_coord").astype(bool)
    anns = set(aa.get_annotation_categories())
    seqfix = (aa.get_annotation("is_motif_atom_with_fixed_seq").astype(bool)
              if "is_motif_atom_with_fixed_seq" in anns else None)

    per_res, seq_res = {}, {}
    for i, (chn, rid, nm, keep) in enumerate(
            zip(aa.chain_id, aa.res_id, aa.atom_name, fixed)):
        if not keep or str(chn) != BINDER:
            continue
        per_res.setdefault(int(rid), []).append(str(nm))
        if seqfix is not None:
            seq_res.setdefault(int(rid), []).append(bool(seqfix[i]))

    # ---- G1: bare guideposts, identified BY ATOM-NAME SET
    gp = {r: FUNC_SETS[frozenset(v)] for r, v in per_res.items()
          if frozenset(v) in FUNC_SETS}
    want_gp = len(prov["bare_guideposts"])
    ck(len(gp) == want_gp, "G1 bare guidepost count (by atom-name set)",
       "%d found, %d declared: %s" % (len(gp), want_gp, prov["bare_guideposts"]))
    ck(all(len(per_res[r]) == 3 for r in gp),
       "G1 each bare guidepost holds exactly 3 atoms (a SUBSET)",
       "sizes %s" % sorted({len(per_res[r]) for r in gp}))

    # ---- G3: elements held as WHOLE residues
    elem_res = {r: v for r, v in per_res.items() if r not in gp}
    want_elem = sum(prov_elem_len(prov, e) for e in prov["elements"])
    ck(len(elem_res) == want_elem, "G3 held element residues",
       "%d fixed residues, %d expected from the element ranges"
       % (len(elem_res), want_elem))
    ck(all(len(v) >= 4 for v in elem_res.values()),
       "G3 element residues are WHOLE (>=4 atoms, never a 3-atom subset)",
       "min %d atoms" % min((len(v) for v in elem_res.values()), default=0))

    # ---- G4: select_unfixed_sequence actually applied
    if seqfix is None:
        ck(False, "G4 is_motif_atom_with_fixed_seq annotation present",
           "ABSENT -- cannot verify the sequence release")
    else:
        el_seqfixed = sum(1 for r in elem_res if any(seq_res.get(r, [])))
        gp_seqfixed = sum(1 for r in gp if any(seq_res.get(r, [])))
        ck(el_seqfixed == 0, "G4 element sequence is RELEASED",
           "%d of %d element residues still sequence-fixed" % (el_seqfixed, len(elem_res)))
        ck(gp_seqfixed == len(gp), "G4 guidepost identity is KEPT",
           "%d of %d guideposts sequence-fixed" % (gp_seqfixed, len(gp)))

    # ---- G2: ori_token consumed as a pure translation, over GUIDEPOST atoms only
    ori = np.asarray(arm["ori_token"], dtype=float)
    deltas = []
    for r in sorted(gp):
        names = set(per_res[r])
        bc = np.mean([aa.coord[i] for i, (chn, rid, nm, keep)
                      in enumerate(zip(aa.chain_id, aa.res_id, aa.atom_name, fixed))
                      if keep and str(chn) == BINDER and int(rid) == r
                      and str(nm) in names], axis=0)
        deltas.append(bc)
    built_c = np.mean(deltas, axis=0)
    # ⛔⛔ THE ATOM SET MUST BE PER RESIDUE TYPE, NOT A UNION. An ARG has a CG, so the union
    # {NE,NH1,NH2,ND1,NE2,CG,NZ} averages FOUR atoms on an arginine instead of three and the
    # pure-translation test reads a spurious 0.2369 A. Resolve the residue NAME first and use
    # exactly its own functional atoms -- the same set select_fixed_atoms named.
    BY_RESN = {"ARG": {"NE", "NH1", "NH2"}, "HIS": {"ND1", "NE2", "CG"}, "LYS": {"NZ"}}
    src = []
    for lab in prov["bare_guideposts"]:
        ch = lab[0]
        resi = resn = None
        for ln in open(pdb_abs):
            if ln.startswith("ATOM") and ln[21] == ch:
                resi, resn = int(ln[22:26]), ln[17:20].strip()
                break
        if resn not in BY_RESN:
            raise SystemExit("guidepost %s is %s -- no functional atom set" % (lab, resn))
        got = pdb_centroid(pdb_abs, ch, resi, BY_RESN[resn])
        if got is None:
            raise SystemExit("guidepost %s: no functional atoms found in the input" % lab)
        src.append(got)
    input_c = np.mean(src, axis=0)
    err = float(np.linalg.norm(built_c - (input_c - ori)))
    ck(err < 0.05, "G2 ori_token consumed as a PURE TRANSLATION (exact)",
       "|delta| = %.4f A over %d guideposts" % (err, len(gp)))

    # ---- G5: realised length inside the declared band
    nres = len({int(r) for c, r in zip(aa.chain_id, aa.res_id) if str(c) == BINDER})
    lo, hi = prov["len_min"], prov["len_max"]
    ck(lo <= nres <= hi, "G5 realised binder length inside the declared band",
       "%d aa in %d-%d" % (nres, lo, hi))
    print("    held %d element residues + %d guideposts | clash vs target %s"
          % (len(elem_res), len(gp), prov["clash_atoms_vs_target"]))
    return nres


def prov_elem_len(prov, e):
    return prov["element_lengths"][str(e)]


def main():
    args = sys.argv[1:]
    if not args:
        sys.exit("usage: j897 <spec.json> ...")
    ck = Check()
    for p in args:
        try:
            one(p, ck)
        except Exception as exc:
            ck(False, "ARM %s built" % os.path.basename(p),
               "%s: %s" % (type(exc).__name__, str(exc).splitlines()[0][:120]))
    print("\n== %d checks, %d FAIL ==" % (ck.n, len(ck.fails)))
    for f in ck.fails:
        print("   FAILED: %s" % f)
    sys.exit(1 if ck.fails else 0)


main()
