#!/usr/bin/env python3
"""Score the r4 co-fold: what does Boltz build from the adopted bias cohort's SEQUENCES alone?

Metrics and their justification are inherited from `s36_score_g9.py` unchanged -- same cutoffs,
same stratification, same two-layout discovery. What differs is that the classes are read from the
RUN'S OWN MANIFEST rather than a hard-coded map, so a design can never be silently misclassified.

⛔⛔ WHAT THIS SCORE CANNOT DO
  - ⛔ It cannot rank step-1 selectivity (`J989` §1): no bound-state interface metric from ANY model
    ranks WT-vs-G12D here, and the crystals invert the same way. Nothing below speaks to it.
  - ⛔ The designs were SELECTED on design-stage `cont_pep` (`s59`), so a high folded peptide-contact
    count is partly selection. The fold tests whether the SEQUENCE rebuilds the pose, not whether
    the cohort engages the peptide.
  - ⛔ Global ipTM is recorded and never used: every input is 4 chains
    ([[global-iptm-is-confounded-by-a-multimeric-target]]).
  - ⛔ n = 3 designs per cell. That is FAR too small for `ADR 0063`'s arm3-vs-pinned falsifier,
    which needs the 8-batch campaign. Any per-cell number here is descriptive only.

⭐ The read, per `J986` §3: peptide-STRATIFIED contact count first, binder<->peptide
`pair_chains_iptm` second, the binder's own pLDDT as the "did it even fold" sanity check, and the
floor set by the cross-paired non-binders folded IN THIS SAME BATCH.

Usage: python3 -I s60_score_r4_cofold.py <batch_out_root> <manifest.tsv> <out.txt>
"""
import collections
import csv
import glob
import json
import os
import sys
import warnings

import numpy as np
from Bio.PDB import PDBParser

warnings.filterwarnings("ignore")

ROOT, MAN, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
IFACE = 5.0      # [M] s36:38 -- binder residue is interface if within this of the target
CONTACT = 4.5    # [M] s36:39 -- residue-residue contact cutoff

P = PDBParser(QUIET=True)


def res_list(chain):
    return [r for r in chain if r.id[0] == " "]


def heavy(r):
    return [a for a in r if a.element != "H"]


def metrics(pdb):
    """Verbatim in behaviour from s36:57-87. Chains are OURS: A=HLA B=b2m C=peptide D=binder."""
    m = P.get_structure("m", pdb)[0]
    ids = [c.id for c in m]
    if "D" not in ids:
        return None
    binder = res_list(m["D"])
    out = {"binder_plddt": float(np.mean([a.get_bfactor() for r in binder for a in r])),
           "binder_len": len(binder)}
    for tgt, label in (("C", "peptide"), ("A", "hla")):
        if tgt not in ids:
            continue
        tres = res_list(m[tgt])
        tat = np.array([a.coord for r in tres for a in heavy(r)])
        iface, contacts = [], 0
        for r in binder:
            ba = np.array([a.coord for a in heavy(r)])
            if np.linalg.norm(tat[:, None, :] - ba[None, :, :], axis=-1).min() <= IFACE:
                iface.append(r)
            for tr in tres:
                ta = np.array([a.coord for a in heavy(tr)])
                if np.linalg.norm(ta[:, None, :] - ba[None, :, :], axis=-1).min() <= CONTACT:
                    contacts += 1
        out[f"n_iface_{label}"] = len(iface)
        out[f"contacts_{label}"] = contacts
        out[f"iface_plddt_{label}"] = (
            float(np.mean([a.get_bfactor() for r in iface for a in r])) if iface else float("nan"))
    return out


def find_prediction_dirs(root):
    """⛔ Boltz writes TWO layouts -- `boltz_results_<input>/predictions/<input>/` for a single input
    and `boltz_results_<dirname>/predictions/<input>/` for a directory. A scorer that knows only one
    SILENTLY FINDS NOTHING and prints a clean, well-formatted empty report (`J986` §6: this happened
    to both s36 and s37 on batches that had written every file). Locate directories that actually
    contain model PDBs. [[a-clean-empty-report-is-a-failure-mode]]"""
    return sorted({os.path.dirname(p) for p in
                   glob.glob(os.path.join(root, "**", "*_model_*.pdb"), recursive=True)})


def pair_iptm(c, a, b):
    """pair_chains_iptm, either index order. Our A,B,C,D -> boltz 0,1,2,3 (s36:132)."""
    p = c.get("pair_chains_iptm", {})
    for x, y in ((a, b), (b, a)):
        if str(x) in p and str(y) in p[str(x)]:
            return p[str(x)][str(y)]
    return float("nan")


man = {r["tag"]: r for r in csv.DictReader(open(MAN), delimiter="\t")}
preds = find_prediction_dirs(ROOT)

L = []
W = L.append
W("# s60 -- Boltz co-fold of the ADOPTED r4 bias cohort (A:-0.495,E:-0.354,L:-0.092).")
W("# single-sequence (msa: empty on every chain), diffusion_samples 5, recycling 3, seed 20261104")
W("# chain ids are OURS: A=HLA B=b2m C=peptide D=binder")
W(f"# interface = binder residue within {IFACE} A of the target chain; contacts at {CONTACT} A")
W("# ⛔ NOT usable for step-1 selectivity (J989 §1). ⛔ Global ipTM recorded, never read.")
W("")

# ---- ⛔ count the files BEFORE reading any number ([[a-clean-empty-report-is-a-failure-mode]]) ----
W(f"## File census -- ⛔ rc=0 is not evidence a fold happened ([[boltz-oom-exits-zero]])")
W(f"  manifest rows      : {len(man)}")
W(f"  prediction dirs    : {len(preds)}")
missing = sorted(set(man) - {os.path.basename(p) for p in preds})
n_pdb = len(glob.glob(os.path.join(ROOT, "**", "*_model_*.pdb"), recursive=True))
W(f"  model PDBs written : {n_pdb}   (expected {len(man) * 5} = {len(man)} inputs x 5 samples)")
if missing:
    W(f"  ⛔⛔ {len(missing)} INPUT(S) PRODUCED NOTHING: {', '.join(missing)}")
if n_pdb != len(man) * 5:
    W(f"  ⛔⛔ PDB COUNT MISMATCH -- a silent OOM writes nothing and still exits 0.")
if not missing and n_pdb == len(man) * 5:
    W("  ✅ complete: every manifest row has a prediction dir and 5/5 samples.")
W("")

rows = []
for pred in preds:
    tag = os.path.basename(pred)
    if tag not in man:
        W(f"  ⚠️ unmanifested prediction dir, skipped: {tag}")
        continue
    mrow = man[tag]
    pdbs = sorted(glob.glob(os.path.join(pred, "*_model_*.pdb")))
    confs = [json.load(open(c)) for c in
             sorted(glob.glob(os.path.join(pred, "confidence_*_model_*.json")))]
    ms = [m for m in (metrics(p) for p in pdbs) if m]
    if not ms:
        W(f"  ⛔ {tag}: PDBs present but no chain D -- not scorable")
        continue
    rows.append(dict(
        tag=tag, klass=mrow["klass"], cell=mrow["cell"], n=len(pdbs),
        cpep=[m.get("contacts_peptide", 0) for m in ms],
        chla=[m.get("contacts_hla", 0) for m in ms],
        npep=[m.get("n_iface_peptide", 0) for m in ms],
        bpl=[m["binder_plddt"] for m in ms],
        ipep=[m.get("iface_plddt_peptide", float("nan")) for m in ms],
        dpep=[pair_iptm(c, 3, 2) for c in confs],
        dhla=[pair_iptm(c, 3, 0) for c in confs],
        iptm=[c.get("iptm", float("nan")) for c in confs]))

if not rows:
    W("⛔ no scorable results under " + ROOT)
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    sys.exit(0)

mean = lambda v: float(np.nanmean(v))
nb = [r for r in rows if r["klass"].startswith("NONBINDER")]
pos = [r for r in rows if r["klass"].startswith("POSITIVE")]
des = [r for r in rows if r["klass"] == "DESIGN"]

# ---- the floor, from the non-binders folded in THIS batch (`J986` §3 rule 5) ----
W("## The floor and the ceiling, both folded in THIS batch")
W(f"  {'tag':<26}{'class':<24}{'bindPLDDT':>10}{'D-C ipTM':>10}{'cont_C':>8}{'cont_A':>8}{'glob':>7}")
for r in sorted(nb + pos, key=lambda r: r["klass"]):
    W(f"  {r['tag']:<26}{r['klass']:<24}{mean(r['bpl']):>10.3f}{mean(r['dpep']):>10.3f}"
      f"{mean(r['cpep']):>8.1f}{mean(r['chla']):>8.1f}{mean(r['iptm']):>7.3f}")
W("")
HAVE_FLOOR = bool(nb)
if not HAVE_FLOOR:
    W("  ⛔⛔ NO NON-BINDER IN THIS BATCH -- every design number below is uninterpretable.")
    FLOOR_C = FLOOR_I = float("nan")
else:
    FLOOR_C = max(mean(r["cpep"]) for r in nb)
    FLOOR_I = max(mean(r["dpep"]) for r in nb)
    LEN_C = min(mean(r["cpep"]) for r in nb)
    LEN_I = min(mean(r["dpep"]) for r in nb)
    W(f"  ⇒ FLOOR (worst case = the HIGHEST non-binder): contacts_peptide {FLOOR_C:.1f}, "
      f"D-C ipTM {FLOOR_I:.3f}")
    W(f"  ⇒ lenient floor (the LOWER non-binder):        contacts_peptide {LEN_C:.1f}, "
      f"D-C ipTM {LEN_I:.3f}")
    # ⛔ With n=2 the floor is POORLY DETERMINED. Quoting one conjunction count would hide which
    # non-binder set the bar, so both are reported and the gap is the honest uncertainty.
    W(f"  ⛔ the two non-binders span {LEN_I:.3f}-{FLOOR_I:.3f} on D-C ipTM ⇒ the ipTM floor is "
      f"POORLY DETERMINED at n=2;")
    W("     ⇒ the conjunction is reported against BOTH below. ⛔ Neither is a calibrated threshold.")
    W("  ⚠️ n=2 non-binders ⇒ this is a DIRECTION, not a calibrated threshold, and real re-fold")
    W("     noise is ~1.56x the within-batch spread ([[within-batch-se-understates-refold-noise]]).")
W("")

# ---- designs, ranked on the CONJUNCTION ([[rank-on-conjunction-not-mean-of-fractions]]) ----
# ⛔ With no floor there is no test. Recording that as `False` would print "0 of N clear", which
# reads as a RESULT when it is a missing control ([[a-clean-empty-report-is-a-failure-mode]]).
# Found by running this scorer's failing branch deliberately, before the real batch landed.
for r in des:
    r["pass_c"] = (mean(r["cpep"]) > FLOOR_C) if HAVE_FLOOR else None
    r["pass_i"] = (mean(r["dpep"]) > FLOOR_I) if HAVE_FLOOR else None
    r["conj"] = (bool(r["pass_c"]) and bool(r["pass_i"])) if HAVE_FLOOR else None

W("## Designs -- ranked on the CONJUNCTION, never a mean of fractions")
W(f"  {'tag':<44}{'cell':<26}{'bPLDDT':>8}{'D-C':>7}{'cont_C':>8}{'cont_A':>8}{'C>f':>5}{'I>f':>5}{'BOTH':>6}")
def mark(x):
    return "?" if x is None else ("Y" if x else ".")


for r in sorted(des, key=lambda r: (-int(bool(r["conj"])), -mean(r["cpep"]), -mean(r["dpep"]))):
    W(f"  {r['tag']:<44}{r['cell']:<26}{mean(r['bpl']):>8.3f}{mean(r['dpep']):>7.3f}"
      f"{mean(r['cpep']):>8.1f}{mean(r['chla']):>8.1f}"
      f"{mark(r['pass_c']):>5}{mark(r['pass_i']):>5}"
      f"{('⭐' if r['conj'] else '') if r['conj'] is not None else '?':>6}")
W("")
if HAVE_FLOOR:
    lenient = sum(1 for r in des if mean(r["cpep"]) > LEN_C and mean(r["dpep"]) > LEN_I)
    W(f"  ⇒ {sum(bool(r['conj']) for r in des)} of {len(des)} designs clear BOTH floors "
      f"(strict, vs the HIGHER non-binder).")
    W(f"  ⇒ {lenient} of {len(des)} clear both against the LENIENT floor (the lower non-binder).")
    W(f"  ⇒ {sum(1 for r in des if mean(r['cpep']) > FLOOR_C)} of {len(des)} clear the "
      f"peptide-CONTACT floor alone; "
      f"{sum(1 for r in des if mean(r['dpep']) > FLOOR_I)} clear the ipTM floor alone.")
else:
    W("  ⛔⛔ NO CONJUNCTION COULD BE EVALUATED -- the floor is missing, so this is not '0 of N'.")
W("  C>f = contacts to peptide above the non-binder floor;  I>f = binder<->peptide ipTM above it")
W("  ⛔ cont_A (HLA) is shown but NEVER ranked on: the highest HLA-contact count in the whole G9")
W("     table belonged to a NON-BINDER (J986 §3 rule 3).")
W("")

# ---- per-cell, descriptive only ----
W("## Per cell -- ⛔ DESCRIPTIVE ONLY, n=3 per cell")
by = collections.defaultdict(list)
for r in des:
    by[r["cell"]].append(r)
W(f"  {'cell':<28}{'n':>3}{'conj':>6}{'cont_C mean':>13}{'D-C mean':>10}{'bPLDDT':>9}")
for cell in sorted(by):
    v = by[cell]
    W(f"  {cell:<28}{len(v):>3}{(sum(bool(r['conj']) for r in v) if HAVE_FLOOR else '?'):>6}"
      f"{np.mean([mean(r['cpep']) for r in v]):>13.1f}"
      f"{np.mean([mean(r['dpep']) for r in v]):>10.3f}"
      f"{np.mean([mean(r['bpl']) for r in v]):>9.3f}")
W("")
W("  ⛔⛔ This table does NOT adjudicate ADR 0063's falsifier (pinned vs arm3). That needs the")
W("     8-batch campaign and a per-backbone paired comparison; n=3 cells here cannot carry it.")
W("")
W("## Sanity: did the binders fold at all?")
# ⛔ The question is "did this binder fold", which is an ABSOLUTE property of the model's own
# confidence -- NOT "did it fold as well as a crystallised minibinder". An earlier version of this
# section compared each design against the worst non-binder and flagged 17 of 18 as suspect; but
# the two non-binder binders are REAL, well-resolved proteins (9MIN's 136-mer, 9O5S's 102-mer),
# chosen precisely because they fold confidently ([[s35 docstring]]), so they sit at the top of the
# scale by construction. A de novo design at ~88 pLDDT has folded. Gate on PLDDT_FOLDED instead and
# report the relative gap separately, as the expected observation it is.
PLDDT_FOLDED = 70.0
lo = [r for r in des if mean(r["bpl"]) < PLDDT_FOLDED]
if lo:
    W(f"  ⚠️ {len(lo)} design(s) below pLDDT {PLDDT_FOLDED:.0f} ⇒ for these a low interface score is")
    W("     ambiguous between 'no interface' and 'the binder did not fold':")
    for r in lo:
        W(f"       {r['tag']}  binder pLDDT {mean(r['bpl']):.3f}")
else:
    W(f"  ✅ every design folds above pLDDT {PLDDT_FOLDED:.0f} "
      f"(range {min(mean(r['bpl']) for r in des):.1f}-{max(mean(r['bpl']) for r in des):.1f})")
    W("     ⇒ a low interface score here means 'NO INTERFACE', not 'the binder failed to fold'.")
    W("     ⭐ That is exactly the ambiguity a scrambled decoy could not have resolved.")
if HAVE_FLOOR:
    W(f"  ⚠️ The two non-binder binders fold MORE confidently "
      f"({min(mean(r['bpl']) for r in nb):.1f}-{max(mean(r['bpl']) for r in nb):.1f}) than our "
      f"designs -- ⛔ EXPECTED, they are crystallised")
    W("     proteins, and ⛔ NOT evidence of a design defect.")

txt = "\n".join(L) + "\n"
open(OUT, "w").write(txt)
print(txt)
