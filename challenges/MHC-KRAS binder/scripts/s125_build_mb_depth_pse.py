"""s125 -- build the ONE .pse for the J1058 MB-depth co-fold round (standing rule: one per arm).

Run under:  pymol -cq s125_build_mb_depth_pse.py -- <extract_root> <s124_per_design.tsv> <manifest.tsv> <out.pse>

⛔ `pymol -cq script.py` does NOT run the script as `__main__`, so a `__main__` guard writes
   NOTHING and exits 0 [[pymol-script-not-main]] ⇒ all code here is top-level.

⭐⭐ WHAT IT SHOWS, AND WHY THIS CONSTRUCTION. The finding of this whole workstream is not that our
   designs look plausible -- it is that **the floor is brutal**: `[M]` cross-paired non-binders,
   which are REAL PUBLISHED BINDERS FOR OTHER pMHCs folded against OUR target, make 8.6-11.2
   peptide contacts and ~0.88-0.89 binder<->peptide ipTM **without being designed for it**, and the
   typical design touches the peptide LESS than they do. A .pse of our candidates alone would hide
   exactly that and would read as a set of handsome complexes. ⇒ **every candidate is shown beside
   the non-binder that defines its own shard's floor**, in the same frame, so the comparison is
   seen rather than asserted.

⛔ THE COHORT IS READ FROM `s124`'s TSV, ⛔ NEVER RETYPED. Designs are selected by the
   pre-registered rule (clears the conjunction) and ordered by CONTACT MARGIN, which is the
   pre-registered ranking statistic -- ⛔ not by ipTM, which may reject but never promote.

⛔ CHAIN LAYOUT IS THE TRAP. In the REFOLD (what Boltz built) the layout is
       A = HLA    B = b2m    C = peptide    D = binder
   ⛔ which is NOT the design-PDB layout (A=binder B=HLA C=b2m D=peptide). Every selection below
   names its chain explicitly; a mix-up would colour the HLA as the binder and still render a
   handsome session [[always-name-the-pdb-chain]].

TRAPS THIS SCRIPT IS WRITTEN AGAINST, all previously paid for:
  - a DISABLED PyMOL group hides its enabled members, so 12 .pse once shipped showing only the bare
    target while every check passed ⇒ members AND groups are enabled explicitly and asserted
    [[pymol-group-disabled-hides-enabled-members]].
  - a scene message with a non-ASCII character SEGFAULTS the saved .pse
    [[pymol-scene-message-must-be-ascii]] ⇒ `ascii_msg()` strips, and every string here is ASCII.
  - `retain_order 1` silently corrupts `cmd.align`: `[M]` a 276-residue fit collapsed to **14 atoms**
    at a NUMERICALLY BETTER rms than the correct 218-atom fit ⇒ ⛔ an RMSD can never validate a
    superposition; `fit()` gates on ATOM COUNT [[pymol-retain-order-corrupts-align]].
"""
import csv
import os
import sys

from pymol import cmd

argv = sys.argv[1:]
if len(argv) != 4:
    sys.exit(__doc__)
ROOTS, TSV, MANIFEST, OUT = argv
# ⭐ `<extract_root>` may be a COLON-SEPARATED LIST, so a cohort spanning several rounds (the
#   validated designs live in three different extract trees) can be shown in ONE session. Each root
#   is searched in order for the tag; ⛔ the FIRST match wins, so order the roots newest-first if a
#   tag could legitimately exist in two rounds.
ROOT_LIST = [r for r in ROOTS.split(":") if r]
for _r in ROOT_LIST:
    if not os.path.isdir(_r):
        sys.exit("REFUSING: extract root %s is not a directory" % _r)
print("extract roots (%d): %s" % (len(ROOT_LIST), ", ".join(os.path.basename(r) for r in ROOT_LIST)))

HLA_LEN = 276
MIN_FIT = 200
MAX_DESIGNS = 6


def ascii_msg(s):
    """A non-ASCII char in a scene message SEGFAULTS the saved .pse. Strip, never trust."""
    return s.encode("ascii", "ignore").decode("ascii")


def need(path):
    if not os.path.exists(path):
        sys.exit("REFUSING: missing input %s" % path)
    return path


def pred_path(shard, tag):
    """Find <tag>'s model_0 across every extract root. ⛔ Returns the first HIT; if the row's own
    shard is absent from every root it falls back to any shard that has the tag, because for a
    rendering the seed is not load-bearing -- ⛔ but it is reported so the session is auditable."""
    for r in ROOT_LIST:
        p = os.path.join(r, "boltz_results_shard%s" % shard, "predictions", tag,
                         "%s_model_0.pdb" % tag)
        if os.path.exists(p):
            return p
    for r in ROOT_LIST:
        for sd in sorted(os.listdir(r)):
            if not sd.startswith("boltz_results_shard"):
                continue
            p = os.path.join(r, sd, "predictions", tag, "%s_model_0.pdb" % tag)
            if os.path.exists(p):
                print("  !! %s not in shard%s; using %s/%s instead (seed differs, rendering only)"
                      % (tag, shard, os.path.basename(r), sd))
                return p
    return os.path.join(ROOT_LIST[0], "boltz_results_shard%s" % shard, "predictions", tag,
                        "%s_model_0.pdb" % tag)          # let need() report the miss


def fit(mobile_sel, target_sel, label):
    """Superpose HLA-on-HLA, then ASSERT the fit actually used the HLA.

    ⛔⛔ Gate on the ATOM COUNT, never the rms: a 14-atom fit on a 276-residue chain reports a
    BETTER rms than the correct 218-atom fit, returns no error, and would ship a .pse whose whole
    visual claim was false with every other check green.
    """
    rms, n = cmd.align(mobile_sel, target_sel)[:2]
    print("  fit %-50s rms %6.3f A over %4d / %d atoms" % (label, rms, n, HLA_LEN))
    if n < MIN_FIT:
        sys.exit("REFUSING: %s fitted only %d of %d HLA CA atoms (floor %d). A low RMSD over a "
                 "handful of atoms is not an alignment." % (label, n, HLA_LEN, MIN_FIT))
    return rms, n


# ── the cohort, READ not retyped ──────────────────────────────────────────────────────────────
rows = list(csv.DictReader(open(TSV), delimiter="\t"))
if not rows:
    sys.exit("REFUSING: %s holds no rows" % TSV)
clear = [r for r in rows if r["conj"] == "True"]
# ⛔⛔ DEDUPLICATE BY MOLECULE, i.e. by (design, mpnn_rank). A control/re-draw tag is folded in
#    EVERY shard, and a backbone can clear at several ranks, so ranking the raw rows shows THE SAME
#    MOLECULE several times at different seeds and crowds the distinct candidates out of the
#    session. One object per distinct molecule, kept at its BEST margin.
best = {}
for r in clear:
    k = (r["design"], r["rank"])
    if k not in best or float(r["margin"]) > float(best[k]["margin"]):
        best[k] = r
clear = sorted(best.values(), key=lambda r: -float(r["margin"]))   # ⛔ margin, never ipTM
picked = clear[:MAX_DESIGNS]
print("distinct clearing molecules: %d (from %d clearing rows across seeds)"
      % (len(clear), sum(1 for r in rows if r["conj"] == "True")))

man = list(csv.DictReader(open(MANIFEST), delimiter="\t"))
nb_tags = sorted({r["tag"] for r in man if r["klass"] == "NONBINDER_crosspaired"})
if len(nb_tags) != 2:
    sys.exit("REFUSING: %d cross-paired non-binders in the manifest, expected 2 -- without the "
             "floor there is nothing to compare the candidates against." % len(nb_tags))

if not picked:
    print("!! NO design cleared the conjunction. Building a FLOOR-ONLY session: that is a real")
    print("   outcome and the .pse must show it rather than silently rendering nothing.")

print("cohort for the .pse: %d design(s) + %d non-binder floor references"
      % (len(picked), len(nb_tags)))

cmd.reinitialize()
cmd.set("assembly", "")
# retain_order is DELIBERATELY NOT SET -- see fit().
cmd.set("cartoon_transparency", 0.0)
cmd.bg_color("white")

anchor = None
loaded_designs, loaded_nb = [], []

# ── load the designs ──────────────────────────────────────────────────────────────────────────
for i, r in enumerate(picked, start=1):
    tag, shard = r["tag"], r["shard"]
    obj = "d%d_%s_m%+.1f" % (i, tag[:28].replace(".", "_"), float(r["margin"]))
    obj = obj.replace("+", "p").replace("-", "m")
    cmd.load(need(pred_path(shard, tag)), obj)
    loaded_designs.append((obj, r))
    if anchor is None:
        anchor = obj
    else:
        fit("%s and chain A and name CA" % obj, "%s and chain A and name CA" % anchor,
            obj + " HLA-on-HLA")

# ── load each shard's own non-binder floor, from the SAME shard as a picked design where possible
nb_shards = sorted({r["shard"] for r in picked}) or sorted({r["shard"] for r in rows})
for tag in nb_tags:
    for shard in nb_shards[:1]:
        p = pred_path(shard, tag)
        if not os.path.exists(p):
            continue
        obj = "floor_%s_s%s" % (tag.split("_")[0], shard)
        cmd.load(p, obj)
        loaded_nb.append(obj)
        if anchor is None:
            anchor = obj
        else:
            fit("%s and chain A and name CA" % obj, "%s and chain A and name CA" % anchor,
                obj + " HLA-on-HLA")
if not loaded_nb:
    sys.exit("REFUSING: no non-binder floor structure could be loaded -- refusing to ship a .pse "
             "that shows candidates with nothing to compare them to.")

# ── styling. ⛔ every selection names its chain explicitly (REFOLD layout A=HLA B=b2m C=pep D=binder)
cmd.hide("everything")
for obj, _r in loaded_designs:
    cmd.show("cartoon", "%s and chain A+B" % obj)
    cmd.color("grey80", "%s and chain A+B" % obj)
    cmd.show("cartoon", "%s and chain C" % obj)
    cmd.color("orange", "%s and chain C" % obj)
    cmd.show("sticks", "%s and chain C" % obj)
    cmd.show("cartoon", "%s and chain D" % obj)
    cmd.color("purple", "%s and chain D" % obj)      # CCNY purple for OUR designs
for obj in loaded_nb:
    cmd.show("cartoon", "%s and chain A+B" % obj)
    cmd.color("grey90", "%s and chain A+B" % obj)
    cmd.show("cartoon", "%s and chain C" % obj)
    cmd.color("orange", "%s and chain C" % obj)
    cmd.show("cartoon", "%s and chain D" % obj)
    cmd.color("firebrick", "%s and chain D" % obj)   # the NON-BINDER, i.e. the floor
    cmd.set("cartoon_transparency", 0.45, obj)

# ── groups. ⛔ a DISABLED group hides its ENABLED members ⇒ enable members AND groups, then assert
G_D, G_F = "designs", "floor_nonbinders"
if loaded_designs:
    cmd.group(G_D, " ".join(o for o, _ in loaded_designs))
cmd.group(G_F, " ".join(loaded_nb))
for o, _r in loaded_designs:
    cmd.enable(o)
for o in loaded_nb:
    cmd.enable(o)
for g in ([G_D] if loaded_designs else []) + [G_F]:
    cmd.enable(g)

# ── scenes, one per candidate plus the floor comparison. ⛔ ASCII ONLY in every message.
cmd.orient("chain C")
cmd.zoom("chain C", 12.0)
cmd.scene("000_all", "store",
          message=ascii_msg("J1058 MB-depth. PURPLE = our designs, FIREBRICK = cross-paired "
                            "NON-binders (the in-batch floor), ORANGE = the G12D peptide. "
                            "The floor references are real published binders for OTHER pMHCs, "
                            "folded against our target without being designed for it."))
for i, (obj, r) in enumerate(loaded_designs, start=1):
    cmd.disable(G_D)
    cmd.enable(G_F)
    cmd.enable(obj)
    cmd.orient("%s and chain C+D" % obj)
    cmd.zoom("%s and chain C+D" % obj, 6.0)
    msg = ("%s  shard %s  MPNN rank %s | cont_pep %.1f, margin %+.1f over its OWN shard floor | "
           "binder-peptide ipTM %.3f | CANDIDATE, not a result."
           % (r["design"], r["shard"], r["rank"], float(r["cpep"]), float(r["margin"]),
              float(r["dpep"])))
    cmd.scene("%03d_%s" % (i, obj[:20]), "store", message=ascii_msg(msg))
cmd.enable(G_D)
cmd.enable(G_F)

# ── assertions BEFORE saving ──────────────────────────────────────────────────────────────────
names = cmd.get_names("objects")
for o, _r in loaded_designs:
    if o not in names:
        sys.exit("REFUSING: object %s vanished before save" % o)
n_shown = 0
for o in [x for x, _ in loaded_designs] + loaded_nb:
    if cmd.count_atoms("%s and chain D" % o) == 0:
        sys.exit("REFUSING: %s has ZERO chain-D (binder) atoms -- this is the 12-empty-.pse "
                 "signature; a session showing only the bare target passes every other check." % o)
    n_shown += 1
print("  asserted: %d objects each carry chain-D binder atoms" % n_shown)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
cmd.scene("000_all", "recall")
cmd.save(OUT)
if not os.path.exists(OUT) or os.path.getsize(OUT) < 10000:
    sys.exit("REFUSING: %s was not written or is implausibly small" % OUT)
print("wrote %s  (%.1f MB, %d objects, %d scenes)"
      % (OUT, os.path.getsize(OUT) / 1e6, len(names), len(cmd.get_scene_list())))
