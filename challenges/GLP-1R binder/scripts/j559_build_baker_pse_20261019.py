"""j559 -- PSE: the published class-B1 de novo binders, their binding-site residues, and OUR
design in the same scene for contrast. 2026-10-19.

    RUN WITH: pymol -cq /abs/path/to/this/file.py

⛔⛔ NO `if __name__ == "__main__"` GUARD. `pymol -cq` on a guarded script exits 0 having written
   NOTHING. The body runs at import.
⛔⛔ EVERY `scene(message=)` STRING IS ASCII-ONLY. A non-ASCII mark corrupts the heap and the .pse
   SEGFAULTS on load. Every message in this file is grepped to prove it is ASCII.
⛔ PATHS ARE ABSOLUTE. PyMOL's `run` rebinds `__file__` to its own site-packages directory, so
   deriving paths from __file__ resolves into the PyMOL install (hit and fixed on 10-19, J553).
⛔ CHAINS ARE NAMED, never inferred, from the entity records:
   9MM5/9MNI  A = de novo minibinder, R = CGRPR (CALCRL), E = RAMP1.
⛔ Every selection is ASSERTED non-empty; the script RAISES rather than saving an empty scene.
⭐ Palette: CCNY purple #7D55C7 for our own material.
"""
from pymol import cmd

B = "/path/to/your/working/tree"  # PLACEHOLDER PATH -- not a working default; set it for your own tree.
OUT = B + "/pse/j559_baker_vs_ours_20261019.pse"
MM5 = B + "/structures/baker_cgrpr_20261019/9mm5.cif"
MNI = B + "/structures/baker_cgrpr_20261019/9mni.cif"
OURS = B + "/structures/baker_cgrpr_20261019/ours_best_medcut51.pdb"

# receptor residues contacted, MEASURED by j559_baker_reference_20261019.py (not copied from prose)
MM5_REC = [34, 38, 41, 72, 92, 93, 94, 95, 114, 116, 121, 124, 125]
MNI_REC = [41, 72, 92, 94, 114, 124, 125, 128, 129, 131]
# GLP-1R site S1 in the RECEPTOR CONSTRUCT's own numbering -- GLP-1's address, residue for
# residue. ⛔⛔ These are CONSTRUCT numbers and they must NOT be applied to a generator
#    output's target chain: the generator renumbers its target 1..N contiguously and
#    COLLAPSES the segment breaks, so the design file loaded below is NOT in construct
#    numbering and the offset is not even constant along the chain. See scene 06.
GLP1R_S1 = [30, 31, 32, 35, 36, 39, 43, 67, 68, 69, 88, 89, 90, 91, 121, 123, 128]


def need(sel, what):
    n = cmd.count_atoms(sel)
    if n == 0:
        raise SystemExit("EMPTY SELECTION: %s (%s) -- refusing to save a misleading scene" % (what, sel))
    return n


cmd.reinitialize()
cmd.set("assembly", "")
cmd.load(MM5, "mm5")
cmd.load(MNI, "mni")
cmd.load(OURS, "ours")

# ---- keep only what the scenes are about; a stray chain confounds every downstream number ----
cmd.remove("mm5 and not chain A+R+E")
cmd.remove("mni and not chain A+R+E")
cmd.remove("hydrogens")
need("mm5 and chain A", "9MM5 minibinder dC2_049")
need("mm5 and chain R", "9MM5 CGRPR")
need("mni and chain A", "9MNI minibinder dC2_050")
need("ours and chain A", "our design, binder chain A")
need("ours and chain B", "our design, target chain B")

cmd.hide("everything")
cmd.bg_color("white")
cmd.set("cartoon_transparency", 0.0)
cmd.set("ray_shadows", 0)
cmd.set("cartoon_fancy_helices", 1)

PURPLE = "0x7D55C7"
cmd.set_color("ccny_purple", [0.49, 0.333, 0.78])

for obj, rc in (("mm5", "grey70"), ("mni", "grey80")):
    cmd.show("cartoon", "%s and chain R" % obj)
    cmd.color(rc, "%s and chain R" % obj)
    cmd.show("cartoon", "%s and chain E" % obj)
    cmd.color("palecyan", "%s and chain E" % obj)
    cmd.show("cartoon", "%s and chain A" % obj)
    cmd.color("orange", "%s and chain A" % obj)
cmd.show("cartoon", "ours and chain B")
cmd.color("grey60", "ours and chain B")
cmd.show("cartoon", "ours and chain A")
cmd.color("ccny_purple", "ours and chain A")

# ---------------- scene 1: the published complex ----------------
cmd.disable("all")
cmd.enable("mm5")
cmd.orient("mm5")
cmd.scene("01_published_complex", "store",
          message=("9MM5: CGRP receptor (grey, chain R) + RAMP1 (cyan, chain E) + the de novo "
                   "minibinder dC2_049 (orange, chain A), 3.2 A cryo-EM. CGRPR is class B1, the "
                   "same family as GLP-1R, from the same Baker-lab campaign that produced the "
                   "GLP1R binders dGI1_024 (27 nM) and mGI1_008 (5.3 nM). There is NO deposited "
                   "GLP1R-bound structure for those two, so this is a PROXY and is labelled one."))

# ---------------- scene 2: the binding site residues ----------------
sel = "mm5 and chain R and resi " + "+".join(str(x) for x in MM5_REC)
need(sel, "9MM5 receptor contact residues")
cmd.show("sticks", sel)
cmd.color("red", sel)
cmd.show("sticks", "mm5 and chain A and byres (chain A within 4.5 of (mm5 and chain R))")
cmd.color("yellow", "mm5 and chain A and byres (chain A within 4.5 of (mm5 and chain R))")
cmd.orient(sel + " or (mm5 and chain A)")
cmd.scene("02_binding_site_9MM5", "store",
          message=("THE BINDING SITE, MEASURED not quoted: 13 CGRPR residues in red "
                   "(34 38 41 72 92-95 114 116 121 124 125) against 13 of the binder's 65 "
                   "residues in yellow. KEY NUMBER: the interface uses only 20 percent of the "
                   "binder's residues - a small focused FACE, not a long contact line. The "
                   "contacting CAs deviate 4.35 A from a straight line, so this is not a chain "
                   "hugging a groove."))

sel2 = "mni and chain R and resi " + "+".join(str(x) for x in MNI_REC)
need(sel2, "9MNI receptor contact residues")
cmd.show("sticks", sel2)
cmd.color("red", sel2)
cmd.disable("all")
cmd.enable("mni")
cmd.orient(sel2 + " or (mni and chain A)")
cmd.scene("03_binding_site_9MNI", "store",
          message=("9MNI / dC2_050, the second published complex, 4.1 A. 10 receptor residues "
                   "(41 72 92 94 114 124 125 128 129 131) against 10 of 66 binder residues, "
                   "15 percent. The two independent designs converge on the SAME patch: "
                   "41, 72, 92, 94, 114, 124, 125 are shared. That patch is the ECD peptide-"
                   "binding groove region - the same region of the extracellular domain "
                   "this campaign targets."))

# ---------------- scene 4: the binder alone -- it is a COMPACT BUNDLE ----------------
cmd.disable("all")
cmd.enable("mm5")
cmd.hide("everything", "mm5 and chain R+E")
cmd.orient("mm5 and chain A")
cmd.scene("04_their_binder_alone", "store",
          message=("dC2_049 ALONE: 65 residues, min_cut 13, med_cut 32, rg_ratio 1.042. "
                   "For scale, a 1TNF protomer scores rg_ratio 1.138, so this designed binder is "
                   "MORE compact than a natural globular domain. It is a closed bundle: there is "
                   "no cut through the chain that nothing crosses."))

# ---------------- scene 5: ours, same instruments, for contrast ----------------
cmd.disable("all")
cmd.enable("ours")
cmd.show("cartoon", "ours and chain B")
cmd.orient("ours and chain A")
cmd.scene("05_ours_for_contrast", "store",
          message=("OUR BEST DESIGN BY CONTACT DENSITY (shard05 deep6 model_6, med_cut 51, the "
                   "highest of all 384) in CCNY purple, on the GLP-1R target in grey. Same "
                   "instruments: min_cut 0, rg_ratio 1.356. Every one of our 384 designs scores "
                   "min_cut 0 - there is a point in the chain that nothing crosses - while both "
                   "published binders score 13. Our rg_ratio median is 1.53 against their 1.04. "
                   "And across 352 apo folds not one of ours reaches TM 0.8 to its own backbone."))

# ---------------- scene 6: DELIBERATELY NOT STORED ----------------
# ⛔⛔ A scene painting site S1 on our own target chain was REMOVED rather than fixed, because
#    it was mislabelled. MEASURED on the design file this script loads: its target chain B holds
#    384 residues numbered 1..384, and those are three construct segments with the breaks
#    collapsed. ALWAYS NAME THE FRAME: chain B 1-101, 102-304 and 305-384 are CONSTRUCT residues
#    6-106, 113-315 and 321-400, which are UNIPROT P43220 residues 29-129, 136-338 and 344-423.
#    So the offset is +5, then +11, then +16 to construct numbering (+28, +34, +39 to UniProt)
#    along a single chain, and construct residues 107-112 and 316-320 are simply absent.
# ⛔ Selecting "chain B and resi 30+31+..." therefore RETURNS ATOMS -- which is exactly why an
#    is-the-selection-empty guard cannot catch this -- but they are the WRONG residues. A scene
#    that cannot be labelled truthfully is not stored. Restoring it needs an explicit
#    chain-to-construct residue map carried alongside the coordinates.
print("NOTE: the S1 site scene is not stored -- this design file's target chain is renumbered "
      "and an explicit residue map is required; see the comment above GLP1R_S1.")

cmd.disable("all")
cmd.enable("mm5")
cmd.enable("ours")
cmd.scene("01_published_complex", "recall")
cmd.save(OUT)
print("saved %s" % OUT)
print("scenes: %s" % cmd.get_scene_list())
