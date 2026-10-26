#!/usr/bin/env python3
"""j641 -- .pse of the best Challenge-3 folds by the Baker filter, with the FIXED RESIDUES marked.

Selection comes from `data/j637_baker_filter_c3_20261024_v2.tsv` (see `J638`/`J640`).

## WHAT "FIXED RESIDUES" ARE HERE, AND HOW THE POSITIONS WERE DERIVED
⛔ Not from memory and ⛔ not hardcoded from a note -- read off the arm's OWN shipped contig
specification. That file is the generator input for the batch and is **not** deposited; the
contig itself is reproduced here, which is what the derivation below actually rests on:

    B1,1-1,C1,3-3,D1,1-1,E1,2-2,F1,2-2,G1,5-5,H1,34-44,/0,A29-38,A88-91,A136-321,A360-423

Walking it: each `X1` is a single-residue pin chain, each `n-n` a free stretch. ⇒ pins land at
binder positions **1, 3, 7, 9, 12, 15, 21**, then a variable 34-44 tail.
⭐ **That is checked, not assumed:** 21 fixed/spacer positions + a 34-44 tail predicts binder
lengths **55-65**, and `J628` measured exactly 55-65. The pins sit BEFORE the variable tail, so
their indices are constant across batches -- which is the only reason indexing them is safe at all
(one contig length is sampled per batch, so n is BATCHES, and a pin index is never
hardcoded blindly).

Held atoms differ per pin (from the same spec) and are shown in the object names:
    B1 TIP · C1 TIP · D1 OH,CZ,CG · E1 TIP · F1 OH,CZ,CG · G1 TIP · H1 NE,CZ,NH1,NH2

## THE REFERENCE OBJECT
`glp1_medoid` = the `j631` medoid of 21 GLP-1 co-folds (`j630` seed20) -- the SAME structure that
became the RFd3 input for the 4,608 `J636` backbones. It carries GLP-1 itself in chain P, so it
shows where a real agonist sits against the same receptor. ⛔ It is a REFERENCE, not a design.

⛔ PyMOL cautions obeyed: scene messages are **pure ASCII** -- a marker character corrupts the
heap and the .pse segfaults on load -- and there is **no `__main__` guard**, because `pymol -cq`
would then exit 0 having written nothing.

run:  pymol -cq j641_build_baker_best_pse_20261024.py
"""
import os

from pymol import cmd

SCRATCH = "/path/to/your/working/tree"  # PLACEHOLDER PATH -- not a working default; set it for your own tree.
J628 = f"{SCRATCH}/j628_unpack"
MEDOID = f"{SCRATCH}/j630_glp1_cofold_21/medoid/glp1r_medoid_F260_withGLP1.pdb"
C3DIR = "/path/to/your/working/tree"  # PLACEHOLDER PATH -- not a working default; set it for your own tree.
OUT = f"{C3DIR}/data/J641_baker_best_20261024.pse"

PINS = [1, 3, 7, 9, 12, 15, 21]
HELD = {1: "B1 TIP", 3: "C1 TIP", 7: "D1 OH,CZ,CG", 9: "E1 TIP",
        12: "F1 OH,CZ,CG", 15: "G1 TIP", 21: "H1 NE,CZ,NH1,NH2"}

# (object name, design, seed-shard, one-line why)  -- from j637's TSV
PICKS = [
    ("pass_both_1", "arm_shard08_v1_deep6_a207_5_model_2", "shard06",
     "clears BOTH surviving cuts: ipTM 0.909, iface PAE min 0.843 A"),
    ("pass_both_2", "arm_shard06_v1_deep6_a207_4_model_2", "shard07",
     "clears BOTH surviving cuts: ipTM 0.910, iface PAE min 0.853 A"),
    ("best_iptm", "arm_shard15_v1_deep6_a207_0_model_3", "shard05",
     "best ipTM of all 960: 0.924 (PAE min 1.682 A, fails that cut)"),
    ("best_pae", "arm_shard04_v1_deep6_a207_1_model_6", "shard03",
     "best iface PAE min of all 960: 0.745 A (ipTM 0.886)"),
    ("best_plddt", "arm_shard12_v1_deep6_a207_4_model_0", "shard06",
     "best complex pLDDT of all 960: 0.787 (the cut is 0.90; no fold clears it)"),
    ("worst_contrast", "arm_shard14_v1_deep6_a207_1_model_3", "shard06",
     "CONTRAST: worst ipTM 0.082, PAE min 27.98 A -- a binder that went elsewhere"),
]

cmd.set("retain_order", 1)
cmd.set("pse_export_version", 1.8)
cmd.set("cartoon_transparency", 0.0)

loaded = []
for obj, design, shard, why in PICKS:
    p = (f"{J628}/{shard}/boltz_results_{shard}/predictions/"
         f"{design}__id1/{design}__id1_model_0.pdb")
    assert os.path.exists(p), f"MISSING fold pdb: {p}"
    cmd.load(p, obj)
    # verify the layout rather than trusting it
    nb = cmd.count_atoms(f"{obj} and chain A and name CA")
    nt = cmd.count_atoms(f"{obj} and chain B and name CA")
    assert nt == 409, f"{obj}: receptor chain B has {nt} CA, expected 409"
    assert 55 <= nb <= 65, f"{obj}: binder chain A has {nb} CA, outside the contig's 55-65"
    loaded.append((obj, design, shard, why, nb))

assert os.path.exists(MEDOID), f"MISSING medoid: {MEDOID}"
cmd.load(MEDOID, "glp1_medoid")

# ---- put everything in one frame: superpose on the RECEPTOR only -------------------------------
# ⛔ align on the receptor, never on the binder -- the binder pose is the thing being looked at.
ref = loaded[0][0]
for obj, *_ in loaded[1:]:
    cmd.align(f"{obj} and chain B and name CA", f"{ref} and chain B and name CA")
cmd.align("glp1_medoid and chain A and name CA", f"{ref} and chain B and name CA")

# ---- styling -----------------------------------------------------------------------------------
cmd.hide("everything")
cmd.bg_color("white")
cmd.show("cartoon")
cmd.color("grey80", "chain B")                       # receptor on every design object
cmd.color("grey70", "glp1_medoid and chain A")       # receptor of the reference

BINDER_COLOURS = ["marine", "teal", "forest", "orange", "purple", "firebrick"]
for (obj, design, shard, why, nb), col in zip(loaded, BINDER_COLOURS):
    cmd.color(col, f"{obj} and chain A")
cmd.color("yellow", "glp1_medoid and chain P")       # GLP-1 itself

# ---- THE FIXED RESIDUES ------------------------------------------------------------------------
pinsel = "+".join(str(i) for i in PINS)
for obj, *_ in loaded:
    sel = f"{obj} and chain A and resi {pinsel}"
    n = cmd.count_atoms(f"{sel} and name CA")
    assert n == len(PINS), f"{obj}: found {n} pin CA, expected {len(PINS)}"
    cmd.show("sticks", f"{sel} and not hydro")
    cmd.color("red", f"{sel} and elem C")
    cmd.show("spheres", f"{sel} and name CA")
    cmd.set("sphere_scale", 0.30, f"{sel} and name CA")
    # a named selection per pin so each is clickable in the GUI, with its held atoms in the name
    for r in PINS:
        tag = HELD[r].split()[0]
        cmd.select(f"{obj}_pin{r}_{tag}", f"{obj} and chain A and resi {r}")
    # receptor residues actually contacted by the pins, MEASURED per structure
    cmd.select(f"{obj}_pin_contacts",
               f"byres ({obj} and chain B within 5.0 of ({sel} and not hydro))")
    cmd.show("sticks", f"{obj}_pin_contacts and not hydro and not name N+C+O")
    cmd.color("salmon", f"{obj}_pin_contacts and elem C")
cmd.deselect()

# ---- scenes. ASCII ONLY in every message -------------------------------------------------------
cmd.set("ray_opaque_background", 1)

def only(obj):
    cmd.disable("all")
    cmd.enable(obj)

cmd.orient(f"{ref} and chain B")
for obj, design, shard, why, nb in loaded:
    only(obj)
    cmd.orient(f"{obj} and chain A")
    cmd.zoom(f"{obj} and chain A", 8)
    msg = (f"{obj}  |  {design}  seed {shard}  |  binder {nb} aa\n"
           f"{why}\n"
           f"RED sticks/spheres = the 7 FIXED pin residues, binder positions "
           f"{', '.join(str(i) for i in PINS)}\n"
           f"held atoms: " + " . ".join(HELD[r] for r in PINS) + "\n"
           f"SALMON = receptor residues within 5.0 A of a pin (measured, not designed)")
    assert msg.isascii(), "scene message must be pure ASCII"
    cmd.scene(obj, "store", message=msg)

# reference scene: where a real agonist actually sits
only("glp1_medoid")
cmd.orient("glp1_medoid")
ref_msg = ("REFERENCE, not a design: j631 medoid of 21 GLP-1 co-folds (j630 seed20).\n"
           "This is the same structure used as the RFd3 input for the 4,608 J636 backbones.\n"
           "YELLOW = GLP-1 itself, bound to the GREY receptor. Compare where our binders landed.")
assert ref_msg.isascii(), "scene message must be pure ASCII"
cmd.scene("ref_glp1", "store", message=ref_msg)

# all-designs overlay
cmd.disable("all")
for obj, *_ in loaded:
    cmd.enable(obj)
cmd.enable("glp1_medoid")
cmd.hide("cartoon", "glp1_medoid and chain A")      # one receptor copy is enough
cmd.orient(f"{ref} and chain B")
all_msg = ("ALL picks superposed on the RECEPTOR (grey). Each binder a different colour;\n"
           "GLP-1 in yellow for reference. RED = the 7 fixed pin residues on every binder.\n"
           "The spread between coloured binders is the pose disagreement J638 is about.")
assert all_msg.isascii(), "scene message must be pure ASCII"
cmd.scene("all_overlay", "store", message=all_msg)

cmd.scene("all_overlay", "recall")
cmd.save(OUT)
print(f"wrote {OUT}")
print(f"objects: {len(loaded)} designs + 1 reference; scenes: {len(loaded) + 2}")
print(f"fixed residues marked at binder positions {PINS}")
