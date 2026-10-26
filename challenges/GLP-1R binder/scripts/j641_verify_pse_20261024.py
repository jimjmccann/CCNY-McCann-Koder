#!/usr/bin/env python3
"""j641 verify -- reload the .pse and recall EVERY scene.

⛔ The failure this guards against is not hypothetical: a non-ASCII character in a
`scene(message=)` corrupts the heap and the .pse SEGFAULTS when opened. Writing the file proves
nothing; only a reload does.
⛔ No `__main__` guard -- `pymol -cq` on a guarded script exits 0 having written nothing.
"""
from pymol import cmd

PSE = "/path/to/your/working/tree/data/J641_baker_best_20261024.pse"  # PLACEHOLDER PATH -- set it for your own tree.
PINS = [1, 3, 7, 9, 12, 15, 21]

cmd.load(PSE)
objs = cmd.get_object_list()
scenes = cmd.get_scene_list()
print(f"reloaded OK: {len(objs)} objects, {len(scenes)} scenes")
print(f"    objects: {objs}")
print(f"    scenes : {scenes}")

for s in scenes:
    cmd.scene(s, "recall")
print(f"all {len(scenes)} scenes recalled without crashing")

bad = 0
for o in objs:
    if o == "glp1_medoid":
        nt = cmd.count_atoms(f"{o} and chain A and name CA")
        npep = cmd.count_atoms(f"{o} and chain P and name CA")
        print(f"    {o:<16} receptor {nt} CA, GLP-1 {npep} CA")
        continue
    nb = cmd.count_atoms(f"{o} and chain A and name CA")
    nt = cmd.count_atoms(f"{o} and chain B and name CA")
    sel = f"{o} and chain A and resi {'+'.join(str(i) for i in PINS)}"
    npin = cmd.count_atoms(f"{sel} and name CA")
    nred = cmd.count_atoms(f"{sel} and elem C")
    print(f"    {o:<16} binder {nb:>3} CA, receptor {nt} CA, pins {npin}/7, pin C atoms {nred}")
    if npin != 7:
        bad += 1

assert bad == 0, f"{bad} object(s) do not carry all 7 fixed residues"
print("OK -- every design object carries all 7 fixed residues")
