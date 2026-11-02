"""j957 probe -- read-only inspection of the j956 session. Run under `pymol -cq`.

Reports objects, scenes, chains, atom/residue counts, saved selections and which
representations are enabled, so the render script is written against the session that
actually exists rather than against its description.

No __main__ guard on purpose: under `pymol -cq` this module is not __main__ and a guarded
script exits 0 having written nothing.
Loads the session read-only; saves nothing, so the .pse cannot be modified.

run: pymol -cq j957_probe_pse_20261102.py -- <in.pse>
"""
import sys

from pymol import cmd

argv = sys.argv[1:]
if "--" in argv:
    argv = argv[argv.index("--") + 1:]
pse = argv[0]

cmd.reinitialize()
cmd.load(pse)

print("== objects ==")
for o in cmd.get_object_list():
    print("  %s  atoms=%d" % (o, cmd.count_atoms(o)))
    for ch in cmd.get_chains(o):
        n_at = cmd.count_atoms("%s and chain %s" % (o, ch))
        n_ca = cmd.count_atoms("%s and chain %s and name CA" % (o, ch))
        print("      chain %s: atoms=%d CA=%d" % (ch, n_at, n_ca))

print("== scenes ==")
print("  %s" % cmd.get_scene_list())

print("== named selections (get_names of type selection) ==")
for s in cmd.get_names("selections", enabled_only=0):
    print("  %s  atoms=%d" % (s, cmd.count_atoms(s)))

print("== per-scene enabled reps ==")
for sc in cmd.get_scene_list():
    cmd.scene(sc, "recall")
    on = [o for o in cmd.get_object_list() if cmd.get_object_state(o) is not None]
    vis = []
    for o in cmd.get_object_list():
        reps = []
        for rep in ("cartoon", "sticks", "surface", "lines", "spheres", "ribbon",
                    "nonbonded", "mesh", "dots"):
            if cmd.count_atoms("%s and rep %s" % (o, rep)) > 0:
                reps.append("%s:%d" % (rep, cmd.count_atoms("%s and rep %s" % (o, rep))))
        vis.append("    %s enabled=%s reps=%s" % (o, o in cmd.get_names("objects", enabled_only=1),
                                                  ",".join(reps) if reps else "none"))
    print("  scene %s" % sc)
    for v in vis:
        print(v)
    print("    view = %s" % (",".join("%.3f" % x for x in cmd.get_view()),))

print("== transparency / global settings of interest ==")
for o in cmd.get_object_list():
    print("  %s transparency=%s cartoon_transparency=%s stick_radius=%s" % (
        o, cmd.get("transparency", o), cmd.get("cartoon_transparency", o),
        cmd.get("stick_radius", o)))
for k in ("bg_rgb", "ray_opaque_background", "antialias", "ray_shadows", "specular",
          "ray_interior_color", "surface_quality", "cartoon_side_chain_helper", "orthoscopic"):
    print("  global %s = %s" % (k, cmd.get(k)))

print("== residue ranges per object/chain ==")
for o in cmd.get_object_list():
    for ch in cmd.get_chains(o):
        resis = []
        cmd.iterate("%s and chain %s and name CA" % (o, ch),
                    "resis.append(int(resi))", space={"resis": resis})
        if resis:
            print("  %s chain %s: %d CA, resi %d..%d" % (o, ch, len(resis), min(resis), max(resis)))

print("== PROBE DONE ==")
