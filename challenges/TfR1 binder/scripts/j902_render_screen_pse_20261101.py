"""j902 -- render the j901 session to PNGs. Run under `pymol -cq`.

⛔ NO `__main__` guard: under `pymol -cq` this module is not __main__, so a guarded script
   exits 0 having written nothing.

TWO DEFECTS IN THE FIRST RENDER PASS, BOTH FIXED HERE AND BOTH VISIBLE IN THE IMAGE
------------------------------------------------------------------------------------
1. ⛔ THE TARGET SURFACE WAS CLIPPED, and the clipped interior rendered as a BLACK WEDGE sitting
   exactly over the binding interface -- i.e. the defect hid the one region the figure exists to
   show. `cmd.orient` sets near/far planes tight to the selection, and a surface crossing the near
   plane is cut open. ⇒ the slab is widened explicitly after every orient, and
   `ray_interior_color` is set so any residual cut face reads as material rather than a void.
2. ⛔ THE SURFACE RENDERED OPAQUE. `transparency` was set BEFORE the objects were loaded, so it
   applied to nothing. ⇒ it is set PER OBJECT after loading, which is also what scene recall
   restores.

⚠️ `ray_interior_color` is a rendering choice, not data: it changes how a cut face LOOKS and
nothing about geometry. Said plainly so nobody reads a grey interior as a measured surface.

run: pymol -cq j902_render_screen_pse_20261101.py -- <in.pse> <manifest.json> <outdir>
"""
import json
import os
import sys

from pymol import cmd

W, H = 2400, 1800


def main(pse, manp, outdir):
    man = json.load(open(manp))
    os.makedirs(outdir, exist_ok=True)
    cmd.reinitialize()
    cmd.load(pse)

    # ---- global render quality
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.2)
    cmd.set("ray_interior_color", "grey70")   # ⚠️ cosmetic only -- see the docstring
    cmd.set("surface_quality", 1)
    cmd.set("cartoon_side_chain_helper", 1)

    objs = cmd.get_object_list()
    for o in objs:
        cmd.set("transparency", 0.45, o)      # ⛔ PER OBJECT, after load
        cmd.set("stick_radius", 0.24, o)

    def shoot(path, focus):
        cmd.orient(focus)
        cmd.zoom(focus, buffer=4.0)
        # ⛔ widen the slab or the surface is cut open across the interface
        cmd.clip("slab", 400)
        cmd.ray(W, H)
        cmd.png(path, dpi=300)
        print("   wrote %s" % os.path.basename(path))

    written = []
    for m in man:
        obj = "%s_b%s_m%s" % (m["opt"], m["design"].split("_")[-3],
                              m["design"].split("_")[-1])
        cmd.scene(obj, "recall")
        cmd.disable("all")
        cmd.enable(obj)
        for o in objs:
            cmd.set("transparency", 0.45, o)
        p = os.path.join(outdir, "%s.png" % obj)
        shoot(p, "%s and chain A" % obj)
        written.append((p, obj, m))

    # ---- overview: all six binders, one target surface
    cmd.scene("00_overview", "recall")
    cmd.enable("all")
    first = "%s_b%s_m%s" % (man[0]["opt"], man[0]["design"].split("_")[-3],
                            man[0]["design"].split("_")[-1])
    for m in man[1:]:
        o = "%s_b%s_m%s" % (m["opt"], m["design"].split("_")[-3],
                            m["design"].split("_")[-1])
        cmd.hide("surface", o)
    cmd.set("transparency", 0.5, first)
    p = os.path.join(outdir, "00_overview.png")
    shoot(p, "all")
    written.append((p, "00_overview", None))

    json.dump([{"png": os.path.basename(p), "obj": o,
                "design": (m or {}).get("design", ""),
                "opt": (m or {}).get("opt", ""),
                "why": (m or {}).get("why", ""),
                "max_step": (m or {}).get("max_step", ""),
                "bb_min": (m or {}).get("bb_min", ""),
                "min_cut": (m or {}).get("min_cut", ""),
                "rg_ratio": (m or {}).get("rg_ratio", ""),
                "n_under_2p5": (m or {}).get("n_under_2p5", ""),
                "n_held": len((m or {}).get("held", [])),
                "n_built": len((m or {}).get("built", [])),
                "n_res": (m or {}).get("n_res", "")}
               for p, o, m in written],
              open(os.path.join(outdir, "j902_renders.json"), "w"), indent=1)
    print("\n   %d PNGs -> %s" % (len(written), outdir))


argv = sys.argv[1:]
if "--" in argv:
    argv = argv[argv.index("--") + 1:]
main(argv[0], argv[1], argv[2])
