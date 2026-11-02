"""j957 -- render the j956 corrected N-terminal top-3 session to raw PNGs. Run `pymol -cq`.

⛔ NO `__main__` guard: under `pymol -cq` this module is not __main__, so a guarded script
   exits 0 having written nothing. (Recorded lesson; see j902_render_screen_pse_20261101.py.)

⛔ READ-ONLY ON THE SESSION. The .pse is loaded and never saved, and no scene is re-stored,
   so there is no path by which a non-ASCII `scene(message=)` could be written here.

THE THREE RECORDED RENDER DEFECTS, AND WHAT IS DONE ABOUT EACH
--------------------------------------------------------------
1. ⛔ CLIPPED TARGET SURFACE RENDERING AS A BLACK WEDGE OVER THE INTERFACE. `orient`/`zoom`
   set the near/far planes tight to the selection and a surface crossing the near plane is
   cut open -- the defect hid the one region the figure exists to show. ⇒ the slab is widened
   explicitly after the camera is set, and `ray_interior_color` is set so any residual cut
   face reads as material rather than as a void.
   ⚠️ `ray_interior_color` is a RENDERING CHOICE, NOT DATA: it changes how a cut face looks
   and nothing about geometry. Said plainly so nobody reads a grey interior as a measurement.
2. ⛔ THE SURFACE RENDERED FULLY OPAQUE because `transparency` was set BEFORE the objects
   were loaded, so it applied to nothing. ⇒ it is set PER OBJECT AFTER loading, and re-set
   after EVERY scene recall (recall restores the stored per-object settings), and then
   ASSERTED per object immediately before each ray. A ray with an opaque surface aborts.
3. ⛔ label overflow is handled in the stamping pass, not here; raw renders carry no text.

COMPARABILITY: the three per-design cameras saved in the session differ from each other, so
using them would show the target in a different orientation in each panel. One COMMON camera
is derived once (the vetted 00_overview rotation, zoomed to hold all three binders) and is
applied to every panel, so the target is rendered identically across panels. The common frame
itself was asserted at build time (target-vs-target RMSD 0.0000 A); this script does not
re-derive it, and says so rather than implying it re-checked.

run: pymol -cq j957_render_top3_pse_20261102.py -- <in.pse> <outdir>
"""
import json
import os
import sys

from pymol import cmd

W, H = 3000, 2250            # exceeds the 2400x1800 of the previous pass
DPI = 300
TRANSP = 0.60                # target surface transparency, per object, after load.
                             # 0.50 was tried first: the grey BUILT segment lies behind the
                             # target surface in this view and read too faintly. 0.60 keeps
                             # the surface clearly present while the built segment shows
                             # through. ⚠️ a presentation choice, not a measurement.
MIN_TRANSP = 0.20            # below this the surface is effectively opaque -> refuse

# rank order is fixed and the numbers are the MEASURED values carried in from the
# scoring pass; this script does not re-measure them and does not claim to.
PANELS = [
    {"obj": "rank1_0_model_1", "rank": 1,
     "design": "arm_shard01_tfr1_nterm25_0_model_1", "batch": "s01b0", "pairs": 165},
    {"obj": "rank2_2_model_0", "rank": 2,
     "design": "arm_shard03_tfr1_nterm25_2_model_0", "batch": "s03b2", "pairs": 150},
    {"obj": "rank3_0_model_4", "rank": 3,
     "design": "arm_shard01_tfr1_nterm25_0_model_4", "batch": "s01b0", "pairs": 141},
]


def set_transparency(objs):
    """⛔ PER OBJECT, AFTER LOAD / AFTER RECALL -- defect 2."""
    for o in objs:
        cmd.set("transparency", TRANSP, o)
        cmd.set("stick_radius", 0.26, o)


def assert_transparency(objs):
    """refuse to ray an opaque surface rather than shipping defect 2 again"""
    bad = []
    for o in objs:
        v = float(cmd.get("transparency", o))
        if v < MIN_TRANSP:
            bad.append((o, v))
    if bad:
        raise RuntimeError("OPAQUE SURFACE REFUSED: %s" % bad)
    return {o: float(cmd.get("transparency", o)) for o in objs}


def main(pse, outdir, transp=None):
    # ⭐ ADDITIVE 2026-11-02: optional surface-transparency override, default unchanged at 0.60.
    # The BUILT segment sits partly INSIDE the target surface, so at 0.60 the thing the figure
    # exists to show reads faintly. Raising it is a PRESENTATION choice and changes no geometry.
    if transp is not None:
        global TRANSP
        TRANSP = float(transp)
        if TRANSP < MIN_TRANSP:
            raise SystemExit("REFUSING: transparency %.2f below MIN_TRANSP %.2f" % (TRANSP, MIN_TRANSP))
        print("== surface transparency overridden to %.2f ==" % TRANSP)
    os.makedirs(outdir, exist_ok=True)
    cmd.reinitialize()
    cmd.load(pse)
    objs = cmd.get_object_list()
    print("== loaded %s: %d objects: %s" % (os.path.basename(pse), len(objs), objs))

    # ---- global render quality
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)      # needed so the white-only autocrop is valid
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.2)
    cmd.set("ray_interior_color", "grey70")  # ⚠️ cosmetic only -- see the docstring
    cmd.set("surface_quality", 1)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("orthoscopic", 1)                # comparable panels: no perspective difference
    set_transparency(objs)                   # ⛔ defect 2: per object, after load

    # ---- establish the reps the scenes carry (recall sets surface/cartoon/sticks per object)
    for sc in cmd.get_scene_list():
        cmd.scene(sc, "recall")
    set_transparency(objs)

    # ---- ONE COMMON CAMERA for every panel, from the vetted overview rotation
    cmd.scene("00_overview", "recall")
    cmd.enable("all")
    cmd.zoom("all", buffer=6.0)
    cmd.clip("slab", 400)                    # ⛔ defect 1: widen the slab
    common = cmd.get_view()
    print("== common view (applied to all panels) ==\n   %s"
          % ",".join("%.4f" % x for x in common))

    def shoot(path, note):
        cmd.set_view(common)
        cmd.clip("slab", 400)                # ⛔ defect 1, again after any camera change
        tr = assert_transparency(objs)       # ⛔ defect 2, asserted not assumed
        cmd.ray(W, H)
        cmd.png(path, dpi=DPI)
        ok = os.path.exists(path) and os.path.getsize(path) > 0
        print("   %-26s %s  %dx%d  %d bytes  transparency=%s"
              % (os.path.basename(path), "OK" if ok else "FAILED", W, H,
                 os.path.getsize(path) if ok else 0,
                 {k: round(v, 2) for k, v in tr.items()}))
        if not ok:
            raise RuntimeError("no PNG written: %s (%s)" % (path, note))

    written = []

    # ---- one panel per design: that object only, so exactly one target surface is shown
    for p in PANELS:
        o = p["obj"]
        cmd.scene(o, "recall")               # colours + reps as built, never recoloured here
        set_transparency(objs)               # ⛔ recall restores stored settings -> re-set
        cmd.disable("all")
        cmd.enable(o)
        cmd.show("surface", "%s and chain B+C" % o)
        png = os.path.join(outdir, "j957_panel_rank%d_20261102.png" % p["rank"])
        shoot(png, o)
        rec = dict(p)
        rec["png"] = os.path.basename(png)
        rec["kind"] = "design"
        rec["n_held"] = cmd.count_atoms("%s_held and name CA" % o)
        rec["n_built"] = cmd.count_atoms("%s_built and name CA" % o)
        rec["n_pin_atoms"] = cmd.count_atoms("%s_pin" % o)
        written.append(rec)

    # ---- overview: all three binders against ONE copy of the target surface
    cmd.scene("00_overview", "recall")
    set_transparency(objs)
    cmd.enable("all")
    keep = PANELS[0]["obj"]
    cmd.show("surface", "%s and chain B+C" % keep)
    for p in PANELS[1:]:
        cmd.hide("surface", p["obj"])
    png = os.path.join(outdir, "j957_panel_overview_20261102.png")
    shoot(png, "overview")
    written.append({"png": os.path.basename(png), "kind": "overview", "rank": 0,
                    "obj": "", "design": "", "batch": "", "pairs": "",
                    "surface_from": keep})

    man = {"pse": os.path.abspath(pse), "w": W, "h": H, "dpi": DPI,
           "transparency": TRANSP, "orthoscopic": 1, "slab": 400,
           "ray_interior_color": "grey70 (cosmetic, not data)",
           "antialias": 2, "ray_shadows": 0, "specular": 0.2,
           "common_view": [round(x, 5) for x in common],
           "panels": written}
    mp = os.path.join(outdir, "j957_renders_20261102.json")
    json.dump(man, open(mp, "w"), indent=1)
    print("\n   %d raw PNGs + manifest -> %s" % (len(written), outdir))


argv = sys.argv[1:]
if "--" in argv:
    argv = argv[argv.index("--") + 1:]
main(*argv[:3])
