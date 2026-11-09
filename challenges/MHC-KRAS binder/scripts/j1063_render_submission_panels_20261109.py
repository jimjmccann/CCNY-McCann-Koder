"""j1063 stage 2/4 -- render the submission .pse to RAW PNGs: one panel per design + an overview.

    pymol -cq j1063_render_submission_panels_20261109.py -- <build.json> <outdir>

⭐ METHOD AND CREDIT: `an earlier challenge's scripts/j531_render_submission_panel_20261018.py` (the
   submission-panel case, per-design colouring) and `TfR1 binder/scripts/
   j957_render_top3_pse_20261102.py` (the one-common-camera + recorded-defect discipline).
   Neither is reusable as-is -- `j531` renders ONE grid image via `grid_mode` and hard-codes a
   TNF epitope, `j957` hard-codes a three-design TfR1 cohort -- so the DISCIPLINE is carried over
   and the cohort comes from stage 1's manifest.

⛔ NO `__main__` GUARD: under `pymol -cq` this module is not `__main__`, so a guarded script exits
   0 having written NOTHING [[pymol-script-not-main]].
⛔ NO `__file__`: PyMOL rebinds it to its own site-packages dir. Every path comes from argv.
⛔ READ-ONLY ON THE SESSION: the .pse is loaded, never saved, and no scene is re-stored, so there
   is no path by which a non-ASCII `scene(message=)` could be written here.
⛔ retain_order is DELIBERATELY NOT SET (it corrupts `cmd.align`); nothing is aligned here anyway.

THE RECORDED RENDER DEFECTS, CARRIED OVER, AND WHAT IS DONE ABOUT EACH
---------------------------------------------------------------------
1. ⛔⛔ A CLIPPED OBJECT RENDERS AS A BLACK WEDGE OVER THE INTERFACE -- i.e. the defect hid the
   one region the figure exists to show. `orient`/`zoom` set the near and far planes tight to the
   selection. ⇒ the slab is widened EXPLICITLY after the camera is set and again after EVERY
   camera change, and `ray_interior_color` makes any residual cut face read as material rather
   than as a void.
   ⚠️ `ray_interior_color` is a RENDERING CHOICE, NOT DATA: it changes how a cut face looks and
   nothing about geometry. Said plainly so nobody reads a grey interior as a measurement.
2. ⛔ TRANSPARENCY APPLIED TO NOTHING because it was set BEFORE the objects existed. ⇒ set PER
   OBJECT after load, re-set after EVERY scene recall (recall restores the stored per-object
   settings), and ASSERTED immediately before each ray.
   ⛔⛔ RE-POINTED, NOT COPIED: `j957` guarded `transparency`, which is the SURFACE setting. This
   session draws CARTOONS AND STICKS AND NO SURFACE AT ALL, so a surface guard here would assert
   a property of zero objects and pass vacuously -- the exact shape of a dead tripwire. ⇒ this
   guards `cartoon_transparency` on the object stage 1 actually set it on, and reports the
   measured value with every panel.
3. ⛔ Label overflow is handled in stage 3, not here; raw renders carry NO text. ⭐ A label also
   cannot live in the image area at all: an automatic labeller that strips the common prefix and
   suffix cannot carry a constant marker, and the TIER WORD on these panels is exactly such a
   constant [[figure-label-cannot-carry-a-constant-marker]].

⛔⛔ AND A CHECK NO COUNT CAN REPLACE. `J1024` section 1: twelve .pse shipped showing a bare
   target because a DISABLED GROUP HIDES ITS ENABLED MEMBERS, and objects 9/9 and scenes 6/6 both
   passed. So before every ray this asserts that the panel's own binder object is enabled AND
   carries atoms, and it counts the non-white pixels of the finished PNG -- but ⛔ even that is
   not enough, and the last step of this job is a human opening every image.

COMPARABILITY: ONE COMMON CAMERA, derived once from the `00_overview` scene and applied to every
panel, so the target is rendered identically across panels and the panels are pixel-comparable.
⛔ The per-scene cameras stored in the session differ; using them would show the target in a
different orientation in each panel.

⛔ NUMBERS ARE CARRIED, NEVER RE-MEASURED. Every figure in the manifest this writes came from
stage 1, which asserted it against the scored tables. This script measures no science.
"""
import hashlib
import json
import os
import sys

from pymol import cmd

# ⛔ 400 dpi: at lower dpi a stamped panel label is unreadable [[read-printed-labels-not-bar-lengths]].
# ⛔⛔ AND THE CANVAS IS PORTRAIT ON PURPOSE. At 3000x2250 landscape this complex -- a pMHC with
#    the binder above the groove -- cropped honestly to 1097x1525, i.e. it used barely a third of
#    the frame width. The stamped band then had to wrap into 1100-1470 px against a 1525 px
#    picture: every label passed its re-measurement and the figure was still mostly text. Matching
#    the canvas to the subject's own 0.72 aspect is what fixes that, NOT shrinking the type
#    towards the refusal floor.
W, H, DPI = 2400, 3200, 400
SLAB = 400                      # ⛔ defect 1
MAX_CARTOON_TRANSP = 0.90       # above this the target is invisible -> refuse
MIN_INK_FRAC = 0.01             # a panel with almost no ink is a blank render dressed as success

argv = [a for a in sys.argv[1:] if a != "--" and not a.endswith(".py")]
if len(argv) != 2:
    sys.exit(__doc__)
BUILD, OUTDIR = argv

man = json.load(open(BUILD))
if not man.get("ok"):
    sys.exit("REFUSING: the stage-1 manifest does not carry ok=true")
picks = man["designs"]
TARGET, PEPTIDE = man["target_object"], man["peptide_object"]


def palette_fingerprint(picks):
    """⭐ MUST equal stage 1's. The palette travels between two programs; this is what catches a
    divergence instead of two figures that silently disagree about which design is which."""
    payload = ";".join("%d=%s:%s" % (d["slot"], d["color_name"], d["color_hex"]) for d in picks)
    return hashlib.md5(payload.encode("ascii")).hexdigest()[:8]


fp = palette_fingerprint(picks)
if fp != man["palette_fingerprint"]:
    sys.exit("REFUSING: palette fingerprint %s != stage 1's %s -- the two stages disagree about "
             "which design is which colour." % (fp, man["palette_fingerprint"]))

os.makedirs(OUTDIR, exist_ok=True)
cmd.reinitialize()
cmd.load(man["pse"])
objs = cmd.get_object_list()
print("== j1063 stage 2/4: render %d panels + overview ==" % len(picks))
print("   loaded %s -- %d objects" % (os.path.basename(man["pse"]), len(objs)))
print("   palette fingerprint %s (matches stage 1)" % fp)

# ---- global render quality
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)      # the white-only autocrop in stage 3 depends on this
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.set("specular", 0.2)
cmd.set("ray_interior_color", "grey70")  # ⚠️ cosmetic only -- see the docstring
cmd.set("cartoon_fancy_helices", 1)
cmd.set("orthoscopic", 1)                # comparable panels: no perspective difference


def set_transp():
    """⛔ defect 2: PER OBJECT, after load AND after every scene recall."""
    cmd.set("cartoon_transparency", man["target_transparency"], TARGET)
    cmd.set("cartoon_transparency", 0.0, PEPTIDE)
    for d in picks:
        cmd.set("cartoon_transparency", 0.0, d["obj"])


def assert_transp():
    """⛔ defect 2 ASSERTED, not assumed. Refuses an invisible target or a faded binder."""
    got = {TARGET: float(cmd.get("cartoon_transparency", TARGET))}
    if got[TARGET] > MAX_CARTOON_TRANSP:
        raise RuntimeError("TARGET EFFECTIVELY INVISIBLE, REFUSED: %s" % got)
    for d in picks:
        v = float(cmd.get("cartoon_transparency", d["obj"]))
        if v > 0.0:
            raise RuntimeError("BINDER %s is %.2f transparent; the binder is the SUBJECT of its "
                               "own panel and must be fully opaque" % (d["obj"], v))
    return got


set_transp()
scenes = cmd.get_scene_list()
if "00_overview" not in scenes:
    sys.exit("REFUSING: no 00_overview scene; scenes are %s" % scenes)
for sc in scenes:                        # establish the reps the scenes carry
    cmd.scene(sc, "recall")
set_transp()                             # ⛔ recall overwrote them

# ---- ONE COMMON CAMERA, from the overview scene
cmd.scene("00_overview", "recall")
for d in picks:
    cmd.enable(d["obj"])
cmd.enable(TARGET)
cmd.enable(PEPTIDE)
cmd.zoom("all", buffer=6.0)
cmd.clip("slab", SLAB)                   # ⛔ defect 1
COMMON = cmd.get_view()
print("   common view applied to every panel:\n     %s"
      % ",".join("%.4f" % x for x in COMMON))

# ---- the target must be the SAME OBJECT IN THE SAME FRAME in every panel. It is literally one
#      object (stage 1 keeps a single copy), so this records its identity, its atom count and its
#      centre of mass ONCE and re-asserts them before every single ray. ⛔ Asserted, not assumed.
TARGET_N = cmd.count_atoms(TARGET)
PEP_N = cmd.count_atoms(PEPTIDE)
TARGET_COM = [round(v, 4) for v in cmd.centerofmass(TARGET)]
PEP_COM = [round(v, 4) for v in cmd.centerofmass(PEPTIDE)]
print("   target identity locked: %s %d atoms com %s | %s %d atoms com %s"
      % (TARGET, TARGET_N, TARGET_COM, PEPTIDE, PEP_N, PEP_COM))


def assert_target_identity(where):
    if cmd.count_atoms(TARGET) != TARGET_N or cmd.count_atoms(PEPTIDE) != PEP_N:
        raise RuntimeError("target atom count changed before %s" % where)
    if [round(v, 4) for v in cmd.centerofmass(TARGET)] != TARGET_COM:
        raise RuntimeError("target moved before %s -- the panels are no longer comparable" % where)
    if [round(v, 4) for v in cmd.centerofmass(PEPTIDE)] != PEP_COM:
        raise RuntimeError("peptide moved before %s" % where)
    for o in (TARGET, PEPTIDE):
        if o not in cmd.get_names("all", enabled_only=1):
            raise RuntimeError("%s is NOT ENABLED for %s" % (o, where))


def ink_fraction(path):
    """Count the non-white pixels of the FINISHED file. ⛔ A clean empty report is a failure mode:
    a ray that produced a white rectangle exits 0 and writes a perfectly valid PNG."""
    try:
        from PIL import Image
    except ImportError:
        return None
    im = Image.open(path).convert("L").resize((300, 225))
    px = im.getdata()
    return sum(1 for v in px if v < 250) / float(len(px))


panels = []


def shoot(scene_name, png, kind, d):
    cmd.scene(scene_name, "recall")      # colours and reps exactly as built; never recoloured here
    set_transp()                         # ⛔ defect 2, after the recall
    cmd.disable("designs")               # clear, then enable exactly what this panel shows
    for x in picks:
        cmd.disable(x["obj"])
    cmd.enable(TARGET)
    cmd.enable(PEPTIDE)
    cmd.enable("target")
    cmd.enable("designs")                # ⛔⛔ the GROUP too: a disabled group hides its members
    shown = [d] if d else picks
    for x in shown:
        cmd.enable(x["obj"])
        if cmd.count_atoms(x["obj"]) == 0:
            raise RuntimeError("%s has zero atoms" % x["obj"])
        if x["obj"] not in cmd.get_names("all", enabled_only=1):
            raise RuntimeError("%s is NOT ENABLED -- this is the bare-target signature" % x["obj"])
    assert_target_identity(png)
    cmd.set_view(COMMON)
    cmd.clip("slab", SLAB)               # ⛔ defect 1, again after the camera change
    tr = assert_transp()                 # ⛔ defect 2, asserted
    path = os.path.join(OUTDIR, png)
    cmd.ray(W, H)
    cmd.png(path, dpi=DPI)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise RuntimeError("no PNG written: %s" % path)
    ink = ink_fraction(path)
    if ink is not None and ink < MIN_INK_FRAC:
        raise RuntimeError("%s is %.4f ink -- effectively blank, refusing to call it a panel"
                           % (png, ink))
    print("   %-34s OK %dx%d %7d B  ink %.3f  enabled %d binder(s)  target_transp %.2f"
          % (png, W, H, os.path.getsize(path), ink if ink is not None else -1,
             len(shown), tr[TARGET]))
    ent = {"png": png, "kind": kind, "scene": scene_name, "ink_fraction": round(ink, 4),
           "n_binders_shown": len(shown)}
    if d:
        ent.update({k: d[k] for k in
                    ("slot", "design", "sub_name", "origin_round", "tag", "shard", "rank",
                     "arm", "tier",
                     "seeds_clear", "seeds_total", "margin_mean", "iptm_mean",
                     "sd_pep_contacts", "failed_limb", "color_name", "color_hex", "binder_len",
                     "fit_rms", "fit_n")})
    panels.append(ent)


shoot("00_overview", "j1063_panel_00_overview.png", "overview", None)
for d in picks:
    shoot("%02d_slot%d" % (d["slot"], d["slot"]),
          "j1063_panel_%02d_%s.png" % (d["slot"], d["design"]), "design", d)

# ⛔ THE VERDICT IS A POSITIVE MARKER IN A FILE, never absent error text on stdout.
out = {"ok": True, "build": os.path.abspath(BUILD), "pse": man["pse"],
       "palette_fingerprint": fp, "W": W, "H": H, "dpi": DPI, "slab": SLAB,
       "common_view": [round(x, 5) for x in COMMON],
       "target_object": TARGET, "peptide_object": PEPTIDE,
       "target_atoms": TARGET_N, "target_com": TARGET_COM,
       "peptide_atoms": PEP_N, "peptide_com": PEP_COM,
       "target_transparency": man["target_transparency"],
       "ray_interior_color": "grey70 (cosmetic, NOT data)",
       "hla_len": man["hla_len"], "min_fit_atoms": man["min_fit_atoms"],
       "panels": panels}
mp = os.path.join(OUTDIR, "j1063_renders_20261109.json")
json.dump(out, open(mp, "w"), indent=1)
print("   wrote %s" % mp)
print("== j1063 stage 2 ok: %d raw PNGs ==" % len(panels))

# ⛔ os._exit skips PyMOL's crash-prone atexit teardown -- ⛔⛔ and ALSO Python's stdout flush,
#    so without this flush every line above is lost on a pipe and the run "succeeds" silently.
sys.stdout.flush()
os._exit(0)
