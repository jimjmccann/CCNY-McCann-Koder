"""j901 builder -- the .pse of the best screen designs. Run under `pymol -cq`.

⛔ NO `if __name__ == "__main__"` GUARD. Under `pymol -cq` the module is NOT __main__, so a
guarded script exits 0 having written NOTHING.

⛔⛔ EVERY `scene(message=...)` STRING IS PURE ASCII. A non-ASCII mark corrupts PyMOL's heap and
the saved session SEGFAULTS on load. No arrows, no marks.

WHAT THE COLOURING SHOWS, which is the point of the session
-----------------------------------------------------------
Per design, three disjoint residue sets, measured by j901 prep rather than assumed:
  HELD   = the 5XOX core elements RFd3 was GIVEN   -> purple
  PINS   = the four fixed epitope-facing residues  -> yellow sticks
  BUILT  = everything RFd3 generated to connect them -> grey
⇒ the grey is the answer to "did RFd3 connect the core to the pins", which is the question the
screen was launched to settle.

run: pymol -cq j901_build_screen_pse_20261101.py -- <indir> <out.pse>
"""
import json
import os
import sys

from pymol import cmd

PURPLE = "0x7D55C7"      # the group's stated CCNY purple, for the GIVEN core
GREY = "0xBFBFBF"        # RFd3-BUILT connector
YELLOW = "0xF2C230"      # the four fixed pins
TEAL = "0x5FB3B3"        # TfR1 target


def main(indir, outpse):
    man = json.load(open(os.path.join(indir, "j901_manifest.json")))
    print("== j901 build: %d designs ==" % len(man))

    cmd.reinitialize()
    cmd.set("assembly", "")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("surface_quality", 1)
    cmd.set("transparency", 0.55)
    cmd.set("stick_radius", 0.22)
    cmd.set("cartoon_oval_length", 0.9)

    ref = None
    for m in man:
        obj = "%s_b%s_m%s" % (m["opt"], m["design"].split("_")[-3],
                              m["design"].split("_")[-1])
        cmd.load(os.path.join(indir, m["cif"]), obj)
        # ⛔ exclude any V pseudo-atom by ELEMENT AND NAME, never a text match. Prep measured 0
        #    in this cohort, but the builder must not depend on that staying true.
        cmd.remove("%s and elem V and name V*" % obj)

        tsel = " or ".join("chain %s" % c for c in m["target_chains"])
        cmd.hide("everything", obj)
        cmd.show("cartoon", "%s and chain A" % obj)
        cmd.show("surface", "%s and (%s)" % (obj, tsel))
        cmd.color(GREY, "%s and chain A" % obj)
        cmd.color(TEAL, "%s and (%s)" % (obj, tsel))

        held = "+".join(str(r) for r in m["held"])
        pins = "+".join(str(r) for r in m["pins"])
        cmd.color(PURPLE, "%s and chain A and resi %s" % (obj, held))
        cmd.color(YELLOW, "%s and chain A and resi %s" % (obj, pins))
        cmd.show("sticks", "%s and chain A and resi %s and not hydro" % (obj, pins))
        cmd.set("cartoon_side_chain_helper", 1, obj)

        # a named selection per design so the viewer can re-find the sets by clicking
        cmd.select("%s_held" % obj, "%s and chain A and resi %s" % (obj, held))
        cmd.select("%s_pins" % obj, "%s and chain A and resi %s" % (obj, pins))
        cmd.deselect()

        # ⭐ a COMMON FRAME is asserted, not assumed: every design carries the same TfR1 crop, so
        #    the target-on-target RMSD to the first object must be ~0. It is PRINTED.
        if ref is None:
            ref = obj
        else:
            try:
                r = cmd.align("%s and (%s) and name CA" % (obj, tsel),
                              "%s and (%s) and name CA" % (ref, tsel),
                              cycles=0, transform=0)
                print("   target-vs-%s RMSD for %s: %.4f A (expect ~0: shared crop)"
                      % (ref, obj, r[0]))
            except Exception as exc:
                print("   NOTE could not compare target frames for %s: %s" % (obj, exc))
        print("   %-22s %-36s held %3d  pins %d  built %3d" %
              (obj, m["design"], len(m["held"]), len(m["pins"]), len(m["built"])))
        m["_obj"] = obj

    # ---- scenes. ASCII ONLY in every message.
    cmd.orient(ref)
    for m in man:
        obj = m["_obj"]
        cmd.disable("all")
        cmd.enable(obj)
        cmd.orient("%s and chain A" % obj)
        msg = ("%s  %s\n"
               "purple = 5XOX core GIVEN to RFd3 (%d res)\n"
               "grey   = RFd3 BUILT this to connect it (%d res)\n"
               "yellow = the 4 fixed pins\n"
               "teal   = TfR1 target\n"
               "max CA-CA step %s A | backbone min to target %s A\n"
               "min_cut %s (one fold if high) | rg_ratio %s | res under 2.5 A: %s\n"
               "%s"
               % (obj, m["design"], len(m["held"]), len(m["built"]),
                  m["max_step"], m["bb_min"], m["min_cut"], m["rg_ratio"],
                  m["n_under_2p5"], m["why"]))
        assert msg.isascii(), "non-ASCII scene message for %s" % obj
        cmd.scene(obj, "store", message=msg)

    # overview: all binders, target once
    cmd.enable("all")
    for m in man[1:]:
        tsel = " or ".join("chain %s" % c for c in m["target_chains"])
        cmd.hide("surface", "%s and (%s)" % (m["_obj"], tsel))
    cmd.orient(ref)
    ov = ("ALL %d designs, one TfR1 surface shown\n"
          "purple = given 5XOX core   grey = RFd3-built connector   yellow = pins\n"
          "Six designs from SIX DISTINCT BATCHES across arms opt05, opt01, opt04.\n"
          "Nothing here is scored for binding: no co-fold, no sequence design." % len(man))
    assert ov.isascii()
    cmd.scene("00_overview", "store", message=ov)

    cmd.scene("00_overview", "recall")
    cmd.save(outpse)
    print("\n   wrote %s (%d B)" % (outpse, os.path.getsize(outpse)))
    print("   scenes: %s" % ",".join(cmd.get_scene_list()))


argv = sys.argv[1:]
if "--" in argv:
    argv = argv[argv.index("--") + 1:]
main(argv[0], argv[1])
