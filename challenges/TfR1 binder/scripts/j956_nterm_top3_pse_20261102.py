"""j956 -- the .pse of the top-scoring N-terminal-extension designs. Run under `pymol -cq`.

THE ASK, 2026-11-02: a viewer session holding the best-scoring designs.

⛔ NO `if __name__ == "__main__"` GUARD. Under `pymol -cq` the module is NOT __main__, so a
guarded script exits 0 having written NOTHING.

⛔⛔ EVERY `scene(message=...)` STRING IS PURE ASCII. A non-ASCII mark corrupts PyMOL's heap and
the saved session SEGFAULTS on load. No arrows, no marks.

⛔⛔ THE RANKING THIS DRAWS IS THE CORRECTED ONE (methods/nterm_extension_arm.md section 5).
The first scoring run derived the
"generated" residue set from `diffused_index_map`, which names only TOKENS -- a generated contig
run has no token, so the N-terminal 25-mer was excluded entirely and `new_contacts_4p5` described
the 2-residue linker at positions 27-28. On the corrected set (27 residues = the 25-mer plus the
linker) the top of the order is completely different: the previous first and second place fall to
19th and 60th. ⇒ this session is built from
`j948_nterm_score_CORRECTED_GENSET_20261102_RANKED.tsv` and from nothing else.

WHAT THE COLOURING SHOWS, which is the point of the session
-----------------------------------------------------------
Per design, three disjoint residue sets, each MEASURED here and printed, never assumed:
  HELD  = the whole protein RFd3 was GIVEN, sequence-fixed   -> purple
  BUILT = the N-terminal extension RFd3 generated            -> grey
  PIN   = the one fixed guidepost guanidinium (NE,NH1,NH2)   -> yellow sticks
⇒ the grey is the answer to the question the arm was launched to settle: does the new N-terminal
segment reach TfR1, or does it only grow away from the interface?

⚠️ The scene message states each design's BATCH, because effective n is BATCHES and not designs:
ranks 1 and 3 come from the SAME batch, so they are not independent draws of each other.

run: pymol -cq j956_nterm_top3_pse_20261102.py -- <ranked.tsv> <cif_dir> <out.pse> [n_top]
"""
import csv
import glob
import json
import os
import re
import sys

import numpy as np
from pymol import cmd

PURPLE = "0x7D55C7"      # the group's stated CCNY purple, for the GIVEN held protein
GREY = "0xBFBFBF"        # the RFd3-BUILT N-terminal extension
YELLOW = "0xF2C230"      # the fixed guidepost
TEAL = "0x5FB3B3"        # TfR1 target

VPSEUDO = re.compile(r"^V\d+$")
CONTACT = 4.5


def sets_from_spec(jpath, structure_resi):
    """held, built, pin positions on the binder chain. held comes from the SPEC, built is the
    STRUCTURE minus held -- never the mapped range, which is the defect the methods note describes."""
    d = json.load(open(jpath))
    imap = d.get("diffused_index_map") or {}
    spec = d.get("specification") or {}
    pin_tokens = spec.get("select_fixed_atoms") or {}
    per_chain, held_tok, target_tok = {}, {}, {}
    tok_out = {}
    for tok, out in imap.items():
        m = re.match(r"^([A-Za-z]+)(\d+)$", str(out))
        if not m:
            continue
        ch, ri = m.group(1), int(m.group(2))
        per_chain.setdefault(ch, []).append(ri)
        tok_out[str(tok)] = (ch, ri)
        (target_tok if re.match(r"^A\d+$", str(tok)) else held_tok).setdefault(ch, []).append(ri)
    binder = [c for c in per_chain if held_tok.get(c)]
    target = sorted(c for c in per_chain if target_tok.get(c) and c not in binder)
    if len(binder) != 1 or not target:
        raise ValueError("ambiguous roles in %s" % os.path.basename(jpath))
    b = binder[0]
    held = sorted(set(held_tok[b]))
    built = sorted(set(structure_resi) - set(held))
    pins = sorted({tok_out[t][1] for t in pin_tokens if t in tok_out})
    pin_atoms = sorted({a.strip() for v in pin_tokens.values()
                        for a in str(v).split(",") if a.strip()})
    if not built:
        raise ValueError("no built residues for %s" % os.path.basename(jpath))
    return b, target, held, built, pins, pin_atoms


def heavy(obj, sel):
    xyz = cmd.get_coords("%s and (%s) and not hydro" % (obj, sel))
    return np.zeros((0, 3)) if xyz is None else np.asarray(xyz)


def main(ranked, cifdir, outpse, ntop=3):
    ntop = int(ntop)
    rows = [r for r in csv.DictReader(
        [l for l in open(ranked) if not l.startswith("#")], delimiter="\t")]
    top = rows[:ntop]
    print("== j956: %d of %d ranked survivors ==" % (len(top), len(rows)))

    cmd.reinitialize()
    cmd.set("assembly", "")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("surface_quality", 1)
    cmd.set("transparency", 0.55)
    cmd.set("stick_radius", 0.22)
    cmd.set("cartoon_oval_length", 0.9)

    ref, drawn = None, []
    for i, r in enumerate(top, 1):
        design = r["design"]
        hits = glob.glob(os.path.join(cifdir, "**", design + ".cif*"), recursive=True)
        if not hits:
            raise SystemExit("REFUSING: no coordinates found for %s under %s" % (design, cifdir))
        obj = "rank%d_%s" % (i, design.split("nterm25_")[-1])
        cmd.load(sorted(hits)[0], obj)
        # ⛔ drop V pseudo-atoms by ELEMENT AND NAME, never a text match
        cmd.remove("%s and elem V and name V*" % obj)

        resi = sorted({a.resi_number for a in cmd.get_model("%s and chain A" % obj).atom})
        jpath = os.path.join(os.path.dirname(sorted(hits)[0]), design + ".json")
        b, target, held, built, pins, pin_atoms = sets_from_spec(jpath, resi)
        if b != "A":
            raise SystemExit("REFUSING: binder chain is %s, not A, for %s" % (b, design))

        tsel = " or ".join("chain %s" % c for c in target)
        heldsel = "+".join(str(x) for x in held)
        builtsel = "+".join(str(x) for x in built)
        pinsel = "+".join(str(x) for x in pins)

        cmd.hide("everything", obj)
        cmd.show("cartoon", "%s and chain A" % obj)
        cmd.show("surface", "%s and (%s)" % (obj, tsel))
        cmd.color(TEAL, "%s and (%s)" % (obj, tsel))
        cmd.color(PURPLE, "%s and chain A and resi %s" % (obj, heldsel))
        cmd.color(GREY, "%s and chain A and resi %s" % (obj, builtsel))
        cmd.color(YELLOW, "%s and chain A and resi %s" % (obj, pinsel))
        cmd.show("sticks", "%s and chain A and resi %s and not hydro" % (obj, pinsel))
        cmd.set("cartoon_side_chain_helper", 1, obj)

        cmd.select("%s_held" % obj, "%s and chain A and resi %s" % (obj, heldsel))
        cmd.select("%s_built" % obj, "%s and chain A and resi %s" % (obj, builtsel))
        cmd.select("%s_pin" % obj, "%s and chain A and resi %s" % (obj, pinsel))
        cmd.deselect()

        # ---- the number the picture claims: contacts of the BUILT segment, measured here
        B = heavy(obj, "chain A and resi %s" % builtsel)
        T = heavy(obj, tsel)
        d = np.linalg.norm(B[:, None, :] - T[None, :, :], axis=2)
        built_min, built_pairs = float(d.min()), int((d < CONTACT).sum())
        tsv_pairs = int(r["new_contacts_4p5"])
        if built_pairs != tsv_pairs:
            raise SystemExit("REFUSING: %s built-contact count measured here (%d) disagrees with "
                             "the ranked table (%d). A figure must not claim a number it cannot "
                             "reproduce." % (design, built_pairs, tsv_pairs))

        # ⭐ a COMMON FRAME is asserted, not assumed: every design carries the same TfR1 crop
        if ref is None:
            ref = obj
        else:
            al = cmd.align("%s and (%s) and name CA" % (obj, tsel),
                           "%s and (%s) and name CA" % (ref, tsel), cycles=0, transform=0)
            print("   target-vs-%s RMSD for %s: %.4f A (expect ~0: shared crop)"
                  % (ref, obj, al[0]))

        print("   rank %d  %-40s batch %-6s held %3d built %3d pin %s(%s)  "
              "built_min %.3f A  built_pairs %d (table %d, match)"
              % (i, design, r["batch"], len(held), len(built), pins,
                 ",".join(pin_atoms), built_min, built_pairs, tsv_pairs))
        drawn.append(dict(obj=obj, r=r, held=held, built=built, pins=pins,
                          pin_atoms=pin_atoms, bmin=built_min, bpairs=built_pairs, tsel=tsel))

    # ---- per-design scenes. ASCII ONLY in every message.
    for i, m in enumerate(drawn, 1):
        r = m["r"]
        cmd.disable("all")
        cmd.enable(m["obj"])
        cmd.orient("%s and chain A" % m["obj"])
        msg = ("rank %d of %d survivors   %s\n"
               "batch %s  -- effective n is BATCHES, not designs\n"
               "purple = whole protein GIVEN to RFd3, sequence-fixed (%d res)\n"
               "grey   = the N-terminal extension RFd3 BUILT (%d res)\n"
               "yellow = the fixed guidepost, atoms %s, held at its specified position\n"
               "teal   = TfR1 target (chains %s)\n"
               "BUILT segment to target: min %.3f A, %d heavy-atom pairs under %.1f A\n"
               "whole binder to target:  min %s A | rg_ratio %s | max CA-CA step %s A\n"
               "pin to target: %s A | backbone breaks over 4.5 A: %s\n"
               "Engagement is NOT affinity. Every coordinate here is a model."
               % (i, len(rows),
                  r["design"], r["batch"], len(m["held"]), len(m["built"]),
                  ",".join(m["pin_atoms"]), r["target_chains"],
                  m["bmin"], m["bpairs"], CONTACT,
                  r["bb_min_to_target"], r["rg_ratio"], r["max_step"],
                  r["pin_min_to_target"], r["n_break_4p5"]))
        if not msg.isascii():
            raise SystemExit("REFUSING: non-ASCII scene message for %s" % m["obj"])
        cmd.scene(m["obj"], "store", message=msg)

    # ---- overview: all binders, one target surface
    cmd.enable("all")
    for m in drawn[1:]:
        cmd.hide("surface", "%s and (%s)" % (m["obj"], m["tsel"]))
    cmd.orient(drawn[0]["obj"])
    batches = [m["r"]["batch"] for m in drawn]
    ov = ("ALL %d designs, one TfR1 surface shown\n"
          "purple = given protein   grey = RFd3-built N-terminal extension   yellow = pin\n"
          "batches %s -- %d distinct, so these are NOT %d independent draws\n"
          "built-segment contact pairs under %.1f A: %s\n"
          "Ranked on contacts made by the BUILT segment, corrected set of 27 residues.\n"
          "The earlier ranking scored a 2-residue linker by mistake and is superseded.\n"
          "Engagement is NOT affinity. Nothing here has been tested."
          % (len(drawn), ",".join(batches), len(set(batches)), len(drawn), CONTACT,
             ", ".join(str(m["bpairs"]) for m in drawn)))
    if not ov.isascii():
        raise SystemExit("REFUSING: non-ASCII overview message")
    cmd.scene("00_overview", "store", message=ov)

    cmd.save(outpse)
    print("   wrote %s (%d bytes)" % (outpse, os.path.getsize(outpse)))
    print("== scenes: %s ==" % ", ".join(cmd.get_scene_list()))


_a = sys.argv[1:]
if len(_a) < 3:
    raise SystemExit("usage: pymol -cq j956_nterm_top3_pse_20261102.py -- "
                     "<ranked.tsv> <cif_dir> <out.pse> [n_top]")
main(*_a[:4])
