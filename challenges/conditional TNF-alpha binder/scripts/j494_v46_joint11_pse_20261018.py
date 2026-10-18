#!/usr/bin/env python
"""j494 -- ONE .pse of the 11 JOINT-passing C1 `v46` backbones.

⭐ FULL, BALANCED COHORT (384/384, 128/arm) -- NOT j489's mid-run partial look. The JOINT gate is
   `continuous` (4.2 A) AND `min_cut >= 19`, reproduced here from `J490`'s own two TSVs: arms come
   out 5 (b1) + 6 (b2) + 0 (b3) = 11, and continuity 15/7/1 = 23, which is `J490` section 3 exactly.

⛔⛔ THE FIVE BUILD TRAPS, carried VERBATIM in substance from this project's earlier .pse builders:
  1. strip to CA at LOAD for bulk objects. (Here n=11, so every binder keeps full backbone and
     the .pse stays small -- the trap is recorded because the next script may not be n=11.)
  2. `scene(message=)` MUST be pure ASCII -- a non-ASCII mark corrupts the heap and the saved
     .pse SEGFAULTS on load: a scene message must be ASCII. Everything goes through
     ascii_msg().
  3. enable member OBJECTS by name, never a group.
  4. `cartoon_trace_atoms 0` + explicit `dss`, or CA-only objects render as garbage.
  5. ⛔ NO `__main__` guard: `pymol -cq` exits 0 writing NOTHING if the
     body is guarded. And ⛔ NO `os._exit()` -- it SKIPS THE STDOUT FLUSH, so the reload check
     prints nothing and reads as success (J475).

⛔⛔ INPUT AND OUTPUT DO NOT SHARE A FRAME -- RFd3 re-centres. Each arm's bundle input is SUPERPOSED
   onto the output target and the RMSD is ASSERTED before its core or pins are drawn; a core drawn
   in the wrong frame looks exactly like a real finding.

⛔ THE LAYOUT IS READ FROM THE SCORED TSV, never assumed, and all 11
   rows must agree before anything is drawn.

⛔ The pins drawn here are the INPUT (spec) pins superposed into the output frame. That is honest
   because pins move <=0.07 A input->output over 384 designs, but it is the ASK, not
   an independent measurement of the output.

usage: pymol -cq j494_v46_joint11_pse_20261018.py -- <extract-dir> <out-dir>
"""
import csv
import os
import re
import statistics
import sys

from pymol import cmd

# PLACEHOLDER PATH -- not a working default; set it for your own tree.
REPO = "/path/to/your/working/tree"
BUNDLE_DIR = f"{REPO}/bundles/rfd3_tnf_c1core_v46"
CONT_TSV = f"{REPO}/data/j490_continuity_c1core_v46_FULL_20261018.tsv"
PATH_TSV = f"{REPO}/data/j490_pathfold_c1core_v46_FULL_20261018.tsv"

EXTRACT = sys.argv[1]
OUT_DIR = sys.argv[2]

WARN = 4.2          # the continuity bar
MIN_CUT_BAR = 19    # the topology bar (LCB1's median)

ASCII_MAP = [("≥", ">="), ("≤", "<="), ("Å", "A"), ("→", "->"),
             ("⭐", "[STAR]"), ("⛔", "[STOP]"), ("⚠", "[WARN]"), ("⇒", "=>")]


def ascii_msg(s):
    for a, b in ASCII_MAP:
        s = s.replace(a, b)
    return s.encode("ascii", "replace").decode("ascii")


def load_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"))


def key(d):
    return d[:-7] if d.endswith(".cif.gz") else d


# ---- per-arm core + input, read from the bundle's own ARMS.tsv rather than hardcoded
ARM_OF = {}
for row in csv.DictReader(open(f"{BUNDLE_DIR}/ARMS.tsv"), delimiter="\t"):
    ARM_OF[row["arm"].lower()] = (row["input_pdb"], int(row["n_core"]), row["core_select"])

# ---- the JOINT gate, recomputed from the FULL TSVs
cont = load_tsv(CONT_TSV)
path = load_tsv(PATH_TSV)
assert len(cont) == 384 and len(path) == 384, f"not the full cohort: {len(cont)}/{len(path)}"
MC = {r["design"]: float(r["min_cut"]) for r in path}
RG = {r["design"]: float(r["rg_ratio"]) for r in path}

winners = []
n_cont_all = 0
for r in cont:
    k = key(r["design"])
    is_cont = r["continuous"] == "1"
    n_cont_all += is_cont
    if is_cont and MC[k] >= MIN_CUT_BAR:
        r["min_cut_f"] = MC[k]
        r["rg_ratio_f"] = RG[k]
        r["max_step_f"] = float(r["max_step"])
        r["arm"] = re.search(r"_c1core_v46_(b\d)_", r["design"]).group(1)
        winners.append(r)

# ⛔⛔ QUOTE THIS COHORT'S OWN DISTRIBUTION, NEVER ANOTHER CAMPAIGN'S. The first build of
#    this script put `J411`'s "our cohort median 1.26" in every per-design message. That figure is
#    the v38-era cohort's; `v46`'s OWN rg_ratio median is 0.969 (n=384), so the quoted number was
#    from a different campaign and would have read as these 11 beating their cohort by a mile.
RG_ALL = sorted(float(r["rg_ratio"]) for r in path)
RG_MED = statistics.median(RG_ALL)
RG_Q1, RG_Q3 = statistics.quantiles(RG_ALL, n=4)[0], statistics.quantiles(RG_ALL, n=4)[2]
MC_MED = statistics.median(float(r["min_cut"]) for r in path)
print(f"[cohort] v46 own rg_ratio median {RG_MED:.3f} (IQR {RG_Q1:.3f}-{RG_Q3:.3f}), "
      f"min_cut median {MC_MED:.1f}")

assert n_cont_all == 23, f"continuity does not reproduce J490: {n_cont_all} vs 23"
assert len(winners) == 11, f"JOINT does not reproduce J490: {len(winners)} vs 11"
by_arm = {a: sum(1 for w in winners if w["arm"] == a) for a in ("b1", "b2", "b3")}
assert by_arm == {"b1": 5, "b2": 6, "b3": 0}, f"arm split wrong: {by_arm}"

# ---- the layout, READ not assumed, and all 11 must agree
binders = {w["binder_chain"] for w in winners}
targets = {w["chains_target"] for w in winners}
assert len(binders) == 1 and len(targets) == 1, f"layout disagrees across winners: {binders} {targets}"
BINDER = binders.pop()
TARGET_OUT = targets.pop().split(",")
TARGET_IN = ["A", "B", "C"]          # the INPUT file's TNFa chains
print(f"[layout] binder chain {BINDER}, target chains {TARGET_OUT} (read from the TSV)")

# best first, by topology then continuity
winners.sort(key=lambda w: (-w["min_cut_f"], w["max_step_f"]))

ARM_COLOR = {"b1": "palegreen", "b2": "skyblue"}
OUT = f"{OUT_DIR}/J494_C1CORE_V46_JOINT11_FULL_20261018.pse"

cmd.set("cartoon_trace_atoms", 0)          # trap 4
cmd.set("assembly", "")
cmd.set("cartoon_transparency", 0.0)

# ---- load the 11 binders
target_made = False
obj_of = {}
for w in winners:
    base = w["design"]
    src = os.path.join(EXTRACT, base)
    assert os.path.exists(src), f"missing extracted design: {src}"
    obj = "d_%s_%s" % (w["arm"], re.search(r"_s(\d+_\d+)_model_0", base).group(1))
    cmd.load(src, obj)
    if not target_made:
        cmd.create("target", f"{obj} and chain {'+'.join(TARGET_OUT)}")
        n_t = cmd.count_atoms("target")
        assert n_t > 2000, f"target selection empty/short ({n_t}) -- OUTPUT layout is wrong"
        print(f"[target] {n_t} atoms from {'+'.join(TARGET_OUT)}")
        target_made = True
    cmd.remove(f"{obj} and not chain {BINDER}")
    n_b = cmd.count_atoms(f"{obj} and name CA")
    assert n_b == int(w["n_res"]), f"{obj}: {n_b} CA vs n_res {w['n_res']}"
    obj_of[base] = obj

# ---- each arm's core + pins, SUPERPOSED into the OUTPUT frame, RMSD ASSERTED
pins_made = False
core_objs = []
for arm in ("b1", "b2"):
    input_pdb, n_core, core_sel = ARM_OF[arm]
    cmd.load(f"{BUNDLE_DIR}/{input_pdb}", "bundle")
    rms = cmd.align(f"bundle and chain {'+'.join(TARGET_IN)} and name CA",
                    "target and name CA")[0]
    print(f"[align] {arm} bundle target -> output target RMSD = {rms:.3f} A")
    assert rms < 1.0, (f"FATAL: {arm} superposition RMSD {rms:.3f} A too high -- the core would be "
                       "drawn in the wrong frame and would look exactly like a real finding")
    co = f"core_{arm}"
    cmd.create(co, "bundle and chain F")
    n_ca = cmd.count_atoms(co)
    assert n_ca == 2 * n_core, f"{co}: expected {2*n_core} atoms (CA+CB x {n_core}), got {n_ca}"
    core_objs.append(co)
    print(f"[core] {co}: N={n_core} {core_sel}, {n_ca} atoms")
    if not pins_made:
        cmd.create("pin_his", "bundle and chain D")
        cmd.create("pin_arg", "bundle and chain E")
        n_h, n_a = cmd.count_atoms("pin_his"), cmd.count_atoms("pin_arg")
        assert n_h == 10, f"pin_his expected 10 atoms (1 HIS), got {n_h}"
        assert n_a == 33, f"pin_arg expected 33 atoms (3 x 11-atom ARG), got {n_a}"
        print(f"[pins] HIS {n_h} atoms + ARG {n_a} atoms = 4 pins")
        pins_made = True
    cmd.delete("bundle")

# ---- representation
cmd.hide("everything")
cmd.show("cartoon", "target")
cmd.color("grey80", "target")
cmd.dss("target")                                             # trap 4
cmd.set("cartoon_transparency", 0.55, "target")

for co in core_objs:
    cmd.show("spheres", co)
    cmd.color("orange", co)
    cmd.set("sphere_scale", 0.35, co)

cmd.show("sticks", "pin_his")
cmd.color("magenta", "pin_his")
cmd.show("sticks", "pin_arg")
cmd.color("red", "pin_arg")

design_objs = []
for w in winners:
    o = obj_of[w["design"]]
    design_objs.append(o)
    cmd.show("cartoon", o)
    cmd.dss(o)
    cmd.color(ARM_COLOR[w["arm"]], o)

SCENE_KW = dict(view=0, color=1, rep=1, active=1)


def only(objs):                                               # trap 3
    cmd.disable("all")
    for o in ["target"] + list(objs):
        cmd.enable(o)


scene_names = []


def scene(name, objs, zoom_sel, zoom_buf, msg):
    only(objs)
    cmd.zoom(zoom_sel, zoom_buf)
    cmd.scene(name, "store", **SCENE_KW, message=ascii_msg(msg))
    scene_names.append(name)


site = core_objs + ["pin_his", "pin_arg"]

scene("F1", design_objs + site, "visible", 3,
      f"1/{len(winners)+2} ALL {len(winners)} JOINT-PASSING C1 v46 BACKBONES, full balanced cohort "
      f"(384/384, 128/arm). JOINT = continuous at {WARN} A AND min_cut >= {MIN_CUT_BAR}: "
      f"{len(winners)}/384 = 2.9%. Green = arm b1 (core N=8 compact, 5 designs), "
      f"blue = arm b2 (core N=8 near_hyd, 6 designs). b3 (N=12) contributed ZERO. "
      f"Orange spheres = the requested core CA+CB. Magenta = the HIS pin, red = the 3 ARG pins. "
      f"Target is transparent grey. "
      f"[WARN] This is the TOP 2.9% of the run, NOT a typical design, and NOTHING here is folded.")

scene("F2", site, "core_b1 or pin_arg", 8,
      f"2/{len(winners)+2} THE ASK ALONE -- no designed chain shown, so the empty space is where "
      f"the binder must be built. 4 pins: HIS (magenta, on the 3-fold axis) + 3 ARG (red, to "
      f"GLU107 of each protomer). Two orange cores overlaid: b1 N=8 compact and b2 N=8 near_hyd. "
      f"[STOP] These are the INPUT pins superposed into the output frame (RMSD asserted <1 A); "
      f"pins move <=0.07 A input->output, so this is the ASK, not a measurement of the output.")

for i, w in enumerate(winners, start=3):
    o = obj_of[w["design"]]
    core_arm = f"core_{w['arm']}"
    scene("D%02d" % (i - 2), [o, core_arm, "pin_his", "pin_arg"], f"{o} or {core_arm}", 4,
          f"{i}/{len(winners)+2} arm {w['arm'].upper()} rank {i-2} of {len(winners)}: "
          f"{w['design'][:52]} | min_cut {w['min_cut_f']:.0f} (bar {MIN_CUT_BAR}, LCB1 median 19) "
          f"| max_step {w['max_step_f']:.3f} A (bar {WARN}) | rg_ratio {w['rg_ratio_f']:.3f} "
          f"(v46's OWN cohort median {RG_MED:.3f}, IQR {RG_Q1:.3f}-{RG_Q3:.3f}) | "
          f"{w['n_res']} residues | "
          f"{w['n_pins']} pins, first {w['first_pin']} last {w['last_pin']}. "
          f"min_cut = can the chain be cut in two with few long-range contacts crossing; "
          f"HIGHER is more one-fold. [STOP] Backbone only -- no sequence, no fold prediction.")

cmd.scene("F1", "recall")
os.makedirs(OUT_DIR, exist_ok=True)
cmd.save(OUT)

# ---- ⛔⛔ SIXTH TRAP, MEASURED HERE 2026-10-18: do NOT reload-verify IN THIS PROCESS.
#      j489's `cmd.reinitialize()` + `cmd.load(OUT)` in the SAME interpreter SEGFAULTS at this
#      object/scene count (rc=139, core dumped, and because Python buffers stdout when it is not a
#      tty the segfault ALSO LOSES EVERY PRINT -- the run looked like it produced nothing at all).
#      The .pse written just before the crash was PERFECT: 16 objects, 13 scenes, rc=0 on load.
#      ⇒ verify in a SEPARATE process (j494_verify), which is the stronger check anyway because a
#      fresh interpreter is exactly what opening the file on another machine does.
#      ⛔ never os._exit() either: it skips the stdout flush (J475).
n_obj = len(cmd.get_names("objects"))
print(f"[built] {OUT}")
print(f"[built] objects {n_obj} | scenes {scene_names}")
print(f"[built] n_joint={len(winners)} b1={by_arm['b1']} b2={by_arm['b2']} b3={by_arm['b3']} "
      f"continuous={n_cont_all}/384")
print(f"[built] size {os.path.getsize(OUT)/1e6:.2f} MB")
print(f"[built] EXPECT_OBJECTS={n_obj}")
print(f"[built] EXPECT_SCENES={','.join(scene_names)}")
sys.stdout.flush()
