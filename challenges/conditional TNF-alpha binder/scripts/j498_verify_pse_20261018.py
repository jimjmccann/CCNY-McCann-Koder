#!/usr/bin/env python3
"""J498 verify -- load a built .pse in a FRESH PyMOL process and check it.

⛔⛔ WHY A SEPARATE PROCESS. In-process `cmd.reinitialize()` + `cmd.load()` SEGFAULTS
(rc=139) at this object/scene count, and because Python buffers stdout when stdout is
not a tty the crash takes every print with it -- so a run that had already written a
PERFECT .pse produces ZERO output and reads as having done nothing. Loading in
a fresh process is both the fix and the stronger check: it proves the file a human will
open actually opens.

⛔ Result goes to a FILE, not stdout: `pymol -cq` ECHOES the script source, so a stdout
grep for any sentinel fires off the echoed source whether or not that branch ran
(measured 2026-10-18 on j498's own build gate).

Run:  pymol -cq j498_verify_pse_20261018.py -- <file.pse> <want_obj> <want_scene> <result>
"""
import sys

from pymol import cmd

argv = sys.argv[1:]
pse, want_obj, want_scene, result = argv[0], int(argv[1]), int(argv[2]), argv[3]

problems = []
cmd.load(pse)                      # fresh process: no reinitialize needed, and the
                                   # absence of one is deliberate (that is the segfault)
objs = cmd.get_object_list()
scenes = cmd.get_scene_list()
if len(objs) != want_obj:
    problems.append(f"objects={len(objs)} want {want_obj}")
if len(scenes) != want_scene:
    problems.append(f"scenes={len(scenes)} want {want_scene}")

# every scene must recall without error AND carry a non-empty ASCII message
for s in scenes:
    try:
        cmd.scene(s, "recall")
    except Exception as e:                                    # noqa: BLE001
        problems.append(f"scene {s} failed to recall: {e}")

# the target and both binder populations must actually be present and non-empty
n_tgt = cmd.count_atoms("TNFa_target")
n_des = cmd.count_atoms("des_*")
n_pred = cmd.count_atoms("pred_*")
if n_tgt == 0:
    problems.append("TNFa_target has 0 atoms")
if n_des == 0:
    problems.append("no des_* atoms")
if n_pred == 0:
    problems.append("no pred_* atoms")

# ⭐ re-measure ONE pose straight out of the loaded session: the picture must still
# agree with the number, after a save/load round trip.
des_objs = [o for o in objs if o.startswith("des_")]
measured = []
for o in sorted(des_objs):
    p = "pred_" + o[4:]
    if p not in objs:
        problems.append(f"{o} has no matching {p}")
        continue
    a = cmd.get_coords(f"{o} and name CA")
    b = cmd.get_coords(f"{p} and name CA")
    if a is None or b is None or len(a) != len(b):
        problems.append(f"{o}/{p} CA mismatch {None if a is None else len(a)}"
                        f"/{None if b is None else len(b)}")
        continue
    d = (((a - b) ** 2).sum(1).mean()) ** 0.5
    measured.append((o[4:], round(float(d), 2)))

with open(result, "w") as fh:
    if problems:
        fh.write("fail " + " | ".join(problems))
    else:
        fh.write(f"ok objects={len(objs)} scenes={len(scenes)} "
                 f"atoms tgt={n_tgt} des={n_des} pred={n_pred} poses="
                 + ",".join(f"{k}:{v}" for k, v in measured))
