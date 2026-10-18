#!/usr/bin/env python
"""j494_verify -- open a built .pse in a FRESH interpreter and assert it survives the load.

⛔⛔ WHY THIS IS A SEPARATE PROCESS: `cmd.reinitialize()` + `cmd.load()` in the SAME interpreter
   that built the scene SEGFAULTS (measured 2026-10-18, rc=139 at 16 objects / 13 scenes) and
   takes the whole buffered stdout with it, so a perfectly good .pse reads as "produced nothing".
   A fresh process is also the stronger check: it is exactly what opening the file elsewhere does.

⛔ NO `__main__` guard and NO `os._exit()`. Prints are flushed
   explicitly, because a crash in PyMOL loses anything still sitting in Python's stdout buffer.

usage: pymol -cq j494_verify_pse_20261018.py -- <pse> <n_objects> <comma,separated,scenes>
"""
import os
import sys

from pymol import cmd

PSE = sys.argv[1]
WANT_OBJ = int(sys.argv[2])
WANT_SCENES = sys.argv[3].split(",")

print(f"[verify] loading {PSE} ({os.path.getsize(PSE)/1e6:.2f} MB)")
sys.stdout.flush()

cmd.load(PSE)

objs = cmd.get_names("objects")
scenes = cmd.get_scene_list()
print(f"[verify] objects {len(objs)}: {objs}")
print(f"[verify] scenes  {len(scenes)}: {scenes}")
sys.stdout.flush()

assert len(objs) == WANT_OBJ, f"objects {len(objs)} != expected {WANT_OBJ}"
assert scenes == WANT_SCENES, f"scenes {scenes} != expected {WANT_SCENES}"

# every scene must RECALL without dying -- a corrupt message string kills it here, not at load
for s in scenes:
    cmd.scene(s, "recall")
print(f"[verify] all {len(scenes)} scenes recalled cleanly")

# the binder/target/core/pin objects must all carry atoms -- an empty object means a selection
# silently matched nothing and the scene would look like a real finding
empty = [o for o in objs if cmd.count_atoms(o) == 0]
assert not empty, f"objects present but EMPTY: {empty}"
print(f"[verify] every object carries atoms (min {min(cmd.count_atoms(o) for o in objs)})")
print("[verify] OK")
sys.stdout.flush()
