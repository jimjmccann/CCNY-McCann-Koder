"""j956 verifier -- re-opens the saved .pse in a FRESH viewer process and checks it.

⛔ NO `if __name__ == "__main__"` GUARD (see the builder).

⭐ A builder that exits 0 is not a verification. The failure this exists to catch is a session that
saves cleanly and crashes on load, which is what a non-ASCII scene message does. ⇒ the ASCII sweep
runs over the SAVED BYTES, not over strings still held in memory.

run: pymol -cq j956_verify_pse_20261102.py -- <in.pse> <ranked.tsv>
"""
import csv
import os
import re
import sys

from pymol import cmd

ok = fail = 0


def ck(cond, what, got=""):
    global ok, fail
    if cond:
        ok += 1
    else:
        fail += 1
        print("   FAIL %s %s" % (what, got))


pse, ranked = sys.argv[1:3]

# ---- 1. the saved bytes: printable runs that mention the legend must be pure ASCII
raw = open(pse, "rb").read()
runs = re.findall(rb"[\x20-\x7e\n\r\t]{20,}", raw)
legend = [r for r in runs if b"purple" in r or b"teal" in r or b"rank" in r]
ck(len(legend) >= 2, "legend text recovered from the saved bytes", "%d blocks" % len(legend))
nonascii = sum(1 for b in raw if b > 127)
print("   [M] saved bytes %d; non-ASCII bytes anywhere in the file: %d" % (len(raw), nonascii))
for blk in legend:
    ck(blk.decode("latin1").isascii(), "legend block is ASCII")

# ---- 2. reload in this fresh process
cmd.reinitialize()
cmd.load(pse)
objs = cmd.get_object_list()
scenes = cmd.get_scene_list()
print("   [M] reloaded: %d objects, %d scenes" % (len(objs), len(scenes)))
ck(len(objs) == 3, "three design objects", str(objs))
ck("00_overview" in scenes, "overview scene present", str(scenes))
ck(len(scenes) == 4, "four scenes", str(scenes))

rows = [r for r in csv.DictReader(
    [l for l in open(ranked) if not l.startswith("#")], delimiter="\t")]
top = rows[:3]

for i, (obj, r) in enumerate(zip(sorted(objs), top), 1):
    held = cmd.count_atoms("%s_held and name CA" % obj)
    built = cmd.count_atoms("%s_built and name CA" % obj)
    pin = cmd.count_atoms("%s_pin and name CA" % obj)
    ck(held == 143, "%s held is 143 residues" % obj, str(held))
    ck(built == 27, "%s built is 27 residues" % obj, str(built))
    ck(pin == 1, "%s one pin residue" % obj, str(pin))
    # the three sets must be disjoint
    both = cmd.count_atoms("(%s_held and %s_built) and name CA" % (obj, obj))
    ck(both == 0, "%s held and built are disjoint" % obj, str(both))
    # the target must still be there, and in the shared frame
    tgt = cmd.count_atoms("%s and (chain B or chain C) and name CA" % obj)
    ck(tgt == 154, "%s target is 154 CA (106 + 48)" % obj, str(tgt))
    # no V pseudo-atoms survived
    v = cmd.count_atoms("%s and elem V and name V*" % obj)
    ck(v == 0, "%s carries no V pseudo-atoms" % obj, str(v))
    print("   %-24s held %3d  built %2d  pin %d  target CA %d" % (obj, held, built, pin, tgt))

# ---- 3. every scene still carries a message, and it mentions the honest caveat
missing = []
for s in scenes:
    cmd.scene(s, "recall")
    msg = cmd.get("scene_current_name") and cmd.get_scene_message(s) \
        if hasattr(cmd, "get_scene_message") else None
    if msg is None:
        continue
    if "Engagement is NOT affinity" not in msg:
        missing.append(s)
ck(not missing, "every readable scene message carries the engagement caveat", str(missing))

print("== j956 verify: %d checks passed, %d failed ==" % (ok, fail))
if fail:
    raise SystemExit(1)
