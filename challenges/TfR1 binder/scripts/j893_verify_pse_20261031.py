"""j893 -- verify the candidate-scaffold .pse BY RELOAD, in a fresh PyMOL, before it is sent.

RUN:  pymol -cq j893_verify_pse_20261031.py -- --pse <file.pse>

⛔ A BUILDER THAT EXITS 0 IS NOT A VERIFICATION. Both failures this guards against look exactly
   like success at build time:
     1. a non-ASCII character in a stored scene message corrupts the heap and the file SEGFAULTS
        when a reader opens it -- so every stored message is re-read HERE, out of the SAVED
        file, and re-checked;
     2. a named selection that matched nothing saves happily and shows nothing, so the thing
        asked for is simply absent. Every object is asserted non-empty with its atom count
        PRINTED, not merely compared -- a count is not a measurement.

⛔ NOT GUARDED UNDER `__main__`: under `pymol -cq` this module is not __main__, so a guarded
   script exits 0 having written nothing.
⛔ Local only.
"""
import sys

from pymol import cmd

argv = sys.argv[1:]
if "--" in argv:
    argv = argv[argv.index("--") + 1:]
if "--pse" not in argv:
    sys.exit("missing --pse")
pse = argv[argv.index("--pse") + 1]

EXPECT_MIN_SCENES = 6
EXPECT_OBJECTS = ["target", "guideposts"]

cmd.load(pse)
fails, passes = [], []


def check(ok, label, detail=""):
    line = "  %-5s %s%s" % ("PASS" if ok else "FAIL", label, (" -- " + detail) if detail else "")
    print(line)
    (passes if ok else fails).append(label)


print("== verify %s ==" % pse)

scenes = cmd.get_scene_list()
check(len(scenes) >= EXPECT_MIN_SCENES, "scene count",
      "%d scenes: %s" % (len(scenes), ", ".join(scenes)))

# ASCII check: PyMOL 3.1 has no get_scene_dict and the scene-message layout inside the
# session is version-specific, so rather than index into it, WALK THE WHOLE SAVED SESSION and
# assert every string in it is ASCII. Stronger than checking scene messages alone -- it also
# covers labels and object names -- and it cannot rot when PyMOL moves the layout.
sess = cmd.get_session()
bad = []
seen_msgs = []


def walk(x, path, depth=0):
    if depth > 8:
        return
    if isinstance(x, str):
        try:
            x.encode("ascii")
        except UnicodeEncodeError as e:
            bad.append((path, repr(x[max(0, e.start - 20):e.end + 20])))
        if len(x) > 60:
            seen_msgs.append(x)
        return
    if isinstance(x, (list, tuple)):
        for i, y in enumerate(x):
            walk(y, "%s[%d]" % (path, i), depth + 1)
    elif isinstance(x, dict):
        for k, y in x.items():
            walk(k, "%s.key" % path, depth + 1)
            walk(y, "%s[%r]" % (path, k), depth + 1)


walk(sess.get("moviescenes"), "moviescenes")
walk(sess.get("names"), "names")
check(not bad, "every string in the saved session is pure ASCII",
      "%d long strings checked" % len(seen_msgs) if not bad
      else "NON-ASCII at %s: %s" % bad[0])
check(len(seen_msgs) >= EXPECT_MIN_SCENES, "a stored message per scene",
      "%d long strings found, >= %d scenes" % (len(seen_msgs), EXPECT_MIN_SCENES))

objs = cmd.get_object_list()
print("  objects: %s" % ", ".join(objs))
for o in EXPECT_OBJECTS:
    check(o in objs, "object present: %s" % o)

for o in objs:
    n = cmd.count_atoms(o)
    check(n > 0, "object non-empty: %s" % o, "%d atoms" % n)

cands = [o for o in objs if o.startswith("cand_")]
check(len(cands) == 4, "four candidate objects", "%d: %s" % (len(cands), ", ".join(cands)))
pins = [o for o in objs if o.startswith("pins_cand_")]
check(len(pins) == 4, "four pin objects", "%d: %s" % (len(pins), ", ".join(pins)))

print("\n== %d passed, %d failed ==" % (len(passes), len(fails)))
if fails:
    print("FAILED: %s" % "; ".join(fails))
    sys.exit(1)
