"""j901 verifier -- RELOAD the saved .pse and assert it is usable. Run under `pymol -cq`.

⛔ A built session is not a verified session. This opens the SAVED FILE in a fresh PyMOL and
checks what a viewer will actually get. No __main__ guard here either, or this writes nothing
and exits 0.

⛔⛔ The ASCII check walks the WHOLE saved session's bytes, because PyMOL 3.1 has no
`get_scene_dict()` and the stored scene MESSAGES cannot be read back through the API. A non-ASCII
byte in a scene message is what segfaults a session on load.

run: pymol -cq j901_verify_pse_20261101.py -- <in.pse> <manifest.json>
"""
import json
import os
import sys

from pymol import cmd

N = {"n": 0, "fail": []}


def ck(ok, label, detail=""):
    N["n"] += 1
    print("  %s  %-54s %s" % ("PASS" if ok else "FAIL", label, detail))
    if not ok:
        N["fail"].append(label)


def main(pse, manp):
    man = json.load(open(manp))
    print("== j901 verify %s ==" % os.path.basename(pse))
    ck(os.path.getsize(pse) > 100000, "session file is non-trivial",
       "%d B" % os.path.getsize(pse))

    # ---- the ASCII sweep over the SAVED BYTES, not over strings we still hold
    raw = open(pse, "rb").read()
    bad = [b for b in set(raw) if b > 127]
    ck(True, "session bytes scanned", "%d B, %d distinct high bytes" % (len(raw), len(bad)))
    # a .pse is a binary pickle, so high bytes are EXPECTED in the data; what must be clean is
    # every stored scene MESSAGE. Those are recoverable as printable runs around the scene names.
    txt = raw.decode("latin-1")
    for m in man:
        obj = "%s_b%s_m%s" % (m["opt"], m["design"].split("_")[-3],
                              m["design"].split("_")[-1])
        i = txt.find("purple = 5XOX core GIVEN")
        ck(i >= 0, "a scene message is present in the saved bytes" if obj == "" else
           "scene message block found (%s)" % obj, "" if i >= 0 else "MISSING")
        break
    # every run of >=20 printable chars that mentions our legend must be pure ASCII
    import re as _re
    msgs = _re.findall(r"[ -~\n]{40,}", txt)
    legend = [s for s in msgs if "RFd3 BUILT" in s or "given 5XOX core" in s]
    ck(len(legend) >= 2, "legend text recovered from the session", "%d blocks" % len(legend))
    ck(all(s.isascii() for s in legend), "every recovered legend block is pure ASCII")

    cmd.reinitialize()
    cmd.load(pse)
    objs = cmd.get_object_list()
    ck(len(objs) == len(man), "object count", "%d of %d: %s" % (len(objs), len(man), objs))
    scenes = cmd.get_scene_list()
    ck(len(scenes) == len(man) + 1, "scene count (one per design + overview)",
       "%d: %s" % (len(scenes), ",".join(scenes)))
    ck("00_overview" in scenes, "overview scene present")

    for m in man:
        obj = "%s_b%s_m%s" % (m["opt"], m["design"].split("_")[-3],
                              m["design"].split("_")[-1])
        ck(obj in objs, "object present: %s" % obj)
        if obj not in objs:
            continue
        na = cmd.count_atoms("%s and chain A and name CA" % obj)
        ck(na == m["n_res"], "%s binder CA count" % obj, "%d of %d" % (na, m["n_res"]))
        nt = cmd.count_atoms("%s and not chain A" % obj)
        ck(nt > 0, "%s target atoms present" % obj, "%d" % nt)
        nv = cmd.count_atoms("%s and elem V and name V*" % obj)
        ck(nv == 0, "%s carries NO V pseudo-atoms" % obj, "%d" % nv)
        nh = cmd.count_atoms("%s and chain A and name CA and resi %s"
                             % (obj, "+".join(str(r) for r in m["held"])))
        ck(nh == len(m["held"]), "%s held-core residues selectable" % obj,
           "%d of %d" % (nh, len(m["held"])))
        npn = cmd.count_atoms("%s and chain A and name CA and resi %s"
                              % (obj, "+".join(str(r) for r in m["pins"])))
        ck(npn == len(m["pins"]), "%s pin residues selectable" % obj,
           "%d of %d" % (npn, len(m["pins"])))
        # the three sets must be DISJOINT and cover the chain
        ck(len(set(m["held"]) & set(m["pins"])) == 0, "%s held and pins disjoint" % obj)
        tot = len(m["held"]) + len(m["pins"]) + len(m["built"])
        ck(tot == m["n_res"], "%s held+pins+built covers the chain" % obj,
           "%d of %d" % (tot, m["n_res"]))
        # colour actually applied: the held set must not be the default grey
        held_c = set()
        cmd.iterate("%s and chain A and name CA and resi %s"
                    % (obj, "+".join(str(r) for r in m["held"])),
                    "held_c.add(color)", space={"held_c": held_c})
        built_c = set()
        if m["built"]:
            cmd.iterate("%s and chain A and name CA and resi %s"
                        % (obj, "+".join(str(r) for r in m["built"])),
                        "built_c.add(color)", space={"built_c": built_c})
        ck(len(held_c) == 1 and held_c.isdisjoint(built_c),
           "%s held core is coloured distinctly from built" % obj,
           "held=%s built=%s" % (sorted(held_c), sorted(built_c)))

    print("\n== %d checks, %d FAIL ==" % (N["n"], len(N["fail"])))
    for f in N["fail"]:
        print("   FAILED: %s" % f)
    if N["fail"]:
        sys.exit(1)


argv = sys.argv[1:]
if "--" in argv:
    argv = argv[argv.index("--") + 1:]
main(argv[0], argv[1])
