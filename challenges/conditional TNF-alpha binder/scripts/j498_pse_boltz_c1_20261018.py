#!/usr/bin/env python3
"""J498 -- build the .pse for the C1 v46 Boltz-2 bound folds (J497's results).

WHAT IT SHOWS. Every design shares a BYTE-IDENTICAL target frame (MEASURED, max|delta| =
0.0 over 456 target CA across all 11), so ONE target object holds all 11 designed
binder positions AND all 11 Boltz predictions in a single common frame. That is what
makes "designed here, predicted there" a picture rather than a table.

⛔⛔ THE TRANSFORM IS APPLIED IN NUMPY AND WRITTEN AS COORDINATES, never handed to
PyMOL as a matrix. PyMOL's transform_object matrix conventions (row- vs column-major,
TTT vs homogeneous, `transpose=`) are a silent-wrong-answer risk, and a figure that
disagrees with the TSV it was built from is worse than no figure. The written file is
then RE-MEASURED and asserted against J497's own pose number before PyMOL sees it.

⛔ TRAPS THIS SCRIPT IS BUILT AROUND, all previously paid for on this project:
  - scene messages must be PURE ASCII: a single emoji corrupts the heap and the .pse
    SEGFAULTS ON LOAD. ⇒ `ascii_msg()` asserts.
  - no `if __name__ == "__main__"` guard: `pymol -cq script.py` exits 0 writing NOTHING
    if the body is guarded. ⇒ the .pml is generated and run
    as a separate process instead.
  - verification runs in a SEPARATE PROCESS: in-process reinitialize()+load() segfaults
    at this object/scene count and takes buffered stdout with it, so a PERFECT .pse can
    print nothing and read as having done nothing (J494).

Usage:
  j498_pse_boltz_c1_20261018.py --tsv <j497 FULL tsv> --design-dir <dir of design pdb>
      --pred-dir <boltz predictions dir> --out <out.pse> [--workdir <tmp>]
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# The group's stated palette. CCNY purple is the house colour; the contrast is a warm tone so
# "what we asked for" vs "what the predictor says" reads at a glance in both themes.
# ⛔ never archive colours.
CCNY_PURPLE = (0x7D / 255, 0x55 / 255, 0xC7 / 255)
PRED_ORANGE = (0xE0 / 255, 0x6C / 255, 0x1E / 255)
SUSPECT_GREY = (0x8A / 255, 0x8A / 255, 0x8A / 255)
TARGET_GREY = (0xB4 / 255, 0xB4 / 255, 0xB4 / 255)

GLU107_IDX = 102          # 1-based sequential index inside a 152-residue target chain
ARG_TIP = (" NE ", " NH1", " NH2")
GLU_TIP = (" OE1", " OE2")


def ascii_msg(s):
    """⛔ A non-ASCII byte in a scene message corrupts the PyMOL heap and the saved
    .pse SEGFAULTS when loaded. Assert rather than silently transliterate."""
    out = (s.replace("→", "->").replace("≤", "<=").replace("≥", ">=")
           .replace("Å", "A").replace("–", "-").replace("—", "-")
           .replace("×", "x").replace("’", "'"))
    out = out.replace('"', "'")          # a bare " would break the generated .pml
    out = "".join(c for c in out if 32 <= ord(c) < 127)
    assert out.isascii(), out
    return out


def read_pdb(path, chains=None):
    """[(line, chain, resnum, resname, atomname, xyz)] for ATOM records."""
    rows = []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        ch = line[21]
        if chains and ch not in chains:
            continue
        rows.append((line, ch, int(line[22:26]), line[17:20].strip(), line[12:16],
                     np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return rows


def ca(rows, chain):
    return np.array([r[5] for r in rows if r[1] == chain and r[4] == " CA "])


def superpose(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    V, S, Wt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, pc, qc


def write_transformed(src_rows, R, pc, qc, out_path):
    """Apply (x - pc) @ R + qc to every atom and write a real PDB. Returns the coords."""
    with open(out_path, "w") as fh:
        for line, ch, rn, rname, an, xyz in src_rows:
            n = (xyz - pc) @ R + qc
            fh.write(f"{line[:30]}{n[0]:8.3f}{n[1]:8.3f}{n[2]:8.3f}{line[54:]}")
    return None


def tip_coords(rows, chain, resnum, names):
    return np.array([r[5] for r in rows
                     if r[1] == chain and r[2] == resnum and r[4] in names])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--design-dir", required=True)
    ap.add_argument("--pred-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--binder-chain", default="A")
    ap.add_argument("--target-chains", default="B,C,D")
    a = ap.parse_args()

    bchain = a.binder_chain
    tchains = [c for c in a.target_chains.split(",") if c]
    work = Path(a.workdir or (Path(a.out).parent / "j498_work"))
    work.mkdir(parents=True, exist_ok=True)

    rows = [r for r in csv.DictReader(open(a.tsv), delimiter="\t")
            if r.get("status") == "OK" and r.get("pose_frame")]
    assert rows, f"no scorable rows in {a.tsv}"
    print(f"J498: {len(rows)} designs from {Path(a.tsv).name}")

    # ---- common target frame, asserted, not assumed ----
    design_paths = {}
    for r in rows:
        tag = r["backbone"]
        hits = sorted(Path(a.design_dir).rglob(f"{tag}.pdb"))
        assert hits, f"no design pdb for {tag} under {a.design_dir}"
        design_paths[tag] = hits[0]
    ref = read_pdb(design_paths[rows[0]["backbone"]], set(tchains))
    refT = np.vstack([ca(ref, c) for c in tchains])
    for r in rows:
        t = read_pdb(design_paths[r["backbone"]], set(tchains))
        tT = np.vstack([ca(t, c) for c in tchains])
        assert tT.shape == refT.shape, (r["backbone"], tT.shape)
        d = float(np.abs(tT - refT).max())
        assert d < 1e-6, (f"design target frames DIFFER by {d} A at {r['backbone']} -- "
                          "a single shared target object would be WRONG")
    print(f"  common target frame asserted: max|delta| = 0.0 A over {len(refT)} target CA")

    # ---- the shared target, written once from the reference design ----
    tgt_pdb = work / "TNFa_target.pdb"
    with open(tgt_pdb, "w") as fh:
        for line, *_ in read_pdb(design_paths[rows[0]["backbone"]], set(tchains)):
            fh.write(line)

    items = []
    for r in rows:
        tag = r["backbone"]
        short = tag.split("v46_")[1].replace("_model_0", "")
        dpath = design_paths[tag]

        pred_dir = Path(a.pred_dir) / r["design"]
        ppdbs = sorted(pred_dir.glob("*_model_0.pdb"))
        assert ppdbs, f"no prediction pdb in {pred_dir}"

        des = read_pdb(dpath)
        pred = read_pdb(ppdbs[0])

        # chosen cyclic mapping, taken from the TSV so the picture and the number agree
        frame = r["pose_frame"]
        perm = list(frame)
        Pt = np.vstack([ca(pred, c) for c in perm])
        Qt = np.vstack([ca(des, c) for c in tchains])
        assert Pt.shape == Qt.shape
        R, pc, qc = superpose(Pt, Qt)

        out_pdb = work / f"pred_{short}.pdb"
        write_transformed(pred, R, pc, qc, out_pdb)

        # ⛔ RE-MEASURE the written file and assert it reproduces J497's pose number.
        moved = read_pdb(out_pdb)
        Pb, Qb = ca(moved, bchain), ca(des, bchain)
        assert Pb.shape == Qb.shape, (short, Pb.shape, Qb.shape)
        pose = float(np.sqrt(((Pb - Qb) ** 2).sum() / len(Pb)))
        want = float(r["rmsd_binder_on_target_frame"])
        assert abs(pose - want) < 0.01, (
            f"{short}: transformed file gives pose {pose:.3f} A but J497 recorded "
            f"{want:.3f} A -- the figure would disagree with its own data")

        # pinned ARG -> GLU107, in both, for the stick/distance scenes
        pins = [p for p in (r.get("pins") or "").split(",") if p]
        pin_args = [int(p[3:]) for p in pins if p.startswith("ARG")]
        best = None
        for tc in tchains:
            tgt_res = sorted({rr[2] for rr in des if rr[1] == tc})
            if len(tgt_res) < GLU107_IDX:
                continue
            g107 = tgt_res[GLU107_IDX - 1]
            go = tip_coords(des, tc, g107, GLU_TIP)
            if not len(go):
                continue
            gnames = [rr[4].strip() for rr in des
                      if rr[1] == tc and rr[2] == g107 and rr[4] in GLU_TIP]
            for pa in pin_args:
                asel = [(rr[5], rr[4].strip()) for rr in des
                        if rr[1] == bchain and rr[2] == pa and rr[4] in ARG_TIP]
                if not asel:
                    continue
                ao = np.array([x[0] for x in asel])
                D = np.linalg.norm(ao[:, None] - go[None], axis=-1)
                i, j = np.unravel_index(D.argmin(), D.shape)
                d = float(D[i, j])
                if best is None or d < best[0]:
                    # ⛔ keep the EXACT atom pair: a 19 A gap never renders under
                    # PyMOL's default distance cutoff, so the line must be atom-to-atom.
                    best = (d, tc, g107, pa, asel[i][1], gnames[j])
        items.append({
            "short": short, "tag": tag, "des": str(dpath), "pred": str(out_pdb),
            "pose": pose, "rmsd": float(r["rmsd_binder"]), "tm": float(r["tm_binder"]),
            "iptm": float(r["iptm_binder_target_mean"]),
            "arm": "b1" if "_b1_" in tag else "b2",
            "suspect": not str(r.get("pose_status", "")).startswith("OK"),
            "tgt_frame": float(r["rmsd_target_frame"]),
            "des_bridge": best[0] if best else None,
            "bridge_tc": best[1] if best else None,
            "bridge_g107": best[2] if best else None,
            "bridge_arg": best[3] if best else None,
            "arg_atom": best[4] if best else None,
            "glu_atom": best[5] if best else None,
            "pred_bridge": min((float(r[f"pred_pin_glu107_{c}"]) for c in tchains
                                if r.get(f"pred_pin_glu107_{c}")), default=None),
        })
        print(f"  {short:<10} pose {pose:6.2f} A  (J497 {want:6.2f})  "
              f"frame {items[-1]['tgt_frame']:5.2f}"
              f"{'  [FRAME_SUSPECT]' if items[-1]['suspect'] else ''}")

    # ⛔ MUST be .pml, not .py: PyMOL compiles a .py as pure Python, so the command
    # lines die with SyntaxError. As .pml it is a command script, and the `python` /
    # `python end` blocks inside it are what make the comma-safe scene() calls work.
    script = work / "j498_build.pml"
    result_file = work / "j498_result.txt"
    if result_file.exists():
        result_file.unlink()
    n_obj, n_scene = write_pymol_script(
        script, tgt_pdb, items, Path(a.out).resolve(), tchains, bchain, result_file)

    pymol = shutil.which("pymol") or "pymol"
    print(f"\nrunning {pymol} -cq {script.name}  (expect {n_obj} objects, {n_scene} scenes)")
    p = subprocess.run([pymol, "-cq", str(script)], capture_output=True, text=True)
    sys.stdout.write(p.stdout[-4000:])
    if p.returncode != 0:
        sys.stderr.write(p.stderr)
        sys.exit(f"FATAL: pymol build exited {p.returncode}")
    # ⛔⛔ pymol RETURNS 0 AFTER A FAILED scene() -- measured 2026-10-18: a comma in a
    # scene message made PyMOL's command parser read the tail as the `view` argument,
    # 12 of 14 scenes never stored, and the run still exited 0 AND wrote a 9.1 MB .pse.
    # An existence-and-size check passes that. So gate on the COUNTS the build printed.
    if "Traceback" in p.stdout or "Traceback" in p.stderr:
        sys.stderr.write((p.stdout + p.stderr)[-3000:])
        sys.exit("FATAL: the pymol build raised an exception (see above)")
    if not result_file.exists():
        sys.exit("FATAL: the build wrote no result file -- it did not reach the end")
    res = result_file.read_text().strip()
    print(f"  build result: {res}")
    if not res.startswith("ok"):
        sys.exit(f"FATAL: build refused to save -- {res}")
    out = Path(a.out)
    if not out.exists() or out.stat().st_size == 0:
        sys.exit(f"FATAL: {out} missing or empty -- pymol returned 0 but wrote nothing")
    print(f"\nwrote {out}  {out.stat().st_size/1e6:.2f} MB")
    print("⇒ NOW VERIFY IN A SEPARATE PROCESS (--verify).")


def write_pymol_script(path, tgt_pdb, items, out_pse, tchains, bchain,
                       result_file):
    """Generate a PyMOL *Python* script, not a .pml.

    ⛔⛔ WHY: in .pml syntax a scene message is parsed by PyMOL's own comma-splitting
    command parser, so ONE COMMA inside the message is read as the next argument
    (`view`) and the scene silently fails -- `ValueError: invalid literal for int()`,
    12 of 14 scenes lost, exit code still 0, .pse still written. Measured 2026-10-18.
    The cmd.* Python API takes the message as a real string argument and cannot
    misparse it. ⛔ No `if __name__ == '__main__'` guard: `pymol -cq script.py` exits 0
    writing NOTHING if the body is guarded.
    """
    L = []
    w = L.append

    def scene(key, msg):
        """⛔ Store a scene from inside a python block, NEVER as a .pml `scene` line.
        In .pml the message is comma-split by PyMOL's command parser, so one comma
        loses the scene while the run still exits 0."""
        w("python")
        w(f"cmd.scene({key!r}, 'store', message={ascii_msg(msg)!r})")
        w("python end")

    w("# J498 -- C1 v46 Boltz-2 bound folds. GENERATED; do not hand-edit.")
    w("set assembly, ''")
    w("set cartoon_transparency, 0.0")
    w("bg_color white")
    w("set ray_opaque_background, 1")
    w("set cartoon_fancy_helices, 1")
    w("set specular, 0.2")
    w(f"set_color ccny_purple, [{CCNY_PURPLE[0]:.4f},{CCNY_PURPLE[1]:.4f},{CCNY_PURPLE[2]:.4f}]")
    w(f"set_color pred_orange, [{PRED_ORANGE[0]:.4f},{PRED_ORANGE[1]:.4f},{PRED_ORANGE[2]:.4f}]")
    w(f"set_color suspect_grey, [{SUSPECT_GREY[0]:.4f},{SUSPECT_GREY[1]:.4f},{SUSPECT_GREY[2]:.4f}]")
    w(f"set_color target_grey, [{TARGET_GREY[0]:.4f},{TARGET_GREY[1]:.4f},{TARGET_GREY[2]:.4f}]")

    w(f'load {tgt_pdb}, TNFa_target')
    w("hide everything, TNFa_target")
    w("show cartoon, TNFa_target")
    w("color target_grey, TNFa_target")
    w("set cartoon_transparency, 0.55, TNFa_target")

    # GLU107 of each protomer, always visible as the thing we aimed at
    for it in items:
        if it["bridge_g107"]:
            g = it["bridge_g107"]
            w(f'select glu107_all, TNFa_target and resi {g} and chain {"+".join(tchains)}')
            w("show sticks, glu107_all")
            w("color firebrick, glu107_all")
            w("set stick_radius, 0.22, glu107_all")
            break

    for it in items:
        s = it["short"]
        w(f'load {it["des"]}, des_{s}')
        w(f'load {it["pred"]}, pred_{s}')
        for o, col in ((f"des_{s}", "ccny_purple"),
                       (f"pred_{s}", "suspect_grey" if it["suspect"] else "pred_orange")):
            w(f"hide everything, {o}")
            w(f"show cartoon, {o} and chain {bchain}")
            w(f"color {col}, {o}")
        w(f"remove des_{s} and not chain {bchain}")
        w(f"remove pred_{s} and not chain {bchain}")

    w("deselect")
    w("set cartoon_side_chain_helper, 1")

    # ---------- overview scenes ----------
    w("hide everything")
    w("show cartoon, TNFa_target")
    w("show sticks, glu107_all")
    for it in items:
        w(f'show cartoon, des_{it["short"]}')
    w("orient TNFa_target")
    w("zoom TNFa_target, 6")
    msg = (
        "DESIGNED: all 11 JOINT v46 binders in their RFd3 positions (CCNY purple). "
        "Target grey, GLU107 of each protomer red. All 11 designs share a "
        "byte-identical target frame, so these are directly comparable.")
    scene("01_designed", msg)

    w("hide everything")
    w("show cartoon, TNFa_target")
    w("show sticks, glu107_all")
    for it in items:
        w(f'show cartoon, pred_{it["short"]}')
    msg = (
        "BOLTZ-2 PREDICTED: the same 11, same frame (orange; grey = FRAME_SUSPECT, "
        "target superposed at ~21.9 A so its pose is not trustworthy). "
        "Median displacement from the designed site 41.7 A.")
    scene("02_predicted", msg)

    w("hide everything")
    w("show cartoon, TNFa_target")
    w("show sticks, glu107_all")
    for it in items:
        w(f'show cartoon, des_{it["short"]}')
        w(f'show cartoon, pred_{it["short"]}')
    msg = (
        "OVERLAY: purple = designed, orange = Boltz. The two populations barely "
        "overlap. 8 of 10 rows have a 1.5 A target frame, so the separation is real "
        "and not a superposition artefact. iptm did NOT show this: global iptm 0.695 "
        "was the target's own trimer (0.963); binder-to-target was 0.256.")
    scene("03_overlay", msg)

    # ---------- per-design scenes ----------
    for i, it in enumerate(sorted(items, key=lambda x: x["pose"]), start=1):
        s = it["short"]
        w("hide everything")
        w("show cartoon, TNFa_target")
        w("show sticks, glu107_all")
        w(f"show cartoon, des_{s}")
        w(f"show cartoon, pred_{s}")
        if it["bridge_arg"]:
            pa = it["bridge_arg"]
            w(f"select arg_des, des_{s} and resi {pa}")
            w(f"select arg_pred, pred_{s} and resi {pa}")
            w("show sticks, arg_des or arg_pred")
            w("color ccny_purple, arg_des")
            w("color pred_orange, arg_pred")
            if it["bridge_tc"] and it["bridge_g107"] and it["arg_atom"]:
                gsel = (f'TNFa_target and chain {it["bridge_tc"]} '
                        f'and resi {it["bridge_g107"]} and name {it["glu_atom"]}')
                for obj, col in ((f"des_{s}", "ccny_purple"),
                                 (f"pred_{s}", "pred_orange")):
                    nm = f"d_{obj}"
                    w(f'distance {nm}, ({obj} and resi {pa} and name {it["arg_atom"]}), '
                      f'({gsel})')
                    w(f"color {col}, {nm}")
                    w(f"set dash_width, 3.0, {nm}")
                    w(f"set label_size, 16, {nm}")
        w(f"orient des_{s} or pred_{s} or TNFa_target")
        db = f"{it['des_bridge']:.2f}" if it["des_bridge"] is not None else "n/a"
        pb = f"{it['pred_bridge']:.2f}" if it["pred_bridge"] is not None else "n/a"
        msg = (
            f"{s} ({it['arm']}): pose {it['pose']:.1f} A off the designed site. "
            f"fold RMSD {it['rmsd']:.2f} A, TM {it['tm']:.3f}. "
            f"pinned ARG{it['bridge_arg']} -> GLU107: designed {db} A, predicted {pb} A. "
            f"binder-to-target iptm {it['iptm']:.3f}. "
            f"target frame {it['tgt_frame']:.2f} A"
            + ("  *** FRAME_SUSPECT: do not read the pose ***" if it["suspect"] else ""))
        scene(f"{i+3:02d}_{s}", msg)

    w("scene 01_designed, recall")
    # ⛔⛔ ASSERT THE COUNTS **BEFORE** SAVING. A failed scene() leaves PyMOL alive and
    # exiting 0, so the only way a lost scene cannot ship is to refuse to write the
    # file. Measured 2026-10-18: 12 of 14 scenes lost, rc=0, 9.1 MB .pse written.
    n_scene = 3 + len(items)
    n_obj = 1 + 2 * len(items)
    # ⛔⛔ THE RESULT GOES TO A FILE, NOT TO STDOUT. `pymol -cq` ECHOES THE SCRIPT
    # SOURCE, so any sentinel string appears in stdout whether or not its branch ran --
    # measured 2026-10-18, a stdout grep for 'FATAL' fired off the echoed else-branch
    # of a build that had in fact SUCCEEDED. A file can only be written by execution.
    w("python")
    w(f"_want_obj, _want_scene = {n_obj}, {n_scene}")
    w("_o, _s = len(cmd.get_object_list()), len(cmd.get_scene_list())")
    w(f"_res = open({str(result_file)!r}, 'w')")
    w("if _o != _want_obj or _s != _want_scene:")
    w("    _res.write('fail objects=%d/%d scenes=%d/%d' % (_o,_want_obj,_s,_want_scene))")
    w("else:")
    w(f"    cmd.save({str(out_pse)!r})")
    w("    _res.write('ok objects=%d scenes=%d' % (_o, _s))")
    w("_res.close()")
    w("python end")
    path.write_text("\n".join(L) + "\n")
    print(f"  wrote {path} ({len(L)} lines)")
    return n_obj, n_scene


main()
