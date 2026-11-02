#!/usr/bin/env python3
"""j901 prep -- select the best screen designs and MEASURE, per design, which residues RFd3 was
GIVEN (held 5XOX elements), which are the fixed pins, and which it BUILT.

Local. ⛔ Runs no model. Writes gunzipped .cif copies plus one JSON manifest for the PyMOL
builder; it does NOT open PyMOL.

WHY THE HELD RESIDUES ARE MEASURED AND NOT READ OFF THE CONTIG
--------------------------------------------------------------
The output .cif carries `is_motif_atom_with_fixed_seq`, which marks ONLY the four bare guideposts
-- the held element SEQUENCE was deliberately released (`select_unfixed_sequence`), so the held
elements are NOT flagged in the file at all. ⛔ And the contig cannot be read off either: the
linker lengths are RANGES drawn per batch, so the output index of a held element differs
design to design.

⇒ The held residues are identified GEOMETRICALLY, which also re-verifies the construct:
`j897` measured `ori_token` to be consumed as a PURE TRANSLATION at exactly 0.0000 A, so a held
residue's output CA must sit at `input_CA - ori_token`. Each match is asserted under MATCH_TOL and
the COUNT is asserted equal to the input's own held-chain residue count. MEASURED: on opt05 the
worst match is 0.0621 A over 43/43 residues, so the motif really is held, not merely requested.

⛔⛔ PSEUDO-ATOMS: raw RFd3 output can carry V0-V8 pseudo-atoms, and a TEXT match returns ZERO and
reads as a false clean. The check here matches
type_symbol == "V" AND label_atom_id == V<digits>, and the count is PRINTED per design rather than
assumed to be zero.

usage: python3 -I j901_prep_pse_inputs_20261101.py --outdir <dir>
"""
import argparse
import gzip
import json
import os
import re
import shutil
import sys

import numpy as np

# ⚠️ SET THESE THREE BEFORE RUNNING. ARMS and GATE resolve inside this directory as staged;
# DESIGNS is the unpacked generator output, which is NOT distributed here (it is tens of
# thousands of coordinate files). See ../methods/README.md, "Reproducibility".
HERE = os.path.dirname(os.path.abspath(__file__))
ARMS = os.path.join(HERE, "..", "data", "core_arm_specs_20261101")
DESIGNS = os.path.join(HERE, "..", "..", "..", "screen_designs")
GATE = os.path.join(HERE, "..", "data", "j900_screen_gate_20261101.tsv")

MATCH_TOL = 0.25          # A; the measured worst case is ~0.06
VPSEUDO = re.compile(r"^V\d+$")
NAME = re.compile(r"arm_shard(\d+)_tfr1_core_(opt\d+)_(\d+)_model_(\d+)$")

# ⭐ SIX designs, chosen from the j900 conjunction survivors across SIX DISTINCT BATCHES and three
#    options, so the session is not six views of one contig draw. Ranked within an option by
#    clearance to target then fold integrity (min_cut).
PICK = [
    ("arm_shard04_tfr1_core_opt05_2_model_6", "opt05", "best arm, best design: cut 44, 0 res under 2.5 A"),
    ("arm_shard04_tfr1_core_opt05_0_model_0", "opt05", "same arm, independent batch 0"),
    ("arm_shard04_tfr1_core_opt05_3_model_4", "opt05", "same arm, independent batch 3; cut 50"),
    ("arm_shard00_tfr1_core_opt01_1_model_5", "opt01", "best clearance in the whole cohort, 3.44 A"),
    ("arm_shard00_tfr1_core_opt01_4_model_5", "opt01", "most compact survivor, rg_ratio 1.056"),
    ("arm_shard03_tfr1_core_opt04_4_model_7", "opt04", "third arm, for contrast"),
]

ATOM_COLS = dict(name=2, resn=4, auth_resi=9, auth_chain=11, elem=1,
                 fixed_seq=17, x=18, y=19, z=20)


def read_gate():
    rows, hdr = {}, None
    for ln in open(GATE):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            continue
        rows[f[0]] = dict(zip(hdr, f))
    return rows


def input_held_ca(opt):
    """-> {5XOX_resi: CA coord} for chain H of the arm's input PDB, plus the spec's ori_token."""
    spec = json.load(open(os.path.join(ARMS, "tfr1_core_%s_spec.json" % opt)))
    _k, arm = next(iter(spec.items()))
    ori = np.asarray(arm["ori_token"], float)
    held = {}
    for ln in open(os.path.join(ARMS, arm["input"])):
        if ln.startswith("ATOM") and ln[21] == "H" and ln[12:16].strip() == "CA":
            held[int(ln[22:26])] = np.array(
                [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return held, ori


def parse_out(path):
    """-> (list of CA records for the binder chain, n_vpseudo, chains seen)"""
    ca, nv, chains = [], 0, set()
    C = ATOM_COLS
    for ln in gzip.open(path, "rt"):
        p = ln.split()
        if len(p) <= C["z"] or p[0] != "ATOM":
            continue
        chains.add(p[C["auth_chain"]])
        # ⛔ element AND name, never a text match on the line
        if p[C["elem"]] == "V" and VPSEUDO.match(p[C["name"]]):
            nv += 1
        if p[C["name"]] == "CA":
            ca.append(dict(resi=int(p[C["auth_resi"]]), chain=p[C["auth_chain"]],
                           resn=p[C["resn"]], fixed=p[C["fixed_seq"]] == "True",
                           xyz=np.array([float(p[C["x"]]), float(p[C["y"]]),
                                         float(p[C["z"]])])))
    return ca, nv, chains


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    gate = read_gate()

    print("== j901 prep: %d designs ==" % len(PICK))
    nfail = 0
    manifest = []
    for design, opt, why in PICK:
        src = os.path.join(DESIGNS, design + ".cif.gz")
        if not os.path.exists(src):
            print("  FAIL  %s: source missing" % design)
            nfail += 1
            continue
        held_in, ori = input_held_ca(opt)
        ca, nv, chains = parse_out(src)
        binder = [r for r in ca if r["chain"] == "A"]
        target = sorted({r["chain"] for r in ca} - {"A"})

        # --- MEASURE the held residues: output CA == input CA - ori_token
        tgt = {hr: (hc - ori) for hr, hc in held_in.items()}
        held_out, devs, used = [], [], set()
        for r in binder:
            best = min(((hr, float(np.linalg.norm(r["xyz"] - hc)))
                        for hr, hc in tgt.items()), key=lambda t: t[1])
            if best[1] < MATCH_TOL:
                held_out.append(r["resi"])
                devs.append(best[1])
                used.add(best[0])
        pins = [r["resi"] for r in binder if r["fixed"]]
        held_set, pin_set = set(held_out), set(pins)
        built = [r["resi"] for r in binder
                 if r["resi"] not in held_set and r["resi"] not in pin_set]

        ok_count = len(used) == len(held_in)
        g = gate.get(design, {})
        print("  %s  %-40s %s" % ("PASS" if ok_count else "FAIL", design, opt))
        print("        held %d/%d matched (worst dev %.4f A) | pins %d | built %d | %d res | V-pseudo %d | chains %s"
              % (len(used), len(held_in), max(devs) if devs else -1, len(pins),
                 len(built), len(binder), nv, ",".join(sorted(chains))))
        if not ok_count:
            nfail += 1
            continue
        if nv:
            print("        ⛔ %d V-pseudo-atoms present -- the builder MUST exclude them" % nv)

        dst = os.path.join(a.outdir, design + ".cif")
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        manifest.append(dict(
            design=design, opt=opt, why=why, cif=os.path.basename(dst),
            held=sorted(held_out), pins=sorted(pins), built=sorted(built),
            target_chains=target, n_res=len(binder), n_vpseudo=nv,
            worst_held_dev=round(max(devs), 4) if devs else None,
            max_step=g.get("max_step", ""), bb_min=g.get("bb_min_to_target", ""),
            min_cut=g.get("min_cut", ""), rg_ratio=g.get("rg_ratio", ""),
            n_under_2p5=g.get("n_res_under_2p5", "")))

    mp = os.path.join(a.outdir, "j901_manifest.json")
    json.dump(manifest, open(mp, "w"), indent=1)
    print("\n== %d of %d prepared, %d FAIL ==" % (len(manifest), len(PICK), nfail))
    print("   manifest: %s" % mp)
    if nfail:
        sys.exit(1)


main()
