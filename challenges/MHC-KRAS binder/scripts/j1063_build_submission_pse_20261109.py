"""j1063 stage 1/4 -- build the Challenge-5 SUBMISSION .pse: 7 designs, one colour each,
against ONE copy of the KRAS G12D / HLA-A*11:01 target.

    pymol -cq j1063_build_submission_pse_20261109.py -- <picks.json> <out.pse>

PIPELINE (the Challenge-2 submission-panel pattern, `an earlier challenge's scripts/j530 -> j531 -> j575
-> j582`, which is the lab's established four-stage shape for exactly this figure):
    1. j1063_build_submission_pse        <- THIS FILE   (.pse, colours, common frame)
    2. j1063_render_submission_panels    (raw PNGs, one per design + overview, ONE camera)
    3. j1063_stamp_panels                (crop to one box, stamp the title band, re-measure)
    4. j1063_combine_panels              (the single combined submission figure)

⭐ WHAT THE FIGURE IS FOR. James McCann, 2026-11-09: "a panel per design and different coloured designs
   for each submission", in the format the earlier challenges were deposited in.
⛔⛔ THESE ARE VIEWS, NOT EVIDENCE. Nothing here ranks, scores, re-folds or re-measures anything,
   and ⛔ A MAINTAINED GEOMETRY IS NEVER EVIDENCE OF AFFINITY. Every number drawn downstream is
   CARRIED from `picks.json`, which carries it from the scored tables; this file's only arithmetic
   is the superposition it must do to put seven independent folds in one frame.

⛔⛔ THE TIER WORDS ARE NOT COSMETIC AND ARE NOT NEGOTIABLE HERE. James McCann putting all seven in the
   submission is a SCOPE decision, which is his; it does not change what any of them is.
     VALIDATED  = clears the pre-registered conjunction at 4 of 4 seeds over 20 draws  (4 designs)
     NOISE      = 3 of 4 seeds = the "INSIDE NOISE, NOT PROMOTED" tier of J1058 section 6 (3)
   ⛔ The word CANDIDATE is NOT used: in that pre-registered scheme it means a 5-draw-only design,
   and none of these are. ⛔ No design here may be called validated, promoted or a result unless
   `picks.json` says VALIDATED.

CHAIN LAYOUT IS THE TRAP AND EVERY SELECTION BELOW NAMES ITS CHAIN:
      A = HLA-A*11:01 (276 res)   B = beta-2 microglobulin (100 res)
      C = KRAS G12D 10-mer VVVGADGVGK   D = the designed binder (66-82 res)
   ⛔ This is the REFOLD layout and it is NOT the design-PDB layout (A=binder B=HLA C=b2m
   D=peptide). A mix-up would colour the HLA as the binder and still render a handsome session
   [[always-name-the-pdb-chain]].

TRAPS THIS FILE IS WRITTEN AGAINST, every one of them already paid for by this lab:
  - ⛔ `pymol -cq script.py` does NOT run the script as `__main__`, so a `__main__` guard writes
    NOTHING and exits 0 ⇒ all code here is top level [[pymol-script-not-main]].
  - ⛔ `__file__` is rebound by PyMOL to its own site-packages dir ⇒ every path comes from argv.
  - ⛔⛔ A DISABLED GROUP HIDES ITS ENABLED MEMBERS -- 12 .pse shipped showing a bare target while
    every count-based check passed (`J1024` section 1) ⇒ members AND groups are enabled and then
    ASSERTED via `cmd.get_type(n) == "object:group"`, and the check REFUSES if it finds no groups
    at all so it cannot pass vacuously [[pymol-group-disabled-hides-enabled-members]].
  - ⛔⛔ A non-ASCII character in `cmd.scene(message=)` SEGFAULTS the saved .pse ⇒ `ascii_msg()`
    strips, and it is used on EVERY message [[pymol-scene-message-must-be-ascii]].
  - ⛔⛔ `set retain_order, 1` silently corrupts `cmd.align`: a 276-residue fit collapsed to 14
    atoms at a NUMERICALLY BETTER rms than the correct 218-atom fit ⇒ retain_order is DELIBERATELY
    LEFT UNSET and `fit()` gates on ATOM COUNT, never rms. ⛔ AN RMSD CAN NEVER VALIDATE A
    SUPERPOSITION [[pymol-retain-order-corrupts-align]].
  - ⛔ `os._exit(0)` at the end skips Python's stdout flush, so on a pipe every print above is
    LOST and the run "succeeds" silently ⇒ `sys.stdout.flush()` first. (Measured here 2026-11-09:
    a probe script printed nothing at all until the flush was added.)

⭐ WHY ONE TARGET COPY. All seven are INDEPENDENT folds, so their HLA frames differ; each is
   superposed HLA-on-HLA onto the slot-1 reference and only the REFERENCE's target is kept. The
   target is then LITERALLY THE SAME OBJECT in every panel rather than seven near-copies, which is
   the strongest available form of the claim the Challenge-4 figure could only assert. The price is
   that each binder's drawn pose carries that design's own fit residual, so the fit rmsd AND the
   atom count are recorded per design and stamped on its panel.
"""
import csv
import hashlib
import json
import os
import sys

from pymol import cmd

HLA_LEN = 276          # chain A CA count, from the target spec
MIN_FIT = 200          # ⛔ the ATOM-COUNT floor. Never an rms floor.
TARGET_TRANSP = 0.25   # cartoon transparency on the pMHC, per object, AFTER load

TARGET_GREY = [0.800, 0.800, 0.815]     # chains A+B, muted so the binder reads first
PEPTIDE_PURPLE = [0x7D / 255.0, 0x55 / 255.0, 0xC7 / 255.0]   # CCNY purple #7D55C7

argv = [a for a in sys.argv[1:] if a != "--" and not a.endswith(".py")]
if len(argv) != 2:
    sys.exit(__doc__)
PICKS, OUT = argv


def ascii_msg(s):
    """⛔ A non-ASCII char in a scene message SEGFAULTS the saved .pse. Strip, never trust."""
    return s.encode("ascii", "ignore").decode("ascii")


def need(path, what):
    if not os.path.exists(path):
        sys.exit("REFUSING: missing %s: %s" % (what, path))
    return path


def fit(mobile, target, label):
    """Superpose HLA-on-HLA, then ASSERT the fit actually used the HLA.

    ⛔⛔ Gate on the ATOM COUNT, never the rms. A 14-atom fit on a 276-residue chain reports a
    BETTER rms than the correct 218-atom fit, returns no error, and would ship a figure whose
    whole visual claim was false with every other check green.
    """
    rms, n = cmd.align(mobile, target)[:2]
    print("   fit %-46s rms %6.3f A over %4d / %d CA" % (label, rms, n, HLA_LEN))
    if n < MIN_FIT:
        sys.exit("REFUSING: %s fitted only %d of %d HLA CA atoms (floor %d). A low RMSD over a "
                 "handful of atoms is not an alignment." % (label, n, HLA_LEN, MIN_FIT))
    return float(rms), int(n)


def palette_fingerprint(picks):
    """Short digest of the resolved slot -> colour assignment.

    ⭐ Carried from `j530`/`j531`: the palette is read by TWO programs (this one and the
    renderer), so this is what catches a divergence instead of two figures that silently disagree
    about which design is which colour. Stage 2 MUST print the same eight characters.
    """
    payload = ";".join("%d=%s:%s" % (d["slot"], d["color_name"], d["color_hex"]) for d in picks)
    return hashlib.md5(payload.encode("ascii")).hexdigest()[:8]


def validate_palette(picks):
    """⛔ Fail loudly rather than emit a figure whose colours mislead. Two designs sharing a
    colour would make the overview unreadable as a legend, which is its whole job."""
    names = [d["color_name"] for d in picks]
    hexes = [d["color_hex"].upper() for d in picks]
    if len(set(names)) != len(names):
        sys.exit("REFUSING: a colour NAME is reused across designs: %s" % names)
    if len(set(hexes)) != len(hexes):
        sys.exit("REFUSING: a colour VALUE is reused across designs: %s" % hexes)
    if "7D55C7" in hexes:
        sys.exit("REFUSING: CCNY purple #7D55C7 is reserved for the KRAS G12D peptide in this "
                 "figure; a design carrying it would make the peptide unreadable.")
    slots = [d["slot"] for d in picks]
    if slots != list(range(1, len(picks) + 1)):
        sys.exit("REFUSING: slots are not 1..n in order: %s" % slots)


def assert_numbers_against_tables(spec):
    """⛔⛔ THE STAMPED NUMBERS ARE ASSERTED AGAINST THE SCORED TABLES, CELL BY CELL.

    This is deliberately an EQUALITY CHECK and not a re-derivation: it re-reads each design's own
    four per-seed rows out of the table they came from and refuses if any margin, ipTM or
    conjunction flag in `picks.json` differs. It computes NOTHING, so it cannot quietly become a
    second, disagreeing scorer -- but a wrong tag, a wrong shard, a wrong rank or a transposed
    seed count cannot survive it either.

    ⛔ It is also what catches the specific error this cohort invites: the design STEM for
    `19_model_0` and `4_model_0` is `arm_shard03_...`, NOT `arm_shard00_...` as the result note's
    shorthand suggests, because those two came from rounds whose shard numbering differs.
    """
    tables = {k: list(csv.DictReader(open(need(v, "scored table")), delimiter="\t"))
              for k, v in spec["tables"].items()}
    nchecked = 0
    for d in spec["designs"]:
        rows = [r for r in tables[d["root"]] if r["tag"] == d["tag"]]
        if len(rows) != d["seeds_total"]:
            sys.exit("REFUSING: %s has %d rows in the %s table, picks.json claims %d seeds"
                     % (d["tag"], len(rows), d["root"], d["seeds_total"]))
        rows.sort(key=lambda r: r["shard"])
        got_design = {r["design"] for r in rows}
        if got_design != {d["design"]}:
            sys.exit("REFUSING: tag %s maps to design %s in the table, picks.json says %s"
                     % (d["tag"], sorted(got_design), d["design"]))
        got_rank = {int(r["rank"]) for r in rows}
        if got_rank != {d["rank"]}:
            sys.exit("REFUSING: %s is MPNN rank %s in the table, picks.json says %d"
                     % (d["tag"], sorted(got_rank), d["rank"]))
        for k, r in enumerate(rows):
            if abs(float(r["margin"]) - d["per_seed_margin"][k]) > 1e-6:
                sys.exit("REFUSING: %s seed %s margin %s != picks.json %s"
                         % (d["tag"], r["shard"], r["margin"], d["per_seed_margin"][k]))
            if abs(float(r["dpep"]) - d["per_seed_iptm"][k]) > 1e-9:
                sys.exit("REFUSING: %s seed %s ipTM %s != picks.json %s"
                         % (d["tag"], r["shard"], r["dpep"], d["per_seed_iptm"][k]))
            if (r["conj"] == "True") != bool(d["per_seed_conj"][k]):
                sys.exit("REFUSING: %s seed %s conj %s != picks.json %s"
                         % (d["tag"], r["shard"], r["conj"], d["per_seed_conj"][k]))
            nchecked += 1
        clear = sum(1 for r in rows if r["conj"] == "True")
        if clear != d["seeds_clear"]:
            sys.exit("REFUSING: %s clears %d of %d seeds in the table, picks.json says %d"
                     % (d["tag"], clear, len(rows), d["seeds_clear"]))
        want = "VALIDATED" if clear == len(rows) else "NOISE"
        if d["tier"] != want:
            sys.exit("REFUSING: %s clears %d of %d seeds, so its tier must be %s, not %s. The "
                     "tier is derived from the pre-registered conjunction, never chosen."
                     % (d["tag"], clear, len(rows), want, d["tier"]))
        print("   %-44s %s  %d/%d seeds  rank %d  [table agrees on %d cells]"
              % (d["design"], d["tier"], clear, len(rows), d["rank"], 3 * len(rows)))
    print("   asserted %d table cells against picks.json; 0 disagreements" % nchecked)


# ── read the picks, check them before touching PyMOL ──────────────────────────────────────────
spec = json.load(open(need(PICKS, "picks.json")))
picks = spec["designs"]
print("== j1063 stage 1/4: build the Challenge-5 submission .pse ==")
print("   %d designs: %d VALIDATED (4/4 seeds) + %d INSIDE NOISE (3/4 seeds, NOT promoted)"
      % (len(picks), sum(1 for d in picks if d["tier"] == "VALIDATED"),
         sum(1 for d in picks if d["tier"] == "NOISE")))
validate_palette(picks)
print("   palette fingerprint %s   (stage 2 must print the same)" % palette_fingerprint(picks))
assert_numbers_against_tables(spec)

for d in picks:
    d["pdb"] = os.path.join(spec["roots"][d["root"]], "boltz_results_shard%s" % d["shard"],
                            "predictions", d["tag"], "%s_model_0.pdb" % d["tag"])
    need(d["pdb"], "co-fold structure for %s" % d["design"])

# ── load ──────────────────────────────────────────────────────────────────────────────────────
cmd.reinitialize()
cmd.set("assembly", "")
# ⛔ retain_order is DELIBERATELY NOT SET -- it corrupts cmd.align. See fit().
cmd.bg_color("white")
cmd.set("cartoon_transparency", 0.0)

raws = []
for d in picks:
    raw = "_raw%d" % d["slot"]
    cmd.load(d["pdb"], raw)
    chains = sorted(cmd.get_chains(raw))
    if chains != ["A", "B", "C", "D"]:
        sys.exit("REFUSING: %s has chains %s, not the expected refold layout A/B/C/D"
                 % (d["design"], chains))
    nd = cmd.count_atoms("%s and chain D" % raw)
    if nd == 0:
        sys.exit("REFUSING: %s has ZERO chain-D (binder) atoms -- this is the 12-empty-.pse "
                 "signature; a session showing only the bare target passes every other check."
                 % d["design"])
    d["binder_len"] = cmd.count_atoms("%s and chain D and name CA" % raw)
    # ⛔⛔ DESIGN BASENAMES ARE NOT UNIQUE ACROSS ROUNDS. `MB_nomotif_broad_11_model_1` exists
    #    under BOTH `arm_shard03_` (2x2 round) and `arm_shard01_` (MB-depth round) -- DIFFERENT
    #    BACKBONES SHARING A STEM -- and `19_model_0` / `4_model_0` each exist under three
    #    prefixes. So a stem match can silently return the WRONG MOLECULE and every other check
    #    here would still pass. The binder length is an independent handle on identity: it is
    #    recorded per design and a mismatch is a STOP, not something to work around.
    if d["binder_len"] != d["expect_binder_len"]:
        sys.exit("REFUSING: %s loaded from %s has a %d-residue chain D, expected %d. The design "
                 "basenames collide across rounds -- this is very likely the WRONG MOLECULE."
                 % (d["design"], d["pdb"], d["binder_len"], d["expect_binder_len"]))
    raws.append((raw, d))
print("   loaded %d structures, all four chains present, all carry chain-D atoms" % len(raws))

# ── one common frame: HLA-on-HLA onto the slot-1 reference ────────────────────────────────────
ref = raws[0][0]
picks[0]["fit_rms"], picks[0]["fit_n"] = 0.0, HLA_LEN      # the reference, by definition
print("   common frame: slot 1 (%s) is the reference" % picks[0]["design"])
for raw, d in raws[1:]:
    d["fit_rms"], d["fit_n"] = fit("%s and chain A and name CA" % raw,
                                   "%s and chain A and name CA" % ref,
                                   "slot %d %s" % (d["slot"], d["design"][:28]))

# ── the ONE target copy, taken from the reference BEFORE anything is deleted ───────────────────
# ⛔ `ref` is itself a design object; taking the target copy inside the loop is how the
#    Challenge-2 panel destroyed its own reference at cell 2 of 4 (`J514` section 3 item 3).
cmd.create("target_pMHC", "%s and chain A+B" % ref)
cmd.create("peptide_G12D", "%s and chain C" % ref)
for d in picks:
    cmd.create("d%d_%s" % (d["slot"], d["design"][10:38]), "_raw%d and chain D" % d["slot"])
    d["obj"] = "d%d_%s" % (d["slot"], d["design"][10:38])
for raw, _d in raws:
    cmd.delete(raw)

npep = cmd.count_atoms("peptide_G12D and name CA")
if npep != 10:
    sys.exit("REFUSING: the peptide object holds %d residues, not the 10-mer KRAS G12D decamer"
             % npep)

# ── A CANONICAL FRAME, AND THE DEFECT THAT FORCED IT ──────────────────────────────────────────
# ⛔⛔ THE FIRST OVERVIEW RENDER PUT THE PEPTIDE UNDER THE BINDERS AND THE FIGURE CONTRADICTED ITS
#    OWN CAPTION. `cmd.orient` on peptide+designs chose a view looking roughly DOWN onto the
#    HLA groove, so all seven binders stacked directly on top of chain C and the KRAS G12D peptide
#    -- the feature the whole gate is about, and the thing the caption calls the obvious one --
#    showed as a few purple flecks. ⭐ Found by OPENING THE RAW PNG, not by any check: every
#    count, every assertion and the ink fraction all passed. Same class as the TfR1 N-terminal
#    figure, where the built segment was the least visible thing in its own panel.
# ⇒ the camera is no longer left to `orient`. A frame is DERIVED FROM THE GEOMETRY so that the
#    groove is seen EDGE-ON: the peptide runs left-right across the picture and the binders sit
#    above it, which is the view in which peptide contact is visible at all.
# ⛔⛔ AND THE COORDINATES ARE ROTATED, NOT THE CAMERA. PyMOL's 9-float view rotation is the
#    TRANSPOSE of the [right, up, back] row convention, and a hand-built camera matrix silently
#    produced sideways figures in `j695`. Rotating coordinates is convention-free and, unlike a
#    camera matrix, it can be CHECKED by re-measuring the geometry afterwards -- which is what the
#    assertions below do.
import numpy as np

def _ca(sel):
    return np.array(cmd.get_coords("%s and name CA" % sel, 1), dtype=float)


_pep = _ca("peptide_G12D")
_bind = np.vstack([_ca(d["obj"]) for d in picks])
_cp = _pep.mean(axis=0)
# x-screen: the peptide's own long axis, from an SVD of its CA coordinates
_u = np.linalg.svd(_pep - _cp)[2][0]
_u /= np.linalg.norm(_u)
# y-screen: from the peptide towards the binders, orthogonalised against the peptide axis.
# ⇒ the binders are ABOVE the peptide on screen BY CONSTRUCTION, not by eye.
_b = _bind.mean(axis=0) - _cp
_n = _b - np.dot(_b, _u) * _u
_n /= np.linalg.norm(_n)
_w = np.cross(_u, _n)
TTT = [_u[0], _u[1], _u[2], 0.0,
       _n[0], _n[1], _n[2], 0.0,
       _w[0], _w[1], _w[2], 0.0,
       0.0, 0.0, 0.0, 1.0]
for _o in cmd.get_object_list():
    cmd.transform_object(_o, TTT)

# ⛔ RE-MEASURE IN THE NEW FRAME rather than trusting the matrix.
_pep2, _bind2 = _ca("peptide_G12D"), np.vstack([_ca(d["obj"]) for d in picks])
_spread = _pep2.max(axis=0) - _pep2.min(axis=0)
print("   canonical frame: peptide spread x/y/z = %.1f / %.1f / %.1f A; binder centroid y %+.1f "
      "vs peptide centroid y %+.1f" % (_spread[0], _spread[1], _spread[2],
                                       _bind2.mean(axis=0)[1], _pep2.mean(axis=0)[1]))
if not (_spread[0] > _spread[1] and _spread[0] > _spread[2]):
    sys.exit("REFUSING: the peptide is not the widest along screen-x after the rotation -- the "
             "frame is wrong and the groove would not be edge-on.")
if not _bind2.mean(axis=0)[1] > _pep2.mean(axis=0)[1]:
    sys.exit("REFUSING: the binders are not above the peptide after the rotation -- the figure "
             "would be upside down.")
print("   one target copy kept: %d HLA CA + %d b2m CA + %d peptide CA"
      % (cmd.count_atoms("target_pMHC and chain A and name CA"),
         cmd.count_atoms("target_pMHC and chain B and name CA"), npep))

# ── styling. ⛔ transparency PER OBJECT and AFTER load, or it applies to nothing ───────────────
cmd.hide("everything")
cmd.set("cartoon_fancy_helices", 1)
cmd.set_color("c5_target_grey", TARGET_GREY)
cmd.set_color("ccny_purple", PEPTIDE_PURPLE)

cmd.show("cartoon", "target_pMHC")
cmd.color("c5_target_grey", "target_pMHC")
cmd.set("cartoon_transparency", TARGET_TRANSP, "target_pMHC")

# ⭐ THE PEPTIDE IS THE SUBJECT, SO IT IS THE LOUDEST THING IN THE PICTURE. Peptide contact is the
#    entire basis of the pre-registered gate, so chain C gets cartoon AND sticks, fully opaque, in
#    CCNY purple -- one colour, the same in every panel.
cmd.show("cartoon", "peptide_G12D")
cmd.show("sticks", "peptide_G12D")
cmd.color("ccny_purple", "peptide_G12D")
cmd.set("cartoon_transparency", 0.0, "peptide_G12D")
cmd.set("stick_radius", 0.30, "peptide_G12D")       # fat enough to read behind a binder helix
cmd.set("cartoon_tube_radius", 0.40, "peptide_G12D")

for d in picks:
    cmd.set_color(d["color_name"], [int(d["color_hex"][i:i + 2], 16) / 255.0 for i in (0, 2, 4)])
    cmd.show("cartoon", d["obj"])
    cmd.color(d["color_name"], d["obj"])
    cmd.set("cartoon_transparency", 0.0, d["obj"])

# ── groups. ⛔⛔ A DISABLED GROUP HIDES ITS ENABLED MEMBERS ⇒ enable BOTH, then assert ──────────
cmd.group("designs", " ".join(d["obj"] for d in picks))
cmd.group("target", "target_pMHC peptide_G12D")
for d in picks:
    cmd.enable(d["obj"])
cmd.enable("target_pMHC")
cmd.enable("peptide_G12D")
cmd.enable("designs")
cmd.enable("target")

groups = [n for n in cmd.get_names("all") if cmd.get_type(n) == "object:group"]
if not groups:
    sys.exit("REFUSING: no groups found. The group check must not be able to pass vacuously -- "
             "that is how the first fix for this defect reported success while doing nothing.")
not_enabled = [g for g in groups if g not in cmd.get_names("all", enabled_only=1)]
if not_enabled:
    sys.exit("REFUSING: groups %s are DISABLED. A disabled group hides its enabled members and "
             "the session would open as a bare target." % not_enabled)
enabled = set(cmd.get_names("all", enabled_only=1))
missing = [d["obj"] for d in picks if d["obj"] not in enabled]
if missing:
    sys.exit("REFUSING: binder objects not enabled: %s" % missing)
print("   visibility ASSERTED: %d groups %s all enabled, %d/%d binder objects enabled"
      % (len(groups), groups, len(picks) - len(missing), len(picks)))

# ── scenes. ⛔ ASCII ONLY in every message ─────────────────────────────────────────────────────
def frame(sel, buf):
    """⛔ IDENTITY ROTATION ON THE CANONICAL FRAME, never `orient`. `orient` is what chose the
    down-the-groove view that hid the peptide; the frame above already encodes which way is up,
    so the camera's only remaining job is to zoom."""
    v = list(cmd.get_view())
    v[0:9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    cmd.set_view(v)
    cmd.zoom(sel, buf)
    cmd.clip("slab", 400)        # ⛔ widen the slab after EVERY camera change


frame("target_pMHC or peptide_G12D or designs", 4.0)
cmd.scene("00_overview", "store", message=ascii_msg(
    "CHALLENGE 5 SUBMISSION SET -- 7 designs vs KRAS G12D 10-mer on HLA-A*11:01. "
    "GREY = HLA + b2m (ONE copy, the slot-1 reference; all 7 superposed HLA-on-HLA onto it). "
    "CCNY PURPLE = the KRAS G12D peptide. One colour per design. "
    "4 designs are VALIDATED (4 of 4 seeds, 20 draws); 3 clear only 3 of 4 seeds and are INSIDE "
    "NOISE, NOT PROMOTED. These are VIEWS, not evidence: a maintained geometry is never "
    "evidence of affinity and nothing here has been tested."))
for d in picks:
    for o in [x["obj"] for x in picks]:
        cmd.disable(o)
    cmd.enable(d["obj"])
    frame("target_pMHC or peptide_G12D or %s" % d["obj"], 4.0)
    tier = ("VALIDATED -- clears the pre-registered conjunction at 4 of 4 seeds over 20 draws"
            if d["tier"] == "VALIDATED" else
            "3 of 4 SEEDS -- INSIDE NOISE, NOT PROMOTED (J1058 s6 tiering); the %s limb failed "
            "at one seed" % d["failed_limb"])
    cmd.scene("%02d_slot%d" % (d["slot"], d["slot"]), "store", message=ascii_msg(
        "%s  (submission name %s) | generated in the %s | arm %s | MPNN rank %d | %d aa | "
        "mean contact margin %s over 4 seeds | mean binder-peptide ipTM %s | "
        "HLA-on-HLA fit %.3f A over %d/%d CA | %s"
        % (d["design"], d["sub_name"], d["origin_round"], d["arm"], d["rank"], d["binder_len"],
           d["margin_mean"], d["iptm_mean"], d["fit_rms"], d["fit_n"], HLA_LEN, tier)))
for d in picks:
    cmd.enable(d["obj"])
cmd.scene("00_overview", "recall")

# ── save ──────────────────────────────────────────────────────────────────────────────────────
os.makedirs(os.path.dirname(OUT), exist_ok=True)
cmd.save(OUT)
if not os.path.exists(OUT) or os.path.getsize(OUT) < 10000:
    sys.exit("REFUSING: %s was not written or is implausibly small" % OUT)

man = {"ok": True,
       "pse": os.path.abspath(OUT),
       "picks": os.path.abspath(PICKS),
       "palette_fingerprint": palette_fingerprint(picks),
       "target_object": "target_pMHC",
       "peptide_object": "peptide_G12D",
       "target_transparency": TARGET_TRANSP,
       "reference_slot": 1,
       "min_fit_atoms": MIN_FIT,
       "hla_len": HLA_LEN,
       "groups": groups,
       "designs": [{k: d[k] for k in
                    ("slot", "design", "sub_name", "origin_round", "tag", "root", "shard",
                     "rank", "arm", "tier", "seeds_clear", "seeds_total", "margin_mean",
                     "iptm_mean", "sd_pep_contacts", "failed_limb", "color_name", "color_hex",
                     "obj", "pdb", "binder_len", "fit_rms", "fit_n")} for d in picks]}
mp = os.path.splitext(OUT)[0] + "_build.json"
json.dump(man, open(mp, "w"), indent=1)
print("   wrote %s  (%.1f MB, %d objects, %d scenes)"
      % (OUT, os.path.getsize(OUT) / 1e6, len(cmd.get_names("objects")),
         len(cmd.get_scene_list())))
print("   wrote %s" % mp)
print("== j1063 stage 1 ok ==")

# ⛔ PyMOL can segfault in its own atexit teardown even after everything passed; os._exit skips
#    it. ⛔⛔ But os._exit ALSO skips Python's stdout flush, so without the flush every line above
#    is lost on a pipe and the run "succeeds" printing nothing. Measured here 2026-11-09.
sys.stdout.flush()
os._exit(0)
