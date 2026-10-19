#!/usr/bin/env python3
"""j532 -- recompute the TARGET FRAME and BINDER POSE over ALL chain permutations.

⛔⛔ THE BUG THIS FIXES, measured 2026-10-18 ~23:1x on the real `c1_holo` cohort.
`j497` superposes the predicted target onto the designed target and keeps the best of
**three CYCLIC** chain mappings (`tchains[k:] + tchains[:k]`), on the stated grounds that
the target is a C3 trimer so a binder at a symmetry-equivalent site must not read as a
failure. That reasoning is right; the enumeration is incomplete. There are **SIX** ways to
map three chains, and Boltz assigns chain IDs per run with no guarantee of cyclic order.

**34 of 88 `c1_holo` designs scored `FRAME_SUSPECT` at 21.78-21.93 A** -- a 0.15 A
spread across 34 designs, which is the signature of ONE systematic alternative, not of 34
bad folds. Under the three ODD permutations those same designs superpose at **~1.5 A**:

    OK       BCD*=1.47  BDC =21.95  CBD =21.92  CDB*=1.49  DBC*=1.46  DCB =21.93
    SUSPECT  BCD*=21.85 BDC = 1.55  CBD = 1.56  CDB*=21.87 DBC*=21.86 DCB = 1.53
    (* = tried by j497)

⇒ **all 34 are fine**; 39% of the C1 pool was being discarded for a bookkeeping reason.
Three independent checks ruled out the alternatives before this conclusion was drawn:
  1. the predicted trimer geometry is the SAME in both groups (centroid distances 21.5 vs
     21.2-21.3 A, Rg 21.2 both) -- Boltz did NOT mis-assemble the target;
  2. every backbone is MIXED (OK and SUSPECT sequences share one design file) -- so it is
     not a property of the design;
  3. all three cyclic maps agree to within 0.05 A on every row -- so the cyclic search was
     not merely picking a poor member.

⛔⛔ WHY THIS IS A SEPARATE SCRIPT AND NOT AN EDIT TO `j497`.
`j521` imports `j497` **under a SHA1 pin** and exits FATAL on drift. Editing `j497` with the
C3 arms minutes from landing would have broken the C3 scorer at the worst moment. `J512`
already settled this exact trade when `j465` needed changing: build alongside, do not edit
the pinned thing mid-deadline. ⛔ `j497` and `j521` are untouched by this file.

⭐ SELECTION RULE, which preserves `j497`'s intent rather than maximising a number:
a permutation is a CANDIDATE only if its target RMSD is under `--frame-max` (the frame must
be real). Among candidates -- which are the genuinely symmetry-equivalent framings -- take
the one that best places the binder. ⛔ Never minimise binder pose across ALL permutations:
that would let a bad frame flatter a badly-placed binder, which is the artefact `j497`'s own
comment warns about.

Emits `*_v2` columns beside the originals and never overwrites them, so the correction is
auditable and the original verdict stays readable.

usage:
  python3 scripts/j532_frame_permutation_fix_20261018.py \
      --in-tsv data/j519_c1holo_selfconsistency_20261018.tsv \
      --pred-dir boltz/.../out_c1_holo/boltz_results_c1_holo/predictions \
      --design-dir mpnn/j513_c1_rebias_20261018/in \
      --binder-chains A --target-chains B,C,D \
      --out data/j532_c1holo_framefix_20261018.tsv
"""
import argparse
import csv
import itertools
import os
import sys

import numpy as np


def load_j497():
    """Reuse j497's own readers so the comparison is like-for-like, not a reimplementation."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "j497_c1_boltz_selfconsistency_20261018.py")
    spec = importlib.util.spec_from_file_location("j497_mod", path)
    mod = importlib.util.module_from_spec(spec)
    saved, sys.argv = sys.argv, ["j497_mod"]
    try:
        spec.loader.exec_module(mod)
    except SystemExit:
        pass
    finally:
        sys.argv = saved
    return mod


def kabsch(P, Q):
    """Proper-rotation superposition. Returns (rmsd, R, Pmean, Qmean)."""
    Pm, Qm = P.mean(0), Q.mean(0)
    Pc, Qc = P - Pm, Q - Qm
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return float(np.sqrt((((Pc @ R) - Qc) ** 2).sum() / len(P))), R, Pm, Qm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-tsv", required=True)
    ap.add_argument("--pred-dir", required=True)
    ap.add_argument("--design-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--binder-chains", default="A")
    ap.add_argument("--target-chains", default="B,C,D")
    ap.add_argument("--design-binder-chains", default=None)
    ap.add_argument("--design-target-chains", default=None)
    ap.add_argument("--frame-max", type=float, default=5.0,
                    help="target RMSD above which a framing is not real (j497 uses 5.0)")
    a = ap.parse_args()

    if os.path.exists(a.out):
        sys.exit(f"REFUSING: {a.out} exists. Fresh --out per run : scripts APPEND "
                 f"or merge silently and an analysis then averages two campaigns.")

    bchains = [c for c in a.binder_chains.split(",") if c]
    tchains = [c for c in a.target_chains.split(",") if c]
    dbchains = ([c for c in a.design_binder_chains.split(",") if c]
                if a.design_binder_chains else list(bchains))
    dtchains = ([c for c in a.design_target_chains.split(",") if c]
                if a.design_target_chains else list(tchains))
    assert len(dtchains) == len(tchains), (dtchains, tchains)
    assert len(dbchains) == len(bchains), (dbchains, bchains)

    j = load_j497()
    rows = list(csv.DictReader(open(a.in_tsv), delimiter="\t"))
    print(f"RESOLVED  binder={','.join(bchains)}  target={','.join(tchains)}  "
          f"design-target={','.join(dtchains)}  frame-max={a.frame_max}")
    print(f"  {len(rows)} rows from {a.in_tsv}")
    perms = list(itertools.permutations(range(len(tchains))))
    print(f"  {len(perms)} permutations per row "
          f"(j497 tried {len(tchains)}: the cyclic ones only)")

    n_rescued = n_still = n_err = 0
    out_rows = []
    for r in rows:
        name = r.get("design")
        r2 = dict(r)
        ppath = os.path.join(a.pred_dir, name, f"{name}_model_0.pdb")
        dpath = os.path.join(a.design_dir, r.get("design_file", ""))
        try:
            pred = j.read_atoms(ppath, bchains + tchains)
            des = j.read_atoms(dpath, dbchains + dtchains)
            Q = np.vstack([des[c]["ca"] for c in dtchains])
            dB = np.vstack([des[c]["ca"] for c in dbchains])
            B = np.vstack([pred[c]["ca"] for c in bchains])
            cands = []
            for pm in perms:
                P = np.vstack([pred[tchains[i]]["ca"] for i in pm])
                if P.shape != Q.shape:
                    continue
                trms, R, Pmean, Qmean = kabsch(P, Q)
                if B.shape != dB.shape:
                    pose = None
                else:
                    pose = float(np.sqrt((((B - Pmean) @ R - (dB - Qmean)) ** 2).sum()
                                         / len(B)))
                cands.append(("".join(tchains[i] for i in pm), trms, pose))
            if not cands:
                raise ValueError("no permutation produced matching atom counts")
            # ⭐ frame validity FIRST, then best binder placement among real framings
            good = [c for c in cands if c[1] <= a.frame_max]
            pool = good if good else cands
            key = (lambda c: (c[2] if c[2] is not None else 9e9)) if good else (lambda c: c[1])
            best = min(pool, key=key)
            r2["pose_frame_v2"] = best[0]
            r2["rmsd_target_frame_v2"] = round(best[1], 3)
            r2["rmsd_binder_on_target_frame_v2"] = (
                round(best[2], 3) if best[2] is not None else "")
            r2["n_valid_frames_v2"] = len(good)
            r2["pose_status_v2"] = ("OK" if good
                                    else f"FRAME_SUSPECT target_rmsd={round(best[1], 3)}")
            r2["frame_was_cyclic_v2"] = str(
                best[0] in {"".join(tchains[k:] + tchains[:k]) for k in range(len(tchains))})
            if str(r.get("pose_status", "")).startswith("FRAME_SUSPECT") and good:
                n_rescued += 1
            elif not good:
                n_still += 1
            out_rows.append(r2)
        except Exception as exc:                                   # noqa: BLE001
            n_err += 1
            r2["pose_status_v2"] = f"ERROR {type(exc).__name__}: {exc}"
            for c in ("pose_frame_v2", "rmsd_target_frame_v2",
                      "rmsd_binder_on_target_frame_v2", "n_valid_frames_v2",
                      "frame_was_cyclic_v2"):
                r2.setdefault(c, "")
            out_rows.append(r2)

    cols = list(rows[0].keys()) + ["pose_status_v2", "pose_frame_v2", "rmsd_target_frame_v2",
                                   "rmsd_binder_on_target_frame_v2", "n_valid_frames_v2",
                                   "frame_was_cyclic_v2"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    was = sum(1 for r in rows if str(r.get("pose_status", "")).startswith("FRAME_SUSPECT"))
    now = sum(1 for r in out_rows if str(r.get("pose_status_v2", "")).startswith("FRAME_SUSPECT"))
    ncyc = sum(1 for r in out_rows if r.get("frame_was_cyclic_v2") == "False")
    print(f"\n  wrote {a.out}  ({len(out_rows)} rows)")
    print(f"  FRAME_SUSPECT  was {was}  ->  now {now}   (RESCUED {n_rescued})")
    print(f"  rows whose winning frame is NON-cyclic (invisible to j497): {ncyc}")
    if n_err:
        print(f"  ⛔ {n_err} rows errored -- see pose_status_v2")
    if now:
        print(f"  ⛔ {now} rows remain FRAME_SUSPECT under ALL {len(perms)} permutations: "
              f"those are real, not bookkeeping.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
