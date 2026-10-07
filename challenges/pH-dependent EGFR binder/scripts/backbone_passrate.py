"""BACKBONE-LEVEL pass rate -- the primary metric.

Not per sequence. For each BACKBONE, pool all of its draws (100 sequences x 3 draws = 300) and
ask how many passed the epitope-contact gate.

⛔ Backbones SPAN shards (foldset.tsv is length-sorted and shards are sequential chunks), so a
backbone's 300 draws can be split across two shards. This aggregates across every shard and
asserts the total, because a backbone silently scored on a partial draw set reads TOO HIGH --
so never rank on partial data. That has already reversed a winner on this project.

GATE: epitope CONTACT is primary; `pose_rmsd` is reported as a diagnostic and is NEVER a cut
(on the 8x measured swing). Thresholds are reported at several levels rather
than one, because the threshold is a choice and should be made in the open.

usage: backbone_passrate.py <contacts.csv> [out.csv]
"""
import sys, csv, collections
import numpy as np

EXPECT_DRAWS = 300          # 100 sequences x 3 draws per backbone


def med(v):
    v = [x for x in v if x == x]
    return float(np.median(v)) if v else float("nan")


def main(inp, outp=None):
    rows = list(csv.DictReader(open(inp)))
    for r in rows:
        r["n_epitope"] = int(r["n_epitope"])
        r["n_epitope_total"] = int(r["n_epitope_total"])
        r["n_epitope_alt"] = int(r["n_epitope_alt"]) if r.get("n_epitope_alt") else None
        for k in ("pose_rmsd", "sc_rmsd", "tgt_fit_rmsd"):
            try:
                r[k] = float(r[k])
            except (TypeError, ValueError):
                r[k] = float("nan")

    by = collections.defaultdict(list)
    for r in rows:
        by[(r["site"], r["design"])].append(r)

    out = []
    for (site, bb), rs in by.items():
        tot = rs[0]["n_epitope_total"]
        n = len(rs)
        ne = [r["n_epitope"] for r in rs]
        # thresholds as fractions of the site's own footprint
        full = sum(1 for x in ne if x >= tot)
        t90 = sum(1 for x in ne if x >= int(np.ceil(0.90 * tot)))
        t80 = sum(1 for x in ne if x >= int(np.ceil(0.80 * tot)))
        t50 = sum(1 for x in ne if x >= int(np.ceil(0.50 * tot)))
        rec = dict(
            site=site, backbone=bb, n_draws=n,
            complete="YES" if n == EXPECT_DRAWS else f"NO({n})",
            shards=",".join(sorted({r["shard"] for r in rs})),
            foot=tot,
            pass_full=full, pct_full=100.0 * full / n,
            pass_90=t90, pct_90=100.0 * t90 / n,
            pass_80=t80, pct_80=100.0 * t80 / n,
            pass_50=t50, pct_50=100.0 * t50 / n,
            med_contacts=med(ne),
            med_pose_rmsd=med([r["pose_rmsd"] for r in rs]),
            med_sc_rmsd=med([r["sc_rmsd"] for r in rs]),
        )
        if site == "B":
            # alt = FP8, the canonical 8 incl. Q28. STRETCH diagnostic only -- the gate is FP7
            # because FP7 is the cut these designs were selected under.
            alt = [r["n_epitope_alt"] for r in rs if r["n_epitope_alt"] is not None]
            rec["pass_FP8_stretch"] = sum(1 for x in alt if x >= 8) if alt else ""
        out.append(rec)

    out.sort(key=lambda r: (-r["pct_90"], -r["pct_full"]))

    print(f"BACKBONE-LEVEL PASS RATE -- {len(out)} backbones, "
          f"{sum(r['n_draws'] for r in out)} draws total\n")
    print(f"{'site':>4} {'backbone':44s} {'draws':>6} {'cmpl':>6} {'foot':>4} "
          f"{'ALL':>11} {'>=90%':>11} {'>=80%':>11} {'medC':>5} {'medPose':>8} {'medSC':>6}")
    for r in out:
        print(f"{r['site']:>4} {r['backbone']:44s} {r['n_draws']:6d} {r['complete']:>6} "
              f"{r['foot']:4d} "
              f"{r['pass_full']:5d} {r['pct_full']:5.1f}% "
              f"{r['pass_90']:5d} {r['pct_90']:5.1f}% "
              f"{r['pass_80']:5d} {r['pct_80']:5.1f}% "
              f"{r['med_contacts']:5.1f} {r['med_pose_rmsd']:8.2f} {r['med_sc_rmsd']:6.2f}")

    inc = [r for r in out if r["complete"] != "YES"]
    if inc:
        print("\n⛔ INCOMPLETE BACKBONES -- do NOT rank these against the complete ones:")
        for r in inc:
            print(f"   {r['site']} {r['backbone']}  {r['n_draws']}/{EXPECT_DRAWS} draws  "
                  f"shards {r['shards']}")
    else:
        print(f"\n✅ every backbone has its full {EXPECT_DRAWS} draws")

    if outp:
        # fieldnames must be the UNION: site-B rows carry pass_FP8_stretch and site-A rows do not,
        # so taking out[0]'s keys drops a column or raises, depending which site sorts first.
        cols = list(out[0].keys())
        for r in out:
            for k in r:
                if k not in cols:
                    cols.append(k)
        with open(outp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, restval="")
            w.writeheader()
            w.writerows(out)
        print(f"\nwrote {outp}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
