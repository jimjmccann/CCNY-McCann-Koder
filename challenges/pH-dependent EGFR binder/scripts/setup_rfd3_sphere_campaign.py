#!/usr/bin/env python3
"""Build the RFd3 SPHERE campaign: both sites, 5 lengths each, 3 shards per length.

⭐⭐ WHY A SPHERE AND NOT A CROP. MEASURED 2026-10-02: RFd3 holds a NON-CONTIGUOUS target RIGID
-- 21 scattered segments kept their mutual geometry to 0.0402 A, 978/978 atoms out, nothing
invented in the gaps, no linking. So the earlier conclusion that no contiguous crop works is not a
dead end: contiguity was never the requirement.
A sphere round the footprint keeps far more real receptor than any contiguous window could.

⛔⛔ THE SPHERE DOES NOT REPLACE THE FULL-RECEPTOR GATE. Trimming preserves coordinates, not
physics: a residue whose neighbours were cut is artificially exposed. The radius is a GUESS and the
gate is the only thing that proves it was generous enough. RFd3 also RECENTERS output by ~67 A, so
the gate MUST superpose, never diff raw coordinates.

RADII (chosen generously): site A 22 A = 231 res; site B 20 A = 230 res. Chosen to match each other
so both arms run at the same rate, and to stay under RFd3's 384-residue training crop with the
binder included (A: 231+80 = 311; B: 230+70 = 300).

LENGTHS. Site A 60-80 by 5; site B 50-70 by 5 (site B is smaller, so its window sits 10 aa lower).
⛔ The length lives in the contig's LEADING INTEGER, not in a `length` field -- a `length` key
would be silently ignored by the thing that matters (make_rfd3_shards.py's standing warning).

⛔ THREE SHARDS PER LENGTH IS ONLY SAFE BECAUSE INFERENCE IS UNSEEDED. Verified in source:
`seed` defaults None in rfd3/engine.py and is only written into output metadata (L302);
`seed_everything` appears ONLY in train.py and testing/, never in the inference path. If it were
seeded, three shards at one length would be byte-identical and we would pay 3x for one shard.

⛔ ONE .pdb PER BUNDLE. The generation runner refuses a bundle whose in/ directory holds more
than one target pdb, because it cannot tell which is the target. => TWO BUNDLES, one per site.

⛔ Hotspots are the PROVEN sets from the earlier single-site arms, not invented. Every one is
ASSERTED BY RESIDUE TYPE against the input file before anything is written, because a wrong
numbering has silently produced a null result before.

CFG stays OFF (refuted: 2.76x null). n_batches/overrides are CLI-only, NOT spec fields.
2026-10-02.
"""
import argparse, json, os, shutil, subprocess, sys
import numpy as np

REPO = "/PATH/TO/CHECKOUT"
MODEL = f"{REPO}/EGFR_receptor_WTnumbering_tethered.pdb"
RFD3_PY = "/PATH/TO/rfd3-venv/bin/python"

SITES = {
    "A": dict(radius=22.0, lengths=(60, 65, 70, 75, 80),
              footprint=[380, 382, 384, 408, 409, 410, 411, 412, 417, 438, 465],
              hotspots={409: ("HIS", "ND1,NE2"),
                        412: ("PHE", "CG,CD1,CD2,CE1,CE2,CZ"),
                        438: ("ILE", "CG2,CD1")}),
    "B": dict(radius=20.0, lengths=(50, 55, 60, 65, 70),
              footprint=[19, 20, 21, 25, 29, 30, 50],
              hotspots={29: ("ARG", "NE,NH1,NH2"),
                        20: ("PHE", "CG,CD1,CD2,CE1,CE2,CZ"),
                        25: ("LEU", "CD1,CD2"),
                        30: ("MET", "CG,SD,CE")}),

    # ======================================================================================
    # ⭐⭐ BRIDGE — ONE binder spanning BOTH epitopes. Motivated by the site-A and site-B
    #    survivors sitting close together on the receptor surface.
    #
    # MEASURED on the tethered model before building any of this:
    #    closest approach site A <-> site B   11.46 A  (A GLY410 O  <->  B THR19 N)
    #    footprint centroid separation        25.25 A
    #    MAX SPAN across the union            45.59 A   <- what a binder must actually cover
    # ⇒ radius 14 A measured to the UNION of both footprints gives 224 residues / 18 segments,
    #   as much real receptor as the single-site spheres (230/231) and leaving 160 for the
    #   binder inside RFd3's 384-residue training crop (224 + 150 = 374).
    #
    # ⭐⭐ WHY THIS IS WELL MOTIVATED AND NOT JUST GEOMETRICALLY POSSIBLE. Site A is domain III,
    #    site B is domain I. They are adjacent ONLY because the receptor is TETHERED; in the
    #    extended, EGF-bound state those domains separate to clamp EGF. So a bridging binder is
    #    CONFORMATION-SELECTIVE for the tethered (autoinhibited) state, and the organiser has
    #    CONFIRMED the assay target is tethered.
    #    ⇒ It would lock the autoinhibited conformation, and EGF could not displace it by
    #    extending, because extension requires separating the very domains it clamps. That turns
    #    site A's EGF problem into an asset: MEASURED, EGF shares 9 of site A's 11 footprint
    #    residues, so a site-A-only binder competes head-on, while a bridge competes with the
    #    CONFORMATIONAL CHANGE instead of for the surface.
    #
    # ⛔⛔ THE RISK, AND IT IS THE M30 LESSON. A spread hotspot set may be unsatisfiable:
    #    `infer_ori_strategy: hotspots` orients on the hotspot CENTROID, and at site B the
    #    requested hotspot M30 -- 12.2 A from the centroid of the rest -- was reached by
    #    0 of 336 designs, with r(distance from M30, contact rate) = +0.967. A 45 A span is a
    #    far more extreme version of the same geometry. ⇒ TWO ARMS, deliberately:
    #      bridge6  all six hotspots (both anchors + both apolar pairs) -- asks for everything
    #      bridge2  ONLY the two anchors H409 + R29 -- the minimal ask, likelier to reach both
    #    If bridge6 returns nothing and bridge2 returns designs, the spread set is the problem
    #    and we learn it for ~1 shard instead of a campaign.
    # ⚠️ M30 is EXCLUDED from both arms: it is not a key contact residue.
    # ⚠️ Sequon N32's tree passes 10.35 A from the A->B centroid line, and only 1 of its ~15-20
    #    native sugars is modelled -- the bridge path is the most glycan-exposed route we have.
    # ⚠️ The earlier length finding -- extra length buys HELICES, not interface, hence site A
    #    60-80 -- was measured on SINGLE-SITE binders and does NOT apply here: a bridge must physically reach
    #    45.6 A, so length is load-bearing.
    # ⭐ LENGTH = 130-150 -- the reach requirement plus margin. Budget check at the top of the ladder: 224 + 150 = 374 <= 384. ✅
    # ⛔⛔ AND THE 384 IS A **TRAINING CROP, NOT A VRAM LIMIT** -- a bigger GPU does NOT extend it.
    #    A larger card buys memory, not validity: past 384 the model runs outside the size
    #    regime it was trained on. So the target/binder trade below is not purchasable.
    # ⚠️⚠️ THE SHELL IS KNOWINGLY SHALLOWER THAN A BINDER REACHES, AND THIS IS THE REAL COST.
    #    MEASURED on the 12 site-A + 23 site-B survivors: binder atoms reach up to **22.1 A**
    #    from their own footprint (medians 19.8 / 17.9). A 14 A shell is less than that, so a
    #    bridge design's far end faces vacuum where real receptor may sit -- artificial exposure,
    #    the same class of error that made the crop round uninterpretable, though far milder.
    #    A genuine 22 A shell round the 45 A union needs ~320 residues, leaving ~60 for the
    #    binder, which cannot span. ⇒ THE BUDGET CANNOT BUY BOTH.
    #    ⇒ MITIGATION, and it is why this is acceptable: `full_receptor_gate.py` catches exactly
    #    this, for free, after generation. At 22 A (site A) the nearest OMITTED atom was
    #    11.57-19.29 A; at 14 A expect that margin to shrink and a real fraction to FAIL. Budget
    #    for a higher gate-rejection rate rather than for a bigger GPU.
    # ======================================================================================
    "bridge6": dict(radius=14.0, lengths=(130, 135, 140, 145, 150),
                    footprint=[19, 20, 21, 25, 28, 29, 50,
                               380, 382, 384, 408, 409, 410, 411, 412, 417, 438, 465],
                    hotspots={409: ("HIS", "ND1,NE2"),
                              412: ("PHE", "CG,CD1,CD2,CE1,CE2,CZ"),
                              438: ("ILE", "CG2,CD1"),
                              29:  ("ARG", "NE,NH1,NH2"),
                              20:  ("PHE", "CG,CD1,CD2,CE1,CE2,CZ"),
                              25:  ("LEU", "CD1,CD2")}),
    "bridge2": dict(radius=14.0, lengths=(130, 135, 140, 145, 150),
                    footprint=[19, 20, 21, 25, 28, 29, 50,
                               380, 382, 384, 408, 409, 410, 411, 412, 417, 438, 465],
                    hotspots={409: ("HIS", "ND1,NE2"),
                              29:  ("ARG", "NE,NH1,NH2")}),
}
REPS = 3                      # shards per length -- safe ONLY because inference is unseeded
# ⭐ RATE AT 374 TOKENS, MEASURED, and why the shard size did NOT change.
#    bridge L150 = 275.92 s / 8 = 34.49 s/design on the slower of the two machines used;
#    2.84x faster on the other, on the identical workload ⇒ ~12.15 s/design there.
#    ⇒ n_batches=42 (336 designs) = 68.0 min/shard, i.e. PAST `0031`'s original "~1 h" bound.
# ⛔ THE BOUND WAS RAISED RATHER THAN THE SHARD SHRUNK. The bound is arbitrary in itself; what it
#    protects is UNPROTECTED WORK PER FAILURE, i.e. how much generation is lost if a shard dies.
#    At 68 min/shard that exposure is just over an hour per lost shard. See `decisions/0035`,
#    which amends the earlier ~1 h bound to ~2 h.
#    ⇒ the bridge arms keep REPS=3 and n_batches=42, so they match sites A and B exactly:
#    3 shards x 5 lengths x 336 = 5,040 backbones per arm, 1,008 per length.
CHAIN = "A"                   # chain id in MODEL


def load_model():
    res, names = {}, {}
    for ln in open(MODEL):
        if ln.startswith("ATOM") and ln[21] == CHAIN:
            r = int(ln[22:26])
            res.setdefault(r, []).append(ln)
            names.setdefault(r, ln[17:20].strip())
    if not res:
        sys.exit(f"⛔ no chain {CHAIN} ATOM records in {MODEL}")
    return res, names


def sphere(res, names, footprint, radius):
    xyz = {k: np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in v])
           for k, v in res.items()}
    miss = [k for k in footprint if k not in xyz]
    if miss:
        sys.exit(f"⛔ footprint residues absent from the model: {miss}")
    fp = np.vstack([xyz[k] for k in footprint])
    keep = sorted(k for k, v in xyz.items()
                  if np.sqrt(((v[:, None, :] - fp[None, :, :]) ** 2).sum(-1)).min() <= radius)
    segs, s, p = [], keep[0], keep[0]
    for k in keep[1:]:
        if k == p + 1:
            p = k; continue
        segs.append((s, p)); s = p = k
    segs.append((s, p))
    return keep, segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--validate", action="store_true", help="check each spec against the INSTALLED schema")
    a = ap.parse_args()
    res, names = load_model()
    summary = {}

    for site, cfg in SITES.items():
        keep, segs = sphere(res, names, cfg["footprint"], cfg["radius"])
        # ⛔ assert every hotspot by TYPE before writing anything
        for r, (want, atoms) in cfg["hotspots"].items():
            if r not in keep:
                sys.exit(f"⛔ site {site}: hotspot {r} is OUTSIDE the {cfg['radius']} A sphere")
            if names[r] != want:
                sys.exit(f"⛔ site {site}: residue {r} is {names[r]}, expected {want} -- NUMBERING")
        print(f"✅ site {site}: all {len(cfg['hotspots'])} hotspots asserted by type, inside the sphere")

        bdir = os.path.join(a.outdir, f"bundle_site{site}")
        sdir = os.path.join(bdir, "shards")
        if os.path.exists(bdir):
            sys.exit(f"⛔ {bdir} exists -- fresh path per experiment")
        os.makedirs(sdir)
        tgt_name = f"egfr_site{site}_sphere{int(cfg['radius'])}A_WT.pdb"
        with open(os.path.join(bdir, tgt_name), "w") as fh:
            n = 0
            for k in keep:
                for l in res[k]:
                    fh.write(l); n += 1
            fh.write("TER\nEND\n")
        # ⛔ MIRROR THE RUNNER'S LOAD-BEARING COPY: the specs name
        #    the target as a BARE filename and live one level BELOW it, so a copy must sit beside
        #    them or the spec cannot load. Doing it here means local validation tests what the worker
        #    will actually do, instead of passing on a layout the worker does not have.
        shutil.copy(os.path.join(bdir, tgt_name), os.path.join(sdir, tgt_name))
        contig_t = ",".join(f"{CHAIN}{x}-{y}" for x, y in segs)
        hs = {f"{CHAIN}{r}": atoms for r, (_, atoms) in cfg["hotspots"].items()}

        idx, mapping = 0, []
        for L in cfg["lengths"]:
            for rep in range(REPS):
                spec = {f"shard{idx:02d}": {
                    "dialect": 2,
                    "infer_ori_strategy": "hotspots",
                    "input": tgt_name,
                    "contig": f"{L},/0,{contig_t}",
                    "select_hotspots": hs,
                    "is_non_loopy": True}}
                with open(os.path.join(sdir, f"arm_shard{idx:02d}.json"), "w") as fh:
                    json.dump(spec, fh, indent=1)
                mapping.append(dict(shard=f"{idx:02d}", site=site, length=L, rep=rep))
                idx += 1
        with open(os.path.join(bdir, "SHARD_MAP.tsv"), "w") as fh:
            fh.write("shard\tsite\tlength\trep\n")
            for m in mapping:
                fh.write(f"{m['shard']}\t{m['site']}\t{m['length']}\t{m['rep']}\n")
        summary[site] = dict(res=len(keep), segs=len(segs), shards=idx, target=tgt_name,
                             bundle=bdir, tokens_max=len(keep) + max(cfg["lengths"]))
        print(f"   site {site}: {len(keep)} res / {len(segs)} segs, {idx} shards "
              f"({len(cfg['lengths'])} lengths x {REPS}), max tokens {len(keep)+max(cfg['lengths'])}")
        if summary[site]["tokens_max"] > 384:
            print(f"   ⚠️ site {site} max tokens {summary[site]['tokens_max']} EXCEEDS the 384 "
                  f"training crop -- it will run, out of distribution")

        if a.validate:
            specs = [os.path.join(sdir, f) for f in sorted(os.listdir(sdir)) if f.endswith(".json")]
            code = ("import json,sys\n"
                    "from rfd3.inference.input_parsing import DesignInputSpecification as D\n"
                    "for p in sys.argv[1:]:\n"
                    "    d=json.load(open(p))\n"
                    "    for k,v in d.items(): D(**v)\n"
                    "print('✅ %d spec(s) validated against the INSTALLED schema'%(len(sys.argv)-1))\n")
            r = subprocess.run([RFD3_PY, "-c", code] + specs, cwd=sdir,
                               capture_output=True, text=True)
            print("   " + (r.stdout.strip() or r.stderr.strip()[-600:]))
            if r.returncode != 0:
                sys.exit("⛔ SPEC VALIDATION FAILED -- not shipping a bundle that cannot load")

    print("\n⭐ TOTALS")
    tot = sum(v["shards"] for v in summary.values())
    print(f"   {tot} shards x 336 designs (n_batches=42) = {tot*336:,} backbones")
    for site, v in summary.items():
        print(f"   site {site}: {v['bundle']}  target {v['target']}")
    print("\n⛔ ONE .pdb PER BUNDLE (runner dies TARGET_AMBIGUOUS otherwise) -- verify:")
    for site, v in summary.items():
        print(f"   site {site}: {len([f for f in os.listdir(v['bundle']) if f.endswith('.pdb')])} pdb(s)")


if __name__ == "__main__":
    main()
