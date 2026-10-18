#!/usr/bin/env python3
"""j484 - build the C1 YIELD bundle, v46. 2026-10-17. Local, no spend, NOT uploaded, NOT launched.

⭐⭐ THE PI'S OBJECTIVE, verbatim 2026-10-17: *"We are now looking to produce backbones that we
   can move on with mpnn and boltz folding"* / *"this should be the best possible try given our
   results and progress"* / *"We can do single variable experiments first though if they are
   needed"* / *"I just want to get the backbone generation job up before tomorrow."*
   ⇒ This is a YIELD run. EXPLOIT the measured config; do not explore. ⛔ Not an experiment.

⛔⛔ WHY NO SINGLE-VARIABLE PRE-EXPERIMENT. Every axis is already decided by committed data:
     pins     4, not 7   -- 7 pins gave 0/96 continuous (`J481` §9c); 4 pins 3/112 (`J476`).
                            ⭐ `J484` also measured 4 pins is better on TOPOLOGY: 66.1% of v44
                            clears LCB1's min_cut vs 18.8% of v45, and 4/112 two-lobe vs 26/96.
     core     N=8        -- ⭐⭐ ALL THREE continuous designs the project has produced are in the
                            two N=8 arms: 3/32 vs 0/80 at N=15, -1.393 A at p=0.0004, and the
                            swap costs NO wall time (72.6/72.5 vs 73.0 s/design).
     n_recycle 3         -- axis CLOSED, p=0.379, 0/16 in every A arm (`J476`).
     length   170-220    -- slack x3.02-3.77, far above `J469`'s x1.69 funded bar.
   ⇒ a sweep would spend the deadline window on knowledge this run does not need. The one axis
     still open -- the PIN SET -- is blocked on new geometry and is the PI's call.

⭐⭐ THE N=12 ARM, AND WHY IT IS FREE RATHER THAN AN EXPERIMENT THAT COSTS TIME.
   `J484` measured a TRADE across core size that nobody had seen, because continuity and fold
   topology had never been scored on the same cohort:
       N=8  (B1+B2)  continuous  9.4%   min_cut>=19  50.0%   JOINT  6.25%
       N=15 (A1-A5)  continuous  0.0%   min_cut>=19  72.5%   JOINT  0.00%
   ⇒ N=8 wins CONTINUITY, N=15 wins TOPOLOGY, neither gets both. ⭐ And `J482` independently
     argues a real core at our length is BIGGER than ours (CA Rg 10.85 vs our 9.01).
   ⇒ TWO independent reasons to expect N=12 in the middle. Since a core swap costs ZERO wall time,
     it rides INSIDE the yield run: 256 designs stay at the proven N=8 config (the PI's stated floor)
     and the third arm is bought for nothing.
   ⛔ STATED AS A HYPOTHESIS, NOT A RESULT: two points do not establish monotonicity, and the
     JOINT rate rests on n=2, which is far too few to rank on.

⛔⛔ SIZING IS THE WHOLE QUESTION, and the rate that matters is the JOINT one, 6.25% -- NOT the
   9.4% continuity figure. `J483` measured continuity and fold topology as ORTHOGONAL (overlap 0),
   so MPNN needs a continuous chain and Boltz needs a real fold, and only the conjunction is
   usable; rank on the conjunction, never on a mean of fractions. At 6.25%, 384 designs project to ~24
   usable backbones -- and that projection rests on 2 of 32, so it is a central estimate with a
   very wide interval, which is itself the argument for the larger n.

⭐ HELD BYTE-IDENTICAL TO v44 (the pins are unchanged, so `ori_token` is unchanged too): CONTIG,
  ORI_TOKEN, the runner, the motif-guard patch, `diffusion_batch_size=1` EXPLICIT at launch.
  ⛔ `diffusion_batch_size=1` is load-bearing: unset = 8 means one contig length per batch and the
  effective n becomes BATCHES not designs.

⛔ ARM TOKENS ARE b1/b2/b3 ON PURPOSE. `j476_arm_aggregate.DESIGN_RE` is
  `_shard(\\d\\d)_.*?_([a-z]\\d)_s(\\d\\d)_`, so an arm token MUST be letter+digit immediately
  before `_s<NN>_`. Verified against this slug before the bundle was written -- that exact class of
  bug surfaced only AFTER the GPU time was spent (`J481` §8j).

Local, no spend. ⛔ Writes a NEW directory; never edits a source bundle.
"""
import hashlib
import json
import os
import re
import shutil
import sys

# PLACEHOLDER PATH -- not a working default; set it for your own tree.
REPO = "/path/to/your/working/tree"
BUN = os.path.join(REPO, "status/competition2/bundles")
SRC = os.path.join(BUN, "rfd3_tnf_c1recyc_v44")          # the bundle that already ran clean
DATA = os.path.join(REPO, "status/competition2/data")

RUNNER = "prod_rfd3_tnf_multi.sh"
RUNNER_SHA1 = "623b3c69924b4e0afb245330d5eb2b6f39df10de"

CORE8C = "j470_4pin_core8_compact_20261017.pdb"
CORE8N = "j470_4pin_core8_nearhyd_20261017.pdb"
CORE12 = "j484_4pin_core12_compact_20261017.pdb"

# ⛔ byte-identical to v44/v41a; pins unchanged (D1,E1,E2,E3) so ori_token is unchanged.
CONTIG = ("25-35,E1,39-48,D1,39-47,E3,43-56,E2,20-30,"
          "/0,A6-157,/0,B6-157,/0,C6-157")
ORI_TOKEN = [0.927, 0.807, 40.226]

SLUG = "c1core_v46"
OUT = os.path.join(BUN, "rfd3_tnf_" + SLUG)
SHARDS_PER_ARM = 8
NBATCH = 16                      # with diffusion_batch_size=1 => 16 designs per shard
BINDER_BAND = "170-220"
N_PINS = 4

# (arm, core pdb, n_core, select, n_recycle, role)
ARMS = [
    ("b1", CORE8C, 8, "compact", 3,
     "YIELD: N=8 compact - the config holding 2 of the project's 3 continuous designs"),
    ("b2", CORE8N, 8, "near_hyd", 3,
     "YIELD: N=8 near_hyd - tied with b1 at p=0.217; held because there is NO basis to choose"),
    ("b3", CORE12, 12, "compact", 3,
     "FREE HYPOTHESIS: N=12 - N=8 wins continuity, N=15 wins topology; test the middle at zero "
     "wall-time cost. J482 also argues a real core at 195 aa is BIGGER than ours"),
]

DESIGN_RE = re.compile(r"_shard(\d\d)_.*?_([a-z]\d)_s(\d\d)_")


def digest(p, algo):
    h = hashlib.new(algo)
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 16), b""):
            h.update(c)
    return h.hexdigest()


def main():
    if os.path.exists(OUT):
        sys.exit("⛔ refusing to overwrite an existing bundle: %s" % OUT)

    # ---- 0. the scoring chain must be able to PARSE these arm names, before anything is written
    for arm, _, _, _, _, _ in ARMS:
        probe = "arm_shard00_%s_%s_s00_0_model_0.cif.gz" % (SLUG, arm)
        m = DESIGN_RE.search(probe)
        if not m:
            sys.exit("⛔⛔ DESIGN_RE does not parse %r -- j476_arm_aggregate would aggregate "
                     "NOTHING and the failure would surface only after the GPU time is spent "
                     "(J481 §8j)" % probe)
        if m.group(2) != arm:
            sys.exit("⛔⛔ DESIGN_RE parses %r as arm %r, not %r -- an arm token that resolves to "
                     "the WRONG arm would invert every contrast" % (probe, m.group(2), arm))
    print("✅ DESIGN_RE parses all %d arm tokens correctly against slug %r" % (len(ARMS), SLUG))

    # ---- 1. the runner must be the one that has already run clean
    got = digest(os.path.join(SRC, RUNNER), "sha1")
    if got != RUNNER_SHA1:
        sys.exit("⛔ runner sha1 %s != pinned %s" % (got, RUNNER_SHA1))
    print("✅ gated runner sha1 %s (from the v44 bundle, which ran 14/14 shards clean)"
          % RUNNER_SHA1[:8])

    # ---- 2. every core must exist, declare its own N, and match its own core.json md5
    for pdb, n, sel in ((CORE8C, 8, "compact"), (CORE8N, 8, "near_hyd"), (CORE12, 12, "compact")):
        src = os.path.join(DATA, pdb)
        if not os.path.exists(src):
            sys.exit("⛔ missing core %s" % src)
        meta = json.load(open(src.replace(".pdb", ".core.json")))
        if meta["n_core"] != n:
            sys.exit("⛔ %s declares n_core=%s, expected %d" % (pdb, meta["n_core"], n))
        if meta["select"] != sel:
            sys.exit("⛔ %s declares select=%r, expected %r" % (pdb, meta["select"], sel))
        if digest(src, "md5") != meta["out_pdb_md5"]:
            sys.exit("⛔ %s md5 does not match its own core.json" % pdb)
        cl = meta["clearances"]
        tnf = cl["chainF_CA_CB  -> TNFa A/B/C heavy"]
        pin = cl["chainF_CA_CB  -> pin heavy"]
        # ⛔ gate on pin HEAVY atoms, never CA -- the pins are whole residues and every atom is held
        if pin < 3.5:
            sys.exit("⛔ %s core->pin heavy %.3f A is inside the 3.5 A gate" % (pdb, pin))
        print("✅ gated %-44s N=%2d %-9s md5 %s  ->TNFa %6.2f  ->pin %5.2f  Rg %4.2f"
              % (pdb, n, sel, meta["out_pdb_md5"][:8], tnf, pin, meta["core_geometry"]["rg"]))

    # ---- 3. write the bundle
    os.makedirs(os.path.join(OUT, "shards"))
    shutil.copy2(os.path.join(SRC, RUNNER), os.path.join(OUT, RUNNER))
    shutil.copy2(os.path.join(SRC, "patch_held_motif_guard.py"),
                 os.path.join(OUT, "patch_held_motif_guard.py"))
    for pdb in (CORE8C, CORE8N, CORE12):
        shutil.copy2(os.path.join(DATA, pdb), os.path.join(OUT, pdb))

    n_shards = len(ARMS) * SHARDS_PER_ARM
    rows = []
    for sid_i in range(n_shards):
        arm, pdb, ncore, sel, nrec, role = ARMS[sid_i % len(ARMS)]   # ⭐ INTERLEAVED
        sid = "%02d" % sid_i
        body = {"dialect": 2, "input": pdb, "contig": CONTIG,
                "unindex": ",".join("F%d" % i for i in range(1, ncore + 1)),
                "select_fixed_atoms": {"F%d" % i: "CA,CB" for i in range(1, ncore + 1)},
                "select_unfixed_sequence": "F1-%d" % ncore,
                "ori_token": list(ORI_TOKEN), "is_non_loopy": True}
        key = "%s_%s_s%s" % (SLUG, arm, sid)
        with open(os.path.join(OUT, "shards", "arm_shard%s.json" % sid), "w") as fo:
            json.dump({key: body}, fo, indent=1)
        with open(os.path.join(OUT, "shards", "arm_shard%s.flags" % sid), "w") as fo:
            fo.write("inference_sampler.n_recycle=%d\n" % nrec)
        rows.append((sid, arm, pdb, ncore, sel, nrec, role))

    note = ("C1 YIELD RUN v46 (PI, 10-17: backbones to move on with MPNN and Boltz, best "
            "possible try, job up before tomorrow) ; NOT an experiment - every axis is already "
            "decided: 4 pins (7 gave 0/96), N=8 core (ALL 3 continuous designs are in the N=8 "
            "arms, 3/32 vs 0/80, p=0.0004), n_recycle=3 (axis closed p=0.379), 170-220 aa "
            "(slack x3.02-3.77) ; b1/b2 are the two N=8 cores, TIED at p=0.217 so both are held "
            "and 96 draws each raises that comparison from n=16 to n=96 FREE ; b3 N=12 is a FREE "
            "hypothesis - J484 measured N=8 wins continuity 9.4% vs 0% while N=15 wins topology "
            "72.5% vs 50%, so neither gets both, and J482 argues a real core at 195 aa is BIGGER "
            "than ours ; sized on the JOINT rate 6.25% (continuous AND min_cut>=19) NOT the 9.4% "
            "continuity rate, because J483 measured the two axes ORTHOGONAL with overlap 0 ; "
            "shards INTERLEAVED (shard NN = arm NN mod 3) so a short run stays BALANCED ; CONTIG "
            "and ori_token byte-identical to v44 ; diffusion_batch_size=1 EXPLICIT at launch")
    with open(os.path.join(OUT, "ARMS.tsv"), "w") as fo:
        fo.write("shard\tslug\tarm\ttest\tn_recycle\tn_core\tcore_select\tinput_pdb\tbinder_aa\t"
                 "pins\trole\tnote\n")
        for sid, arm, pdb, ncore, sel, nrec, role in rows:
            test = "T_YIELD" if arm in ("b1", "b2") else "T_CORE12"
            fo.write("%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%s\n"
                     % (sid, SLUG, arm.upper(), test, nrec, ncore, sel, pdb,
                        BINDER_BAND, N_PINS, role, note))

    # ---- 4. verify what was WRITTEN, not what was intended
    seen = {}
    for sid_i in range(n_shards):
        sid = "%02d" % sid_i
        d = json.load(open(os.path.join(OUT, "shards", "arm_shard%s.json" % sid)))
        (k,) = list(d)
        b = d[k]
        arm, pdb, ncore, sel, nrec, _ = ARMS[sid_i % len(ARMS)]
        if arm not in k or not k.endswith("_s%s" % sid):
            sys.exit("⛔ shard%s key %r does not name arm %s" % (sid, k, arm))
        if b["contig"] != CONTIG or b["ori_token"] != ORI_TOKEN:
            sys.exit("⛔ shard%s drifted on a HELD field" % sid)
        if b["input"] != pdb or not os.path.exists(os.path.join(OUT, pdb)):
            sys.exit("⛔ shard%s names input %r not in the bundle" % (sid, b["input"]))
        if len(b["select_fixed_atoms"]) != ncore or b["select_unfixed_sequence"] != "F1-%d" % ncore:
            sys.exit("⛔ shard%s core-atom bookkeeping is not N=%d" % (sid, ncore))
        f = open(os.path.join(OUT, "shards", "arm_shard%s.flags" % sid)).read().strip()
        if f != "inference_sampler.n_recycle=%d" % nrec:
            sys.exit("⛔ shard%s flags %r != n_recycle=%d" % (sid, f, nrec))
        seen.setdefault(arm, []).append(sid)
    if len(seen) != len(ARMS) or any(len(v) != SHARDS_PER_ARM for v in seen.values()):
        sys.exit("⛔ allocation is not %d x %d: %s"
                 % (len(ARMS), SHARDS_PER_ARM, {k: len(v) for k, v in seen.items()}))

    # the three arms must differ ONLY in the core (input + the F-atom bookkeeping)
    base = {}
    for arm in seen:
        s = seen[arm][0]
        b = list(json.load(open(os.path.join(OUT, "shards/arm_shard%s.json" % s))).values())[0]
        base[arm] = {k: v for k, v in b.items()
                     if k not in ("input", "unindex", "select_fixed_atoms",
                                  "select_unfixed_sequence")}
    if len({json.dumps(v, sort_keys=True) for v in base.values()}) != 1:
        sys.exit("⛔ the arms differ in something OTHER than the core -- not single-variable")
    print("✅ all %d arms are identical outside the core (input + F-atom bookkeeping)" % len(ARMS))

    # cross-check ARMS.tsv against the design keys the runner will emit
    import csv
    for row in csv.DictReader(open(os.path.join(OUT, "ARMS.tsv")), delimiter="\t"):
        sid, arm = row["shard"], row["arm"].lower()
        probe = "arm_shard%s_%s_%s_s%s_0_model_0.cif.gz" % (sid, SLUG, arm, sid)
        m = DESIGN_RE.search(probe)
        if not m or m.group(1) != sid or m.group(2) != arm or m.group(3) != sid:
            sys.exit("⛔ ARMS.tsv row shard=%s arm=%s will not round-trip through DESIGN_RE"
                     % (sid, arm))
    print("✅ every ARMS.tsv row round-trips through the scorer's DESIGN_RE")

    nf = sum(len(fs) for _, _, fs in os.walk(OUT))
    print("\n⇒ %s" % OUT)
    print("   %d files | %d shards x %d = %d designs | %d arms x %d designs"
          % (nf, n_shards, NBATCH, n_shards * NBATCH, len(ARMS), SHARDS_PER_ARM * NBATCH))
    print("   arm -> shards: %s" % {k: ",".join(v) for k, v in sorted(seen.items())})
    print("   LAUNCH KNOBS: RFD3_SHARDS=%d RFD3_NBATCH=%d  (diffusion_batch_size=1 => NBATCH IS "
          "designs/shard; ⛔ the launcher's own 'nb*8' echo is EGFR arithmetic and is WRONG here)"
          % (n_shards, NBATCH))
    print("   per-shard generation: %d x 72.6 s = %.2f h  (inside the shard wall-time ceiling)"
          % (NBATCH, NBATCH * 72.6 / 3600.0))
    print("   projected wall time: %.2f h  (72.6 s/design + ~85 s/shard overhead)"
          % ((n_shards * NBATCH * 72.6 + n_shards * 85) / 3600.0,))
    print("⛔ NOT uploaded. ⛔ NOT launched.")


if __name__ == "__main__":
    main()
