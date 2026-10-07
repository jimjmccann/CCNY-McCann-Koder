#!/usr/bin/env python3
"""Build the Scheme A shard set for the RFd3 EGFR production run.

⭐⭐ SCHEME A: 31 shards x 333 designs, ONE FIXED LENGTH EACH, 70-100 aa.
   70..100 inclusive is exactly 31 values, so the mapping is one length per shard with no
   choice left to make -- shard NN gets length 70+NN.

⭐ WHY FIXED LENGTH PER SHARD rather than a range per shard: it makes novelty-vs-length FREE.
   Every design in an output archive has a known length, so the novelty rate can be read
   against length with no extra bookkeeping. `0c` predicts longer = more novel; this is what
   lets us check that on landing instead of assuming it.

⛔ THE LENGTH LIVES IN `contig`, NOT IN A `length` FIELD. The base arm.json carries
   "contig": "70,/0,B370-475" -- the leading integer is the binder chain length. Writing a
   `length` key alongside it would be silently ignored by the one that matters, so this script
   rewrites the contig and refuses if it cannot find that leading integer.

⛔ `DesignInputSpecification` sets extra="forbid", so a stray key fails at load. This script
   therefore COPIES a base arm that already validates and changes exactly one thing.

⛔ NOT set here (CLI-only, they are NOT spec fields): n_batches, inference_sampler.*,
   use_classifier_free_guidance, cfg_scale. Production is CFG OFF -- both expensive sampler
   knobs were REFUTED (CFG 2.76x null, n_recycle 2.55x null-to-worse).

Usage:
  python3 make_rfd3_shards.py --base ../rfd3_arms/armR1_two_island.json --out <bundle>/shards
"""
import argparse, json, os, re, shutil, subprocess, sys

RFD3_PY = "/PATH/TO/rfd3-venv/bin/python"


def rewrite_contig(contig: str, length: int) -> str:
    """Replace the leading binder length in a contig, e.g. '70,/0,B370-475' -> '85,/0,B370-475'."""
    m = re.match(r"^(\d+)(,.*)$", contig)
    if not m:
        raise SystemExit(f"⛔ contig {contig!r} does not start with an integer length -- "
                         "refusing to guess where the length lives")
    return f"{length}{m.group(2)}"


def validate(paths):
    """Validate each arm.json against the INSTALLED schema, not against our idea of it.

    ⛔ A JSON that parses is not a spec that loads: extra='forbid' means a typo'd key is only
    caught here. Skipped (loudly) if the venv is absent, never silently.
    """
    if not os.path.exists(RFD3_PY):
        print(f"⚠️  rfd3 venv not found at {RFD3_PY} -- SCHEMA VALIDATION SKIPPED, not passed")
        return None
    code = (
        "import json,sys\n"
        "from rfd3.inference.input_parsing import DesignInputSpecification as D\n"
        "bad=0\n"
        "for p in sys.argv[1:]:\n"
        "    d=json.load(open(p))\n"
        "    for name,spec in d.items():\n"
        "        try: D(**spec)\n"
        "        except Exception as e: print('FAIL',p,name,e); bad+=1\n"
        "print('VALIDATED',len(sys.argv)-1,'file(s), failures',bad)\n"
        "sys.exit(1 if bad else 0)\n")
    # ⛔ CWD MUST BE THE BUNDLE DIR. DesignInputSpecification resolves `input` against the
    #    PROCESS CWD, not against the arm.json naming it (MEASURED -- the
    #    same trap that would have killed every worker with FileNotFoundError). Validating from
    #    anywhere else reports 31 spurious failures and hides any real one.
    cwd = os.path.dirname(os.path.abspath(paths[0]))          # <bundle>/shards
    cwd = os.path.dirname(cwd)                                # <bundle>, where the pdb is staged
    r = subprocess.run([RFD3_PY, "-c", code, *paths], capture_output=True, text=True, cwd=cwd)
    print(r.stdout.strip() or r.stderr.strip()[-2000:])
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="arm.json to clone (the conditioning)")
    ap.add_argument("--out", required=True, help="output shards/ directory")
    ap.add_argument("--lmin", type=int, default=70)
    ap.add_argument("--lmax", type=int, default=100)
    ap.add_argument("--target", default=None,
                    help="target pdb to copy next to the shards (default: the base's `input`)")
    a = ap.parse_args()

    base = json.load(open(a.base))
    if len(base) != 1:
        raise SystemExit(f"⛔ {a.base} holds {len(base)} arms; expected exactly 1")
    base_name, spec = next(iter(base.items()))
    if "contig" not in spec:
        raise SystemExit("⛔ base arm has no `contig` -- cannot set a length")

    lengths = list(range(a.lmin, a.lmax + 1))
    os.makedirs(a.out, exist_ok=True)
    written = []
    for i, L in enumerate(lengths):
        s = dict(spec)
        s["contig"] = rewrite_contig(spec["contig"], L)
        name = f"shard{i:02d}_L{L}"
        p = os.path.join(a.out, f"arm_shard{i:02d}.json")
        json.dump({name: s}, open(p, "w"), indent=1)
        written.append(p)

    # ⭐ Stage the target pdb beside the shards: the worker resolves `input` against its CWD.
    tgt = a.target or spec.get("input")
    if tgt:
        src = os.path.join(os.path.dirname(os.path.abspath(a.base)), tgt)
        dst = os.path.join(os.path.dirname(os.path.abspath(a.out)), os.path.basename(tgt))
        if os.path.exists(src):
            shutil.copy2(src, dst); print(f"staged target -> {dst}")
        else:
            print(f"⚠️  target {src} NOT FOUND -- stage it by hand before uploading the bundle")

    print(f"\nwrote {len(written)} shards into {a.out}")
    print(f"  base arm      : {base_name}  ({a.base})")
    print(f"  lengths       : {a.lmin}..{a.lmax}  ({len(lengths)} shards, one length each)")
    print(f"  contig example: {json.load(open(written[0]))[f'shard00_L{a.lmin}']['contig']!r}"
          f" .. {json.load(open(written[-1]))[f'shard{len(lengths)-1:02d}_L{a.lmax}']['contig']!r}")
    ok = validate(written)
    if ok is False:
        raise SystemExit("⛔ schema validation FAILED -- do not upload this bundle")
    print("⭐ shards ready" + ("" if ok else " (UNVALIDATED -- venv missing)"))


if __name__ == "__main__":
    main()
