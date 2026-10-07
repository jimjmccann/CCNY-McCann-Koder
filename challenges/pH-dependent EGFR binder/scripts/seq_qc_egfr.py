#!/usr/bin/env python3
"""Stage 3b — the SEQUENCE-LEVEL half of the design filter, on SolubleMPNN output.

⛔ THE FILTER IS SPLIT IN TWO, because NCPR and the Chan patch need
different inputs. THIS script is the sequence half and runs on fasta alone:
    - guidepost identity held (the real check that --fixed_residues worked)
    - cysteine count == 0            (`--omit_AA C` compliance)
    - N-glycosylation sequon N-X-[ST], X != P   (a free filter; NOT covered by omitting C)
    - NCPR at pH 7.4, against the brief's band
⛔ THE CHAN LARGEST-POSITIVE-PATCH METRIC IS NOT HERE AND CANNOT BE: it needs per-atom SASA on a
folded structure. It lives in `cycle2_F/charge_gate.py` and runs AFTER the fold-back, not here.

NCPR = (#R + #K + f*#H - #D - #E) / N over the BINDER CHAIN ONLY (before the ':').
⚠️ f = 0.5 for His at pH 7.4, matching `charge_gate.py`'s convention so the two agree.
BRIEF BAND (`decisions/0026`): net -0.071/residue, cap -0.094. Gate used here: -0.13 <= NCPR <= -0.02.
⚠️ The band is the TUNING TARGET for `--bias_AA`, whose magnitude is deliberately left unset.
"""
import argparse, glob, os, re, sys
from collections import Counter

POS, NEG = set("RK"), set("DE")
SEQUON = re.compile(r"N[^P][ST]")
LO, HI = -0.13, -0.02


def ncpr(s, f_his=0.5):
    c = Counter(s)
    return (c["R"] + c["K"] + f_his * c["H"] - c["D"] - c["E"]) / len(s)


def read_fastas(d):
    """yield (design, idx, binder_seq). ⛔ SKIPS entry 0: LigandMPNN writes the INPUT sequence
    first (it has no `id=` field), and counting it as a design inflates every statistic."""
    for fa in sorted(glob.glob(os.path.join(d, "*.fa"))):
        name, hdr = os.path.basename(fa)[:-3], None
        for line in open(fa):
            line = line.strip()
            if line.startswith(">"):
                hdr = line; continue
            if not line:
                continue
            if hdr and "id=" not in hdr:
                continue                      # the native/input row
            yield name, hdr, line.split(":")[0]


#: the 20 LigandMPNN designs; anything else in chain A (RFd3's `UNK` placeholder) is SKIPPED by
#: LigandMPNN and therefore absent from its output sequence.
_STD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
           "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
_POSMAP_CACHE = {}


def resid_to_seqpos(pdb_name, guideposts_path):
    """{chain-A PDB residue number -> 1-based position in LigandMPNN's output sequence}.

    Built from the design's own PDB, dropping residues LigandMPNN cannot design, so a leading
    `UNK` shifts nothing. Returns {} if the PDB cannot be found, in which case the caller falls
    back to the (naive) residue number and the check degrades to its old behaviour.
    ⭐ The PDB dir is `pdb/` beside guideposts.json -- the layout stage_egfr_fold_arm.py writes.
    """
    if not guideposts_path:
        return {}
    if pdb_name in _POSMAP_CACHE:
        return _POSMAP_CACHE[pdb_name]
    path = os.path.join(os.path.dirname(os.path.abspath(guideposts_path)), "pdb", pdb_name)
    out, seen = {}, set()
    if os.path.exists(path):
        pos = 0
        for line in open(path):
            if not line.startswith("ATOM") or line[21] != "A":
                continue
            rid = int(line[22:26])
            if rid in seen:
                continue
            seen.add(rid)
            if line[17:20].strip() not in _STD_AA:
                continue                      # UNK / placeholder: LigandMPNN omits it
            pos += 1
            out[rid] = pos
    _POSMAP_CACHE[pdb_name] = out
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("seqdir", help="<outdir>/mpnn/seqs")
    ap.add_argument("--guideposts", default="", help="guideposts.json from mpnn_egfr.py")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()

    gp = {}
    if a.guideposts:
        import json
        gp = json.load(open(a.guideposts))

    rows, n_bad_gp = [], 0
    for name, hdr, s in read_fastas(a.seqdir):
        c = Counter(s)
        # guidepost identity: the ONLY direct evidence that --fixed_residues held
        gok = ""
        key = name + ".pdb"
        if key in gp:
            checks = []
            # ⛔⛔ NEVER INDEX THE SEQUENCE BY PDB RESIDUE NUMBER. LigandMPNN SKIPS residues it
            # cannot design (RFd3 emits a leading `UNK` placeholder with V0-V8 virtual atoms on
            # some backbones), so its output is SHORTER than chain A and every position after the
            # skip is shifted. On 2026-10-01 that made one site-B design
            # (arm_shard30_shard30_L80_13_model_1, UNK at A1, 79-aa sequence vs 80 residues) report
            # "⛔⛔ MPNN MUTATED A GUIDEPOST" for all 100 of its sequences when the Asp was held
            # perfectly -- the check read position 16 (the next residue, L) instead of 15.
            # Site A had 0/244 such backbones, which is why this never fired there.
            posmap = resid_to_seqpos(key, a.guideposts)
            for label, (resid, _d) in gp[key].items():
                want = "D" if label.startswith("p1") else "Y"
                idx = posmap.get(resid, resid) if posmap else resid
                got = s[idx - 1] if 0 < idx <= len(s) else "?"
                checks.append(f"{label}:{resid}{got}{'' if got == want else '!=' + want}")
                if got != want:
                    n_bad_gp += 1
            gok = ";".join(checks)
        rows.append(dict(design=name, n=len(s), ncpr=ncpr(s), cys=c["C"],
                         sequon=len(SEQUON.findall(s)), D=c["D"], E=c["E"], K=c["K"], R=c["R"],
                         guidepost=gok))

    if not rows:
        raise SystemExit(f"⛔ no sequences parsed from {a.seqdir}")
    N = len(rows)
    v = sorted(r["ncpr"] for r in rows)
    inband = sum(LO <= r["ncpr"] <= HI for r in rows)
    print(f"n = {N} designed sequences  ({len({r['design'] for r in rows})} backbones)")
    print(f"  NCPR(7.4)   min {v[0]:+.4f}  p25 {v[N//4]:+.4f}  median {v[N//2]:+.4f}  "
          f"p75 {v[3*N//4]:+.4f}  max {v[-1]:+.4f}")
    print(f"  ⭐ IN BAND [{LO}, {HI}] : {inband}/{N} ({100*inband/N:.1f}%)   "
          f"too POSITIVE {sum(r['ncpr'] > HI for r in rows)}   "
          f"too NEGATIVE {sum(r['ncpr'] < LO for r in rows)}")
    print(f"  cysteine-free       : {sum(r['cys'] == 0 for r in rows)}/{N}"
          + ("   ⛔ --omit_AA C DID NOT HOLD" if any(r["cys"] for r in rows) else "   ✅"))
    print(f"  N-X-[ST] sequon-free: {sum(r['sequon'] == 0 for r in rows)}/{N}"
          f"   (designs carrying >=1 sequon: {sum(r['sequon'] > 0 for r in rows)})")
    if gp:
        print(f"  guidepost identity held: {N - n_bad_gp}/{N}"
              + ("   ⛔⛔ MPNN MUTATED A GUIDEPOST" if n_bad_gp else "   ✅ --fixed_residues held"))
    else:
        print("  ⚠️ no guideposts.json given -- guidepost identity NOT checked")

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print(f"  wrote {a.csv}")


if __name__ == "__main__":
    main()
