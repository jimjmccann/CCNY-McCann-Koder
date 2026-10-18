#!/usr/bin/env python3
"""j168_mpnn_tnf.py -- SolubleMPNN sequence design on the apex TNF binder backbones.

A TNF-specific driver rather than a reuse of challenge 1's deposited
`../../pH-dependent EGFR binder/scripts/mpnn_egfr.py`: that script
requires a `--wt2out` WT->output residue map and asserts an EGFR anchor type per design, neither
of which exists for this target. The MPNN invocation itself is copied from it verbatim (ADR 0030
flags), so the settled decisions carry over.

SETTLED FLAGS KEPT FROM ADR 0030 / mpnn_egfr.py:
  model_type=soluble_mpnn, checkpoint solublempnn_v_48_020.pt
  --omit_AA C   ⛔⛔ NON-NEGOTIABLE: cysteines break the competition's expression system
                (the PI). Free thiols are an oxidation and aggregation liability.
  P PERMITTED   (do not omit; proline caps helices and turns)
  --chains_to_design A   the target is CONTEXT, never designed.

⛔⛔ THIS IS THE **C1** (ONE-BINDER) DRIVER. It REFUSES a multi-copy file — see guard 1.

⛔ WHAT THIS SCRIPT GUARDS, because getting it wrong is silent:
  1. **CHAIN IDENTITY, PER DESIGN, EVERY CHAIN.** — chain selection
     has silently changed a conclusion four times here. ⚠️ FIXED 2026-10-12 (`J242` §1a): the old
     check read only chains A and B, so a C3 file (binder A,C,E @222; target B,D,F @152) matched
     BOTH and sailed through, after which `--chains_to_design A` designed ONE copy of three while
     holding the other two as fixed context carrying RFd3's own mostly-Ala sequence ⇒ a sequence
     optimised against poly-Ala neighbours, no self-interface, **no error**. Roles are now read
     from the file's own `is_motif_atom_with_fixed_seq` (`J248` §4's discriminator, not length),
     every chain is inspected, and >1 binder chain REFUSES and points at the tied C3 driver.
  2. **THE PINNED MOTIF IS PASSED AS `--fixed_residues_multi`.** ⚠️ ADDED 2026-10-12 (`J242` §1b):
     the driver passed NO fixed residues, so MPNN would overwrite HIS/ARG/ARG/ARG/GLU — the motif
     required on every design — silently. ⛔ Positions are read PER DESIGN from the file, never
     hardcoded: contig ranges are re-sampled per parser call and the two C3 batches landed on
     different indices (61/89/97/124/154 vs 65/91/99/126/155). A zero count REFUSES.
  3. **NO CYSTEINE IN THE DESIGNED CHAIN.** Checked after the run on the DESIGNED chain's FASTA
     field, not field 1 — the FASTA holds ALL chains colon-separated
     (MPNN's homo-oligomer mode copies the target).
Local, no spend, seconds per backbone on a single workstation GPU.

usage: j168_mpnn_tnf.py <indir-with-cif.gz> <outdir> [--seqs 8] [--run]
"""
import argparse
import glob
import gzip
import json
import os
import subprocess
import sys

LMPNN = os.environ.get("LIGANDMPNN_DIR", "/path/to/LigandMPNN")  # PLACEHOLDER default
PY = os.environ.get("MPNN_PY", "/path/to/your/python")  # PLACEHOLDER default


def load(path):
    from biotite.structure.io.pdbx import CIFFile, get_structure
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        cf = CIFFile.read(fh)
    st = get_structure(cf, model=1)
    return st


TRUTHY = ("True", "true", "TRUE", "1", "Y", "yes")


def roles_from_cif(path):
    """Classify chains as binder/target and return guidepost residues, from the file itself.

    ⛔ THE DISCRIMINATOR IS NOT LENGTH (`J248` §4): it is the per-chain fraction of residues
    carrying `is_motif_atom_with_fixed_seq` -- the TARGET is ~100% (fixed by definition), a
    designed BINDER ~2% (its pinned guideposts only). A 152-aa target and a 152-aa binder copy
    are indistinguishable by length. Same logic as `j261_detect_layout.py`.
    ⛔ ALLOWLIST, never a denylist: the column holds the STRING "False", and a denylist of falsy
    spellings read it as TRUE and called every chain TARGET (`J262`).
    """
    from collections import defaultdict
    cols, rows, in_loop = [], [], False
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1])
                in_loop = True
                continue
            if in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
                f = s.split()
                if len(f) == len(cols):
                    rows.append(dict(zip(cols, f)))
            elif in_loop and rows and (s.startswith("#") or s.startswith("loop_")):
                break
    fix_key = next((c for c in cols if "fixed_seq" in c), None)
    if not rows or not fix_key:
        raise SystemExit(
            f"⛔⛔ {os.path.basename(path)}: no atom_site rows parsed, or NO *fixed_seq* column "
            f"({len(rows)} rows, column present={bool(fix_key)}). Refusing: without that column "
            f"the binder/target call would fall back to LENGTH, which `J248` §4 measured is not a "
            f"sound discriminator, and the guidepost set would be unknown.")
    ch_key = "auth_asym_id" if "auth_asym_id" in cols else "label_asym_id"
    seq_key = "auth_seq_id" if "auth_seq_id" in cols else "label_seq_id"
    comp_key = "auth_comp_id" if "auth_comp_id" in cols else "label_comp_id"

    res, fixed, names = defaultdict(set), defaultdict(set), {}
    for r in rows:
        ch, ri = r[ch_key], r[seq_key]
        res[ch].add(ri)
        names[(ch, ri)] = r[comp_key]
        if r.get(fix_key, "") in TRUTHY:
            fixed[ch].add(ri)
    binder, target, amb = [], [], []
    for c in sorted(res):
        pct = 100.0 * len(fixed.get(c, ())) / max(len(res[c]), 1)
        (target if pct > 90 else binder if pct < 20 else amb).append(c)
    gp = {c: sorted(((ri, names[(c, ri)]) for ri in fixed.get(c, ())), key=lambda t: int(t[0]))
          for c in binder}
    return binder, target, amb, gp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", help="directory of RFd3 *.cif.gz")
    ap.add_argument("outdir")
    ap.add_argument("--seqs", type=int, default=8, help="sequences per backbone (batch_size)")
    ap.add_argument("--temperature", type=float, default=0.1)
    ap.add_argument("--omit-aa", default="C", help="⛔ keep C unless the PI has agreed in writing")
    ap.add_argument("--bias-raw", default="", help="raw --bias_AA string, e.g. 'A:-1.0'")
    ap.add_argument("--binder-len", type=int, default=130)
    ap.add_argument("--target-len", type=int, default=77)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--run", action="store_true", help="actually run; otherwise print the command")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    pdbdir = os.path.join(a.outdir, "pdb")
    os.makedirs(pdbdir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(a.indir, "*.cif.gz")))
    if a.limit:
        files = files[:a.limit]
    if not files:
        raise SystemExit(f"⛔ no *.cif.gz in {a.indir}")

    from biotite.structure.io.pdb import PDBFile
    import numpy as np

    paths = {}
    BLEN = {}          # design stem -> binder length, read PER DESIGN (see the AUTO guard below)
    fixed_multi = {}
    for f in files:
        st = load(f)
        # ⛔⛔ VALIDATE **EVERY** CHAIN, not A and B. The old check read only chains A and B and
        #    on a C3 file (binder A,C,E @222; target B,D,F @152) BOTH matched, so it sailed
        #    through and then ran `--chains_to_design A`, designing ONE copy of three against two
        #    poly-Ala neighbours -- no error, plausible-looking output (`J242` §1a).
        chains = [c for c in dict.fromkeys(st.chain_id.tolist())]
        binder, target, amb, gp = roles_from_cif(f)
        nres = {c: len(np.unique(st.res_id[st.chain_id == c])) for c in chains}
        if amb:
            raise SystemExit(
                f"⛔⛔ {os.path.basename(f)}: chains {amb} are AMBIGUOUS (fixed-seq fraction "
                f"20-90%). Refusing -- resolve the layout before designing.")
        if len(binder) != 1:
            raise SystemExit(
                f"⛔⛔ {os.path.basename(f)}: detected {len(binder)} BINDER chains {binder} "
                f"(target {target}; residues {nres}).\n"
                f"   THIS DRIVER IS THE C1 (one-binder) PATH and would design only {binder[:1]}, "
                f"holding the others as FIXED CONTEXT carrying RFd3's own mostly-Ala sequence -- "
                f"a sequence optimised against poly-Ala neighbours, with no self-interface and no "
                f"error (`J242` §1a, the most dangerous hole found).\n"
                f"   ⇒ Use the TIED C3 driver for a multi-copy file. Refusing.")
        bch, tch = binder[0], target
        # ⭐ `--binder-len 0` = AUTO: take the binder length PER DESIGN from the file instead of
        #    asserting one global constant. Needed from 2026-10-18: the C1 `v46` JOINT cohort is
        #    187-204 aa across 11 designs, so no single value can gate it, and the old default of
        #    130 (an earlier campaign's length) REFUSED every one of them. ⛔ AUTO still refuses a
        #    binder whose length collides with a target chain's -- see the FASTA guard below, which
        #    identifies the designed chain BY LENGTH and would otherwise pick the wrong field.
        if a.binder_len and nres[bch] != a.binder_len:
            raise SystemExit(
                f"⛔ {os.path.basename(f)}: binder chain {bch} has {nres[bch]} residues "
                f"(expected {a.binder_len}). Refusing -- the designed chain may not be the binder.")
        if any(nres[c] == nres[bch] for c in tch):
            raise SystemExit(
                f"⛔⛔ {os.path.basename(f)}: binder chain {bch} is {nres[bch]} aa and a TARGET "
                f"chain has the same length ({ {c: nres[c] for c in tch} }). The FASTA post-check "
                f"identifies the designed chain BY LENGTH and would be ambiguous. Refusing "
                f"(MPNN's homo-oligomer mode copies the target).")
        BLEN[os.path.basename(f).replace(".cif.gz", "")] = nres[bch]
        ntgt = sum(nres[c] for c in tch)
        if ntgt != a.target_len:
            raise SystemExit(
                f"⛔ {os.path.basename(f)}: target chains {tch} total {ntgt} residues "
                f"(expected {a.target_len}). Refusing.")
        if bch != "A":
            raise SystemExit(
                f"⛔ {os.path.basename(f)}: the binder is chain {bch}, not A, but the invocation "
                f"below hardcodes `--chains_to_design A`. Refusing rather than designing the "
                f"wrong chain -- the chain must always be named explicitly.")

        out = os.path.join(pdbdir, os.path.basename(f).replace(".cif.gz", ".pdb"))
        pf = PDBFile()
        pf.set_structure(st)
        pf.write(out)
        paths[out] = ""

        # ⛔⛔ HOLE (b), `J242` §1b: the old driver passed NO --fixed_residues, so MPNN redesigned
        #    the entire pinned motif -- the Arg/His cluster required on EVERY design -- with no
        #    error. `mpnn_egfr.py` has this guard and says so in its header; this driver was
        #    written fresh and LOST it. ⛔ Positions are read PER DESIGN from the file's own
        #    `is_motif_atom_with_fixed_seq`, never hardcoded: contig ranges are re-sampled per
        #    parser call, and the two C3 batches really did land on different indices (61/89/97/
        #    124/154 vs 65/91/99/126/155), so a hardcoded list is wrong for half a cohort.
        #    Key format is chain+resnum, matching `run.py:212` (`chain_letters[i] + R_idx + icode`).
        gplist = [f"{bch}{ri}" for ri, _ in gp.get(bch, [])]
        if not gplist:
            raise SystemExit(
                f"⛔⛔ {os.path.basename(f)}: ZERO fixed-seq (guidepost) residues found on binder "
                f"chain {bch}. Refusing: with no --fixed_residues MPNN would overwrite the motif "
                f"silently, and a zero count is exactly the measurement that must not be trusted "
                f"(`J262`: a denylist read the string \"False\" as true and reported 0 of 15).")
        fixed_multi[out] = " ".join(gplist)  # ⛔ WHITESPACE-separated (run.py:131), unlike
        #                                       --chains_to_design which is COMMA-separated (:281)

    pj = os.path.join(a.outdir, "pdb_path_multi.json")
    json.dump(paths, open(pj, "w"), indent=1)
    fj = os.path.join(a.outdir, "fixed_residues_multi.json")
    json.dump(fixed_multi, open(fj, "w"), indent=1)
    ngp = sorted({len(v.split()) for v in fixed_multi.values()})
    print(f"⭐ {len(paths)} backbones converted; ONE binder chain + target validated in ALL of "
          f"them, every chain inspected.")
    print(f"⭐ guidepost residues PINNED per design (counts {ngp}); e.g. "
          f"{next(iter(fixed_multi.values()))}")

    cmd = [PY, os.path.join(LMPNN, "run.py"),
           "--model_type", "soluble_mpnn",
           "--checkpoint_soluble_mpnn",
           os.path.join(LMPNN, "model_params", "solublempnn_v_48_020.pt"),
           "--pdb_path_multi", pj,
           # ⛔⛔ WITHOUT THIS, MPNN REDESIGNS THE WHOLE Arg/His MOTIF, silently (`J242` §1b).
           #    `_multi` because positions differ per design; the file maps pdb path -> a
           #    WHITESPACE-separated residue list (`run.py:126-129`). The mask composes as
           #    `chain_mask * fixed_positions` (`run.py:299`), i.e. designable = in
           #    chains_to_design AND not fixed -- a plain multiplication, no conflict.
           "--fixed_residues_multi", fj,
           "--chains_to_design", "A",
           "--omit_AA", a.omit_aa,
           "--out_folder", os.path.join(a.outdir, "mpnn"),
           "--batch_size", str(a.seqs),
           "--number_of_batches", "1",
           "--temperature", str(a.temperature),
           "--save_stats", "1",
           "--seed", "37"]
    if a.bias_raw:
        cmd += ["--bias_AA", a.bias_raw]

    print(" ".join(cmd))
    if not a.run:
        print("(dry run -- pass --run to execute)")
        return 0
    rc = subprocess.call(cmd)
    if rc != 0:
        raise SystemExit(f"⛔ LigandMPNN exited {rc}")

    # ---- post-check: NO CYSTEINE in the DESIGNED chain ----
    fas = sorted(glob.glob(os.path.join(a.outdir, "mpnn", "seqs", "*.fa")))
    nseq, ncys, lens = 0, 0, set()
    nat_cys = 0
    for p in fas:
        hdr = ""
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    hdr = line
                    continue
                if not line:
                    continue
                # ⛔⛔ THE FIRST RECORD OF EVERY LigandMPNN FASTA IS THE **INPUT** SEQUENCE, not a
                #    design: only the designed records carry `id=` in the header. Counting the
                #    native record made this gate report "CYSTEINES PRESENT" on a run where all
                #    64 designs were clean (measured 2026-10-10) -- and RFd3's own emitted partial
                #    sequence really does contain C, so the native count is not zero and the false
                #    alarm looks exactly like a true one.
                # ⭐ the expected length is resolved PER DESIGN (BLEN, from the file itself), not
                #    from one global constant -- the v46 JOINT cohort spans 187-204 aa.
                blen = BLEN.get(os.path.basename(p).replace(".fa", ""), a.binder_len)
                if not blen:
                    raise SystemExit(f"⛔ {os.path.basename(p)}: no recorded binder length for "
                                     f"this design -- cannot identify the designed chain.")
                if "id=" not in hdr:
                    f = [s for s in line.split(":") if len(s) == blen]
                    if f:
                        nat_cys += f[0].count("C")
                    continue
                # ⛔ the FASTA holds ALL chains colon-separated; the DESIGNED chain is the one
                #    whose length matches the binder. Never read field 1 blindly.
                fields = line.split(":")
                designed = [s for s in fields if len(s) == blen]
                if not designed:
                    raise SystemExit(f"⛔ {os.path.basename(p)}: no {blen}-aa field in "
                                     f"{[len(s) for s in fields]} -- cannot identify the "
                                     f"designed chain.")
                if len(designed) > 1:
                    raise SystemExit(f"⛔⛔ {os.path.basename(p)}: {len(designed)} fields are "
                                     f"{blen} aa -- the designed chain is AMBIGUOUS by length. "
                                     f"Refusing rather than scoring the wrong chain.")
                s = designed[0]
                nseq += 1
                lens.add(len(s))
                ncys += s.count("C")
    print(f"⭐ {nseq} DESIGNED sequences over {len(fas)} backbones; designed-chain lengths "
          f"{sorted(lens)}; CYSTEINES in designs: {ncys}")
    print(f"   (for reference, RFd3's own input sequences carry {nat_cys} Cys across those "
          f"{len(fas)} backbones -- excluded, they are not designs)")
    if ncys:
        raise SystemExit("⛔⛔ CYSTEINES PRESENT -- these cannot be submitted.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
