#!/usr/bin/env python3
"""j496 -- build Boltz-2 BOUND inputs for the C1 v46 JOINT-11 SolubleMPNN sequences.

One FASTA per designed sequence: the designed binder as chain A, plus the WHOLE TNFa TRIMER as
chains B/C/D (the PI's instruction, 2026-10-18: fold them with Boltz-2 against the whole
TNF-alpha trimer).

⛔⛔ CHAIN NAMING MATCHES THE DESIGN FILE so the self-consistency RMSD needs no remap: the v46
   output is binder A + target B,C,D, and so is this. The chain is always named explicitly.

⭐ MSA CONVENTION, copied from the proven j457 v41c inputs rather than invented:
     binder  ->  `>A|protein|empty`                 (de novo; an MSA is meaningless and
                                                     language models do not transfer)
     target  ->  `>B|protein|<tnfa_msa.csv>`        (the real TNFa MSA, 2198 rows)

⭐ THE TARGET SEQUENCE IS NOT RETYPED. It is read from the design .cif's own chain B and ASSERTED
   byte-identical to the 152-aa chain the proven v41c bound inputs used. Two independent sources
   agreeing is what settles identity here (the SEQADV records warn that
   geometry comes from the complex but IDENTITY must not be taken from it alone).

⛔ The DESIGNED chain is pulled from the MPNN FASTA BY LENGTH, and the length is the one recorded
   for that backbone -- the FASTA holds ALL chains colon-separated
   (MPNN's homo-oligomer mode copies the target). A record without `id=` is the INPUT sequence, not a
   design, and is skipped.

⛔ NO CYSTEINE is re-asserted here on the sequence that actually goes to Boltz, not only on MPNN's
   output, because this is the last point before the fold, and cysteines are excluded from every design.

usage: j496_make_boltz_inputs_20261018.py <mpnn-outdir> <design-cif-dir> <out-root>
"""
import glob
import gzip
import os
import sys

from biotite.sequence import ProteinSequence
from biotite.structure import get_residues
from biotite.structure.io.pdbx import CIFFile, get_structure

MPNN_DIR, CIF_DIR, OUT_ROOT = sys.argv[1], sys.argv[2], sys.argv[3]

# PLACEHOLDER PATH -- not a working default; set it for your own tree.
REPO = "/path/to/your/working/tree"
MSA = f"{REPO}/boltz/j457_v41c_fold_20261016/tnfa_msa.csv"
# the 152-aa chain the PROVEN v41c bound inputs used, read back from that file
PROVEN_FASTA = f"{REPO}/boltz/j457_v41c_fold_20261016/inputs/bound/v41c_s11_seq01_bound.fasta"

assert os.path.exists(MSA), f"missing TNFa MSA: {MSA}"


def proven_target():
    with open(PROVEN_FASTA) as fh:
        lines = [l.strip() for l in fh]
    for i, l in enumerate(lines):
        if l.startswith(">D|"):
            return lines[i + 1]
    raise SystemExit("could not read chain D from the proven v41c bound FASTA")


def cif_chain_seqs(path, chains):
    with gzip.open(path, "rt") as fh:
        st = get_structure(CIFFile.read(fh), model=1)
    out = {}
    for ch in chains:
        sub = st[st.chain_id == ch]
        _, names = get_residues(sub)
        out[ch] = "".join(ProteinSequence.convert_letter_3to1(n) for n in names)
    return out


TGT = proven_target()
assert len(TGT) == 152, f"proven target is {len(TGT)} aa, expected 152"

os.makedirs(OUT_ROOT, exist_ok=True)
inp_dir = os.path.join(OUT_ROOT, "inputs")
os.makedirs(inp_dir, exist_ok=True)

fas = sorted(glob.glob(os.path.join(MPNN_DIR, "mpnn", "seqs", "*.fa")))
assert fas, f"no MPNN FASTAs under {MPNN_DIR}/mpnn/seqs"

n_written = 0
manifest = []
for fa in fas:
    stem = os.path.basename(fa).replace(".fa", "")
    cif = os.path.join(CIF_DIR, stem + ".cif.gz")
    assert os.path.exists(cif), f"no design cif for {stem}: {cif}"

    seqs = cif_chain_seqs(cif, ["A", "B", "C", "D"])
    blen = len(seqs["A"])
    # ⛔ the target identity is asserted against the proven file, per chain, not assumed
    for ch in ("B", "C", "D"):
        assert seqs[ch] == TGT, (
            f"{stem} chain {ch} target sequence differs from the proven v41c target -- "
            f"refusing to fold against an unverified target")
    assert blen != 152, f"{stem}: binder is 152 aa, same as a target chain -- ambiguous by length"

    # ---- pull the designed sequences out of the MPNN FASTA
    designs = []
    hdr = ""
    with open(fa) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                hdr = line
                continue
            if not line:
                continue
            if "id=" not in hdr:
                continue                      # the INPUT sequence, not a design
            fields = line.split(":")
            cand = [s for s in fields if len(s) == blen]
            assert len(cand) == 1, (
                f"{stem}: {len(cand)} fields are {blen} aa in {[len(s) for s in fields]} -- "
                f"the designed chain is ambiguous")
            designs.append(cand[0])

    assert designs, f"{stem}: no designed sequences found"
    for i, s in enumerate(designs, start=1):
        assert "C" not in s, f"{stem} seq{i:02d}: CYSTEINE in the sequence going to Boltz"
        tag = f"{stem}_seq{i:02d}"
        out = os.path.join(inp_dir, tag + "_bound.fasta")
        with open(out, "w") as fh:
            fh.write(f">A|protein|empty\n{s}\n")
            for ch in ("B", "C", "D"):
                fh.write(f">{ch}|protein|{MSA}\n{TGT}\n")
        manifest.append((tag, blen, blen + 3 * 152))
        n_written += 1

with open(os.path.join(OUT_ROOT, "manifest.tsv"), "w") as fh:
    fh.write("tag\tbinder_aa\ttotal_res\n")
    for t, b, tot in manifest:
        fh.write(f"{t}\t{b}\t{tot}\n")

tot = {t for _, _, t in manifest}
print(f"[j496] {n_written} bound FASTAs over {len(fas)} backbones -> {inp_dir}")
print(f"[j496] target 152 aa x3 ASSERTED byte-identical to the proven v41c target, on every chain")
print(f"[j496] binder lengths {sorted({b for _, b, _ in manifest})}")
print(f"[j496] total residues per complex {sorted(tot)}  (binder + 456)")
print(f"[j496] ZERO cysteines asserted on every sequence written")
