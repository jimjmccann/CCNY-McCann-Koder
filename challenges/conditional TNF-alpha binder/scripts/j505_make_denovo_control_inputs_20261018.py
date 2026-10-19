#!/usr/bin/env python3
"""j505 -- build the DE NOVO positive control Boltz-2 input: LCB1 + its own target (SARS-CoV-2 RBD).

WHY. `J502` showed Boltz-2 places the natural receptor TNFR2 on the real TNFα epitope under exactly
our conditions (recall 0.781, 6.25× over chance), which removed "the instrument is blind to this
target" as an explanation for `J497`'s 0/110. ⛔⛔ BUT IT LEFT ONE ALTERNATIVE STANDING, and `j501`'s
header pre-registered it: `3ALQ` was published in 2010 and is very likely in Boltz-2's training set,
so the pass may be **recall of a memorised complex** rather than an ability to place a **novel**
chain. TNFR2 is also a *natural* sequence whose family the model has seen, while our binders are
de novo and have no evolutionary context at all.

⇒ THIS ARM TESTS THE DE NOVO REGIME SPECIFICALLY. **LCB1** is a computationally designed
miniprotein (Cao et al. 2020) with **measured picomolar-range affinity** for the SARS-CoV-2 RBD and
a solved complex, `7JZU`. It has no natural homologs, so it is handicapped the way our designs are
handicapped — a designed sequence, empty MSA, no family. If Boltz places LCB1 on the true RBD
epitope, then "the method cannot dock a de novo designed chain" is refuted.

⭐⭐ AND IT CARRIES A GROUND TRUTH `J502` COULD NOT HAVE: the complex is solved, so the predicted
binder pose can be compared to the CRYSTAL pose directly (superpose the predicted target onto the
crystal target, then measure binder RMSD). For our own designs no crystal exists, which is why
`J497` had to measure displacement from the *designed* site instead. ⇒ scored by `j506`.

## ⛔⛔ THE CAVEAT THAT SURVIVES EVEN THIS ARM -- state it, do not let a PASS over-claim
`7JZU` was published in **2020** and may ALSO be in Boltz-2's training set. ⇒ a pass here shows the
model can place a **de novo designed chain** (which TNFR2 did not test), but it does **NOT** prove
generalisation to a chain the model has never seen. ⛔ Only a de novo binder with measured affinity
and **no deposited complex** closes that, and we have exactly that locally in Adaptyv's
`status/competition_egfr/labels_402.csv` -- see `J505` §5. ⇒ this arm narrows the gap; it does not
shut it.

## THE TWO ARMS
⭐ **POSITIVE (de novo):** LCB1, `7JZU` chain `A`, 55 aa -- designed binder, measured affinity.
⚠️ **NEGATIVE:** the TNFα protomer, 152 aa, read from the same proven v41c FASTA our own folds used.
   A real, well-folded natural protein with **no RBD affinity**. It calibrates what "no interface"
   looks like against this target, which `J502` §3 showed is essential -- there the negative
   enriched 4.96× purely from SITE BIAS, so an absolute recall number means nothing without it.
   ⛔ It is **not length-matched** to LCB1 (152 vs 55 aa); `j501` could length-match its negative
   and this one cannot, because the comparison here is fixed by what LCB1 is. Stated, not hidden.

## ⛔ TWO DELIBERATE DEVIATIONS FROM THE `j496`/`j501` PATH, BOTH IN THE CONSERVATIVE DIRECTION
1. **The target gets an `empty` MSA, not a real one.** We hold a 2198-row TNFα MSA but have no RBD
   MSA locally, and building one is not a $0 local step. Our own designs' target had a real MSA, so
   this arm is run with **LESS** information than they got. ⇒ a PASS here is conservative and the
   conclusion is safe; a FAIL would be **confounded** by the missing target MSA and would need the
   MSA before it could be read as a de novo failure. ⛔ Do not read a fail as decisive.
2. **One target chain, not three.** TNFα is a trimer and the RBD is a monomer; LCB1's target simply
   is a single chain. ⇒ the 2-of-3-protomer groove geometry of the TNFα arms has no analogue here,
   and `iptm` has only one binder↔target pair (which incidentally removes the multimeric confound
   of).
Everything else is byte-identical to the proven path: binder on `empty`, same Boltz call, same
flags, same seed.

usage: j505_make_denovo_control_inputs_20261018.py <out-dir>
"""
import os
import sys

REPO = "<placeholder>"
PROVEN_FASTA = f"{REPO}/boltz/j457_v41c_fold_20261016/inputs/bound/v41c_s11_seq01_bound.fasta"
SEVENJZU = f"{REPO}/data/J130_inputs_20261010_7jzu.pdb"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}


def chain_seq(path, chain):
    out, nums = [], []
    for line in open(path):
        if line.startswith("ATOM") and line[21] == chain and line[12:16] == " CA ":
            out.append(AA3.get(line[17:20].strip(), "X"))
            nums.append(int(line[22:26]))
    assert out, f"no CA for chain {chain} in {path}"
    gaps = [(nums[i], nums[i + 1]) for i in range(len(nums) - 1) if nums[i + 1] != nums[i] + 1]
    return "".join(out), nums, gaps


def proven_tnfa():
    """The 152-aa TNFα string our own folds used, identified by MSA header, never by letter."""
    seqs, order, cur, is_tgt = {}, [], None, {}
    for line in open(PROVEN_FASTA):
        line = line.strip()
        if line.startswith(">"):
            p = line[1:].split("|")
            cur = p[0]
            seqs[cur], is_tgt[cur] = "", (len(p) > 2 and p[2] != "empty")
            order.append(cur)
        elif cur:
            seqs[cur] += line
    tch = [c for c in order if is_tgt[c]]
    assert len(tch) == 3, f"expected 3 MSA-bearing chains, got {tch}"
    s = seqs[tch[0]]
    assert len(s) == 152, f"proven target {len(s)} aa, expected 152"
    return s


def main():
    out_dir = sys.argv[1]
    os.makedirs(out_dir, exist_ok=True)

    lcb1, lnums, lgaps = chain_seq(SEVENJZU, "A")
    rbd, rnums, rgaps = chain_seq(SEVENJZU, "B")
    tnfa = proven_tnfa()
    assert not lgaps, f"LCB1 numbering gaps {lgaps} -- a gapped chain is not a clean control"
    assert not rgaps, f"RBD numbering gaps {rgaps}"
    assert "X" not in lcb1 + rbd, "unparsed residue in 7JZU"
    assert len(lcb1) == 55, f"LCB1 is {len(lcb1)} aa, expected 55 (J417 measured 55)"
    assert len(rbd) == 193, f"RBD is {len(rbd)} aa, expected 193 (j501 measured 193)"

    print(f"target  RBD        {len(rbd):>4} aa  resnum {rnums[0]}-{rnums[-1]}  (7JZU chain B)")
    print(f"⛔ target MSA = empty (no RBD MSA locally) -- LESS information than our own folds got")
    print()
    arms = [
        ("pos_lcb1_denovo", lcb1, "POSITIVE de novo: LCB1, designed, measured affinity (7JZU ch A)"),
        ("neg_tnfa_vs_rbd", tnfa, "NEGATIVE: TNFα protomer, folded, no RBD affinity (proven FASTA)"),
    ]
    for tag, binder, label in arms:
        fa = os.path.join(out_dir, f"{tag}_bound.fasta")
        with open(fa, "w") as fh:
            fh.write(f">A|protein|empty\n{binder}\n")
            fh.write(f">B|protein|empty\n{rbd}\n")
        print(f"  {tag:<18} {len(binder):>4} aa binder + {len(rbd)} target = "
              f"{len(binder)+len(rbd)} res   cys={binder.count('C')}")
        print(f"                     {label}")
        print(f"                     -> {fa}")
    print()
    print("⛔ Verify an OUTPUT COUNT, never rc -- a Boltz OOM exits 0.")
    print("   Expect 2 of 2 *_model_0.pdb files.")


main()
