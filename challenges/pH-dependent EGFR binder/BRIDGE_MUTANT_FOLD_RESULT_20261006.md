# The bridge mutants work — and they do not change the novelty level

Measured 2026-10-06.

The question this answered: can the sequences be modified so that they pass? The sequence half was
modified successfully. **The level did not move, and this is a measurement rather than an argument.**

## 1. What was built and folded

SolubleMPNN with the 7OBI residue **forbidden per position** (`--omit_AA_per_residue`) at 35–44
non-interface positions, giving variant B mutants of `id85`, `id1104` and `id47`.

Folded **monomer-only**, `boltz 2.2.1`, `--diffusion_samples 21 --recycling_steps 3 --no_kernels`,
**3 seeds (1001/2002/3003) giving 63 draws per molecule**. The parents were folded under identical
settings as the control: the pre-existing parent structures came from complexes, and comparing those
to fresh monomer folds would have confounded the mutation with the context.

**Verified by count, never by exit code:** 378/378 PDBs, 18/18 prediction directories, and three
reports of "Number of failed examples: 0".

## 2. The result

| molecule | draws | TM median | TM p10 | frac ≥0.80 | seq-id | coverage | structural | level |
|---|---|---|---|---|---|---|---|---|
| `id85` parent | 63 | 0.973 | 0.968 | **1.00** | 64.2% | 0.81 | HIGH | 2 |
| `id85` **mutB** | 63 | **0.975** | 0.967 | **1.00** | **29.4%** | 0.77 | HIGH | **2** |
| `id1104` parent | 63 | 0.975 | 0.971 | **1.00** | 59.5% | 0.77 | HIGH | 2 |
| `id1104` **mutB** | 63 | **0.970** | 0.964 | **1.00** | **30.0%** | 0.77 | HIGH | **2** |
| `id47` parent | 63 | 0.972 | 0.966 | **1.00** | 59.0% | 0.81 | HIGH | 2 |
| `id47` **mutB** | 63 | **0.941** | 0.918 | **1.00** | **30.7%** | 0.77 | HIGH | **2** |

**The sequence half moved, exactly as designed:** 59–64% to **29.4–30.7%** (foldseek local; MMseqs2
global on the same mutants is ~19%). Zero interface residues changed.

**The structural half did not move:** ΔTM = **+0.002 / −0.005 / −0.030**, and every one of the
**189 mutant draws** is ≥0.80. The fold is untouched.

"High structural similarity" is an independent sufficient condition for Level 2, so the mutants
remain Level 2 at 19–30% sequence identity. This is the rule working as written, not a failure of the
mutations: we asked MPNN to change the sequence *while keeping the fold*, and it did exactly that.
**The design goal and the novelty rule are in direct opposition here.**

## 3. The one lever left, and its arithmetic

Both **high** (TM ≥ 0.8) and **moderate** (TM ≥ 0.5) require **more than 70% of the sequence
covered**. Below 70% coverage the structural term is "less than moderate" — and that is the only
route to Level 3 or 4 for a molecule whose fold is a known one.

- **Measured coverage: 0.77** for the mutants, 0.77–0.81 for the parents. The matched region is
  **~112 of 145 residues**.
- To reach below 0.70 while keeping the TPR core: 112/L < 0.70 requires L > 160, so **16–20 aa** of
  added sequence that does not match a known structure. That addition is distal and need not touch
  the paratope, unlike every other option available.

**This is unverified and is not a plan.** Three things could defeat it:

1. Coverage may be segmented with a three-predictor domain consensus (TED) and measured **by
   domains**. An appended segment recognised as a known domain would count, and coverage would not
   fall. The appendage must be structurally novel *and* survive domain segmentation as separate.
2. A ~20 aa appendage on a 145 aa binder is a **new molecule with no binding evidence**, and nothing
   in the remaining time re-establishes that.
3. Our coverage proxy is `alnlen/qlen` from a whole-chain foldseek search, **not** a per-domain
   weighted figure. The 0.77 is indicative, not directly comparable.

## 4. What this means

- **The nine submitted site-A designs need none of this.** They are Level 3 on both halves already: structure
  moderate rather than high, and sequence clear — see `SEQUENCE_NOVELTY_SWISSPROT_20261006.md`, which
  reports zero significant hits in Swiss-Prot or the PDB.
  - **The range "TM 0.62–0.76" is retracted, corrected on the monomer fold.** It is supported by
    neither measurement. The complex-derived table spans **0.567–0.885** — and that span includes the
    one design flagged as high similarity, which the retracted sentence silently dropped — while the
    monomer medians span **0.574–0.713**. See `NOVELTY_CORRECTION_MONOMER_20261006.md`. The
    **conclusion** that the designs are moderate rather than high survives; only the range used to
    state it was never measured.
- **The three bridges (`egfr-bridge-id85`, `egfr-bridge-id1104`, `egfr-bridge-id47`) are Level 2
  both before and after mutation.** Before, they failed both halves
  (`id85` against `1na0_A`, CTPR3, at E = 4.4e-33); after, they fail only the structural half. That is
  a real improvement in the honest description of the molecule, and it is **not** a level change.
- The mutants are worth keeping as designs regardless: a bridge at 19–30% identity to the CTPR series
  is defensibly our own molecule in a way the parent is not.

## 5. Files

`data/bridge_vs_ctpr_20261006/` holds `bridge_mutants_20261006.csv` (6 mutants, B and C variants),
`mutable_sites_20261006.csv`, `per_residue_identity_vs_7obi.csv`, `bridges_vs_7obi_20261006.pse` and
the renders. The 378 fold PDBs and the 328,256-row foldseek hit table are **not deposited** — they
are too large to distribute — and the per-draw numbers they support are in the table in §2 above.
