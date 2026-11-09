# s118 -- J1056 2x2 REFOLD, scored PER SHARD with each shard's own in-batch floor

## ⭐ THE RESOLVED CONFIG, recorded here because `s60`'s own header is STALE
- `[M]` sequence design: SolubleMPNN via `s53`, **NO `--bias_AA`** (the unbiased measurement), `--omit_AA C`, T=0.1, 8 seqs/backbone, best by `overall_confidence`
- `[M]` fold: Boltz, `msa: empty` on every chain, `diffusion_samples 5`, recycling 3, **seeds 01-04, one per shard**
- `[M]` chains: A=**HLA-A\*11:01**(276) B=b2m(100) C=peptide G12D `VVVGADGVGK`(10) D=binder ⇒ binder<->peptide is `pair_chains_iptm(3,2)`
- ⛔⛔ **THE ALLELE LABEL IN THIS ROUND'S OWN `manifest.json` IS WRONG.** `[M]` `the pinned fold runner`:381 hard-codes *"A=HLA-A\*03:01"*, and this round reused that runner **byte-for-byte** (deliberately — its sha1 is pinned and pods refuse on `FATAL_SHA1_MISMATCH`), so all four shard manifests carry the wrong allele. ⭐ The true allele is **A\*11:01** (`J960`:1; `s35`:15-16). ⛔ The SEQUENCES are correct and unaffected — they are lifted from `s35`'s source by `ast` and `s116` asserted per design that the folded HLA/b2m/peptide chains match it exactly — so this is a FALSE LABEL in an artefact, ⛔ not a scientific error. Same class as `s83`'s hard-coded banner; ⛔ do not "fix" the runner, its bytes are load-bearing.
- ⛔ `s60`:107-108 prints *"ADOPTED r4 bias cohort"* and *"seed 20261104"*. ⛔⛔ BOTH FALSE here. The per-shard `s60` reports are scratch intermediates and must not be cited.
- ⭐ `s115` picked **one design per diffusion batch**, so the 20 designs per arm are 20 INDEPENDENT contigs ⇒ design-level == batch-level.

## 0 · PER-SHARD FLOORS — ⛔ never pooled (the seed IS the shard)

| shard | seed | designs | non-binders | STRICT floor cont_pep | STRICT floor D-C ipTM | lenient cont_pep | lenient ipTM |
|---|---|---|---|---|---|---|---|
| 01 | 1 | 20 | 2 | 10.4 | 0.888 | 9.8 | 0.779 |
| 02 | 2 | 20 | 2 | 11.2 | 0.882 | 10.0 | 0.788 |
| 03 | 3 | 20 | 2 | 11.0 | 0.888 | 10.8 | 0.782 |
| 04 | 4 | 20 | 2 | 10.8 | 0.892 | 8.6 | 0.795 |

- `[M]` **80 designs** aggregated over 4 shards, each floored by its OWN 2 non-binders.
- ⚠️ At n=2 per shard a floor is POORLY DETERMINED. The STRICT floor (higher non-binder) is the only headline; the lenient column is a sensitivity line and ⛔ was fixed as such before the counts below were seen.

## 1 · ⭐⭐ THE PRE-REGISTERED CONJUNCTION (`J1055` §3)

| arm | motif | breadth | n | **conj (STRICT)** | contacts alone | ipTM alone | conj (lenient) |
|---|---|---|---|---|---|---|---|
| `B_motif_broad` | held | broad | 20 | **0** | 2 | 0 | 0 |
| `MB_nomotif_broad` | none | broad | 20 | **2** | 6 | 2 | 2 |
| `M_nomotif_narrow` | none | narrow | 20 | **0** | 4 | 0 | 1 |
| `R_motif_narrow` | held | narrow | 20 | **0** | 1 | 0 | 0 |

### ⇒ **2 of 80 designs clear the pre-registered conjunction** against their own in-batch floor.

## 2 · THE 2x2 ON THE ENDPOINT, and on the CONTACT MARGIN

| | narrow | broad |
|---|---|---|
| **motif held** | 0/20 | 0/20 |
| **motif none** | 0/20 | 2/20 |

| factor | level A | level B | mean contact margin A | B | Δ | mean ipTM margin A | B |
|---|---|---|---|---|---|---|---|
| motif | held (n=40) | none (n=40) | -6.51 | -4.20 | -2.31 | -0.449 | -0.372 |
| breadth | broad (n=40) | narrow (n=40) | -5.02 | -5.69 | +0.67 | -0.399 | -0.422 |

## 3 · ⭐ THE PROMOTION LIST — ranked on CONTACT MARGIN, ipTM as TIE-BREAK ONLY

⛔ Ranking on `pair(3,2)` would rank on the confounded quantity. The confound only INFLATES, so ipTM can REJECT and never PROMOTE ⇒ it breaks ties and nothing more.
⭐ Clearing the floor makes a design a **CANDIDATE, never a result** (`J1055` §3).

| # | tag | arm | shard | cont_pep | floor | **margin** | D-C ipTM | floor | binder pLDDT |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `D_s03_MB_nomotif_broad_2_model_3` | MB_nomotif_broad | 01 | 16.4 | 10.4 | **+6.0** | 0.957 | 0.888 | 97.5 |
| 2 | `D_s03_MB_nomotif_broad_4_model_0` | MB_nomotif_broad | 03 | 11.6 | 11.0 | **+0.6** | 0.909 | 0.888 | 93.3 |

- wrote `t151_2x2_refold_conjunction_20261108_per_design.tsv` — 80 per-design rows
- ⛔ the four `s60_shard*_SCRATCH.txt` files carry `s60`'s STALE header and are NOT the record.
