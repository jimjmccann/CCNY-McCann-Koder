# Sequence-half novelty measurement — Adaptyv/Proteinbase Level-3 gate
**Measured 2026-10-06, on local CPU.**

The structural half is already measured by us: every one of the **16 queries** has **at least
moderate** structural similarity to a known protein. ⛔ **"16" here counts QUERIES, not submitted
designs** — it is the 12-design verdict cohort plus 4 additional candidates, so it is neither the
submitted set nor the 12. The word "designs" was used loosely in an earlier draft and is corrected
throughout this file to "queries". So the Proteinbase level turns entirely
on one number per design:

| sequence similarity | + moderate structure | level | gate (>= 3) |
|---|---|---|---|
| <= 30% | moderate | **3** | PASS |
| > 30% | moderate | **2** | FAIL |

## ⇒ THE ANSWER — the 9 **novelty-screened site-A candidates** (the "PASS_L3*" set)

⛔ **This nine is NOT the same nine site-A designs that were submitted.** The two sets differ by
exactly one member and are easy to confuse, because both are naturally described as "nine site-A
designs". Enumerated here so no reader has to infer which is meant:

| set | members |
|---|---|
| **PASS_L3\* set** — what this table screens | `id2015` `id3017` `id3020` `id3035` `id3036` `id89` `id2046` **`cand-siteA-id2`** `id86` |
| **submitted site-A nine** — what was actually sent | `id2015` `id3017` `id3020` `id3035` `id3036` `id89` `id2046` **`id2048`** `id86` |

⇒ `cand-siteA-id2` was screened here but **never submitted**; `id2048` **was** submitted but is
not in this table. ⛔ **A claim that "the nine site-A designs passed" is therefore true of the
screening cohort and NOT of the submitted set.** Wherever this document says "the 9", it means the
PASS_L3* set above.

**All 9 of the PASS_L3\* set are BELOW the 30% threshold and all 9 still pass the Level-3 gate on
the sequence axis.** None of them has a single statistically significant hit in Swiss-Prot or the PDB, by
either of two independent search tools, at any sensitivity we could push.

| design | struct_tier | significant hits (MMseqs2 SP/PDB, phmmer) | max GLOBAL id among significant hits | sequence axis |
|---|---|---|---|---|
| `egfr-siteA-id2015` | PASS_L3 | 0 / 0, 0 | none exist | ✅ **PASS** |
| `egfr-siteA-id3017` | PASS_L3 | 0 / 0, 0 | none exist (but see below) | ✅ **PASS** |
| `egfr-siteA-id3020` | PASS_L3 | 0 / 0, 0 | none exist | ✅ **PASS** |
| `egfr-siteA-id3035` | PASS_L3 | 0 / 0, 0 | none exist (but see below) | ✅ **PASS** |
| `egfr-siteA-id3036` | PASS_L3 | 0 / 0, 0 | none exist | ✅ **PASS** |
| `egfr-siteA-id89` | PASS_L3 | 0 / 0, 0 | none exist | ✅ **PASS** |
| `egfr-siteA-id2046` | PASS_L3_borderline | 0 / 0, 0 | none exist | ✅ **PASS** |
| `cand-siteA-id2` | PASS_L3_candidate | 0 / 0, 0 | none exist | ✅ **PASS** |
| `egfr-siteA-id86` | PASS_L3_seq_at_risk | 0 / 0, 0 | none exist | ✅ **PASS** |

The best *insignificant* hit for these designs sits at **25.3% global identity** (MMseqs2,
E = 0.8 to 134). The highest identity reachable at **any** E-value, including pure noise out to
E = 10000, is **33.3-40.0% global** (30.0% for id86) — see the three caveats below, which is
where the real risk lives.

**The 4 bridges fail the sequence axis outright**, independently of the structural verdict:
41.4-46.2% **global** identity (and 61-68% local) at E = 1e-28 to 1e-33 against *designed
consensus TPR* structures — `1na0_A`/`1na3_A` (CTPR3/CTPR2), `7obi_A` (CTPR-rv4), `2fo7_A`,
`8chy_A`, `9qhb_A`. Both identity conventions clear 30% by a wide margin, so there is no
convention under which the bridges pass. This **confirms and strengthens** the existing
`FAIL_L2_bridge` tag.

---

## 1. Exactly what was built

| component | value |
|---|---|
| primary tool | **MMseqs2 18.8cc5c** (bioconda), in a dedicated conda env |
| cross-check tool | **HMMER 3.4 `phmmer`** (Aug 2023) |
| aligner cross-check | **Biopython 1.84** `PairwiseAligner`, Smith-Waterman |
| DB 1 | **UniProtKB/Swiss-Prot release 2026_03** (02-Sep-2026) — **575,748** sequences |
| DB 2 | **PDB `pdb_seqres`** (wwPDB, file dated 25-Sep-2026) — 1,169,269 chains, **1,103,516** after `mol:protein` filter |

MMseqs2 was installed into a **new** env; the pre-existing `mmseqs` env was not modified
(it happens to hold the same version, 18.8cc5c). Raw DB files live under `$DESIGN_DATA_ROOT`
(`novelty_db_20261006/`), on a scratch volume rather than the root filesystem.

### Search parameters

MMseqs2 `easy-search`, both DBs:
```
- s 7.5 -e 10000 --max-seqs 5000 --alignment-mode 3 -a -c 0 --min-seq-id 0
- -max-seq-len 100000 --threads 8 --format-mode 4
```
BLOSUM62, `--gap-open aa:11 --gap-extend aa:1`, `--comp-bias-corr 1`. `-s 7.5` is maximum
sensitivity; `-e 10000` / `-c 0` / `--min-seq-id 0` make this simultaneously the sensitive
search *and* the no-threshold control, so a "no hits" answer can never be an artifact of an
insensitive search.

phmmer, both DBs:
```
- -max -E 1000 --domE 1000 --incE 1000 --incdomE 1000 --cpu 16 --notextw
```
`--max` disables all three heuristic filters (MSV/Viterbi/Forward) for maximal sensitivity.

Independent realignment: Biopython `PairwiseAligner`, `mode="local"`, BLOSUM62,
`open_gap_score=-11`, `extend_gap_score=-1` (matching MMseqs2's matrix and gap costs).

### Sensitivity actually achieved — "no hits" is not the finding

Every one of the 16 queries returned **hundreds to thousands** of permissive hits in both
DBs (MMseqs2: 19-1452 per query; phmmer: 3,401-11,534). The finding is never "the search returned
nothing" — it is "the search returned plenty, and none of it is significant or above 30%".

| | MMseqs2 hit rows | phmmer domain rows |
|---|---|---|
| Swiss-Prot | 9,474 | 86,464 |
| PDB seqres | 8,567 | 106,454 |

### Two methodological traps hit during this work

1. **`-a` (backtrace) is MANDATORY in MMseqs2.** Without it, `fident` and `nident` are
   silently reported as **0.000 / 0** for every hit instead of erroring — which reads as
   "0% identity to everything". The first full run produced exactly that and was discarded.
   Measured 2026-10-06; the fix is recorded in the script's own comments.
2. **`nident/qlen` must be derived from the LOCAL alignment.** Forcing a global
   Needleman-Wunsch alignment against a longer target scatters identities across the whole
   target and inflates the global figure: `id3017` vs `7udk_A` reads **30.7%** global from the
   local alignment but **40.0%** from an NW alignment. The local-derived convention is the one
   that matches the worked example in the brief (67.6% over 65/145 → ~30% global), and it is
   what is used throughout this report.

### Verification (a count is not a measurement)

- Downloads are real data, not HTML error pages: `file` reports gzip; first Swiss-Prot header
  is `sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 ...`; first PDB
  header is `100d_A mol:na length:10`. Sequence counts as tabled above.
- **All 18,041 hits had `nident` and `alnlen` recomputed from the `qaln`/`taln` alignment
  strings: 0 mismatches on both.** MMseqs2's identity columns reproduce exactly.
- Independent Smith-Waterman vs MMseqs2 local identity over 320 shared query-target pairs:
  median difference **1.52 pp**, p90 15.6 pp, max 41.0 pp. The large tail is expected — SW
  picks its own optimal local alignment and sometimes selects a different, shorter/longer
  window than MMseqs2's seeded alignment; it never changed a gate verdict.
- No `2>/dev/null` anywhere in any script. stderr was shown for every command.

---

## 2. Which identity is which

For every hit both numbers are reported:

- **local identity** = `nident / alnlen` — identical residues over the *aligned* region.
  This is MMseqs2's native `fident` column (and `pident` = the same × 100).
- **global identity** = `nident / qlen` — identical residues over the *full query*, taking
  `nident` from the local alignment.

phmmer's tabular output has **no identity column at all**, so every phmmer identity in this
report was recomputed by Smith-Waterman realignment.

The gap between the two conventions is decisive here: `cand-bridge-id76` vs `1na3_A` reads
**61.3% local but 31.7% global**; `egfr-siteA-id89` vs `9k7w_F` reads **47.5% local but 25.3%
global**. Unfiltered "max local identity" is worthless on its own — it is dominated by 9-16
residue fragments at 80-100% identity (e.g. `egfr-bridge-id85` hits `8fee_A` at **100.0% local
over 9 residues** = 6.2% global, E = 6e3). Every local figure quoted below carries a
`qcov >= 50%` floor for that reason.

---

## 3. Full results over all 16 queries

### Table 1 — significant hits (E <= 1e-3) and the gate verdict

| design | len | struct_tier | MMseqs2 sig hits (SP / PDB) | phmmer domain-sig hits | max LOCAL id (sig) | max GLOBAL id (sig) | sequence axis |
|---|---|---|---|---|---|---|---|
| `cand-siteA-id19-sh09` | 75 | FAIL_L2_both | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `cand-bridge-id76` | 145 | FAIL_L2_bridge | 187 / 391 | 2 | 61.3% | 46.2% | **FAIL** (46.2% > 30%) |
| `egfr-bridge-id1104` | 145 | FAIL_L2_bridge | 218 / 329 | 2 | 64.6% | 44.1% | **FAIL** (44.1% > 30%) |
| `egfr-bridge-id47` | 145 | FAIL_L2_bridge | 181 / 308 | 2 | 64.0% | 41.4% | **FAIL** (41.4% > 30%) |
| `egfr-bridge-id85` | 145 | FAIL_L2_bridge | 312 / 370 | 2 | 67.5% | 46.2% | **FAIL** (46.2% > 30%) |
| `cand-siteA-id73` | 80 | FAIL_L2_seq | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id2048` | 75 | FAIL_L2_struct_coinflip | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id2015` | 75 | PASS_L3 | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id3017` | 75 | PASS_L3 | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id3020` | 75 | PASS_L3 | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id3035` | 75 | PASS_L3 | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id3036` | 75 | PASS_L3 | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id89` | 75 | PASS_L3 | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id2046` | 75 | PASS_L3_borderline | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `cand-siteA-id2` | 75 | PASS_L3_candidate | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |
| `egfr-siteA-id86` | 80 | PASS_L3_seq_at_risk | 0 / 0 | 0 | n/a | n/a | **PASS** (no significant hit) |

### Table 2 — no-threshold CONTROL: the single worst hit at any E-value (up to E=10000)

| design | max LOCAL id (qcov>=50%) | its hit / E | max GLOBAL id | its hit / E | above 30% global? |
|---|---|---|---|---|---|
| `cand-siteA-id19-sh09` | 42.3% | 6mfv_A / 9e+02 | 50.7% | 8qai_A / 0.0029 (ph) | yes |
| `cand-bridge-id76` | 61.3% | 1na3_A / 7.3e-22 | 46.2% | 7obi_A / 3.3e-30 | yes |
| `egfr-bridge-id1104` | 64.6% | 8ur5_A / 1.3e-27 | 44.1% | 1na0_A / 3.2e-33 | yes |
| `egfr-bridge-id47` | 64.0% | 8ur5_A / 3.8e-25 | 41.4% | 7obi_A / 6.1e-45 (ph) | yes |
| `egfr-bridge-id85` | 67.5% | 1na3_A / 1.1e-25 | 46.2% | 9qhb_A / 2.2e-186 (ph) | yes |
| `cand-siteA-id73` | 38.6% | B1YAJ4 / 1.7e+03 | 35.0% | A5HYY0 / 6.1e+03 | yes |
| `egfr-siteA-id2048` | 42.0% | Q5JYT7 / 1.1e+03 | 40.0% | sp|A6VAK1|CMOA_PSEP7 / 24 (ph) | yes |
| `egfr-siteA-id2015` | 44.0% | Q5JYT7 / 4.3e+02 | 40.0% | sp|A6VAK1|CMOA_PSEP7 / 6.9 (ph) | yes |
| `egfr-siteA-id3017` | 41.7% | Q24UT2 / 7.3e+03 | 38.7% | sp|Q04QS6|SYI_LEPBJ / 9.6 (ph) | yes |
| `egfr-siteA-id3020` | 42.0% | Q5JYT7 / 1.1e+03 | 40.0% | sp|A6VAK1|CMOA_PSEP7 / 10 (ph) | yes |
| `egfr-siteA-id3035` | 41.7% | Q24UT2 / 7.3e+03 | 37.3% | sp|Q04QS6|SYI_LEPBJ / 11 (ph) | yes |
| `egfr-siteA-id3036` | 43.2% | Q8EEE7 / 3.1e+02 | 40.0% | sp|A6VAK1|CMOA_PSEP7 / 27 (ph) | yes |
| `egfr-siteA-id89` | 47.5% | 9k7w_F / 0.82 | 33.3% | sp|A4SPN7|CMOA_AERS4 / 0.61 (ph) | yes |
| `egfr-siteA-id2046` | 42.0% | Q5JYT7 / 1.1e+03 | 40.0% | sp|A6VAK1|CMOA_PSEP7 / 10 (ph) | yes |
| `cand-siteA-id2` | 47.6% | 23wn_A / 4.8e+02 | 36.0% | sp|A9WDL9|THIM_CHLAA / 2.7 (ph) | yes |
| `egfr-siteA-id86` | 39.5% | Q8A9W1 / 8.3e+03 | 30.0% | sp|Q9CEC9|UPP_LACLA / 22 (ph) | no |

### Table 3 — top hits per design (MMseqs2, by bitscore, both DBs)

| design | db | hit | protein | local id | global id | nident/alnlen | qcov | E-value |
|---|---|---|---|---|---|---|---|---|
| `cand-siteA-id19-sh09` | swissprot | Q58743 | Uncharacterized protein MJ1348 | 50.0% | 18.7% | 14/28 | 37.3% | 2.61 |
| `cand-siteA-id19-sh09` | swissprot | Q9USI5 | Heat shock protein sti1 homolog | 48.1% | 17.3% | 13/27 | 36.0% | 33.5 |
| `cand-siteA-id19-sh09` | swissprot | P15705 | Heat shock protein STI1 | 35.3% | 16.0% | 12/34 | 45.3% | 164 |
| `cand-bridge-id76` | swissprot | P56558 | UDP-N-acetylglucosamine--peptide N-acetylg | 42.2% | 24.1% | 35/83 | 57.2% | 7.9e-12 |
| `cand-bridge-id76` | swissprot | Q8CGY8 | UDP-N-acetylglucosamine--peptide N-acetylg | 42.2% | 24.1% | 35/83 | 57.2% | 7.9e-12 |
| `cand-bridge-id76` | swissprot | Q27HV0 | UDP-N-acetylglucosamine--peptide N-acetylg | 42.2% | 24.1% | 35/83 | 57.2% | 7.9e-12 |
| `egfr-bridge-id1104` | swissprot | P56558 | UDP-N-acetylglucosamine--peptide N-acetylg | 43.4% | 31.7% | 46/106 | 73.1% | 5.34e-17 |
| `egfr-bridge-id1104` | swissprot | Q8CGY8 | UDP-N-acetylglucosamine--peptide N-acetylg | 43.4% | 31.7% | 46/106 | 73.1% | 5.34e-17 |
| `egfr-bridge-id1104` | swissprot | Q27HV0 | UDP-N-acetylglucosamine--peptide N-acetylg | 43.4% | 31.7% | 46/106 | 73.1% | 5.34e-17 |
| `egfr-bridge-id47` | swissprot | Q58208 | TPR repeat-containing protein MJ0798 | 39.0% | 28.3% | 41/105 | 72.4% | 1.48e-11 |
| `egfr-bridge-id47` | swissprot | P56558 | UDP-N-acetylglucosamine--peptide N-acetylg | 37.5% | 26.9% | 39/104 | 71.7% | 2.76e-11 |
| `egfr-bridge-id47` | swissprot | Q8CGY8 | UDP-N-acetylglucosamine--peptide N-acetylg | 37.5% | 26.9% | 39/104 | 71.7% | 2.76e-11 |
| `egfr-bridge-id85` | swissprot | P56558 | UDP-N-acetylglucosamine--peptide N-acetylg | 41.5% | 30.3% | 44/106 | 73.1% | 3.87e-14 |
| `egfr-bridge-id85` | swissprot | Q8CGY8 | UDP-N-acetylglucosamine--peptide N-acetylg | 41.5% | 30.3% | 44/106 | 73.1% | 3.87e-14 |
| `egfr-bridge-id85` | swissprot | Q27HV0 | UDP-N-acetylglucosamine--peptide N-acetylg | 41.5% | 30.3% | 44/106 | 73.1% | 3.87e-14 |
| `cand-siteA-id73` | swissprot | O32156 | Uncharacterized ABC transporter extracellu | 43.2% | 20.0% | 16/37 | 46.2% | 0.339 |
| `cand-siteA-id73` | swissprot | Q5UQ72 | Putative replication factor C small subuni | 29.2% | 26.2% | 21/72 | 90.0% | 28.9 |
| `cand-siteA-id73` | swissprot | P74745 | Serine/threonine-protein kinase C | 40.7% | 13.8% | 11/27 | 33.8% | 102 |
| `egfr-siteA-id2048` | swissprot | Q4QLP6 | Magnesium transport protein CorA | 37.3% | 25.3% | 19/51 | 68.0% | 6.8 |
| `egfr-siteA-id2048` | swissprot | P44998 | Magnesium transport protein CorA | 37.3% | 25.3% | 19/51 | 68.0% | 6.8 |
| `egfr-siteA-id2048` | swissprot | Q15NL2 | Carboxy-S-adenosyl-L-methionine synthase | 40.9% | 24.0% | 18/44 | 58.7% | 12.9 |
| `egfr-siteA-id2015` | swissprot | Q15NL2 | Carboxy-S-adenosyl-L-methionine synthase | 43.2% | 25.3% | 19/44 | 58.7% | 6.8 |
| `egfr-siteA-id2015` | swissprot | Q4QLP6 | Magnesium transport protein CorA | 38.0% | 25.3% | 19/50 | 66.7% | 12.9 |
| `egfr-siteA-id2015` | swissprot | P44998 | Magnesium transport protein CorA | 38.0% | 25.3% | 19/50 | 66.7% | 12.9 |
| `egfr-siteA-id3017` | swissprot | Q4QLP6 | Magnesium transport protein CorA | 37.3% | 25.3% | 19/51 | 68.0% | 1.9 |
| `egfr-siteA-id3017` | swissprot | P44998 | Magnesium transport protein CorA | 37.3% | 25.3% | 19/51 | 68.0% | 1.9 |
| `egfr-siteA-id3017` | swissprot | Q0ZQ41 | CMP-5'-(3-aminopropyl)phosphonate synthase | 58.6% | 22.7% | 17/29 | 38.7% | 33.5 |
| `egfr-siteA-id3020` | swissprot | Q4QLP6 | Magnesium transport protein CorA | 38.0% | 25.3% | 19/50 | 66.7% | 4.95 |
| `egfr-siteA-id3020` | swissprot | P44998 | Magnesium transport protein CorA | 38.0% | 25.3% | 19/50 | 66.7% | 4.95 |
| `egfr-siteA-id3020` | swissprot | Q15NL2 | Carboxy-S-adenosyl-L-methionine synthase | 40.9% | 24.0% | 18/44 | 58.7% | 12.9 |
| `egfr-siteA-id3035` | swissprot | Q4QLP6 | Magnesium transport protein CorA | 37.3% | 25.3% | 19/51 | 68.0% | 1.9 |
| `egfr-siteA-id3035` | swissprot | P44998 | Magnesium transport protein CorA | 37.3% | 25.3% | 19/51 | 68.0% | 1.9 |
| `egfr-siteA-id3035` | swissprot | Q0ZQ41 | CMP-5'-(3-aminopropyl)phosphonate synthase | 58.6% | 22.7% | 17/29 | 38.7% | 33.5 |
| `egfr-siteA-id3036` | swissprot | Q4QLP6 | Magnesium transport protein CorA | 37.3% | 25.3% | 19/51 | 68.0% | 2.61 |
| `egfr-siteA-id3036` | swissprot | P44998 | Magnesium transport protein CorA | 37.3% | 25.3% | 19/51 | 68.0% | 2.61 |
| `egfr-siteA-id3036` | swissprot | Q15NL2 | Carboxy-S-adenosyl-L-methionine synthase | 38.6% | 22.7% | 17/44 | 58.7% | 33.5 |
| `egfr-siteA-id89` | swissprot | Q31JJ7 | Carboxy-S-adenosyl-L-methionine synthase | 51.4% | 25.3% | 19/37 | 49.3% | 2.61 |
| `egfr-siteA-id89` | swissprot | A0KIF6 | Carboxy-S-adenosyl-L-methionine synthase | 43.2% | 25.3% | 19/44 | 58.7% | 3.6 |
| `egfr-siteA-id89` | swissprot | A4SPN7 | Carboxy-S-adenosyl-L-methionine synthase | 53.6% | 20.0% | 15/28 | 37.3% | 9.36 |
| `egfr-siteA-id2046` | swissprot | Q4QLP6 | Magnesium transport protein CorA | 38.0% | 25.3% | 19/50 | 66.7% | 9.36 |
| `egfr-siteA-id2046` | swissprot | P44998 | Magnesium transport protein CorA | 38.0% | 25.3% | 19/50 | 66.7% | 9.36 |
| `egfr-siteA-id2046` | swissprot | Q15NL2 | Carboxy-S-adenosyl-L-methionine synthase | 40.9% | 24.0% | 18/44 | 58.7% | 12.9 |
| `cand-siteA-id2` | swissprot | P13445 | RNA polymerase sigma factor RpoS | 41.9% | 17.3% | 13/31 | 41.3% | 24.3 |
| `cand-siteA-id2` | swissprot | P0A2E7 | RNA polymerase sigma factor RpoS | 41.9% | 17.3% | 13/31 | 41.3% | 24.3 |
| `cand-siteA-id2` | swissprot | D0ZVL4 | RNA polymerase sigma factor RpoS | 41.9% | 17.3% | 13/31 | 41.3% | 24.3 |
| `egfr-siteA-id86` | swissprot | O67118 | Chaperone protein DnaK | 35.7% | 25.0% | 20/56 | 70.0% | 28.9 |
| `egfr-siteA-id86` | swissprot | O58969 | Uncharacterized ABC transporter extracellu | 50.0% | 17.5% | 14/28 | 35.0% | 39.6 |
| `egfr-siteA-id86` | swissprot | A0A509AL23 | Surface phospholipase | 45.5% | 18.8% | 15/33 | 41.2% | 39.6 |
| `cand-siteA-id19-sh09` | pdbseqres | 7obi_A | CTPR-rv4 | 34.9% | 38.7% | 29/83 | 110.7% | 10.5 |
| `cand-siteA-id19-sh09` | pdbseqres | 7obi_B | CTPR-rv4 | 34.9% | 38.7% | 29/83 | 110.7% | 10.5 |
| `cand-siteA-id19-sh09` | pdbseqres | 6deh_A | TPR repeat protein, protein-protein intera | 84.6% | 14.7% | 11/13 | 17.3% | 134 |
| `cand-bridge-id76` | pdbseqres | 8chy_A | Crystal structure of an 8-repeat consensus | 46.9% | 41.4% | 60/128 | 88.3% | 3.64e-31 |
| `cand-bridge-id76` | pdbseqres | 8bu0_A | Consensus tetratricopeptide repeat protein | 46.9% | 41.4% | 60/128 | 88.3% | 3.64e-31 |
| `cand-bridge-id76` | pdbseqres | 8cig_A | Consensus tetratricopeptide repeat protein | 46.9% | 41.4% | 60/128 | 88.3% | 3.64e-31 |
| `egfr-bridge-id1104` | pdbseqres | 2fo7_A | SYNTHETIC CONSENSUS TPR PROTEIN | 55.8% | 43.4% | 63/113 | 77.9% | 1.71e-33 |
| `egfr-bridge-id1104` | pdbseqres | 2hyz_A | SYNTHETIC CONSENSUS TPR PROTEIN | 55.8% | 43.4% | 63/113 | 77.9% | 1.71e-33 |
| `egfr-bridge-id1104` | pdbseqres | 8ch0_A | Consensus tetratricopeptide repeat protein | 55.8% | 43.4% | 63/113 | 77.9% | 1.71e-33 |
| `egfr-bridge-id47` | pdbseqres | 1na0_A | designed protein CTPR3 | 57.7% | 41.4% | 60/104 | 71.7% | 5.12e-28 |
| `egfr-bridge-id47` | pdbseqres | 1na0_B | designed protein CTPR3 | 57.7% | 41.4% | 60/104 | 71.7% | 5.12e-28 |
| `egfr-bridge-id47` | pdbseqres | 7obi_A | CTPR-rv4 | 56.1% | 41.4% | 60/107 | 73.8% | 9.62e-28 |
| `egfr-bridge-id85` | pdbseqres | 1na0_A | designed protein CTPR3 | 63.5% | 45.5% | 66/104 | 71.7% | 4.4e-33 |
| `egfr-bridge-id85` | pdbseqres | 1na0_B | designed protein CTPR3 | 63.5% | 45.5% | 66/104 | 71.7% | 4.4e-33 |
| `egfr-bridge-id85` | pdbseqres | 2fo7_A | SYNTHETIC CONSENSUS TPR PROTEIN | 55.2% | 44.1% | 64/116 | 80.0% | 1.13e-32 |
| `cand-siteA-id73` | pdbseqres | 5f7v_A | Lmo0181 protein | 52.6% | 12.5% | 10/19 | 23.8% | 44.5 |
| `cand-siteA-id73` | pdbseqres | 4zza_A | Sugar binding protein of ABC transporter s | 66.7% | 12.5% | 10/15 | 18.8% | 115 |
| `cand-siteA-id73` | pdbseqres | 4zza_B | Sugar binding protein of ABC transporter s | 66.7% | 12.5% | 10/15 | 18.8% | 115 |
| `egfr-siteA-id2048` | pdbseqres | 9dtd_A | C2-01018 | 48.1% | 17.3% | 13/27 | 36.0% | 70.9 |
| `egfr-siteA-id2048` | pdbseqres | 9dtd_B | C2-01018 | 48.1% | 17.3% | 13/27 | 36.0% | 70.9 |
| `egfr-siteA-id2048` | pdbseqres | 4iwn_A | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 253 |
| `egfr-siteA-id2015` | pdbseqres | 9dtd_A | C2-01018 | 48.1% | 17.3% | 13/27 | 36.0% | 134 |
| `egfr-siteA-id2015` | pdbseqres | 9dtd_B | C2-01018 | 48.1% | 17.3% | 13/27 | 36.0% | 134 |
| `egfr-siteA-id2015` | pdbseqres | 4iwn_A | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 348 |
| `egfr-siteA-id3017` | pdbseqres | 4iwn_A | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 253 |
| `egfr-siteA-id3017` | pdbseqres | 4iwn_B | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 253 |
| `egfr-siteA-id3017` | pdbseqres | 3t6g_B | Breast cancer anti-estrogen resistance pro | 56.0% | 18.7% | 14/25 | 33.3% | 253 |
| `egfr-siteA-id3020` | pdbseqres | 9dtd_A | C2-01018 | 51.9% | 18.7% | 14/27 | 36.0% | 51.6 |
| `egfr-siteA-id3020` | pdbseqres | 9dtd_B | C2-01018 | 51.9% | 18.7% | 14/27 | 36.0% | 51.6 |
| `egfr-siteA-id3020` | pdbseqres | 4iwn_A | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 253 |
| `egfr-siteA-id3035` | pdbseqres | 4iwn_A | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 253 |
| `egfr-siteA-id3035` | pdbseqres | 4iwn_B | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 253 |
| `egfr-siteA-id3035` | pdbseqres | 3t6g_B | Breast cancer anti-estrogen resistance pro | 56.0% | 18.7% | 14/25 | 33.3% | 253 |
| `egfr-siteA-id3036` | pdbseqres | 9dtd_A | C2-01018 | 51.9% | 18.7% | 14/27 | 36.0% | 27.3 |
| `egfr-siteA-id3036` | pdbseqres | 9dtd_B | C2-01018 | 51.9% | 18.7% | 14/27 | 36.0% | 27.3 |
| `egfr-siteA-id3036` | pdbseqres | 4iwn_A | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 253 |
| `egfr-siteA-id89` | pdbseqres | 9k7w_F | C3-Ni1-HH*3-18 | 47.5% | 25.3% | 19/40 | 53.3% | 0.818 |
| `egfr-siteA-id89` | pdbseqres | 9k7w_A | C3-Ni1-HH*3-18 | 47.5% | 25.3% | 19/40 | 53.3% | 0.818 |
| `egfr-siteA-id89` | pdbseqres | 9k7w_B | C3-Ni1-HH*3-18 | 47.5% | 25.3% | 19/40 | 53.3% | 0.818 |
| `egfr-siteA-id2046` | pdbseqres | 9dtd_A | C2-01018 | 51.9% | 18.7% | 14/27 | 36.0% | 51.6 |
| `egfr-siteA-id2046` | pdbseqres | 9dtd_B | C2-01018 | 51.9% | 18.7% | 14/27 | 36.0% | 51.6 |
| `egfr-siteA-id2046` | pdbseqres | 4iwn_A | tRNA (cmo5U34)-methyltransferase | 58.8% | 13.3% | 10/17 | 22.7% | 253 |
| `cand-siteA-id2` | pdbseqres | 6omf_F | RNA polymerase sigma factor RpoS | 41.9% | 17.3% | 13/31 | 41.3% | 37.5 |
| `cand-siteA-id2` | pdbseqres | 6kj6_F | RNA polymerase sigma factor RpoS | 41.9% | 17.3% | 13/31 | 41.3% | 37.5 |
| `cand-siteA-id2` | pdbseqres | 6uu4_FFF | RNA polymerase sigma factor RpoS | 41.9% | 17.3% | 13/31 | 41.3% | 37.5 |
| `egfr-siteA-id86` | pdbseqres | 8bfy_A | Putative secreted cellobiose-binding (Tran | 57.1% | 15.0% | 12/21 | 26.2% | 32.4 |
| `egfr-siteA-id86` | pdbseqres | 4zza_A | Sugar binding protein of ABC transporter s | 68.8% | 13.8% | 11/16 | 20.0% | 61.1 |
| `egfr-siteA-id86` | pdbseqres | 4zza_B | Sugar binding protein of ABC transporter s | 68.8% | 13.8% | 11/16 | 20.0% | 61.1 |

### PASS_L3 cohort summary

- `egfr-siteA-id2015` (PASS_L3): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 40.0%  -> **PASS** (no significant hit)
- `egfr-siteA-id3017` (PASS_L3): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 38.7%  -> **PASS** (no significant hit)
- `egfr-siteA-id3020` (PASS_L3): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 40.0%  -> **PASS** (no significant hit)
- `egfr-siteA-id3035` (PASS_L3): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 37.3%  -> **PASS** (no significant hit)
- `egfr-siteA-id3036` (PASS_L3): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 40.0%  -> **PASS** (no significant hit)
- `egfr-siteA-id89` (PASS_L3): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 33.3%  -> **PASS** (no significant hit)
- `egfr-siteA-id2046` (PASS_L3_borderline): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 40.0%  -> **PASS** (no significant hit)
- `cand-siteA-id2` (PASS_L3_candidate): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 36.0%  -> **PASS** (no significant hit)
- `egfr-siteA-id86` (PASS_L3_seq_at_risk): significant hits = NONE; sig GLOBAL = n/a; control worst GLOBAL = 30.0%  -> **PASS** (no significant hit)

---

## 4. Three caveats — the things that could still bite

### Caveat 1 — `id3017` and `id3035` touch 30.7% against a designed helical repeat protein

These two designs are the only PASS designs with *any* nominally significant hit, and it lands
**0.7 points above** the threshold:

| | value |
|---|---|
| hit | `7udk_A` / `7udo_A` — **"Designed helical repeat protein (DHR) RPB_LRP2_R4"**, 172 aa |
| identity | **39.7% local / 30.7% global** (23 identical / 58 aligned, 23/75 of query) |
| phmmer full-sequence E | **1.7e-06** (significant) |
| phmmer best *domain* E | **36** (NOT significant; the four repeat domains score 36-130) |
| MMseqs2 E | **253** (not significant; MMseqs2 does not call this a hit at all) |

**Why this is judged not to be a real homology hit:** the full-sequence E-value of 1.7e-06 is
assembled from four individually *insignificant* repeat alignments against a repetitive target.
That is the classic HMMER repeat-accumulation artifact — no single alignment is significant, but
summing four weak ones crosses the full-sequence threshold. MMseqs2, which does not accumulate
this way, puts the same pair at E = 253. The same artifact inflates `cand-siteA-id19-sh09`
(`9qhb_A`, 564 aa engineered TPR, full E = 3.1e-08 from 12+ repeat matches each at domE 68-1100).

**But it is flagged because:** 30.7% > 30%, `7udk_A` is a Baker-lab *designed* helical repeat
protein and our designs are de novo helical bundles, so the resemblance is to the idealized
designed-helix sequence idiom rather than to a natural protein — and that is exactly the class
of similarity a novelty checker built for de novo designs might be tuned to catch. If
Proteinbase reports a full-sequence E-value without a per-domain check, `id3017` and `id3035`
are the two designs that could come back at 31% instead of 25%.

### Caveat 2 — with no significance filter at all, 8 of 9 exceed 30% global

Taken to the absurd limit — the best hit at *any* E-value out to 10000 — the PASS designs reach
**33.3-40.0% global identity** (only `id86` stays at 30.0%). The targets responsible are
compositionally biased, not homologous: **Carboxy-S-adenosyl-L-methionine synthase (CmoA)**,
**Golgin subfamily A member 8B** (a long low-complexity coiled-coil), **isoleucyl-tRNA
synthetase**, **Magnesium transport protein CorA**. Our designs are Ala/Glu/Lys/Leu-rich
idealized helices, which is why they align against charged coiled-coils at twilight-zone
identity. E-values run 0.6 to 2.4e4.

This matters only if Proteinbase applies **no** E-value threshold. No sane pipeline does, and
MMseqs2's own default is 1e-3 — six to seven orders of magnitude tighter than these hits. But
it is the honest worst case and it is the reason the margin is not described here as comfortable.
The single most significant of these is `egfr-siteA-id89` at **33.3% global vs `A4SPN7`
(CmoA)**, phmmer full E = **0.61**, domE = 1.6 — still insignificant, but closer to the line
than anything else in the cohort.

### Caveat 3 — these are not 9 independent measurements

Six of the nine PASS designs are **94.7-98.7% identical to each other** — they are 1-4 point
mutants of a single sequence:

```
id3035 <-> id3036  98.7%     id3020 <-> id3036  98.7%     id3020 <-> id2046  98.7%
id3017 <-> id3035  98.7%     id2015 <-> id2046  98.7%     id3017 <-> id3036  97.3%
id2015 <-> id3020  97.3%     id3035 <-> id2046  96.0%     id3017 <-> id3020  96.0%
id2015 <-> id3017  94.7%  ... (21 pairs at >= 94.7% within this cluster)
```

`egfr-siteA-id89` sits at 72-77% to that cluster; `cand-siteA-id2` at 60-63%; only
`egfr-siteA-id86` (80 aa) is genuinely separate, and it is 75.0% identical to
`cand-siteA-id73`. **The 9 PASS designs are effectively 2-3 distinct molecules.** Two
consequences for the submission decision:

1. Their novelty verdicts are **perfectly correlated** — they pass together and they would
   fail together. There is no diversification benefit from submitting 6 near-clones.
2. The challenge's own hard criterion asks for *"adequate sequence- and structural-diversity
   from known proteins"*. A slate of six 98%-identical sequences is a diversity exposure on a
   different axis from the one measured here, and nothing in this measurement addresses it.

### One status change vs the structural-only expectation

`cand-siteA-id19-sh09` is tagged **`FAIL_L2_both`**, i.e. expected to fail the sequence half
too. **It does not fail the sequence half on this measurement**: 0 significant MMseqs2 hits in
either DB, 0 domain-significant phmmer hits, best MMseqs2 E = 2.61 (Swiss-Prot) / 10.5 (PDB).
Its no-threshold worst case is high (**44.0-50.7% global** against `7obi_A` CTPR-rv4 and
`8qai_A`) and phmmer's repeat-accumulated full E reaches 3.1e-08 against an engineered TPR, so
this one is genuinely ambiguous rather than clean. It is not in the submission cohort, so the
stakes are low — but the `_both` half of its tag is not reproduced here.

---

## 5. The unexpected findings, and what could not be determined

**Unexpected**

- **The bridges fail on sequence, not just structure.** At 41-46% *global* identity to designed
  consensus TPRs at E ≈ 1e-30, they fail under every convention. The prior read was that the
  bridges' problem was the structural half; the sequence half is independently fatal.
- **The nearest neighbours are other people's *designed* proteins**, not natural ones — CTPR2/
  CTPR3/CTPR-rv4 for the bridges, a Baker-lab DHR for `id3017`/`id3035`, `C2-01018`/`9exk_A`
  "De novo designed protein K12"/DARPin for the site-A designs. Our designs are novel against
  nature and much less novel against the de novo design literature. If Proteinbase's databases
  keep growing with deposited designed proteins, this margin erodes over time.
- **phmmer is markedly more sensitive than MMseqs2 on repetitive targets** and disagreed on 3
  of 16 queries (`id19-sh09`, `id3017`, `id3035`). Running both was load-bearing, not
  decorative — a single-tool answer would have missed Caveat 1 entirely.
- The `-a` flag trap: a silent 0% identity for every hit, with exit code 0.

**Could not be determined**

- **Proteinbase's exact identity convention** (local vs global vs some coverage-weighted
  variant) and **its E-value / coverage thresholds.** This is the single largest residual
  uncertainty. On our numbers the PASS designs land at 25.3% (best significant-tier hit),
  30.7% (Caveat 1), or up to 40.0% (no filter at all) — the convention is worth up to ~15
  points and the gate is at 30.
- **What "moderate" vs "high" structural similarity means numerically** to Proteinbase. The
  brief's statement that all 16 queries have at least moderate structure was taken as given and
  not re-derived.
- Whether Proteinbase deduplicates near-identical submissions (Caveat 3) or scores each
  independently.

---

## 6. Limitations — databases searched vs not searched

Proteinbase runs MMseqs2 against **five** databases. Two were reproduced here.

| database | searched? | notes |
|---|---|---|
| **Swiss-Prot** | ✅ yes | release 2026_03, 575,748 seqs — full reproduction |
| **PDB sequences** | ✅ yes | 1,103,516 protein chains, `pdb_seqres` 25-Sep-2026 |
| **patent sequences** | ❌ **no** | not downloaded |
| **therapeutic antibody DB** | ❌ **no** | not downloaded |
| **PLAbDab** | ❌ **no** | not downloaded |

**Direction of the bias: every number in this report is a LOWER bound on the maximum identity
Proteinbase will find, and therefore an UPPER bound on novelty.** Adding databases can only
raise the maximum identity, never lower it. A design that passes here could still fail on a
database not searched here; a design that fails here cannot be rescued by one.

Assessment of the residual risk from the three unsearched databases:

- **Therapeutic antibody DB and PLAbDab are low risk.** Both are antibody-specific (paired
  V-region sequences). Our designs are non-antibody de novo helical bundles with no
  immunoglobulin fold, no CDR structure and no framework homology. A >30% global hit from
  either would be surprising.
- **Patent sequences are the real unsearched risk.** Patent corpora are large, redundant, and
  contain a great many *designed* and engineered proteins — including scaffold patents covering
  consensus TPRs, DARPins and designed helical repeat proteins. Given that our nearest
  neighbours in the two databases that *were* searched are already other people's designed
  proteins (CTPR, DHR, DARPin, "de novo designed protein K12"), the patent set is where an
  additional hit would be most likely. It cannot be bounded from this work.
- Also not searched: **TrEMBL** (unreviewed UniProt). Proteinbase specifies Swiss-Prot, so this
  is a match to their stated method rather than a gap — but it means a design could have a
  close unreviewed relative that neither they nor we would see.

Other limitations:

- Only the Swiss-Prot and PDB **sequence** sets were searched; no profile-vs-profile
  (HHblits) or structure-based search was run here — the structural half is a separate,
  already-completed measurement.
- `--max-seqs 5000` caps hits per query. The cap was reached for some queries, but it bounds
  only how many *insignificant* hits are listed, not whether the best hit was found (hits are
  returned in prefilter-score order).
- The `qcov >= 50%` floor on local identity is our choice, not Proteinbase's. Without it the
  local numbers are meaningless (100% over 9 residues); with a different floor they would shift.

---

## 7. Files

**Scripts** (`scripts/`)

| file | purpose |
|---|---|
| `novelty_seq_search_20261006.sh` | DB prep + MMseqs2 sensitive search, both DBs |
| `novelty_seq_phmmer_20261006.sh` | phmmer `--max` cross-check, both DBs |
| `novelty_seq_analyze_20261006.py` | local/global identity by E-value tier + alignment-string verification |
| `novelty_seq_crosscheck_20261006.py` | MMseqs2 vs phmmer significance + independent SW realignment |
| `novelty_seq_phmmer_identity_20261006.py` | identities for phmmer hits (phmmer emits none) |
| `novelty_seq_report_20261006.py` | generates the tables in section 3 |

**Data** (`data/novelty_20261006/`)

| file | contents |
|---|---|
| `mmseqs_swissprot_hits_20261006.tsv` | 9,474 hits, all columns incl. `nident`/`alnlen`/`qlen` |
| `mmseqs_pdbseqres_hits_20261006.tsv` | 8,567 hits |
| `phmmer_{swissprot,pdbseqres}_domtbl_E10_20261006.txt` | phmmer domain tables, E <= 10 rows only (3.7 / 5.0 MB) |
| `phmmer_{swissprot,pdbseqres}_domtbl_20261006.txt` | **full** phmmer tables, 25 MB each — **deliberately NOT deposited** (size). They are retained privately; the `_E10_` subsets above cover every significant hit, so nothing load-bearing is missing from this deposit. |
| `novelty_seq_summary_20261006.csv` | per-query best hit by E-value tier and convention |
| `novelty_seq_top5_20261006.tsv` | top-5 hits per query per DB, both identities |
| `novelty_seq_crosscheck_20261006.csv` | tool-vs-tool agreement + SW identities |
| `novelty_seq_phmmer_identity_20261006.csv` | phmmer hit identities by significance tier |

**Not deposited:** the search databases (Swiss-Prot + PDB FASTA) and the full hit tables *with* the
`qaln`/`taln` alignment strings and the phmmer raw output. The deposited copies drop only those two
bulky alignment-string columns; `nident` + `alnlen` + `qlen` fully determine both identity
conventions.

