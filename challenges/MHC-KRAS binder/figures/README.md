# Challenge 5 — submission renderings, 7 designs, one colour each

Nine ray-traced PNGs of the Challenge-5 submission set: **one panel per design**, an
`00_overview.png` holding all seven against one copy of the target, and a combined figure. They
are views of `../structures/j1063_c5_submission_7designs_20261109.pse`, and they exist because what they show is
easier to see than to read: **where each binder sits relative to the KRAS G12D peptide**, and how
differently the seven do it.

⛔⛔ **These are VIEWS, not evidence.** Nothing in them ranks, scores, re-folds or re-measures
anything, and ⛔ **a maintained geometry is never evidence of affinity.** Every number stamped on
a panel is *carried* from the scored tables (`t152` / `t154` / `t155`) and is re-asserted against
them cell by cell at build time — 28 cells, 0 disagreements. Nothing here has been expressed or
assayed.

---

## ⛔⛔ The 4 / 3 split, which is the thing this figure must not let you misread

James McCann put **all seven** designs in the submission on 2026-11-09. **That is a scope decision and it
does not change what any of them is.** Seven panels side by side is exactly the image that could
be taken for seven equally supported designs, so every panel states its own tier on its own face
and the combined figure says it in the header.

| tier | n | what it means |
|---|---|---|
| **VALIDATED** | **4** | clears the pre-registered conjunction at **4 of 4 seeds over 20 draws** |
| **3 of 4 SEEDS — INSIDE NOISE, NOT PROMOTED** | **3** | clears at 3 of 4 seeds; the *J1058* §6 tier |

⛔ The word **CANDIDATE** is deliberately **not used anywhere** in these figures. In that
pre-registered scheme it means a *5-draw-only* design, and none of these are — all seven have been
to 20 draws. ⛔ None of the three is called validated, promoted, or a result, on any panel, in any
filename, or here.

⭐ Every demoted design still cleared **3 of 4** seeds and none collapsed to 0 of 4, so the gate
separates **robust from marginal**, not real from fake.

---

## What is here

| file | what it shows |
|---|---|
| `00_overview.png` | all seven binders at once against **one** copy of the target — ⭐ this is the legend for the set |
| `01_arm_shard03_MB_nomotif_broad_2_model_3.png` | **VALIDATED** · rank 1 · margin **+6.0** · ipTM 0.958 · 82 aa |
| `02_arm_shard01_MB_nomotif_broad_11_model_1.png` | **VALIDATED** · rank 2 · margin **+5.2** · ipTM 0.937 · 75 aa |
| `03_arm_shard02_M_nomotif_narrow_8_model_2.png` | **VALIDATED** · rank 5 · margin **+3.4** · ipTM 0.931 · 81 aa |
| `04_arm_shard01_MB_nomotif_broad_6_model_3.png` | **VALIDATED** · rank 4 · margin **+1.1** · ipTM 0.931 · 82 aa |
| `05_arm_shard00_MB_nomotif_broad_11_model_2.png` | **3 of 4 seeds** · rank 3 · margin +6.4 · ipTM 0.894 · 66 aa · the **ipTM** limb failed at one seed |
| `06_arm_shard03_MB_nomotif_broad_19_model_0.png` | **3 of 4 seeds** · rank 2 · margin +2.7 · ipTM 0.907 · 76 aa · the **ipTM** limb failed at one seed |
| `07_arm_shard03_MB_nomotif_broad_4_model_0.png` | **3 of 4 seeds** · rank 1 · margin +0.9 · ipTM 0.916 · 72 aa · the **peptide-contact** limb failed at one seed |
| `j1063_c5_submission_7designs_combined_20261109.png` | all eight panels in one 4 × 2 figure, each keeping its own band |

Deposited alongside these images, and all regenerable from `../scripts/`: the session
`../structures/j1063_c5_submission_7designs_20261109.pse`, its build manifest
`../data/j1063_c5_submission_7designs_20261109_build.json` (which structure file each panel came
from, and the alignment atom counts), and `../data/j1063_labels_20261109.json` (every stamped line,
for the re-measurement pass).

⚠️ **The MPNN rank is not 1 for most of them, and that is correct.** The re-draw re-folds *the
sequence that actually cleared*, which is that design's own best-margin rank — rank 5 for slot 3,
rank 4 for slot 4. Ranking the rows instead would have re-drawn a different molecule.

---

## The colouring, which is what was asked for

| colour | meaning |
|---|---|
| **grey** | HLA-A\*11:01 (chain A, 276 res) + beta-2 microglobulin (chain B, 100 res) — **one shared copy**, 25 % cartoon transparency |
| **CCNY purple #7D55C7** | the **KRAS G12D 10-mer `VVVGADGVGK`** (chain C), cartoon **and** sticks, fully opaque |
| **one colour per design** | each binder (chain D), the same colour in its own panel and in the overview |

The seven design colours are the **Okabe–Ito colour-blind-safe qualitative set**, chosen so the
figure survives deuteranopia and protanopia and does not rely on red/green contrast:

| slot | colour | hex |
|---|---|---|
| 1 | blue | `#0072B2` |
| 2 | bluish green | `#009E73` |
| 3 | orange | `#E69F00` |
| 4 | sky blue | `#56B4E9` |
| 5 | vermillion | `#D55E00` |
| 6 | reddish purple | `#CC79A7` |
| 7 | black | `#000000` |

⚠️ **One deliberate departure from the Challenge-2 convention, stated because it is a real
difference.** `j530`/`j531` reserve CCNY purple **#7D55C7 for the rank-1 design**. Here purple is
the **peptide** instead, following the newest Challenge-5 instrument (`s94`, 2026-11-06, whose
legend already keys purple to the KRAS G12D peptide). The reason is that **peptide contact is the
entire basis of the pre-registered gate**, so the peptide has to be the visually obvious feature,
and the lab accent is the strongest thing in the palette. The build program **refuses** to assign
`#7D55C7` to any design, so the two cannot collide by accident.

⚠️ Okabe–Ito's yellow `#F0E442` is too pale to read as a cartoon on white, so **black** takes that
slot — black is part of the published Okabe–Ito set. It renders flat (no shading), which costs
some depth cue on slot 7 while keeping it unmistakable.

---

## Each image is self-describing, deliberately

Every panel carries a stamped **title band above the image area** holding the slot, the full design
name, **the tier sentence**, the public submission name, the generation round, arm, MPNN rank,
binder length, that design's own carried numbers, the fit it was placed by, the colour legend, and
two caveat lines.

⭐ **A label goes in the band and never into the image area.** A PyMOL per-cell labeller strips
whatever is *common* across the label set — `[M]` three times on 2026-10-19, `DES_s07_3_s05` drew
as `s07_3_s05` — so **a constant marker can never survive inside a per-cell label**, and the tier
word is exactly such a constant.

⚠️ **One target copy, and the price of it.** The seven are independent folds, so their HLA frames
differ. Each is superposed **HLA-on-HLA** onto the slot-1 reference and only the reference's target
is kept, so the target is **literally the same object** in every panel. The price is that each
binder's drawn pose carries that design's own fit residual, so **the fit rmsd and the atom count
are stamped on every panel**: 0.168–0.251 Å over 227–258 of 276 CA.

⚠️ **The shown draw is `model_0` of fold shard 01 — that seed's own top-ranked model of 5, not a
random draw**, and it is labelled as such on every panel. One model is not a distribution; the
distribution is the 20-draw conjunction the tier sentence refers to.

---

## Four render defects were found and fixed rather than shipped

1. ⛔⛔ **The overview put the peptide underneath the binders and the figure contradicted its own
   caption.** `cmd.orient` on peptide + designs chose a view looking roughly **down** onto the HLA
   groove, so all seven binders stacked directly on top of chain C and the KRAS G12D peptide — the
   feature the whole gate is about — showed as a few purple flecks. ⭐ **Found by opening the raw
   PNG. Every automated check passed**: object counts, group-enable assertions, per-panel atom
   counts and the ink fraction. Same class as the TfR1 N-terminal figure, where the built segment
   was the least visible thing in its own panel. ⇒ the camera is no longer left to `orient`: a
   frame is **derived from the geometry** so the groove is edge-on, with the peptide's own long
   axis on screen-x and the binders above it by construction. ⛔ And the **coordinates** are
   rotated, not the camera — PyMOL's 9-float view rotation is the transpose of the
   `[right, up, back]` convention and a hand-built camera matrix silently produced sideways figures
   in `j695`. Rotating coordinates can be **re-measured** afterwards, and it is: peptide spread
   22.2 / 8.1 / 3.4 Å, binder centroid y +24.1 against peptide centroid y +13.8, both asserted.
2. ⛔ **A 38-character design name is one unsplittable token and it stopped the run.** The honest
   crop of these renders is 1097 px wide, and `arm_shard03_MB_nomotif_broad_2_model_3` does not fit
   on one line even at the 22 px refusal floor — so `fit()` **refused**, correctly. Shrinking
   further would have produced a label nobody can read on a printed panel. ⇒ lines **wrap**, and
   the font shrinks until the longest single token fits before wrapping. ⛔ Nothing is ever
   truncated and nothing is drawn below the floor.
3. ⛔ **The band was nearly as tall as the picture.** At 3000 × 2250 landscape this complex cropped
   to 1097 × 1525 — it used barely a third of the frame width — and the band wrapped to 1100–1470 px
   against a 1525 px image. Every label passed its re-measurement and the figure was still mostly
   text. ⇒ the canvas is **portrait, matched to the subject's own 0.72 aspect**, and the type
   ratios follow the deposited precedent rather than being retuned by eye. The band is now
   469–596 px against a 2142 px image.
4. ⛔ **A nine-key legend cannot sit on one row.** The overview keys the target, the peptide and
   all seven designs. A text-only fit is how a six-key row once ran 258 px off the right edge once
   the chips, gaps and separators were laid down. ⇒ the legend fit **includes the chips and the
   gaps**, and the row **wraps** over as many rows as it needs rather than overflowing. It uses
   three.

⚠️ **The crop can only remove border.** The box is the union of the white-background content boxes
of all eight renders, so it cannot clip structure in any of them, and it is the same box for every
panel so the panels stay pixel-comparable. Narrow crops are then **padded** out to 1900 px and
centred — that adds border, it never removes structure — and the before/after sizes are printed.

**Verification, printed rather than assumed:** **0 overflowing labels of 104**, re-measured on the
*finished* PNGs by a separate pass that re-opens each file, re-creates each recorded font at its
recorded size and re-measures each recorded string; **0 images with ink touching a left or right
edge**, checked independently of the text metrics; and **0 of 4** header lines overflowing on the
combined figure.

---

## ⛔⛔ And the check no checker can do

Every defect above except the first was found by a program. **The first was found by opening the
image**, and it is the one that mattered: a label checker measures labels, and it cannot measure
whether the picture shows what the caption claims. All nine finished PNGs were opened and read
individually before this file was written. In each one the binder is plainly visible, in its own
colour, above the groove, and the purple peptide is visible in the groove beneath it.

⚠️ **What the overview cannot tell you:** the seven binders overlap heavily in one frame, so slot 5
(vermillion) is largely hidden behind the others. The overview is a **legend and a spread**, not a
per-design view; the per-design panels are.

---

## How to regenerate them

Four stages — build the session, render, stamp, combine — following the shape this group has used
for submission panels in earlier challenges. All four scripts are deposited in `../scripts/`, and
the inputs they take are in `../data/`. **Local, CPU only, no GPU and no rented compute.**

```
pymol  -cq j1063_build_submission_pse_20261109.py     -- j1063_picks_20261109.json <out.pse>
pymol  -cq j1063_render_submission_panels_20261109.py -- <out>_build.json "$LAB_SCRATCH/c5_j1063_raw_20261109"
python3 -I j1063_stamp_panels_20261109.py                "$LAB_SCRATCH/c5_j1063_raw_20261109" <this dir>
python3 -I j1063_combine_panels_20261109.py              <this dir> <this dir>/..._combined_20261109.png
```

⛔ The two PyMOL programs have **no `__main__` guard** on purpose: under `pymol -cq` the module is
not `__main__`, so a guarded script exits 0 having written nothing. ⛔ They end in `os._exit(0)` to
skip PyMOL's crash-prone teardown — **and flush stdout first**, because `os._exit` skips Python's
flush and on a pipe the whole log is otherwise lost while the run "succeeds" silently (measured
here: a probe printed nothing at all until the flush was added).

⛔ Raw renders go to a large scratch volume (`$LAB_SCRATCH` in our environment), not the root drive —
a full-resolution render set will fill a system disk and take running jobs down with it.

`../data/j1063_picks_20261109.json` is the input and the only place a number is written by hand; stage 1
asserts all of it against the scored tables and **refuses** on any disagreement.

---

## Other guards that did not fire, recorded because they are the ones that have cost real work

- ⛔⛔ **A disabled group hides its enabled members** — twelve `.pse` once shipped with no binders
  while objects 9/9 and scenes 6/6 both passed. Members **and** groups are enabled explicitly and
  asserted via `cmd.get_type(n) == "object:group"`, and the check **refuses when it finds no
  groups at all** so it cannot pass vacuously, which is how the first fix for this defect reported
  success while doing nothing.
- ⛔⛔ **`retain_order 1` corrupts `cmd.align`** — a 276-residue fit once collapsed to 14 atoms at a
  *numerically better* rms than the correct 218-atom fit. `retain_order` is left **unset** and
  every fit is gated on **atom count**, never rms: floor 200 of 276, achieved 227–258.
- ⛔⛔ **A non-ASCII character in `cmd.scene(message=)` segfaults the saved `.pse`** — every message
  goes through an ASCII strip.
- ⛔⛔ **Design basenames are not unique across rounds.** `MB_nomotif_broad_11_model_1` exists under
  **both** `arm_shard03_` and `arm_shard01_` — different backbones, same stem — and `19_model_0`
  and `4_model_0` each exist under three prefixes, so matching by stem can silently load the wrong
  molecule. Every structure is matched on its **full** name, and the **binder length is asserted
  per design** (82 / 75 / 81 / 82 / 66 / 76 / 72 aa) as an independent handle on identity. A
  mismatch is a stop, not something to work around.
- ⛔ **A blank render exits 0 and writes a valid PNG.** Each raw panel's non-white **ink fraction**
  is counted on the finished file before it is accepted (0.070–0.151 here).
- ⚠️ `ray_interior_color grey70` is a **rendering choice, not data**: it changes how a cut face
  looks and nothing about geometry. The slab is widened explicitly after every camera change, so
  no surface or cartoon is cut open across the interface.

---

## What the numbers mean, and what they do not

**Contact margin** is peptide contacts above that design's **own shard floor**, and the floor is
set by **cross-paired non-binders** — real published binders for *other* pMHCs, folded against this
target without being designed for it. ⛔ **That floor is brutal**: they make 9.6–12.6 peptide
contacts and 0.861–0.895 binder–peptide ipTM by accident, and the mean contact margin is
**negative at every rank**, i.e. the typical design touches the peptide *less* than a non-binder
does. Clearing it makes a design worth re-drawing; it does not make it a binder.

⚠️ The three demoted designs are not demoted for being obviously worse. Slot 5 has the **highest**
mean contact margin in the whole set (+6.4) and the **lowest** ipTM (0.894, inside the non-binder
floor band): plenty of peptide contacts, an interface no more confident than a cross-paired
non-binder's. ⛔ Ranking on contacts alone would have put it top of the workstream. That is the
whole reason the gate is a conjunction, and it is visible in these panels: slot 5's binder sits
close against the peptide and still did not hold.
