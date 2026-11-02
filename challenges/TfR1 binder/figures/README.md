# Figures

Seven ray-traced renders of the scaffolded-design screen, 2400 × 1800 before cropping, one per
design shown plus an overview. They are views of the session file in `../structures/`, and they
exist because what they show is easier to see than to read: **how much of each binder the generator
built**, and whether the chain stays clear of the receptor.

⛔⛔ **These are VIEWS, not evidence.** Nothing in them ranks, scores or re-measures anything, and ⛔
**a maintained geometry is never evidence of affinity.** The clash figures stamped on them are raw
generator output, which clashes normally; read `../methods/core_scaffolded_design.md` §5 for what the
numbers mean and what they do not.

## What is here

| file | what it shows |
|---|---|
| `00_overview.png` | all six designs at once, against one copy of the receptor surface |
| `opt05_b2_m6.png` | best arm, best design — highest fold integrity, no residue inside 2.5 Å of the receptor |
| `opt05_b0_m0.png` | same arm, an **independent** batch |
| `opt05_b3_m4.png` | same arm, a third independent batch |
| `opt01_b1_m5.png` | runner-up arm; the best receptor clearance in the whole cohort, 3.44 Å |
| `opt01_b4_m5.png` | runner-up arm; the most compact survivor |
| `opt04_b4_m7.png` | third arm, for contrast |

**Six designs from six distinct batches across three arms**, so the set is not six views of one
contig draw. Batch membership matters here: one contig length is drawn per batch, so designs within
a batch are not independent (`../methods/core_scaffolded_design.md` §5.2).

## The colouring, which is the whole point

| colour | meaning |
|---|---|
| **purple** | the kept 5XOX core elements — what the generator was **given** |
| **grey** | what the generator **built** to connect them |
| **yellow** | the four fixed guidepost residues |
| **teal** | the receptor surface |

⭐ The grey is most of every chain: **85 to 98 residues of a 141–145 residue binder.** The held and
built sets are **measured geometrically** per design, not read off the input or the contig, and the
reason they have to be is in `../methods/core_scaffolded_design.md` §6.

## Each image is self-describing, deliberately

Every render carries a stamped **title band** above the image area holding the design identifier, the
colour legend, and that design's own measured numbers (largest Cα–Cα step, backbone minimum distance
to the receptor, minimum cut, radius-of-gyration ratio, residues within 2.5 Å). ⇒ a file that is
detached from this directory still says what it is.

⭐ **A label goes in the band and never into the image area.** An automatic labeller that strips a
common prefix and suffix cannot carry a constant marker, so the marker must live outside the picture.

## Three render defects were found and fixed rather than shipped

Recorded because the first one hid the only region the figure exists to show.

1. ⛔⛔ **The receptor surface was being clipped, and the clipped interior rendered as a black wedge
   sitting exactly over the binding interface.** Cause: the auto-orient command sets the near and far
   clipping planes tight to the selection, and a surface crossing the near plane is cut open. Fix:
   widen the slab explicitly after every orient, and set the ray-tracer's interior colour so any
   residual cut face reads as material rather than as a void. ⚠️ That colour is a **rendering choice,
   not data** — it changes how a cut face looks and nothing about geometry, and it is labelled as such
   in the program so nobody reads a grey interior as a measured surface.
2. ⛔ **The surface rendered fully opaque**, because transparency had been set **before** the objects
   were loaded, so it applied to nothing. Fix: set it per object after loading.
3. ⛔ **A stamped metrics line ran off the right edge on the narrower crops and truncated a column
   name.** ⇒ a **wrong** label, not a cosmetic issue. Fix: every title, subtitle, metrics line and
   legend row is now fitted to the available width, with a refusal below a minimum size, and a
   separate pass re-measures every stamped line against its canvas — **0 overflowing labels across
   all six designs.**

⚠️ **The crop can only remove border.** It is computed from the uniform white background alone and
applied with a fixed pad, so it cannot clip structure, and the before and after sizes are printed so
a bad crop would be visible rather than silent.

## How to regenerate them

`../scripts/j902_render_screen_pse_20261101.py` renders the session under a headless viewer, and
`../scripts/j902_stamp_renders_20261101.py` crops and stamps. ⛔ The render program has no
`__main__` guard on purpose: under that invocation the module is not `__main__`, so a guarded script
exits successfully having written nothing.
