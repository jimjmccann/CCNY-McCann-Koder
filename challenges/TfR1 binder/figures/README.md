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

---

# `nterm_top3_corrected_20261102.png` — the extension arm's corrected top three

Three panels plus an overview, built from
`../structures/j956_nterm_top3_corrected_20261102.pse` and documented in
`../methods/nterm_extension_arm.md`. Colouring: **purple** the whole protein given to the generator
and held sequence-fixed, **grey** the N-terminal extension the generator built, **yellow** the held
guidepost shown as sticks, **teal** the receptor as a surface. The receptor is the same object, in
the same frame, in every panel — asserted, not assumed.

⚠️ **The figure states on its own face that ranks 1 and 3 come from the same generation batch**, so
the three panels are **two** independent draws and not three.

## ⛔⛔ A fourth render defect, and the only one no check could catch

The three defects above were each found by a program. This one was found by **opening the image**.

At the surface transparency the earlier figures used, the **grey built segment — the entire subject
of the figure — was the least visible thing in it.** It read as a few faint patches behind the
receptor surface, which a reader would reasonably take as *the extension does not engage*, while the
measurement printed in the same panel says **141 to 165 contact pairs** under 4.5 Å. ⇒ a figure that
contradicted its own caption, with **0 overflowing labels** and every automated check passing.

⭐ Fix: the transparency is now raised through an **additive** override that leaves the program's
default unchanged, and the segment is plainly visible threading into the receptor. The original
renders are kept, not replaced.

⭐⭐ **And raising it exposed something worth stating rather than hiding:** the built segment sits
**inside** the receptor's solvent-excluded surface, which is what a **2.43–2.55 Å** minimum
separation looks like (`[M]` over the three panels; the figure's own caption brackets it more
loosely as 2.1–2.6 Å). These are raw generator backbones gated at a 2.0 Å floor with no relaxation.
⇒ a third caption line now says so **on the figure itself**, so a reader cannot mistake
interpenetration for good packing. Re-verified after the change: **0 overflowing labels of 47**, and
**0** images with ink touching an edge.

⚠️ ⇒ **The lesson generalises past this figure:** a label checker measures labels. It cannot measure
whether the picture shows what the caption claims, and that is the failure mode that matters most in
a published figure.
