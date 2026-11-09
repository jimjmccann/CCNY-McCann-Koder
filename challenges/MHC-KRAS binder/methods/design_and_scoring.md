# Design and scoring

Read `README.md` in this directory first for the target and the four design arms.

## 1 · What gets scored

A designed sequence is re-folded **as the full four-chain complex** with the target unchanged, and
then scored **only** on the part of the result that answers the actual question: does this molecule
engage *the peptide*?

Two quantities, both measured per fold:

- **`cont_pep`** — contacts between the binder and the **peptide** (chain C), computed over a
  consistent atom set across arms. ⛔ Not contacts to the complex, which the MHC surface dominates.
- **binder↔peptide interface confidence** — the **pairwise** ipTM for exactly that chain pair.
  ⛔⛔ **Not the global ipTM**, which on a four-chain target is dominated by the HLA/β2m interface
  that is identical in every design and tells you nothing about the binder.

## 2 · ⛔⛔ The floor is other people's binders, not a null model

The comparison set is **real, published binders for *other* peptide–MHC complexes, folded against
our target**. They are good molecules for this class of problem with no reason to recognise this
particular peptide. Measured, they set:

| | range across rounds |
|---|---|
| peptide contacts | **9.6 – 12.6** |
| binder↔peptide interface confidence | **0.861 – 0.895** |

⇒ ⭐ **The mean contact margin over the whole screen is negative at every sequence-design rank.** The
typical design touches the peptide *less* than one of these cross-paired non-binders does by
accident. Any screen that reported most designs as promising against this target would be measuring
something other than peptide engagement.

⭐ **The floor is replicated into every fold shard** and each shard is scored against **its own**
floor. The shards differ enough that pooling them would smear the comparison, and a design is never
compared to a floor measured in a different shard.

## 3 · ⭐⭐ The gate: a conjunction, pre-registered

A design clears if and only if **both** limbs clear its own shard's floor:

```
cont_pep > floor_cont_pep   AND   iptm(binder, peptide) > floor_iptm
```

⛔ **Never one limb alone. Never a ratio.**

- A **ratio** or a weighted sum pays partial credit, so a molecule that is excellent on one limb and
  fails the other can outrank one that clears both. Averaging fractions hides exactly the failure the
  conjunction exists to catch.
- **Contacts alone** reward a binder draped across the groove without a confident interface.
- **Interface confidence alone** rewards a confident grip on the shared MHC surface.

### The conjunction earned its keep, on the design we had publicly called best

The design with the **highest peptide contacts in the whole cohort** cleared the contact limb at all
four seeds — 15.6, 16.8, 17.2, 17.4 — and failed **only** the interface limb, at a single seed, with
a value **inside the non-binder band**.

⇒ ⭐⭐ **Ranked on contacts alone it would have been the top design of the workstream.** It had been
named in writing, before the data existed, as the strongest object in its round and the first thing
to re-test. It did not survive. Plenty of peptide contact; an interface no more confident than a
cross-paired non-binder's.

⚠️ The direction of this is not a constant. In one earlier round the interface limb rejected most
contact-clearing designs; in another it was the permissive limb. ⇒ the interface metric is recorded
as **stringent-or-confounded depending on the comparison**, ⛔ not as a property of the metric.

## 4 · Gates applied before anything is scored

In order, and each one refuses rather than down-weights:

1. ⛔⛔ **Backbone continuity first.** A clash-based floor *rewards* a broken chain, because breaking
   the backbone is the cheapest way to relieve a clash. Continuity is therefore gated **before** any
   clash metric. `[M]` In this campaign continuity passed 160/160 in both independent rounds — it has
   never been the binding constraint here, which is worth stating so nobody adds a gate that does no
   work.
2. **Raw-clash survival**, on unrelaxed generator output. ⛔ Raw backbones clash normally; a tight
   clash gate applied to them measures the generator's output convention, not design quality.
   Measured survival replicated across rounds at **90.6%** and **91.9%**.
3. ⛔ **No cysteines**, enforced at sequence design (`--omit_AA C`) and re-asserted when the
   submission file is built.
4. **Chain-length and composition checks** against the target's known chemistry — **276 / 99 / 10**
   in the design-stage frame — so a mis-assembled input is caught rather than folded. ⚠️ Those are
   the design-stage numbers; the co-folded complex carries β2m as 100. See the two-convention table
   in `README.md`, which is a live trap and not a footnote.

## 5 · What a cleared design is, and is not

⛔ **Clearing this gate makes a design a candidate. It is not evidence of affinity, and it is not a
result.** The gate is a comparison between two sets of model outputs. Nothing in this directory was
measured at a bench, and a maintained geometry in a figure is never evidence of binding.

See `robustness_and_tiering.md` for the step that separates a one-fold clear from a claim.
