# Structures — PyMOL sessions

Two sessions, built by `../scripts/s125_build_mb_depth_pse.py`.

| file | what is in it |
|---|---|
| `j1062_validated_cohort.pse` | the **submitted cohort's validated designs** beside the cross-paired non-binders that define the floor |
| `j1058_mb_depth_cohort.pse` | the full MB-only round, for context on what did **not** clear |

⭐ **The floor controls are in the session on purpose.** The interesting comparison is not a design
against empty space, it is a design against a real published binder for a *different* peptide–MHC
that was folded against this target. Opening the validated cohort shows both, in the same frame.

## How they are aligned, and why that is stated

Every object is superposed on the **HLA heavy chain**, so the target sits in one frame and the
binders can be compared directly. Two things about that are worth knowing:

- ⭐ The fit is gated on **atom count**, not RMSD: each superposition is required to land a large
  fraction of the 276 HLA residues, and in these sessions every fit landed **224–260 atoms of 276**.
  A low RMSD over a handful of atoms is the failure mode this check exists to catch.
- The designed binder is **chain D**, and its presence was asserted on every object in the session
  rather than assumed — a session can otherwise look complete while the thing it is about is absent.

⛔ **A disabled group hides its enabled members in PyMOL.** Groups and their members are enabled
explicitly in the builder and the visibility asserted, because sessions have previously shipped from
this group with the binders invisible.

## ⛔ What these sessions are not

⛔ **Every coordinate is a model, and the binder coordinates are predictions of a designed sequence,
not structures of a molecule anyone has made.** Nothing here was measured at a bench.

⛔ **A maintained geometry is never evidence of affinity.** A design can sit convincingly in the
groove and bind nothing. What these sessions are good for is seeing *where* a binder sits and
*whether it engages the peptide or only the shared MHC surface* — which is the whole question of this
challenge, and is easier to see than to read off a table.

⚠️ Raw generated backbones interpenetrate the target's solvent-excluded surface at the separations
used here; that is what an unrelaxed design at a low minimum-distance floor looks like, ⛔ not good
packing.
