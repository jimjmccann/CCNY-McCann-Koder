# How the 20 submitted designs were selected

This note records the selection procedure for the 20 sequences in `selected_20_with_scores.tsv`.
It is a description of what was done, not a claim that any of these designs binds: none has been
assayed, and every number here is computed.

## The starting cohort

Backbones were generated against the GLP-1 receptor ectodomain and the receptor core, using
guidepost residues taken from a GLP-1 co-fold peptide. Nineteen generation arms were run, differing
in which guideposts were used, how widely they were spread, and the direction in which the chain
traverses them. Sequences were designed onto the resulting backbones with a soluble variant of
ProteinMPNN, cysteines excluded, and screened for N-linked glycosylation sequons.

Each sequence was co-folded with the full 409-residue receptor in five independent runs. The five
runs per sequence are the unit of evidence: re-fold scatter is large on this target (median
standard deviation 0.142 against a median value of 0.352), so a single prediction is not
informative, and 45 per cent of sequences disagree with themselves across their own five runs. Only
sequences whose five runs agree to within 0.15 were considered. Collapsing to one representative
sequence per backbone gives the 160 backbones scored in `scoring/`.

## The axes, and what each is worth

Five computed axes were used. They are not equally trustworthy, and the weighting reflects that.

- **Interface contact count and buried surface area per residue.** These are the only axes with
  external evidence of correlating with measured binding in published work on comparable designs.
  Buried area is expressed per residue so that a longer binder does not score well on size alone.
- **Interface predicted aligned error.** Lower is better. This cohort is poor on this axis
  throughout, with a median of 18.0 A and a best decile of 11.8 A, so it separates designs from one
  another but certifies none of them.
- **Absolute Gaussian linking number**, a measure of topological entanglement between the designed
  chain and the receptor. Two experimentally validated de novo binders of a related class-B1
  receptor score 0.0941 and 0.0995, so 0.0941 is used as a reference line rather than a tuned
  threshold.
- **Co-fold confidence**, taken as the lower of the two directions of the pairwise interface
  confidence rather than the global value. This axis was deliberately given the least weight: on an
  earlier set of 248 designs with measured binding outcomes, the same confidence panel achieved a
  best area under the ROC curve of 0.601 with a confidence interval spanning 0.5, so it is not a
  demonstrated predictor of binding.

Designs were ordered on the **minimum** of their percentile ranks across these axes rather than the
mean. An average rewards a design that is excellent on one axis and disqualifying on another; the
minimum requires a design to be acceptable on all of them.

## Expert visual review, and the screen derived from it

The 40 highest-ranked designs were inspected individually. Eleven were judged to engage the
receptor at the residues the designs were built to reach; twenty-nine were not, including two that
scored well but were judged on inspection to be unsuitable. Those eleven occupy slots 1 to 11 of
the submitted set and were not re-scored, because a computed proxy should not override a direct
assessment of the property it is a proxy for.

The contrast between the two groups was then used to derive a screen. Sixteen receptor residues are
contacted by most of the designs judged suitable and few of the others; designs judged suitable
contact a mean of 13.7 of these sixteen, against 5.3 for the rest. The set includes the receptor
positions the guideposts were constructed against, which it recovered without being given them. The
distance from the designed chain's C-terminus to this site also separates the two groups
(26.9 A against 32.9 A), while the same measurement taken from the N-terminus does not (29.0 A
against 27.5 A), indicating a directional effect rather than an artefact of chain length.

⚠️ This screen rests on eleven positively labelled examples. A single cut at twelve contacted
residues reproduces 88 per cent of the visual assessment on the designs it was derived from, but
only 72 per cent under leave-one-out cross-validation. It is treated here as a hypothesis used to
prioritise, not as a validated filter.

The remaining nine slots were filled from the 120 backbones not inspected visually, requiring at
least twelve contacted site residues — the floor of the visually approved group — and applying the
conjunction described above. No more than three designs come from any one generation arm, because
the arms target partly disjoint surfaces and twenty designs from one arm would constitute a single
hypothesis tested twenty times.

## A tension that is recorded rather than resolved

Topological cleanliness and target engagement are anti-correlated in this cohort, to the point of
being close to mutually exclusive. Every design in the uninspected pool that beats the 0.0941
reference line contacts **none** of the sixteen site residues, and the best linking number among
designs that do engage the site is 0.272. The submitted set therefore prioritises engagement, and
nine of its twenty members carry linking numbers between 0.33 and 0.46. This is a property of the
generation method, not of the ranking, and it is the clearest target for improvement in the next
round.

## What is not claimed

No design here has been experimentally tested. The receptor conformation is itself predicted and
varies between runs, by up to 9.5 A between independent predictions of the same receptor, so
comparisons of where different designs sit on the receptor are approximate. Binding and receptor
activation are separate properties and nothing in this selection addresses activation.
