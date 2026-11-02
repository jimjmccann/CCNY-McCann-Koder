# Transferable findings from a structural-mining line that closed

McCann / Koder groups, The City College of New York. Written 2026-10-31; **§12–§18 and the additions
marked 2026-11-01 were appended in a second pass** covering the scaffolded-design line that followed
the mining line.

⚠️ **The new findings are appended, not interleaved.** §1–§11 are ordered by how transferable they
are; §12 onward are not inserted into that order, because renumbering would invalidate every
citation of this file. Read §12–§18 as a second list of the same kind.

This is the part of the TfR1 work most likely to be useful to someone else. It is not about TfR1,
and none of it depends on the particular binding site involved. It is a list of things that went
wrong, or that had to be built to stop something going wrong, in a two-day effort to search large
structure databases for chains that already present a designed arrangement of side-chain functional
groups.

**What the line was:** define a binding site as a small set of typed functional groups at fixed
relative positions; search two structure databases for chains whose own functional groups reproduce
that arrangement within an RMSD bar; then test whether any match could actually be used as a
scaffold. The method and all the numbers are in
`methods/structural_mining.md`. **Net yield of usable scaffolds: zero.**

The findings below are ordered by how transferable they are, not by when they happened.

---

## 1. ⛔⛔ A distance pre-filter must be admissible against the criterion it pre-filters for, and the floor is computable

**The error.** Candidate matches were enumerated under a pairwise-distance pre-filter: a pair of
database points may serve a pair of query points only if their separations agree to within `tol_d`.
Acceptance was then an RMSD test, `RMSD ≤ tol_rmsd`. `tol_d` was set to 1.5 Å and `tol_rmsd` to
2.0 Å.

⛔ **A pre-filter tighter than its own criterion discards candidates the criterion would accept.**
The published result — "the six-group arrangement is matched by 0 chains, 0.00%, the geometry is
unsatisfiable" — was a statement about the filter. At an admissible tolerance roughly **1%** of
chains in each of two databases present it.

**The floor, derived.** For a Kabsch superposition of *k* points with RMSD ≤ *r*, the squared
per-point deviations satisfy `Σδ² ≤ k·r²`. A pairwise distance can deviate by at most `δ_i + δ_j`.
Maximising `δ_i + δ_j` subject to `δ_i² + δ_j² ≤ k·r²` gives `δ_i = δ_j = r·√(k/2)`, so

```
max | Δd_ij |  =  2 · r · √(k/2)  =  √(2k) · r
```

⇒ **`tol_d ≥ √(2k) · tol_rmsd`.** At `tol_rmsd` = 2.0 Å that is 4.899 / 5.657 / 6.325 / 6.928 Å at
*k* = 3 / 4 / 5 / 6.

⛔ **It is per-*k*, never one constant across *k*.** A single pre-filter constant applied to a sweep
over *k* is wrong at every *k* but one.

### 1.1 ⭐⭐ The invariant that tells you which of your conclusions are at risk, before you re-run anything

Because `tol_d` gates only a pre-filter, raising it can **only admit** candidates previously
discarded. It can never reject one it previously kept.

⇒ ⭐⭐ **Hit counts are monotone non-decreasing in the pre-filter tolerance. Any count measured at a
too-tight tolerance is a LOWER BOUND.**

That single fact triages an entire damaged result set for free:

- conclusions of the form **"this is common"** are safe by construction — they can only get stronger;
- conclusions of the form **"this is absent / rare / unsatisfiable"** are the ones that can reverse.

In this case the two low-*k* rows of a four-row curve needed no re-litigation at all, and the two
high-*k* rows reversed. Both predictions were written down before the corrected run, and the
monotonicity check — run on every one of 42 subsets rather than on the best one — came back with 0
violations, which is also how you know the new code's wiring is sound.

⇒ ⭐ **Make monotonicity in your own parameters an assertion, not an expectation.** It costs nothing
and it distinguishes "the result changed because the parameter was wrong" from "the result changed
because the new code is wrong".

### 1.2 ⛔ Widening an inadmissible filter adds loose matches; it does not produce a better one

Worth knowing before celebrating a recovered result. Per *k*, the **best** RMSD achieved was
identical at the inadmissible and admissible tolerances. The recovered population sat entirely
between 1.07 Å and the 2.0 Å bar, and vanished completely at a 1.0 Å bar.

⇒ the corrected finding was *"this is satisfiable, but only loosely"* — a different and much weaker
statement than *"this is satisfiable"*, and the difference is the whole practical question.

⚠️ At the level of an individual cell the floor **can** be better — it was in 3 of 84 — because at a
too-tight `tol_d` a match is only *guaranteed* to pass if `r ≤ tol_d/√(2k)`. Never worse. State such
a claim at the granularity you measured it.

---

## 2. ⛔⛔ Publish the identities, not just the counts. A statistic is not a result if nobody can look at a single instance

Two days of searching produced curves whose output columns were
`k / subset / tolerance / hits / chains / percentage / best_rmsd` — and **no chain identity.** Not
one of the ~1,300 matching chains was known by name.

The identities had been **computed and thrown away**: the matcher already returned a chain index per
hit, and a chain-name file sat beside each database. The driver summed them into a count.

**Emitting the named list took eleven seconds of local CPU once it was asked for.** It immediately
produced the two most important findings of the whole line (§4), neither of which is visible in a
count.

⇒ ⭐ **If a search's purpose is to find candidates, its output row is a candidate, not a tally.**
Aggregate afterwards; you can always sum a list of names, and you can never un-sum a count.

⇒ ⭐ A cheap diagnostic for this failure mode: ask what the **next action** on the output would be.
If the answer is "run something else to find out which ones they are", the output is incomplete.

---

## 3. ⛔⛔ A control that shares a parameter with the anomaly it explains proves nothing

This is the most general finding here, and it appeared twice in two different forms.

### 3.1 The rescue and the thing it rescued were one artefact

The "0 chains match" result of §1 had been explained away by a positive control: take a **real**
protein's own six functional groups as the query, and it *also* returns 0 hits in 20,000 other
chains — therefore the zero is a property of the regime, not of our design, and "such an arrangement
is essentially unique in nature."

⛔ **That control ran at the same inadmissible tolerance as the zero it was written to explain.** At
an admissible tolerance, four real control arrangements return 286 / 14 / 1,242 / 156 matching
chains. The uniqueness conclusion is refuted, and the honest version is *"our filter hid these
arrangements — for our design and for nature alike."*

⇒ ⛔⛔ **When an anomaly and the control that normalises it share a parameter, the control is
invariant to that parameter and therefore silent about it.** List what your control shares with your
measurement before you trust it. Everything shared is outside its reach.

### 3.2 ⛔ A cross-database replication cannot validate a parameter both arms use

The same curve was replicated in a second, independently sourced database: same best subset at 4 of 4
point counts, rank correlation 0.9904 across all 42 subsets, percentages within 0.88–1.01× at every
point count. That replication was real, and it was **never** evidence that the tolerance was right —
both arms shared the bad parameter, so the agreement was invariant to it.

⇒ ⭐ **Replication tests whether a method is consistent. Algebra tests whether its parameters are
admissible.** They are different checks and one does not substitute for the other. What caught the
error here was a derivation, not a replication.

### 3.3 ⚠️ And independence is usually weaker than the two database names suggest

One store was PDB-derived and carried crystallographic redundancy — homologues, re-depositions,
multiple chains of one assembly. The other was predicted monomers built from sequences, many of which
have PDB relatives. ⇒ the two sample **overlapping biology by different routes**.

⇒ ⭐ That is a **strong cross-check on the method** and a **weak one on nature's frequencies.** ⛔ A
percentage from either store is not a frequency in nature, and a high correlation between them is not
evidence that it is. What would make it one is a store that is *compositionally* different, not merely
differently sourced.

---

## 4. ⛔⛔ Selecting on one axis can anti-select another, and two oppositions can close a line

The named hit list produced two findings of the same shape, and a third measurement finished the job.

**Size against fit.** Of 795 PDB-derived chains presenting the six groups, **only 2 were under 120
residues** and 91% were over 250; across both databases the shortest chain presenting the
constellation at all was **116 residues**. And the small ones were the **worst** fits: best match
at ≤ 250 residues 1.4662 Å, against the list's overall best of 1.0680 Å.

⛔ **A measurement lesson sits inside that one.** The first pass took chain length to be the highest
residue *number* among a chain's matched groups. That is not a length and not a bound on one: a
domain numbered 296–411 reads as 407 and is 116 residues long, so the error has no fixed sign. The
conclusion survived recomputation from true residue counts, but a specific published claim — *"not
one is under 120"* — did not. ⇒ **Never let a residue identifier stand in for a residue count.**

**Accessibility against fit.** Of 1,314 matching chains, only **31 (2.4%)** presented all six groups
on a face a partner could reach from one direction; **1,140 (86.8%)** presented them where nothing
could reach them. And again the best fits were exactly the ones that failed: best RMSD among failures
**1.0680 Å**, best among survivors **1.5400 Å**. The entire tight head of the ranked list was buried.

⇒ ⭐⭐ **Both oppositions have one reading:** a wide arrangement of several charged and aromatic
groups is, in nature, usually an **interior or active-site** arrangement rather than a surface patch.
Rank order on geometric fit was therefore an *inverse* proxy for usability.

⇒ ⭐ **The transferable habit: whenever a ranked list is produced, measure the property you actually
need against the rank.** Here, two independent usability properties were both negatively correlated
with the ranking score, which no amount of looking at the top of the list would have revealed.

---

## 5. ⛔⛔ "Can it be reached?" is not "does the host body stay out of the way?" — and the second is strictly stronger

The accessibility gate above cast rays outward from each matched group and asked whether a partner
could reach all six from one direction. 31 chains passed.

Then the survivors were **actually placed** on the arrangement in the target's own coordinate frame.
Each of the four distinct small candidates put **679 / 1,271 / 2,469 / 1,694** heavy atoms within
2.5 Å of the target, with closest approaches of **0.40 / 0.15 / 0.15 / 0.07 Å**. Not grazing
contacts: the host protein's body sitting where the target is.

⇒ ⛔⛔ **A reachability test is about lines of sight from the site. A clearance test is about the
volume the host occupies.** The first does not imply the second, and in this case **all four
survivors failed the second.** The gate had been designed to the weaker requirement, and the gap was
not noticed until the placements were actually built and looked at.

⚠️ **The warning had been visible for a day and was read as a nuisance instead of a requirement.** A
clash filter had been retired earlier for having no discriminating power, and the single number that
explained its behaviour was that the arrangement's centroid sits **1.91 Å from the nearest target
heavy atom** — *inside* the target's van der Waals envelope. Any protein whose groups match the
arrangement therefore has its body centred on the target surface, and roughly half of it lands inside
the target at **any** orientation. That is the body-clearance finding, a day early, mis-filed as an
instrument problem.

⇒ ⭐ **When an instrument fails because "everything overlaps", ask whether the overlap is the
finding.**

---

## 6. ⛔⛔ Match the thing that was specified, not the thing your builder emitted

The first two searches matched on backbone **Cα**, then **Cα + Cβ**, and both failed. The cause:

- The design specification constrained **only functional tip atoms** — guanidinium nitrogens,
  imidazole ring atoms. The backbone was deliberately left free.
- But the coordinate file fed to the search carried backbone positions anyway, from **one idealised
  rotamer** the builder program emitted: every arginine at χ1 = −176.4°, χ2 = 177.2°, every histidine
  at χ1 = −66.9°, χ2 = 90.2°, spread under **0.07°**.

⇒ ⛔⛔ **Matching on those atoms matched the builder's default, not the design.** And the penalty is
quantified: over 13,058 arginines in real structures, superimposing a guanidinium group *exactly*
leaves the two Cα atoms agreeing to within 1.0 Å only **4.5%** of the time. ⇒ a chemically perfect
match fails a Cα gate essentially always.

⇒ ⭐ **A degree of freedom that nothing constrained will be filled in by whatever wrote the file, and
it will look like data.** Before searching on a coordinate, check that something actually specified
it. A tell-tale sign: zero variance across instances where real structures have a broad
distribution.

⭐ **The same reframing made the search three to five orders of magnitude cheaper.** A 250-residue
chain holds ~35 functional groups but ~250 residues, so enumerating typed functional groups instead
of backbone positions shrinks the subset space enormously. A full-database search dropped from hours
to seconds.

---

## 7. ⛔ A point-set search is not selective without context, and the standard fix does not always apply

Measured, before the functional-group reframing: a **Cα-only, identity-free five-point search found
the query arrangement in 73% of all protein chains**, and three control arrangements drawn from real
chains matched at the same rate or higher, at every RMSD cutoff down to 0.25 Å. Tightening the cutoff
never inverted the ordering in the query's favour.

The combinatorics force it: a 170-residue chain holds C(170,5) ≈ 1.1 × 10⁹ five-subsets, and ten
pairwise distance constraints cannot cut that to a small number.

The standard cure is to match **segments** rather than isolated points — requiring i−2…i+2 around
each matched position turns 5 constraints into 25 — which is what established substructure-search
tools do.

⛔ **It could not be used here, because the query had no backbone context to match.** The specified
positions were isolated single residues with no neighbours; there was no i±2 to constrain.

⇒ ⭐ **Check that your query can express the constraint your tool needs before adopting the tool.**
The right instrument for a specification of isolated functional groups is a typed functional-group
search, not a backbone substructure search — and the absence of selectivity was the diagnostic that
pointed at it.

---

## 8. Correctness machinery that earned its place

Each of these caught something real, or was cheap enough that the one that did not still paid for
itself.

### 8.1 ⛔⛔ Never summarise a sorted file's prefix

A set of clash statistics was computed on **the first 25 rows of a 58-row file that was sorted
ascending by the very column being summarised.** Published: min 52, median 92, mean 93. Correct over
all 58 rows: min 52, median **129.0**, mean **131.86**.

⇒ ⛔⛔ **The conclusion's sign reversed.** The text concluded that matched placements clash *slightly
less* than a random control (92 against ~116). The correct comparison is 129 against ~116: they
clash **more**. The correct number was printed nine lines above the wrong one, in the same file.

⇒ ⭐ And the corrected sign had a mechanism, which is why it is the believable one: 65.5% of matched
chains had their centre of mass on the target side of the arrangement, median +3.67 Å. ⇒ **a sign
that reverses and then acquires a mechanism is usually the right sign.**

### 8.2 ⭐⭐ A vacuity control: a gate that can only return FAIL is not a gate

Before trusting a 2.4% pass rate, the accessibility gate was re-run in a mode that replaced the
matched residues with the **six most-exposed** cationic and imidazole residues of the same chain.
**11 of 12 returned a pass**, with both channels maxed and minimum exposure 0.675–0.969.

⇒ ⭐ **The pass branch is reachable, so every failure is a verdict rather than an artefact of an
unsatisfiable rule.** Without this, "almost everything fails" and "the rule is impossible" are the
same observation. Cheap, and it should be standard for any new gate.

⚠️ One detail: a control run still prints the main run's prediction lines, where they are meaningless.
Label control output as control output, or someone will quote it.

### 8.3 ⭐⭐ Assert identity, not just presence: the wrong-target guard

Every downstream step looked residues up by number in the real structure. The gate asserted that the
residue found at a recorded number **must** match the recorded group type — a cation record must find
an arginine or lysine, an imidazole record a histidine — and excluded any mismatch from every
statistic rather than measuring it silently.

A rule was fixed in advance: at a 5% mismatch rate or higher, the results must not be published at
all, because the mapping between the search store and the real structures would be broken.

**Measured: 0 mismatches in 1,314.** The guard fired zero times, which is the outcome you want and
not a reason to have skipped it: it is the difference between "the numbers are right" and "the
numbers were never checked".

⇒ ⭐ This is the guard that catches the hardest class of bug in any pipeline that re-reads source
data by index — **the right command run against the wrong target.** It produces plausible numbers and
nothing else detects it.

### 8.4 ⭐⭐ Reproduce the old numbers before believing the new ones — and don't over-specify the gate

A corrected run must reproduce the *original* numbers at the *original* parameters, exactly, before
any corrected number is reported. Here that meant counts, best RMSD, **and** the identity of the best
subset at each point count — reproduced to the digit, in both databases. ⇒ the arithmetic and the
matcher were never wrong; **one parameter was**, and that distinction is what made the correction
publishable rather than a rewrite.

⚠️ **The trap inside it.** The reproduction gate was written to require a best RMSD of "0.5720 Å",
because the published source printed **0.572** and that was extended to the driver's four decimals.
The true value is **0.5718**. Enforcing the invented fourth decimal would have failed the gate on a
rounding artefact, and the next step would have been hunting a database mismatch that did not exist.

⇒ ⛔ **Never state a reproduction gate to more precision than its source carries.** The bar itself was
the rounded thing.

### 8.5 ⭐ If you must copy an instrument to instrument it, prove the copy equivalent

The matcher could not report one diagnostic without modification. Rather than edit the reference
implementation, its core routine was copied verbatim with a flag added, and the copy was then proved
equivalent: 4 subsets × 2 tolerances × 400 chains = **3,200 paired comparisons**, 0 mismatches in
both hit count and best RMSD. The original was imported and never edited, so it remained the
reference against which everything else is checked.

### 8.6 ⛔⛔ An overflow that returns the same value as "no result" must be counted separately

The matcher abandons a chain whose intermediate candidate set exceeds a cap by returning **an empty
list** — the same value it returns for *no match*. ⇒ a resource overflow and a true zero are
indistinguishable in the output.

⇒ ⛔ **Count them separately and report the count per cell.** This matters most in exactly the runs
meant to replace a zero, because those runs are the ones with inflated intermediate sets.

⚠️ **Measured here: zero overflows in 84 cells of each database plus 12 control cells.** So the
defect was **latent, not active**, and it must not be claimed to have corrupted anything. Reporting a
prediction that failed in the direction of "the suspected bug did not bite" is part of the same
discipline as reporting one that succeeded.

### 8.7 ⭐ Verify a saved artefact by reloading it, in a fresh process

A builder that exits successfully has not verified anything. Two failure modes here look exactly like
success at build time: a non-ASCII character stored in a viewer session's metadata crashes the file
when it is next opened, and a named selection that matched nothing saves happily and displays nothing.
The verifier reloads the saved file in a fresh process, asserts every object non-empty with its atom
count, and **walks the entire saved session asserting every stored string is ASCII** — which is
stronger than checking the specific fields, because it also covers labels and object names and does
not rot when the viewer moves its internal layout.

⭐ A related rule for the same toolchain: under a headless viewer invocation the script's module is
**not** `__main__`, so a conventional `if __name__ == "__main__":` guard makes the program exit
successfully having written nothing.

### 8.8 ⭐⭐ Gate a picture on the number it claims to depict

The superposition builder **refuses to write output** unless the RMSD it achieves equals the
published value for that chain to 1e-3.

⇒ ⭐⭐ **A superposition that does not reproduce the published number is a different alignment**, and
a figure built on it shows a measurement nobody made. All four candidates passed to the digit. This
is a two-line check that makes a figure auditable, and it generalises to any visual built from a
scored table.

### 8.9 ⛔⛔ An inapplicable check is not a passing check — and this cost nothing only because it was caught twice (2026-11-01)

Two instances, in two different programs, both of which **reported success**:

1. A generic input validator ran over a new construct and its summary read as clean. Two of its
   checks had not passed: one asserts that every fixed residue holds exactly three atoms, which is
   false by design once whole backbone segments are held (sizes run 3 to 14), and the other's
   arithmetic test was **skipped** because its precondition — that the fixed set contains only
   guidepost atoms — did not hold. A skipped check is not a satisfied one. ⇒ both underlying
   properties were re-verified directly, and the guideposts are now identified by their **atom-name
   set** rather than by "has three fixed atoms", which is a coincidence test.
2. A sequence-design driver reconciles a repaired coordinate file against the raw file, because an
   earlier cohort's files lost an annotation during repair. This cohort needed **no repair**, so the
   same file was passed as both arguments and the comparison compared a file with itself. It
   reported a pass.

⇒ ⭐⭐ **Record such a check as VACUOUS, by name, rather than counting it among the passes.** The
failure mode is not a wrong answer; it is a true answer to a question nobody asked, sitting in a
column that looks like evidence. ⚠️ This is the complement of 8.2: there, a gate that could only
return FAIL; here, a check that could only return PASS.

### 8.10 ⭐⭐ Record the configuration the run *resolved*, not the one you passed (2026-11-01)

A rate-measuring run on one arm wrote out its **resolved** sampler configuration, and the full
screen's ten shards were compared against it: identical in all ten.

⇒ ⭐⭐ **That comparison is the only thing that licenses comparing the test to the screen.** Flags are
passed through layers that have their own defaults, so "the flag was passed" and "the run used the
value" are different claims. ⭐ The cheap discipline: have the run print what it actually resolved,
to a file beside the output, and diff it between runs before comparing their numbers.

⚠️ The same run exposed the opposite case: its launch record carried no field for which accelerator
was used, so the hardware a result was produced on is **not auditable from the record** — see §18.

### 8.11 ⭐⭐ A control built to fail is how you prove a column reads what it claims (2026-11-01)

One arm of a ten-arm screen was built around a fragment already **known** to sit inside the target's
volume, and its input was measured and recorded as carrying **379** close atom pairs against the target
*before anything ran*. In the result it scores **0 of 40** on clearance with a median backbone distance to
target of 0.379 Å — while its continuity is **perfect**, 8 of 8 in every batch.

⇒ ⭐⭐ Two things at once, and the second is the transferable one: the element is confirmed as the
part to drop, **and the clearance column is confirmed to be reading clearance**, because something
designed to fail on exactly that column did. ⭐ A screen in which nothing fails has not demonstrated
that its discriminating column works. ⚠️ And note it passes the *other* columns, which is what makes
it a clean control rather than simply a bad design.

---

## 9. Counting lessons, which sound trivial and each cost something

- ⛔ **`wc -l` is the wrong instrument for a record count.** A store's chain count read from its own
  offset array was 85,699; read by counting newlines it was 85,698, because the file had no trailing
  newline. The error was in the direction that would have been published.
- ⛔ **A store may count entries, not unique names.** A 795-row hit list held **786 distinct chain
  names** — nine chains entered twice from two overlapping source subsets. This invalidated nothing,
  because every curve counted entries the same way, but "795 unique chains" was wrong and had been
  written.
- ⛔ **A key count is not a coverage count.** Taking distinct identifiers from a database's lookup
  table gave 224,856; the number actually searchable, by joining the index to the names, was
  **151,880**. The lookup listed entries the database did not hold.

⇒ ⭐ The shared shape: **the thing counted was not the datum.** Count the datum, from the structure
that defines it, and say which structure that was.

**Two more of the same kind, 2026-11-01:**

- ⛔⛔ **Concatenating files that lack a trailing newline undercounts by exactly `n_files − 1`.** A
  tool wrote one output file per input, each ending **without** a newline. Concatenating them and
  counting record headers reported **129** where the true per-file sum was **144** — short by exactly
  the 15 file boundaries, because each file's last data line fuses onto the next file's header and
  that header stops being at the start of a line. ⇒ ⭐ **parse each file on its own, and assert the
  per-file count against what the run was asked to produce.** The undercount is silent, it is
  proportional to how many files you have, and it looks like a plausible number.
- ⛔ **The first record of each of those files was the *input*, not an output.** It carries no output
  identifier, which is the only thing distinguishing it. Counting it as an output would have inflated
  every yield by one ninth. ⇒ ⭐ **label it and exclude it, rather than trimming it away** — a
  discarded row cannot be checked, and a labelled one can.

---

## 10. ⛔ Infrastructure: a mirror script that deletes its own log

The public archive mirror was fetched with the distributor's own template script. That template sets
its log file **inside** the mirrored directory, while the transfer runs with deletion of
extra-to-source files enabled. The remote tree has no such entry, so the run **deletes its own log
while writing it.**

**Measured:** after a completed run the log was absent and the directory held no stray top-level
files. ⇒ the process exit status was the **only** surviving evidence of completion, and the mirror had
to be verified against the filesystem directly instead — 87 GiB, 1,294 hash directories, 261,584
compressed entries, integrity spot-checked. The data were fine; the defect cost the log, not the
mirror.

⇒ ⭐ **Never write a transfer's log inside the tree the transfer is allowed to delete from.** And ⭐
**when the only evidence of success is an exit code, go and measure the result.**

---

## 11. ⛔⛔ And the process finding, which is the one the group should keep

The specification problem this line was meant to help with — a designed binder too short to span the
site it was asked to present — had been **identified and written down three days before the mining
line ended**, with options priced. None was taken. The three days went to control arms and tolerance
re-derivation: work that characterised the *instrument* and moved nothing toward a better structure.

The mining line itself then produced, in order: a negative, a null, a retracted filter, a corrected
curve, a named list, two oppositions and a clash result. Every one of those is a real measurement.
**None of them is a design.**

⇒ ⭐⭐ **Characterising an instrument is not progress on the object.** Both are necessary and they are
not interchangeable, and a line of work can be rigorous, well controlled, honestly reported, and still
be the wrong thing to have been doing. ⇒ the useful question to ask of a measurement, before running
it, is not *"is this sound?"* but *"what would be done differently depending on the answer?"* — and, here,
a count of matching chains had no action attached to either outcome until the identities were emitted.

---
---

# Appended 2026-11-01 — from the scaffolded-design line

The work these come from is in `methods/core_scaffolded_design.md`: keep
subsets of a real protein's secondary-structure elements as a rigid core, generate the rest around
them with the target present, screen ten such subsets, design sequences on the survivors. ⛔ None of
it is a binding result, which is itself the point of §18.

## 12. ⛔⛔ Before you trim, find out *where* the problem is — a middle is not a tail

The line in §5 ended with candidate scaffolds whose bodies interpenetrate the target: the best of
them put **679** heavy atoms within 2.5 Å. The obvious repair is to shorten the chain.

**MEASURED:** trimming **both termini** removes **none** of the 679. Holding **subsets of
secondary-structure elements** instead, each re-placed on its own, takes it to **zero** for nine of
ten subsets tried, and where a closest-approach distance was recorded alongside, the best of them sit
at 2.54 Å — ordinary packing between two chains rather than interpenetration.

⇒ ⭐⭐ **The collision was in the middle of the chain, so no amount of trimming could reach it.** The
useful generalisation is not "decompose your scaffold"; it is that **"too big" and "wrong in the
middle" are different diagnoses with the same summary statistic**, and one cheap measurement — the
per-residue or per-element contribution to the total — separates them before any repair is attempted.
Attributing a total to its parts costs minutes and chooses the repair for you.

⚠️ **And the decomposition has a cost that the total also hides.** Of the six functional groups the
site specification named, two sit inside the kept backbone and four do not; the two inside are inside
a segment whose *sequence* was deliberately released, so only **four** remain
identity-constrained downstream. ⇒ ⭐ **a construct that fixes a geometric problem can quietly relax a
chemical constraint**, and the only way to notice is to re-derive what is still constrained **from the
deposited input**, rather than from the intent.

## 13. ⛔⛔ Count the degrees of freedom before you report a fit, and never compare fits over different point counts

The kept fragment's placement is set by superimposing the specified points it **contains**. That
number varies with the element set, and it changes the kind of answer available:

| contained points | what a least-squares fit determines |
|---|---|
| ≥ 3 | the rigid body, fully |
| **= 2** | ⛔ **everything except one rotation about the axis through them.** That rotation is real and must be chosen on some stated criterion — here, minimise collision with the target, swept in 5° steps |
| ≤ 1 | nothing. The placement is arbitrary and the option was **rejected** |

⇒ ⭐⭐ **A two-point "fit" is not a fit until the surviving rotation is disposed of, and how it was
disposed of is part of the result.**

⛔⛔ **And the resulting number is not comparable to a fit over more points.** The two-point
placements here report 0.9187 Å and the three-point control 0.8644 Å, against **1.5400 Å** for the
six-point superposition of the same chain. The small numbers are not better alignments; fewer
constraints simply cost less to satisfy. ⇒ ⭐ **always print the point count next to an RMSD**, and
treat a fit that improves when you drop constraints as arithmetic rather than progress. This project
has now been bitten by three mutually incomparable radii in one analysis and by this; the rule that
covers both is that **a distance is meaningless without the set it was computed over.**

### 13.1 ⭐⭐ The impossible number is the cheapest bug detector you will ever build

The defect that invalidated the first run of the option designer was found this way. A **four**-point
least-squares fit was reported as **worse** than the **six**-point fit of the same chain — 7.101 Å
against 1.5400 Å. That cannot happen: the six-point problem contains the four-point problem as a
relaxation, so the sub-fit can only be better or equal.

The actual cause: the point correspondence was hard-coded as the residue numbers **sorted
numerically**, with a comment asserting that this was the query order. It was not. Three of six
points were therefore fitted onto the wrong targets.

⇒ ⭐⭐ **Look for a quantity whose ordering is guaranteed by mathematics rather than by your code, and
assert it.** Monotonicity under relaxation, a subset count never exceeding its superset, a
sub-problem's optimum never beating the full problem's — each is free, each needs no ground truth,
and each catches a whole class of silent correspondence errors. ⚠️ Note what did **not** catch it:
every program in this line asserts that a matched residue's **type** agrees with the specified type,
and two of the three mismatches were type mismatches — but this one program had omitted that
assertion. ⇒ ⛔ **a guard is only a guard where it runs**; a convention that holds in "the pipeline"
and not in a particular file is a convention with a hole in it.

## 14. ⛔⛔ Gate first on the failure that your next metric rewards

Designs were scored for (a) backbone continuity — the largest step between consecutive residues —
and (b) minimum distance from the backbone to the target.

⛔ **Read in that order (b) is corrupt**, because a chain that has broken into pieces can let its
pieces retreat from the target in a way an intact chain cannot. A clearance floor therefore **rewards
a break**. Reading clearance first would have put the broken designs at the top of the ranking.

⇒ ⭐⭐ **Order your gates so that each one is measured on a population where it means what you think.**
The general test is cheap to apply: for each metric, ask *what does an outright failure of some other
kind do to this number?* — and if the answer is "improves it", that other kind gets gated first.

⚠️ **Related, and it is why these two were never collapsed:** continuity and fold topology are
**orthogonal** here. One arm scored 38 of 40 on continuity and **0 of 40** on whether the chain is a
single fold — it is continuous and it is two lobes on a tether. Overlap between the two failure sets
was essentially nil. ⇒ ⛔ **a chain can be unbroken and still not be one protein**, and no single
backbone number detects both.

### 14.1 ⭐⭐ A threshold is admissible only where the counts stop moving

Before applying a binary continuity cut, the gate re-measured its own bar ladder on **this** cohort:
303 / 346 / 350 / 351 / 353 / 371 designs pass at 4.0 / 4.1 / 4.5 / 5.0 / 7.0 / 10.0 Å.

⇒ the counts at 5, 7 and 10 Å **differ**, so the population genuinely occupies the space between
"ideal" and "broken". The rule the gate enforces is that a cut may be used only where moving the bar
over a plausible range does not move the count — and here it moves. ⛔ **So no verdict column was
written at all**, and the bar quoted in the write-up is labelled a *reporting* bar.

⇒ ⭐⭐ **Measure your threshold's sensitivity on the cohort you are about to apply it to, and if the
answer is "sensitive", publish a distribution instead of a verdict.** A threshold is a claim about
where a gap in the data is. If there is no gap, the threshold is inventing one, and every count
downstream inherits the invention.

## 15. ⛔⛔ `n` is the number of independent draws, not the number of files

Ten arms produced 40 design files each. But the generator draws **one chain layout per batch** and
emits eight models from it, so the eight share a draw. Effective `n` per arm is **5**, not 40.

⇒ ⛔ **"40 of 40" is "5 of 5".** Every count in the write-up is reported with its per-batch
breakdown, for exactly this reason, and the aggregation program prints the per-batch table so a 40/40
cannot read as strong evidence on its own.

⭐ **What the per-batch view buys, concretely:** it distinguishes a consistent effect from one lucky
draw. Three arms are 8/8 in **every** one of their five batches, while the worst runs 4/8, 2/8,
**0/8**, 3/8, 6/8 — the same mean would have hidden that. ⇒ **the batch table is not a caveat, it is
the evidence.**

⚠️ **And the headline result of that screen is still only 4 arms against 6.** The separation is clean
— the worst arm holding one particular large helix beats the best arm without it, with no overlap —
and it is *still* a comparison of ten groups of five. ⇒ ⭐ it was written up as **"strongly suggested
and worth designing around"**, explicitly not as established. ⛔ A clean separation at small `n` is a
hypothesis with a good shape, not a result.

⚠️ **One confound was named rather than resolved:** the arms holding the large helix also hold
**fewer and larger** segments overall, and nothing measured here separates the helix from the
segmentation. Naming a confound you cannot remove is part of reporting the effect.

## 16. ⛔⛔ Do not inherit a settled parameter across a changed construct — and a small response to a large push is a finding, not a failure

A composition bias had been fitted, over nine rounds, on an earlier cohort of this project, and was
the settled production setting. The new cohort differed: **141–148 residue chains built around a held
fragment of a real protein**, against **52–79 residue de novo bundles**.

**MEASURED on the new cohort before anything was inherited:** its dominant compositional excess is
**glycine** (13.18% against 7.07% natural). The settled bias's largest term is **alanine**, which
here sits at 11.24% against that cohort's 21.03%; its glycine term is the **smallest** of its four;
and two residues it suppresses are already **below** natural here, so it would have pushed them the
wrong way. The unbiased composition of the new cohort was already closer to natural than the old
cohort reached after all nine of its rounds.

⇒ ⭐⭐ **A parameter fitted to a pathology is only valid where that pathology is.** The standing rule
the earlier work itself recorded — *measure this cohort rather than assume the last one's* — is the
transferable part, and it is worth more than the fitted values it governs. ⛔ Inheriting the numbers
and inheriting the **method** are different acts; only the second one transferred here.

### 16.1 ⭐⭐ Response per unit of pressure separates "mis-tuned" from "constrained"

After refitting on the new cohort: **glycine moved 0.60 percentage points under a cumulative
adjustment that, at the sampling temperature in use, is an effective shift of −3.07 in the logit.
Alanine moved 1.71 points under a shift less than half as large.**

⇒ ⭐⭐ **The residue that barely moves under the larger push is the one the structure is demanding, not
the one the sampler is over-producing.** That ratio — change in outcome per unit of applied pressure
— is a diagnostic available from the refit itself, with no extra run: a knob that pushes hard and
achieves little is telling you it is attached to the wrong variable. ⭐ Here the reading is that the
backbones genuinely need glycine, which is consistent with most of each chain being generated
connector between held elements, and turns needing glycine.

⚠️ **Stated as the leading explanation and not as a measured cause.** The test would be glycine
content against loop fraction per design, and it **was not run.** ⇒ a mechanism that explains a
measurement is still a hypothesis until it predicts one.

⚠️ **A carried limit, kept because it bounds the whole exercise:** on another target in this project
the same kind of bias was **null for design quality** (paired test, p = 0.657), and here its cost in
model confidence was −0.0007, i.e. nothing. ⇒ ⛔ **it earns its place on composition plausibility
alone.** Do not expect it to buy foldability, and do not justify it by diversity.

## 17. ⛔⛔ A composite score defined over a variable support is not comparable across runs that change the support

The composition work reports an aggregate "deviation from natural" — the summed absolute difference
over the **designable** residue set. That set is defined by which residues the run was told to omit.

⇒ ⛔ **Changing the omission list changes the support, the renormalisation of the observed
frequencies, and the number of terms in the sum — so two such scores computed under different
omission lists are not the same quantity.** Re-measuring the same cohort under the two omission lists
used in this work gives **32.3** and **30.9**; the figures quoted in the source material pair
different supports against each other, which makes a before-and-after look like an improvement partly
driven by the definition. ⚠️ A third figure in that material, **30.4**, matches neither
re-measurement and most likely belongs to a **separate** unbiased run under the second support whose
sequences were not kept — in which case it is a different cohort rather than a different calculation,
and it cannot be reproduced from anything deposited. ⇒ ⛔ **none of the three is quoted as a
headline**; the per-component values are published instead.

⇒ ⭐⭐ **Publish the per-component values, which are comparable, and quote a composite only with its
support stated.** Here the per-residue percentages carry the entire finding and the composite carries
none of it. ⭐ The general rule: a scalar built by summing over a set is a statement about that set,
and the set belongs in the label.

## 18. ⛔⛔ A quantity you never sampled is not a reassurance — it is a gap with the shape of one

The generative runs in this line were deliberately placed on hardware large enough that memory could
not be a confound, **on the assumption that the runner would report memory use.** It does not sample
it at any point. ⇒ **no memory figure exists for any of these runs.**

The consequence is precise and it is not "probably fine": it is **unmeasured** whether the same work
fits on smaller hardware, so **no sizing recommendation can be made from a completed, successful
campaign.** A run that succeeds tells you the resource was sufficient and nothing about how much of
it was needed.

⇒ ⭐⭐ **Decide what you will need to have measured *before* the run, because a successful run cannot
be interrogated afterwards for a quantity it did not record.** Two further instances from the same
runs, both one-line fixes that were not in place: the execution record carried no field for which
accelerator was used, so a result's hardware is not auditable from it; and its timestamp is when the
record was **written**, not when the work started, understating duration by minutes. ⇒ ⛔ **an
execution record that cannot answer "on what, and for how long" is a log, not provenance.**

⚠️ **And the rate that *was* measured comes with a comparability warning**, which is the same shape
as §13: an earlier rate from this project differs in **both** chain length and hardware, so the two
numbers cannot be differenced. ⇒ quote them separately or not at all.
