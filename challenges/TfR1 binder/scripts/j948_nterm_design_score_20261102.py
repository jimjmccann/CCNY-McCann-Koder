#!/usr/bin/env python3
"""j948 -- score the N-terminal-extension arm's RFd3 designs. Local, reads .cif.gz on disk.

⛔ NO network, NO GPU, NO spend, nothing launched. Reads each design's `.cif.gz` and its sibling
`.json`, writes ONE TSV. ⛔ Applies NO verdict column and NO binary cut.

=== WHY THIS EXISTS RATHER THAN THE STANDARD GATE-CHAIN DRIVER ===
`[M]` 2026-11-02: that driver REFUSES this arm. Its G0 discriminator assigns a chain by fixed-sequence
fraction (`>90%` target, `<20%` binder, else AMBIGUOUS). This arm holds 143 of 170 binder residues
fixed, so the binder chain reads **84.1% fixed** and lands in the dead zone:
    `AMBIGUOUS chain(s) ['A'] at fixed% [84.1] -- NOT scored`  ⇒ 0 of 2 resolved.
⭐ The driver did the RIGHT thing -- it refused rather than silently calling the binder a target, and
printed "absent from the TSV is NOT a pass". ⛔ Do NOT widen its thresholds: a band that admits
84% would admit a genuine target too. This is the documented failure mode
[[partial-diffusion-output-breaks-chain-role-detection]], whose prescription is: **read the role
from the SPEC.**

⇒ ⭐⭐ ROLES COME FROM EACH DESIGN'S OWN `.json`, never from a flag and never from a chain letter
([[always-name-the-pdb-chain]]: a SCORER must never name one). `specification.contig` plus
⭐ CORRECTED 2026-11-02 (methods/nterm_extension_arm.md section 5): the generated set is the binder chain's OWN residues minus the
spec's held residues. ⛔ NOT derived from `diffused_index_map`, which names only TOKENS -- a
generated run has no token, so the old derivation excluded an N-terminal extension entirely and
measured the 2-residue linker instead.
`diffused_index_map` map every contig token to an output chain: the token stream's generated (`NN-NN`)
and held-protein (`H*`/`E*`) entries land on the BINDER chain, the `A<resi>` target-crop entries land
on the TARGET chains. The resolved layout is PRINTED BEFORE ANY ROW IS WRITTEN.

=== WHAT IS MEASURED, AND THE ORDER IS DELIBERATE ===
1. CONTINUITY FIRST, before any clash or interface number
   [[gate-backbone-continuity-before-any-clash-metric]]. Measured reason: a clash floor REWARDS a
   break -- a chain with half of itself 56 A away has fewer close contacts. Reported as a
   DISTRIBUTION (max_step, n_break_4p5), ⛔ never pass/fail.
2. COMPACTNESS as `rg_ratio` = Rg / (2.2 * n^0.38), the only length-safe column
   [[continuity-and-fold-topology-are-orthogonal]] -- continuity and topology are orthogonal and
   must not be collapsed.
3. CLASH, binder-vs-target and intra-binder SEPARATELY, as a continuous minimum heavy-atom
   distance. ⛔ Raw RFd3 output clashes normally [[raw-rfd3-output-clashes-normally]], so this is
   reported, not gated.
4. ⭐⭐ THE QUESTION THIS ARM WAS BUILT TO ANSWER: does the NEW 25-mer reach the target, or grow
   away from it? The arm's own prediction (its builder's note) is that residues 1-3 of the parent
   design already sit 3.46 / 4.09 / 4.13 A from TfR1, so "the only direction with room is away from
   the interface". ⇒ `new_contacts_4p5` and `new_min_dist` are measured over the GENERATED residues
   only, and they are what can falsify it.
5. The held arginine's guidepost clearance, on its PIN HEAVY ATOMS (NE, NH1, NH2) and ⛔ never on CA
   [[gate-guidepost-clearance-on-pin-heavy-atoms-not-ca]] -- CA-gating once passed a 1.67 A clash.

=== THREE TRAPS THIS HANDLES EXPLICITLY ===
⛔⛔ PSEUDO-ATOMS. Raw RFd3 output carries V0-V8 pseudo-atoms
   [[rfd3-raw-output-carries-pseudoatoms]]. They are dropped by matching element `V` AND a name
   `V<digits>`, and the dropped COUNT is emitted per design. ⛔ A text match on the whole file
   returns zero and reads as a FALSE CLEAN.
⛔⛔ n IS BATCHES, NOT DESIGNS [[rfd3-designs-are-clustered-by-batch]]. One contig length is drawn
   per BATCH, so the models inside a batch share a draw and are NOT independent. The batch id is
   emitted per row so a 24/24 that is really 3/3 cannot read as strong evidence.
⛔ A reported value is not a measurement: RFd3 writes its OWN `metrics` block into each json
   (`n_chainbreaks`, `max_ca_deviation`, `radius_of_gyration`). Those are carried through with an
   `rfd3_` prefix for COMPARISON ONLY, next to the columns computed here from the coordinates.
   ⛔ They are the model's claims about itself and are never used in place of a measurement.

usage:
  j948_nterm_design_score_20261102.py --outdir <dir> '<glob>' ['<glob>' ...] [--limit N]
"""

import argparse
import glob as globmod
import gzip
import json
import math
import os
import re
import sys
import time
from collections import defaultdict

import numpy as np

BB = ("N", "CA", "C", "O")
PSEUDO = re.compile(r"^V\d+$")
CONTACT = 4.5          # heavy-atom contact bar, A
BREAK = 4.5            # consecutive CA-CA step above this = a break, A
IDEAL_CA = 3.8


def openmaybe(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def read_cif_atoms(path):
    """Minimal mmCIF _atom_site loop reader. Returns (atoms, n_pseudo_dropped).

    atoms: list of (chain, resi, resname, atomname, element, xyz)
    ⛔ Column order is READ FROM THE HEADER, never assumed -- an index guess is how a chain gets
    silently mislabelled.
    """
    cols, rows, in_loop, header_done = [], [], False, False
    with openmaybe(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                in_loop, cols = True, cols + [s.split(".", 1)[1]]
                continue
            if in_loop and not header_done:
                if s.startswith(("ATOM", "HETATM")):
                    header_done = True
                else:
                    if s.startswith("#") or s.startswith("loop_") or not s:
                        continue
            if header_done:
                if s.startswith(("ATOM", "HETATM")):
                    rows.append(s.split())
                elif s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
                    break
    if not cols or not rows:
        raise ValueError(f"no _atom_site loop parsed from {path}")
    ix = {c: i for i, c in enumerate(cols)}

    def pick(*names):
        for n in names:
            if n in ix:
                return ix[n]
        raise KeyError(f"none of {names} in {sorted(ix)}")

    i_ch = pick("auth_asym_id", "label_asym_id")
    i_ri = pick("auth_seq_id", "label_seq_id")
    i_rn = pick("auth_comp_id", "label_comp_id")
    i_an = pick("auth_atom_id", "label_atom_id")
    i_el = pick("type_symbol")
    i_x, i_y, i_z = ix["Cartn_x"], ix["Cartn_y"], ix["Cartn_z"]

    atoms, dropped = [], 0
    for r in rows:
        if len(r) <= max(i_ch, i_ri, i_rn, i_an, i_el, i_z):
            continue
        an, el = r[i_an].strip('"'), r[i_el]
        # ⛔⛔ pseudo-atoms: element AND name, both required
        if el.upper() == "V" and PSEUDO.match(an):
            dropped += 1
            continue
        if el.upper() == "H":
            continue
        try:
            resi = int(r[i_ri])
        except ValueError:
            continue
        atoms.append((r[i_ch], resi, r[i_rn], an, el,
                      np.array([float(r[i_x]), float(r[i_y]), float(r[i_z])])))
    return atoms, dropped


def roles_from_json(jpath):
    """⭐ Resolve binder/target chains and the GENERATED residue set from the design's own spec.

    Returns (binder_chain, target_chains, generated_resi, pin_tokens).
    ⛔ Raises if the spec cannot settle it -- an unresolved role is never guessed.
    """
    d = json.load(open(jpath))
    imap = d.get("diffused_index_map") or {}
    if not imap:
        raise ValueError("no diffused_index_map")
    spec = d.get("specification") or {}
    pin_tokens = spec.get("select_fixed_atoms") or {}

    # token -> "A26" ; a token naming a target crop position looks like A<resi>
    per_chain = defaultdict(list)
    target_tok, held_tok = defaultdict(list), defaultdict(list)
    for tok, out in imap.items():
        m = re.match(r"^([A-Za-z]+)(\d+)$", str(out))
        if not m:
            continue
        ch, ri = m.group(1), int(m.group(2))
        per_chain[ch].append(ri)
        if re.match(r"^A\d+$", tok):        # target-crop token, e.g. A193
            target_tok[ch].append(ri)
        else:                                # held protein / held pin token, e.g. H4, E1
            held_tok[ch].append(ri)

    binder = [c for c in per_chain if held_tok.get(c)]
    target = sorted(c for c in per_chain if target_tok.get(c) and c not in binder)
    if len(binder) != 1 or not target:
        raise ValueError(f"ambiguous roles from spec: held={sorted(held_tok)} "
                         f"target={sorted(target_tok)}")
    b = binder[0]
    # ⛔⛔ THE GENERATED SET CANNOT BE DERIVED FROM THIS MAP. `diffused_index_map` maps TOKENS,
    # and a generated run (`25-25`, `2-2`) has no token name, so it has NO entry here. The old
    # derivation took `range(min(mapped), max(mapped)+1) - held`, whose lower bound is the first
    # HELD residue -- so an N-terminal extension sitting BELOW it was silently excluded and
    # `generated` collapsed to the inter-element linker alone (see the arm methods note). The honest answer is
    # STRUCTURE minus HELD, and only the caller has the structure.
    held_resi = set(held_tok[b])
    return b, target, held_resi, pin_tokens


def rg_of(coords):
    c = coords.mean(axis=0)
    return float(np.sqrt(((coords - c) ** 2).sum(axis=1).mean()))


def mindist(a, b):
    if len(a) == 0 or len(b) == 0:
        return float("nan"), 0
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return float(d.min()), int((d < CONTACT).sum())


def score_one(cif, jpath):
    atoms, dropped = read_cif_atoms(cif)
    b, targets, held_resi, pin_tokens = roles_from_json(jpath)

    by_chain = defaultdict(list)
    for a in atoms:
        by_chain[a[0]].append(a)
    if b not in by_chain:
        raise ValueError(f"spec binder chain {b} absent from structure chains {sorted(by_chain)}")

    binder = by_chain[b]
    tgt = [a for c in targets for a in by_chain.get(c, [])]

    # ⭐ generated = what the STRUCTURE has on the binder chain minus what the spec HOLDS (the corrected derivation)
    generated = sorted({a[1] for a in binder} - held_resi)
    if not generated:
        raise ValueError(f"no generated residues on chain {b}: structure has "
                         f"{len({a[1] for a in binder})} residues and the spec holds "
                         f"{len(held_resi)} of them")

    # ---- 1. CONTINUITY FIRST, on the binder chain only
    ca = {}
    for ch, ri, rn, an, el, xyz in binder:
        if an == "CA":
            ca[ri] = xyz
    order = sorted(ca)
    steps = [float(np.linalg.norm(ca[order[i + 1]] - ca[order[i]]))
             for i in range(len(order) - 1) if order[i + 1] == order[i] + 1]
    max_step = max(steps) if steps else float("nan")
    med_step = float(np.median(steps)) if steps else float("nan")
    n_break = sum(1 for s in steps if s > BREAK)
    max_dev = max((abs(s - IDEAL_CA) for s in steps), default=float("nan"))

    # ---- 2. COMPACTNESS, length-safe
    bca = np.array([ca[r] for r in order]) if order else np.zeros((0, 3))
    n_res = len(order)
    rg = rg_of(bca) if n_res else float("nan")
    rg_exp = 2.2 * (n_res ** 0.38) if n_res else float("nan")
    rg_ratio = rg / rg_exp if n_res else float("nan")

    # ---- 3. CLASH, separately
    bh = np.array([a[5] for a in binder])
    th = np.array([a[5] for a in tgt]) if tgt else np.zeros((0, 3))
    bb_min, n_contact = mindist(bh, th)
    # intra-binder: non-adjacent residues only, backbone atoms
    bbb = defaultdict(list)
    for ch, ri, rn, an, el, xyz in binder:
        if an in BB:
            bbb[ri].append(xyz)
    ris = sorted(bbb)
    intra = float("inf")
    for i in range(len(ris)):
        for j in range(i + 3, len(ris)):
            A = np.array(bbb[ris[i]]); B = np.array(bbb[ris[j]])
            d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2).min()
            if d < intra:
                intra = d
    intra_min = float(intra) if intra != float("inf") else float("nan")

    # ---- 4. ⭐ THE ARM'S OWN QUESTION: does the NEW segment reach the target?
    gen = set(generated)
    gh = np.array([a[5] for a in binder if a[1] in gen])
    new_min, new_contacts = mindist(gh, th)
    n_gen_res = len({a[1] for a in binder if a[1] in gen})

    # ---- 5. the held pin's clearance, on PIN HEAVY ATOMS, never CA
    pin_names = set()
    for v in pin_tokens.values():
        pin_names |= {x.strip() for x in str(v).split(",") if x.strip()}
    ph = np.array([a[5] for a in binder if a[3] in pin_names and a[3] not in BB]) \
        if pin_names else np.zeros((0, 3))
    pin_min, _ = mindist(ph, th)

    # ---- carried, for COMPARISON ONLY: RFd3's own claims about itself
    d = json.load(open(jpath))
    m = d.get("metrics") or {}
    nclash = m.get("n_clashing") or {}

    name = os.path.basename(cif).replace(".cif.gz", "").replace(".cif", "")
    mb = re.search(r"_(\d+)_model_(\d+)$", name)
    # ⛔⛔ THE BATCH KEY MUST CARRY THE SHARD. `[M]` 2026-11-02: keying on the trailing index alone
    #    collapsed shard00's batch 0 with shard01's batch 0 and reported "3 batches" for a 4-shard
    #    cohort whose true independent count is 12. One contig length is drawn per batch PER SHARD
    #    [[rfd3-designs-are-clustered-by-batch]], and shards are separate processes on separate
    #    pods, so their batch indices are unrelated. ⇒ the key is shard+index.
    ms = re.search(r"_shard(\d+)_", name)
    shard = ms.group(1) if ms else "?"
    batch = f"s{shard}b{mb.group(1)}" if mb else "?"
    model = mb.group(2) if mb else "?"

    return {
        "design": name, "batch": batch, "model": model,
        "binder_chain": b, "target_chains": ",".join(targets),
        "n_res_binder": n_res, "n_gen_res": n_gen_res, "n_pseudo_dropped": dropped,
        "max_step": round(max_step, 4), "median_step": round(med_step, 4),
        "n_break_4p5": n_break, "max_ca_dev": round(max_dev, 4),
        "rg": round(rg, 4), "rg_ratio": round(rg_ratio, 4),
        "bb_min_to_target": round(bb_min, 4), "n_contact_4p5": n_contact,
        "intra_min_nonadj": round(intra_min, 4),
        "new_min_dist": round(new_min, 4), "new_contacts_4p5": new_contacts,
        "pin_min_to_target": round(pin_min, 4),
        "rfd3_n_chainbreaks": m.get("n_chainbreaks"),
        "rfd3_max_ca_deviation": m.get("max_ca_deviation"),
        "rfd3_rg": m.get("radius_of_gyration"),
        "rfd3_helix_fraction": m.get("helix_fraction"),
        "rfd3_clash_bb": nclash.get("interresidue_clashes_w_backbone"),
        "rfd3_clash_sc": nclash.get("interresidue_clashes_w_sidechain"),
    }


COLS = ["design", "batch", "model", "binder_chain", "target_chains",
        "n_res_binder", "n_gen_res", "n_pseudo_dropped",
        "max_step", "median_step", "n_break_4p5", "max_ca_dev",
        "rg", "rg_ratio", "bb_min_to_target", "n_contact_4p5", "intra_min_nonadj",
        "new_min_dist", "new_contacts_4p5", "pin_min_to_target",
        "rfd3_n_chainbreaks", "rfd3_max_ca_deviation", "rfd3_rg",
        "rfd3_helix_fraction", "rfd3_clash_bb", "rfd3_clash_sc"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--limit", type=int, default=0, help="probe mode: first N files only")
    ap.add_argument("globs", nargs="+")
    a = ap.parse_args()

    files = []
    for g in a.globs:
        files.extend(sorted(globmod.glob(g)))
    files = [f for f in files if f.endswith((".cif", ".cif.gz"))]
    if a.limit:
        files = files[:a.limit]
    if not files:
        sys.exit("⛔ no .cif/.cif.gz matched -- nothing to score (⛔ not a result of 0%)")

    os.makedirs(a.outdir, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    out = os.path.join(a.outdir, f"j948_nterm_score_{stamp}.tsv")

    print("=== j948 N-term arm design score ===")
    print(f"  files            : {len(files)}")
    print(f"  contact bar      : {CONTACT} A heavy-atom;  break bar: {BREAK} A CA-CA")
    print(f"  roles            : ⭐ FROM EACH DESIGN'S OWN .json (never a flag, never a chain letter)")
    print(f"  out              : {out}")

    rows, skipped, layouts = [], [], defaultdict(int)
    for f in files:
        j = f.replace(".cif.gz", ".json").replace(".cif", ".json")
        if not os.path.exists(j):
            skipped.append((f, "no sibling .json -- roles unresolvable"))
            continue
        try:
            r = score_one(f, j)
            rows.append(r)
            layouts[f"binder={r['binder_chain']} target={r['target_chains']}"] += 1
        except Exception as e:
            skipped.append((f, f"{type(e).__name__}: {e}"))

    print("\n  RESOLVED LAYOUT (⛔ from the spec, printed before any row is written):")
    for k, v in sorted(layouts.items()):
        print(f"    {k}   on {v} design(s)")
    if skipped:
        print(f"\n  ⚠️ {len(skipped)} file(s) NOT scored -- ⛔ absent from the TSV is NOT a pass:")
        for f, why in skipped[:12]:
            print(f"    {os.path.basename(f)}: {why}")
    if not rows:
        sys.exit("\n⛔⛔ ZERO designs scored -- that is a FAILURE, not a result of 0% anything.")

    with open(out, "w") as fh:
        fh.write(f"# j948 N-term arm design score, {stamp}\n")
        fh.write(f"# roles resolved FROM EACH DESIGN'S OWN .json spec (diffused_index_map)\n")
        fh.write(f"# ⛔ NO verdict column, NO binary cut. n is BATCHES not designs.\n")
        fh.write(f"# rfd3_* columns are the MODEL'S OWN claims, carried for comparison only\n")
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")

    nb = defaultdict(int)
    for r in rows:
        nb[r["batch"]] += 1
    print(f"\n  wrote {len(rows)} rows over {len(nb)} batch(es): "
          f"{dict(sorted(nb.items()))}")
    print(f"  ⛔ effective n is {len(nb)} BATCHES, not {len(rows)} designs")

    # distributions, never a cut
    def dist(col):
        v = [r[col] for r in rows if isinstance(r[col], (int, float))
             and not (isinstance(r[col], float) and math.isnan(r[col]))]
        if not v:
            return "n/a"
        v = sorted(v)
        return (f"min {v[0]:.3f}  med {v[len(v)//2]:.3f}  max {v[-1]:.3f}")
    for c in ["max_step", "n_break_4p5", "rg_ratio", "bb_min_to_target",
              "new_min_dist", "new_contacts_4p5", "pin_min_to_target"]:
        print(f"    {c:20s} {dist(c)}")
    print(f"\n  pseudo-atoms dropped, total: {sum(r['n_pseudo_dropped'] for r in rows)}")
    print(f"  ⇒ {out}")


if __name__ == "__main__":
    main()
