#!/usr/bin/env python3
"""Boltz-2 YAMLs: the 402 COMPETITION binders vs the FULL EGFR ECTODOMAIN, +/- glycans.

⭐ WHY. Every previous fold of these designs used a 106 aa (site A) or 77 aa (site B) CROP, and on
2026-10-02 that crop was found to delete receptor the binders occupy -- 70%/100% of our own designs
within 2 A of deleted backbone (measured 2026-10-02). The crop also makes the
measured-binding result uninterpretable in one specific way its own caveat names: a design whose
real epitope is elsewhere HAS NOWHERE ELSE TO DOCK, so "engages the footprint" is partly forced by
the construct. Folding against the WHOLE ectodomain removes both problems at once: nothing is
deleted, and every design is free to choose its own face.

⭐ THE TARGET IS NOT A MODEL CHOICE -- IT IS READ OUT OF THE ASSAY DATA. The 621-residue target is
the constant suffix of the `full_sequence` column in the competition's released table: a single
variant shared by 400 of 402 rows, so it is the assay construct AS PUBLISHED. Against our
structural model it has been checked residue-for-residue at the same WT index -- 604 of 604
MODELLED residues match, zero mismatches, which is 97.3% of the target. The remaining 17 residues
are NOT PRESENT in the model (WT 1-3, 295-299, 307, 577, 615-621) and are therefore unverified by
that comparison. Whether the published sequence matches the physical protein in the tube cannot be
verified here and is not claimed. (WT = mature = UniProt P00533 - 24.)
⛔ CORRECTED 2026-10-06, and the retraction is kept rather than the text silently rewritten: an
earlier version of this paragraph said the 621-mer "matches all 604 residues ... with ZERO
mismatches", which reads as whole-sequence verification. It is not -- 17 of the 621 residues are
absent from the model and were never checked.

⭐ GLYCANS. Nobody has ever folded these designs glycosylated, and the organiser confirmed
HEK293 / glycosylated / tethered. Occupancy is taken from MS on
NATURALLY EXPRESSED receptor (Zhen 2003) --
not from a crystal, which resolves only 1-7 sugars of a native 15-20-residue tree.
  --glycans none   aglycosyl control
  --glycans nag    one NAG per occupied sequon (minimal: attachment + proximal volume)
  --glycans core5  Man3GlcNAc2 core per occupied sequon (NAG-NAG-BMA-(MAN)(MAN))
⚠️ BOTH glycan modes UNDER-represent native trees. They are a lower bound on glycan interference,
never an upper one. GlycoSHIELD remains the right tool for the real exclusion volume.

⛔ MSA IS MANDATORY ON THE TARGET AND SETTLED BY MEASUREMENT: target CA-RMSD to 8UKW/B is
0.47-0.52 A WITH an MSA vs 1.18-1.48 A without (measured 2026-09-29). The binder keeps `msa: empty`
-- de novo sequences have no homologs (de-novo-proteins-break-plms). The a3m path written here
is WORKER-ABSOLUTE (`/workspace/msa/target.a3m`) to match the bundle builder's assert.

2026-10-02.
"""
import argparse, collections, csv, os, re, sys

AA = set("ACDEFGHIKLMNPQRSTVWY")

# ⭐ Occupancy by MS on naturally expressed receptor (Zhen 2003), WT numbering.
#    FULLY occupied -- the defensible "best guess" set.
SEQUONS_FULL = [151, 328, 337, 389, 420, 504, 544, 599]
#    ⛔ Deliberately EXCLUDED and why:
#      104, 172 -- measured NOT glycosylated.
#      579      -- partial.
#      32       -- atypical (NNC, not N-X-S/T) and only ~80% occupied => a MIXTURE, so no single
#                  model is right for all molecules. It sits 3 A from the site-B footprint, which
#                  makes it the most consequential one to get wrong; it gets its own arm, not a
#                  silent inclusion here.
SEQUONS_PARTIAL = [32, 579]

# N-glycan core linkages. (donor_residue_index, donor_atom, acceptor_residue_index, acceptor_atom)
# Asn ND2 -> NAG1 C1 is the protein link; the rest build Man3GlcNAc2.
CORE5 = [("NAG", None), ("NAG", (1, "O4")), ("BMA", (2, "O4")),
         ("MAN", (3, "O3")), ("MAN", (3, "O6"))]

# PDB has a ONE-CHARACTER chain column and the runner writes --output_format pdb, so stay inside a
# single-character pool rather than assuming the writer invents something sane.
POOL = [c for c in "CDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"]


def load_target(tsv):
    """The constant suffix of `full_sequence`. ⛔ ASSERTED unique, not assumed."""
    suf = collections.Counter()
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        fs = (r.get("full_sequence") or "").strip().upper()
        sq = (r.get("sequence") or "").strip().upper().split(":")[0]
        if fs and sq and fs.startswith(sq):
            suf[fs[len(sq):]] += 1
    if not suf:
        sys.exit("⛔ no design whose full_sequence starts with its sequence -- cannot derive target")
    (tgt, n), = suf.most_common(1)
    if len(suf) != 1:
        sys.exit(f"⛔ {len(suf)} DISTINCT constant suffixes, not 1 -- the assay target is not "
                 f"unique and must not be guessed: {[(len(k), v) for k, v in suf.most_common()]}")
    print(f"✅ assay target derived from {n} designs: {len(tgt)} aa, one variant only")
    return tgt


def check_pdb_agreement(tgt, pdb):
    """⛔ The structure and the assay construct must be the same protein in the same numbering."""
    three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
             "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
             "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    d = {}
    for ln in open(pdb):
        if ln.startswith("ATOM"):
            d[int(ln[22:26])] = three.get(ln[17:20].strip(), "X")
    mis = [(k, v, tgt[k - 1]) for k, v in sorted(d.items()) if k <= len(tgt) and tgt[k - 1] != v]
    if mis:
        sys.exit(f"⛔ {len(mis)} residue(s) disagree between {os.path.basename(pdb)} and the assay "
                 f"construct, first {mis[0]} -- refusing to fold against a mismatched numbering")
    print(f"✅ {len(d)} structure residues agree with the assay construct at the same WT index")


def glycan_block(tgt, mode, sequons):
    """(ligand yaml lines, constraint yaml lines, n_sugars). Asserts every sequon is ASN."""
    if mode == "none":
        return [], [], 0
    bad = [k for k in sequons if tgt[k - 1] != "N"]
    if bad:
        sys.exit(f"⛔ sequon(s) {bad} are not ASN in the assay construct -- wrong numbering")
    tree = [("NAG", None)] if mode == "nag" else CORE5
    lig, con, ids, n = [], [], iter(POOL), 0
    for k in sequons:
        local = []
        for resname, link in tree:
            try:
                cid = next(ids)
            except StopIteration:
                sys.exit("⛔ out of single-character chain ids -- PDB output cannot hold this many "
                         "glycan chains; reduce the tree or switch the runner to mmCIF")
            local.append(cid)
            lig += [f"- ligand:", f"    id: {cid}", f"    ccd: {resname}"]
            n += 1
            if link is None:
                con += [f"- bond:", f"    atom1: [B, {k}, ND2]", f"    atom2: [{cid}, 1, C1]"]
            else:
                parent, atom = link
                con += [f"- bond:", f"    atom1: [{local[parent - 1]}, 1, {atom}]",
                        f"    atom2: [{cid}, 1, C1]"]
    return lig, con, n


def read_binders(tsv, minlen, maxlen):
    out, skipped = [], []
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        name = re.sub(r"[^A-Za-z0-9_.-]", "_", (r.get("name") or "").strip())
        seq = (r.get("sequence") or "").strip().upper().split(":")[0]
        if not name or not seq:
            skipped.append((name, "empty")); continue
        if set(seq) - AA:
            skipped.append((name, f"nonstandard {sorted(set(seq)-AA)}")); continue
        if not (minlen <= len(seq) <= maxlen):
            skipped.append((name, f"len {len(seq)}")); continue
        out.append((name, seq))
    return out, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="/PATH/TO/assay_table.tsv",
                    help="the competition's released table; needs a full_sequence column")
    ap.add_argument("--pdb", default="/PATH/TO/"
                                    "EGFR_receptor_glycosylated_WTnumbering_tethered.pdb")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--glycans", choices=("none", "nag", "core5"), required=True)
    # ⛔⛔ THE TEMPLATE IS NOT OPTIONAL POLISH -- WITHOUT IT BOLTZ BUILDS THE WRONG STATE.
    # MEASURED 2026-10-03: untemplated, Boltz folds the 621-mer with domains I and III CLAMPED
    # TOGETHER (H409<->R29 CA-CA 10.4 A vs 29.0 A in the tethered crystal, consistent across 12
    # independent folds, spread only 3.3 A). Domain I then sits 1.8-4.4 A from EVERY site-A binder
    # and domain III from every site-B binder -- so ipTM/pdockq/contacts are measured against a
    # complex that cannot exist in the tethered assay target. ⛔ A SOFT template does NOT fix it
    # (tgt_fit 17.5 A, no better than none). ⭐ force+threshold+--use_potentials DOES:
    # tgt_fit 0.69 A, H409<->R29 28.49 A, zero intrusion, binder unharmed (sc_rmsd 0.40 vs 0.38).
    # ⚠️ Costs ~1.36x wall time and leaves ~6 mild backbone stretches (4.2-4.9 A) in the target --
    # all 33-83 A from the binder, i.e. harmless here. Known upstream, boltz issue #541, closed
    # "not planned". ⛔ template_id is boltz's SUBCHAIN name ("Axp"), NOT the author chain "A".
    ap.add_argument("--template", default="",
                    help="mmCIF of the TETHERED receptor. Must be gemmi-written with the entity "
                         "sequence spanning the FULL 621-mer, or boltz's parse_polymer raises.")
    ap.add_argument("--template-chain", default="Axp",
                    help="chain name INSIDE the template as boltz parses it (subchain id)")
    ap.add_argument("--template-threshold", type=float, default=2.0,
                    help="A of allowed deviation; force is always on when --template is given")
    ap.add_argument("--msa", default="/workspace/msa/target.a3m",
                    help="worker-absolute a3m for chain B. ⛔ msa:empty on a REAL protein is REFUTED.")
    ap.add_argument("--per-shard", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-len", type=int, default=25)
    ap.add_argument("--max-len", type=int, default=300)
    ap.add_argument("--include-partial-sequons", action="store_true",
                    help="also glycosylate WT 32 and 579 (partial occupancy -- a MIXTURE)")
    a = ap.parse_args()

    tgt = load_target(a.tsv)
    check_pdb_agreement(tgt, a.pdb)
    sequons = SEQUONS_FULL + (SEQUONS_PARTIAL if a.include_partial_sequons else [])
    lig, con, nsug = glycan_block(tgt, a.glycans, sorted(sequons))

    binders, skipped = read_binders(a.tsv, a.min_len, a.max_len)
    if a.limit:
        binders = binders[:a.limit]
    if not binders:
        sys.exit("⛔ no binder sequences survived filtering")

    per = a.per_shard or len(binders)
    nsh = (len(binders) + per - 1) // per
    for i, (name, seq) in enumerate(binders):
        d = os.path.join(a.outdir, f"shard{i // per:02d}")
        os.makedirs(d, exist_ok=True)
        body = ["version: 1", "sequences:",
                "- protein:", "    id: A", f"    sequence: {seq}", "    msa: empty",
                "- protein:", "    id: B", f"    sequence: {tgt}", f"    msa: {a.msa}"]
        body += ["  " + l if False else l for l in lig]
        if con:
            body += ["constraints:"] + con
        if a.template:
            body += ["templates:",
                     f"    - cif: {a.template}",
                     "      chain_id: B",
                     f"      template_id: {a.template_chain}",
                     "      force: true",
                     f"      threshold: {a.template_threshold}"]
        with open(os.path.join(d, f"{name}.yaml"), "w") as fh:
            fh.write("\n".join(body) + "\n")

    L = sorted(len(s) for _, s in binders)
    print(f"\n⭐ {len(binders)} YAMLs over {nsh} shard(s) -> {a.outdir}")
    print(f"   binder {L[0]}-{L[-1]} aa (median {L[len(L)//2]})   target {len(tgt)} aa")
    print(f"   glycans: {a.glycans}  ({nsug} sugar residues over {len(sequons)} sequons "
          f"{sorted(sequons)})")
    print(f"   chain B MSA: {a.msa}")
    if a.template:
        print(f"   ⭐ TEMPLATE (forced, threshold {a.template_threshold} A): {a.template}")
        print(f"      chain_id B <- template chain {a.template_chain}"
              f"   ⛔ the worker MUST pass --use_potentials or force is ignored")
    else:
        print("   ⛔⛔ NO TEMPLATE: boltz will build the EXTENDED state (domains I/III clamped, "
              "H409<->R29 10.4 A vs 29.0 A tethered). Interface metrics will be INVALID.")
    print(f"   approx tokens/complex: {len(tgt)} + {L[len(L)//2]} protein + ~{nsug*13} ligand atoms")
    if skipped:
        print(f"   skipped {len(skipped)}: {skipped[:4]}")
    if a.glycans != "none":
        print("   ⚠️ glycan trees UNDER-represent native HEK293 glycans (15-20 residues, "
              "sialylated): a LOWER bound on interference.")


if __name__ == "__main__":
    main()
