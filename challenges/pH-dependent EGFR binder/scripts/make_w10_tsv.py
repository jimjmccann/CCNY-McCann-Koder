#!/usr/bin/env python3
"""Build the 100-sequence TSV for the egfrW10 arm: the WINNING site-A backbone
(arm_shard11_shard11_24_model_2) at 10 draws.

Source: the SAME SolubleMPNN fasta egfrP14 and egfrA21 drew from (ids 1-100), so the
existing 3-draw results map onto these sequences by id with no re-derivation.

⛔ Sequence length is ASSERTED against the skeleton's chain A. A length mismatch would
mean a different backbone's sequences and every downstream contact number would be wrong.

usage: make_w10_tsv.py <fasta> <skeleton.yaml> <out.tsv>
"""
import re, sys

fa, skel, out = sys.argv[1], sys.argv[2], sys.argv[3]
DESIGN = "arm_shard11_shard11_24_model_2"

# --- skeleton chain A length, to assert against ---
lines = open(skel).read().splitlines()
want = None
for i, ln in enumerate(lines):
    if ln.strip() == "id: A":
        for j in range(i, min(i + 4, len(lines))):
            if lines[j].strip().startswith("sequence:"):
                want = len(lines[j].split("sequence:", 1)[1].strip())
if want is None:
    sys.exit("⛔ could not read chain A length from the skeleton")
print(f"skeleton chain A: {want} aa")

# --- parse the fasta: only records carrying an explicit id= ---
recs, cur = {}, None
for ln in open(fa):
    ln = ln.rstrip("\n")
    if ln.startswith(">"):
        m = re.search(r"\bid=(\d+)", ln)
        cur = int(m.group(1)) if m else None
    elif cur is not None and ln.strip():
        # ⛔ LigandMPNN writes BINDER:TARGET_CONTEXT on one line, ':'-separated. Chain A is
        # field 0 (75 aa); the 231-aa remainder is the RFd3 sphere-cut target CONTEXT and is
        # irrelevant here -- chain B comes from the skeleton (full ecto621, 621 aa). Taking
        # the whole line would write a 307-aa chain A and every contact number would be wrong.
        recs.setdefault(cur, "")
        recs[cur] += ln.strip().split(":")[0]

print(f"parsed {len(recs)} numbered sequences, ids {min(recs)}-{max(recs)}")
missing = [i for i in range(1, 101) if i not in recs]
if missing:
    sys.exit(f"⛔ ids absent from the fasta: {missing}")

bad = {i: len(s) for i, s in recs.items() if len(s) != want}
if bad:
    sys.exit(f"⛔ LENGTH MISMATCH vs skeleton chain A ({want} aa): {bad}")
print(f"✅ all 100 sequences are {want} aa, matching the skeleton's chain A")

aa = set("ACDEFGHIKLMNPQRSTVWY")
for i, s in recs.items():
    odd = set(s) - aa
    if odd:
        sys.exit(f"⛔ id{i} has non-standard residues {sorted(odd)}")
print("✅ all sequences are standard residues only")

# --- CROSS-CHECK against the egfrA21 input, which actually folded and landed -------------
# If our extraction matches that TSV for the ids it contains, the extraction is provably the
# same as the one that produced the 21-draw structures. This is cheap and it is the only
# independent check available.
import os
A21 = "/PATH/TO/SCRATCH/siteA_top3x10.tsv"
if os.path.exists(A21):
    n_ck = 0
    for ln in open(A21):
        name, seq = ln.rstrip("\n").split("\t")
        if not name.startswith(f"A_{DESIGN}_id"):
            continue
        i = int(name.rsplit("_id", 1)[1])
        if recs.get(i) != seq:
            sys.exit(f"⛔⛔ MISMATCH vs the PROVEN egfrA21 input at id{i}:\n"
                     f"   ours {recs.get(i)}\n   A21  {seq}")
        n_ck += 1
    print(f"✅ cross-check: {n_ck} ids match the egfrA21 TSV byte-for-byte "
          f"(the input that actually folded)")
else:
    print(f"⚠️  egfrA21 TSV absent at {A21} -- cross-check SKIPPED, not passed")

with open(out, "w") as fh:
    for i in range(1, 101):
        fh.write(f"A_{DESIGN}_id{i}\t{recs[i]}\n")
print(f"✅ wrote {out}: 100 rows")
print(f"   uniq sequences: {len({recs[i] for i in range(1,101)})} of 100")
