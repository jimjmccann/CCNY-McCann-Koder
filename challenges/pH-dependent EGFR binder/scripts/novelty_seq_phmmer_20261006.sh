#!/usr/bin/env bash
# phmmer cross-check of the MMseqs2 sequence-novelty search.
# MMseqs2 is the PRIMARY (it is what Proteinbase runs); phmmer is an
# independent, profile-based second opinion on whether any SIGNIFICANT
# homolog exists at all.  2026-10-06.  No 2>/dev/null anywhere.
set -euo pipefail

# $DESIGN_DATA_ROOT: the same scratch volume the MMseqs2 script uses; it holds the
# sequence databases and the full (~25 MB per DB) phmmer tables.
: "${DESIGN_DATA_ROOT:?set DESIGN_DATA_ROOT to a scratch volume with room for the sequence DBs}"
DB="$DESIGN_DATA_ROOT/novelty_db_20261006"
WORK="$DESIGN_DATA_ROOT/novelty_search_20261006"
REPO=/PATH/TO/CHECKOUT
Q="$REPO/data/novelty_20261006/novelty_seq_query_20261006.fasta"
OUT="$REPO/data/novelty_20261006"
PH=phmmer   # a phmmer (HMMER 3.4+) on your PATH
mkdir -p "$WORK"

echo "### phmmer version"; "$PH" -h | sed -n '2p'

run_ph () {   # $1=name  $2=db fasta
  local name="$1" dbf="$2"
  # --max turns OFF all three heuristic filters (MSV/Vit/Fwd): maximal
  # sensitivity, which is the point -- a de novo sequence with no real homolog
  # must not be reported as "no hits" merely because a filter dropped it.
  echo "### phmmer --max vs $name  (-E 1000 --domE 1000 --incE 1000, 16 cpu)"
  "$PH" --max -E 1000 --domE 1000 --incE 1000 --incdomE 1000 \
        --cpu 16 --notextw \
        --tblout "$WORK/phmmer_${name}_tbl_20261006.txt" \
        --domtblout "$WORK/phmmer_${name}_domtbl_20261006.txt" \
        -o "$WORK/phmmer_${name}_full_20261006.txt" \
        "$Q" "$dbf"
  # repo copy keeps only the E<=10 rows: the full table is ~25 MB of mostly
  # noise (E up to 1000) and the full version stays under $DESIGN_DATA_ROOT.
  head -3 "$WORK/phmmer_${name}_domtbl_20261006.txt" \
      > "$OUT/phmmer_${name}_domtbl_E10_20261006.txt"
  awk '!/^#/ && $7+0<=10' "$WORK/phmmer_${name}_domtbl_20261006.txt" \
      >> "$OUT/phmmer_${name}_domtbl_E10_20261006.txt"
  echo "### $name domain rows: $(grep -vc '^#' "$WORK/phmmer_${name}_domtbl_20261006.txt")"
}

run_ph swissprot "$DB/uniprot_sprot.fasta"
run_ph pdbseqres "$DB/pdb_seqres_prot.fasta"
echo "### DONE"
