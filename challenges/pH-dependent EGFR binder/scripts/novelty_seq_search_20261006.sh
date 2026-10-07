#!/usr/bin/env bash
# Sequence-half novelty measurement for the Adaptyv/Proteinbase Level-3/4 gate.
# Reproduces the MMseqs2-vs-SwissProt(+PDB) part of what Proteinbase runs.
# 2026-10-06.  NO 2>/dev/null anywhere: stderr is always shown.
set -euo pipefail

# $DESIGN_DATA_ROOT must point at a volume with room for the sequence databases
# (Swiss-Prot + pdb_seqres, a few GB uncompressed) and the raw hit tables.
# Set it however you like before running.
: "${DESIGN_DATA_ROOT:?set DESIGN_DATA_ROOT to a scratch volume with room for the sequence DBs}"
DB="$DESIGN_DATA_ROOT/novelty_db_20261006"
WORK="$DESIGN_DATA_ROOT/novelty_search_20261006"
REPO=/PATH/TO/CHECKOUT
Q="$REPO/data/novelty_20261006/novelty_seq_query_20261006.fasta"
OUT="$REPO/data/novelty_20261006"
mkdir -p "$WORK" "$OUT"

export PATH="/PATH/TO/mmseqs2-env/bin:$PATH"   # MMseqs2 18.8cc5c (bioconda)
echo "### mmseqs version"; mmseqs version

# ---- protein-only PDB seqres --------------------------------------------
if [ ! -s "$DB/pdb_seqres_prot.fasta" ]; then
  echo "### filtering pdb_seqres to mol:protein"
  awk '/^>/{keep=/mol:protein/} keep' "$DB/pdb_seqres.txt" > "$DB/pdb_seqres_prot.fasta"
fi
echo "### pdb protein chains: $(grep -c '^>' "$DB/pdb_seqres_prot.fasta")"
echo "### swissprot seqs:     $(grep -c '^>' "$DB/uniprot_sprot.fasta")"

# ---- search -------------------------------------------------------------
# fident = nident/alnlen  (LOCAL identity, mmseqs native column)
# we also emit nident + qlen so GLOBAL identity = nident/qlen is computable.
# qaln/taln are emitted so identity can be RECOMPUTED from the alignment and
# checked against mmseqs' own columns (see the analyzer's verification pass).
#
# ⛔ `-a` (backtrace) IS MANDATORY.  Without it mmseqs SILENTLY reports
#    fident=0.000 and nident=0 for every hit instead of erroring, which reads
#    as "0% identity to everything".  Measured 2026-10-06.
FMT='query,target,pident,fident,nident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits,qaln,taln,theader'

run_search () {   # $1=name  $2=db fasta
  local name="$1" dbf="$2"
  echo "### mmseqs easy-search vs $name  (-s 7.5, -e 10000, --max-seqs 5000, --alignment-mode 3, -a, -c 0, --min-seq-id 0)"
  rm -rf "$WORK/tmp_$name"
  local full="$WORK/mmseqs_${name}_hits_full_20261006.tsv"
  mmseqs easy-search "$Q" "$dbf" "$full" "$WORK/tmp_$name" \
    -s 7.5 -e 10000 --max-seqs 5000 --alignment-mode 3 -a -c 0 --min-seq-id 0 \
    --max-seq-len 100000 --threads 8 --format-mode 4 --format-output "$FMT"
  # repo copy drops the two bulky alignment-string columns (15,16); nident +
  # alnlen + qlen fully determine both identity conventions anyway.
  cut --complement -f15,16 "$full" > "$OUT/mmseqs_${name}_hits_20261006.tsv"
  echo "### $name hit rows: $(( $(wc -l < "$full") - 1 ))   (full+aln: $full)"
  echo "### $name sanity: rows with nident==0 ->" \
       "$(awk -F'\t' 'NR>1 && $5==0' "$full" | wc -l)"
}

run_search swissprot "$DB/uniprot_sprot.fasta"
run_search pdbseqres "$DB/pdb_seqres_prot.fasta"
echo "### DONE"
