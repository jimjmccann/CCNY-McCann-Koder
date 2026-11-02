#!/bin/sh

############################################################################
#
# Script for mirroring PDB FTP archive using rsync
#
############################################################################
#
# This is the wwPDB-provided template, https://files.wwpdb.org/pub/pdb/software/rsyncPDB.sh
# (author Thomas Solomon, 2002, last updated 2014-05-14), configured as its own instructions
# direct: the three local lines filled in, one SERVER/PORT pair uncommented, and exactly one
# rsync target uncommented.
#
# Local configuration, 2026-10-30:
#   - SERVER = PDBe (rsync.ebi.ac.uk), port 873. MEASURED: the RCSB server rsync.wwpdb.org::ftp
#     was UNREACHABLE from the machine used here -- rc=143 (timeout) on port 33444 AND on the
#     default port.
#     PDBe is one of the three servers the template itself offers, so this is a choice the
#     template sanctions, not a deviation from it.
#   - TARGET = the mmCIF coordinate tree only. The template calls it "Approximately 24 GB";
#     MEASURED with --dry-run --stats against the live server on 2026-10-30, the tree is
#     several times that. The template's figure is from 2014. Trust the measurement.
#   - MIRRORDIR must be on a volume with room for the whole tree, and MUST NOT be inside
#     a directory you cannot afford to have --delete applied to. See the defect note below.
#
# Bandwidth only. No GPU and no accelerator.
############################################################################

############################################################################
# DEFECT FOUND 2026-10-31, NOT FIXED. READ BEFORE ANY RE-MIRROR OR REUSE.
#
# LOGFILE (below) is set INSIDE MIRRORDIR, while the rsync at the bottom of this
# file runs --delete against that same directory. The remote tree has no `logs`
# entry, so --delete removes it: the run deletes its own log while writing it.
#
# MEASURED, not inferred: after the 2026-10-30 run completed (rc=0, finished
# ~12:05 on 10-31) the log was ABSENT from MIRRORDIR, and the directory held 0
# stray top-level files. So the exit status was the ONLY completion evidence that
# survived, and the mirror had to be verified against the filesystem instead:
# 87 GiB, 1,294 two-char hash dirs, 261,584 *.cif.gz, gzip -t OK on a spot check.
# The mirror itself is GOOD -- the defect cost the log, not the data.
# INFERRED (the --delete mechanism above): consistent with both observations, but
# not separately tested by re-running.
#
# The defect is in the wwPDB template AS VENDORED, not in our three edited lines.
# FIX is to put LOGFILE outside MIRRORDIR -- deliberately NOT applied here: this
# file is the record of a run that already happened, and editing it would make the
# record disagree with the run. Anyone re-mirroring should move LOGFILE first.
############################################################################

# ⚠️ SET THESE TWO BEFORE RUNNING. The absolute paths used for the run this file records
# were machine-specific and have been replaced with placeholders for release; the rest of
# the file is unchanged.
MIRRORDIR=/path/to/your/pdb_mmcif_mirror
LOGFILE=/path/to/your/pdb_mmcif_mirror/logs
# ^^ DEFECT, REPRODUCED AS IT WAS RUN: this path is inside MIRRORDIR, which the rsync below
# --deletes, so the run destroys its own log. See the header. Move it outside MIRRORDIR.
RSYNC=/usr/bin/rsync

SERVER=rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated     # PDBe server name
PORT=873                                                      # port PDBe server is using

mkdir -p "$MIRRORDIR"

############################################################################
# Rsync only the mmCIF format coordinates  /pub/pdb/data/structures/divided/mmCIF
############################################################################
${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/data/structures/divided/mmCIF/ $MIRRORDIR > $LOGFILE 2>/dev/null
