#!/usr/bin/env bash
# s113 -- upload the option-③ designability 2x2 RFd3 bundle to B2 and VERIFY it round-trip.
#
# Usage: bash s113_upload_designability_2x2_bundle.sh [--upload]
#   Without --upload: pre-flight only. ⛔ NOTHING LEAVES THIS BOX.
#
# ⛔⛔ EVERY GUARD CHECKED HERE IS COPIED FROM `RFD3_START` IN `the pod launcher`,
#    i.e. from THE LAUNCHER THAT WILL RUN THIS, ⛔ not from a remembered list
#    [[check-a-bundle-against-its-runners-own-guards]]. A bundle that fails one of these does not
#    fail cheaply: the pod boots, dies, and `die()` SLEEPS 1800 s AT FULL BILLING RATE.
#      FATAL_NO_RUNNER     <- /workspace/in/<runner> missing or empty
#      FATAL_NO_ARM_SPEC   <- neither arm.json nor shards/ present
#      FATAL_NO_TARGET     <- no *.pdb at the bundle ROOT (glob, not recursive)
#      FATAL_SHA1_MISMATCH <- staged runner sha1 != RFD3_SCRIPT_SHA1 passed at launch
#
# ⛔ `rclone`'s exit status is ALWAYS 0 [[rclone-exit-status-is-always-zero]] ⇒ every check gates on
#    STDOUT. ⛔ B2 stores SHA1, not MD5. ⛔ There is NO rclone.conf on this host; remotes come from
#    env files sourced IN A SUBSHELL. ⛔ FILENAMES ONLY -- never read, print or quote a credential.
set -uo pipefail

BUNDLE="rfd3_c5_designability_2x2_v1"
SRC="/PATH/TO/bundles/${BUNDLE}"
RUNNER="prod_rfd3_multi.sh"
WANT_SHA1="623b3c69924b4e0afb245330d5eb2b6f39df10de"
CREDR="/PATH/TO/CREDENTIALS/<READ_CREDS_ENV>"
CREDW="/PATH/TO/CREDENTIALS/<WRITE_CREDS_ENV>"
REMOTE_PATH="OBJECT-STORE-BUCKET/bundles/${BUNDLE}"
DO_UPLOAD="${1:-}"

[ -d "$SRC" ] || { echo "FATAL: no bundle at $SRC"; exit 1; }

echo "=== s113 pre-flight: $BUNDLE ==="

# ---- guard 1: the runner exists, is non-empty, and its sha1 is the one we will pass ----
[ -s "$SRC/$RUNNER" ] || { echo "⛔ FATAL_NO_RUNNER would fire: $SRC/$RUNNER missing/empty"; exit 1; }
got=$(sha1sum "$SRC/$RUNNER" | cut -c1-40)
if [ "$got" != "$WANT_SHA1" ]; then
    echo "⛔ FATAL_SHA1_MISMATCH would fire: staged $got != $WANT_SHA1"
    echo "   ⇒ pass RFD3_SCRIPT_SHA1=$got at launch, or re-copy a reviewed runner."
    exit 1
fi
echo "  ✅ runner $RUNNER  sha1 $got  (matches RFD3_SCRIPT_SHA1 to pass)"

# ---- guard 2: the multi-shard spec shape ----
[ -d "$SRC/shards" ] || { echo "⛔ FATAL_NO_ARM_SPEC would fire: no shards/ and no arm.json"; exit 1; }
nj=$(find "$SRC/shards" -maxdepth 1 -name 'arm_shard*.json' | wc -l)
nf=$(find "$SRC/shards" -maxdepth 1 -name 'arm_shard*.flags' | wc -l)
[ "$nj" -gt 0 ] || { echo "⛔ FATAL: shards/ holds no arm_shardNN.json"; exit 1; }
echo "  ✅ shards/: $nj spec(s) + $nf sidecar(s)"
# ⛔ contiguously numbered from 00, or RFD3_SHARDS="0-3" expands onto a shard that is not there.
i=0
while [ "$i" -lt "$nj" ]; do
    p=$(printf '%s/shards/arm_shard%02d.json' "$SRC" "$i")
    [ -s "$p" ] || { echo "⛔ FATAL: $p missing -- shards must be contiguous from 00"; exit 1; }
    i=$((i + 1))
done
echo "  ✅ shard numbering contiguous 00..$(printf '%02d' $((nj - 1)))"

# ---- guard 3: at least one .pdb AT THE ROOT. ⛔ the launcher uses a GLOB, not a recursive find,
#      so a target pdb tucked inside a subdirectory reads as absent and dies FATAL_NO_TARGET.
npdb=$(find "$SRC" -maxdepth 1 -name '*.pdb' | wc -l)
[ "$npdb" -gt 0 ] || { echo "⛔ FATAL_NO_TARGET would fire: no *.pdb at the bundle ROOT"; exit 1; }
echo "  ✅ $npdb target pdb(s) at the bundle root"

# ---- guard 4: every spec's `input` is a BARE BASENAME and that file is present at the root ----
for p in "$SRC"/shards/arm_shard*.json; do
    inp=$(python3 -I -c 'import json,sys
d=json.load(open(sys.argv[1]))
for k,v in d.items(): print(v["input"])' "$p")
    case "$inp" in
        */*) echo "⛔ FATAL: $(basename "$p") input '$inp' is a PATH; must be a bare basename"; exit 1;;
    esac
    [ -s "$SRC/$inp" ] || { echo "⛔ FATAL: $(basename "$p") names input '$inp' which is not at the root"; exit 1; }
done
echo "  ✅ every spec's input is a bare basename present at the root"

# ---- guard 5: the mandatory engine patch ----
[ -s "$SRC/patch_held_motif_guard.py" ] || {
    echo "⛔ FATAL: patch_held_motif_guard.py absent -- the runner DIES if the patch refuses,"
    echo "   and it is NOT optional alongside a held motif."; exit 1; }
echo "  ✅ patch_held_motif_guard.py present"

nall=$(find "$SRC" -type f | wc -l)
echo "=== pre-flight PASSED: $nall files ==="

if [ "$DO_UPLOAD" != "--upload" ] && [ "$DO_UPLOAD" != "--verify" ]; then
    echo "⛔ NOT uploaded (pass --upload). Nothing has left this box."
    exit 0
fi

[ -f "$CREDW" ] || { echo "FATAL: no write env file $CREDW"; exit 1; }
[ -f "$CREDR" ] || { echo "FATAL: no read env file $CREDR"; exit 1; }

# ⛔⛔ REFUSE A NON-EMPTY PREFIX. A bundle prefix is written ONCE; overwriting one another round is
#    still using is how a pod silently generates against the wrong specs. Gate on STDOUT, never rc.
# ⭐ `--verify` re-runs only the round-trip checks below, for when the bytes are already up and it
#   was the CHECK that was wrong. ⛔ It never writes, so it cannot be a way around this refusal.
if [ "$DO_UPLOAD" = "--upload" ]; then
    pre=$( ( set +u; . "$CREDR"; rclone lsf "<REMOTE_RO>:${REMOTE_PATH}/" 2>/dev/null ) | wc -l )
    if [ "$pre" -ne 0 ]; then
        echo "⛔⛔ REFUSING: <REMOTE_RO>:${REMOTE_PATH}/ already holds $pre entries. Bump the version in"
        echo "   BUNDLE= here and in RFD3_BUNDLE_DIR at launch."
        exit 1
    fi
    echo "target prefix verified EMPTY (0 entries)"
    ( set +u; . "$CREDW"; rclone copy "$SRC" "<REMOTE_RW>:${REMOTE_PATH}/" )
else
    echo "⭐ --verify: round-trip checks only, no write."
fi

# ---- round-trip through the POD'S OWN bundles-scoped read key, gating on stdout ----
echo "=== round-trip via b2r (the key the pod itself uses) ==="
# ⛔ `rclone lsf -R` EMITS DIRECTORY ENTRIES TOO (a trailing "/"), so a naive line count compares
#   local FILES against remote files+dirs and reports a phantom surplus -- measured here on the
#   first run: 14 "objects" against 13 files, the difference being exactly `shards/`.
#   ⇒ filter them out and count FILES, which is what the pod actually needs present.
listing=$( ( set +u; . "$CREDR"; rclone lsf -R --files-only "<REMOTE_RO>:${REMOTE_PATH}/" 2>/dev/null ) | sort )
echo "$listing" | sed 's/^/  /'
nrem=$(echo "$listing" | grep -c .)
echo "  objects in B2: $nrem   local files: $nall"
[ "$nrem" -eq "$nall" ] || { echo "⛔⛔ FATAL: $nrem of $nall objects present -- do NOT launch"; exit 1; }

# ⛔ the runner is the one object whose bytes the pod HASHES, so compare its sha1 explicitly.
rem=$( ( set +u; . "$CREDR"; rclone hashsum sha1 "<REMOTE_RO>:${REMOTE_PATH}/${RUNNER}" 2>/dev/null ) | awk '{print $1}')
[ "$rem" = "$WANT_SHA1" ] || { echo "⛔⛔ FATAL: remote runner sha1 ${rem:0:12}... != ${WANT_SHA1:0:12}... -- do NOT launch"; exit 1; }
echo "  ✅ remote runner sha1 matches (the pod's FATAL_SHA1_MISMATCH guard will pass)"

# ⛔ and the shards/ subdirectory must have survived the copy, or FATAL_NO_ARM_SPEC fires.
echo "$listing" | grep -q '^shards/arm_shard00\.json$' || {
    echo "⛔⛔ FATAL: shards/ did not survive the upload -- FATAL_NO_ARM_SPEC would fire"; exit 1; }
echo "  ✅ shards/ subdirectory survived the upload"

echo "=== ✅ bundle $BUNDLE uploaded and VERIFIED round-trip ($nall objects). ==="
echo "⛔ NEXT, IN THIS ORDER:"
echo "   1. confirm the SHORT-deadline reaper from the J1054 fold round is GONE (END-ANCHORED pgrep)"
echo "   2. arm a NEW reaper sized to the RFd3 PESSIMISTIC WHOLE-POD wall, and verify its own log"
echo "   3. write the dated SPEND_LEDGER line naming the pod id -- 'No line, no campaign'"
echo "   4. THEN launch, with ALL SIX defaults overridden explicitly (POD_GPU especially)"
