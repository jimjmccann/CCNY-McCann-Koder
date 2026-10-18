#!/usr/bin/env bash
# j496 stage 2 -- the remaining 99 Boltz folds, chained to start when stage 1's boltz exits.
# ⛔ Chained rather than run in parallel: each fold holds ~9.4 GiB and this card is 16 GiB, so two
#    at once would OOM -- and a Boltz OOM exits 0 and looks fine.
# ⛔ Every stage verifies an OUTPUT COUNT, never an exit status, for the same reason.
# ⭐ Run it as a visible, interruptible job -- never as a detached background process, so it
#    can always be stopped.
set -uo pipefail
# PLACEHOLDER PATH -- not a working default; set it for your own tree.
D=/path/to/the/fold/round
BOLTZ=/path/to/your/boltz
S1_PID="${1:?need the stage1 boltz PID}"

cd "$D" || exit 1

echo "[chain] waiting for stage1 boltz PID $S1_PID to exit"
while [ -d "/proc/$S1_PID" ]; do sleep 30; done
echo "[chain] stage1 process gone at $(date '+%H:%M:%S')"

N1=$(find stage1_out -name '*_model_0.pdb' 2>/dev/null | wc -l)
echo "[chain] stage1 OUTPUT COUNT = $N1 / 10 (exit status is not evidence: a Boltz OOM exits 0)"
if [ "$N1" -lt 10 ]; then
  echo "[chain] ⛔ stage1 is SHORT ($N1/10). NOT starting stage2 -- diagnose first."
  exit 1
fi

echo "[chain] starting stage2: 99 folds, ~1m55s each, about 3.2 h"
"$BOLTZ" predict stage2 --out_dir stage2_out \
  --diffusion_samples 1 --recycling_steps 3 --devices 1 --accelerator gpu \
  --seed 20261018 --no_kernels --output_format pdb --write_full_pae 2>&1 | tail -4

N2=$(find stage2_out -name '*_model_0.pdb' 2>/dev/null | wc -l)
NP=$(find probe_out -name '*_model_0.pdb' 2>/dev/null | wc -l)
TOT=$((N1 + N2 + NP))
echo "[chain] stage2 OUTPUT COUNT = $N2 / 99"
echo "[chain] TOTAL FOLDED = $TOT / 110  (probe $NP + stage1 $N1 + stage2 $N2)"
[ "$TOT" -eq 110 ] && echo "[chain] ✅ all 110 folded" || echo "[chain] ⛔ SHORT: $TOT of 110"
