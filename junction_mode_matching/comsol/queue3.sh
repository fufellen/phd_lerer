#!/bin/bash
CB="/c/Program Files/COMSOL/COMSOL62/Multiphysics/bin/win64/comsolbatch.exe"
cd /c/workspace/COMSOL/queue_2026/R_reverse || exit 1
run_case () {
  local name="$1"; shift
  if tasklist 2>/dev/null | grep -qi comsolbatch; then echo "$(date +%T) another comsolbatch is running, abort" >> queue.log; exit 2; fi
  echo "$(date +%T) START $name" >> queue.log
  env "$@" timeout 5400 "$CB" -inputfile RunReverse.class -batchlog "log_$name.log" -nosave > "out_$name.txt" 2>&1
  echo "$(date +%T) END $name exit=$?" >> queue.log
}
run_case abrupt_rev_a_lp500 MR_ONLY=abrupt MR_DIR=rev MR_PCM=2.712 MR_LPCM=5.0
run_case sync405_rev_a_lp500 MR_ONLY=sync405 MR_DIR=rev MR_PCM=2.712 MR_LPCM=5.0
echo "$(date +%T) QUEUE3 DONE" >> queue.log
