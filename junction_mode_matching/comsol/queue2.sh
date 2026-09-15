#!/bin/bash
# Second COMSOL queue: mesh refinement of the sync405 pair, started only after queue.sh has finished.
CB="/c/Program Files/COMSOL/COMSOL62/Multiphysics/bin/win64/comsolbatch.exe"
cd /c/workspace/COMSOL/queue_2026/R_reverse || exit 1
until grep -q "QUEUE DONE" queue.log 2>/dev/null; do sleep 30; done
run_case () {
  local name="$1"; shift
  if tasklist 2>/dev/null | grep -qi comsolbatch; then echo "$(date +%T) another comsolbatch is running, abort" >> queue.log; exit 2; fi
  echo "$(date +%T) START $name" >> queue.log
  env "$@" timeout 5400 "$CB" -inputfile RunReverse.class -batchlog "log_$name.log" -nosave > "out_$name.txt" 2>&1
  echo "$(date +%T) END $name exit=$?" >> queue.log
}
run_case sync405_fwd_a_m80 MR_ONLY=sync405 MR_DIR=fwd MR_PCM=2.712 MR_MESH=0.8
run_case sync405_fwd_a_m70 MR_ONLY=sync405 MR_DIR=fwd MR_PCM=2.712 MR_MESH=0.7
run_case sync405_rev_a_m70 MR_ONLY=sync405 MR_DIR=rev MR_PCM=2.712 MR_MESH=0.7
run_case sync405_fwd_c_m70 MR_ONLY=sync405 MR_DIR=fwd MR_PCM=3.308 MR_MESH=0.7
run_case sync405_rev_c_m70 MR_ONLY=sync405 MR_DIR=rev MR_PCM=3.308 MR_MESH=0.7
run_case abrupt_fwd_a_m70 MR_ONLY=abrupt MR_DIR=fwd MR_PCM=2.712 MR_MESH=0.7
run_case shared450_fwd_a_m70 MR_ONLY=shared450 MR_DIR=fwd MR_PCM=2.712 MR_MESH=0.7
echo "$(date +%T) QUEUE2 DONE" >> queue.log
