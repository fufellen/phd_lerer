#!/bin/bash
# Serialized COMSOL queue: one comsolbatch at a time, each bounded by a wall-clock limit.
CB="/c/Program Files/COMSOL/COMSOL62/Multiphysics/bin/win64/comsolbatch.exe"
cd /c/workspace/COMSOL/queue_2026/R_reverse || exit 1
run_case () {
  local name="$1"; shift
  if tasklist 2>/dev/null | grep -qi comsolbatch; then echo "$(date +%T) another comsolbatch is running, abort" >> queue.log; exit 2; fi
  echo "$(date +%T) START $name" >> queue.log
  env "$@" timeout 5400 "$CB" -inputfile RunReverse.class -batchlog "log_$name.log" -nosave > "out_$name.txt" 2>&1
  echo "$(date +%T) END $name exit=$?" >> queue.log
}
run_case sync405_fwd_a  MR_ONLY=sync405   MR_DIR=fwd MR_PCM=2.712 MR_IMAGE=1
run_case sync405_rev_a  MR_ONLY=sync405   MR_DIR=rev MR_PCM=2.712 MR_IMAGE=1
run_case sync405_fwd_c  MR_ONLY=sync405   MR_DIR=fwd MR_PCM=3.308
run_case sync405_rev_c  MR_ONLY=sync405   MR_DIR=rev MR_PCM=3.308 MR_IMAGE=1
run_case shared450_fwd_a MR_ONLY=shared450 MR_DIR=fwd MR_PCM=2.712
run_case shared450_rev_a MR_ONLY=shared450 MR_DIR=rev MR_PCM=2.712
run_case abrupt_fwd_a   MR_ONLY=abrupt    MR_DIR=fwd MR_PCM=2.712 MR_IMAGE=1
run_case abrupt_rev_a   MR_ONLY=abrupt    MR_DIR=rev MR_PCM=2.712 MR_IMAGE=1
run_case sync405_rev_a_m80 MR_ONLY=sync405 MR_DIR=rev MR_PCM=2.712 MR_MESH=0.8
run_case sync405_rev_a_lf250 MR_ONLY=sync405 MR_DIR=rev MR_PCM=2.712 MR_LFEED=2.5 MR_IMAGE=1
run_case sync405_fwd_a_lf250 MR_ONLY=sync405 MR_DIR=fwd MR_PCM=2.712 MR_LFEED=2.5
echo "$(date +%T) QUEUE DONE" >> queue.log
