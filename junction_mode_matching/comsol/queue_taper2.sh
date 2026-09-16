#!/bin/bash
CB="/c/Program Files/COMSOL/COMSOL62/Multiphysics/bin/win64/comsolbatch.exe"
cd /c/workspace/COMSOL/queue_2026/R_reverse || exit 1
run_case () {
  local name="$1"; shift
  if tasklist 2>/dev/null | grep -qi comsolbatch; then echo "$(date +%T) another comsolbatch is running, abort" >> queue_taper.log; exit 2; fi
  echo "$(date +%T) START $name" >> queue_taper.log
  env "$@" timeout 3600 "$CB" -inputfile RunTaper.class -batchlog "log_$name.log" -nosave > "out_$name.txt" 2>&1
  echo "$(date +%T) END $name exit=$?" >> queue_taper.log
}
for sh in quad ellipsec expo sqrt ellipse klop step cubic; do
  img=""; if [ "$sh" = "ellipsec" ]; then img="MR_IMAGE=1"; fi
  run_case tp2_${sh}_L50_fwd  MT_SHAPE=$sh MT_SHAPE_L=0.5 MR_DIR=fwd $img
  run_case tp2_${sh}_L50_rev  MT_SHAPE=$sh MT_SHAPE_L=0.5 MR_DIR=rev $img
done
for sh in ellipsec expo quad linear; do
  run_case tp2_${sh}_L25_fwd MT_SHAPE=$sh MT_SHAPE_L=0.25 MR_DIR=fwd
  run_case tp2_${sh}_L25_rev MT_SHAPE=$sh MT_SHAPE_L=0.25 MR_DIR=rev
done
for sh in ellipsec expo linear; do
  run_case tp2_${sh}_L100_fwd MT_SHAPE=$sh MT_SHAPE_L=1.0 MR_DIR=fwd
  run_case tp2_${sh}_L100_rev MT_SHAPE=$sh MT_SHAPE_L=1.0 MR_DIR=rev
done
echo "$(date +%T) TAPER QUEUE2 DONE" >> queue_taper.log
