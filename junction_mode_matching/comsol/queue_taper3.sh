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
T25="0.0000:0.0000;0.0049:0.2864;0.8997:0.3223;0.9789:0.5444;1.0000:1.0000"
T50="0.0000:0.0000;0.0331:0.2976;0.9658:0.3012;0.9872:0.7008;1.0000:1.0000"
T100="0.0000:0.0000;0.0700:0.1728;0.9835:0.3305;0.9860:0.5229;1.0000:1.0000"
run_case tp3_table_L50_fwd  MT_SHAPE=table MT_TABLE="$T50"  MT_SHAPE_L=0.5  MR_DIR=fwd MR_IMAGE=1
run_case tp3_table_L50_rev  MT_SHAPE=table MT_TABLE="$T50"  MT_SHAPE_L=0.5  MR_DIR=rev
run_case tp3_table_L25_fwd  MT_SHAPE=table MT_TABLE="$T25"  MT_SHAPE_L=0.25 MR_DIR=fwd
run_case tp3_table_L25_rev  MT_SHAPE=table MT_TABLE="$T25"  MT_SHAPE_L=0.25 MR_DIR=rev
run_case tp3_table_L100_fwd MT_SHAPE=table MT_TABLE="$T100" MT_SHAPE_L=1.0  MR_DIR=fwd
run_case tp3_table_L100_rev MT_SHAPE=table MT_TABLE="$T100" MT_SHAPE_L=1.0  MR_DIR=rev
echo "$(date +%T) TAPER QUEUE3 DONE" >> queue_taper.log
