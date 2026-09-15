#!/bin/bash
cd /c/workspace/femwell_mm || exit 1
run () {
  echo "$(date +%T) START $*" >> mm_queue.log
  python -u mode_matching_junction.py "$@" > "run_$(echo $* | tr ' ' '_').txt" 2>&1
  echo "$(date +%T) END $* exit=$?" >> mm_queue.log
}
run sync405 a A 300
run sync405 a B 200
run sync405 a A 200 0.7
run sync405 c A 200
run shared450 a A 200
run abrupt a A 200
echo "$(date +%T) MM QUEUE DONE" >> mm_queue.log
