#!/usr/bin/env bash
# Lane ordstarts, build round: dev/ordstarts-after.R on the eight
# random-effect designs (40 seeds), the four designs without random
# effects (40 seeds) and the K = 3 design (20 seeds), in two queues (two
# R processes at once).
#   bash dev/ordstarts-after.sh
cd "$(dirname "$0")/.." || exit 1
export NOT_CRAN=true
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/ordstarts-out
one() {
  local tag=$1
  shift
  rm -f "$O/$tag.tsv"
  local t0=$(date +%s)
  "$R" dev/ordstarts-after.R "$@" "$O/$tag.tsv" > "$O/$tag.log" 2>&1
  echo "JOB $tag exit $? wall $(( $(date +%s) - t0 ))s" >> "$O/after-jobs.txt"
}
q1() {
  one after-cum_sratio-int-1000 cum_sratio int 1000 50 1:40
  one after-cum_cum-int-1000 cum_cum int 1000 50 1:40
  one after-cum_sratio-int-300 cum_sratio int 300 20 1:40
  one after-cum_cum-int-300 cum_cum int 300 20 1:40
  one after-cum_sratio-none-1000 cum_sratio none 1000 1 1:40
  one after-cum_cum-none-1000 cum_cum none 1000 1 1:40
  one after-k3-none-3000 k3 none 3000 50 1:20
}
q2() {
  one after-cum_sratio-slope-1000 cum_sratio slope 1000 50 1:40
  one after-cum_cum-slope-1000 cum_cum slope 1000 50 1:40
  one after-cum_sratio-slope-300 cum_sratio slope 300 20 1:40
  one after-cum_cum-slope-300 cum_cum slope 300 20 1:40
  one after-cum_sratio-none-300 cum_sratio none 300 1 1:40
  one after-cum_cum-none-300 cum_cum none 300 1 1:40
}
rm -f "$O/after-jobs.txt"
q1 &
q2 &
wait
echo "ALL DONE" >> "$O/after-jobs.txt"
