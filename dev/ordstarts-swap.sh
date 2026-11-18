#!/usr/bin/env bash
# Lane ordstarts fix round: dev/ordstarts-swap.R on the eight K = 2
# designs (40 seeds) and the k3 design (20 seeds, lean), in two queues
# (two R processes at once).
#   bash dev/ordstarts-swap.sh
cd "$(dirname "$0")/.." || exit 1
export NOT_CRAN=true
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/ordstarts-out
one() {
  local tag=$1
  shift
  rm -f "$O/$tag.tsv"
  local t0=$(date +%s)
  "$R" dev/ordstarts-swap.R "$@" "$O/$tag.tsv" $LEAN > "$O/$tag.log" 2>&1
  echo "JOB $tag exit $? wall $(( $(date +%s) - t0 ))s" >> "$O/swap-jobs.txt"
}
q1() {
  one swap-cum_cum-int-1000 cum_cum int 1000 50 1:40
  one swap-cum_sratio-int-1000 cum_sratio int 1000 50 1:40
  one swap-cum_cum-int-300 cum_cum int 300 20 1:40
  one swap-cum_sratio-int-300 cum_sratio int 300 20 1:40
  LEAN=lean one swapk3-a k3 int 1000 50 1:10
}
q2() {
  one swap-cum_cum-slope-1000 cum_cum slope 1000 50 1:40
  one swap-cum_sratio-slope-1000 cum_sratio slope 1000 50 1:40
  one swap-cum_cum-slope-300 cum_cum slope 300 20 1:40
  one swap-cum_sratio-slope-300 cum_sratio slope 300 20 1:40
  LEAN=lean one swapk3-b k3 int 1000 50 11:20
}
rm -f "$O/swap-jobs.txt"
q1 &
q2 &
wait
echo "ALL DONE" >> "$O/swap-jobs.txt"
