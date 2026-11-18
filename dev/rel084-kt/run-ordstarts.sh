#!/usr/bin/env bash
# rel084-kt: ordstarts key tables on rellib-r23. The lane's
# dev/ordstarts-after.sh ran two queues at once; here the same 13 jobs run
# in one queue (one R process at a time). Output: dev/rel084-log/kt/ordstarts-out/.
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
O=dev/rel084-log/kt/ordstarts-out
mkdir -p $O
one() {
  local tag=$1; shift
  rm -f "$O/$tag.tsv"
  local t0=$(date +%s)
  "$RS" $K/ordstarts-after.R "$@" "$O/$tag.tsv" > "$O/$tag.log" 2>&1
  echo "JOB $tag exit $? wall $(( $(date +%s) - t0 ))s" >> "$O/after-jobs.txt"
}
rm -f "$O/after-jobs.txt"
one after-cum_sratio-int-1000 cum_sratio int 1000 50 1:40
one after-cum_cum-int-1000 cum_cum int 1000 50 1:40
one after-cum_sratio-int-300 cum_sratio int 300 20 1:40
one after-cum_cum-int-300 cum_cum int 300 20 1:40
one after-cum_sratio-none-1000 cum_sratio none 1000 1 1:40
one after-cum_cum-none-1000 cum_cum none 1000 1 1:40
one after-k3-none-3000 k3 none 3000 50 1:20
one after-cum_sratio-slope-1000 cum_sratio slope 1000 50 1:40
one after-cum_cum-slope-1000 cum_cum slope 1000 50 1:40
one after-cum_sratio-slope-300 cum_sratio slope 300 20 1:40
one after-cum_cum-slope-300 cum_cum slope 300 20 1:40
one after-cum_sratio-none-300 cum_sratio none 300 1 1:40
one after-cum_cum-none-300 cum_cum none 300 1 1:40
echo "ALL DONE" >> "$O/after-jobs.txt"
krun ordstarts-sum6.log $K/ordstarts-sum6.R
krun ordstarts-ulps2-weights.log $K/ordstarts-ulps2.R weights
krun ordstarts-short-p9.log $K/ordstarts-short-p9.R
echo ALLDONE > $KLOG/ordstarts-done.txt
