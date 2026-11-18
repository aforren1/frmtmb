#!/usr/bin/env bash
# Reviewer short check of lane ordstarts: two queues, one R process each.
#   bash dev/reviews/ordstarts-short-run.sh <a|b>
cd "$(dirname "$0")/../.." || exit 1
export NOT_CRAN=true
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/reviews
J=$O/ordstarts-short-jobs.txt
if [ "$1" = a ]; then
  for des in cs1000 cc300; do
    f=$O/ordstarts-short-rep-lane-$des.tsv
    rm -f "$f"
    "$R" $O/ordstarts-short-rep.R lane $des 1:10 "$f" > $O/ordstarts-short-rep-lane-$des.log 2>&1
    echo "JOB rep $des exit $?" >> $J
  done
  f=$O/ordstarts-short-p4-lane.tsv
  rm -f "$f"
  "$R" $O/ordstarts-short-p4.R lane compat2 400 1 "$f" > $O/ordstarts-short-p4-lane.log 2>&1
  for cfg in "cc 300" "cs 300" "cs 1000"; do
    set -- $cfg
    "$R" $O/ordstarts-short-p4.R lane $1 $2 1:40 "$f" >> $O/ordstarts-short-p4-lane.log 2>&1
  done
  echo "JOB p4 done" >> $J
  echo "QUEUE a done" >> $J
else
  for p in p1 p2 p6 p7 p8; do
    "$R" $O/ordstarts-short-$p.R > $O/ordstarts-short-$p.txt 2>&1
    echo "JOB $p exit $?" >> $J
  done
  echo "QUEUE b done" >> $J
fi
