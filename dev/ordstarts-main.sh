#!/usr/bin/env bash
# Lane ordstarts: the main measurement, eight designs of 40 seeds in two
# queues (two R processes at once), each design written to
# dev/ordstarts-out/main-<family>-<re>-<n>.tsv by dev/ordstarts-run.R.
#   bash dev/ordstarts-main.sh
cd "$(dirname "$0")/.." || exit 1
export NOT_CRAN=true
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/ordstarts-out
one() {
  local fam=$1 re=$2 n=$3 G=$4
  local tag="main-$fam-$re-$n"
  rm -f "$O/$tag.tsv"
  local t0=$(date +%s)
  "$R" dev/ordstarts-run.R "$fam" "$re" "$n" "$G" 1:40 "$O/$tag.tsv" \
    > "$O/$tag.log" 2>&1
  echo "JOB $tag exit $? wall $(( $(date +%s) - t0 ))s" >> "$O/main-jobs.txt"
}
q1() {
  one cum_cum int 1000 50
  one cum_sratio slope 1000 50
  one cum_cum int 300 20
  one cum_sratio slope 300 20
}
q2() {
  one cum_sratio int 1000 50
  one cum_cum slope 1000 50
  one cum_sratio int 300 20
  one cum_cum slope 300 20
}
rm -f "$O/main-jobs.txt"
q1 &
q2 &
wait
echo "ALL DONE" >> "$O/main-jobs.txt"
