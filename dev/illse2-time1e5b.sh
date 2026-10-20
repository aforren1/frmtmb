#!/usr/bin/env bash
# Lane illse2, item 5: the reviewer's n = 1e5 cost runs
# (dev/illse-rev-cost.R's design through dev/illse2-cost.R, no
# profiler), every arm at once so that they share the machine's load:
# 0.71.0 (base) and the lane, with and without (1 | g); rellib-r8
# (0.70.0) no longer exists.
#   bash dev/illse2-time1e5.sh <tag>
set -u
tag="$1"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
D=/c/Users/adf44/source/r/frmtmb-wt-illse2
for a in base lane; do
  for re in 0 1; do
    ( env ILLSE2_ARM=$a "$RS" "$D/dev/illse2-cost.R" $re 1e5 0 \
        > "$D/dev/illse2-log/t1e5-$tag-$a-$re.txt" 2>&1 ) &
  done
done
wait
grep -h "^arm .* re " "$D"/dev/illse2-log/t1e5-$tag-*.txt
