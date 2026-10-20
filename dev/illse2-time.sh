#!/usr/bin/env bash
# Lane illse2, item 5: wall and CPU time of the wide design
# (dev/illse2-cost.R, no profiler), one fit at a time, arms interleaved,
# REPS rounds. Arms: r8 (0.70.0), base (0.71.0), lane with I2_RICH
# unset, lane with I2_RICH=single.
#   bash dev/illse2-time.sh <re 0|1> <n> <reps> <tag> [arms]
set -u
re="$1"; n="$2"; reps="$3"; tag="$4"; arms="${5:-r8 base lane}"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
D=/c/Users/adf44/source/r/frmtmb-wt-illse2
out="$D/dev/illse2-log/time-$tag-$re-$n.txt"
: > "$out"
for r in $(seq 1 "$reps"); do
  for a in $arms; do
    case "$a" in
      laneE) env ILLSE2_ARM=lane I2_RICH=single "$RS" "$D/dev/illse2-cost.R" "$re" "$n" 0 2>&1 | grep "^arm .* re " | sed "s/^/rep $r $a | /" >> "$out" ;;
      *) env ILLSE2_ARM="$a" "$RS" "$D/dev/illse2-cost.R" "$re" "$n" 0 2>&1 | grep "^arm .* re " | sed "s/^/rep $r $a | /" >> "$out" ;;
    esac
  done
done
cat "$out"
