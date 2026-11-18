#!/usr/bin/env bash
# Lane nlband: compare the survey's ROW lines between two arms, job by
# job, with the wall-time fields dropped (fresh, fresh2: field 7; scan:
# field 9; limit: field 5), and sum the COST lines.
# Usage: bash dev/nlband-accsum.sh <before tag> <after tag>
b="$1"; a="$2"
cd "$(dirname "$0")/nlband-out" || exit 1
drop() { # file, field to drop
  grep -a "^ROW " "$1" | awk -v k="$2" '{ $k = ""; print }'
}
tot_rows=0; tot_diff=0
for set in fresh fresh2 scan limit; do
  case $set in fresh|fresh2) k=7;; scan) k=9;; limit) k=5;; esac
  for fa in "$set-$a"/*; do
    fb="$set-$b/$(basename "$fa")"
    [ -f "$fb" ] || { echo "MISSING $fb"; continue; }
    n=$(drop "$fa" $k | wc -l)
    nb=$(drop "$fb" $k | wc -l)
    dd=$(diff <(drop "$fb" $k) <(drop "$fa" $k) | grep -c "^>")
    tot_rows=$((tot_rows + n)); tot_diff=$((tot_diff + dd))
    echo "JOB $set/$(basename "$fa") rows after $n before $nb differ $dd"
    if [ "$dd" -gt 0 ]; then
      diff <(drop "$fb" $k) <(drop "$fa" $k) | grep "^[<>]" | cut -c1-220 |
        sed 's/^/    /'
    fi
  done
done
echo "ALL rows $tot_rows differ $tot_diff"
for arm in "$b" "$a"; do
  cat fresh-"$arm"/* fresh2-"$arm"/* scan-"$arm"/* limit-"$arm"/* |
    grep -a "^COST" | awk -v arm="$arm" '
      { for (i = 1; i <= NF; i++) {
          if ($i == "escapes") e += $(i + 1)
          if ($i == "holes") h += $(i + 1)
          if ($i == "past") p += $(i + 1)
          if ($i == "evals") v += $(i + 1)
          if ($i == "escape") ve += $(i + 1)
          if ($i == "secs") s += $(i + 1)
        } n++ }
      END { printf "COSTSUM %s jobs %d escapes %d holes %d past %d evals %d in escape %d secs %.1f\n", arm, n, e, h, p, v, ve, s }'
done
