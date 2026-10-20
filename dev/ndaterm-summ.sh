#!/usr/bin/env bash
# Lane ndaterm: rebuild RESULTS.txt of a dev/ndaterm-par.sh run from its
# per-file logs, and list every file that failed, errored or has no
# RESULT line, and every lib: line that is not the lane library.
#   bash dev/ndaterm-summ.sh <tag>
OUT=/c/Users/adf44/source/r/frmtmb-wt-ndaterm/dev/ndaterm-suite-$1
list="$OUT/jobs.txt"; res="$OUT/RESULTS.txt"; : > "$res"
while read -r p f; do
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" 2>/dev/null | tail -1)
  if [ -n "$hit" ]; then echo "$p $hit" >> "$res"
  else echo "$p NO RESULT $(basename "$f")" >> "$res"; fi
  grep -a "^lib:" "$o" 2>/dev/null | grep -v "wt-ndaterm-lib/frmtmb " |
    sed "s|^|$p WRONG LIB |" >> "$res"
done < "$list"
echo "ran $(grep -c RESULT "$res") of $(wc -l < "$list")" >> "$res"
awk '{for (i = 1; i <= NF; i++) { split($i, kv, "=");
       if (kv[1] in s) s[kv[1]] += kv[2]; else if (kv[2] != "") s[kv[1]] = kv[2] } }
     END { printf "TOTAL pass=%d fail=%d error=%d skip=%d warn=%d\n",
           s["pass"], s["fail"], s["error"], s["skip"], s["warn"] }' \
  <(grep RESULT "$res") >> "$res"
grep -E "fail=[1-9]|error=[1-9]|NO RESULT|WRONG LIB|^ran|^TOTAL" "$res"
