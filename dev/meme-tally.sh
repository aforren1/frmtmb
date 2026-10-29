#!/bin/bash
# Lane meme: per-package totals of a run-par.sh tier from its per-file
# logs. Usage: bash dev/meme-tally.sh <tier files dir>
D="$1"
for f in "$D"/*.txt; do
  p=$(basename "$f" | sed 's/--.*//')
  r=$(grep -a "^RESULT" "$f" | tail -1)
  [ -z "$r" ] && { echo "$p NORESULT $(basename $f)"; continue; }
  echo "$p $r"
done | awk '{p=$1; if ($2=="NORESULT") {nr[p]++; next}
  for (i=3;i<=NF;i++) {split($i,a,"="); s[p,a[1]]+=a[2]}; n[p]++}
  END {for (p in n) printf "%s files=%d pass=%d fail=%d err=%d skip=%d warn=%d noresult=%d\n",
    p, n[p], s[p,"pass"], s[p,"fail"], s[p,"err"], s[p,"skip"], s[p,"warn"], nr[p]}' | sort
