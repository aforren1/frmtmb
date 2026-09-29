#!/usr/bin/env bash
# Reviewer's rerun of the 40 files the worker lists, plus the guard files,
# ONE per R process. ARM=lane|base, JOBS=n, OUT=dir.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-thresrefit
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
ARM=${ARM:-lane}
OUT=${OUT:-"$ROOT/dev/thresrefit-rev-testlog-$ARM"}
JOBS=${JOBS:-5}
LIST=${LIST:-"$ROOT/dev/thresrefit-rev-files.txt"}
mkdir -p "$OUT"

n=0
while read -r f; do
  [ -z "$f" ] && continue
  case "$f" in \#*) continue ;; esac
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done
  n=$((n + 1))
  ( cd "$ROOT" && FRMTMB_LIB="$ARM" "$RS" dev/thresrefit-rev-runtest.R "$f" \
      > "$OUT/$f.log" 2>&1 ) &
done < "$LIST"
wait
echo "launched $n files  arm=$ARM"
grep -h "^REVRESULT" "$OUT"/*.log | sort
echo "--- files with no REVRESULT line (aborted) ---"
while read -r f; do
  [ -z "$f" ] && continue
  case "$f" in \#*) continue ;; esac
  grep -q "^REVRESULT" "$OUT/$f.log" 2>/dev/null || echo "$f"
done < "$LIST"
echo "--- totals ---"
grep -h "^REVRESULT" "$OUT"/*.log | awk '
  { for (i = 1; i <= NF; i++) {
      split($i, a, "=");
      if (a[1] == "pass") p += a[2];
      if (a[1] == "fail") fl += a[2];
      if (a[1] == "err") e += a[2];
      if (a[1] == "skip") s += a[2]; } ; c++ }
  END { print "files=" c " pass=" p " fail=" fl " err=" e " skip=" s }'
