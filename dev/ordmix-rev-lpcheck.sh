#!/usr/bin/env bash
# Reviewer of lane ordmix: run cases of dev/ordmix-rev-lpcheck.R in
# parallel. Logs: dev/ordmix-rev-lpcheck-log/<case>.txt
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/ordmix-rev-lpcheck-log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export TMPDIR="$TMP"
mkdir -p "$OUT"
for c in "$@"; do
  ( "$RS" "$ROOT/dev/ordmix-rev-lpcheck.R" "$c" lane > "$OUT/$c.txt" 2> "$OUT/$c.err" ) &
done
wait
for c in "$@"; do
  grep -aE "^(LP|GRAD|GRADFRM) " "$OUT/$c.txt"
  grep -aE "^Error" "$OUT/$c.err" | sed "s/^/ERR $c : /"
done
