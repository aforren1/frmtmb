#!/usr/bin/env bash
# Every mode of dev/samplegen-check.R, one process each, two at a time.
#   bash dev/samplegen-checkall.sh <LIB> <label>
LIB="$1"; LABEL="$2"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd /c/Users/adf44/source/r/frmtmb-wt-samplegen || exit 1
OUT=dev/samplegen-out/check-$LABEL
mkdir -p "$OUT"
set -- S T U N P Q G R I D
while [ $# -gt 0 ]; do
  "$RS" dev/samplegen-check.R "$LIB" "$1" > "$OUT/$1.txt" 2>&1 &
  if [ $# -gt 1 ]; then
    "$RS" dev/samplegen-check.R "$LIB" "$2" > "$OUT/$2.txt" 2>&1 &
    shift
  fi
  shift
  wait
done
echo DONE
