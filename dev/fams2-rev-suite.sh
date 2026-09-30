#!/bin/bash
# Reviewer: run a list of test files, one per R process, P at a time.
#   bash dev/fams2-rev-suite.sh <arm> <pkg> <outdir> <gated> <P> file...
arm=$1; pkg=$2; out=$3; gated=$4; P=$5; shift 5
cd /c/Users/adf44/source/r/frmtmb-wt-fams2
mkdir -p "$out"
export TMP="C:\Users\adf44\AppData\Local\Temp\1" TEMP="C:\Users\adf44\AppData\Local\Temp\1"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export FRMTMB_BRMS_FIT_TESTS=$gated
n=0
for b in "$@"; do
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/fams2-rev-run-test.R "$arm" "$pkg" "$b" > "$out/$b.log" 2>&1 &
  n=$((n + 1))
  if [ $n -ge $P ]; then wait -n; n=$((n - 1)); fi
done
wait
grep -h "^RESULT" "$out"/*.log > "$out/RESULTS.txt"
echo "done $(wc -l < $out/RESULTS.txt) of $#"
