#!/bin/bash
# Run every test file of a package, one file per R process, P at a time.
#   bash dev/fams2-suite.sh <pkg> <outdir> [P] [gated] [pattern]
pkg=$1; out=$2; P=${3:-12}; gated=${4:-false}; pat=${5:-test-*.R}
cd /c/Users/adf44/source/r/frmtmb-wt-fams2
mkdir -p "$out"
if [ "$pkg" = "frmtmb" ]; then tdir=tests/testthat; else tdir=extensions/$pkg/tests/testthat; fi
export TMP="C:\Users\adf44\AppData\Local\Temp\1" TEMP="C:\Users\adf44\AppData\Local\Temp\1"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export FRMTMB_BRMS_FIT_TESTS=$gated
n=0
for f in $tdir/$pat; do
  b=$(basename "$f")
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/fams2-run-test.R "$pkg" "$b" > "$out/$b.log" 2>&1 &
  n=$((n + 1))
  if [ $n -ge $P ]; then wait -n; n=$((n - 1)); fi
done
wait
grep -h "^RESULT" "$out"/*.log > "$out/RESULTS.txt"
echo "files: $(ls $tdir/$pat | wc -l)  results: $(wc -l < $out/RESULTS.txt)"
