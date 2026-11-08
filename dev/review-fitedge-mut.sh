#!/usr/bin/env bash
# Review of lane fitedge: mutants, one test file per R process, sequentially.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-fitedge
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$WT/dev/review-fitedge-out"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_PARALLEL_PORT=11834
run() {
  "$RS" "$WT/dev/review-fitedge-mutrun.R" "$1" "$2" > "$OUT/mut-$1-$2.txt" 2>&1
  grep -h "^lib\|^RESULT\|^DETAIL" "$OUT/mut-$1-$2.txt"
}
run nogn test-gp-start.R
run gn1 test-gp-start.R
run noothers test-gp-start.R
run kappa12 test-gp-start.R
run nokappa test-gp-start.R
run novalue test-gp-start.R
run miback test-latres.R
run minever test-latres.R
run noridge test-nl-family.R
run ridgekept test-nl-family.R
echo MUT-DONE
