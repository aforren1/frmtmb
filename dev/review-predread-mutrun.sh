#!/usr/bin/env bash
# Build one mutant, run the lane's audit on the shapes it reaches, then the
# lane's own test files with the mutant ahead of the lane's library.
#   bash dev/review-predread-mutrun.sh <name> <audit shapes, comma separated>
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
nm="$1"; shapes="$2"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
O=dev/review-predread-out/mut-$nm
mkdir -p "$O"
LIB=C:/Users/adf44/source/r/rev-predread-lib/$nm
"$RS" dev/review-predread-mutants.R "$nm" > "$O/build.txt" 2>&1
sed 's#dev/predread-log/audit-%s.tsv#dev/review-predread-out/audit-%s.tsv#' dev/predread-audit.R \
  > dev/review-predread-audit-copy.R.tmp.$nm
"$RS" dev/review-predread-audit-copy.R.tmp.$nm "$LIB" "mut-$nm" "$shapes" > "$O/audit.txt" 2>&1 &
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-predread/dev/stan-cache
export FRMTMB_RELLIB="$LIB;C:/Users/adf44/source/r/wt-predread-lib"
for f in tests/testthat/test-linpred-layers.R tests/testthat/test-linpred-layers-brms.R \
  tests/testthat/test-me-newdata.R tests/testthat/test-me-newdata-brms.R tests/testthat/test-resp-read.R; do
  "$RS" dev/review-predread-runone.R frmtmb "$f" "$LIB" > "$O/$(basename "$f" .R).txt" 2>&1 &
done
"$RS" dev/review-predread-runone.R frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-readers-resp-me.R "$LIB" \
  > "$O/test-readers-resp-me.txt" 2>&1 &
wait
rm -f dev/review-predread-audit-copy.R.tmp.$nm
echo "=== $nm" > "$O/summary.txt"
grep -h "^lib:\|^cells:" "$O/audit.txt" >> "$O/summary.txt"
grep -ah "^RESULT" "$O"/test-*.txt >> "$O/summary.txt"
cat "$O/summary.txt"
