#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
"$RS" review-impsamp-fr1-fixtau.R
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/rev-impsamp-lib/fix_tau --no-multiarch --no-docs --no-test-load C:/Users/adf44/source/r/rev-impsamp-src/fr1-fix_tau 2>&1 | tail -1
REV_ARM=mutant REV_MUTLIB=C:/Users/adf44/source/r/rev-impsamp-lib/fix_tau "$RS" review-impsamp-fr1-ordfix.R 2>&1 | grep -v "^Boundary" > review-impsamp-out/fr1/ordfix-tau.txt &
"$RS" review-impsamp-fr1-mut.R test fix_tau > review-impsamp-out/fr1/mut/test-fix_tau.log 2>&1 &
wait
echo done
