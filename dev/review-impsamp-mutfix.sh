#!/usr/bin/env bash
# the scale-free step: two edits, so this one is built by hand
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
L=review-impsamp-out/mutants
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
"$RS" review-impsamp-mutants.R make fix_h > $L/make-fix_h.log 2>&1
"$RS" -e 'p <- "C:/Users/adf44/source/r/rev-impsamp-src/fix_h/R/importance.R"; s <- readLines(p); i <- grep("^  for \(j in seq_along\(fix\)\) \{$", s); stopifnot(length(i) == 1L); s <- append(s, c("  djoint <- as.numeric(Matrix::diag(env$spHess(full, random = FALSE)))[fix]"), i - 1L); writeLines(s, p)' >> $L/make-fix_h.log 2>&1
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/rev-impsamp-lib/fix_h --no-multiarch --no-docs --no-test-load C:/Users/adf44/source/r/rev-impsamp-src/fix_h >> $L/make-fix_h.log 2>&1
tail -1 $L/make-fix_h.log
for fam in poisson bernoulli; do
  REV_ARM=mutant REV_MUTLIB=C:/Users/adf44/source/r/rev-impsamp-lib/fix_h "$RS" review-impsamp-scale.R $fam 1 > review-impsamp-out/scale-fixh-$fam.log 2>&1 &
done
wait
"$RS" review-impsamp-mutants.R test fix_h test-importance.R > $L/test-fix_h.log 2>&1
echo done
