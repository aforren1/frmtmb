#!/bin/bash
# The review's trial merge (dev/review-predread-merge.sh) on this round's hunks,
# in the session scratchpad: three-way over d5d0056d with the worktrees of latres
# and nlfix as they stand, installed under wt-predread-lib/merge, the lane's own.
# Runs the me() and predict tests on it, for M5 and m7 on the merged tree.
#   bash dev/predread-merge-check.sh <scratch dir>
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
S="$1"
L=C:/Users/adf44/source/r/wt-predread-lib/merge
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
mkdir -p "$S" "$L"
sed "s#^M=.*#M=$S/merge#; s#^LOG=.*#LOG=dev/predread-log/merge.txt#" \
  dev/review-predread-merge.sh > dev/predread-log/merge-copy.sh
bash dev/predread-log/merge-copy.sh > /dev/null 2>&1
grep -c "<<<<<<<" "$S"/merge/R/*.R | grep -v ":0" > dev/predread-log/merge-markers-before.txt
sed "s#^M <- .*#M <- \"$S/merge\"#" dev/review-predread-merge-resolve.R \
  > dev/predread-log/merge-resolve-copy.R
# this round's two more hunks, keeping both sides: latres's mean and this lane's
# in get_predict(), latres's shift and this lane's cov = TRUE frame in idx_at()
{
  echo 'res("R/interop.R", function(ours, theirs) c(theirs, ours))'
  echo 'res("extensions/frmtmb.sample/R/methods-draws.R", function(ours, theirs) c(theirs, ours[-1L]))'
  echo 'left <- system2("grep", c("-l", shQuote("^<<<<<<< "), "-r", shQuote(file.path(M, "R")), shQuote(file.path(M, "extensions/frmtmb.sample/R"))), stdout = TRUE)'
  echo 'cat("markers left after this round:", length(left), "\n")'
} >> dev/predread-log/merge-resolve-copy.R
"$RS" dev/predread-log/merge-resolve-copy.R > dev/predread-log/merge-resolve.txt 2>&1
export R_LIBS="$L;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$RC" CMD INSTALL --library="$L" --no-multiarch --no-docs "$S/merge" \
  > dev/predread-log/merge-install.txt 2>&1
"$RC" CMD INSTALL --library="$L" --no-multiarch --no-docs \
  "$S/merge/extensions/frmtmb.sample" >> dev/predread-log/merge-install.txt 2>&1
tail -1 dev/predread-log/merge-install.txt
printf '%s\n' 'frmtmb tests/testthat/test-me-newdata.R' \
  'frmtmb tests/testthat/test-me-newdata-brms.R' \
  'frmtmb tests/testthat/test-linpred-layers.R' \
  'frmtmb tests/testthat/test-resp-read.R' \
  'frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-readers-resp-me.R' \
  > dev/predread-log/jobs-merge.txt
NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_RELLIB="$L" \
  FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-predread/dev/stan-cache \
  bash dev/release/run-par.sh predread-merge dev/predread-log/jobs-merge.txt 5 \
  > dev/predread-log/run-merge.txt 2>&1
echo merge check done
