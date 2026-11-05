#!/usr/bin/env bash
# Short check, item 7: install the refreshed trial merge and run core's and frmtmb.sample's
# whole suites on it, one file per process, gated, against the union of the three lanes'
# Stan caches, first the two files that never had a reading.
set -u
cd /c/Users/adf44/source/r/rev-predread-src/merge
P=/c/Users/adf44/source/r/frmtmb-wt-predread
O=$P/dev/review-predread-out/r2/merge-tests
mkdir -p "$O"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
LIB=C:/Users/adf44/source/r/rev-predread-lib/merge
rm -rf "$LIB"; mkdir -p "$LIB"
cp -rn "$P/extensions/." extensions/
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
: > "$O/../merge-install.txt"
"$RC" CMD INSTALL --library="$LIB" --no-multiarch --no-docs --no-test-load . >> "$O/../merge-install.txt" 2>&1
for e in sample eam learn survival; do
  "$RC" CMD INSTALL --library="$LIB" --no-multiarch --no-docs --no-test-load "extensions/frmtmb.$e" >> "$O/../merge-install.txt" 2>&1
done
grep -a "DONE\|ERROR" "$O/../merge-install.txt" > "$O/../merge-install.sum"
unset R_LIBS
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/rev-predread-src/stan-cache
N=${1:-6}
: > "$O/../merge-jobs.txt"
echo "frmtmb tests/testthat/test-latres-brms.R" >> "$O/../merge-jobs.txt"
echo "frmtmb tests/testthat/test-nl-brms.R" >> "$O/../merge-jobs.txt"
for f in tests/testthat/test-*.R; do
  case "$f" in */test-latres-brms.R|*/test-nl-brms.R) ;; *) echo "frmtmb $f" >> "$O/../merge-jobs.txt" ;; esac
done
for f in extensions/frmtmb.sample/tests/testthat/test-*.R; do echo "frmtmb.sample $f" >> "$O/../merge-jobs.txt"; done
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  ( "$RS" "$P/dev/review-predread-runone.R" "$p" "$f" "$LIB" > "$O/$p--$(basename "$f" .R).txt" 2>&1 ) &
done < "$O/../merge-jobs.txt"
wait
echo finished > "$O/../merge-tests.done"
