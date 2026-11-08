#!/usr/bin/env bash
# Release review 0.80.0: frmtmb.ode 0.7.3 and frmtmb.survival 0.2.2 (copies
# of rellib-r19's builds) over rellib-r18, every test file, one per process.
#   bash dev/review-rel-floor.sh
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/review-rel-out/floor; mkdir -p $O
for p in ode survival; do
  for f in extensions/frmtmb.$p/tests/testthat/test-*.R; do
    o=$O/$p-$(basename $f .R).log
    "$R" dev/review-rel-floor-run.R C:/Users/adf44/source/r/rev-rel-lib/floor frmtmb.$p $f > $o 2>&1
    echo "$(grep -h '^RESULT' $o | tail -1) | $(grep -h '^lib:' $o | head -1)"
  done
done
echo "FLOOR END"
