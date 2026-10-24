#!/usr/bin/env bash
# Release review of 0.74.0: one-session forward runs with the covstruct
# registry added to the LEAK snapshot (dev/review-rel074-run-session.R).
cd /c/Users/adf44/source/r/frmtmb-wt-release
export NOT_CRAN=true R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win FRMTMB_STAN_CACHE=$PWD/dev/stan-cache
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/review-rel074-log/cssess; mkdir -p $L
T=$PWD
"$RS" dev/review-rel074-run-session.R frmtmb forward $T/tests/testthat > $L/frmtmb.txt 2>&1 &
for p in frmtmb.spatial frmtmb.survival frmtmb.latent frmtmb.sample; do
  "$RS" dev/review-rel074-run-session.R $p forward $T/extensions/$p/tests/testthat > $L/$p.txt 2>&1 &
done
wait
grep -h "^SESSION" $L/*.txt
grep -h "^LEAK" $L/*.txt | grep -v "with load" || echo "no untagged LEAK"
grep -h "^LEAK.*covstruct" $L/*.txt || echo "no covstruct LEAK"
