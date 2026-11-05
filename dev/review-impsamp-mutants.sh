#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
L=review-impsamp-out/mutants
mkdir -p $L
rm -f $L/*.log
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
cat > $L/one.sh <<'EOS'
m=$1
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev
L=review-impsamp-out/mutants
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" review-impsamp-mutants.R make "$m" > $L/make-$m.log 2>&1
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/rev-impsamp-lib/$m --no-multiarch --no-docs --no-test-load C:/Users/adf44/source/r/rev-impsamp-src/$m >> $L/make-$m.log 2>&1
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" review-impsamp-mutants.R test "$m" test-importance.R > $L/test-$m.log 2>&1
EOS
for m in ${MUTS:-partial_rate no_rate no_cross frozen no_jacobian half_jacobian map_recycle no_prior_mode no_refusal one_round grad_always no_cut}; do echo $m; done | xargs -P 4 -L 1 bash $L/one.sh
echo done
