#!/usr/bin/env bash
# Run the lane's mutants on test-importance.R, 4 at a time, each stopped at 25
# minutes (a mutant with a wrong density climbs the refusal's ladder on every
# fit): bash dev/impsamp-mutrun.sh [mutants]
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
L=dev/impsamp-out/mut
mkdir -p $L
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export R_LIBS="C:/Users/adf44/source/r/wt-impsamp-lib;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
one() {
  m=$1
  "$RS" dev/impsamp-mutants.R make "$m" > $L/make-$m.log 2>&1
  timeout 1500 "$RS" dev/impsamp-mutants.R test "$m" > $L/test-$m.log 2>&1 ||
    echo "ENDED exit $? $m" >> $L/test-$m.log
}
n=0
for m in ${1:-partial_rate abs_step beta_only no_rate no_cross frozen no_jacobian half_jacobian map_recycle no_prior_mode no_refusal one_round grad_always no_cut capped_settles no_prior}; do
  one "$m" &
  n=$((n + 1))
  if [ $((n % 4)) -eq 0 ]; then wait; fi
done
wait
grep -h "RESULT\|FAILED\|ENDED" $L/test-*.log
echo mutants-done
