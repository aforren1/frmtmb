#!/usr/bin/env bash
# Review of lane hygiene, item 3042: the ported bin-1 brms suite and the
# helper-copy test, gated, on the final tree (the lane's plain tier ran
# them with 0 expectations: skip_unless_brms_suite() needs the gate).
#   bash dev/reviews/hygiene-review-suite.sh <mine|base> > dev/reviews/hygiene-review-suite-<arm>.log
arm=${1:-mine}
cd /c/Users/adf44/source/r/frmtmb-wt-hygiene || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=/c/Users/adf44/source/r/frmtmb-wt-hygiene/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
{
  for f in tests/testthat/test-brms-suite-*.R; do echo "frmtmb $f"; done
  for f in extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R; do
    echo "frmtmb.sample $f"
  done
} | awk '{print NR, $0}' | xargs -P 3 -L 1 bash -c '
  port=$((12010 + $0))
  t0=$(date +%s)
  out=$(HYG_ARM='"$arm"' R_PARALLEL_PORT=$port "'"$RS"'" \
        dev/hygiene-run-tests.R $1 $2 2>&1)
  echo "$out" | grep -E "^(RESULT|DETAIL|lib:)" | sed "s|^|[$1] |"
  echo "[$1] TIME $(basename $2) $(( $(date +%s) - t0 )) s"
'
