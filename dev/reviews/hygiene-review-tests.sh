#!/usr/bin/env bash
# Review of lane hygiene: rerun the lane's changed and new test files on
# the final tree, three R processes at most, one file per process.
#   bash dev/reviews/hygiene-review-tests.sh > dev/reviews/hygiene-review-tests.log
cd /c/Users/adf44/source/r/frmtmb-wt-hygiene || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
T=tests/testthat
E=extensions
# arm bare package file
jobs="base true frmtmb $T/test-compat-register.R
base true frmtmb $T/test-structure.R
base true frmtmb.sample $E/frmtmb.sample/tests/testthat/test-loo.R
base true frmtmb.ode $E/frmtmb.ode/tests/testthat/test-ode.R
base false frmtmb $T/test-register-once.R
base false frmtmb $T/test-registry-helper.R
mine false frmtmb $T/test-register-once.R
mine false frmtmb $T/test-registry-helper.R
mine false frmtmb $T/test-weights-constant.R
mine false frmtmb $T/test-codemeta.R
mine false frmtmb $T/test-port-harness.R
mine false frmtmb $T/test-guardrail1.R
mine false frmtmb $T/test-compat-register.R
mine false frmtmb $T/test-structure.R
mine false frmtmb $T/test-prior-table-brms.R
mine false frmtmb.sample $E/frmtmb.sample/tests/testthat/test-loo.R
mine false frmtmb.ode $E/frmtmb.ode/tests/testthat/test-ode.R
mine false frmtmb.sample $E/frmtmb.sample/tests/testthat/test-brms-suite-helper-copy.R"
i=0
echo "$jobs" | while read -r arm bare p f; do
  i=$((i + 1))
  echo "$i $arm $bare $p $f"
done | xargs -P 3 -L 1 bash -c '
  port=$((11950 + $0))
  out=$(HYG_ARM=$1 HYG_BARE=$2 R_PARALLEL_PORT=$port "'"$RS"'" \
        dev/hygiene-run-tests.R $3 $4 2>&1)
  echo "$out" | grep -E "^(RESULT|DETAIL|lib:)" | sed "s|^|[$1 bare=$2] |"
'
