#!/usr/bin/env bash
# Short check of lane hygiene's fix round: the six changed or new test
# files on the lane library, one file per process, two at a time (the
# latent session is the third process), then the review's regonce
# reproducer on both arms.
#   bash dev/reviews/hygiene-short-tests.sh > dev/reviews/hygiene-short-tests.log
cd /c/Users/adf44/source/r/frmtmb-wt-hygiene || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
T=tests/testthat
jobs="mine false frmtmb $T/test-register-once.R
mine false frmtmb $T/test-registry-helper.R
mine false frmtmb $T/test-structure.R
mine false frmtmb $T/test-compat-register.R
mine false frmtmb $T/test-weights-constant.R
mine false frmtmb $T/test-codemeta.R"
i=0
echo "$jobs" | while read -r arm bare p f; do
  i=$((i + 1))
  echo "$i $arm $bare $p $f"
done | xargs -P 2 -L 1 bash -c '
  port=$((12430 + $0))
  out=$(HYG_ARM=$1 HYG_BARE=$2 R_PARALLEL_PORT=$port "'"$RS"'" \
        dev/hygiene-run-tests.R $3 $4 2>&1)
  echo "$out" | grep -E "^(RESULT|DETAIL|lib:)" | sed "s|^|[$1] |"
'
for arm in base mine; do
  "$RS" dev/reviews/hygiene-review-regonce.R $arm \
    > dev/reviews/hygiene-short-regonce-$arm.log 2>&1
  echo "regonce $arm exit $?"
done
