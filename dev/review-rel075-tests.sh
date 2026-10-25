#!/usr/bin/env bash
# Review: the brief's test files, one file per R process, rellib-r13.
REL=/c/Users/adf44/source/r/frmtmb-wt-rel075
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$REL/dev/review-rel075-log/tests
mkdir -p "$OUT"
export NOT_CRAN=true
cat > "$OUT/jobs" <<J
frmtmb tests/testthat/test-guardrail1.R
frmtmb tests/testthat/test-extension-home-mixed.R
frmtmb tests/testthat/test-brms-names-specials.R
frmtmb tests/testthat/test-autocor-seam.R
frmtmb tests/testthat/test-fcor.R
frmtmb tests/testthat/test-aterm-option.R
frmtmb tests/testthat/test-simplex-pars.R
frmtmb tests/testthat/test-compat.R
frmtmb tests/testthat/test-message-uniqueness.R
frmtmb.survival extensions/frmtmb.survival/tests/testthat/test-cox-bhaz.R
frmtmb.spatial extensions/frmtmb.spatial/tests/testthat/test-sar.R
frmtmb.spatial extensions/frmtmb.spatial/tests/testthat/test-guardrail.R
frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-autocor-seam-draws.R
frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-draws-names-specials.R
J
run1() {
  p=$1; f=$2
  log="$OUT/$p--$(basename "$f" .R).log"
  "$RS" "$REL/dev/ciharden-run1.R" base "$p" "$REL/$f" > "$log" 2>&1
  echo "$(grep '^RESULT' "$log" || echo "NO RESULT $f") | $(grep '^lib:' "$log")"
}
export -f run1; export OUT REL RS
xargs -P 5 -L 1 bash -c 'run1 "$0" "$1"' < "$OUT/jobs" > "$OUT/summary.txt"
echo ALLDONE >> "$OUT/summary.txt"
