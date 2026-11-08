#!/bin/sh
# Review of lane mfxse: build mutants of R/interop.R into
# rev-mfxse-lib/<mutant> (source copies in rev-mfxse-src/<mutant>), then run
# test-marginaleffects-re.R on each, one process per file, 2 at once.
#   sh dev/review-mfxse-mutants.sh
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
W=/c/Users/adf44/source/r/frmtmb-wt-mfxse
SRC=/c/Users/adf44/source/r/rev-mfxse-src
LIB=/c/Users/adf44/source/r/rev-mfxse-lib
OUT=$W/dev/review-mfxse-out/mut
mkdir -p "$SRC" "$LIB" "$OUT"
export R_LIBS=/c/Users/adf44/source/r/rellib-r18\;/c/Users/adf44/AppData/Local/R/win-library/4.6
mk() {
  name=$1; sedexpr=$2
  d=$SRC/$name
  rm -rf "$d"; mkdir -p "$d"
  cp -r $W/DESCRIPTION $W/NAMESPACE $W/R $W/inst $W/man $W/data "$d"/
  sed -i "$sedexpr" "$d/R/interop.R"
  if cmp -s "$d/R/interop.R" "$W/R/interop.R"; then
    echo "MUTANT $name: sed changed nothing"; return 1
  fi
  diff "$W/R/interop.R" "$d/R/interop.R" > "$OUT/$name.diff"
  mkdir -p "$LIB/$name"
  "$R" CMD INSTALL --no-multiarch --no-docs --no-html --library="$LIB/$name" "$d" \
    > "$OUT/$name-install.log" 2>&1 || echo "INSTALL FAILED $name"
}
mk oneblock '723s/b_idx <- c(b_idx, bk\[\["b_idx"\]\])/if (!length(b_idx)) b_idx <- c(b_idx, bk[["b_idx"]])/'
mk notheta '752s/return(TRUE)/return(FALSE)/'
mk cross0 '877s/C <- M %\*%/C <- 0 * M %*%/'
mk nohyp '903s/if ("hypothesis_is_formula" %in% names(dots)) {/if (FALSE) {/'
mk sandown '1234s/interop_vcov_re(model, V, dots, Vf = Vf)/interop_vcov_re(model, V, dots)/'
run() {
  "$RS" $W/dev/review-mfxse-run-test.R $W/tests/testthat/test-marginaleffects-re.R \
    "$LIB/$1" > "$OUT/$1-test.log" 2>&1
}
run oneblock & run notheta & wait
run cross0 & run nohyp & wait
run sandown
for m in oneblock notheta cross0 nohyp sandown; do
  echo "$m: $(grep -h -E '^lib:' $OUT/$m-test.log) $(grep -h RESULT $OUT/$m-test.log)"
done
