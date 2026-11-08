#!/bin/sh
# Review of lane mfxse, short check: the oneblock and sandown mutants rebuilt
# on the fix-round source (R/interop.R of 13:49), into
# rev-mfxse-lib/sc-<mutant>, then test-marginaleffects-re.R on each, one at
# a time.
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
W=/c/Users/adf44/source/r/frmtmb-wt-mfxse
SRC=/c/Users/adf44/source/r/rev-mfxse-src
LIB=/c/Users/adf44/source/r/rev-mfxse-lib
OUT=$W/dev/review-mfxse-out/sc/mut
mkdir -p "$OUT"
export R_LIBS=/c/Users/adf44/source/r/rellib-r18\;/c/Users/adf44/AppData/Local/R/win-library/4.6
mk() {
  name=sc-$1; sedexpr=$2
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
  "$RS" $W/dev/review-mfxse-run-test.R $W/tests/testthat/test-marginaleffects-re.R \
    "$LIB/$name" > "$OUT/$name-test.log" 2>&1
  echo "$name: $(grep -h -E '^lib:' $OUT/$name-test.log) $(grep -h RESULT $OUT/$name-test.log)"
  grep -E "^── (Failure|Error)" "$OUT/$name-test.log"
}
mk oneblock '749s/b_idx <- c(b_idx, bk\[\["b_idx"\]\])/if (!length(b_idx)) b_idx <- c(b_idx, bk[["b_idx"]])/'
mk sandown '1307s/interop_vcov_re(model, V, dots, Vf = Vf)/interop_vcov_re(model, V, dots)/'
