#!/bin/sh
# Does a frmtmb.sample test file leave Rplots.pdf behind? Deletes the
# file BEFORE each run, so a leftover cannot read as a hit, runs each
# file in its own process with the gated environment on (the superset:
# a drawing call in a gated block draws only there), and checks.
#   sh dev/simnewdata-rplots.sh <label> [file ...]
# SIMNEWDATA_LIB=base measures the base build.
ROOT=C:/Users/adf44/source/r/frmtmb-wt-simnewdata
cd "$ROOT" || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
DIR=extensions/frmtmb.sample/tests/testthat
PDF=$DIR/Rplots.pdf
label=$1
shift
LOG=dev/simnewdata-log/rplots-$label.txt
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
: > "$LOG"
echo "SIMNEWDATA_LIB=${SIMNEWDATA_LIB:-lane}" >> "$LOG"
for f in "$@"; do
  rm -f "$PDF"
  if [ -e "$PDF" ]; then echo "could not clear $PDF" >> "$LOG"; exit 1; fi
  res=$("$R" dev/simnewdata-runtest.R "$DIR/$f" frmtmb.sample 2>&1 |
        grep "^RESULT")
  [ -n "$res" ] || res="NO RESULT LINE"
  if [ -e "$PDF" ]; then
    echo "DRAWS  $f  $res" >> "$LOG"
  else
    echo "clean  $f  $res" >> "$LOG"
  fi
done
rm -f "$PDF"
cat "$LOG"
