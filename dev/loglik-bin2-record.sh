#!/bin/sh
# Lane loglik, fix round 1: record every generated bin-2 file on the
# lane's build into dev/loglik-bin2-log/, for a dry run of the bin-2
# ledger (FRMTMB_BIN2_REC_DIR). The runner stages per package and clears
# the stage first, so the frmtmb files run one after another; the one
# frmtmb.sample file runs beside them.
cd C:/Users/adf44/source/r/frmtmb-wt-loglik || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export FRMTMB_PORT_ROOT="$PWD"
export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/wt-loglik-lib
export R_LIBS_SITE=C:/Users/adf44/source/r/rellib-r13
L=dev/loglik-bin2-log
mkdir -p $L
one() {
  pkg=$(basename "$(dirname "$1")")
  topic=$(basename "$1" .R | sed 's/^test-brms-suite-bin2-//')
  "$R" dev/brmsport-bin2-run.R "$pkg" "$1" "$L/rec-$pkg-$topic.tsv" \
    > "$L/run-$pkg-$topic.txt" 2>&1
}
one dev/brmsport-bin2/frmtmb.sample/test-brms-suite-bin2-methods.R &
for f in dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-*.R; do one "$f"; done
wait
grep -h "^lib:\|^RESULT" $L/run-*.txt
