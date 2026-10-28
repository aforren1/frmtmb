#!/bin/sh
# Record every generated bin-2 file on the lane build into
# dev/sigmainn-bin2-log/, for a dry run of the bin-2 ledger
# (FRMTMB_BIN2_REC_DIR), as dev/loglik-bin2-record.sh did. The runner
# stages per package and clears the stage first, so the frmtmb files
# run one after another; the frmtmb.sample file runs beside them.
cd C:/Users/adf44/source/r/frmtmb-wt-sigmainn || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export FRMTMB_PORT_ROOT="$PWD"
export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/wt-sigmainn-lib
export R_LIBS_SITE=C:/Users/adf44/source/r/rellib-r14
L=dev/sigmainn-bin2-log
mkdir -p $L
mode=${1:-record}
one() {
  pkg=$(basename "$(dirname "$1")")
  topic=$(basename "$1" .R | sed 's/^test-brms-suite-bin2-//')
  if [ "$mode" = tier ]; then
    "$R" dev/brmsport-bin2-run.R "$pkg" "$1" > "$L/tier-$pkg-$topic.txt" 2>&1
  else
    "$R" dev/brmsport-bin2-run.R "$pkg" "$1" "$L/rec-$pkg-$topic.tsv" \
      > "$L/run-$pkg-$topic.txt" 2>&1
  fi
}
one dev/brmsport-bin2/frmtmb.sample/test-brms-suite-bin2-methods.R &
for f in dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-*.R; do one "$f"; done
wait
grep -h "^lib:\|^RESULT" $L/run-*.txt $L/tier-*.txt 2>/dev/null
