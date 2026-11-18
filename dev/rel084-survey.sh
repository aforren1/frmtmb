#!/usr/bin/env bash
# The recorded survey of the 0.84.0 release. First the compat frames of
# rellib-r22 (0.83.0) and of this build, and their diff (frm_compat() rows
# that moved this round), then every pair with RECORD=1, one pair at a
# time: each pair runs a live-row child, so it holds 2 R processes, which
# leaves one for the work beside it. Adapted from dev/rel081-survey.sh.
#
#   bash dev/rel084-survey.sh [frames|survey|all]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel084 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
G=dev/rel084-log
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
case "${1:-all}" in
  frames|all)
    "$RS" dev/equiv-dump-compat.R "$L/rellib-r22" $G/frame-r22 > $G/frame-r22.txt 2>&1
    "$RS" dev/equiv-dump-compat.R "$L/rellib-r23" $G/frame-r23 > $G/frame-r23.txt 2>&1
    "$RS" dev/compat-frame-diff.R $G/frame-r22/compat.tsv $G/frame-r23/compat.tsv \
      > $G/frame-diff-r22-r23.txt 2>&1
    echo "FRAMES END";;
esac
case "${1:-all}" in
  survey|all)
    rm -rf dev/rel084-equiv-out
    FRMTMB_EQUIV_LIB="$L/rellib-r23;$L/rellib-r22" RECORD=1 bash dev/equiv-run.sh 1 dev/rel084-equiv-out \
      > $G/survey.txt 2>&1
    echo "SURVEY END $?" >> $G/survey.txt;;
esac
