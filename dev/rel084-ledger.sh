#!/usr/bin/env bash
# The ported brms suite, regenerated once on the final build, as at 0.78.0
# (dev/round-20261029.md, "The ported brms suite"): bin 1 in record mode, the ledger,
# the generated test files, the guards, then bin 2. Each step's log is
# dev/rel084-log/<step>.txt. The two ledgers before are kept for the diffs.
#
#   bash dev/rel084-ledger.sh [bin1|ledger|bin2]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel084 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/rel084-log
export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r23 FRMTMB_PORT_EXTENSIONS=all
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
[ -f $L/ledger-before.tsv ] || git show HEAD:dev/brmsport-ledger.tsv | tr -d '\r' > $L/ledger-before.tsv
[ -f $L/bin2-ledger-before.tsv ] || git show HEAD:dev/brmsport-bin2-ledger.tsv | tr -d '\r' > $L/bin2-ledger-before.tsv
case "${1:-all}" in
  bin1|all)
    bash dev/rel084-record.sh brmsformula brmsterms families helper-copy methods brm brmsfit-helpers > $L/record-a.txt 2>&1
    bash dev/rel084-record.sh data-helpers emmeans priors standata > $L/record-b.txt 2>&1
    echo "bin1 done: $(cat $L/record-a.txt $L/record-b.txt | grep -c '^RESULT')";;
esac
case "${1:-all}" in
  ledger|all)
    "$RS" dev/brmsport-ledger.R > $L/ledger-1.txt 2>&1; echo "ledger exit $?"
    tail -3 $L/ledger-1.txt | cut -c1-200;;
esac
case "${1:-all}" in
  gen|all)
    "$RS" dev/brmsport-gen.R > $L/gen.txt 2>&1; echo "gen exit $?"
    "$RS" dev/brmsport-guards.R > $L/guards.txt 2>&1; echo "guards exit $?"
    tail -1 $L/guards.txt | cut -c1-200
    # bin 2 itself: dev/rel084-bin2.sh, which copies the bin-2 files first;;
esac
echo "LEDGER STEP END"
