#!/bin/sh
# The whole bin-2 port of brms's suite against one build, in order:
# survey, generate, record, verdicts, ledger, regenerate, and the tier
# run that asserts the verdicts. The consolidation runs it on the
# merged build after dev/brmsport-ledger.R has regenerated bin 1, so the
# composed totals read the regenerated bin-1 ledger.
#
#   FRMTMB_PORT_LIB=<lib> sh dev/brmsport-bin2-all.sh
#
# FRMTMB_PORT_LIB defaults to rellib-r9 (dev/brmsport-bin2-run.R). The
# ledger stops on a stale or missing verdict, so a run on another build
# names every row that moved.
set -e
cd "$(dirname "$0")/.."
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export FRMTMB_PORT_ROOT="$PWD"
"$R" dev/brmsport-bin2-survey.R > dev/brmsport-bin2-log/survey.txt 2>&1 ||
  { cat dev/brmsport-bin2-log/survey.txt; exit 1; }
"$R" dev/brmsport-bin2-gen.R
sh dev/brmsport-bin2-record.sh
"$R" dev/brmsport-bin2-verdicts-2t.R
if [ -f dev/brmsport-bin2-verdicts-2p.R ]; then
  "$R" dev/brmsport-bin2-verdicts-2p.R
fi
"$R" dev/brmsport-bin2-ledger.R > dev/brmsport-bin2-log/ledger.txt 2>&1 ||
  { cat dev/brmsport-bin2-log/ledger.txt; exit 1; }
"$R" dev/brmsport-bin2-gen.R
sh dev/brmsport-bin2-record.sh tier
tail -8 dev/brmsport-bin2-summary.md
