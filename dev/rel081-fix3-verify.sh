#!/usr/bin/env bash
# The reruns after dev/rel081-release-fix3.R, which changed frmtmb.sample's
# code and tests (and NEWS text only elsewhere), on rellib-r20 with
# frmtmb.sample reinstalled (dev/rel081-install-sample.sh). One step at a
# time, at most three R processes:
#   1. roxygen of frmtmb.sample again (diff must be empty), installed code
#      against the tree (all nine), check_pkgdown() (all nine);
#   2. every frmtmb.sample test file, one per process, plain then gated;
#   3. the frmtmb.sample sessions, three orders, plain then gated;
#   4. R CMD check --as-cran of frmtmb.sample.
# Progress lines go to dev/rel081-log/chain.txt as FIX3.
#   bash dev/rel081-fix3-verify.sh > dev/rel081-log/fix3-verify.txt 2>&1
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r20
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/rel081-log
note() { echo "$(date '+%F %T') FIX3 $*" >> $L/chain.txt; }
note "verify start"
## 1
snap() { find extensions/frmtmb.sample/man extensions/frmtmb.sample/NAMESPACE -type f | sort | xargs md5sum; }
snap > $L/fix3-doc-md5-1.txt
REL081_LIB=$ROUND_BASE "$RS" dev/rel081-rox.R frmtmb.sample > $L/fix3-rox.txt 2>&1
snap > $L/fix3-doc-md5-2.txt
echo "sample files whose roxygen output moved: $(diff $L/fix3-doc-md5-1.txt $L/fix3-doc-md5-2.txt | grep -c '^>')"
"$RS" dev/rel081-installed.R > $L/fix3-installed.txt 2>&1
cut -c1-110 $L/fix3-installed.txt
# check_pkgdown() reads the docs through pandoc (the first run of this
# script left it off PATH and read "Pandoc not available" nine times)
( export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
  export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:$PATH"
  "$RS" dev/rel081-pkgdown.R > $L/fix3-pkgdown.txt 2>&1 )
grep PKGDOWN $L/fix3-pkgdown.txt
note "docs end"
## 2
: > $L/fix3-sample.jobs
for f in extensions/frmtmb.sample/tests/testthat/test-*.R; do
  echo "frmtmb.sample $f" >> $L/fix3-sample.jobs
done
echo "sample test files: $(grep -c '' $L/fix3-sample.jobs)"
CIHARDEN_UNGATED=true bash dev/ciharden-scan.sh r081-fix3-plain ref base 3 $L/fix3-sample.jobs
note "plain end"
bash dev/ciharden-scan.sh r081-fix3-ref ref base 3 $L/fix3-sample.jobs
note "gated end"
## 3
NOT_CRAN=true bash dev/release/run-sessions.sh rel0.81.0s 3 frmtmb.sample > $L/fix3-sessions-plain.txt 2>&1
bash dev/rel081-sessions-sum.sh rel0.81.0s > $L/fix3-sessions-plain-sum.txt 2>&1
NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
  FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.81.0sg 3 frmtmb.sample > $L/fix3-sessions-gated.txt 2>&1
bash dev/rel081-sessions-sum.sh rel0.81.0sg > $L/fix3-sessions-gated-sum.txt 2>&1
note "sessions end"
## 4
bash dev/rel081-check.sh frmtmb.sample > $L/fix3-check.txt 2>&1
cat $L/fix3-check.txt
mkdir -p $L/checklogs
cp dev/rel081-check/frmtmb.sample/check.log $L/checklogs/check-frmtmb.sample-fix3.log
cp dev/rel081-check/frmtmb.sample/frmtmb.sample.Rcheck/00check.log $L/checklogs/00check-frmtmb.sample-fix3.log
rm -rf dev/rel081-check
note "check end"
note "verify END"
