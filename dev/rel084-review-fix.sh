#!/usr/bin/env bash
# After the release review (dev/reviews/2026-11-18-rel084.md, F1 to F4):
# roxygen of core, checkRd on the Rd files it moved, core reinstalled into
# rellib-r23, the touched test files (one R process each; test-compat.R
# without and with every extension), build-docs and check_pkgdown for core.
# One R process at a time. Logs in dev/rel084-log/fix/.
#   bash dev/rel084-review-fix.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-env.sh
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel084
cd "$ROOT" || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LIB=C:/Users/adf44/source/r/rellib-r23
L=dev/rel084-log/fix
mkdir -p $L
export REL084_LIB=$LIB NOT_CRAN=true
snap() { find man NAMESPACE -type f | sort | xargs md5sum | sed 's/\*//'; }
snap > $L/md5-before.txt
"$RS" dev/rel084-rox.R frmtmb > $L/rox.txt 2>&1; echo "rox exit $?"
snap > $L/md5-after.txt
moved=$(diff $L/md5-before.txt $L/md5-after.txt | awk '/^>/ {print $3}')
echo "Rd files roxygen moved: ${moved:-none}"
for f in $moved man/quantile_elf.Rd; do
  "$RS" -e "r <- tools::checkRd('$f'); cat('checkRd $f:', length(r), 'problems\n'); print(r)"
done 2>&1 | grep -v '^character(0)$'
bash dev/rel084-install.sh $LIB core > $L/install.txt 2>&1; tail -2 $L/install.txt
"$RS" dev/rel084-installed.R > $L/installed.txt 2>&1; grep '^frmtmb ' $L/installed.txt | cut -c1-110
printf 'frmtmb tests/testthat/test-families-quantile.R\nfrmtmb tests/testthat/test-compat.R\nfrmtmb tests/testthat/test-ald-fixed.R\n' > $L/fix.jobs
ROUND_BASE=$LIB bash dev/ciharden-scan.sh r084-fix ref base 1 $L/fix.jobs > $L/scan.txt 2>&1
cat dev/ciharden-log/r084-fix.sum
grep -h '^lib:' dev/ciharden-log/r084-fix-suite/*.txt | awk '{print $4}' | sort | uniq -c
for a in none all; do
  "$RS" dev/rel084-kt/compat-ext.R $a > $L/compat-$a.txt 2>&1
  grep -h '^lib\|^RESULT' $L/compat-$a.txt
done
export R_LIBS="$LIB;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:$PATH"
"$RS" -e 'pkgdown::check_pkgdown("C:/Users/adf44/source/r/frmtmb-wt-rel084"); cat("PKGDOWN frmtmb OK\n")' > $L/pkgdown.txt 2>&1
grep PKGDOWN $L/pkgdown.txt || tail -5 $L/pkgdown.txt
rm -rf dev/docsci-site-rel084
"$RS" dev/release/build-docs.R --dest=C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/docsci-site-rel084 --pkgs=frmtmb --min-sites=1 > $L/build-docs.txt 2>&1
echo "build-docs exit $?"; tail -4 $L/build-docs.txt | cut -c1-140
rm -rf dev/docsci-site-rel084
echo REVIEW-FIX-END
