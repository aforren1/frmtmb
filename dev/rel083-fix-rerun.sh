#!/usr/bin/env bash
# After the release review's F2 and F3 and the test-sugar.R seed: roxygen
# of core, a core reinstall into rellib-r22, a second roxygen run whose
# moves must be 0, check_pkgdown() and checkRd of core, installed code
# against the tree, the regenerated guards, and the changed test files
# (gated, one per process, 3 at once).
#   bash dev/rel083-fix-rerun.sh > dev/rel083-log/fix-rerun.txt 2>&1
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel083 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export REL083_LIB=C:/Users/adf44/source/r/rellib-r22 ROUND_BASE=C:/Users/adf44/source/r/rellib-r22
export RSTUDIO_PANDOC="C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:$PATH"
snap() { find man NAMESPACE -type f | sort | xargs md5sum | sed 's/\*//'; }
snap > dev/rel083-log/fix-md5-0.txt
"$RS" dev/rel083-rox.R frmtmb > dev/rel083-log/fix-rox-1.txt 2>&1
snap > dev/rel083-log/fix-md5-1.txt
echo "roxygen 1 moved: $(diff dev/rel083-log/fix-md5-0.txt dev/rel083-log/fix-md5-1.txt | grep '^>' | awk '{print $2}' | tr '\n' ' ')"
bash dev/rel083-install.sh C:/Users/adf44/source/r/rellib-r22 core
"$RS" dev/rel083-rox.R frmtmb > dev/rel083-log/fix-rox-2.txt 2>&1
snap > dev/rel083-log/fix-md5-2.txt
echo "roxygen 2 moved: $(diff dev/rel083-log/fix-md5-1.txt dev/rel083-log/fix-md5-2.txt | grep -c '^>')"
"$RS" -e '.libPaths(c("C:/Users/adf44/source/r/rellib-r22","C:/Users/adf44/AppData/Local/R/win-library/4.6")); cat("checkRd problems:", length(tools::checkRd("man/frmtmb_control.Rd")), "\n"); r <- tryCatch({pkgdown::check_pkgdown("."); "OK"}, error = function(e) conditionMessage(e)); cat("PKGDOWN frmtmb", r, "\n")' 2>&1 | grep -E 'checkRd|PKGDOWN'
"$RS" dev/rel083-installed.R 2>&1 | cut -c1-110
FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r22 FRMTMB_PORT_EXTENSIONS=all "$RS" dev/brmsport-guards.R > dev/rel083-log/guards-2.txt 2>&1
echo "guards exit $?: $(tail -1 dev/rel083-log/guards-2.txt | cut -c1-120)"
printf '%s\n' "frmtmb tests/testthat/test-nl-starts-control.R" "frmtmb tests/testthat/test-sugar.R" \
  "frmtmb tests/testthat/test-brms-suite-methods.R" > dev/rel083-log/fix.jobs
bash dev/ciharden-scan.sh r083-fix ref base 3 dev/rel083-log/fix.jobs
cat dev/ciharden-log/r083-fix.sum
grep -h '^lib' dev/ciharden-log/r083-fix-suite/*.txt | awk '{print $4, $7}' | sort | uniq -c
echo FIX RERUN END
