#!/usr/bin/env bash
# The docs at the end: roxygen once more over all nine packages, whose diff must be
# empty; installed code against the tree; check_pkgdown() for the reference indexes;
# the nine sites built into a gitignored folder.
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/rel079-env.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel079
cd "$ROOT" || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/rel079-log
LIB=C:/Users/adf44/source/r/rellib-r17
export R_LIBS="$LIB;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:$PATH"
export REL079_LIB="$LIB" NOT_CRAN=true
snap() { find man NAMESPACE extensions/*/man extensions/*/NAMESPACE -type f | sort | xargs md5sum | sed 's/\*//'; }
snap > $L/doc-md5-1.txt
"$RS" dev/rel079-rox.R frmtmb frmtmb.eam frmtmb.latent frmtmb.ode frmtmb.survival frmtmb.spatial frmtmb.coupling frmtmb.sample frmtmb.learn > $L/rox-again.txt 2>&1
snap > $L/doc-md5-2.txt
echo "files whose roxygen output moved: $(diff $L/doc-md5-1.txt $L/doc-md5-2.txt | grep -c '^>')"
diff $L/doc-md5-1.txt $L/doc-md5-2.txt | head -6
"$RS" dev/rel079-installed.R > $L/installed.txt 2>&1
cut -c1-110 $L/installed.txt
"$RS" dev/rel079-pkgdown.R > $L/pkgdown.txt 2>&1
grep PKGDOWN $L/pkgdown.txt | cut -c1-140
rm -rf dev/docsci-site-rel079
"$RS" dev/release/build-docs.R --dest=C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/docsci-site-rel079 > $L/build-docs.txt 2>&1
echo "build-docs exit $?"
tail -5 $L/build-docs.txt | cut -c1-140
echo DOCS-END
