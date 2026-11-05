#!/usr/bin/env bash
# After the chain: the last text fix on core (a link in ?frmtmb-autocor) and one test
# bound, so core is roxygenised and installed again, its check rerun, and the one
# test file run on the reference BLAS and both OpenBLAS copies. No function body
# changed after the chain: dev/rel079-installed.R counts that.
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/rel079-env.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel079
cd "$ROOT" || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
L=dev/rel079-log
LIB=C:/Users/adf44/source/r/rellib-r17
export R_LIBS="$LIB;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win REL079_LIB="$LIB"
"$RS" dev/rel079-rox.R frmtmb > $L/rox-final.txt 2>&1
"$RC" CMD INSTALL --library="$LIB" --no-multiarch "$ROOT" > $L/install-final.txt 2>&1
tail -1 $L/install-final.txt
printf 'frmtmb tests/testthat/test-latres-brms.R\nfrmtmb tests/testthat/test-latres.R\n' > $L/final.jobs
export ROUND_BASE=$LIB
bash dev/ciharden-scan.sh r079-final ref base 2 $L/final.jobs > /dev/null 2>&1
bash dev/ciharden-scan.sh r079-final-ob26 ob0.3.26-lapack-t4 base 2 $L/final.jobs > /dev/null 2>&1
bash dev/ciharden-scan.sh r079-final-ob32 ob0.3.32-lapack-t4 base 2 $L/final.jobs > /dev/null 2>&1
bash dev/rel079-check.sh frmtmb > $L/check-2.txt 2>&1
for r in r079-final r079-final-ob26 r079-final-ob32; do cut -c1-90 dev/ciharden-log/$r.sum; done
grep -E "Status:" $L/check-2.txt
bash dev/rel079-rerun.sh latres review-latres-tk.R fresh f_geom 7001 7030 $ROOT/$L/rr/latres/tk/fresh-f_geom.tsv
bash dev/rel079-rerun.sh latres review-latres-tk.R m4b blog_w 2001 2040 $ROOT/$L/rr/latres/tk/m4b-blog_w.tsv
bash dev/rel079-rerun.sh nlfix review-nlfix-d7.R cratio
echo FINAL-END
