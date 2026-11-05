#!/usr/bin/env bash
# The targeted rerun after the short check's edits, on the install of 18:36:20 of core.
# Each step writes a line to dev/rel079-log/short.txt, and the script ends by printing
# "SHORT END". Do not stop this shell.
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/rel079-env.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel079
cd "$ROOT" || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/rel079-log
S=$L/short.txt
LIB=C:/Users/adf44/source/r/rellib-r17
export ROUND_BASE=$LIB
note() { echo "$(date '+%F %T') $*" >> "$S"; }
: > "$S"
note "install stamp $(stat -c %y ../rellib-r17/frmtmb/DESCRIPTION | cut -c1-19)"
cat > $L/short.jobs <<'EOF'
frmtmb tests/testthat/test-latres-newrow.R
frmtmb tests/testthat/test-latres-newrow-brms.R
frmtmb tests/testthat/test-latres.R
frmtmb tests/testthat/test-latres-brms.R
frmtmb tests/testthat/test-me-newdata.R
frmtmb tests/testthat/test-marginaleffects-args.R
frmtmb tests/testthat/test-conditions-census.R
frmtmb tests/testthat/test-message-uniqueness.R
frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-latres.R
EOF
bash dev/ciharden-scan.sh r079-short ref base 8 $L/short.jobs > /dev/null 2>&1
note "gated end"
printf 'frmtmb tests/testthat/test-latres-newrow.R\n' > $L/short-ob.jobs
bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/openblas-build2.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack >> $L/openblas-build2.txt 2>&1
bash dev/ciharden-scan.sh r079-short-ob26 ob0.3.26-lapack-t4 base 1 $L/short-ob.jobs > /dev/null 2>&1
bash dev/ciharden-scan.sh r079-short-ob32 ob0.3.32-lapack-t4 base 1 $L/short-ob.jobs > /dev/null 2>&1
note "openblas end"
for k in "3000 block" "3001 block" "20001 block"; do
  set -- $k
  NOT_CRAN=true "$RS" dev/rel079-probe-limits.R $1 $2 > $L/limits-short-$1-$2.txt 2>&1
done
note "limits end"
export R_LIBS="$LIB;C:/Users/adf44/AppData/Local/R/win-library/4.6" REL079_LIB="$LIB" NOT_CRAN=true
snap() { find man NAMESPACE extensions/*/man extensions/*/NAMESPACE -type f | sort | xargs md5sum | sed 's/\*//'; }
snap > $L/doc-md5-sc1.txt
"$RS" dev/rel079-rox.R frmtmb frmtmb.eam frmtmb.latent frmtmb.ode frmtmb.survival frmtmb.spatial frmtmb.coupling frmtmb.sample frmtmb.learn > $L/rox-sc2.txt 2>&1
snap > $L/doc-md5-sc2.txt
note "roxygen again moved: $(diff $L/doc-md5-sc1.txt $L/doc-md5-sc2.txt | grep -c '^>')"
note "against the last chain's roxygen: $(diff $L/doc-md5-2.txt $L/doc-md5-sc2.txt | grep '^>' | awk '{print $3}' | tr '\n' ' ')"
"$RS" dev/rel079-installed.R > $L/installed-sc.txt 2>&1
note "installed end"
bash dev/rel079-check.sh frmtmb > $L/check-sc.txt 2>&1
note "check end: $(grep -h 'Status:' $L/check-sc.txt | tr '\n' ' ')"
rm -rf dev/ciharden-out
note "openblas copies deleted"
note "SHORT END"
echo "SHORT END"
