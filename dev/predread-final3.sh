#!/bin/bash
# Fix round 1's final pass, in order: install, the new tests on the base build,
# every suite one file per process with the gates on, OpenBLAS 0.3.26 and 0.3.32 on
# the numerical files, the audits, R CMD check of the two changed packages, the
# one-session orders of all nine suites, and all 172 equivalence pairs.
#   bash dev/predread-final3.sh
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
L=C:/Users/adf44/source/r/wt-predread-lib
B=C:/Users/adf44/source/r/rellib-r16
U=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-predread/dev/stan-cache
log=dev/predread-log/final3.txt
: > "$log"
step() { echo "$(date +%H:%M:%S) $*" >> "$log"; }

step install
bash dev/predread-install.sh core frmtmb.sample frmtmb.coupling frmtmb.eam \
  frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.spatial frmtmb.survival \
  > dev/predread-log/install-final3.txt 2>&1

step seen
NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true bash dev/release/run-par.sh \
  predread-seen dev/predread-log/jobs-seen.txt 4 > /dev/null 2>&1
NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true bash dev/release/run-par.sh \
  predread-seen2 dev/predread-log/jobs-seen2.txt 2 > /dev/null 2>&1

step suite
ls tests/testthat/test-*.R | sed 's/^/frmtmb /' > dev/predread-log/jobs-all.txt
for p in frmtmb.sample frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn \
         frmtmb.ode frmtmb.spatial frmtmb.survival; do
  ls extensions/$p/tests/testthat/test-*.R | sed "s/^/$p /" >> dev/predread-log/jobs-all.txt
done
NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true \
  FRMTMB_RELLIB="$L" bash dev/release/run-par.sh predread-all \
  dev/predread-log/jobs-all.txt 8 > dev/predread-log/run-all.txt 2>&1
"$RS" dev/predread-suite-table.R predread-all > dev/predread-log/suite-all.md 2>&1

step openblas
sed -i 's#^frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-incl-thres-draws.R.*#&#' dev/predread-log/jobs-ob.txt
printf '%s\n' 'frmtmb tests/testthat/test-ordinal-fitted.R' \
  'frmtmb tests/testthat/test-review-v29.R' \
  'frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-sampling-ported.R' \
  >> dev/predread-log/jobs-ob.txt
sort -u dev/predread-log/jobs-ob.txt -o dev/predread-log/jobs-ob.txt
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/predread-log/ob26-build.txt 2>&1
bash dev/ciharden-scan.sh predread-ob26 ob0.3.26-lapack-t4 "$L" 2 \
  dev/predread-log/jobs-ob.txt > dev/predread-log/ob26-run.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > dev/predread-log/ob32-build.txt 2>&1
bash dev/ciharden-scan.sh predread-ob32 ob0.3.32-lapack-t4 "$L" 2 \
  dev/predread-log/jobs-ob.txt > dev/predread-log/ob32-run.txt 2>&1
rm -rf dev/ciharden-out

step audits
"$RS" dev/predread-audit.R "$L" after > dev/predread-log/audit-after.txt 2>&1 &
"$RS" dev/predread-audit-sample.R "$L" after \
  > dev/predread-log/audit-sample-after.txt 2>&1 &
"$RS" dev/predread-audit2.R "$L" base > dev/predread-log/audit2-base.txt 2>&1 &
wait
"$RS" dev/predread-audit-table.R > dev/predread-log/audit-table.md 2>&1

step checks
(
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
export R_LIBS="$L;$B;$U"
ck=dev/predread-check
rm -rf "$ck"
mkdir -p "$ck/core" "$ck/sample" "$ck/latent"
(cd "$ck/core" && "$RC" CMD build ../../.. > build.txt 2>&1 &&
   "$RC" CMD check --as-cran frmtmb_*.tar.gz > check.txt 2>&1) &
(cd "$ck/sample" && "$RC" CMD build ../../../extensions/frmtmb.sample \
   > build.txt 2>&1 &&
   "$RC" CMD check --as-cran frmtmb.sample_*.tar.gz > check.txt 2>&1) &
# fix round 2 changed one test of frmtmb.latent
(cd "$ck/latent" && "$RC" CMD build ../../../extensions/frmtmb.latent \
   > build.txt 2>&1 &&
   "$RC" CMD check --as-cran frmtmb.latent_*.tar.gz > check.txt 2>&1) &
wait
)

step sessions
NOT_CRAN=true FRMTMB_RELLIB="$L" bash dev/release/run-sessions.sh \
  predread-sess 8 > dev/predread-log/run-sessions.txt 2>&1

step equiv
rm -f dev/predread-equiv-out/p*.tsv dev/predread-equiv-out/p*.rds
FRMTMB_EQUIV_LIB="$L;$B" bash dev/equiv-run.sh 8 dev/predread-equiv-out \
  > dev/predread-log/equiv-run3.txt 2>&1
"$RS" dev/predread-equiv-diff.R > dev/predread-log/equiv-diff.txt 2>&1
step done
