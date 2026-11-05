#!/usr/bin/env bash
# Final pass of fix round 2, in order. Nothing is edited while it runs.
set -u
FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
  bash dev/gpgap-par.sh dev/gpgap-out/core-gated frmtmb > dev/gpgap-out/core-gated.log 2>&1
FRMTMB_BRMS_FIT_TESTS=true bash dev/gpgap-par.sh dev/gpgap-out/ext-gated \
  frmtmb.sample frmtmb.spatial frmtmb.survival frmtmb.coupling frmtmb.latent \
  frmtmb.learn frmtmb.ode frmtmb.eam > dev/gpgap-out/ext-gated.log 2>&1
NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-gpgap-lib \
  bash dev/release/run-sessions.sh gpgapsess 8 > dev/gpgap-out/sessions.log 2>&1
cp dev/release/gpgapsess.log dev/gpgap-out/gpgapsess.log
rm -rf dev/release/gpgapsess-sessions dev/release/gpgapsess.log
bash dev/gpgap-ob.sh
rm -rf dev/ciharden-out tests/testthat/_problems extensions/frmtmb.sample/tests/testthat/_problems
rm -f dev/gpgap-equiv-out/p022.tsv dev/gpgap-equiv-out/p023.tsv dev/gpgap-equiv-out/p094.tsv \
  dev/gpgap-equiv-out/p109.tsv dev/gpgap-equiv-out/p136.tsv
FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-gpgap-lib;C:/Users/adf44/source/r/rellib-r16" \
  bash dev/equiv-run.sh 5 dev/gpgap-equiv-out > dev/gpgap-out/equiv-run.txt 2>&1
for p in . extensions/frmtmb.spatial extensions/frmtmb.sample; do
  bash dev/gpgap-check.sh "$p" > "dev/gpgap-out/checkrun-$(basename "$p").txt" 2>&1
done
for p in frmtmb frmtmb.spatial frmtmb.sample; do
  s=$(grep -a "^Status" "dev/gpgap-check/$p/check.txt")
  n=$(grep -a "[.][.][.] .*NOTE$\|[.][.][.] NOTE$\|WARNING$\|ERROR$" "dev/gpgap-check/$p/check.txt" | cut -c1-70 | tr '\n' ';')
  echo "$s; $n" > "dev/gpgap-out/check-$p.txt"
done
echo FINAL2-DONE
