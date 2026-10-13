#!/usr/bin/env bash
# Reviewer of lane seflags (check of punch 2c): the rev4 attack scripts
# (whole, cum, glm, illpred, illid, routes, bitwise) and the B1 to B3
# and degree-7 scripts on the lane (and base where a reference run is
# needed), reference BLAS and OpenBLAS 0.3.26 LAPACK, test-se-flags.R
# both BLAS, and the false-alarm battery. Waits for every job.
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
L=dev/seflags-rev5-log
for c in ref ob; do
  R="$RS"; [ $c = ob ] && R="$OB"
  for s in whole cum illpred routes; do
    "$R" dev/seflags-rev4-$s.R lane < /dev/null > $L/$s-$c-lane.txt 2>&1 &
  done
  "$R" dev/seflags-rev4-illid.R lane $L/illid-$c-lane.tsv < /dev/null > $L/illid-$c-lane.txt 2>&1 &
  "$R" dev/seflags-rev4-glm.R lane $L/glm-$c-lane.tsv < /dev/null > $L/glm-$c-lane.txt 2>&1 &
  "$R" dev/seflags-rev4-bitwise.R lane $L/bitwise-$c-lane.rds < /dev/null > $L/bitwise-$c-lane.txt 2>&1 &
  "$R" dev/seflags-rev2-join.R lane < /dev/null > $L/join-$c-lane.txt 2>&1 &
  "$R" dev/seflags-rev2-deg7.R lane < /dev/null > $L/deg7-$c-lane.txt 2>&1 &
  "$R" dev/seflags-rev3-pred3x.R lane < /dev/null > $L/pred3x-$c-lane.txt 2>&1 &
  "$R" dev/seflags-b1ref.R < /dev/null > $L/b1ref-$c-lane.txt 2>&1 &
  "$R" dev/seflags-testfile.R lane tests/testthat/test-se-flags.R < /dev/null > $L/test-se-flags-$c.txt 2>&1 &
  wait
done
"$RS" dev/seflags-rev2-sd2.R lane < /dev/null > $L/sd2-ref-lane.txt 2>&1 &
"$RS" dev/seflags-rev3-deg7exact.R lane < /dev/null > $L/deg7exact-ref-lane.txt 2>&1 &
"$RS" dev/seflags-rev3-deg7pred.R lane < /dev/null > $L/deg7pred-ref-lane.txt 2>&1 &
"$RS" dev/seflags-rev4-whole.R base < /dev/null > $L/whole-ref-base.txt 2>&1 &
"$RS" dev/seflags-rev4-cum.R base < /dev/null > $L/cum-ref-base.txt 2>&1 &
wait
for c in ref ob; do
  "$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-$c-lane.rds dev/seflags-rev4-log/bitwise-$c-base.rds > $L/bitcmp-$c.txt 2>&1
done
echo all done
