#!/usr/bin/env bash
# Reviewer of lane seflags (re-check): the round-0 review scripts on the
# punch build (lane), the base and the trial merge with optrap's current
# tree (op2), reference BLAS and OpenBLAS 0.3.26 LAPACK. Waits for every
# job before it exits.
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
L=dev/seflags-rev2-log/A
mkdir -p $L
run() { "$1" dev/seflags-rev2-$2.R $3 > $L/$2-$4-$3.txt 2>&1; }
for a in lane base op2; do
  run "$OB" pred3 $a ob &
  run "$OB" fx1theta $a ob &
  run "$RS" join $a ref &
  run "$OB" join $a ob &
  run "$RS" sd2 $a ref &
  run "$RS" deg7 $a ref &
  run "$OB" deg7 $a ob &
done
for a in lane base; do
  run "$RS" sd $a ref &
  run "$RS" pred $a ref &
  run "$RS" pred2 $a ref &
  run "$OB" pred2 $a ob &
done
run "$RS" neg lane ref &
"$RS" dev/seflags-b1ref.R > $L/b1ref-ref-lane.txt 2>&1 &
"$RS" dev/seflags-rev2-b1ref-op2.R > $L/b1ref-ref-op2.txt 2>&1 &
"$OB" dev/seflags-rev2-b1ref-op2.R > $L/b1ref-ob-op2.txt 2>&1 &
for lib in C:/Users/adf44/source/r/wt-seflags-lib C:/Users/adf44/source/r/seflags-rev-lib/op2 C:/Users/adf44/source/r/rellib-r7; do
  t=$(basename $lib)
  "$RS" dev/seflags-illscan.R $lib > $L/illscan-ref-$t.txt 2>&1 &
done
"$OB" dev/seflags-illscan.R C:/Users/adf44/source/r/seflags-rev-lib/op2 > $L/illscan-ob-op2.txt 2>&1 &
wait
echo all done
