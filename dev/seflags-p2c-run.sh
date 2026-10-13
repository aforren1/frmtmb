#!/usr/bin/env bash
# Lane seflags, punch round 2c: the review's dev/seflags-rev4-* scripts
# (illpred, illid, whole, cum, routes, gain) on the lane, the 0.69.0 base
# and the lane's merge with optrap (copies whose "lane" library is the
# merge build), reference BLAS and OpenBLAS 0.3.26 LAPACK, into
# dev/seflags-log/p2c/rev4. Waits for every job.
#   bash dev/seflags-p2c-run.sh
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
L=dev/seflags-log/p2c/rev4
S=dev/seflags-out/p2c-scripts
SW=C:/Users/adf44/source/r/frmtmb-wt-seflags/$S
mkdir -p $L $S
sed "s#C:/Users/adf44/source/r/wt-seflags-lib#C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-merge-out/lib#" \
  dev/seflags-rev4-common.R > $S/common.R
for s in illpred illid whole cum routes gain; do
  sed "s#C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev4-common.R#$SW/common.R#" \
    dev/seflags-rev4-$s.R > $S/$s.R
done
for c in ref ob; do
  R="$RS"; [ $c = ob ] && R="$OB"
  for a in lane base m9; do
    if [ $a = m9 ]; then P=$S/; arm=lane; else P=dev/seflags-rev4-; arm=$a; fi
    "$R" ${P}illpred.R $arm > $L/illpred-$c-$a.txt 2>&1 &
    "$R" ${P}illid.R $arm $L/illid-$c-$a.tsv > $L/illid-$c-$a.txt 2>&1 &
    "$R" ${P}whole.R $arm > $L/whole-$c-$a.txt 2>&1 &
    "$R" ${P}cum.R $arm > $L/cum-$c-$a.txt 2>&1 &
    "$R" ${P}routes.R $arm > $L/routes-$c-$a.txt 2>&1 &
    [ $a = base ] || "$R" ${P}gain.R $arm $L/gain-$c-$a.tsv > $L/gain-$c-$a.txt 2>&1 &
  done
  wait
done
echo all done
