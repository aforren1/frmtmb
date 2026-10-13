#!/usr/bin/env bash
# Lane seflags, punch round 1: the B1 to B3 and degree-7/8 checks on two
# arms, this lane alone ("lane") and the trial merge with lane optrap
# ("merge", dev/seflags-merge.sh), under one R build. The review's
# scripts are copied with their "lane" library pointed at the merge for
# the merge arm. Logs: dev/seflags-log/merge-<tag>/<script>-<arm>.txt.
#   bash dev/seflags-merge-run.sh <tag> <ref|ob26>
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
tag="$1"; cfg="$2"
case "$cfg" in
  ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe" ;;
  ob26) RS="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
        export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 ;;
esac
ML=C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-merge-out/lib
LL=C:/Users/adf44/source/r/wt-seflags-lib
S=$W/dev/seflags-merge-out/scripts
# the same directory as R reads it
SW=C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-merge-out/scripts
L=$W/dev/seflags-log/merge-$tag
mkdir -p "$S" "$L"
sed "s#C:/Users/adf44/source/r/wt-seflags-lib#$ML#" dev/seflags-rev-common.R \
  > "$S/common.R"
for s in deg7 pred3 join sd2 sd fx1theta; do
  sed "s#C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R#$SW/common.R#" \
    dev/seflags-rev-$s.R > "$S/rev-$s.R"
done
sed "s#C:/Users/adf44/source/r/wt-seflags-lib#$ML#" dev/seflags-b1ref.R \
  > "$S/b1ref.R"
# the re-check's prediction script, its "lane" arm pointed at the merge
sed "s#C:/Users/adf44/source/r/wt-seflags-lib#$ML#" dev/seflags-rev2-common.R \
  > "$S/common2.R"
sed "s#C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev2-common.R#$SW/common2.R#" \
  dev/seflags-rev2-illpred.R > "$S/illpred.R"
for arm in lane merge; do
  if [ $arm = lane ]; then
    pre=dev/seflags-rev-; lib=$LL; b1=dev/seflags-b1ref.R
    ip=dev/seflags-rev2-illpred.R
  else
    pre=$S/rev-; lib=$ML; b1=$S/b1ref.R; ip=$S/illpred.R
  fi
  "$RS" "$ip" lane > "$L/illpred-$arm.txt" 2>&1 &
  "$RS" dev/seflags-illscan.R "$lib" > "$L/illscan-$arm.txt" 2>&1 &
  for d in 7 8; do
    "$RS" ${pre}deg7.R lane $d > "$L/deg$d-$arm.txt" 2>&1 &
  done
  for s in pred3 join sd2 sd fx1theta; do
    "$RS" ${pre}$s.R lane > "$L/$s-$arm.txt" 2>&1 &
  done
  "$RS" "$b1" > "$L/b1ref-$arm.txt" 2>&1 &
  "$RS" dev/seflags-testfile.R $arm tests/testthat/test-se-flags.R \
    > "$L/flags-$arm.txt" 2>&1 &
  wait
done
echo MERGE RUN DONE
