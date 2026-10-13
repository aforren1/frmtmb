#!/usr/bin/env bash
# Lane seflags, punch round 2b: the reruns on the final build, in turn.
# Waits for every job it starts.
#   bash dev/seflags-p2b-batch.sh
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
P=dev/seflags-log/p2b
echo "start $(date)"
bash dev/seflags-p2b-falarm.sh &
for a in lane m9; do
  "$RS" dev/seflags-p2b-scan.R $a > $P/run-$a/scan-ref-$a.txt 2>&1 &
  OPENBLAS_NUM_THREADS=2 "$OB" dev/seflags-p2b-scan.R $a > $P/run-$a/scan-ob-$a.txt 2>&1 &
  "$RS" dev/seflags-p2b-paths.R $a > $P/paths-ref-$a.txt 2>&1 &
  OPENBLAS_NUM_THREADS=2 "$OB" dev/seflags-p2b-paths.R $a > $P/paths-ob-$a.txt 2>&1 &
done
wait
echo "falarm, scan, paths done $(date)"
bash dev/seflags-revrun.sh p2b ref > $P/revrun-ref.txt 2>&1
bash dev/seflags-revrun.sh p2bob ob26 > $P/revrun-ob.txt 2>&1
echo "revrun done $(date)"
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp}"
bash dev/seflags-attacks.sh p2b ref 12 > $P/attacks-ref.txt 2>&1
bash dev/seflags-attacks.sh p2bob ob26 12 > $P/attacks-ob.txt 2>&1
echo "attacks done $(date)"
