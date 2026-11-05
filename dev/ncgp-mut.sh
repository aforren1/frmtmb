#!/usr/bin/env bash
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/ncgp-out/mut
for m in control absolute nugget6 nonull noscan nowarm nounion nowork quiet_all quiet_headroom floorless nopolish norowcheck noorder sorted nothcols fd_dtheta; do
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 1; done
  ( "$RS" dev/ncgp-mut.R $m > dev/ncgp-out/mut/$m.txt 2>&1 ) &
done
wait
grep -h NCGPM dev/ncgp-out/mut/*.txt > dev/ncgp-out/mut.txt
echo MUT-DONE
