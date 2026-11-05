#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
{
for N in 64 200 1000; do
  echo "sym20 1 $N"; echo "symt 1 $N"; echo "bigsd_bern 3 $N"; echo "t_nu5 1 $N";
  echo "re_sigma 1 $N"; echo "nearsep 2 $N"; echo "smallpois 2 $N"; echo "wtbin 1 $N";
  echo "t_outlier 2 $N"; echo "unbal_bern 1 $N"; echo "cloglog 5 $N"; echo "gaussloglink 5 $N";
done
} | xargs -P 4 -L 1 bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" review-impsamp-mcse.R "$0" "$1" "$2" 100 > review-impsamp-out/mcse-log/"$0"-"$1"-"$2".log 2>&1'
echo done
