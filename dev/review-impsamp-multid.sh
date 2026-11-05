#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
xargs -P 4 -L 1 bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" review-impsamp-multid.R "$0" "$1" > review-impsamp-out/multid-log/"$0"-"$1".log 2>&1' < review-impsamp-out/multid-jobs.txt
echo done
