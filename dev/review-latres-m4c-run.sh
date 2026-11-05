#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for d in gid_l bid_l ginv_l pid_l pid_20 glog_w blog_w; do for s in 2001 2021; do e=$((s+19)); echo "$d $s $e dev/review-latres-out/m4b/$d-$s.tsv"; done; done | xargs -P 4 -L 1 "$R" dev/review-latres-m4b.R >> dev/review-latres-out/m4b.log 2>&1
