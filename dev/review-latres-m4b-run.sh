#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/review-latres-out/m4b2
{
for d in ginv_b1 ginv_b2 ginv_b3 ginv_s5 gident_b gid_s5 ig_b1 ig_b2 ig_s5 pid_b pid_20 pid_8 bid_b glog_c glog_s5 glog_w blog_w gid_l bid_l ginv_l pid_l bb_T zoib_T; do for s in 2001 2021; do e=$((s+19)); echo "$d $s $e dev/review-latres-out/m4b2/$d-$s.tsv"; done; done
for d in gpost_b wpost_b epost_b igpost_b bepost_b psqrt_b bb_post_lo bb_post_b zoib_post_lo zoib_post_b; do for s in 2001 2011; do e=$((s+9)); echo "$d $s $e dev/review-latres-out/m4b2/$d-$s.tsv"; done; done
} | xargs -P 7 -L 1 "$R" dev/review-latres-m4b.R > dev/review-latres-out/m4b2.log 2>&1
