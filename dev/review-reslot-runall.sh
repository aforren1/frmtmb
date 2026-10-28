#!/bin/bash
# Reviewer's queue: one test file per R process, 4 at a time.
# Each line of the queue file: arm package path gated(0/1)
W=C:/Users/adf44/source/r/frmtmb-wt-reslot
mkdir -p $W/dev/review-reslot-out
xargs -P 4 -L 1 bash $W/dev/review-reslot-run1.sh < "$1"
