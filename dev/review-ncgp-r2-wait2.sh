cd /c/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-out
until grep -q DONE r2-remedyU1.txt && grep -q DONE r2-remedyU3.txt && grep -q DONE r2-sampcmp-lane-1.txt && grep -q DONE r2-sampcmp-lane-2.txt && grep -q DONE r2-equi3.txt && grep -q Status: /c/Users/adf44/source/r/rev-ncgp-src/check/frmtmb/check.txt && grep -q DONE ngbrms-bern.txt; do sleep 30; done; echo ok
