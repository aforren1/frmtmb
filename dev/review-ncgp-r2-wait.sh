cd /c/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-out
until grep -q DONE r2-ng-lane.txt && grep -q DONE r2-remedy3.txt && grep -q DONE r2-fc.txt && grep -q DONE floor-40.txt && grep -q DONE floor-100.txt && grep -q DONE r2-edge2-lane.txt && grep -q DONE r2-edge2-base.txt && grep -q "RME\|ERR" mem-lane-300.txt; do sleep 30; done; echo ok
