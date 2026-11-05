cd /c/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-out
until [ -f r3-queue.done ] && grep -q DONE r3-mtest.txt && grep -q DONE r3-ng-lane-a.txt && grep -q DONE r3-ng-base-a.txt; do sleep 30; done; echo ok
