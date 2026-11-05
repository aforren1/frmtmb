cd /c/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-out
until [ -f r3-queue.done ] && grep -q DONE r3-dist6.txt && grep -q DONE r3-logs.txt; do sleep 30; done; echo ok
