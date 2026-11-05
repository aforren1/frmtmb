#!/usr/bin/env bash
# Round 2's short check: reruns and the rule study, 6 at a time.
O=dev/review-refuse-out/f2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
{
for s in f1-classify f1-classify2 misc4 f1-misc f1-idx f1-names f2-items; do echo "$s x"; done
for s in f1-vm2 f1-vm; do echo "$s x"; done
for f in Gamma exponential beta exgaussian student lognormal von_mises gaussian; do echo "f2-rule table1 $f"; done
for s in 0 0.02 0.05 0.1 0.2 0.4 0.7; do echo "f2-rule table2 $s"; done
echo "f2-rule vm"
echo "f2-rule cost"
} > $O/jobs.txt
cat $O/jobs.txt | xargs -P 6 -L1 bash -c 'n=$(echo "$0-$1-$2" | tr -d " "); "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/review-refuse-$0.R $1 $2 > dev/review-refuse-out/f2/$n.log 2>&1'
touch $O/DONE
