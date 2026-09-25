#!/bin/bash
# Reviewer 2: fresh seeds for item 6, after the seeds-1.. run's last two
# fits (PIDs given) have exited, so that at most two fits run at once.
cd /c/Users/adf44/source/r/frmtmb-wt-phase3b
for p in "$@"; do while kill -0 $p 2>/dev/null; do sleep 15; done; done
echo "$(date +%T) orphans done, launching fresh seeds" >> dev/phase3b-review2/cov-launch2.log
bash dev/phase3b-review2/r2-cov-driver.sh base C $(seq 2001 2 2069) &
sleep 30
bash dev/phase3b-review2/r2-cov-driver.sh base D $(seq 2002 2 2070) &
wait
