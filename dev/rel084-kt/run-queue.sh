#!/usr/bin/env bash
# rel084-kt: the consolidation's one-process queue beside the chain. Waits
# for the ordstarts driver (left running by the key-table agent) to write
# its done file, then runs the remaining lanes one R process at a time,
# then merge step 6 (test-compat.R without and with every extension).
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
D=$ROOT/dev/rel084-kt
until [ -f "$KLOG/ordstarts-done.txt" ]; do sleep 60; done
echo "QUEUE ordstarts done $(date '+%F %T')" >> "$KLOG/queue.txt"
bash $D/run-nbrobust.sh
echo "QUEUE nbrobust done $(date '+%F %T')" >> "$KLOG/queue.txt"
bash $D/run-rest.sh > "$KLOG/run-rest.txt" 2>&1
echo "QUEUE rest done $(date '+%F %T')" >> "$KLOG/queue.txt"
krun compat-none.log $D/compat-ext.R none
krun compat-all.log $D/compat-ext.R all
echo "QUEUE END $(date '+%F %T')" >> "$KLOG/queue.txt"
