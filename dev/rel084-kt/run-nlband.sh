#!/usr/bin/env bash
# rel084-kt: nlband key tables on rellib-r23, one R process at a time.
# The survey's 25 jobs run serially here (the lane ran dev/nlband-queue.sh
# with N processes).
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
export NLBAND_ARM=after
krun nlband-fix1.log $K/nlband-fix1.R sep,stud,band
krun nlband-repro.log $K/nlband-repro.R
krun nlband-studg.log $K/nlband-short-studg.R x
krun nlband-sweep.log $K/nlband-sweep.R
while IFS='|' read -r log envs cmd; do
  [ -z "$log" ] && continue
  mkdir -p "$(dirname "$log")"
  t0=$(date +%s)
  env $envs "$RS" $cmd > "$log" 2>&1
  echo "JOB END $log exit $? wall $(( $(date +%s) - t0 ))s" >> $KLOG/nlband-survey-jobs.txt
done < $K/nlband-acc-afterrel.jobs
bash $K/nlband-accsum.sh before afterrel > $KLOG/nlband-accsum-before.txt 2>&1
bash $K/nlband-accsum.sh after4 afterrel > $KLOG/nlband-accsum-after4.txt 2>&1
echo ALLDONE > $KLOG/nlband-done.txt
