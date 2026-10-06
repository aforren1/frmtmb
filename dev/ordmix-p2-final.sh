#!/usr/bin/env bash
# Punch round 2, the final runs on the final lane build, in turn: the
# core and frmtmb.sample suites (one file per process, every gate on),
# then every warning driver (dev/ordmix-p2-drivers.sh).
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
export TMP="C:/Users/adf44/AppData/Local/Temp/1" TEMP="C:/Users/adf44/AppData/Local/Temp/1"
bash "$ROOT/dev/ordmix-run-par.sh" lane p2final "$ROOT/dev/ordmix-p2-jobs-suites.txt" 20
bash "$ROOT/dev/ordmix-p2-drivers.sh"
echo FINAL DONE
