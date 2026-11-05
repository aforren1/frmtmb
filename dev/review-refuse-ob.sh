#!/usr/bin/env bash
# OpenBLAS 0.3.32 with LAPACK, on the test files the changed code reaches.
set -u
W="$(pwd)"
OB="$W/dev/review-refuse-out/ob"
bash "$OB/dev/ciharden-openblas.sh" 0.3.32 lapack > "$OB/build.log" 2>&1
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp}"
export TREE="$W"
export FRMTMB_STAN_CACHE="$W/dev/review-refuse-out/stan-cache"
bash "$OB/dev/ciharden-scan.sh" reviewob ob0.3.32-lapack-t4 C:/Users/adf44/source/r/wt-refuse-lib 3 "$OB/list.txt" > "$OB/scan.log" 2>&1
