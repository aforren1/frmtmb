#!/usr/bin/env bash
# The SE false-alarm battery (lane illse's 41 models, seeds 1:20, and its
# bitwise set) on the merged build against rellib-r10, under the
# reference BLAS and OpenBLAS 0.3.26 and 0.3.32 (BLAS and LAPACK, 4
# threads), one configuration after another (dev/rel073-sf-run.sh).
#   bash dev/rel073-sf-all.sh
cd /c/Users/adf44/source/r/frmtmb-wt-release || exit 1
TAG=ref CFGS=ref N=8 bash dev/rel073-sf-run.sh
TAG=ob26 CFGS=ob OBV=0.3.26 N=8 bash dev/rel073-sf-run.sh
TAG=ob32 CFGS=ob OBV=0.3.32 N=8 bash dev/rel073-sf-run.sh
echo "SF ALL DONE $(date +%T)"
