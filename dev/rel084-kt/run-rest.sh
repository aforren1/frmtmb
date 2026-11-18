#!/usr/bin/env bash
# rel084-kt: the remaining lanes, one after another (each one R process at
# a time).
D=/c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt
bash $D/run-gpselimit.sh
bash $D/run-slopetrig.sh
bash $D/run-squant.sh
bash $D/run-sparsese.sh
bash $D/run-quadesc.sh
bash $D/run-nlband.sh
bash $D/run-drawdir.sh
echo RESTDONE
