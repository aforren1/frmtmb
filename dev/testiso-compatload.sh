#!/usr/bin/env bash
# compatload over every preload set; $1 = tag, TESTISO_DIR/TESTISO_LIB/TESTISO_ATTACH from the caller
cd "$(dirname "$0")/.."
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=dev/testiso-log/compatload-$1; mkdir -p $out
all="frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline"
"$RS" dev/testiso-compatload.R > $out/none.txt 2>&1 &
for p in $all; do "$RS" dev/testiso-compatload.R $p > $out/$p.txt 2>&1 & done
"$RS" dev/testiso-compatload.R $all > $out/all.txt 2>&1 &
wait
for f in $out/*.txt; do echo "== $(basename $f .txt)"; grep -a -E "^test-compat|^registries" $f; done
