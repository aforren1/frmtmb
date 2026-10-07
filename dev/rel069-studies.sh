#!/usr/bin/env bash
# The reruns the four reviews ask for on the merged tree, against the
# 0.69.0 release library (rellib-r7), reference BLAS. Logs in
# dev/rel069-log/studies/. The attack set has its own driver
# (dev/rel069-attacks.sh).
#
#   bash dev/rel069-studies.sh
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
cd "$REL"
L=dev/rel069-log/studies
mkdir -p "$L"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R7=C:/Users/adf44/source/r/rellib-r7
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$REL/dev/stan-cache
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
# ciharden: test-perf.R's headroom and its guard
"$RS" dev/rel069-perfseeds.R glm 1001 15 > "$L/perf-glm-1001.txt" 2>&1 &
"$RS" dev/rel069-perfseeds.R glm 2001 15 > "$L/perf-glm-2001.txt" 2>&1 &
"$RS" dev/rel069-perfseeds.R glmm 3001 10 > "$L/perf-glmm-3001.txt" 2>&1 &
for v in asis noprofmem bytesNA bytesZero noloop; do
  "$RS" dev/rel069-perfguard.R "$v" > "$L/perfguard-$v.txt" 2>&1 &
done
# optima: the mo() study, 200 seeds in four chunks, and check_laplace()
for c in 1 2 3 4; do
  a=$(( (c - 1) * 50 + 1 )); b=$(( c * 50 ))
  "$RS" dev/rel069-mo-study.R rel "$a:$b" "$L/mo-rel-$c.tsv" \
    > "$L/mo-rel-$c.log" 2>&1 &
done
"$RS" dev/rel069-laplace.R > "$L/laplace.txt" 2>&1 &
# setier: singular study, separation study, the quadrature trap and
# quadrature boundary checks
"$RS" dev/setier-singular.R "$R7" "$L/sing-rel.tsv" > "$L/sing-rel.log" 2>&1 &
"$RS" dev/setier-sep.R "$R7" "$L/sep-rel.tsv" > "$L/sep-rel.log" 2>&1 &
"$RS" dev/setier-rev3-trapq.R "$R7" > "$L/trapq-rel.txt" 2>&1 &
"$RS" dev/setier-rev3-quad.R "$R7" > "$L/quad-rel.txt" 2>&1 &
wait
echo STUDIES DONE
