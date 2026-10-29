#!/usr/bin/env bash
# The floors of the 0.78.0 release, measured. Adapted from
# dev/rel077-floor.sh (the successor of dev/rel073-floor.sh), with its
# method: nothing is installed. Each arm loads a COPY of an installed
# build (none of these packages has compiled code) in front of the core
# it is measured on.
#
# frmtmb.sample (floor raised to frmtmb (>= 0.78.0)):
#   old: rellib-r16's build, its floor lowered to 0.61.0 in DESCRIPTION
#        and Meta/package.rds, in front of rellib-r15 (frmtmb 0.77.0);
#        if it loads, every test file of its suite runs there (the
#        release tree's files), one per process, every gate on. On 0.77.0
#        core it cannot call predict_rescor_draw() or
#        rescor_row_loglik(resp =), and set_prior() refuses the me()
#        classes its default priors write.
#   new: rellib-r16 as built.
#   converse: rellib-r15's frmtmb.sample 0.24.0 in front of rellib-r16.
# The seven other extensions (no file of theirs changed): every test file
# of each, from the release tree, on rellib-r15 alone (frmtmb 0.77.0 with
# the extension's own 0.77.0 build), every gate on. Their counts on the
# new core are the gated tier's; equal counts mean no bump.
#
#   bash dev/rel078-floor.sh [N] [part]    part: sample (default all), others
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/12}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/12}"
export GIT_OPTIONAL_LOCKS=0
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel078
SCR=C:/Users/adf44/source/r/rel078-floor-lib
R15=C:/Users/adf44/source/r/rellib-r15
R16=C:/Users/adf44/source/r/rellib-r16
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
N="${1:-8}"
part="${2:-all}"
cd "$ROOT" || exit 1
mkdir -p "$SCR"
echo "== ac4906bd (0.77.0): the seams the new floor rests on"
# predict_rescor_draw() existed as an internal at 0.77.0; the floor rests
# on its export (floor-1.txt printed the definition count, 1 and 1)
printf '  predict_rescor_draw exported in NAMESPACE: ac4906bd %s, this tree %s\n' \
  "$(git show ac4906bd:NAMESPACE | grep -c 'predict_rescor_draw')" \
  "$(grep -c 'predict_rescor_draw' NAMESPACE)"
printf '  rescor_row_loglik(fit, dpv, resp = NULL): ac4906bd %s, this tree %s\n' \
  "$(git show ac4906bd:R/predict-brms.R | grep -c 'rescor_row_loglik <- function(fit, dpv, resp')" \
  "$(grep -c 'rescor_row_loglik <- function(fit, dpv, resp' R/predict-brms.R)"
printf '  "sdme" in R/priors.R: ac4906bd %s, this tree %s\n' \
  "$(git show ac4906bd:R/priors.R | grep -c '"sdme"')" \
  "$(grep -c '"sdme"' R/priors.R)"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
lower_floor() {   # <installed package dir>
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$1/DESCRIPTION"
  "$RS" -e "f <- '$1/Meta/package.rds'; m <- readRDS(f); for (k in c('Depends', 'Imports')) m[[k]] <- lapply(m[[k]], function(d) { if (identical(d\$name, 'frmtmb') && !is.null(d\$version)) d\$version <- package_version('0.61.0'); d }); saveRDS(m, f); cat('  floor in Meta/package.rds:', paste(vapply(c(m\$Depends, m\$Imports), function(d) paste(d\$name, format(d\$version)), ''), collapse = ', '), '\n')"
}
if [ "$part" = all ] || [ "$part" = sample ]; then
  p=frmtmb.sample
  v=$(grep '^Version:' "extensions/$p/DESCRIPTION" | cut -d' ' -f2 | tr -d '\r')
  fl=$(grep -o 'frmtmb (>= [0-9.]*)' "extensions/$p/DESCRIPTION" | head -1)
  echo "== $p $v, floor $fl"
  FL=$SCR/old-$p; rm -rf "${SCR:?}/old-$p"; mkdir -p "$FL"
  cp -r "$R16/$p" "$FL/"
  echo "-- old: rellib-r16's build over rellib-r15 (floor lowered)"
  lower_floor "$FL/$p"
  out=$(load_line "$FL" "$R15" "$p"); echo "$out"
  if echo "$out" | grep -q "LOADED"; then
    : > "dev/rel078-log/floor-$p.jobs"
    for f in extensions/$p/tests/testthat/test-*.R; do
      echo "$p $f" >> "dev/rel078-log/floor-$p.jobs"
    done
    ROUND_BASE=$R15 /usr/bin/bash dev/ciharden-scan.sh \
      "r078-floor-$p" ref "$FL" "$N" "dev/rel078-log/floor-$p.jobs"
  fi
  echo "-- new: rellib-r16 as built"
  load_line "$R16" "$R16" "$p"
  FC=$SCR/conv-$p; rm -rf "${SCR:?}/conv-$p"; mkdir -p "$FC"
  cp -r "$R15/$p" "$FC/"
  echo "-- converse: $p $(grep '^Version:' "$FC/$p/DESCRIPTION" | cut -d' ' -f2) from rellib-r15 over rellib-r16"
  load_line "$FC" "$R16" "$p"
fi
if [ "$part" = all ] || [ "$part" = others ]; then
  echo "== the seven other extensions: the release tree's tests on rellib-r15"
  : > dev/rel078-log/floor-others.jobs
  for e in coupling eam latent learn ode spatial survival; do
    for f in extensions/frmtmb.$e/tests/testthat/test-*.R; do
      echo "frmtmb.$e $f" >> dev/rel078-log/floor-others.jobs
    done
  done
  ROUND_BASE=$R15 /usr/bin/bash dev/ciharden-scan.sh r078-floor-others ref \
    base "$N" dev/rel078-log/floor-others.jobs
fi
rm -rf "$SCR"
echo "FLOOR DONE"
