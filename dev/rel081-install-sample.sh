#!/usr/bin/env bash
# Roxygenise and reinstall frmtmb.sample alone into rellib-r20, after the
# release review's fixes (dev/rel081-release-fix3.R), which touched only
# frmtmb.sample's code. Core and the other extensions keep their install
# of dev/rel081-install.sh.
#   bash dev/rel081-install-sample.sh > dev/rel081-log/install-sample.txt 2>&1
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel081
LIB=C:/Users/adf44/source/r/rellib-r20
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export R_LIBS="$LIB;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export REL081_LIB="$LIB"
"$RS" "$ROOT/dev/rel081-rox.R" frmtmb.sample 2>&1 | grep -E '^== rox|ROX DONE|Error|Warn|Writing|Deleting'
echo "== INSTALL frmtmb.sample $(date +%T)"
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$ROOT/extensions/frmtmb.sample" 2>&1 |
  grep -E 'DONE|ERROR|Error|WARN' | head -8
"$RS" -e ".libPaths(c('$LIB')); cat('installed', format(packageVersion('frmtmb.sample')), format(file.mtime(system.file('DESCRIPTION', package = 'frmtmb.sample'))), '\n')"
echo "INSTALL-SAMPLE END"
