#!/usr/bin/env bash
# md5 of every generated file of the eight packages (Rd, NAMESPACE and
# frmtmb.sample's copy of helper-brms-suite.R), to show that a second
# roxygen pass changes none.
#   bash dev/rel070-rox-md5.sh <out>
cd /c/Users/adf44/source/r/frmtmb-wt-release || exit 1
{ ls man/*.Rd NAMESPACE extensions/*/man/*.Rd extensions/*/NAMESPACE \
     extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R; } |
  sort | xargs md5sum > "$1"
wc -l < "$1"
