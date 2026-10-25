#!/bin/sh
# Guardrail-1 classification of every brms name, three arms: base
# (rellib-r12) core alone, the lane core alone, the lane with
# frmtmb.spatial loaded. dev/spatial-guard-report.R does the work.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd C:/Users/adf44/source/r/frmtmb-wt-sarfcor
SPATIAL_LIBS="C:/Users/adf44/source/r/rellib-r12" SPATIAL_PKGS="" "$RS" dev/spatial-guard-report.R > dev/sarfcor-out/guard-base.txt 2>&1
SPATIAL_LIBS="C:/Users/adf44/source/r/wt-sarfcor-lib;C:/Users/adf44/source/r/rellib-r12" SPATIAL_PKGS="" "$RS" dev/spatial-guard-report.R > dev/sarfcor-out/guard-lane.txt 2>&1
SPATIAL_LIBS="C:/Users/adf44/source/r/wt-sarfcor-lib;C:/Users/adf44/source/r/rellib-r12" SPATIAL_PKGS="frmtmb.spatial" "$RS" dev/spatial-guard-report.R > dev/sarfcor-out/guard-lane-spatial.txt 2>&1
