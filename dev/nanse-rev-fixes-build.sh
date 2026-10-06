#!/bin/sh
# Reviewer: build dev/nanse-rev-fixes-table.R from the header and the
# body (line 28 on) of lane fixes' dev/fixes-p2-table.R.
W=C:/Users/adf44/source/r/frmtmb-wt-nanse
cat $W/dev/nanse-rev-fixes-head.R > $W/dev/nanse-rev-fixes-table.R
tail -n +28 C:/Users/adf44/source/r/frmtmb-wt-fixes/dev/fixes-p2-table.R >> $W/dev/nanse-rev-fixes-table.R
