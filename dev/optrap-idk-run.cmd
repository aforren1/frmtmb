@echo off
rem Lane optrap, item 3: one run of dev/optrap-idk.R pinned to a core mask.
rem   optrap-idk-run.cmd <hexmask> <label> [lib]
cd /d C:\Users\adf44\source\r\frmtmb-wt-optrap
start "" /b /wait /affinity %1 "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev/optrap-idk.R "%3" %2 > dev\optrap-log\idk-%2.txt 2>&1
