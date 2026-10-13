@echo off
rem Reviewer of lane optrap, re-check: one Rscript pinned to a core mask.
rem   optrap-rev2-pin.cmd <hexmask> <log name> <script> [args...]
cd /d C:\Users\adf44\source\r\frmtmb-wt-optrap
set NOT_CRAN=true
set MASK=%1
set NAME=%2
shift
shift
start "" /b /wait /affinity %MASK% "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" %1 %2 %3 %4 %5 %6 %7 %8 > dev\optrap-rev2-log\%NAME%.txt 2>&1
