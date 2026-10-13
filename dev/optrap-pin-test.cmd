@echo off
rem Lane optrap, item 3: one test file in one R process pinned to a core
rem mask.   optrap-pin-test.cmd <hexmask> <label> <lib> <test file>
cd /d C:\Users\adf44\source\r\frmtmb-wt-optrap
set NOT_CRAN=true
start "" /b /wait /affinity %1 "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev/optrap-testfile.R %3 frmtmb %4 > dev\optrap-log\pin-%2.txt 2>&1
