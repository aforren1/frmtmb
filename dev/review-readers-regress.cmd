@echo off
rem Review of lane readers: both arms of dev/review-readers-regress.R pinned
rem to cores 0-1 (P-cores), dev/rtmb-pitfalls.md item 21.
cd /d C:\Users\adf44\source\r\frmtmb-wt-readers
set R_PARALLEL_PORT=21820
set REV_ARM=before
start "" /affinity 3 /wait /b "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev/review-readers-regress.R > dev\review-readers-out\regress-before.txt 2>&1
set R_PARALLEL_PORT=21821
set REV_ARM=after
start "" /affinity 3 /wait /b "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev/review-readers-regress.R > dev\review-readers-out\regress-after.txt 2>&1
"C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev/review-readers-regress-cmp.R > dev\review-readers-out\regress-cmp.txt 2>&1
