@echo off
rem Short check of lane quadseam2's fix round, item 2: SYMSMALL seed 1
rem (dev/quadseam2-symsmall.R) under OpenBLAS 0.3.26 with LAPACK,
rem OPENBLAS_CORETYPE=Haswell, on rellib-r21 and on the lane library,
rem each pinned to the same P-core mask (rtmb-pitfalls.md item 21).
rem Run from the worktree after dev/ciharden-openblas.sh 0.3.26 lapack.
set OPENBLAS_CORETYPE=Haswell
set NOT_CRAN=true
set RO=dev\ciharden-out\Rob0.3.26-lapack\bin\Rscript.exe
start "" /b /wait /affinity 0C00 %RO% dev/quadseam2-symsmall.R base "1" > dev\reviews\quadseam2-short-log\ob-symsmall-base.txt 2>&1
start "" /b /wait /affinity 0C00 %RO% dev/quadseam2-symsmall.R C:/Users/adf44/source/r/wt-quadseam2-lib "1" > dev\reviews\quadseam2-short-log\ob-symsmall-lib.txt 2>&1
start "" /b /wait /affinity 0C00 %RO% dev/quadseam2-symsmall.R base "1" > dev\reviews\quadseam2-short-log\ob-symsmall-base2.txt 2>&1
echo OB DONE
