@echo off
start "" /affinity 1 /wait cmd /c ""C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev\review-priordraws-warmcmp.R old > dev\review-priordraws-log\warmcmp-old.txt 2>&1"
start "" /affinity 1 /wait cmd /c ""C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev\review-priordraws-warmcmp.R new > dev\review-priordraws-log\warmcmp-new.txt 2>&1"
