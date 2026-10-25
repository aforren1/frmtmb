@echo off
rem The smallfix review's bitwise sample (r1-bitwise.R, run from the
rem lane's tree, unedited) on rellib-r12 and rellib-r13, each arm twice,
rem every process pinned to CPUs 0-1 from its start (rtmb-pitfalls item
rem 21), as the review ran it. Called by dev/rel075-equality.sh.
set NOT_CRAN=true
set RS=C:\Program Files\R\R-4.6.1\bin\Rscript.exe
set S=C:\Users\adf44\source\r\frmtmb-wt-smallfix\dev\review-smallfix-out\scripts\r1-bitwise.R
set O=C:\Users\adf44\source\r\frmtmb-wt-rel075\dev\rel075-log\eq\smallfix
set BASE=C:/Users/adf44/source/r/rellib-r12
set MERGED=C:/Users/adf44/source/r/rellib-r13
for %%i in (1 2) do (
  start "" /affinity 3 /wait /b cmd /c ""%RS%" "%S%" "%BASE%" "%O%\base%%i.rds" > "%O%\base%%i.txt" 2>&1"
  start "" /affinity 3 /wait /b cmd /c ""%RS%" "%S%" "%MERGED%" "%O%\merged%%i.rds" > "%O%\merged%%i.txt" 2>&1"
)
