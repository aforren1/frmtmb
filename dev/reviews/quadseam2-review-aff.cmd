@echo off
rem Rscript pinned to the cores of mask %1 (rtmb-pitfalls.md item 21), rest
rem of the arguments passed on. Review of lane quadseam2.
set MASK=%1
shift
start "" /b /wait /affinity %MASK% "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" %1 %2 %3 %4 %5 %6
