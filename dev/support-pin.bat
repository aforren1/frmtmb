@echo off
rem Lane support: run one R script pinned to core 0 (a P-core), so two
rem builds compare on one core type (dev/rtmb-pitfalls.md, item 21).
rem   dev\support-pin.bat <script> <log>   with FRMTMB_PORT_LIB set
start "" /affinity 1 /wait /b "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" %1 > %2 2>&1
