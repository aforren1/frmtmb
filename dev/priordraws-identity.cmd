@echo off
rem Six runs of dev/priordraws-identity.R, each pinned to logical core 0
cd /d C:\Users\adf44\source\r\frmtmb-wt-priordraws
for %%A in (base lane) do (
  for %%R in (1 2 3) do (
    start "" /wait /affinity 1 "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev/priordraws-identity.R %%A %%R
  )
)
