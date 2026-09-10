# Reviewer, item 6: does the hazard scanner REACH the new R/ file, and
# is roxygen in sync with the sources?
#
# The lane proves the first by planting `est$b` in its own worktree
# source and taking it out again. A reviewer must not edit the lane's
# files, so the same plant is made in a COPY, installed to a throwaway
# library, and test-bracket-access.R is run against that.
#
#   powershell -File dev/rev-latent-hazard.ps1

$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-latent"
$work = "C:\Users\adf44\source\r\rev-latent-haz"
$lib  = "C:\Users\adf44\source\r\rev-latent-hazlib"
$env:R_LIBS = "C:/Users/adf44/source/r/rev-latent-hazlib;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
Set-Location $tree

if (Test-Path $work) { Remove-Item -Recurse -Force $work }
if (Test-Path $lib)  { Remove-Item -Recurse -Force $lib }
New-Item -ItemType Directory -Force $work | Out-Null
New-Item -ItemType Directory -Force $lib  | Out-Null
Copy-Item -Recurse "$tree\extensions\frmtmb.latent\*" $work

Write-Output "======== hazard ABSENT (the shipped source) ========"
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL --library=$lib --no-multiarch $work 2>&1 | Select-String -Pattern "^\* DONE|ERROR"
Rscript "$tree\dev\rev-latent-runtests.R" "$work\tests\testthat\test-bracket-access.R" $lib 2>&1 |
  Select-String -Pattern "^ASSERT|^FAIL|^ERROR|^  \[|^      "
Write-Output "  RUN exit=$LASTEXITCODE"

Write-Output "======== hazard PLANTED in hmm-starts.R ========"
Rscript "$tree\dev\rev-latent-plant.R" $work
if ($LASTEXITCODE -ne 0) { Write-Output "PLANT FAILED"; exit 1 }
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL --library=$lib --no-multiarch $work 2>&1 | Select-String -Pattern "^\* DONE|ERROR"
Rscript "$tree\dev\rev-latent-runtests.R" "$work\tests\testthat\test-bracket-access.R" $lib 2>&1 |
  Select-String -Pattern "^ASSERT|^FAIL|^ERROR|^  \[|^      "
Write-Output "  RUN exit=$LASTEXITCODE  (non-zero means the scanner reaches hmm-starts.R)"

Remove-Item -Recurse -Force $work
Remove-Item -Recurse -Force $lib
Write-Output "======== done ========"
