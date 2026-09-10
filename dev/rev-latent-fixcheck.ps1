# Reviewer: does the proposed fix for BLOCKER 1 actually separate the
# two cases, and does it leave the real detection alone?
#
# Mutant F1 gives the two untoleranced comparisons in hmm_starts() the
# same relative tolerance hmm_starts_modes() already uses. A fix that
# has only been reasoned about is not a fix, so it is built, installed
# and run against both ends: the six unimodal fits that produced the
# false alarm, and the d4 probe the function exists for. The whole
# test file runs too, since a fix that turns the suite red is not a fix
# either.
#
#   powershell -File dev/rev-latent-fixcheck.ps1

$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-latent"
$work = "C:\Users\adf44\source\r\rev-latent-fix"
$lib  = "C:\Users\adf44\source\r\rev-latent-fixlib"
$env:R_LIBS = "C:/Users/adf44/source/r/rev-latent-fixlib;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
Set-Location $tree

New-Item -ItemType Directory -Force $work | Out-Null
New-Item -ItemType Directory -Force $lib  | Out-Null
Copy-Item -Recurse -Force "$tree\extensions\frmtmb.latent\*" $work
Rscript "$tree\dev\rev-latent-mutate.R" $work F1
if ($LASTEXITCODE -ne 0) { Write-Output "FIX PATCH FAILED"; exit 1 }
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL --library=$lib --no-multiarch $work 2>&1 | Select-String -Pattern "^\* DONE|ERROR"
Write-Output "install exit=$LASTEXITCODE"

Write-Output "======== the six unimodal fits, with the fix ========"
Rscript "$tree\dev\rev-latent-fixcheck.R" $lib 2>&1

Write-Output "======== test-hmm-starts.R, with the fix ========"
Rscript "$tree\dev\rev-latent-runtests.R" "$work\tests\testthat\test-hmm-starts.R" $lib 2>&1 |
  Select-String -Pattern "^BLOCKS|^ASSERT|^FAIL|^ERROR|^SKIP|^  \[|^      "
Write-Output "  RUN exit=$LASTEXITCODE"
Write-Output "======== done ========"
