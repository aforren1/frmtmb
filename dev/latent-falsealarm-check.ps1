# Lane `latent`: see each guard test FAIL, on the code it pins.
#
# A test that has only been seen passing is not evidence. This rebuilds
# the package on a COPY with one shipped decision undone at a time, into
# a scratch library that is not this lane's, and runs the test file that
# pins it. The copies and the scratch library are thrown away; nothing
# here touches the worktree's own sources or `latent-lib`.
#
#   notol  round 1's fix removed: the two comparisons carry no
#          tolerance, and a unimodal fit is told it found a local
#          optimum.
#   refs   round 2's fix removed: one tolerance FUNCTION called with a
#          different reference at each site, so a band of `grad_tol`
#          makes the printout contradict itself.
#   w05    the shrink weight moved to the low edge of the usable
#          interval, where it loses seed 20270476 by 256 units.
#   w099   the shrink weight moved to the high edge, where it loses the
#          same seed from the other side.
#
#   powershell -File dev/latent-falsealarm-check.ps1
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-latent"
$work = "C:\Users\adf44\source\r\latent-absent"
$lib  = "C:\Users\adf44\source\r\latent-absent-lib"
$env:R_LIBS = "C:/Users/adf44/source/r/latent-absent-lib;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$modes = @{ notol = "test-hmm-starts.R"; refs = "test-hmm-starts.R";
            w05 = "test-lca.R"; w099 = "test-lca.R" }
Set-Location $tree
foreach ($mode in @("notol", "refs", "w05", "w099")) {
  Write-Output "======== mode $mode, $($modes[$mode]) ========"
  Remove-Item -Recurse -Force $work, $lib -ErrorAction SilentlyContinue
  New-Item -ItemType Directory -Force $work | Out-Null
  New-Item -ItemType Directory -Force $lib  | Out-Null
  Copy-Item -Recurse -Force "$tree\extensions\frmtmb.latent\*" $work
  Rscript dev/latent-falsealarm-revert.R $work $mode
  if ($LASTEXITCODE -ne 0) { Write-Output "REVERT FAILED"; continue }
  & "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL --library=$lib --no-multiarch $work 2>&1 | Select-String -Pattern "^\* DONE|ERROR"
  Write-Output "install exit=$LASTEXITCODE"
  Rscript dev/latent-runtests.R "$work\tests\testthat\$($modes[$mode])" $lib 2>&1 |
    Select-String -Pattern "^BLOCKS|^ASSERT|^PASS|^FAIL|^ERROR|^SKIP|^  \[|^      "
  Write-Output "RUN exit=$LASTEXITCODE (non-zero is the point)"
}
Remove-Item -Recurse -Force $work, $lib -ErrorAction SilentlyContinue
