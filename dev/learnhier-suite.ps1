# Lane `learnhier`: the whole frmtmb.learn suite, ONE TEST FILE PER R
# PROCESS, with counts that cannot be mistaken for a clean run.
#
#   powershell -File dev/learnhier-suite.ps1            # ungated
#   powershell -File dev/learnhier-suite.ps1 -gates stan
#
# Two traps this avoids, both of which have bitten this project: `$?` is
# false the moment a native command writes to stderr, so every gate here
# reads `$LASTEXITCODE`; and PowerShell flattens nested arrays, so the
# file list is an ArrayList and each name is read into a local before it
# reaches the child process.
#
# dev/learnhier-testfile.R exits non-zero on a FAILURE OR AN ERROR, not
# on failures alone: a runner that sums `failed` and not `error` prints
# a clean line for a file that aborted halfway.
param([string]$gates = "")
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-learnhier"
Set-Location $tree
$files = New-Object System.Collections.ArrayList
$dir = "$tree\extensions\frmtmb.learn\tests\testthat"
foreach ($f in Get-ChildItem $dir -Filter "test-*.R" | Sort-Object Name) {
  [void]$files.Add($f.Name)
}
Write-Output "FILES TO RUN: $($files.Count)  gates: '$gates'"
$log = "$tree\dev\learnhier-suite-log.txt"
if (Test-Path $log) { Remove-Item $log }
$bad = 0
foreach ($nm in $files) {
  Write-Output "---- $nm"
  if ($gates -eq "") {
    Rscript dev/learnhier-testfile.R $nm 2>&1 |
      Tee-Object -FilePath $log -Append |
      Select-String -Pattern "^FILE|^PASS|^Error|Failure"
  } else {
    Rscript dev/learnhier-testfile.R $nm $gates 2>&1 |
      Tee-Object -FilePath $log -Append |
      Select-String -Pattern "^FILE|^PASS|^Error|Failure"
  }
  if ($LASTEXITCODE -ne 0) {
    $bad = $bad + 1
    Write-Output "  ** $nm EXITED $LASTEXITCODE **"
  }
}
Write-Output "FILES WITH FAILURES OR ERRORS: $bad of $($files.Count)"
