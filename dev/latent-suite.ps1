# Lane `latent`: the whole frmtmb.latent suite, ONE TEST FILE PER R
# PROCESS, with counts that cannot be mistaken for a clean run.
#
# Two traps this avoids, both of which have bitten this project:
#   * `$?` is false the moment a native command writes to stderr, so
#     every gate here reads $LASTEXITCODE;
#   * PowerShell flattens nested arrays, so the file list is built with
#     an ArrayList and each name is read into a local before it reaches
#     the child process.
#
#   powershell -File dev/latent-suite.ps1
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-latent"
Set-Location $tree
$files = New-Object System.Collections.ArrayList
foreach ($f in Get-ChildItem "$tree\extensions\frmtmb.latent\tests\testthat" -Filter "test-*.R" | Sort-Object Name) {
  [void]$files.Add($f.FullName)
}
Write-Output "FILES TO RUN: $($files.Count)"
$log = "$tree\dev\latent-suite-log.txt"
if (Test-Path $log) { Remove-Item $log }
$bad = 0
foreach ($p in $files) {
  $nm = Split-Path $p -Leaf
  Write-Output "---- $nm"
  Rscript dev/latent-runtests.R $p 2>&1 | Tee-Object -FilePath $log -Append | Select-String -Pattern "^BLOCKS|^ASSERT|^PASS|^FAIL|^ERROR|^SKIP|^WARN|^SECONDS|^  \["
  if ($LASTEXITCODE -ne 0) { $bad = $bad + 1; Write-Output "  ** $nm EXITED $LASTEXITCODE **" }
}
Write-Output "FILES WITH FAILURES OR ERRORS: $bad of $($files.Count)"
