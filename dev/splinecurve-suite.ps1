# Lane splinecurve: the whole frmtmb.spline suite, one test file per R
# process, all at once, against the lane's private library. Each file
# writes its own log; the summary is read from the RESULT lines, and a
# file with no RESULT line is reported as MISSING rather than skipped.
param([string]$Lib = "C:/Users/adf44/source/r/wt-splinecurve-lib",
      [string]$Tag = "after")
. C:\Users\adf44\source\r\frmtmb-wt-release\dev\splinecurve-env.ps1 -Lib $Lib
$ErrorActionPreference = "Continue"
$root = "C:\Users\adf44\source\r\frmtmb-wt-release"
$runner = "$root\dev\splinecurve-run-tests.R"
$logdir = "$root\dev\splinecurve-suite-$Tag"
New-Item -ItemType Directory -Force $logdir | Out-Null
$files = Get-ChildItem "$root\extensions\frmtmb.spline\tests\testthat\test-*.R"
$procs = New-Object System.Collections.ArrayList
foreach ($f in $files) {
  $nm = $f.Name
  $log = "$logdir\$nm.log"
  $p = Start-Process -FilePath "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" `
    -ArgumentList @($runner, $Lib, "frmtmb.spline", $f.FullName) `
    -RedirectStandardOutput $log -RedirectStandardError "$log.err" `
    -NoNewWindow -PassThru
  [void]$procs.Add($p)
}
foreach ($p in $procs) { $p.WaitForExit() }
$ran = 0
foreach ($f in $files) {
  $nm = $f.Name
  $line = Select-String -Path "$logdir\$nm.log" -Pattern "^RESULT " |
    Select-Object -Last 1
  if ($null -eq $line) { "MISSING $nm" } else { $ran++; $line.Line }
  Select-String -Path "$logdir\$nm.log" -Pattern "^FAILED " |
    ForEach-Object { "  " + $_.Line }
}
"SUITE ran $ran of $($files.Count)"
