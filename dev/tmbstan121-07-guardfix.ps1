# The one test change this lane made, in every arm it has to hold in.
#   C, A, B  : the edited file in the worktree
#   A-mut    : edited file, fixture that does NOT restore state (must FAIL)
#   A-mutorig: original file, same mutant (the old assertion, for scale)
$WT  = "C:/Users/adf44/source/r/frmtmb-wt-tmbstan121"
$LIB = "C:/Users/adf44/source/r/tmbstan121-lib"
$R   = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S   = "$WT/dev/tmbstan121-05-run-tests.R"
$F   = "test-tmbstan-build-guard.R"
$env:NOT_CRAN = "true"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$jobs = New-Object System.Collections.ArrayList
[void]$jobs.Add(@{ tag = "C"; arm = "C"; f = "$WT/extensions/frmtmb.sample/tests/testthat/$F" })
[void]$jobs.Add(@{ tag = "A"; arm = "A"; f = "$WT/extensions/frmtmb.sample/tests/testthat/$F" })
[void]$jobs.Add(@{ tag = "B"; arm = "B"; f = "$WT/extensions/frmtmb.sample/tests/testthat/$F" })
[void]$jobs.Add(@{ tag = "A-mut"; arm = "A"; f = "$LIB/mut/$F" })
[void]$jobs.Add(@{ tag = "A-mutorig"; arm = "A"; f = "$LIB/mut-orig/$F" })
$procs = New-Object System.Collections.ArrayList
foreach ($j in $jobs) {
  $env:FRMTMB_STAN_CACHE = "$LIB/stan-cache-" + $j.arm
  $log = "$WT/dev/tmbstan121-guardfix-" + $j.tag + ".log"
  $tf = $j.f
  $p = Start-Process -FilePath $R -ArgumentList @($S, $j.arm, "frmtmb.sample", $tf) -RedirectStandardOutput $log -RedirectStandardError ($log + ".err") -NoNewWindow -PassThru
  [void]$procs.Add($p)
}
foreach ($p in $procs) { $p.WaitForExit() }
foreach ($j in $jobs) {
  Write-Output ("== " + $j.tag)
  Get-Content ("$WT/dev/tmbstan121-guardfix-" + $j.tag + ".log") | Select-String -Pattern "^(RESULT|DETAIL)" | ForEach-Object { $_.Line }
}
