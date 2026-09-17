# The runs the ungated suite cannot make:
#  - the control arm C (tmbstan 1.2.0 built against StanHeaders 2.39.1)
#    through the suite's own guard file and the file that broke CI;
#  - the brms-fit tier, whose two blocks compile FRESH brms programs,
#    in A, in B, and in B with the C++17 override.
# One R process per file; jobs run in parallel, each with its own log.
$WT  = "C:/Users/adf44/source/r/frmtmb-wt-tmbstan121"
$LIB = "C:/Users/adf44/source/r/tmbstan121-lib"
$R   = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S   = "$WT/dev/tmbstan121-05-run-tests.R"
$DIR = "$WT/extensions/frmtmb.sample/tests/testthat"
$env:NOT_CRAN = "true"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

$jobs = New-Object System.Collections.ArrayList
[void]$jobs.Add(@{ tag = "C-guard"; arm = "C"; f = "test-tmbstan-build-guard.R"; fit = ""; mv = "" })
[void]$jobs.Add(@{ tag = "C-ce"; arm = "C"; f = "test-conditional-effects-draws.R"; fit = ""; mv = "" })
foreach ($f in @("test-loo.R", "test-sampling-ported.R")) {
  [void]$jobs.Add(@{ tag = "A-fit-$f"; arm = "A"; f = $f; fit = "true"; mv = "" })
  [void]$jobs.Add(@{ tag = "B-fit-$f"; arm = "B"; f = $f; fit = "true"; mv = "" })
  [void]$jobs.Add(@{ tag = "Bcxx17-fit-$f"; arm = "B"; f = $f; fit = "true"; mv = "$LIB/makevars/Makevars.win" })
}
$procs = New-Object System.Collections.ArrayList
foreach ($j in $jobs) {
  $env:FRMTMB_BRMS_FIT_TESTS = $j.fit
  if ($j.mv) { $env:R_MAKEVARS_USER = $j.mv } else { Remove-Item Env:R_MAKEVARS_USER -ErrorAction SilentlyContinue }
  $env:FRMTMB_STAN_CACHE = "$LIB/stan-cache-" + $j.arm
  $log = "$WT/dev/tmbstan121-extra-" + $j.tag + ".log"
  # Built first: inside @(), the comma binds before +, which once
  # passed the directory instead of the file.
  $tf = "$DIR/" + $j.f
  $p = Start-Process -FilePath $R -ArgumentList @($S, $j.arm, "frmtmb.sample", $tf) -RedirectStandardOutput $log -RedirectStandardError ($log + ".err") -NoNewWindow -PassThru
  [void]$procs.Add($p)
}
foreach ($p in $procs) { $p.WaitForExit() }
foreach ($j in $jobs) {
  $log = "$WT/dev/tmbstan121-extra-" + $j.tag + ".log"
  Write-Output ("== " + $j.tag)
  Get-Content $log | Select-String -Pattern "^(RESULT|DETAIL)" | ForEach-Object { $_.Line }
}
