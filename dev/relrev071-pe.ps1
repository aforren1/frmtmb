# Release reviewer of 0.71.0: the changed and the SE/optimizer test files
# on rellib-r9 under ref pinned to the P-cores (C03C03), to the E-cores
# (3FC), and OpenBLAS 0.3.26 LAPACK at 4 threads, one file per process
# through dev/ciharden-run1.R, every gate on, 8 at once.
#   powershell -File dev/relrev071-pe.ps1 <jobs file> <tag>
param([string]$jobsFile, [string]$tag)
$wt = "C:\Users\adf44\source\r\frmtmb-wt-release"
$log = "$wt\dev\relrev071-log\pe-$tag"
$ref = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$ob26 = "$wt\dev\relrev071-out\Rob0.3.26-lapack\bin\x64\Rscript.exe"
$run1 = "$wt\dev\ciharden-run1.R"
$env:NOT_CRAN = "true"
$env:FRMTMB_BRMS_FIT_TESTS = "true"
$env:FRMTMB_DRMTMB_FIT_TESTS = "true"
$env:FRMTMB_FUZZ = "true"
$env:FRMTMB_SCALE_TESTS = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "$wt\dev\relrev071-out\stan-cache"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
New-Item -ItemType Directory -Force "$log\refP", "$log\refE", "$log\ob26" | Out-Null
$cfgs = @(@("refP", $ref, "C03C03"), @("refE", $ref, "3FC"), @("ob26", $ob26, ""))
$lines = Get-Content $jobsFile | Where-Object { $_ -ne "" }
$procs = New-Object System.Collections.ArrayList
foreach ($c in $cfgs) {
  foreach ($ln in $lines) {
    $p, $f = $ln -split " ", 2
    $b = [System.IO.Path]::GetFileNameWithoutExtension($f)
    $out = "$log\$($c[0])\$p--$b.txt"
    $ff = "$wt\$f"
    if ($c[0] -eq "ob26") {
      $cmd = "set OPENBLAS_NUM_THREADS=4&& set OMP_NUM_THREADS=4&& `"$($c[1])`" `"$run1`" base $p `"$ff`" > `"$out`" 2>&1"
    } else {
      $cmd = "start `"`" /affinity $($c[2]) /b /wait `"$($c[1])`" `"$run1`" base $p `"$ff`" > `"$out`" 2>&1"
    }
    while (@($procs | Where-Object { -not $_.HasExited }).Count -ge 8) {
      Start-Sleep -Milliseconds 500
    }
    [void]$procs.Add((Start-Process -FilePath cmd.exe -ArgumentList "/c", "`"$cmd`"" -WindowStyle Hidden -PassThru))
  }
}
$procs | Wait-Process
"done" | Out-File "$log\done.txt"
