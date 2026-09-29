# Lane splinecurve: R CMD build, then R CMD check --as-cran ONCE, on the
# final frmtmb.spline source, inside dev/splinecurve-check/final.
. C:\Users\adf44\source\r\frmtmb-wt-release\dev\splinecurve-env.ps1
$ErrorActionPreference = "Continue"
$root = "C:\Users\adf44\source\r\frmtmb-wt-release"
$dir = "$root\dev\splinecurve-check\final"
New-Item -ItemType Directory -Force $dir | Out-Null
Set-Location $dir
$R = "C:\Program Files\R\R-4.6.1\bin\R.exe"
& $R CMD build "$root\extensions\frmtmb.spline" *> "$root\dev\splinecurve-build-after.log"
"build exit $LASTEXITCODE"
$tgz = Get-ChildItem "$dir\frmtmb.spline_*.tar.gz" | Select-Object -First 1
if ($null -eq $tgz) { "NO TARBALL"; exit 1 }
& $R CMD check --as-cran $tgz.FullName *> "$root\dev\splinecurve-check-after.log"
"check exit $LASTEXITCODE"
Select-String -Path "$root\dev\splinecurve-check-after.log" -Pattern "^Status"
