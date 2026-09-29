# Lane splinecurve, round 2: R CMD build, then R CMD check --as-cran ONCE,
# on the final frmtmb.spline source, inside dev/splinecurve-check/final-r2.
# The tarball is named by the DESCRIPTION version, not globbed, so a
# tarball left by an earlier round cannot be the one checked.
. C:\Users\adf44\source\r\frmtmb-wt-release\dev\splinecurve-env.ps1
$ErrorActionPreference = "Continue"
$root = "C:\Users\adf44\source\r\frmtmb-wt-release"
$pkg = "$root\extensions\frmtmb.spline"
$dir = "$root\dev\splinecurve-check\final-r2"
New-Item -ItemType Directory -Force $dir | Out-Null
Set-Location $dir
$ver = ((Select-String -Path "$pkg\DESCRIPTION" -Pattern "^Version:").Line -replace "^Version:\s*", "").Trim()
$tgz = "$dir\frmtmb.spline_$ver.tar.gz"
if (Test-Path $tgz) { Remove-Item $tgz -Confirm:$false }
$R = "C:\Program Files\R\R-4.6.1\bin\R.exe"
& $R CMD build $pkg *> "$root\dev\splinecurve-build-r2.log"
"build exit $LASTEXITCODE, version $ver"
if (-not (Test-Path $tgz)) { "NO TARBALL $tgz"; exit 1 }
& $R CMD check --as-cran $tgz *> "$root\dev\splinecurve-check-r2.log"
"check exit $LASTEXITCODE"
Select-String -Path "$root\dev\splinecurve-check-r2.log" -Pattern "^Status"
