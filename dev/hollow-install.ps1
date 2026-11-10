# Lane hollow: roxygenise and install core and frmtmb.sample into the
# lane library. Usage: powershell -File dev/hollow-install.ps1 [core|sample|all]
param([string]$what = "all")
$ErrorActionPreference = "Continue"
$wt = "C:/Users/adf44/source/r/frmtmb-wt-hollow"
$lib = "C:/Users/adf44/source/r/wt-hollow-lib"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$r = "C:\Program Files\R\R-4.6.1\bin\R.exe"
$env:PATH = $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$env:R_LIBS = "$lib;C:/Users/adf44/source/r/rellib-r19;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$pkgs = @()
if ($what -eq "core" -or $what -eq "all") { $pkgs += $wt }
if ($what -eq "sample" -or $what -eq "all") { $pkgs += "$wt/extensions/frmtmb.sample" }
foreach ($p in $pkgs) {
  & $rs -e "roxygen2::roxygenise('$p')"
  Write-Output "ROXYGENISE EXIT $LASTEXITCODE $p"
  & $r CMD INSTALL "--library=$lib" --no-multiarch $p
  Write-Output "INSTALL EXIT $LASTEXITCODE $p"
}
