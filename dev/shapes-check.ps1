# R CMD check --as-cran on the two packages this lane changed the API
# of, built WITH vignettes, against the lane's private library.
# Modeled on dev/release/run-check.ps1, and it keeps that file's two
# rules: no --no-build-vignettes (a tarball with vignette sources and no
# inst/doc manufactures two WARNINGs and a NOTE) and no --no-manual (it
# skips the HTML-manual check, which is where this machine's V8
# math-rendering NOTE and an unescaped % in Rd would show).
#
#   powershell -File dev/shapes-check.ps1 [-Only core|sample|spline]

param([string]$Only = "")

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-shapes"
$R    = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$LOG  = "$ROOT/dev/shapes-log/check.log"
$OUT  = "C:/Users/adf44/source/r/shapes-check"

if (-not (Test-Path $R)) { throw "R missing: $R" }

$env:R_LIBS = "C:/Users/adf44/source/r/shapes-lib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "$ROOT/dev/stan-cache"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH

if (Test-Path $OUT) { Remove-Item -Recurse -Force $OUT }
New-Item -ItemType Directory -Force $OUT | Out-Null
Remove-Item $LOG -ErrorAction SilentlyContinue

$pkgs = New-Object System.Collections.ArrayList
if ($Only -eq "" -or $Only -eq "core") { [void]$pkgs.Add($ROOT) }
if ($Only -eq "" -or $Only -eq "sample") {
  [void]$pkgs.Add("$ROOT/extensions/frmtmb.sample")
}
# punch round 3 changed ?royston_parmar, so frmtmb.spline is checked too
if ($Only -eq "" -or $Only -eq "spline") {
  [void]$pkgs.Add("$ROOT/extensions/frmtmb.spline")
}
foreach ($p in $pkgs) {
  if (-not (Test-Path "$p/DESCRIPTION")) { throw "not a package: $p" }
}

$ErrorActionPreference = "Continue"
$done = 0
foreach ($p in $pkgs) {
  Add-Content -Path $LOG -Value ("===== CHECK " + $p + " =====")
  Push-Location $OUT
  & $R CMD build $p 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
  $tgz = Get-ChildItem -Path $OUT -Filter "*.tar.gz" | Sort-Object LastWriteTime | Select-Object -Last 1
  if ($tgz) {
    & $R CMD check --as-cran $tgz.FullName 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
    $done = $done + 1
  } else {
    Add-Content -Path $LOG -Value "NO TARBALL BUILT"
  }
  Pop-Location
}
Add-Content -Path $LOG -Value ("CHECKED " + $done + " of " + $pkgs.Count)
Write-Output ("CHECKED " + $done + " of " + $pkgs.Count)
