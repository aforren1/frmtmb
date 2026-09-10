# R CMD check --as-cran on all eight, built WITH vignettes.
#
# Do not pass --no-build-vignettes to R CMD build: a tarball with
# vignette sources and no inst/doc manufactures two WARNINGs and a NOTE
# on every package, which at 0.55.2 looked like a finding on four
# packages the round had not touched.
#
# Do not pass --no-manual either. It skips the HTML-manual check, which
# is where this machine's environmental V8 math-rendering NOTE comes
# from. Skipping it makes the release look cleaner than the lanes'
# own checks and is not comparable to them.

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb"
$R    = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$LOG  = "$ROOT/dev/release/check.log"
$OUT  = "C:/Users/adf44/source/r/checkout-r3"

if (-not (Test-Path $R)) { throw "R missing: $R" }

$env:R_LIBS = "C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH

if (Test-Path $OUT) { Remove-Item -Recurse -Force $OUT }
New-Item -ItemType Directory -Force $OUT | Out-Null
Remove-Item $LOG -ErrorAction SilentlyContinue

$pkgs = New-Object System.Collections.ArrayList
[void]$pkgs.Add($ROOT)
foreach ($n in @("frmtmb.eam", "frmtmb.latent", "frmtmb.ode",
                 "frmtmb.spline", "frmtmb.coupling", "frmtmb.sample",
                 "frmtmb.learn")) {
  [void]$pkgs.Add("$ROOT/extensions/$n")
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
