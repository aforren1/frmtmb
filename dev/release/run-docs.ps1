# Rebuild the eight pkgdown sites. Core goes to docs/, each extension
# to docs/<pkg>/ via the destination in its own _pkgdown.yml.
#
# Counts sites that actually built. A driver reporting success over
# zero work has happened three times in this project.

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$LOG  = "$ROOT/dev/release/docs.log"

if (-not (Test-Path $R)) { throw "Rscript missing: $R" }

$env:R_LIBS = "C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:FRMTMB_STAN_CACHE = "$ROOT/dev/stan-cache"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;" + $env:PATH

Remove-Item $LOG -ErrorAction SilentlyContinue

$pkgs = New-Object System.Collections.ArrayList
[void]$pkgs.Add($ROOT)
foreach ($n in @("frmtmb.eam", "frmtmb.latent", "frmtmb.ode",
                 "frmtmb.spline", "frmtmb.coupling", "frmtmb.sample",
                 "frmtmb.learn")) {
  [void]$pkgs.Add("$ROOT/extensions/$n")
}
foreach ($p in $pkgs) {
  if (-not (Test-Path "$p/_pkgdown.yml")) { throw "no _pkgdown.yml: $p" }
}

$ErrorActionPreference = "Continue"
$ok = 0
foreach ($p in $pkgs) {
  Add-Content -Path $LOG -Value ("===== SITE " + $p + " =====")
  $expr = "pkgdown::build_site(pkg = '" + $p + "', preview = FALSE, install = FALSE, new_process = FALSE, lazy = FALSE)"
  $out = & $R -e $expr 2>&1
  Add-Content -Path $LOG -Value $out
  if ($LASTEXITCODE -eq 0) {
    $ok = $ok + 1
    Add-Content -Path $LOG -Value "===== ok ====="
  } else {
    Add-Content -Path $LOG -Value ("===== FAILED exit " + $LASTEXITCODE + " =====")
  }
}
Add-Content -Path $LOG -Value ("DOCS built " + $ok + " of " + $pkgs.Count)
Write-Output ("DOCS built " + $ok + " of " + $pkgs.Count)
