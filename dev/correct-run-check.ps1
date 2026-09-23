# Lane wt-correct copy of dev/release/run-check.ps1: the lane tree, the
# lane library in front, and the two packages the lane touched. The
# release script itself is not edited.
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
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-correct"
$R    = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$LOG  = "$ROOT/dev/correct-log/check.log"
# CORRECT_LIB and CORRECT_CHECK let a second worker on this lane check
# against its own private library and its own scratch directory.
$LIB  = if ($env:CORRECT_LIB) { $env:CORRECT_LIB } else { "C:/Users/adf44/source/r/correct-lib" }
$OUT  = if ($env:CORRECT_CHECK) { $env:CORRECT_CHECK } else { "C:/Users/adf44/source/r/correct-check" }

if (-not (Test-Path $R)) { throw "R missing: $R" }

$env:R_LIBS = "$LIB;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
# StanHeaders 2.39.1 compiles only with the user Makevars C++17 flag, and
# HOME depends on the launcher, so name the file (dev/tmbstan121-findings.md)
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH

if (Test-Path $OUT) { Remove-Item -Recurse -Force $OUT }
New-Item -ItemType Directory -Force $OUT | Out-Null
Remove-Item $LOG -ErrorAction SilentlyContinue

$pkgs = New-Object System.Collections.ArrayList
[void]$pkgs.Add($ROOT)
# the lane touched frmtmb and frmtmb.sample only
foreach ($n in @("frmtmb.sample")) {
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
    # A passing "checking tests ... OK" does NOT echo the testthat
    # summary into the check log: the FAIL/PASS/SKIP counts live only in
    # the .Rcheck tree, and the low-disk rule on this machine says to
    # delete that tree as soon as the Status line is read. A record that
    # quotes in-check counts and cites the check log is therefore citing
    # a file that does not hold them, which is how punch round 1 got one
    # wrong. The file is a few KB, so it is kept beside the log instead.
    $nm = ($tgz.Name -replace "_.*$", "")
    $dest = Join-Path $ROOT "dev/correct-log/$nm-testthat.Rout"
    # CLEAR FIRST. A -Force copy onto a fixed name leaves the PREVIOUS
    # run's file behind when this run produces none, and a reader who
    # quotes the file rather than the log line then gets last run's
    # counts. That is the same failure A2 was, so the destination goes
    # before the source is looked for.
    Remove-Item $dest -Force -ErrorAction SilentlyContinue
    $rout = Join-Path $OUT "$nm.Rcheck/tests/testthat.Rout"
    if (-not (Test-Path $rout)) { $rout = "$rout.fail" }
    if (Test-Path $rout) {
      Copy-Item $rout $dest -Force
    } else {
      Add-Content -Path $LOG -Value "NO testthat.Rout FOR $nm"
    }
  } else {
    Add-Content -Path $LOG -Value "NO TARBALL BUILT"
  }
  Pop-Location
}
Add-Content -Path $LOG -Value ("CHECKED " + $done + " of " + $pkgs.Count)
Write-Output ("CHECKED " + $done + " of " + $pkgs.Count)
