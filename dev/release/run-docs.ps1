# A thin wrapper over dev/release/build-docs.R for the Windows
# development box. It sets the paths and libraries this machine needs
# and does nothing else: the build, the count guard and the layout live
# in the R script, so this file and the GitHub Actions job run the same
# code.
#
# The published site is built by .github/workflows/pkgdown.yaml, not
# here. Use this to see a change before pushing it.
#
# By default it builds into a SCRATCH tree, dev/docsci-site, and not
# into docs/. docs/ is on its way out of the repository; pass
# -Dest <path> to put the site somewhere else.
#
#   run-docs.ps1                              all eight sites
#   run-docs.ps1 -Pkgs frmtmb.eam -MinSites 1  one site
#
# -MinSites is not lowered for you when -Pkgs is given. A partial tree
# has dead navbar links, so the count guard is meant to fail on it
# unless you say you know.

param(
  [string]$Dest = "",
  [string]$Pkgs = "",
  [int]$MinSites = 8
)

$ErrorActionPreference = "Stop"
$ROOT = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"

if (-not (Test-Path $R)) { throw "Rscript missing: $R" }
if ($Dest -eq "") { $Dest = "$ROOT/dev/docsci-site" }

$env:R_LIBS = "C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
# StanHeaders 2.39.1 compiles only with the user Makevars C++17 flag, and
# HOME depends on the launcher, so name the file (dev/tmbstan121-findings.md)
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;" + $env:PATH

# An ArrayList rather than an array literal: PowerShell flattens nested
# arrays, which has silently dropped arguments in this project before.
$argv = New-Object System.Collections.ArrayList
[void]$argv.Add("$ROOT/dev/release/build-docs.R")
[void]$argv.Add("--dest=$Dest")
[void]$argv.Add("--log=$ROOT/dev/release/docs.log")
[void]$argv.Add("--min-sites=$MinSites")
if ($Pkgs -ne "") { [void]$argv.Add("--pkgs=$Pkgs") }

$ErrorActionPreference = "Continue"
& $R $argv.ToArray()
$code = $LASTEXITCODE
Write-Output ("DOCS exit " + $code)
exit $code
