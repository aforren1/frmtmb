# R CMD check --as-cran on frmtmb, ONCE, on the final pass (lane
# wt-famlink). Adapted from dev/argspell-check.ps1, which records the
# reason for each setting. The tarball and check directory live outside
# the worktree so neither is picked up as package source.
$ErrorActionPreference = "Continue"
$root = Split-Path -Parent $PSScriptRoot
$pandoc = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:RSTUDIO_PANDOC = $pandoc
$env:PATH = "$pandoc;" +
  "C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/famlink-lib;" +
  "C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;" +
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:NOT_CRAN = ""
$env:FRMTMB_STAN_CACHE = Join-Path $root "dev/stan-cache"
$work = "C:\Users\adf44\source\r\famlink-check"
if (-not (Test-Path $work)) { New-Item -ItemType Directory $work | Out-Null }
Get-ChildItem -Path $work -Filter "frmtmb_*.tar.gz" |
  Remove-Item -Confirm:$false
Set-Location $work
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build $root
$tar = Get-ChildItem -Path $work -Filter "frmtmb_*.tar.gz" |
  Sort-Object LastWriteTime | Select-Object -Last 1
Write-Output "built $($tar.Name)"
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran `
  --output=$work $tar.FullName
