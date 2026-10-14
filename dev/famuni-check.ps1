# R CMD build and check --as-cran of one package of the worktree, inside
# dev/famuni-check/<pkg>/. Usage:
#   powershell -File dev/famuni-check.ps1 <path relative to the worktree>
param([string]$Rel)
$WT = "C:\Users\adf44\source\r\frmtmb-wt-famuni"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-famuni-lib;C:/Users/adf44/source/r/rellib-r8;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:NOT_CRAN = "true"
$name = Split-Path $Rel -Leaf
if ($name -eq "frmtmb-wt-famuni" -or $Rel -eq ".") { $name = "frmtmb" }
$dir = Join-Path $WT "dev\famuni-check\$name"
New-Item -ItemType Directory -Force $dir | Out-Null
Set-Location $dir
$src = Join-Path $WT $Rel
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build $src *> build.log
$tgz = Get-ChildItem -Filter "$name*.tar.gz" | Sort-Object LastWriteTime | Select-Object -Last 1
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran $tgz.Name *> check.log
Select-String -Path check.log -Pattern "^Status" | ForEach-Object { $_.Line }
