# Lane parity2: R CMD build and R CMD check --as-cran of one package from
# the worktree, inside dev/parity2-check/.
#   powershell -File dev/parity2-check.ps1 <package dir relative to the worktree>
param([string]$Pkg = ".")
$ErrorActionPreference = "Continue"
$wt = "C:/Users/adf44/source/r/frmtmb-wt-parity2"
$out = "$wt/dev/parity2-check"
New-Item -ItemType Directory -Force $out | Out-Null
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-parity2-lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:NOT_CRAN = ""
$r = "C:\Program Files\R\R-4.6.1\bin\R.exe"
Set-Location $out
& $r CMD build "$wt/$Pkg" 2>&1 | Select-Object -Last 3
$name = (Get-Content "$wt/$Pkg/DESCRIPTION" | Select-String "^Package:").Line.Split(" ")[-1]
$tar = (Get-ChildItem "$out/${name}_*.tar.gz" | Sort-Object LastWriteTime | Select-Object -Last 1).FullName
& $r CMD check --as-cran $tar 2>&1 | Out-File -Encoding utf8 "$out/$name-check.log"
Get-Content "$out/$name-check.log" | Select-String "Status:|WARNING|ERROR|NOTE" | Select-Object -Last 15
