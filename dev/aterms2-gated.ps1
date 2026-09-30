# Run one test file, gated tier, one process. Usage:
#   powershell -File dev/aterms2-gated.ps1 <pkg> <dir> <file> <arm> <log>
# <arm> is "lane" (this lane's library first) or "base" (rellib-r3).
param([string]$Pkg, [string]$Dir, [string]$File, [string]$Arm,
      [string]$Log)
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:TMP = "C:\Users\adf44\AppData\Local\Temp\1"
$env:TEMP = "C:\Users\adf44\AppData\Local\Temp\1"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/stan-cache"
$env:NOT_CRAN = "true"
$env:FRMTMB_BRMS_FIT_TESTS = "true"
$ErrorActionPreference = "Continue"
Set-Location $Dir
$runner = "C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/aterms2-testfile.R"
$a = @($runner, $Pkg, $File)
if ($Arm -eq "base") { $a += "base" }
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" @a 2>&1 |
  ForEach-Object { "$_" } |
  Out-File -FilePath $Log -Encoding utf8
