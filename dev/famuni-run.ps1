# Run one R script with the toolchain on PATH, stdout and stderr to
# dev/famuni-out/<name>.txt and .err. Usage:
#   powershell -File dev/famuni-run.ps1 <script.R> <name> [args...]
param([string]$Script, [string]$Name)
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:NOT_CRAN = "true"
$env:FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/stan-cache"
Set-Location C:\Users\adf44\source\r\frmtmb-wt-famuni
$argl = @($Script) + $args
$p = Start-Process -FilePath "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" `
  -ArgumentList $argl -RedirectStandardOutput "dev/famuni-out/$Name.txt" `
  -RedirectStandardError "dev/famuni-out/$Name.err" -NoNewWindow -Wait -PassThru
Write-Output "exit $($p.ExitCode)"
