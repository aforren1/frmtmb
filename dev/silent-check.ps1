# Builds and checks one package of the worktree with --as-cran, inside
# dev/silent-check/. Usage: powershell -File dev/silent-check.ps1 <dir>
param([string]$pkgdir)
$root = "C:/Users/adf44/source/r/frmtmb-wt-silent"
$out = "$root/dev/silent-check"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;" +
  "C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH +
  ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-silent-lib;" +
  "C:/Users/adf44/source/r/rellib-r20;" +
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "$root/dev/stan-cache"
Set-Location $out
$R = "C:\Program Files\R\R-4.6.1\bin\R.exe"
& $R CMD build "$root/$pkgdir" 2>&1 | Out-File -Encoding utf8 "$out/build-$(Split-Path $pkgdir -Leaf).log"
$name = (Get-Content "$root/$pkgdir/DESCRIPTION" | Select-String "^Package:").ToString().Split(" ")[-1]
$ver = (Get-Content "$root/$pkgdir/DESCRIPTION" | Select-String "^Version:").ToString().Split(" ")[-1]
& $R CMD check --as-cran "$name`_$ver.tar.gz" 2>&1 | Out-File -Encoding utf8 "$out/check-$name.log"
Get-Content "$out/check-$name.log" | Select-String "Status:|WARNING|ERROR|NOTE"
