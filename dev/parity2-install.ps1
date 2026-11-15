# Lane parity2: roxygenise and install core (and, with -Ext, the named
# extensions) from the worktree into the lane library.
#   powershell -File dev/parity2-install.ps1 [-Ext frmtmb.sample,...]
param([string[]]$Ext = @(), [switch]$NoCore)
$ErrorActionPreference = "Continue"
$wt = "C:/Users/adf44/source/r/frmtmb-wt-parity2"
$lib = "C:/Users/adf44/source/r/wt-parity2-lib"
$env:PATH = $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$env:R_LIBS = "$lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$r = "C:\Program Files\R\R-4.6.1\bin\R.exe"
$dirs = New-Object System.Collections.ArrayList
if (-not $NoCore) { [void]$dirs.Add($wt) }
foreach ($e in $Ext) { [void]$dirs.Add("$wt/extensions/$e") }
foreach ($d in $dirs) {
  & $rs -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$d')" 2>&1 | Select-String -Pattern "Error|Warning|Writing" | Select-Object -First 20
  & $r CMD INSTALL --library=$lib --no-multiarch $d 2>&1 | Select-String -Pattern "ERROR|DONE|Error" | Select-Object -First 10
}
