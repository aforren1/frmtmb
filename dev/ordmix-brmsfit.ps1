# Run dev/ordmix-brmsfit.R for every case at once, one R process each.
# Logs: dev/ordmix-brmsfit-log/<case>.txt
$wt = "C:\Users\adf44\source\r\frmtmb-wt-ordmix"
$out = "$wt\dev\ordmix-brmsfit-log"
New-Item -ItemType Directory -Force $out | Out-Null
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$procs = New-Object System.Collections.ArrayList
foreach ($c in @("none", "mu", "hgr", "hcs")) {
  $p = Start-Process -FilePath "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" `
    -ArgumentList @("$wt\dev\ordmix-brmsfit.R", $c) `
    -RedirectStandardOutput "$out\$c.txt" -RedirectStandardError "$out\$c.err" `
    -NoNewWindow -PassThru
  [void]$procs.Add($p)
}
foreach ($p in $procs) { $p.WaitForExit() }
foreach ($c in @("none", "mu", "hgr", "hcs")) {
  Get-Content "$out\$c.txt" | Select-String -Pattern "^(ROWS|SUMMARY) "
  Get-Content "$out\$c.err" | Select-String -Pattern "^Error" | ForEach-Object { "ERR $c : $_" }
}
