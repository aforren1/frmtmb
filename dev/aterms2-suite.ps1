# Run test files one per R process, $P at a time, against this lane's
# library (lane) or the base build (base). One log per file in
# dev/aterms2-suite-<tag>/, each ending in a RESULT line.
# Usage: powershell -File dev/aterms2-suite.ps1 <tag> <arm> <gated> <list> [P]
# <list> holds lines "<pkg> <dir> <file>"; <gated> is 0 or 1.
param([string]$Tag, [string]$Arm, [string]$Gated, [string]$List,
      [int]$P = 12)
$WT = "C:\Users\adf44\source\r\frmtmb-wt-aterms2"
$OUT = "$WT\dev\aterms2-suite-$Tag"
New-Item -ItemType Directory -Force $OUT | Out-Null
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:TMP = "C:\Users\adf44\AppData\Local\Temp\1"
$env:TEMP = $env:TMP
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/stan-cache"
$env:NOT_CRAN = "true"
if ($Gated -eq "1") { $env:FRMTMB_BRMS_FIT_TESTS = "true" }
$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$runner = "$WT\dev\aterms2-testfile.R"
$lines = Get-Content $List | Where-Object { $_.Trim() -ne "" }
$running = New-Object System.Collections.ArrayList
foreach ($ln in $lines) {
  $parts = $ln.Split(" ")
  $pkg = $parts[0]; $dir = $parts[1]; $f = $parts[2]
  while ($running.Count -ge $P) {
    Start-Sleep -Milliseconds 500
    foreach ($pr in @($running)) { if ($pr.HasExited) { [void]$running.Remove($pr) } }
  }
  $log = "$OUT\" + $pkg + "__" + $f.Replace(".R", "") + ".txt"
  $argl = @("`"$runner`"", $pkg, $f)
  if ($Arm -eq "base") { $argl += "base" }
  $pr = Start-Process -FilePath $R -ArgumentList $argl -WorkingDirectory $dir `
    -RedirectStandardOutput $log -RedirectStandardError "$log.err" `
    -NoNewWindow -PassThru
  [void]$running.Add($pr)
}
foreach ($pr in @($running)) { $pr.WaitForExit() }
$got = @(Get-ChildItem "$OUT\*.txt" | Where-Object {
  Select-String -Path $_.FullName -Pattern "^RESULT" -Quiet }).Count
"SUITE $Tag ran $got of $($lines.Count)"
Get-ChildItem "$OUT\*.txt" | ForEach-Object {
  Select-String -Path $_.FullName -Pattern "^RESULT" | ForEach-Object { $_.Line }
} | Out-File -FilePath "$WT\dev\aterms2-suite-$Tag-summary.txt" -Encoding utf8
