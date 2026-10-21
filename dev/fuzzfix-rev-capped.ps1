# Reviewer: run one Rscript with a memory cap; kill it (and only it) if
# its working set passes the cap. Prints the peak working set.
#   powershell -File dev/fuzzfix-rev-capped.ps1 <capMB> <outfile> <args...>
param([int]$capMB, [string]$out, [Parameter(ValueFromRemainingArguments)]$rest)
$rs = "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
$argline = ($rest | ForEach-Object { '"' + $_ + '"' }) -join ' '
$p = Start-Process -FilePath $rs -ArgumentList $argline -NoNewWindow `
  -PassThru -RedirectStandardOutput $out -RedirectStandardError ($out + ".err")
$peak = 0; $killed = $false
while (-not $p.HasExited) {
  Start-Sleep -Milliseconds 500
  try { $p.Refresh(); $ws = [int]($p.WorkingSet64 / 1MB) } catch { $ws = 0 }
  if ($ws -gt $peak) { $peak = $ws }
  if ($ws -gt $capMB) { Stop-Process -Id $p.Id -Force; $killed = $true; break }
}
Add-Content -Path $out -Value ("CAPPED peakMB {0} cap {1} killed {2}" -f $peak, $capMB, $killed)
