# Runs one Rscript under a working-set ceiling, because a dense Hessian tape
# can take every free gigabyte of a machine other lanes share. Stops only the
# process it started. Prints the peak working set.
#   powershell -File dev/ncgp-guard.ps1 <limitMB> <timeoutS> <script> <args...>
param([int]$LimitMB, [int]$TimeoutS, [Parameter(ValueFromRemainingArguments=$true)]$Rest)
$rs = "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
$out = [System.IO.Path]::GetTempFileName()
$p = Start-Process -FilePath $rs -ArgumentList $Rest -NoNewWindow -PassThru `
  -RedirectStandardOutput $out -RedirectStandardError "$out.err"
$h = $p.Handle
$peak = 0L
$t0 = Get-Date
$why = "done"
while (-not $p.HasExited) {
  Start-Sleep -Milliseconds 500
  try { $p.Refresh(); $ws = $p.PeakWorkingSet64 } catch { $ws = 0L }
  if ($ws -gt $peak) { $peak = $ws }
  if ($peak -gt $LimitMB * 1MB) { $why = "memlimit"; Stop-Process -Id $p.Id -Force; break }
  if (((Get-Date) - $t0).TotalSeconds -gt $TimeoutS) { $why = "timeout"; Stop-Process -Id $p.Id -Force; break }
}
$el = ((Get-Date) - $t0).TotalSeconds
Get-Content $out | Where-Object { $_ -match '^(NCGP|Error)' -or $_ -match 'rror' }
Get-Content "$out.err" | Where-Object { $_ -match 'rror|bad_alloc' } | Select-Object -First 5
"NCGPGUARD $why peakMB=$([int]($peak / 1MB)) elapsed=$([int]$el)s args=$($Rest -join ' ')"
Remove-Item $out, "$out.err" -ErrorAction SilentlyContinue
