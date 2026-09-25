# Run one recovery seed and report the R child's peak working set.
# Usage: powershell -File dev/phase3b-peakmem.ps1 <lib> <arm> <seed>
param([string]$lib, [string]$arm, [string]$seed)
Set-Location C:/Users/adf44/source/r/frmtmb-wt-phase3b
$env:P3B_LIB = $lib
$p = Start-Process -FilePath "C:/Program Files/R/R-4.6.1/bin/x64/Rscript.exe" `
  -ArgumentList "dev/phase3b-eam-recovery7.R", $arm, $seed -PassThru -NoNewWindow `
  -RedirectStandardOutput "dev/phase3b-log/peakmem-$arm-$seed.out" `
  -RedirectStandardError "dev/phase3b-log/peakmem-$arm-$seed.err"
$peak = 0
while (-not $p.HasExited) {
  try { $p.Refresh(); if ($p.PeakWorkingSet64 -gt $peak) { $peak = $p.PeakWorkingSet64 } } catch {}
  Start-Sleep -Milliseconds 500
}
"{0} {1}: peak working set {2:N2} GB" -f $arm, $seed, ($peak / 1GB)
Get-Content "dev/phase3b-log/peakmem-$arm-$seed.out"
