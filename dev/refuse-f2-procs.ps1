# This lane's R processes, by their script: list them, or stop one by a pattern
# that must name a dev/refuse-f2 script. Usage: -Stop <pattern>
param([string]$Stop = "")
$p = Get-CimInstance Win32_Process -Filter "Name='Rterm.exe' or Name='Rscript.exe'" |
  Where-Object { $_.CommandLine -match 'dev/refuse-f2' }
foreach ($x in $p) {
  $m = [regex]::Match($x.CommandLine, 'dev/refuse-f2[^ ]*( [A-Za-z_.]+)?').Value
  if ($Stop -ne "" -and $x.CommandLine -match $Stop) {
    Stop-Process -Id $x.ProcessId -Force
    "stopped $($x.ProcessId) $m"
  } else {
    "$($x.ProcessId) $m"
  }
}
