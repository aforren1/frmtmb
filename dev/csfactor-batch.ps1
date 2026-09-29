# Run one test file per R process, up to $MaxJobs at a time, and append
# each RESULT line to one log. Usage:
#   powershell -File dev/csfactor-batch.ps1 -Lib <lib> -Log <log> -Files a,b
param(
  [string]$Lib = "C:/Users/adf44/source/r/wt-csfactor-lib",
  [string]$Log = "C:/Users/adf44/source/r/frmtmb-wt-csfactor/dev/csfactor-log/batch.txt",
  [string[]]$Files,
  [int]$MaxJobs = 7,
  [switch]$Gated
)
$ErrorActionPreference = "Continue"
$WT = "C:/Users/adf44/source/r/frmtmb-wt-csfactor"
$Rscript = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$outDir = Split-Path $Log -Parent
if (-not (Test-Path $outDir)) { New-Item -ItemType Directory $outDir | Out-Null }
[System.IO.File]::WriteAllText($Log, "", (New-Object System.Text.UTF8Encoding($false)))
$jobs = New-Object System.Collections.ArrayList
foreach ($f in $Files) {
  while ((@(Get-Job -State Running)).Count -ge $MaxJobs) { Start-Sleep -Milliseconds 500 }
  $j = Start-Job -ScriptBlock {
    param($Rscript, $WT, $Lib, $f, $Gated)
    $env:NOT_CRAN = "true"
    if ($Gated) {
      $env:FRMTMB_BRMS_FIT_TESTS = "true"
      $env:FRMTMB_STAN_CACHE = "$WT/dev/csfactor-stan-cache"
      $env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
    }
    $o = & $Rscript "$WT/dev/csfactor-run-tests.R" $Lib $f 2>&1
    $o -join "`n"
  } -ArgumentList $Rscript, $WT, $Lib, $f, $Gated.IsPresent
  [void]$jobs.Add(@{ job = $j; name = $f })
}
$ran = 0
foreach ($e in $jobs) {
  $null = Wait-Job $e.job
  $txt = Receive-Job $e.job
  Remove-Job $e.job
  $lines = ($txt -split "`n") | Where-Object { $_ -match "^(RESULT|  BAD|  SKIP|     )" }
  if (-not ($lines | Where-Object { $_ -match "^RESULT" })) {
    Add-Content -Path $Log -Value ("RESULT " + $e.name + " NO-RESULT-LINE")
    Add-Content -Path $Log -Value (($txt -split "`n" | Select-Object -Last 25) -join "`n")
  } else {
    Add-Content -Path $Log -Value ($lines -join "`n")
    $ran = $ran + 1
  }
}
Add-Content -Path $Log -Value ("BATCH produced a RESULT line for " + $ran + " of " + $Files.Count + " files")
