# Lane `latent`: run the hmm_starts() cost measurement as soon as the
# replicate arms are done, and not before. The measurement's own
# control assumes a quiet machine, so it waits for the arm driver's
# "ARMS DONE" marker rather than racing it.
#
#   powershell -File dev/latent-2p3-cost.ps1
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
Set-Location C:\Users\adf44\source\r\frmtmb-wt-latent
while (-not (Select-String -Path dev\latent-2p3-run-log.txt -Pattern "ARMS DONE" -Quiet -ErrorAction SilentlyContinue)) {
  Start-Sleep -Seconds 30
}
Write-Output "arms done, starting the cost measurement at $(Get-Date -Format HH:mm)"
Rscript dev/latent-2p3-starts-cost.R 2 3
Write-Output "COST exit=$LASTEXITCODE at $(Get-Date -Format HH:mm)"
Write-Output "COST DONE"
