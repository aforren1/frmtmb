# Lane `latent`, punch round 2: the weight evidence, in order.
#
#   1. 200 FRESH seeds with their own poLCA arm, and the whole grid;
#   2. the lane's own 200 seeds at the three weights the round-one
#      sweep did not reach.
#
# Sequential, so that neither run's counts depend on what else the
# machine was doing. Nothing here is a timing claim.
#
#   powershell -File dev/latent-2p4-weight.ps1
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
Set-Location C:\Users\adf44\source\r\frmtmb-wt-latent
Write-Output "OUT-OF-SAMPLE: 200 seeds from 20270401, $(Get-Date -Format HH:mm)"
Rscript dev/latent-2p4-oos.R 20270401 200 dev/latent-2p4-oos.tsv
Write-Output "oos exit=$LASTEXITCODE at $(Get-Date -Format HH:mm)"
Write-Output "IN-SAMPLE right end: w = 0.95, 0.99, 1"
Rscript dev/latent-2p4-shrink2.R 200
Write-Output "shrink2 exit=$LASTEXITCODE at $(Get-Date -Format HH:mm)"
Write-Output "WEIGHT EVIDENCE DONE"
