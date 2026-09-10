# Lane `latent`, item 2.3: the two replicate arms, SEQUENTIALLY.
#
# Sequential and not parallel. Two of these processes on the same
# machine made each replicate roughly five times slower than the single
# replicate they were priced from, and a shared machine also makes
# every second in the output tables meaningless. The timing CLAIM this
# row makes comes from dev/latent-2p3-starts-cost.R, which runs alone
# and carries its own control; the seconds in these tables are
# indicative and are labelled as such.
#
# Arm A's identity replicates were measured separately, in
# dev/latent-2p3-hmm-A-id.tsv, because depmixS4's random-start EM cost
# 18, 263 and 645 seconds on the first three replicates at
# tol = 1e-12 and maxit = 5000. That variance is the reference's, not
# the model's, and paying it forty times buys nothing: an identity is
# a property of the two codes.
#
#   powershell -File dev/latent-2p3-run.ps1
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
Set-Location C:\Users\adf44\source\r\frmtmb-wt-latent
Write-Output "A recovery: 40 replicates, no third-party arm, $(Get-Date -Format HH:mm)"
Rscript dev/latent-2p3-hmm.R A 40 dev/latent-2p3-hmm-A.tsv 0
Write-Output "A exit=$LASTEXITCODE at $(Get-Date -Format HH:mm)"
Write-Output "B: 15 replicates, hmmTMB on the first 6, $(Get-Date -Format HH:mm)"
Rscript dev/latent-2p3-hmm.R B 15 dev/latent-2p3-hmm-B.tsv 6
Write-Output "B exit=$LASTEXITCODE at $(Get-Date -Format HH:mm)"
Write-Output "ARMS DONE"
