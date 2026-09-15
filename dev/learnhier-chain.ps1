# Lane `learnhier`: run the remaining arms back to back, at the cap,
# without a person in the loop between them.
#
#   dev/learnhier-chain.ps1
#
# ORDER IS BY PRIORITY, so that an interruption loses the least. The
# 60-replicate main arm is item 2.2's mandate and is already running;
# this waits for it rather than competing with it. Then Test 2, which is
# the one experiment that separates the whole remaining candidate set
# for the `sd(ndt)` bias from the alternative, because any bias driven
# by approximation or by information shrinks with the trial count and
# only a structural one survives it. Then the importance arm, which has
# a partial answer already and gains only precision.
#
# The trial-count arms write to their OWN directories. A record is named
# `rlddm-<seed>.rds` with no trial count in it, so three arms sharing a
# directory would silently overwrite each other, and the skip check
# would then refuse to run the ones that matter.
#
# The arms are PAIRED on their learners: lh_rlddm_data() draws the block
# from `set.seed(seed + 2)` before it simulates anything, so the same
# seed gives the same hundred learners at every trial count. Only the
# data they produce, and therefore their floors, differ.
$tree = "C:\Users\adf44\source\r\frmtmb-wt-learnhier"
Set-Location $tree
$rec = "$tree\dev\learnhier-rec"

function Wait-Arm([string]$dir, [string]$pat, [int]$want) {
  while (@(Get-ChildItem $dir -Filter $pat -EA SilentlyContinue).Count -lt $want) {
    Start-Sleep -Seconds 60
  }
}

Write-Output ("CHAIN waiting for the main arm at " + (Get-Date -Format "HH:mm"))
Wait-Arm $rec "rlddm-*.rds" 60
Write-Output ("CHAIN main arm complete at " + (Get-Date -Format "HH:mm"))

# Test 2, cheap point first: fewer trials per learner, where the
# committed prediction is a LARGER relative bias in sd(ndt).
& powershell -NoProfile -File "$tree\dev\learnhier-launch2.ps1" `
  -kind run -design rlddm -n 12 -outdir "$tree\dev\learnhier-rec-nt100" `
  -cap 6 -stagger 30 -nt 100 -runner "dev/learnhier-run2.R"
Write-Output ("CHAIN nt=100 done at " + (Get-Date -Format "HH:mm"))

# Test 2's 400-trial arm is DROPPED, and named with its cost: 12
# replicates, about 2.8 hours. Its purpose was to see whether the
# ESTIMATE converges on a better target as trials rise, and
# dev/learnhier-mvsn.R has since shown the TARGET does not move at all
# across 100, 200 and 400 trials (nuisance share 107.9, 106.7, 108.2
# percent, paired difference +0.0023, p = 0.895), so there is nothing
# for an estimate to converge on. The 100-trial arm above still gives
# the estimate-side comparison at a quarter of the cost.

# The importance arm, to the count the item's success rate needs.
& powershell -NoProfile -File "$tree\dev\learnhier-launch2.ps1" `
  -kind imp -design bandit -n 24 -draws 100 `
  -outdir "$tree\dev\learnhier-imp" -cap 6 -stagger 20
Write-Output ("CHAIN importance done at " + (Get-Date -Format "HH:mm"))

# The five-parameter block, which is the other reading of "every
# parameter": the same data with `bias` in the block too, whose truth
# is a variance of exactly zero. A probe on two seeds, not an arm.
& powershell -NoProfile -File "$tree\dev\learnhier-launch2.ps1" `
  -kind run -design rlddm5 -n 2 -outdir "$tree\dev\learnhier-rec5" `
  -cap 6 -stagger 60 -runner "dev/learnhier-run2.R"
Write-Output ("CHAIN ALL ARMS DONE at " + (Get-Date -Format "HH:mm"))
