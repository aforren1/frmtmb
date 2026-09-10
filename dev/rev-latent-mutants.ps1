# Reviewer, item 6: has test-hmm-starts.R been SEEN TO FAIL?
#
# The lane records the not-converged and error BRANCHES being reached by
# construction, which is not the same claim. A new function has no
# "unfixed code" to run the tests against, so the equivalent evidence is
# a mutant: break one guarantee in hmm_starts() and check the file goes
# red on that guarantee and not merely on a missing symbol.
#
# Four mutants, each a single edit to a COPY of the package. The lane's
# own tree is never touched.
#
#   M1  the jitter is silently dropped (relist returns the estimates
#       unchanged), which is the guard failing open
#   M2  a refit that did not converge may win
#   M3  the converged spread is taken over every finished refit
#   M4  the incumbent is exempted from the convergence test
#
#   powershell -File dev/rev-latent-mutants.ps1

$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-latent"
$work = "C:\Users\adf44\source\r\rev-latent-mut"
$lib  = "C:\Users\adf44\source\r\rev-latent-mutlib"
$env:R_LIBS = "C:/Users/adf44/source/r/rev-latent-mutlib;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
Set-Location $tree

foreach ($m in @("M1", "M2", "M3", "M4")) {
  Write-Output "================ MUTANT $m ================"
  if (Test-Path $work) { Remove-Item -Recurse -Force $work }
  if (Test-Path $lib)  { Remove-Item -Recurse -Force $lib }
  New-Item -ItemType Directory -Force $work | Out-Null
  New-Item -ItemType Directory -Force $lib  | Out-Null
  Copy-Item -Recurse "$tree\extensions\frmtmb.latent\*" $work
  Rscript "$tree\dev\rev-latent-mutate.R" $work $m
  if ($LASTEXITCODE -ne 0) { Write-Output "  MUTATION FAILED"; continue }
  & "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL --library=$lib --no-multiarch $work 2>&1 | Select-String -Pattern "^\* DONE|ERROR"
  Write-Output "  install exit=$LASTEXITCODE"
  if ($LASTEXITCODE -ne 0) { continue }
  Rscript "$tree\dev\rev-latent-runtests.R" "$work\tests\testthat\test-hmm-starts.R" $lib 2>&1 |
    Select-String -Pattern "^frmtmb.latent from|^BLOCKS|^ASSERT|^FAIL|^ERROR|^SKIP|^  \[|^      "
  Write-Output "  RUN exit=$LASTEXITCODE  (non-zero means the suite CAUGHT this mutant)"
}
if (Test-Path $work) { Remove-Item -Recurse -Force $work }
if (Test-Path $lib)  { Remove-Item -Recurse -Force $lib }
Write-Output "================ done ================"
