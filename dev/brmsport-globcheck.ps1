# The brms-suite glob of dev/release/run-gated.ps1, lifted out and run on
# the worktree (both packages must match) and on an empty directory
# (must throw), so the fail-closed branch is seen firing and not firing.
#   powershell -File dev/brmsport-globcheck.ps1
$ErrorActionPreference = "Stop"
$ROOT = (Get-Location).Path -replace "\\", "/"
$text = Get-Content "$ROOT/dev/release/run-gated.ps1" -Raw
$start = $text.IndexOf('$suite = @(')
$end = $text.IndexOf('[void]$jobs.Add(@{ n = "frmtmb.learn"')
if ($start -lt 0 -or $end -lt 0) { throw "glob block not found" }
$block = $text.Substring($start, $end - $start)

function Run-Glob($core, $root) {
  $jobs = New-Object System.Collections.ArrayList
  $ROOT = $root
  Invoke-Expression $block
  return $jobs.Count
}

$n = Run-Glob "$ROOT/tests/testthat" $ROOT
Write-Output ("worktree: " + $n + " jobs")
$empty = "$env:TEMP/brmsport-emptyglob"
New-Item -ItemType Directory -Force $empty | Out-Null
New-Item -ItemType Directory -Force "$empty/extensions/frmtmb.sample/tests/testthat" | Out-Null
try {
  $m = Run-Glob $empty $empty
  Write-Output ("empty: " + $m + " jobs, DID NOT THROW")
} catch {
  Write-Output ("empty: threw: " + $_.Exception.Message)
}
