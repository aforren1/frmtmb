# The full local validation, in the order the cases need.
#
# 1. the eight-site build with the FINAL dev/release/build-docs.R, so
#    that the numbers in dev/docsci-findings.md come from the script that
#    is being delivered and not from an earlier draft;
# 2. the guard cases, which bend the _pkgdown.yml files and so must not
#    run beside a build that reads them.
$ErrorActionPreference = "Continue"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-docsci"

Write-Output "########## FULL BUILD"
& "$ROOT/dev/docsci-run-local.ps1"
Write-Output "########## FULL BUILD done"

Write-Output "########## GUARD CASES"
& "$ROOT/dev/docsci-guard-tests.ps1"
Write-Output "########## GUARD CASES done"
