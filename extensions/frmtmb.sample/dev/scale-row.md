# Scale row

This package has a row in the Phase 0 measurement tier of
`dev/extension-gaps-plan.md`. The design, the numbers and what the row
decided are in the repository root's `dev/scale-findings.md`; how the
run was made is in `dev/scale-lane-notes.md`.

| row | design | what it decides |
|---|---|---|
| `sample` | a 2000-row binomial GLMM with two crossed factors, 4 chains x 2000, then `posterior_epred()`, `posterior_predict()`, `ranef()`, `log_lik()` and `loo()` | how slow the per-draw R loops are, and whether the plan's item 4.6 caching is needed before anything else in this package |

`posterior_epred()` is timed at the full draw count and at a tenth of
it, because a per-draw loop and a cached design differ in their SLOPE
in the number of draws and not only in their level.

Run it with:

    $env:NOT_CRAN = "true"
    $env:FRMTMB_SCALE_TESTS = "true"
    $env:FRMTMB_SCALE_ROW = "<row>"
    $env:FRMTMB_SCALE_OUT = "<a file to append to>"

then `testthat::test_file("tests/testthat/test-scale.R")`. Without
`FRMTMB_SCALE_TESTS` the file skips. `FRMTMB_SCALE_SMALL=true` runs
every design at a size that checks the code path in about a minute and
marks each recorded row `small=TRUE`.
