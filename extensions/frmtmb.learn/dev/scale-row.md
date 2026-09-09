# Scale row

This package has a row in the Phase 0 measurement tier of
`dev/extension-gaps-plan.md`. The design, the numbers and what the row
decided are in the repository root's `dev/scale-findings.md`; how the
run was made is in `dev/scale-lane-notes.md`.

| row | design | what it decides |
|---|---|---|
| `learn` | `bandit2arm_delta()`, 100 subjects x 200 trials, `(1 \| p \| id)` on both parameters | whether the per-trial recursion scales to the population designs the literature runs |
| `learn-rlddm` | `rlddm()` on the same frame, `(1 \| p \| id)` on `alpha`, `drift`, `bs` and `ndt` | the same question with a diffusion choice rule |

`bias` is held at 0.5 in the `rlddm` row, which is what the literature
does with a start point between two arms; `?rlddm` says so.

Run it with:

    $env:NOT_CRAN = "true"
    $env:FRMTMB_SCALE_TESTS = "true"
    $env:FRMTMB_SCALE_ROW = "<row>"
    $env:FRMTMB_SCALE_OUT = "<a file to append to>"

then `testthat::test_file("tests/testthat/test-scale.R")`. Without
`FRMTMB_SCALE_TESTS` the file skips. `FRMTMB_SCALE_SMALL=true` runs
every design at a size that checks the code path in about a minute and
marks each recorded row `small=TRUE`.
