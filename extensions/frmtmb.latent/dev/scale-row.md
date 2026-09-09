# Scale row

This package has a row in the Phase 0 measurement tier of
`dev/extension-gaps-plan.md`. The design, the numbers and what the row
decided are in the repository root's `dev/scale-findings.md`; how the
run was made is in `dev/scale-lane-notes.md`.

| row | design | what it decides |
|---|---|---|
| `hmm` | K = 3 gaussian, 50 sequences x 500 steps, `tr12 ~ (1 \| id)` | the cost of the post-fit passes, `hmm_probs()` and `hmm_viterbi()`, which are R loops over sequences |
| `lca` | K = 4, n = 2000, 10 binary items, two covariates on membership | included for the table's completeness |

The `hmm` row uses `init = "uniform"`. The default stationary initial
distribution is refused when a transition carries a predictor, and the
row's simulator draws the first state uniformly, so uniform is the
correct model here and not only the allowed one.

Run it with:

    $env:NOT_CRAN = "true"
    $env:FRMTMB_SCALE_TESTS = "true"
    $env:FRMTMB_SCALE_ROW = "<row>"
    $env:FRMTMB_SCALE_OUT = "<a file to append to>"

then `testthat::test_file("tests/testthat/test-scale.R")`. Without
`FRMTMB_SCALE_TESTS` the file skips. `FRMTMB_SCALE_SMALL=true` runs
every design at a size that checks the code path in about a minute and
marks each recorded row `small=TRUE`.
