# Scale row

This package has a row in the Phase 0 measurement tier of
`dev/extension-gaps-plan.md`. The design, the numbers and what the row
decided are in the repository root's `dev/scale-findings.md`; how the
run was made is in `dev/scale-lane-notes.md`.

| row | design | what it decides |
|---|---|---|
| `ode` | 100 subjects x 8 samples, the depot and central system, twice-daily dosing for 7 days through `ii`/`addl`, one `ss` row | whether the segmented sensitivity solve is usable at population scale, or whether the closed-form path (the plan's item 5.4) is a prerequisite rather than an option |

This is the row the plan singles out: "Whether 5.4 is Phase 5 or Phase
1 is what Phase 0 decides."

Run it with:

    $env:NOT_CRAN = "true"
    $env:FRMTMB_SCALE_TESTS = "true"
    $env:FRMTMB_SCALE_ROW = "<row>"
    $env:FRMTMB_SCALE_OUT = "<a file to append to>"

then `testthat::test_file("tests/testthat/test-scale.R")`. Without
`FRMTMB_SCALE_TESTS` the file skips. `FRMTMB_SCALE_SMALL=true` runs
every design at a size that checks the code path in about a minute and
marks each recorded row `small=TRUE`.
