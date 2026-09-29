# Scale row

This package has a row in the Phase 0 measurement tier of
`dev/extension-gaps-plan.md`. The design, the numbers and what the row
decided are in the repository root's `dev/scale-findings.md`; how the
run was made is in `dev/scale-lane-notes.md`.

| row | design | what it decides |
|---|---|---|
| `spline` | `s(t) + s(t, subject, bs = "fs")` at 40 subjects x 200 points, then `frm_curve()` and `frm_curve_feature()` | the cost of the feature search, which calls `predict()` once per Newton step |

The memoized joint-precision solve behind `frm_curve()` was already
measured at 6.9 s on 8006 random coefficients, which is why this row
times the curve and the feature separately from the fit.

That memoization is also the reason the row measures a COLD call
before it measures anything else. The solve is cached on the fit
object and paid once, by whichever curve call comes first, so every
later curve call on that fit is cheap and any "first call" measured
after another curve call is a warm number. The row therefore calls
`frm_curve_feature()` on the untouched fit to price the solve, and
interleaves the rest afterwards. Measured: the solve 0.48 s, then the
simultaneous band 0.071 s, the feature search 0.0098 s and the
pointwise band 0.0031 s.

The Royston and Parmar design in the plan's "Realistic scale" table is
not in this tier: the Phase 0 table's spline row names the curves
design only, and the survival design belongs to the plan's item 2.5.

Run it with:

    $env:NOT_CRAN = "true"
    $env:FRMTMB_SCALE_TESTS = "true"
    $env:FRMTMB_SCALE_ROW = "<row>"
    $env:FRMTMB_SCALE_OUT = "<a file to append to>"

then `testthat::test_file("tests/testthat/test-scale.R")`. Without
`FRMTMB_SCALE_TESTS` the file skips. `FRMTMB_SCALE_SMALL=true` runs
every design at a size that checks the code path in about a minute and
marks each recorded row `small=TRUE`.
