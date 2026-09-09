# Scale row

This package has a row in the Phase 0 measurement tier of
`dev/extension-gaps-plan.md`. The design, the numbers and what the row
decided are in the repository root's `dev/scale-findings.md`; how the
run was made is in `dev/scale-lane-notes.md`.

| row | design | what it decides |
|---|---|---|
| `coupling-coh-intercept` | `coh ~ 1` | the floor |
| `coupling-coh-cond` | `coh ~ cond` | what a condition effect costs |
| `coupling-coh-smooth` | `coh ~ cond + s(freq, by = cond)` | what the smooth costs |
| `coupling-coh-id` | the same `+ (1 \| id)` | the survey's 52.5 s model, and the one its 0.77 contrast came from |
| `coupling-coh-full` | the same `+ (1 \| id:cond)` | the survey's 28.9 s model, and the one the plan's item 2.6 names |

All five run on one design of 40 subjects x 2 conditions x 60
frequencies. Only the top rung is asserted to recover the condition
contrast, and only at the realistic size: the rungs below it omit
random effects the truth has, and at the small size the assertion has
no power.

The plan asked for "the four scratchpad benchmarks from the survey" to
be promoted here. Their constructions are in the survey's session
transcript, and the top two rungs are the survey's third and fourth
fits, so the promotion is done. The survey ran 0.53.0 under
`load_all()` and this runs 0.55.0 installed, so no speedup is claimed
from the pairs. `dev/scale-lane-notes.md` records how the first draft
of this file got that wrong.

Run it with:

    $env:NOT_CRAN = "true"
    $env:FRMTMB_SCALE_TESTS = "true"
    $env:FRMTMB_SCALE_ROW = "<row>"
    $env:FRMTMB_SCALE_OUT = "<a file to append to>"

then `testthat::test_file("tests/testthat/test-scale.R")`. Without
`FRMTMB_SCALE_TESTS` the file skips. `FRMTMB_SCALE_SMALL=true` runs
every design at a size that checks the code path in about a minute and
marks each recorded row `small=TRUE`.
