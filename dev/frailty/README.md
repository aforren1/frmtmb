# The scripts behind dev/frailty-findings.md

Item 2.5 of `dev/extension-gaps-plan.md`: the Royston-Parmar frailty
against rstpm2's log-normal frailty, and a random effect on `gamma1`.

They live in the repository rather than in a session scratchpad because
a number without its construction gets re-derived by whoever reads it
next, and that has twice produced a confident and wrong claim that a
shipped figure was false. `dev/lane-rules.md` records the same lesson
about a process document that was cleaned mid-round.

## To run one

Every script `source()`s its neighbours by a BARE name, so run from
inside this directory:

    cd dev/frailty
    Rscript frailty-recover.R

`frailty-lib.R` sets `.libPaths()` and carries the paths of the machine
the round was run on, and `frailty-runtest.R` carries the path of the
worktree. Edit both before running anything. Nothing here installs a
package.

Each script writes its raw rows next to itself as a `.tsv` of
`key=value` pairs separated by tabs, and its `-report.R` companion
turns that into the table the findings document quotes.

The rows of the three replicate runs are KEPT here rather than
regenerated: `frailty-recover.tsv` (200 replicates, 26 minutes),
`frailty-gamma1.tsv` (60) and `frailty-sweep.tsv` (30 cells). Rerunning
a script overwrites its own file, so copy one aside first if you want
to compare. A `-report.R` script reads the `.tsv` and fits nothing, so
any table in `dev/frailty-findings.md` can be checked in a second.

## What each one produces

| script | produces |
|---|---|
| `frailty-lib.R` | the library order. Sourced by every other script |
| `frailty-common.R` | the simulator, the Royston-Parmar basis written out by hand, the exact per-cluster marginal likelihood, and the change of basis from rstpm2's `nsx` coordinates |
| `frailty-ident.R` | the first side-by-side fit, the three integration rules at one parameter point, and rstpm2 at 3, 9, 21 and 51 quadrature nodes |
| `frailty-sweep.R`, `frailty-sweep-report.R` | how the agreement moves with cluster size and with the frailty's own size |
| `frailty-recover.R`, `frailty-recover-report.R` | the 200-replicate recovery table and the agreement with rstpm2 at the plan's realistic design |
| `frailty-gamma1-recover.R`, `frailty-gamma1-report.R` | the 60-replicate table for a random effect on `gamma1`, with the component switched off as the second arm |
| `frailty-gamma1-units2.R` | the time-unit measurement: a slope-only block moves, a paired block does not |
| `frailty-falsify.R` | every guard in `test-frailty.R` run against the case it exists to catch |
| `frailty-floor.R` | 15 fits on a design built to reach the monotonicity floor, which did not reach it |
| `frailty-floor2.R` | the constructed floor case, which pins `rp_floored()` to the centre deviations rather than the population `gamma1` |
| `frailty-barrier.R` | the size of the monotonicity penalty per event, from `sp_floor_pos()`'s own arithmetic |
| `frailty-floor3.R` | the all-censored centre: the case this lane asserted was safe and the punch round falsified, with its control |
| `frailty-parorder.R` | the parameter-ordering guard that `frailty-sweep.R` and `frailty-ident.R` now carry, run against three wrong orderings |
| `frailty-check.sh` | `R CMD build` then `R CMD check --as-cran`, which a lane runs ONCE, on its final pass |
| `frailty-runtest.R`, `frailty-suite.sh` | the test runner, one test file per R process. `frailty-runtest.R` evaluates each file in `testthat::test_env("frmtmb.spline")`, which is what `test_check()` does; a bare `test_file()` cannot see the package's internal functions and reports an ERROR that belongs to the runner |
