# Scale row

This package has a row in the Phase 0 measurement tier of
`dev/extension-gaps-plan.md`. The design, the numbers and what the row
decided are in the repository root's `dev/scale-findings.md`; how the
run was made is in `dev/scale-lane-notes.md`.

| row | design | what it decides |
|---|---|---|
| `eam` | `wiener()`, 30 subjects x 400 trials, 2 conditions, random effects on `mu`, `bs` and `ndt` | whether a hierarchical DDM is minutes or hours |
| `eam-sv` | the same with `wiener(variability = "sv")` | what across-trial drift variability adds |
| `eam-unbounded` | the same data with `wiener(max_ndt = 0.45, allow_unreachable = TRUE)` | whether the non-decision-time BOUND is what the `eam` row's trouble is |

The second question the `eam` row answers is whether the bounded
non-decision-time link tolerates a random effect at all, which is the
trouble item 2.1 of the plan expects. The row records the fitted
between-subject spread of `ndt` on the natural scale beside the
spread the simulator used, so that a collapsed component is visible
rather than inferred.

The `eam-unbounded` arm exists because "the fit did not converge" is
not a cause. It fits the SAME data with the bound raised above every
subject's true non-decision time, so the two arms differ in the bound
and in nothing else.

**It is the diagnosis, not the fix. Do not copy it into a script.**
`allow_unreachable = TRUE` does not add a per-subject bound; it lifts
the refusal that stops a model admitting rows it cannot reach, and
under it an unreachable row is evaluated at a decision time of
`1e-9 * min(y)` rather than refused. This arm is valid only because its
optimum landed where no row is unreachable, which was checked (the
smallest `y - ndt` over all 12,000 rows is 0.0219 s and no row is
floored) and is a property of this fit rather than a guarantee. Item
1.0a is what would make it a guarantee.

The timing convention this tier follows is the one in
`dev/dev-findings.md`, "How this package times things".

Run it with:

    $env:NOT_CRAN = "true"
    $env:FRMTMB_SCALE_TESTS = "true"
    $env:FRMTMB_SCALE_ROW = "<row>"
    $env:FRMTMB_SCALE_OUT = "<a file to append to>"

then `testthat::test_file("tests/testthat/test-scale.R")`. Without
`FRMTMB_SCALE_TESTS` the file skips. `FRMTMB_SCALE_SMALL=true` runs
every design at a size that checks the code path in about a minute and
marks each recorded row `small=TRUE`.

## Disclosure to ship with the next release

Phase 0 found a defect a user hits today, so this is written now rather
than when item 1.0a lands. It is a NEWS bullet for this package,
drafted to be lifted verbatim into the next version's section. This
lane did not put it in `NEWS.md` itself, because that file is organized
by released version and choosing the version is not a measurement
lane's call.

> `wiener()` bounds the non-decision time by `min(rt)`, the GLOBAL
> fastest response in the data, through a scaled logit. With a random
> effect on `ndt` this is the wrong constraint: the information about a
> subject's non-decision time is that subject's own fastest response,
> and a subject whose true `ndt` is above the global minimum cannot be
> represented at any value of the random effect. At 30 subjects x 400
> trials with a between-subject `ndt` spread of 26 ms on a mean of
> 250 ms, 20 of the 30 subjects are in that position while NONE is
> inconsistent with its own data. Without `variability` the fit does
> not converge and says so: maximum absolute gradient 1.3e11, Hessian
> not positive definite, all seven standard errors NaN. With
> `variability = "sv"` it CONVERGES: code 0, maximum gradient 7.5e-05,
> positive definite Hessian, no bad standard errors, and `diagnose()`
> reports nothing. The population non-decision time then comes back
> pinned at the bound, 0.2236 against a truth of 0.25 and 0.42 of a
> standard error below `min(rt)`, with a delta-method standard error of
> 7.2e-06 on it, because the scaled logit's derivative vanishes where
> the estimate has been pushed. Do not put a random effect on `ndt`
> until the per-subject bound lands; the same defect reaches `rlddm()`
> in frmtmb.learn, which takes its diffusion parameterization from this
> package. Fitting the same data with the bound raised above every
> subject's truth recovers everything and finds a log likelihood 121.4
> units higher at the same parameter count.
