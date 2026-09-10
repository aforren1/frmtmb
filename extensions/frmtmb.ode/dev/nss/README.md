# The steady-state run-in measurement scripts

Every number in `dev/nss-findings.md` (repository root) was produced by
a script here. They are RUN scripts, not tests: they print tables
rather than assert, and several take minutes. The assertions that guard
the same behavior in the shipped package are in
`tests/testthat/test-ode-events.R`, under "the steady-state tail".

`prelude.R` sets `.libPaths()` and is sourced by all of them. With
`NSS_ARM=ref` in the environment it puts the round's shared reference
build of the base commit first instead of the lane's library, so the
same script measures the BEFORE arm. `frm_ode()` passes unknown
arguments to the integrator through `...`, so the reference build
accepts `ss_extrapolate` without complaint and ignores it; `nss-09`
relies on that, while `nss-12` and `nss-13` withhold the argument on
that arm rather than depend on it.

## What each one does

| script | what it measures |
|---|---|
| `nss-00-env.R` | the library order and the versions, run before believing anything |
| `nss-01-getvalues.R` | whether a tape-time numeric value can be read off an advector, which is what `n_ss = "auto"` would need |
| `nss-02-accessors.R` | which RTMB accessors exist for the same purpose |
| `nss-03-ratio.R` | the contraction ratio against `exp(-lambda_z * ii)` on exact matrix exponentials, the `n_ss` "auto" would pick, and what extrapolation buys with and without solver noise |
| `nss-04-hazards.R` | four constructions where an unguarded correction is worse than truncating |
| `nss-05-guard.R` | eight candidate guards over nineteen systems, scored on the worst row |
| `nss-06-rule.R` | the stand-down rule against a hard cap, per state against shared, which settles the two constants |
| `nss-07-adops.R` | which of the guard's operations RTMB can tape, and their derivatives |
| `nss-08-rmove.R` | how far the contraction ratio moves between the starting values and the optimum: 54 fits, three designs x three starts x six seeds. This is the measurement that rejects `n_ss = "auto"` |
| `nss-09-verify.R` | the change through `frm_ode()` itself: the compartment grid, the guard constructions, bit-for-bit back compatibility against the reference build, and the solve count |
| `nss-10-accept.R` | the acceptance design fitted through `frm_ode()` itself. Killed three times by the session before finishing; `nss-19-objgrad.R` answers the same question for the cost of one tape build |
| `nss-11-grad.R` | the gradient through the correction, against a finite difference of its own numeric path and against the exact steady state |
| `nss-12-warn.R` | the false-alarm rate of the warning on 22 ordinary population pharmacokinetic schedules, classified against the exact limit |
| `nss-13-warntext.R` | what the warning SAYS, base against lane, on the rows the review measured its understatement on |
| `nss-14-absent.R` | the bars the new tests assert, computed with each piece of the guard removed |
| `nss-15-seefail.sh` | builds four variants of the package with one property each removed and runs `test-ode-events.R` against them, so every new test has been seen to fail. Uses `nss-15-variant.py` and `nss-15-run.R` |
| `nss-16-nonconsec.R` | a defect found and not fixed: a derivative written `0 * y[k]` cannot be solved |
| `nss-17-accept-lin.R` | the acceptance criterion at a 107 h terminal half-life, six seeds |
| `nss-18-accept-265.R` | the same at 264.6 h, which is the design the review's `k21` factor of 3.16 was actually measured on |
| `nss-19-objgrad.R` | `frm_ode()` against `frm_lincmt()` on the acceptance dataset's whole objective and gradient, which licenses the closed form as the vehicle in `nss-17` |
| `nss-20-failopen.R` | punch round 1: the reviewer's three fail-open constructions, run against the fix on both arms |
| `nss-21-tape.R` | punch round 1: AD tape node counts, both arms, against the base commit |
| `nss-22-kink.R` | punch round 1: whether the objective's one-sided derivatives converge across the stand-down point, and what the smoothing costs |
| `nss-23-flipflop.R` | punch round 1: where extrapolation loses to truncation, and the bound on it |
| `nss-24-backcompat.R` | punch round 1: `ss_extrapolate = FALSE` against the base commit on 26 awkward schedules |
| `nss-25-dominance.R` | punch round 1: a floored `n_ss = "auto"` against the correction, error and `lsoda` entries side by side |
| `nss-26-corners.R` | punch round 2: the smoothness of the correction FACTOR at its four junctions, measured on the shape rather than on an objective, with the curvature control that separates a corner from ordinary curvature |
| `nss-27-largenss.R` | punch round 2: whether the warning still fires at `n_ss` in the hundreds, inside the stand-down band, which is where `?frm_ode` sends a long-half-life user |
| `nss-15-one.sh` | one see-it-fail variant with its own work directory and library. Use this rather than starting `nss-15-seefail.sh` twice: two copies share `$WORK` and deadlock |
| `nss-suite.sh` | one test file per R process, with `nss-one.R`; `nss-check.sh` runs `R CMD check --as-cran` once |
