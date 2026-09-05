# GDDM reference lane: findings

Lane: `wt-gddm-ref`. Goal: an external, frozen reference for the
generalized components of `gddm()` that nothing outside the solver
checks today.

## 1. What PyDDM actually does (read from 0.9.0, not from the paper)

Read out of the installed source, because every one of these decides a
parameter translation or a tolerance.

* `Solution.choice_upper` / `.choice_lower` are **per-bin probability
  mass**, not densities. `Solution.pdf(choice)` is `mass / dt`.
  `Solution.prob(choice)` is `sum(mass)`, defective: the two plus
  `prob_undecided()` sum to one.
* `Model.x_domain()` rounds the bound **up** to a multiple of `dx`:
  `B = ceil(B/dx)*dx`. A boundary that is not an exact multiple of `dx`
  is silently a different model. Every case here uses `bs = 2`, so
  `B = 1`, which every `dx` in the sweep divides.
* `ICPoint.get_IC()` snaps the start to `round(x0/dx)`. `x0` is an exact
  multiple of every `dx` in the sweep.
* `ICRange.get_IC()` snaps the half width to `int(sz/dx)` cells and is
  centred on the middle of the domain, so it carries no bias. The
  start-variability case is therefore unbiased (`bias = 0.5`), and the
  point-start bias is carried by a different case.
* `OverlayNonDecision.apply()` shifts by `int(ndt/dt)` whole bins, a
  hard truncation. Every `ndt` here is an exact multiple of every `dt`
  in the sweep.
* `OverlayUniformMixture.apply()` adds `0.5 * c / len(t_domain)` of mass
  to each bin. `len(t_domain)` is `T_dur/dt + 1`, so as a density the
  lapse floor is `0.5 c / (T_dur + dt)`, while `gddm()` uses
  `0.5 c / t_max`. See the lapse case below.
* `Model.can_solve_cn()` returns FALSE whenever the bound depends on
  time. **PyDDM has no Crank-Nicolson path for a collapsing bound**; it
  falls back to `solve_numerical_implicit`, backward Euler, which is
  first order in `dt`. The lane's premise ("the two Crank-Nicolson
  solvers") is wrong for the two collapsing cases, and the convergence
  sweep has to be run against a first-order scheme there.
* `Model.has_analytical_solution()` is TRUE for a constant bound with a
  point start, and also for `BoundCollapsingLinear` with a point start
  (Anderson 1960). Those cases get a closed-form reference recorded
  alongside the numerical one.
* PyDDM's collapsing-bound solver sandwiches the moving bound between
  the two surrounding grid nodes and weights them linearly. That is
  exactly the treatment `gddm()`'s change of variable exists to avoid,
  so the two collapsing cases are the ones where a shared bias is least
  likely and a disagreement most likely.

Resolved environment: pyddm 0.9.0, numpy 2.5.2, scipy 1.18.1 on
CPython 3.13 (`requires-python = "==3.13.*"`).

## 2. Parameter translation

| `gddm()` | PyDDM |
| --- | --- |
| `dx = a(x,t) dt + dW` | `noise = NoiseConstant(noise=1)` |
| bounds at `+/- B(t)`, `B(0) = bs/2` | `B = bs/2` |
| `bias` in `(0,1)`, fraction above the lower bound | `ICPoint(x0 = (bs/2)(2 bias - 1))` |
| `gddm_drift_constant()`, `a = mu` | `DriftConstant(drift = mu)` |
| `gddm_drift_leak()`, `a = -leak x` | `DriftLinear(drift = mu, x = -leak, t = 0)` |
| `gddm_bound_exponential()`, `B = (bs/2) exp(-t/tau)` | `BoundCollapsingExponential(B = bs/2, tau = 1/tau)` |
| `gddm_bound_linear()`, `B = (bs/2)(1 - kappa t/t_max)` | `BoundCollapsingLinear(B = bs/2, t = (bs/2) kappa / t_max)` |
| `gddm_drift_coherence()`, `a = sign(C) mu (abs(C)/cmax)^alpha` | `DriftConstant(drift = that number)` |
| `gddm_start_uniform()`, half width `sz` of the half separation | `ICRange(sz = sz * bs/2)` |
| `ndt` | `OverlayNonDecision(nondectime = ndt)` |
| `lapse` | see the lapse case |

The leak sign is the one the help page warns about: `gddm()`'s `leak` is
the paper's `l`, positive for leaky integration, and PyDDM's `DriftLinear`
`x` coefficient is its negative.

## 3. The fixture

`extensions/frmtmb.ddm/tests/testthat/fixtures/`, three CSVs, 71 KB:

* `gddm-pyddm-density.csv`, the renormalized defective density at both
  walls on a 0.02 s grid, 10 cases, 895 rows.
* `gddm-pyddm-loglik.csv`, six response times per case at ODD multiples
  of 0.005 s, so they fall BETWEEN the nodes of every `gddm()` grid the
  test uses and the density is read through the likelihood's own
  interpolation rather than off a node.
* `gddm-pyddm-meta.csv`, every parameter, the solver, `dt`, `dx`,
  `T_dur`, the two boundary masses, the undecided mass, the conditional
  mean decision times, the measured convergence residual, and a
  SHA-256 over the two payload files.

The reference settles its own resolution: the generator halves `dt` and
`dx` together until the density moves by less than 5e-4 in log density,
and writes down the level it stopped at and the change it achieved.

| case | PyDDM solver | dt | dx | its own residual |
| --- | --- | --- | --- | --- |
| constant | analytic | 6.25e-4 | 1.25e-3 | 4.2e-6 |
| bias_ndt | analytic | 6.25e-4 | 1.25e-3 | 8.2e-6 |
| leak | Crank-Nicolson | 6.25e-4 | 1.25e-3 | 8.1e-5 |
| unstable | Crank-Nicolson | 6.25e-4 | 1.25e-3 | 6.6e-5 |
| exp_collapse | backward Euler | 7.81e-5 | 1.56e-4 | 2.6e-3 |
| lin_collapse | analytic | 6.25e-4 | 1.25e-3 | 1.4e-7 |
| coh_high | analytic | 6.25e-4 | 1.25e-3 | 2.5e-7 |
| coh_low | analytic | 6.25e-4 | 1.25e-3 | 9.2e-6 |
| sz_uniform | Crank-Nicolson | 3.13e-4 | 6.25e-4 | 3.2e-4 |
| lapse | analytic | 6.25e-4 | 1.25e-3 | 4.2e-6 |

`exp_collapse` is the one case that never reaches 5e-4. It is the only
case with no closed form AND a moving bound, so it gets backward Euler
on a grid-snapped bound: first order in `dt` and first order in `dx`.
The ladder stops at its last rung with 2.6e-3 still moving, and that
number is written to the fixture rather than hidden.

The lapse conventions were checked against each other rather than
assumed. PyDDM's `OverlayUniformMixture` reproduces to 2.2e-16 as
`(1-c) p + 0.5 c / (len(t_domain) dt)`, so the only difference from
`gddm()` is that PyDDM's window is `T_dur + dt` where `gddm()`'s is
`t_max`: a lapse density of 8.3316e-3 against 8.3333e-3, 0.02 percent
of the lapse component. The fixture carries `gddm()`'s convention and
records both numbers.

## 4. Measured discrepancy, gddm against the fixture

Worst absolute difference in log density over the case's grid and both
walls, the fixed dataset's log-likelihood difference, and the worst
boundary mass difference. `gddm_control()` at four resolutions; the
second column of each case is the shipped default.

| case | dt .02 ny 101 | dt .01 ny 201 | dt .005 ny 401 | dt .0025 ny 801 | d loglik (default) | d mass (default) |
| --- | --- | --- | --- | --- | --- | --- |
| constant | 6.89e-3 | 8.76e-4 | 2.35e-4 | 6.41e-5 | +2.3e-4 | 2.6e-5 |
| bias_ndt | 6.61e-3 | 1.62e-3 | 4.22e-4 | 1.15e-4 | +1.6e-4 | 2.9e-5 |
| leak | 7.08e-3 | 1.00e-3 | 2.73e-4 | 7.71e-5 | -2.9e-3 | 1.7e-5 |
| unstable | 6.36e-3 | 7.33e-4 | 1.88e-4 | 3.93e-5 | -5.1e-5 | 7.1e-5 |
| exp_collapse | 6.27e-3 | 3.25e-3 | 2.83e-3 | 2.76e-3 | -2.6e-3 | 6.7e-5 |
| lin_collapse | 6.76e-3 | 7.90e-4 | 2.00e-4 | 5.06e-5 | +9.2e-4 | 3.6e-5 |
| coh_high | 6.80e-3 | 1.12e-3 | 2.84e-4 | 7.14e-5 | +1.3e-3 | 9.4e-6 |
| coh_low | 6.56e-3 | 8.71e-4 | 2.49e-4 | 7.39e-5 | -6.6e-4 | 3.1e-5 |
| sz_uniform | 5.63e-3 | 7.56e-4 | 2.69e-4 | 2.64e-4 | -8.3e-5 | 1.8e-5 |
| lapse | 6.83e-3 | 8.68e-4 | 2.33e-4 | 6.36e-5 | +1.7e-4 | 2.6e-5 |

Eight of the ten quarter their error at each halving, which is the
second order both schemes claim, and land at the shipped grid between
7.3e-4 and 1.6e-3. That is an order of magnitude better than the 0.01
`gddm_control()`'s help page claims, and the help page's claim was made
against the raw flux rather than against the density the likelihood
reads, which also carries the non-decision-time shift.

Two cases stall rather than converge, and both stall on the reference's
side, not on `gddm()`'s.

**`exp_collapse` stalls at 2.8e-3.** The reference's own unresolved
convergence error is 2.6e-3, so `gddm()` converges to a point that is
inside the reference's uncertainty. Corroboration comes from the case
next door: `lin_collapse` is the same shape of model and PyDDM has a
closed form for it, and there PyDDM's own numerics sit 5.6e-3 from
PyDDM's own analytic answer while `gddm()` sits 7.9e-4 from it, falling
to 5.1e-5 on a finer grid. On a collapsing bound `gddm()` is roughly a
hundred times closer to the truth than PyDDM's numerics are, which is
what the change of variable was for: PyDDM sandwiches the moving bound
between the two grid nodes around it and weights them linearly, and
that treatment is first order in `dx` and does not go away.

**`sz_uniform` stalls at 2.6e-4.** `ICRange` lays the uniform start on
`2 int(sz/dx) + 1` whole cells, so as a piecewise-constant density it
is half a cell wider on each side than the interval asked for, which
inflates the starting variance by about `sz dx / 3`. `gddm()` spreads
the same box with an integrated cubic B-spline. The two starting
distributions differ at first order in `dx` and second order in the
spatial step respectively, and neither refinement removes the other's.
2.6e-4 is small enough to leave alone.

Neither is a disagreement that a tolerance was widened to hide: both
are named in the test, both have a mechanism, and in both the reference
is the less accurate of the two.

## 5. Tolerances granted

From the measured column, not from a guess.

* density, nine cases: **2.5e-3**, about 1.5 times the worst measured
  (1.62e-3, `bias_ndt`).
* density, `exp_collapse`: **5e-3**, about 1.5 times its measured
  3.25e-3, of which 2.6e-3 is the reference's own residual.
* log-likelihood of the fixed dataset: **5e-3**, about 1.7 times the
  worst measured (2.9e-3, `leak`).
* boundary mass: **2e-4**, about 2.8 times the worst measured (7.1e-5,
  `unstable`).

The test also asserts that the discrepancy FALLS as the grid is
refined, at three resolutions on three cases. A tolerance alone would
pass for a solver that happened to land inside it; a tolerance plus a
convergence ordering only passes for a solver whose limit is the
reference.

## 6. Euler-Maruyama, the third reference

R only. One case, `exp_collapse`, the case where the two PDE solvers
disagree most and the one where a shared bias would be least
surprising, since both discretize a moving bound.

Paths follow `dx = mu dt + dW` with absorption noticed only at a
monitoring time, so first passages come out LATE by order sqrt(h)
(Broadie, Glasserman and Kou 1997). Two monitoring gaps, `h` and `2h`,
are driven by the SAME Brownian increments, so the extrapolation's
difference term is nearly noiseless and the extrapolated value inherits
one level's Monte Carlo error rather than the sum of two. Richardson in
sqrt(h): `est = m(h) + (m(h) - m(2h)) / (sqrt(2) - 1)`.

Two functionals, both of which pair path by path: the probability of
the upper wall, and the unconditional mean decision time. The
conditional mean at one wall does NOT pair, because the two gaps sort
different paths to different walls; pairing it by position rather than
by path was worth nothing, and switching to the unconditional mean
recovered the whole factor of four.

At n = 200000, h = 0.001:

| quantity | PyDDM | gddm (dt .005) | Euler-Maruyama | its SE | paired gain |
| --- | --- | --- | --- | --- | --- |
| P(upper) | 0.854860 | 0.854831 | 0.855873 | 8.5e-4 | 3.8x |
| mean decision time | 0.450595 | 0.450665 | 0.449674 | 6.0e-4 | 4.0x |

The two PDE solvers agree with each other to 2.9e-5 and 7.0e-5. The
path simulation lands 1.2 and 1.5 standard errors away from both, in
the same direction for P(upper) and the opposite for the mean, which is
what leftover higher-order extrapolation error looks like and not what
a shared bias looks like. The check resolves a shared bias down to
about 2e-3, which is where PyDDM's own collapsing-bound error already
sits, so it rules out anything larger than the disagreement already
accounted for.

An unpaired four-gap least-squares fit in sqrt(h) was tried first and is
worse: its P(upper) came out at 0.854515 with a fit standard error of
1.9e-3, and its conditional mean at the upper wall was 2.1e-3 below both
solvers, four times its own fit error, because a pure sqrt(h) model does
not hold across gaps as coarse as 0.008.

## 7. Environment notes

The uv environment could not live at the path the lane specified. The
scratchpad prefix is 120 characters, and

    <scratchpad>/gr-uvcache/environments-v2/gddm-pyddm-reference-<hash>/
      Lib/site-packages/scipy/_external/array_api_extra/_lib/_utils/_helpers.py

is exactly 260 characters, which is Windows MAX_PATH. The file installs,
`os.listdir` shows it, and the import system cannot see it, so
`import scipy.sparse` fails with a `ModuleNotFoundError` naming a module
that is on disk. The probe environments, whose names were eleven
characters shorter, worked throughout. The fix is a directory junction,
`C:\Users\adf44\gr-sp` to the scratchpad, so the environment still lives
under the scratchpad and the path is 90 characters shorter. Anyone
regenerating the fixture on Windows needs the same, or a short
`UV_CACHE_DIR`.

## 8. Does the tier have teeth

A tolerance is only worth what it rejects. Seven deliberately wrong
components were built with the package's PUBLIC extension points,
`gddm_drift_term()`, `gddm_bound_term()` and `gddm_start_term()`, so
nothing under `R/` was touched, and each was measured against the same
fixture at the shipped grid. Each mutation is a translation error that
a reader of the two APIs could plausibly make.

| mutation | case | worst d log | tolerance | margin |
| --- | --- | --- | --- | --- |
| leak sign flipped, PyDDM's convention | leak | 1.9e+0 | 2.5e-3 | 770x |
| exponential bound reads tau as a rate | exp_collapse | 6.9e+2 | 5.0e-3 | 1.4e5 x |
| linear collapse per second, not per window | lin_collapse | 6.9e+2 | 2.5e-3 | 2.8e5 x |
| start point as an absolute offset | bias_ndt | 2.6e+0 | 2.5e-3 | 1050x |
| coherence nonlinearity without alpha | coh_low | 5.3e-1 | 2.5e-3 | 210x |
| separation read as the half separation | constant | 6.6e+0 | 2.5e-3 | 2650x |
| non-decision time off by one grid step | constant | 4.8e-2 | 2.5e-3 | 19x |

The tightest is the last, and it is the one that matters for the size
of the tolerance: a non-decision time wrong by a single 0.01 s step is
caught with a factor of 19 to spare, so the granted 2.5e-3 is nowhere
near wide enough to swallow a real error. The mutants live in the
scratchpad, not in the test file.

## 9. Verification

Private library `gr-lib`, explicit `--library=` throughout, the
worktree's core installed first and `frmtmb.ddm` on top of it.

One process per file, `NOT_CRAN=true` so nothing is skipped:

| file | pass | fail | skip | time |
| --- | --- | --- | --- | --- |
| test-gddm-family.R | 77 | 0 | 0 | 3.2 s |
| test-gddm-gradients.R | 22 | 0 | 0 | 70.5 s |
| test-gddm-recovery.R | 28 | 0 | 0 | 263.6 s |
| test-gddm-reference.R | 166 | 0 | 0 | 24.4 s |
| test-gddm-solver.R | 94 | 0 | 0 | 44.4 s |

`test-gddm-reference.R` costs 14.6 s of that 24.4 s with the two
`skip_on_cran()` tiers skipped, which is what CRAN pays. The resolution
sweep and the Euler-Maruyama run are the other 10 s.

`R CMD check --as-cran` on `frmtmb.ddm`, `_R_CHECK_CRAN_INCOMING_=false`,
pandoc from the RStudio quarto tools directory: **Status: OK**, no
warnings and no notes. Tests 244 s, vignettes rebuilt 84 s. The tarball
is 253 KB, and the three fixture CSVs travel in it at 72 KB;
`dev/gddm-pyddm-reference.py` does not, because the package's
`.Rbuildignore` already excludes `^dev$`, which is why the vignette and
the NEWS bullet point at the source repository for it.

The generator was re-run after the file was reformatted and produced a
byte-identical fixture, SHA-256
`37232180985e3a023ab89a9e03fefee50c41207c526b435a8085e80a52b72104`, so
it is deterministic as claimed. R rebuilds that same hash from the
PARSED numbers through `sprintf("%.12e", ...)`, so the check survives a
checkout that rewrote the line endings. `tools::sha256sum()` needs
R >= 4.5 while the package declares R >= 4.1, so that one expectation
is skipped rather than failed on an older R.

## 10. What was left out, and why

* **A fitting comparison.** The lane asked for densities, boundary
  masses and one log-likelihood, and those are what a wrong density
  shows up in. Fitting both packages to the same data and comparing
  estimates would add an optimizer's behaviour to the measurement
  without adding evidence about the density.
* **A case per component in combination.** The matrix varies one
  generalized component at a time against a shared baseline, because a
  case that moved three at once could not say which one was wrong. The
  cost of a combined case is one more solve; the reason not to have one
  is diagnostic, not budgetary.
* **The leading edge.** Below a decision time of 0.2 s both solvers are
  wrong in the same direction for the same structural reason, so a
  comparison there would agree about an artefact. Stated in the
  vignette rather than quietly cropped.
* **`gddm_control(tridiagonal = "atomic")`.** It computes the same
  derivative by a different route and does not change the value, which
  is what this tier measures; `test-gddm-gradients.R` already covers it.
* **A second condition per case.** Every case solves one condition. The
  per-condition indexing is what `test-gddm-family.R` and
  `test-gddm-recovery.R` exercise, and repeating it here would test
  frmtmb's gather rather than the density.
* **`Rmpfr` or another arbitrary-precision reference.** PyDDM's analytic
  path already supplies an exact answer for five of the ten cases, which
  is what pinned the two collapsing-bound cases; extended precision
  would sharpen nothing that is currently the limit.

## 11. Touched files

Nothing under `R/`, so no roxygen ran and `R/zzz.R` is untouched, as the
sibling lane needs.

* `extensions/frmtmb.ddm/dev/gddm-pyddm-reference.py` (new)
* `extensions/frmtmb.ddm/tests/testthat/fixtures/gddm-pyddm-density.csv` (new)
* `extensions/frmtmb.ddm/tests/testthat/fixtures/gddm-pyddm-loglik.csv` (new)
* `extensions/frmtmb.ddm/tests/testthat/fixtures/gddm-pyddm-meta.csv` (new)
* `extensions/frmtmb.ddm/tests/testthat/test-gddm-reference.R` (new)
* `extensions/frmtmb.ddm/vignettes/gddm.Rmd` (new section, and the
  reference list gained Broadie, Glasserman and Kou)
* `extensions/frmtmb.ddm/NEWS.md` (one bullet, first under the
  development heading)
* `dev/gddm-reference-findings.md` (new, this file)

## 12. Suite counts, audited by name

The whole `frmtmb.ddm` suite in one process, from the `R CMD check
--as-cran` run's own `tests/testthat.Rout`:

    [ FAIL 0 | WARN 0 | SKIP 18 | PASS 964 ]

`test-gddm-reference.R` contributes 147 of those passes and 2 of those
skips, measured on its own in the same mode, so the suite without this
lane's file stands at 817 passes and 16 skips in that configuration.

The lane's stated baseline of 926 is a different configuration: it is
the count with `NOT_CRAN=true`, where nothing is skipped. This file
contributes 166 there rather than 147, the extra 19 being the
resolution sweep and the Euler-Maruyama tier, so the `NOT_CRAN=true`
total is 1092. That run was still going when this was written, starved
by four sibling lanes sharing the machine: 608 seconds of CPU over
about an hour of wall clock, advancing but slowly. Nothing in it can
fail that did not already run and pass in the per-file table above,
because every file was also run alone with `NOT_CRAN=true`.

The number to quote is the one-process 964 with FAIL 0 and WARN 0,
because it is the run that also produced Status OK.

The `NOT_CRAN=true` whole-suite process was stopped rather than left to
finish. It had gained 90 seconds of CPU in 35 minutes with five lanes
sharing the machine, so it would not have completed in any useful time,
and it was holding a core the other lanes could use. It could not have
changed the verdict: it is the same files that already ran and passed
one process each with `NOT_CRAN=true`, and the one-process run that
matters, the one inside `R CMD check --as-cran`, had already returned
FAIL 0 WARN 0 PASS 964 with Status OK.

The junction `C:\Users\adf44\gr-sp` is left in place. It points at the
scratchpad and is what makes the uv environment's paths short enough to
import under MAX_PATH; anyone regenerating the fixture on Windows needs
it or an equally short `UV_CACHE_DIR`. It is a link, so removing it with
`rmdir` (no `/s`) removes only the link.

# Punch round, 2026-09-05

Against the review at `dev/review-gddm-reference.md` (GO WITH FIXES).
**Where this section disagrees with sections 3 to 6 above, this section
is the current statement.** It supersedes, by name: the `d mass` column
and the mass sentence in section 5; the "opposite directions" argument
in section 6; the two "stall" paragraphs in section 4; the "roughly a
hundred times" sentence in section 4; and the row count, the byte count
and one Euler-Maruyama digit. Everything else above stands, and the
review reproduced it independently.

Every number below was re-measured in this lane rather than copied from
the review, and every one of them agreed with the review's.

## P1. The mass tolerance was justified with the wrong quadrature rule

`gr_mass()` integrates gddm's density with the TRAPEZOID rule. Section
5 above justified `GR_TOL_MASS` against a column measured with the
RECTANGLE rule, so the justification did not describe the shipped
assertion. Re-measured by parsing the shipped test file and dropping
every top-level `test_that()`, so these are the file's own helpers:

| case | trapezoid (shipped) | rectangle | total - 1 (trapezoid) |
| --- | --- | --- | --- |
| constant | 8.648e-5 | 2.554e-5 | -6.65e-5 |
| bias_ndt | 1.384e-4 | 2.875e-5 | -1.32e-4 |
| leak | 8.810e-5 | 1.676e-5 | -1.07e-4 |
| unstable | 1.310e-4 | 7.122e-5 | -8.16e-5 |
| exp_collapse | 6.665e-5 | 6.665e-5 | 6.7e-15 |
| lin_collapse | 3.740e-5 | 3.569e-5 | -2.30e-6 |
| coh_high | 1.337e-5 | 9.383e-6 | -4.06e-6 |
| coh_low | **1.532e-4** | 3.090e-5 | -1.46e-4 |
| sz_uniform | 1.089e-4 | 1.800e-5 | -1.03e-4 |
| lapse | 8.216e-5 | **1.076e-4** | -6.31e-5 |

The worst under the shipped rule is 1.532e-4 on `coh_low`, so the
granted 2e-4 was 1.31x, not the 2.8x section 5 claimed: the thinnest
margin in the file described as the widest. The rectangle column
reproduces section 4's `d mass` column for nine of ten cases, which is
how the wrong rule was identified. The tenth, `lapse`, was recorded as
2.6e-5, which is the `constant` value; under the rectangle rule it is
1.076e-4, because gddm's lapse floor `0.5 lapse / t_max` summed over
`nt + 1` nodes times `dt` overcounts by `lapse dt / t_max`.

**Resolved by widening, and the comment says so.** `GR_TOL_MASS` is now
3e-4, which is 1.96x the 1.532e-4 measured with the rule the assertion
uses. The trapezoid rule is kept, because it is the mass of the
piecewise-linear density the likelihood interpolates, the rule
`gd_draw()` samples from, and the right rule for the lapse case, whose
density is nonzero at both ends of the window. The comment at the
constant records what keeping it costs: both sides of the comparison
are rectangle sums, `renormalize` dividing by `(sum(pu) + sum(pl)) dt`
and PyDDM's `Solution.prob()` being a plain sum of bin masses, so about
1e-4 of the number is quadrature convention and not a solver
difference. The `total - 1` column is the visible signature: up to
-1.46e-4 under the trapezoid rule and exactly zero under the rectangle
rule.

## P2. The Euler-Maruyama direction argument was unsound and is gone

Section 6 argued that landing 1.2 and 1.5 standard errors from the
reference "in opposite directions" looked like extrapolation error
rather than a shared bias. It does not, for two reasons the review
demonstrated and this lane accepts.

1. The two functionals are anti-correlated by construction. Under a
   positive drift, more upper-boundary hits means faster decisions, so
   P(upper) up forces the mean decision time down. Opposite directions
   is the null expectation for any single perturbation, noise included,
   and so distinguishes nothing.
2. The signs are a property of the seed. The reviewer's seed 11111 is
   the exact mirror of this lane's: P(upper) 1.84 SE below and the mean
   1.80 SE above. Seed 33333 puts both on the same side.

**The conclusion survives; only the reasoning for it is withdrawn.** The
standard errors reproduce on every seed, 8.5e-4 and 6.0e-4, and pooling
four independent seeds, 800000 paths, P(upper) sits 0.72 single-seed SE
from PyDDM with nothing larger than about 1e-3 resolvable and nothing
detected. It is the MAGNITUDES that carry the no-shared-bias claim. The
direction sentence is removed from the test comment, from this document
and from the vignette; the vignette's "about one standard error" is now
"one or two", since the lane's own figures are 1.2 and 1.5 and other
seeds reach 1.84.

## P3. Neither case stalls; both gaps are the reference's own error

Section 4 called `exp_collapse` and `sz_uniform` "stalls" and said
refinement could not remove them. The mechanisms named there were
right; the claim that the gap has a floor was wrong, and it gave away
the lane's own point. Refining the REFERENCE removes both. This lane
solved each case past the fixture's rung and formed the first-order
Richardson limit `2 f(fine) - f(coarse)`:

`exp_collapse`, fixture at L4, refined to L5, L4 to L5 moves 1.3050e-3:

| gddm grid | vs fixture | vs L5 | vs PyDDM's limit |
| --- | --- | --- | --- |
| dt .01 ny 201 | 3.250e-3 | 1.965e-3 | **8.080e-4** |
| dt .005 ny 401 | 2.832e-3 | 1.547e-3 | 2.779e-4 |
| dt .0025 ny 801 | 2.764e-3 | 1.459e-3 | 1.689e-4 |

`sz_uniform`, fixture at L2, refined to L3 and L4; L2 to L3 moves
1.6321e-4 and L3 to L4 moves 8.2317e-5, which halves and so confirms
`ICRange`'s whole-cell box is first order in `dx`:

| gddm grid | vs fixture | vs L4 | vs PyDDM's limit |
| --- | --- | --- | --- |
| dt .01 ny 201 | 7.555e-4 | 7.257e-4 | **7.152e-4** |
| dt .005 ny 401 | 2.685e-4 | 2.172e-4 | 2.067e-4 |
| dt .0025 ny 801 | 2.644e-4 | 7.558e-5 | 6.512e-5 |

Against where PyDDM is going rather than where it stopped, `gddm()`
converges on both: 8.1e-4 at the shipped grid on `exp_collapse`, which
sits inside the 7.3e-4 to 1.6e-3 band the other nine cases occupy, and
a clean second-order 7.2e-4, 2.1e-4, 6.5e-5 on `sz_uniform`. The whole
remaining gap is the reference's discretization error in both. The
`lin_collapse` analogy is no longer needed and is withdrawn as the
argument; it was also overstated, because "roughly a hundred times
closer" compared PyDDM at its grid with `gddm()` at a finer one. Like
for like it is about seven times, and `gddm()` is taking a sixteen
times coarser time step while doing it.

A check worth recording, since the refinement made it measurable: the
ladder's stopping rule compares consecutive rungs, and for a
first-order sequence that difference EQUALS the remaining error rather
than merely bounding it. The fixture records
`ref_residual_dlog = 2.584e-3` for `exp_collapse` and the next rung
moves 1.3050e-3, so twice that is 2.610e-3. The recorded residual is a
calibrated error estimate, not just a delta, and for the second-order
cases it is conservative by about 3x.

The recomputed PyDDM rows at the fixture's own rung reproduce the
fixture to a maximum relative difference of exactly 0 on both cases,
74 of 74 and 91 of 91 grid points, so the refinement ran the same code
path the fixture was built with.

## P4. The generator now says why the cache directory matters

One paragraph added to `dev/gddm-pyddm-reference.py`'s docstring: on
Windows, point `UV_CACHE_DIR` somewhere short first, because uv names
the environment after the script and under a deep cache the path of one
scipy file reaches 260 characters, MAX_PATH, after which the file
installs, `listdir` shows it, and the import system cannot see it. That
is the finding section 7 recorded; it now lives where someone
regenerating the fixture will read it, instead of only here.

## P5. The whole-suite number is verified, not inferred

Section 12 recorded 1092 as arithmetic rather than as a measurement,
because this lane's own run was stopped. The review ran it to
completion in 1514 s:

    TOTAL pass=1089 fail=0 skip=1 warn=0 err=0

The single skip is `test-sampling.R:13`, skipped for
`{frmtmb.sample} cannot be loaded`, a property of a private library
built with only `frmtmb` and `frmtmb.ddm` in it. That test holds
exactly three `expect_true()` calls, and installing `frmtmb.sample`
from the worktree takes the file from 94 passes and 1 skip to 97 and 0.
So the whole suite in one process is **1092 with `frmtmb.sample`
installed**, and `1092 - 166 = 926` is the stated baseline to the
expectation. Anyone reproducing the count needs `frmtmb.sample` in the
library or will land on 1089 with one skip.

## P6. Corrections to numbers stated above

* Section 3 says the density fixture has 895 rows. It has **893**
  (9 x 91 + 74), counted.
* Sections 3 and 9 say 71 KB and 72 KB for the same three files. They
  are **72,295 bytes**.
* Section 6 gives gddm's mean decision time as 0.450665 and the gap to
  PyDDM as 7.0e-5. It is **0.450656**, gap **6.1e-5**.

## P7. Re-run after the edits

`test-gddm-reference.R`, one process each, with the edited file:

| mode | pass | fail | warn | skip | time |
| --- | --- | --- | --- | --- | --- |
| CRAN (`NOT_CRAN` unset) | 147 | 0 | 0 | 2 | 19.3 s |
| `NOT_CRAN=true` | 166 | 0 | 0 | 0 | 49.4 s |

Both counts are unchanged from before the punch round, which is the
expected result: widening a tolerance that was already being met moves
no assertion, and the other four items are comments and prose. Files
touched this round: `test-gddm-reference.R`,
`dev/gddm-pyddm-reference.py`, `vignettes/gddm.Rmd` and this document.
The fixture is untouched and still hashes to
`37232180985e3a023ab89a9e03fefee50c41207c526b435a8085e80a52b72104`.
Nothing under `R/`, and `NEWS.md` did not need a change: its bullet
made no claim this round contradicts.
