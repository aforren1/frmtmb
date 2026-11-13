# Lane wcheck: the unbounded-weights check (round of 2026-11-11)

Base: e8abed0a (frmtmb 0.81.0), branch wt-r82-wcheck, library
`C:/Users/adf44/source/r/wt-wcheck-lib`. "Before" is rellib-r20.
Brief: `round-briefs/r082-wcheck.md` (items 3447, 3450/3608, 3482,
3619, 3435, 2972; 3466 is quadseam's).

Status: done, including the fix round after the review
(`dev/reviews/2026-11-11-wcheck.md`, verdict FIX); the section "Fix
round" at the end supersedes what it corrects above it. Nothing
committed.

## What changed

`R/weights-unbounded.R` only, in core. No extension changed.

1. **Each response is checked (3447).** `check_weights_unbounded()`
   loops over the responses and calls the new
   `weights_unbounded_response()` for each, on that response's rows (its
   `Z`, weights, `cens()` and, under `subset()`, its own row count). The
   refusal names the response ("in the mu of response ly2") and the
   term as that response writes it (`block_term_in_response()` strips
   the frame's "y2 " label prefix, and takes the response's own
   component of an ID block). `rescor = TRUE` stays out: its joint
   density is not a product over responses, and `frame()` refuses
   `weights()` with `rescor = TRUE`, so no rescor model has an excess
   (unweighted, `sum(w - 1) = 0`; rescor is gaussian or student only, so
   no Laplace-linear family reaches it).
2. **The censored allowance (3450, 3608).** Under Laplace, a censored
   row costs 1/2 only when no constant mean brings its probability to 1
   with the others. When the mean's design spans the constant and there
   is no `trunc()`, the cost is `n_cens - max(n_right, n_left)` halves:
   the cheaper side's rows are free. Otherwise every censored row still
   costs 1/2 (the 0.81.0 rule, a proved upper bound, see the release
   review's section 3 in `dev/reviews/2026-10-21-release.md`). The exact
   route (quadrature, full sampling) is unchanged: censored rows cost 0.
3. **One cell through one coefficient (3435, first two cases).** For a
   family with one spike dpar, each column of its design with a single
   nonzero value on some rows but not all (a treatment dummy of
   `sigma ~ 1 + f` or a cell of `sigma ~ 0 + f`) is a direction that
   moves those rows alone. If that coefficient's class `b` prior is
   free (`prior_column_verdict()`) and the cell's own excess is
   positive, the model is refused, and the message names the
   coefficient. `prior_cellmeans_verdict()` also returns "free" when
   every column is free (the review's fix 1).
4. **Power-law tails (3435, third case).** `unbounded_spike_rate` gives
   the width's power per spike dpar (1/2 for Gamma, Beta,
   inverse.gaussian and von_mises; 1 for weibull, frechet and the
   scales). The prior verdict gets `rate * excess`. On the natural
   scale toward infinity, student_t (`nu + 1`), inv_gamma (`shape + 1`)
   and flat (0) tails are compared with it (MAP: bounds above, free
   below, unknown at equality; sample: bounds when `p - excess > 1`).
   On the log link, exponential, gamma (by `rate`) and logistic (by
   `1 / scale`) are decided the same way where they were "unknown".

The message gains `resp`, `coef` and `cost`; the censored clause now
says how many rows take 1/2 ("the 8 that a constant mean cannot bring
to 1 with the rest"), and the power reads `sum(w) - n - 8/2 = 4`.

## Items

| item | verdict | evidence |
|---|---|---|
| 3447 | FIXED | part mv below: 8 runaways refused, 11 absent cases unchanged |
| 3450, 3608 | FIXED for one-sided censoring | part cens; path slopes below |
| 3435 | FIXED, all three cases | part prior |
| 2972 | neither half reproduces on rellib-r20; closed | section 2972 |
| 3482 | not this rule; boundfam's | excess 0 is bounded and not refused |
| 3619 | not this rule; filed analysis | section 3619 |
| 3466 | quadseam's; partly caught here | mv quadrature now refused |

### Verdict table (generated)

`Rscript dev/wcheck-cases.R <lib|base> <mv|cens|prior>` (logs
`dev/wcheck-out/cases-<part>-<base|lane>.txt`), table by
`Rscript dev/wcheck-table.R` (`dev/wcheck-out/table.txt`), pasted
verbatim. Seeds: part mv and prior `set.seed(20261021)`, n = 80, data
of `dev/fuzzfix-rev-wu.R` and `dev/fuzzfix-rev-r2-edges.R`; part cens
`set.seed(20261021)`, n = 80, data of `dev/rel073-wcens-probe.R`.
"betad" is the dispersion coefficients (log sigma, log shape, ...).

Part mv (19 arms)

| arm | rellib-r20 | lane |
|---|---|---|
| mv gaussian x2 (1 \| obs), w | fit 836.78, betad -37.49 -21.39, conv 1 | refused |
| mv gaussian + poisson, w | fit 385.491, betad -37.06, conv 1 | refused |
| mv poisson + Gamma, w | fit 15491.5, betad 38.08, conv 1 | refused |
| mv gaussian x2 (1 \| p \| obs), w | fit 1100.18, betad -36.85 -35.38, conv 1 | refused |
| mv gaussian x2, w on y2 only | fit 320.025, betad 0.2171 -37.67, conv 1 | refused |
| mv gaussian, subset(sub) on y, w | fit 235.698, betad -36.51 0.4569, conv 1 | refused |
| mv gaussian x2, lognormal prior sigma resp y | fit 858.328, betad -23.61 -36.78, conv 1 | refused |
| mv gaussian x2 (1 \| obs), w, quadrature | fit -179.772, betad -3.782 -3.537, conv 0 | refused |
| mv gaussian x2 (1 \| obs), w7 (sum < n) | fit -241.407, betad 0.1548 0.4356, conv 0 | fit -241.407, betad 0.1548 0.4356, conv 0 |
| mv gaussian x2 (1 \| g), w | fit -343.211, betad 0.1068 0.4061, conv 0 | fit -343.211, betad 0.1068 0.4061, conv 0 |
| mv gaussian x2, lognormal prior sigma both | fit 460.073, betad -19.8 -19.8, conv 0 | fit 460.073, betad -19.8 -19.8, conv 0 |
| poisson (1 \| obs), overdispersed, unweighted | fit -172.495, betad , conv 0 | fit -172.495, betad , conv 0 |
| poisson (1 \| obs), overdispersed, w | fit -202.307, betad , conv 0 | fit -202.307, betad , conv 0 |
| mv poisson x2 (1 \| obs), overdispersed | fit -299.43, betad , conv 0 | fit -299.43, betad , conv 0 |
| mv poisson x2 (1 \| obs), overdispersed, w | fit -356.07, betad , conv 0 | fit -356.07, betad , conv 0 |
| mv poisson x2 (1 \| p \| obs), overdispersed | fit -299.325, betad , conv 0 | fit -299.325, betad , conv 0 |
| mv gaussian + poisson (1 \| obs), unweighted | fit -297.527, betad 0.144, conv 0 | fit -297.527, betad 0.144, conv 0 |
| mv gaussian x2 (1 \| obs), unweighted | fit -274.795, betad -0.1678 0.244, conv 0 | fit -274.795, betad -0.1678 0.244, conv 0 |
| mv gaussian x2 (1 \| obs), w7, quadrature | fit -241.407, betad 0.1548 0.4356, conv 0 | fit -241.407, betad 0.1548 0.4356, conv 0 |

The "lognormal prior sigma both" arm is bounded (a prior on each
response's sigma) and is not refused; its estimate is the degenerate
balance the message describes (log sigma -19.8). The weighted and
unweighted overdispersed poisson arms, one and two responses, are the
brief's identified observation-level effects: they fit with conv 0 and
no warning on both builds.

Part cens (19 arms)

| arm | rellib-r20 | lane |
|---|---|---|
| gaussian right E 0 Laplace | fit -107.111, betad 0.02634, conv 0 | fit -107.111, betad 0.02634, conv 0 |
| gaussian right E 0.5 Laplace | fit -107.731, betad -0.2226, conv 0 | refused |
| gaussian right E 2 Laplace | fit -108.975, betad -0.8617, conv 0 | refused |
| gaussian right E 4 Laplace | fit -108.949, betad -1.923, conv 0 | refused |
| gaussian right E 7.5 Laplace | fit -72.2017, betad -14.16, conv 1 | refused |
| gaussian right E 8 Laplace | fit -48.2781, betad -13.64, conv 1 | refused |
| gaussian right E 16 Laplace | refused | refused |
| gaussian right E 4 quadrature | refused | refused |
| gaussian right E 4 (1 \| g4) | fit -112.473, betad 0.01757, conv 0 | fit -112.473, betad 0.01757, conv 0 |
| gaussian both E 2 Laplace | fit -104.885, betad -0.5761, conv 0 | fit -104.885, betad -0.5761, conv 0 |
| gaussian both E 3.5 Laplace | fit -105.744, betad -0.985, conv 0 | fit -105.744, betad -0.985, conv 0 |
| gaussian both E 4.5 Laplace | fit -105.987, betad -1.264, conv 0 | refused |
| gaussian both E 8 Laplace | fit -104.262, betad -2.336, conv 0 | refused |
| asym_laplace  0 of 80 right censored | refused | refused |
| asym_laplace  8 of 80 right censored | refused | refused |
| asym_laplace 40 of 80 right censored | refused | refused |
| asym_laplace 72 of 80 right censored | fit 36.5069, betad -7.241 -2.137, conv 1 | refused |
| asym_laplace 80 of 80 right censored | fit -1.42109e-14, betad 0 -37.24, conv 1 | fit -1.42109e-14, betad 0 -37.24, conv 1 |
| asym_laplace 72 right censored (1 \| g4) | fit -0, betad 0 -38.61, conv 1 | fit -0, betad 0 -38.61, conv 1 |

"right": 16 of 80 rows right censored at the 80th percentile, the 64
uncensored rows at weight 1 + E / 64. "both": also the 8 rows below the
10th percentile left censored, 56 uncensored rows at 1 + E / 56; the 8
left rows cost 1/2 each, so E above 4 is refused and E up to 4 is not
(not proved bounded either: the check fails open there).

Part prior (20 arms; excess 23.868, cell a 11.587, cell b 12.281)

| arm | rellib-r20 | lane |
|---|---|---|
| cellmeans, no prior | refused | refused |
| cellmeans, normal(0,1) on b | fit 153.938, betad -11.59 -12.28, conv 0 | fit 153.938, betad -11.59 -12.28, conv 0 |
| cellmeans, normal on faca only | fit 377.179, betad -5.12 -36.89, conv 1 | refused |
| cellmeans, normal on facb only | fit 469.509, betad -36.41 -14.67, conv 1 | refused |
| cellmeans, student_t on b | fit 346.267, betad -7.632 -37.54, conv 1 | refused |
| cellmeans, normal faca + student_t facb | fit 352.358, betad -5.122 -37.62, conv 1 | refused |
| fac + x, normal(0,1) on b | fit 157.529, betad -11.59 -12.28 -1.895, conv 0 | fit 157.529, betad -11.59 -12.28 -1.895, conv 0 |
| fac + x, normal on faca + facb | fit 179.043, betad -11.89 -12.94 -6.708, conv 0 | fit 179.043, betad -11.89 -12.94 -6.708, conv 0 |
| 1 + fac, normal on Intercept only | fit 437.856, betad -11.89 -23.24, conv 1 | refused |
| 1 + fac, normal on Intercept and b | fit 438.651, betad -23.69 -0.347, conv 0 | fit 438.651, betad -23.69 -0.347, conv 0 |
| Gamma, student_t(3,0,2.5) shape | fit 1.77265e+10, betad 51.68, conv 1 | refused |
| Gamma, student_t(30,0,2.5) shape | fit -136.548, betad 2.484, conv 0 | fit -136.548, betad 2.484, conv 0 |
| Gamma, inv_gamma(2,1) shape | fit 2677.62, betad 37.23, conv 1 | refused |
| Gamma, gamma(0.01,0.01) shape | fit -81.0358, betad 7.086, conv 0 | fit -81.0358, betad 7.086, conv 0 |
| Beta, student_t(3,0,2.5) phi | fit 3.73091e+12, betad 57, conv 1 | refused |
| Beta, gamma(0.01,0.01) phi | fit 85.3553, betad 6.963, conv 0 | fit 85.3553, betad 6.963, conv 0 |
| Gamma cellmeans, normal(0,1) on b | fit -94.2093, betad 5.812 6.145, conv 0 | fit -94.2093, betad 5.812 6.145, conv 0 |
| Gamma cellmeans, normal on faca only | fit 106040, betad 13.7 40.06, conv 1 | refused |
| Gamma cellmeans, gamma(2,1) on b | error: The optimizer failed on this model (Gamm | refused |
| Gamma cellmeans, gamma(2,1) on b, unweighted | error: The optimizer failed on this model (Gamm | error: The optimizer failed on this model (Gamm |

Every arm that changes ran away on rellib-r20 (conv 1, a log dispersion
beyond 13 in size, or a log-likelihood of 1e4 to 1e12), except three
kinds: the censored arms that stopped inside (see the path below), the
mv quadrature arm (conv 0 at the spurious maximum of 3466), and the
3486 arm (NA/NaN gradient). WITHDRAWN in the fix round: "No arm whose
model is bounded is refused" was false. The arm "1 + fac, normal on
Intercept only" is bounded (the review's finding 1: the prior sits on
the centered intercept, which the contrast moves); its refusal was
wrong, and it fits again after the fix.
Student_t(30) on the shape: `nu + 1 = 31` against `23.868 / 2 = 11.9`,
so bounded, and it fits as before. The last two arms are backlog 3486
("gamma(2, 1) on a log-scale b of shape ~ 0 + fac"): with weights it
is now refused, correctly (the gamma tail `exp(-b)` against a growth
of `exp(6.1 b)` per cell); unweighted it still ends at "NA/NaN
gradient", so 3486 stays open.

### The censored rows: the objective along a constant mean

`Rscript dev/wcheck-cens-path.R base` (log
`dev/wcheck-out/cens-path-base.txt`): `dry_run = "objective"` (no
check), mean slope 0 and intercept 2 above the bound, the effects' log
sd re-optimized at each step, log sigma from -2 to -30. The slope per
unit of `-log sigma`:

- gaussian, 16 right censored, uncensored excess 0.5, 2, 4, 7.5: 0.5,
  2, 4 and 7.5 at every step (9 steps each); E = 0, the control: at
  most 3.1e-12 in size to log sigma -18, then -1.03e-08 and -3.06e-05.
- gaussian, 16 right and 8 left, mean above: E = 3.5: -0.62 to -1.00
  (not refused); E = 4.5: 0.20 to 0.38 to log sigma -14, then -0.41
  from -14 to -18 (refused); E = 8: 3.70 to 3.88 to -18 (refused). The
  8 left rows cost about 1/2 each, as the allowance says, less a term
  that grows like log log (1 / sigma) (the release review's section 3);
  at E = 4.5 the margin is 0.5, which that term eats at this range, so
  this arm's refusal rests on the proof, not on the measured slope.
  Past log sigma -18 the inner optimizer jumps (29.3, 24.9 per unit),
  which I read as numerical, not as a value.
- asym_laplace, 72 of 80 right censored, quantile held at 0.5: 8.00 per
  unit from log sigma -4 on, with the mean 2 above the bound and 2
  below it alike; the effects' log sd ran to 10.1 and 12.8. With the
  quantile free, the local re-optimization drifted the quantile and the
  path flattened at -fn 19.18; that arm is not a bound.

So the 0.81.0 refusals that the allowance missed (E 0.5 to 8 with 16
right censored, asym_laplace 72 of 80) are unbounded, which the
construction proves (each right-censored row's probability is 1 at its
effect's prior mean, so its effect has no data curvature) and the
slopes show. That the asym_laplace slope is 8 with the mean BELOW the
bound too is not explained by the derivation at a fixed sd; it is filed
(the allowance there errs open, the safe side).

### Tests

New blocks in `tests/testthat/test-weights-unbounded.R`, seen to fail
on rellib-r20 (`dev/release/run-tests.R`, log
`dev/wcheck-out/tf-wu-base.txt`) and pass on the lane
(`tf-wu-lane.txt`):

    base: RESULT test-weights-unbounded.R pass=83 fail=3 err=3 skip=0 warn=7
    lane: RESULT test-weights-unbounded.R pass=119 fail=0 err=0 skip=0 warn=0

The base failures are behavioral, not missing symbols: each new
`expect_match(err, wu_refusal)` receives a `frmtmb_fit` (the model was
fitted, not refused), the student_t(3) and inv_gamma(2, 1) Gamma fits
did not throw, and the sample-route check did not throw; the 7 warnings
are the runaways' convergence warnings. The censored block's path
assertion is a ratio to the excess (`abs(slope - E) < 1e-3 * E`) with a
unit-weight control on the same path.

Rounding (`Rscript dev/wcheck-perturb.R <lib>`, log
`dev/wcheck-out/perturb-lane.txt`): the five guard-absent fits of the
new blocks, and two of their refusals, on the test data with x and the
response scaled by `1 + 2e-15 * k`, k = 0..24:

    mv_pois        code 0, 0 other warnings, 0 messages x25
    mv_pois_w      code 0, 0 other warnings, 0 messages x25
    mv_gauss_w07   code 0, 0 other warnings, 0 messages x25
    cells_xs       code 0, 0 other warnings, 0 messages x25
    gamma_t30      code 0, 0 other warnings, 0 messages x25
    mv_refused     REFUSED x25
    cell_refused   REFUSED x25

### 2972

Two halves, rerun on rellib-r20 and the lane:

- Weights: `Rscript dev/wcheck-2972-fuzz.R <lib|base>` fits the 7
  specs of `fuzz_plan(20260901, 300)` with `re = ar1` and `aterm =
  weights`, with the halved weights the tier now draws and with the
  pre-0.73.0 draw (`w * 2`, sum(w) - n from 15.7 to 28.9). All six
  table families (gaussian, student, weibull, Gamma, lognormal, Beta)
  are refused at `w * 2` on both builds; zero_inflated_poisson is not
  in the table and fits at conv 0 (-288.435). The -1e15 objectives were
  these runaways, refused since 0.73.0.
- huber: `Rscript dev/wcheck-2972-huber.R base` (the battery's hir2
  design, 400 rows, seeds 1 to 10): every fit ends between -621.7 and
  -564.4, within 0.31 of the fit from the true values; none is near
  +1.39e5. The degenerate stop was 0.70.0's start, replaced at 0.71.0.
  8 of 10 default fits still stop at code 1 (huber's kinked objective,
  warned): cluster 6 of the triage, not this lane's.

### 3482 (boundfam's)

Weights summing to exactly n: the excess is 0 (tolerance `1e-8 *
rows`), the likelihood is bounded, and this check does not and should
not refuse it; arm "gaussian right E 0" fits at conv 0 on both builds.
The boundary message the item asks for belongs to the edge-verdict
rule of lane boundfam.

### 3619: extension families (not fixed)

`Rscript dev/rel073-rev-fp.R C:/Users/adf44/source/r/rellib-r20 ext`
(log `dev/wcheck-out/3619-ext-base.txt`): frmtmb.eam 0.13.1, wiener
`cens()` with weights 1.5 and `(1 | obs)`: logLik 1782.5521, conv 1,
betad 19.565 1.731 0.000; rdm(2): 766.2075, conv 1, betad -2.922 9.279
3.864. Still runaways. It is not this rule. For wiener the spike is the
boundary separation growing while each row's drift effect grows with it
(`v_i = a / (t_i - ndt)`, so the first-passage sd `a^(1/2) v^(-3/2)`
shrinks like `1 / a`), so condition (iii) fails: the effects do not
stay at a finite point, the growth needs the effects' sd to grow with
`a`, and a prior on that sd can bound it. Each family needs its own
proof and a seam through which frmtmb.eam states it (core-boundary
rule: the seam in core, the entries in the extension). I recommend
that as its own item; the backlog row stays open.

### 3466 (quadseam's)

With several responses the check now runs under `quadrature = TRUE`
too (the exact route): "mv gaussian x2 (1 | obs), w, quadrature"
stopped at conv 0 at log sigma -3.78 and -3.54 on rellib-r20 (the
spurious maximum of 3466) and is refused on the lane. What stays for
quadseam: an unbounded quadrature fit the check does not reach (`se()`,
autocorrelation, `mi()`, a family outside the table, a structure that
does not state its rank) still gets the spurious maximum. The check
reads quadrature as exact; it does not detect the rule's own stop.

## Defects found, not fixed (filed under "Filed in the round of
2026-11-11", "From lane wcheck")

- The Laplace allowance still takes 1/2 per censored row on the other
  side, which errs open; asym_laplace's measured cost is 0 there (slope
  8.00), unexplained.
- The observation-level warning names an mvbf term with its response
  ("`y 1 | obs` gives every observation ..."), in `check_olre` (not
  this file).
- asym_laplace, 72 of 80 right censored, `(1 | g4)`: logLik -0 at log
  sigma -38.61, conv 1, after "Starting value 0 for sigma is -Inf
  through its log link", on both builds.
- 3486 unweighted: unchanged (NA/NaN gradient).

## Not done, and why

- 3619: a different rule per extension family (above).
- 3482: boundfam's.
- Two-sided censoring with a linear program over the mean's design:
  the one-sided rule is what the filed cases need and is provable from
  the design spanning the constant alone; a mixed design would need the
  location map of each family and link.
- Cells for a family with two spike dpars (exgaussian): both dpars
  would have to reach their limit on the same rows; the check keeps the
  constant direction there.

## Files touched

- `R/weights-unbounded.R` (the whole change).
- `tests/testthat/test-weights-unbounded.R` (four new blocks and the
  helper `wu_fit_all()`).
- `NEWS.md` (a development-version section, four BREAKING bullets).
- `dev/test-backlog.md` (CLOSED lines under 2972, 3435, 3447, 3450,
  3608; a new section "Filed in the round of 2026-11-11").
- New scripts: `dev/wcheck-cases.R`, `dev/wcheck-cens-path.R`,
  `dev/wcheck-table.R`, `dev/wcheck-perturb.R`, `dev/wcheck-messages.R`,
  `dev/wcheck-probe-frame.R`, `dev/wcheck-2972-fuzz.R`,
  `dev/wcheck-2972-huber.R`, `dev/wcheck-install.sh`,
  `dev/wcheck-check.sh`.

No compat row, ledger row or bin-2 row is touched (the change is a
refusal on weighted models with an observation-level effect; the
ported-suite verdicts are regenerated at consolidation anyway).
Version: minor (behavior changes marked BREAKING); no new dependency,
same floors.

## Suites

One file per process, every suite, 3 at a time (`bash
dev/release/run-par.sh wcheck-suite-par dev/wcheck-out/jobs-all.txt 3`
with `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true` and `FRMTMB_RELLIB`
set to the lane library; log `dev/release/wcheck-suite-par.log`,
"WCHECK-SUITE-PAR ran 510 of 510"). Summed from its RESULT lines:

    frmtmb files=317 pass=26581 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1854 fail=0 err=0 skip=3 warn=0
    frmtmb.latent files=12 pass=687 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=545 fail=0 err=0 skip=2 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=87 pass=4414 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial files=8 pass=561 fail=0 err=0 skip=0 warn=0
    frmtmb.survival files=15 pass=556 fail=0 err=0 skip=0 warn=0

All 510 `lib:` lines load frmtmb from `wt-wcheck-lib`. The 29 skips are
gates not set here: 15 "set FRMTMB_SCALE_TESTS=true", 13 "set
FRMTMB_DRMTMB_FIT_TESTS=true", 1 "set FRMTMB_FUZZ=true".

OpenBLAS 0.3.26 with LAPACK (`bash dev/ciharden-openblas.sh 0.3.26
lapack`; its own check printed LAPACK 3.12.0):
`test-weights-unbounded.R` pass=119 fail=0 err=0 skip=0 warn=0
(`dev/wcheck-out/tf-wu-openblas.txt`). The change is a refusal before
any tape, and its one numeric assertion (the path slope) is a ratio to
the excess, so 0.3.32 was not run.

One session per suite, three orders (`bash dev/release/run-sessions.sh
wcheck-suite-sess 2`, same environment; log
`dev/release/wcheck-suite-sess.log`, "WCHECK-SUITE-SESS ran 27 of 27"):
all 27 sessions fail=0 err=0 warn=0, 194 LEAK lines, every one tagged
"[with load of ...]" (none untagged, none from
test-weights-unbounded.R). Core: forward pass=26581 skip=15, reverse
26746 skip=15, preload 26849 skip=17 (the two extra skips are the
"names frmtmb.survival / frmtmb.eam when it is not loaded" tests, by
design). Every extension's three sessions give the same counts as its
one-file-per-process run above. The `lib:` lines load frmtmb from
`wt-wcheck-lib` and each extension from rellib-r20.

R CMD check --as-cran, core (`bash dev/wcheck-check.sh`, built with
vignettes in `dev/wcheck-check/`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`,
pandoc and TinyTeX on PATH, `R_LIBS` = lane library, rellib-r20, user
library; logs kept as `dev/wcheck-out/check-frmtmb.txt`, the check tree
deleted):

    * checking examples ... [150s] NOTE
    frm_hazard_reads 3.39   0.11    6.47
    * checking HTML version of manual ... [37s] NOTE
    Skipping checking math rendering: package 'V8' unavailable
      Running 'testthat.R' [77m]
     [77m] OK
    Status: 2 NOTEs

The V8 NOTE is the expected one. The examples NOTE is elapsed time
under load (3.5 s of CPU, 6.47 s elapsed, while two one-session suites
ran beside it) on an example this lane does not touch; lane-rules.md
"An examples-timing NOTE on this box measures LOAD". The first attempt
failed at `R CMD build` ("copying to build directory failed") because
a half-deleted `dev/ciharden-out/` was in the tree; it was rerun after
the deletion finished. No extension changed, so no extension check.

## Fix round (after `dev/reviews/2026-11-11-wcheck.md`, verdict FIX)

The rule: a refusal must be a proof. Where the check cannot prove the
objective unbounded it says "unknown" and the fit runs as on 0.81.0.
Each finding was reproduced on the first lane build before the fix
(logs `dev/wcheck-out/fr/*-lane0.txt`), then rerun on the fixed build
(`dev/wcheck-out/fr/*-lane.txt`). Base numbers are the reviewer's logs
(`dev/reviews/wcheck-review-*-base.txt`).

### What changed in `R/weights-unbounded.R`

1. **Centered intercept (BLOCKER 1).** `prior_column_verdict()` returns
   "unknown" when the dispersion design has an intercept, is centered
   (`lp$center` not FALSE), and an `Intercept` entry's verdict is not
   "free". frmtmb evaluates that prior at `b0 + means_X'b`, which the
   contrast moves too. With center = FALSE the intercept is a class `b`
   coefficient, and `prior_spike_verdict()` now reads its own column's
   `b` entry (`prior_b_column_verdict()`) instead of "no Intercept
   entry, free".
2. **resp spelling (BLOCKER 2).** `prior_resp_matches()` accepts no
   `resp`, the raw name, or `brms_stan_name(r)`; every prior lookup of
   the check uses it.
3. **The shift needs the prior's support (BLOCKER 3).**
   `unbounded_shift_sides()` allows the mean up (down) only when the
   design spans the constant and every prior on what moves (the
   centered `Intercept` entry, or every class `b` entry of a design
   without one) has support unbounded on that side
   (`prior_support_sides()`: normal, student_t, logistic and flat both
   sides; exponential, gamma, inv_gamma and lognormal up only; any
   other density, or a finite `lb` or `ub`, blocks that side).
   Otherwise the 0.81.0 costing (1/2 per censored row) applies.
4. **The Jacobian (MAJOR 4).** `prior_entry_verdict()` now uses
   frmtmb's MAP objective for a natural-scale prior, which carries
   `+ log x`: toward 0 the power is `order + 1 - excess`, toward
   infinity `excess + 1 - p`. Bounded when the margin is positive,
   free when negative; at 0 the MAP route is "unknown" and the sample
   route "free". This is the condition the sample route already used,
   so both routes now read one margin.
5. **lognormal on a log-scale coefficient (MINOR 5)** toward +Inf is
   "free" (tail `exp(-(log eta)^2 / 2)` against `exp(excess eta)`).

### Each finding, reproduced and fixed

| finding | case | rellib-r20 | first lane | fixed lane |
|---|---|---|---|---|
| 1 | R1: `sigma ~ 1 + fac`, normal(0, 1) on Intercept, w 0.8 / 1.05 | logLik -94.2889, conv 0 | refused | logLik -94.2889, conv 0 |
| 1 | R2: the same, no prior (slope 2.000 = E_b) | -94.2888, conv 0 (local) | refused | refused |
| 1 | R3: w 0.3 / 1.5, with the prior (unbounded) | 612.7526, conv 1 | refused | 612.7526, conv 1 (fails open, as 0.81.0) |
| 2 | mvbf `y_a` / `y.b`, lognormal on sigma, resp "ya" / "yb" | 835.458, conv 0 | refused | 835.458, conv 0 |
| 2 | the same, resp "y_a" / "y.b" | 835.458, conv 0 | 835.458, conv 0 | 835.458, conv 0 |
| 3 | uniform(-5, 0.5) on the mean's Intercept, 16 right censored, E 2 | -112.216, conv 0 | refused | -112.216, conv 0 |
| 3 | the same, E 6 | -111.436, conv 0 | refused | -111.436, conv 0 |
| 4 | gaussian E 0.5: normal(0, 1), student_t(3, 0, 2.5), exponential(1) on sigma | refused x3 | refused x3 | -112.358, -112.523, -112.184, conv 0 |
| 4 | gaussian E 0.9: the same three | refused x3 | refused x3 | -112.061, -112.767, -111.024, conv 0 |

Rows of finding 1 from `dev/reviews/wcheck-review-center2.R`, 2 from
`-respname.R`, 3 from `-cens2.R`, 4 from `-misc.R jac`; seeds 20261111
(center, cens2, jac) and 20261021 (respname), n = 80.

### Both sides of every boundary of the corrected tail rule

`Rscript dev/wcheck-fr-boundary.R <lib>` (log
`dev/wcheck-out/fr/boundary-lane.txt`; data of
`dev/reviews/wcheck-review-rates.R`, seed 20261111, n = 80). Each arm:
the check's verdict, then the slope of the objective's profile along
the log dispersion (the rest re-optimized) across the grid, against
the predicted `growth + 1 - p` (per unit of `-log sigma`: `excess -
(order + 1)`):

| arm | check | slopes | predicted |
|---|---|---|---|
| Gamma E 12 (growth 6), student_t(5.5) on shape | REFUSED | 0.567 0.508 0.501 0.500 0.500 0.500 | 0.5 |
| Gamma E 12, student_t(6.5) | fit | -0.426 -0.492 -0.499 -0.500 -0.500 -0.500 | -0.5 |
| Gamma E 12, inv_gamma(5.5, 1) | REFUSED | 0.557 0.509 0.501 0.500 0.500 0.500 | 0.5 |
| Gamma E 12, inv_gamma(6.5, 1) | fit | -0.443 -0.491 -0.499 -0.500 -0.500 -0.500 | -0.5 |
| weibull E 6 (growth 6), student_t(5.5) | REFUSED | 0.530 0.501 0.500 0.500 0.500 | 0.5 |
| weibull E 6, student_t(6.5) | fit | -0.424 -0.498, then -5.0e9 (numeric breakdown past log shape 7) | -0.5 |
| gaussian E 0.5, normal(0, 1) on sigma | fit | -0.496 -0.500 -0.500 -0.500 -0.500 | -0.5 |
| gaussian E 1.5, normal(0, 1) on sigma | REFUSED | 0.504 0.500 0.500 0.500 0.500 | 0.5 |
| gaussian E 2, gamma(1.5, 1) on sigma | REFUSED | 0.559 0.508 0.501 0.500 0.500 | 0.5 |
| gaussian E 2, gamma(2.5, 1) on sigma | fit | -0.441 -0.492 -0.499 -0.500 -0.500 | -0.5 |
| Gamma `shape ~ 1` E 12, lognormal(0, 1) on Intercept | REFUSED | 5.52 5.59 5.65 5.69 5.73 5.75 | toward 6 |

The log-link rules by rate (exponential, gamma, logistic) are
unchanged; the review measured both sides of each
(`dev/reviews/wcheck-review-rates-prior-lane.txt`).

### Filed, not fixed (MINOR 6 and 7, and the open part of 1)

Under "Filed in the round of 2026-11-11 / From lane wcheck" in
`dev/test-backlog.md`: the asym_laplace other-side case with the
review's numbers (36 left, 36 right, 8 uncensored: let through, grows
at 8.00 per unit, warns at conv 1); hurdle_lognormal and
zero_inflated_beta (logLik 968.909 and 2.79e10, unrefused; each needs
a proof per family); and the undecided centered-intercept cell
direction (R3 is unbounded and now let through, as on 0.81.0). MINOR
8: NEWS bullet 2 now says "spans the constant (an intercept, or a full
set of cell means)" and names the prior condition. Bullets 3 and 4 and
a new fifth bullet (the Jacobian false refusals, now fitted) are
rewritten; the CLOSED line under 3435 now reads "decided, with
Jacobian".

### The 58-case table and the 24 guard-absent models

`dev/wcheck-cases.R` (all three parts) rerun on the fixed build;
`Rscript dev/wcheck-table.R` > `dev/wcheck-out/fr/table.txt`. Against
the first lane's table exactly one row changes: "1 + fac, normal on
Intercept only" goes from refused to fit 437.856, betad -11.89 -23.24,
conv 1, the rellib-r20 outcome (bounded per the review's finding 1,
with the maximum at a centered intercept near -24.6). The other 57
rows are identical.

The reviewer's 24 guard-absent models (`wcheck-review-common.R`, log
`dev/wcheck-out/fr/common-lane.txt`) match the reviewer's base log line
for line apart from the `lib:` line (`diff` empty).

Rounding, 25 draws each (`dev/wcheck-perturb.R` rerun, and
`dev/wcheck-fr-perturb.R` for the new arms; logs in
`dev/wcheck-out/fr/`):

    mv_pois        code 0, 0 other warnings, 0 messages x25
    mv_pois_w      code 0, 0 other warnings, 0 messages x25
    mv_gauss_w07   code 0, 0 other warnings, 0 messages x25
    cells_xs       code 0, 0 other warnings, 0 messages x25
    gamma_t30      code 0, 0 other warnings, 0 messages x25
    mv_refused     REFUSED x25
    cell_refused   REFUSED x25
    center_R1      code 0, 0 other warnings, 0 messages x25
    center_R2      REFUSED x25
    cens_unif_E2   code 0, 1 other warnings, 0 messages x25
    cens_unif_E6   code 0, 1 other warnings, 0 messages x25

The one warning of the bounded-mean fits is "Standard errors are not
available for 1 of 3 parameters. (Intercept): a bound holds it", the
intercept at its `ub` of 0.5; the test asserts only conv 0 and a finite
logLik.

### Tests

`tests/testthat/test-weights-unbounded.R`:
- the cell block drops the `sigma ~ 1 + fac` refusal under an
  Intercept prior; it now refuses R2 (no prior) and fits R1 at conv 0;
- the sampler block: at an excess of 0.5 under student_t on sigma the
  MAP check returns NULL too (it refused); at 1.8 it refuses;
- the tail block: guard `E / 2 < 30`; student_t at `nu = E/2 - 0.5`
  refused and at `E/2 + 0.5` fitted;
- new blocks: "a prior on sigma counts its Jacobian" (normal,
  student_t, exponential at E 0.5 fit and at 1.5 refused; gamma(2.5,
  1) fit and gamma(1.5, 1) refused at E 2); "a prior that bounds the
  mean keeps censored rows' cost" (the reviewer's uniform case at E 2
  and 6, and its no-prior refusal); "a prior names a response in
  brms's spelling too" (`y_a` / `y.b`, resp "ya" / "yb" and raw).

    rellib-r20:  RESULT test-weights-unbounded.R pass=93 fail=11 err=4 skip=0 warn=10
    fixed lane:  RESULT test-weights-unbounded.R pass=143 fail=0 err=0 skip=0 warn=0

(`dev/wcheck-out/fr/tf-wu-base.txt`, `tf-wu-lane.txt`.) On rellib-r20
the new refusals fail as behavior (the student_t at `E/2 - 0.5` and
the earlier blocks get a fit, not a refusal), and the new guard-absent
arms of the Jacobian block fail too (rellib-r20 refuses them).

### Runs

The fixed build, installed into `wt-wcheck-lib` (only
`R/weights-unbounded.R` changed in code; no extension changed).

- Core, one file per process, 3 at a time, `NOT_CRAN=true`,
  `FRMTMB_BRMS_FIT_TESTS=true` (`dev/release/wcheck-suite-fr.log`,
  "WCHECK-SUITE-FR ran 317 of 317"; all 317 `lib:` lines name
  `wt-wcheck-lib`): pass=26605 fail=0 err=0 skip=15 warn=0. The first
  round had pass=26581; the 24 more are the new assertions (119 to 143
  in `test-weights-unbounded.R`). The 15 skips are the scale-tier gate.
- Core in one session, three orders (`bash dev/release/run-sessions.sh
  wcheck-suite-frsess 2 frmtmb`; `dev/release/wcheck-suite-frsess.log`,
  "ran 3 of 3"): forward pass=26605, reverse 26770, preload 26873 (skip
  17: the two "when it is not loaded" tests), all fail=0 err=0 warn=0;
  135 LEAK lines, all tagged "[with load of ...]", none from
  `test-weights-unbounded.R`. The full session covers that file and
  every file that fits a weighted model.
- frmtmb.sample's `test-weights-unbounded.R` (the sample route reads
  the changed verdicts) against the fixed core: pass=5 fail=0 err=0
  skip=0 warn=0 (`dev/wcheck-out/fr/tf-sample-wu-lane.txt`).
- OpenBLAS 0.3.26 with LAPACK (`bash dev/ciharden-openblas.sh 0.3.26
  lapack`, LAPACK 3.12.0): `test-weights-unbounded.R` pass=143 fail=0
  err=0 skip=0 warn=0 (`dev/wcheck-out/fr/tf-wu-openblas.txt`);
  `dev/ciharden-out/` deleted afterwards.
- R CMD check --as-cran, core (`bash dev/wcheck-check.sh`, log kept as
  `dev/wcheck-out/fr/check-frmtmb.txt`, check tree deleted): tests
  "[39m] OK", `Status: 1 NOTE` (the expected V8 NOTE).

The extension suites were not rerun in this round: no extension
changed, and the first round's 193 extension files ran against the same
core API; only frmtmb.sample reads the changed verdicts, and its
weights test file was rerun above.

### Fix round 2 (the short check's S1: the Jacobian follows the link)

The first fix added the log link's Jacobian power to every link.
frmtmb adds `log |dx / d eta|` (`natural_log_jacobian()`, R/priors.R),
whose power of x depends on the link. Derived from each inverse link
x(eta):

| link | dx/deta | toward 0 | toward infinity |
|---|---|---|---|
| log | x | 1 | 1 |
| identity | 1 | 0 | 0 |
| softplus, x = log(1 + exp(eta)) | 1 - exp(-x): about x near 0, 1 at infinity | 1 | 0 |
| squareplus, x = (eta + sqrt(eta^2 + 4)) / 2 | x / sqrt(eta^2 + 4): x is about 1/abs(eta) near 0, about eta at infinity | 2 | 0 |

A new `natural_jacobian_power(link, zero)` returns these, and NA for
any other link, where `prior_entry_verdict()` now says "unknown". The
early "bounds" returns (a positive `lb`, a finite `ub`, a faster than
polynomial tail, an infinite order) do not depend on the link and stay.

Confirmed by profile on both sides of every boundary
(`Rscript dev/wcheck-fr-links.R <lib>`, log
`dev/wcheck-out/fr/links-lane.txt`; seed 20261111, n = 80; grid even
in log x, mapped to the link scale; slopes per unit of log(1/sigma) or
log shape):

| arm | check | slopes | predicted |
|---|---|---|---|
| gaussian log, normal on sigma, E 0.5 | fit -112.358, conv 0 | -0.496 -0.500 -0.500 -0.500 | -0.5 |
| gaussian log, E 1.5 | REFUSED | 0.504 0.500 0.500 0.500 | 0.5 |
| gaussian identity, gamma(2, 1) on sigma, E 0.5 | fit -112.184, conv 0 | -0.441 -0.492 -0.499 -0.500 | -0.5 |
| gaussian identity, gamma(2, 1), E 1.5 | REFUSED | 0.559 0.508 0.501 0.500 | 0.5 |
| gaussian softplus, normal, E 0.5 | fit -112.199, conv 0 | -0.467 -0.496 -0.499 -0.500 | -0.5 |
| gaussian softplus, normal, E 1.5 | REFUSED | 0.533 0.504 0.501 0.500 | 0.5 |
| gaussian squareplus, normal, E 1.5 | fit -112.683, conv 0 | -0.487 -0.500 -0.500 -0.500 | -0.5 |
| gaussian squareplus, normal, E 2.5 | REFUSED | 0.513 0.500 0.500 0.500 | 0.5 |
| Gamma log shape, E 12, student_t(5.5) | REFUSED | 0.565 0.508 0.501 0.500 | 0.5 |
| Gamma log shape, student_t(6.5) | fit -85.0966, conv 0 | -0.428 -0.492 -0.499 -0.500 | -0.5 |
| Gamma identity shape, student_t(4.5) | REFUSED | 0.559 0.508 0.501 0.500 | 0.5 |
| Gamma identity shape, student_t(5.5) | fit -86.1921, conv 0 | -0.435 -0.492 -0.499 -0.500 | -0.5 |
| Gamma softplus shape, student_t(4.5) | REFUSED | 0.559 0.508 0.501 0.500 | 0.5 |
| Gamma softplus shape, student_t(5.5) | fit -86.1921, conv 0 | -0.435 -0.492 -0.499 -0.500 | -0.5 |
| Gamma squareplus shape, student_t(4.5) | REFUSED | 0.559 0.508 0.501 0.500 | 0.5 |
| Gamma squareplus shape, student_t(5.5) | fit -86.156, conv 0 | -0.435 -0.492 -0.499 -0.500 | -0.5 |

The reviewer's cases on the fixed build
(`dev/wcheck-out/fr/short-shapelinks-lane.txt`, `short-links-lane.txt`):
Gamma E 12 student_t(5.5) on the shape fits under identity and
softplus (logLik -95.1383, betad 23.3, conv 0, as on rellib-r20) and is
refused under log; gaussian `normal(0, 1)` on sigma: squareplus fits at
E 0.5 (-112.437) and E 1.5 (-112.683, conv 0; rellib-r20 refuses it),
identity is refused at E 0.5 and 1.5 (rellib-r20 refuses both too), log
fits at 0.5 and is refused at 1.5.

Tests: a new block "the Jacobian a natural-scale prior carries follows
the link" (Gamma student_t(5.5) under identity and softplus fit at
conv 0, under log refused; gaussian squareplus E 1.5 fits at conv 0;
identity E 0.5 refused).

    rellib-r20:  RESULT test-weights-unbounded.R pass=97 fail=12 err=5 skip=0 warn=12
    fixed lane:  RESULT test-weights-unbounded.R pass=151 fail=0 err=0 skip=0 warn=0
    OpenBLAS 0.3.26 + LAPACK: pass=151 fail=0 err=0 skip=0 warn=0

(`dev/wcheck-out/fr/tf2-wu-base.txt`, `tf2-wu-lane.txt`,
`dev/wcheck-out/fr2/tf-wu-openblas.txt`.) On rellib-r20 the new block
fails as behavior: the log-link Gamma model is fitted, not refused, and
the squareplus model is refused (an error where a fit is expected).
`dev/ciharden-out/` is deleted again.

NEWS bullets 4 and 5 now give the power per link (1 for log; 0 for
identity, softplus and squareplus toward infinity; 1, 2 and 0 at 0 for
softplus, squareplus and identity) and say that any other link is not
decided. Filed in the backlog: the reviewer's pre-existing MINOR, an
invalid `class = "b", coef = "Intercept"` prior on a centered
dispersion meets the weights refusal before its own error.

Reruns: the 58-case table (`dev/wcheck-out/fr2/table.txt`) is identical
to fix round 1's; the reviewer's 24 guard-absent models
(`dev/wcheck-out/fr/common2-lane.txt`) match the base log apart from
`lib:`.

Core, one file per process, 3 at a time, `NOT_CRAN=true`,
`FRMTMB_BRMS_FIT_TESTS=true` (`dev/release/wcheck-suite-fr2.log`,
"WCHECK-SUITE-FR2 ran 317 of 317"; all 317 `lib:` lines name
`wt-wcheck-lib`): pass=26613 fail=0 err=0 skip=15 warn=0 (8 more than
fix round 1: the new block's assertions, 143 to 151). R CMD check
--as-cran, core (`bash dev/wcheck-check.sh`, log kept as
`dev/wcheck-out/fr2/check-frmtmb.txt`, check tree deleted): tests
"[28m] OK", `Status: 1 NOTE` (the expected V8 NOTE).

### Fix round 3 (short check 2's S2: the sample route reads power 1)

The MAP objective is a density in eta, so its power of x carries the
link's Jacobian power j (`natural_jacobian_power()`). A posterior's
integrability does not: with the Jacobian, the integral over eta is the
integral of p(x) L(x) dx under every link. Fix round 2 applied j to
both routes; now `prior_entry_verdict()` uses j = 1 on
`route = "sample"` for every link (0.81.0's sample conditions:
`order + 1 - E > 0` toward 0, `p - 1 - E > 0` toward infinity) and keeps
`natural_jacobian_power()` for MAP. "free" at equality on the sample
route stays: the integrand in x tends to a positive constant times 1/x.

MINOR: frechet's `nu` uses `logm1` (nu = 1 + exp(eta), dnu/deta =
nu - 1), so j = 1 as nu grows; `natural_jacobian_power()` now returns 1
there (NA toward 0, which nu never reaches). Confirmed
(`Rscript dev/wcheck-fr3-frechet.R <lib|base>`, log
`dev/wcheck-out/fr3/frechet-lane.txt`; frechet data of
`dev/reviews/wcheck-review-rates.R`, seed 20261111, n = 80, E 6, rate
1, profile along log nu 4 to 12):

| prior on nu | MAP | sample | fit | slope | predicted |
|---|---|---|---|---|---|
| student_t(5.5) | REFUSED | REFUSED | refused (rellib-r20: logLik -102.901, conv 0, a stop on an unbounded objective) | 0.503 0.501 0.500 0.500 | 0.5 |
| student_t(6.5) | pass | pass | logLik -102.907, conv 0 (as rellib-r20) | -0.491 -0.499, then numeric breakdown past log nu 8 | -0.5 |

The reviewer's route table on the fixed build
(`dev/reviews/wcheck-short2-routes.R`, log
`dev/wcheck-out/fr3/routes-lane.txt`): every sample verdict now agrees
with the posterior column. gaussian `normal(0, 1)` on sigma, E 0.5:
sample passes under log, identity, softplus and squareplus (proper);
E 1.5: sample refuses under all four (improper), squareplus included.
Gamma E 12, student_t(5.5) on the shape: sample refuses under all four
links (improper); student_t(6.5) and (7.5) pass. MAP stays as fix round
2 left it (identity E 0.5 refused, squareplus E 1.5 passed, Gamma
student_t(5.5) refused under log only). frechet's student_t(3) at E 8
is now refused on both routes, student_t(12) passes both.

Tests: the link block gains the frechet MAP verdicts (student_t(5.5)
refused, (6.5) passes); a new block "the sampler's verdict does not
depend on the link" (identity E 0.5: sample NULL, MAP refused;
squareplus E 1.5: sample improper, MAP NULL; Gamma E 12 student_t(5.5)
under identity, softplus and squareplus: sample improper).

    rellib-r20:  RESULT test-weights-unbounded.R pass=100 fail=12 err=6 skip=0 warn=12
    fixed lane:  RESULT test-weights-unbounded.R pass=160 fail=0 err=0 skip=0 warn=0
    OpenBLAS 0.3.26 + LAPACK: pass=160 fail=0 err=0 skip=0 warn=0
    frmtmb.sample test-weights-unbounded.R: pass=5 fail=0 err=0 skip=0 warn=0

(`dev/wcheck-out/fr3/tf-wu-base.txt`, `tf-wu-lane.txt`,
`tf-wu-openblas.txt`, `tf-sample-wu-lane.txt`.) On rellib-r20 the new
block errors where 0.81.0 refuses the bounded squareplus MAP objective;
the link block fails as before. `dev/ciharden-out/` is deleted.

NEWS: bullet 4 names logm1 beside log; bullet 5 says the powers are
those of `frm()`'s maximum and that `frm_sample()`'s check reads 1
under every link, as before.

Reruns: the 58-case table (`dev/wcheck-out/fr3/table.txt`) is identical
to fix round 2's; the reviewer's 24 guard-absent models
(`dev/wcheck-out/fr3/common-lane.txt`) match the base log apart from
`lib:`.

R CMD check --as-cran, core (`bash dev/wcheck-check.sh`, log kept as
`dev/wcheck-out/fr3/check-frmtmb.txt`, check tree deleted): tests
"[30m] OK", `Status: 1 NOTE` (the expected V8 NOTE). The R CMD check's
test step runs the whole core suite in one session. A separate
one-file-per-process core suite was not rerun in this round; the change
since fix round 2's 317-of-317 run is the sample route's power and one
link entry, inside `prior_entry_verdict()` and
`natural_jacobian_power()`, and every caller of those runs in
`test-weights-unbounded.R` and frmtmb.sample's weights test above.

### For the merge

- Lane priortab inserts `drop_unbounding_slots()` just before the
  check's call in `R/fit.R`; rerun `test-weights-unbounded.R` after the
  merge, since the check then sees priortab's filtered prior list.
- Lane quadseam moves its runaway test to an `mi()` case; the check
  skips `mi()` by design and keeps doing so.
