## Filed at the 0.70.0 release (2026-10-13)

What the three lanes of the round opened 2026-10-07 (`moparity`,
`optrap`, `seflags`) found and did not fix, and what their reviews and
the consolidation listed. Each item names its lane and the record with
the measurement: `dev/<lane>-findings.md` and
`dev/reviews/2026-10-08-<lane>.md`. `dev/round-20261013.md` is the
round's record. Lanes moparity and seflags added their own sections;
they follow this one as the lanes wrote them, one heading level down.
Lane optrap added none; its filings are below. Each item this round
closed is marked in place above with "CLOSED at 0.70.0".

### Closed by this round (were open at 0.69.1)

- The gaussian random-intercept fit short of the maximum at units
  1e3; the restart that raised and discarded a finite first run (fuzz
  seed 20379118 under OpenBLAS LAPACK); the objective's state after a
  discarded run (`obj_settle()` on every path); `mo()` seed 194; the
  test of `test-id-kron.R`, which no longer depends on which of its two
  points a core type reaches (optrap).
- The estimable predictions flagged on fixture 1 (`data-helpers:7`)
  under OpenBLAS LAPACK; the boundary message on an sd the likelihood
  does not read; curved ridges above 1e-9; the degree-7 raw polynomial
  with `(1 | g)` (seflags). The flat simplex direction off the
  coordinates is much reduced but not closed (Open - low, below).
- `mo(x) * f` and `(mo(x) | g)`; `simo` and `sbhaz` Dirichlet priors;
  `update(prior =)` merging as brms; the `cox()` baseline's plateau;
  the dense blocks' `+Inf` and `NaN` log density at an underflowed log
  sd; `parnames()` deciding from the binding; the exact `gp()` field
  sampled centered; `as_tmbstan()`'s message on `mo()` fits
  (moparity).

### Open - high

- **The wide high-R^2 design, seed 1, stops 9107.66 below lmer at
  code 0 with no warning**, on 0.69.0 and on 0.70.0 alike
  (`dev/optrap-p2-wide1.R`, `p2/wide1-r7.txt`; n 101000, 50 columns).
  3 iterations, 9 evaluations; the group sd ends at 6.5e-16 against
  lmer's 0.486 and sigma at 1.11 against 0.999. The start is the
  cause: sigma at `log(sd(y))` = 9.21 and every slope at 0, so the
  first steps put the group sd at the boundary, where the log-scale
  gradient vanishes. The response units do not catch it (the residual
  spread is 1, no unit), and `newton_polish()` does not run (the
  design is not near collinear, so its random-effect steps do not
  apply). Measured candidate (`dev/optrap-p2-wide1-starts.R`): sigma
  started at the residual MAD and the mu coefficients at least squares
  reach lmer at code 0 (-1.3e-8, 9 iterations). That changes the
  starting values of every location-scale fit whose response is far
  from its residual spread, across ten families' sigma conventions, so
  it needs its own fit diff over the suites: its own round. Lane
  optrap, punch round 2; its review's consolidation item 2.
- **Accurate SEs on ill-conditioned designs** is its own lane: see
  lane seflags' section below ("Open - high priority"). 0.70.0 keeps
  0.69.0's NaN for the coefficients an ill-conditioned direction
  carries, with a message that names the design, after the punch
  round 2b attempt gave wrong finite SEs through emmeans,
  marginaleffects, `vcov()` read whole and the ordinal category route
  (review B6 and B7).

### Found at the consolidation (the merged build), fixed before the release

Two lanes can each be right and their union wrong. The reruns of the
three reviews on the first merged build found these three
(`dev/round-20261013.md`, "The reruns the reviewers asked for" and
"Found at the consolidation, fixed"); lane optrap fixed them in its
consolidation fixes (`dev/optrap-findings.md`, "Consolidation fixes",
install 32), and the release tree took the fixes by a rerun of
`dev/rel070-merge.sh`.

- **`newton_polish()` accepted a step whose actual gain was a tiny
  fraction of the predicted one.** CLOSED at 0.70.0 (lane optrap,
  consolidation fixes). On the false-alarm battery's
  `negbin_large_shape` (negative binomial data near the Poisson limit,
  where the shape's MLE is infinite), seed 2 under the reference BLAS:
  the Hessian's eigenvalues were 700.7, 310.6 and 9.9e-7, the actual
  gain 2.965e-8 of a prediction above 1e-6, and the step took the log
  shape from 19.98 to 66.76, where the gradient is 565; every SE was
  NaN, `b_Intercept` and `b_x` included, behind a large-gradient
  warning (`dev/rel070-seflags-attrib.R`). Over 200 seeds
  (`dev/rel070-polishgate.R negbin`): 0.69.1 warned on 66 and lost no
  coefficient SE; the first merged build polished 11, warned on 74 and
  lost both coefficient SEs on 3. The fix: a step is kept only if it
  gains at least 0.1 of its prediction (the consolidation's candidate
  0.25 refused a needed first step of a degree-7 raw polynomial, seed
  17, whose actual gain is 0.203 of the prediction), and steps that
  end where the Hessian has no Cholesky factor are undone.
- **Raw polynomials with `(1 | g)` stopped short of the maximum in
  silence.** CLOSED at 0.70.0 (lane optrap, consolidation fixes). On
  the first merged build, against an exact profile maximum
  (`dev/rel070-optrap-diag-repoly.R`, 20 seeds per degree): degree 6,
  13 of 20 short (to 1.1), all silent (0.69.1: 5, to 0.0229, 4
  silent); degree 7, 20 of 20 short (to 3.1), all silent, where 0.69.1
  warned through the SE check's false "flat" verdict, which lane
  seflags removed. The cause is lane optrap's coefficient units (the
  lane's ablations: without them these fits are 0.69.1's). The fix:
  Newton steps on fits with random effects too, on a Hessian by
  central differences of the exact gradient (`stats::optimHess()` in
  natural units), when the fixed design is near collinear
  (`predictor_design_kappa(frame) > 1e5`: 3e5 at degree 5, 4e6 at
  degree 6, 3.4 for sleepstudy), so other fits pay one QR of their
  design and are otherwise unchanged.
- **The polish could end where the convergence check rejects the
  gradient, so the SE check did not run** (medium). CLOSED at 0.70.0
  (lane optrap, consolidation fixes): on seflags' B5 design without
  random effects (degree 8, seeds 4 and 5) every SE was NaN, z, wb and
  sigma included; the steps ended where the Hessian does not factor
  and are now undone, so the fit is lane seflags' again.
- Not a defect, recorded: ported fixture 1 under OpenBLAS LAPACK stops
  at code 1 (lane optrap's unit 4 on the mu of `count`), so seflags'
  B1 premise there no longer arises; with `predictor_units()` off it is
  code 0 with seflags' lost set (`dev/rel070-seflags-fx1.R`). The B1
  constructions that hold under both BLAS (`b1ref`) keep their test.

### Open - medium

- **A group sd 1000 times the residual sd gives "false convergence"
  at lmer's maximum** on 2 of 20 fits (seeds 6 and 14, unit 1020;
  `re1000-lane.txt` of the review), where 0.69.0 gives 0. The MAD of
  the residuals holds the group variance (unit 1024), so the group sd
  starts right, but sigma starts at the family's own `sd(y)`, 1000
  times its value. The warning is right about the optimizer and wrong
  about the fit. The sigma start of the item above would remove it.
  Lane optrap, review r2.
- **The convergence-warning policy**: the SE check does not run on a
  fit with a convergence warning, so on a ridge the reported
  uncertainty flips with rounding (lane ci691; lane seflags' section
  below, "Pre-existing"). It covered seflags' B5 seeds 4 and 5
  without random effects on the review's merge build (fixed above, by
  undoing the polish). A policy decision for the user before a lane.
- **The ported row `brmsfit-methods:326` fails again**: fixture 1
  (student, `mo()`, `s()`, `arma`, `sigma ~ Trt`) stops at "false
  convergence (8)" with an sd at -723 to -1253 on the log scale, on an
  objective that is unbounded within reach (from one of 41 jittered
  starts lane optrap ran it to -1308.3, `dev/optrap-rev-fix1.R`). At
  the stop the sigma coefficients' `sdreport()` variances are
  negative, so `fitted(fit1, dpar = "sigma")` has NaN errors. 0.69.0's
  pass was luck of the stop point: with 40 jittered starts (sd 0.3) the
  row holds from 34 of 41 on rellib-r7 and fails from the other 7, and
  lane optrap's unit (8 on this fixture) changes which stop the default
  start reaches (`dev/optrap-fix1.R`; optrap review, item 6 and
  consolidation item 1). The ledger records it as cannot transfer
  again (fixture). The remedy is the fixture's degeneracy, not the
  test.

### Open - low

- **A residue of the exactly flat simplex direction off the
  coordinates** (the consolidation review, m1; reopens the 0.69.0 item
  above). Over 40 seeds of four designs (`dev/relrev070-segap.R`, the
  review's `segap/` logs), the fits that are silent with every zeta SE
  finite: 64 of 160 on 0.69.1, 7 with lane seflags alone, 71 with lane
  optrap alone, 8 on the release under either BLAS (`ls ~ mo(inc)` 2,
  `ls ~ mo(inc) + age` 2, `ls ~ mo(inc) * age` 0, `y ~ mo(x2)` with n
  300 4). The release names the direction on 148 (reference) and 149
  (OpenBLAS) of 160; the others warn otherwise. A silent fit reports,
  for example, zeta SEs 547.9 and 949.0 (seed 10 of `y ~ mo(x2)`) on a
  direction the likelihood does not read. Pre-existing and much
  reduced; every fit is at its maximum. Lane seflags' ground.
- **The least-squares subsample of the unit reads evenly spaced
  rows**, so on rows sorted by a covariate with strong
  heteroscedasticity (sd `exp(1.5 x)`) and a subsample of 100 rows the
  unit is more than one power of two off on 8 of 20 seeds (2 of 20
  unsorted; `dev/optrap-rev3-sub2.R`). At the cap's real sizes (5e6 /
  p rows, 5e4 or more for up to 100 columns) neither arm moves, and
  the worst seed's fit with the subsample's unit still ends within
  4.4e-10 of `lm()`. Harmless; a seeded random subsample would remove
  the dependence on row order. Lane optrap, review c3.
- **`y + 1e9`**: a response with a huge offset gets no unit from its
  residual spread and fits as before (lane optrap, m4).
- **Mixture components get no unit** (lane optrap, m3, measured).
- **huber()'s Laplace objective is discontinuous** and stops on either
  side of a jump, on every build (lane optrap, m6).
- **Raw polynomials of degree 8 and up** are singular in double
  precision; Newton steps cannot help, and every such fit is warned
  (lane optrap; see seflags' ill-conditioned item). At degree 7
  without a group term 1 of 20 stays 4.0e-4 short in silence (lane
  optrap's consolidation fixes, seed 4).
- **The one-level group with a bound** (`y ~ x + (1 | g1)`, g1 one
  level, `ub = 0.1` on x): nlminb zigzags to its evaluation limit on 5
  of 40 seeds on optrap's round-0 build, 0 of 40 on rellib-r7; with the
  residual unit the test's own data have unit 1 (lane optrap, round 0
  and punch round 1). The model is degenerate and frmtmb warns so.
- **fit_loss2 and F11 seed 8** stay short (0.032, 0.1418), as lane
  optima filed; their causes are not the response's units or the
  restart (lane optrap).
- **The objective keeps the constant `n log u`**, so the invariance of
  a fit to the response's units is 1e-7 in the log-likelihood rather
  than exact; removing it made relative convergence stricter and broke
  converged `mo()` fits (lane optrap).
- **`frmtmb.eam` `test-defects.R:92`'s mixing weight at 6e-11** is
  called "flat", where "a mixing parameter at its limit", in the
  boundary message that exists, explains the fit better (lane
  seflags' section below; review consolidation list).
- **Read every `||` term as brms does** (moparity's r2, in its section
  below). A core parsing change; its own lane.
