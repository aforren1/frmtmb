# The classic TMB examples: what frmtmb replicates today

Audit date: 2026-09-03. Branch `wt-tmbex`, at v0.41.0 (2619816).

Revised 2026-09-03 on branch `wt-sparsear1`: roadmap item 1 landed, so
`sdv_multi` and `sdv_multi_compact` move from AWKWARD to REPLICATES at
their published size and the `transform` rows are no longer capped. The
rows and counts below are the revised ones.

## What was audited

The universe is the `tmb_examples` directory of two upstream repositories,
enumerated from the repositories themselves rather than from any list:

* `kaskr/RTMB/tmb_examples` holds 19 example models as plain R. They run
  in the installed RTMB with no compilation, so they give exact reference
  fits. A vendored copy is in `dev/tmb-examples/`.
* `kaskr/adcomp/tmb_examples` holds 31 C++ examples. 15 of them are the
  originals of the RTMB ports. The other 16 have no port, so they are the
  fallback set. A vendored copy is in `dev/tmb-examples-adcomp/`.

Total: 35 distinct examples. Both vendored copies are gitignored.

Every verdict is evidenced. For a REPLICATES row the reference model is
fitted through its own `MakeADFun`/`nlminb` and `frm()` is fitted in the
same R process, and the two maximized log-likelihoods are compared. Both
sides tape the same objective through RTMB, so agreement is expected at
1e-6 or better. Rows that agree less closely, and rows where the
reference's own optimizer stops short of its optimum, are explained
below; none of them is a model difference.

The harness is `dev/tmb-examples-check.R` (one R subprocess per case,
because the reference scripts `source()` their data files into the global
environment). The Kalman prototype is `dev/tmbex-kalman-prototype.R`.
The REPLICATES rows are also regression tests in
`tests/testthat/test-tmb-examples.R`, which downloads nothing: the data
generation and the reference negative log-likelihood are both inline.

## Verdict counts

| verdict | count |
| --- | --- |
| (a) REPLICATES | 15 |
| (b) AWKWARD | 1 |
| (c) OUT OF SCOPE | 19 |

Of the 19 out-of-scope rows, 3 are TMB API demonstrations whose
underlying model is another row that replicates, 3 are OpenMP parallel
builds of models that already replicate, and 13 need a capability frmtmb
does not have. (A fourth API demonstration, `dataeval`, is counted under
REPLICATES, because its model fits.)

## (a) REPLICATES

`|diff|` is the absolute difference between the reference maximized
log-likelihood and `logLik(frm(...))`.

Three rows carry a `+`. On those the reference's own cold-start `nlminb`
stops short of its optimum, so the number quoted is the reference
objective CROSS-EVALUATED at frmtmb's estimates against
`logLik(frm(...))`. That is a statement about the two likelihoods, which
is what the audit is about; agreement of the two optima would only be a
statement about two optimizers. Each `+` row says below by how much the
reference's stopping point falls short, and in every case `frm()`
reaches the better objective.

| example | source | `frm()` spelling | `|diff|` |
| --- | --- | --- | --- |
| linreg | RTMB | `bf(Y ~ x)`, `family = gaussian()` | 4.19e-13 |
| dataeval | RTMB | `bf(y ~ x)`, `family = gaussian()` | 1.15e-11 |
| tweedie | RTMB | `bf(y ~ 1)`, `family = tweedie()` | 5.46e-12 |
| compois | RTMB | `bf(x ~ 1)`, `family = compois()` | 1.31e-10 |
| spatial | RTMB | `bf(y ~ 1 + x2 + exp(pos + 0 \| grp))`, `family = poisson()` | 1.01e-08 |
| spde | RTMB | `bf(time \| cens(cens) ~ sex + age + wbc + tpi + spde(fem, gr = node))`, `family = weibull()` | 9.35e-09 |
| adaptive_integration | RTMB | `bf(x \| trials(n) ~ 0 + c1 + c2 + c3 + (1 \| obs))`, `family = binomial()`, `quadrature = TRUE` | 1.36e-03 |
| transform | RTMB | `bf(y ~ qgamma(pnorm(z), shape, scale), z ~ 0 + ar1(tim + 0 \| g), shape ~ 1, scale ~ 1, nl = TRUE)` | 3.93e-09 |
| transform2 | RTMB | `bf(y ~ qbeta(pnorm(z), shape1, shape2), z ~ 0 + ar1(tim + 0 \| g), shape1 ~ 1, shape2 ~ 1, nl = TRUE)` | 1.82e-12 + |
| sdv_multi | RTMB | `mvbf(bf(x1 ~ 0, sigma ~ 1 + ar1(tim + 0 \| g)) + gaussian(), ... , rescor = TRUE)` | 8.50e-09 + |
| sdv_multi_compact | RTMB | the same three-response `mvbf()` spelling | 8.50e-09 + |
| orange_big | adcomp | `bf(y ~ a0 / (1 + exp(-(t - a1) / a2)), a0 ~ 1 + (1 \| tree), a1 ~ 1, a2 ~ 1, nl = TRUE)` | 2.33e-08 |
| socatt | adcomp | `bf(y ~ x1 + x2 + x3 + (1 \| g))`, `family = cumulative()` | 1.35e-08 |
| lr_test | adcomp | `bf(obs ~ 0 + g, sigma ~ 0 + g)` and its two restrictions | 1.4e-09 |
| longlinreg | adcomp | `bf(Y ~ x)`, `family = gaussian()` | 1.86e-09 |

### Rows that need a word

**dataeval.** The example demonstrates `DataEval()` and `Tape$atomic()`:
one tape is built for a chunk of data and reused for ten chunks, with the
chunk index smuggled in as a parameter. The fitted model is one pooled
linear regression across all 100 rows, and that is what `frm()` fits. The
API being demonstrated has no `frm()` surface and does not need one, but
it does reappear below: `Tape$atomic()` is the reason the Kalman node
needs no hand-written adjoint.

**adaptive_integration, at 1.36e-03.** This is the one loose row, and the
looseness is a quadrature rule, not a model difference. The reference
marginalizes each scalar random effect with `stats::integrate()`, which
is adaptive QUADPACK over the whole real line. `frm(quadrature = TRUE)`
uses TMB's experimental `integrate=` with the spec
`list(dim = 1, adaptive = FALSE, method = "marginal_gk")`, a fixed
Gauss-Kronrod rule. Over 1000 observations the accumulated rule
difference is 1.36e-03 on an objective of -1719.39, a relative difference
of 7.9e-07. Both routes are exact marginalizations of the same model to
within their own rule error, so this row is a REPLICATES with an honest
tolerance rather than an AWKWARD. Turning `adaptive` on in that spec is
the way to close it, and it is a one-line change in `R/fit.R` rather than
a new capability. It is listed in the roadmap below.

**transform and transform2, and the pin they need.** The reference draws
the latent field with `dautoreg(u, phi = phi)`, whose scale is fixed at 1,
because the field is pushed through `pnorm()` and must be marginally
standard normal for `qgamma`/`qbeta` to be the intended quantile
transform. frmtmb's `ar1()` block carries a free marginal standard
deviation, so the plain spelling fits a strict superset with one extra
parameter and lands above the pinned fit: 0.066 log-likelihood units at
n = 1000, 0.099 at n = 200. Pinning that parameter with `lower = c(theta_1 = 0)` and
`upper = c(theta_1 = 0)` is the `map=` equivalent and recovers the
reference exactly.

Both rows now run at the example's own n = 1000. Through v0.42 they were
capped at n = 200, because `ar1()` was a dense covariance structure whose
cost is cubic in the series length; the O(d) density of v0.43 removes the
cap. Nothing about the pin story moved with it. At n = 1000:

| row | reference | `frm()` pinned | pinned rho | reference phi |
| --- | --- | --- | --- | --- |
| transform | -2829.4068542850 | -2829.4068542811 | 0.542386 | 0.542383 |
| transform2 | 685.3439060286 | 745.1140235944 | 0.554231 | 0.516259 |

`transform2` is one of the `+` rows of the REPLICATES table. Its reference script
says in its own comments that the Laplace approximation is delicate for
that model (the noise sd is 0.005), and the reference's cold-start
`nlminb` stops 59.77 log-likelihood units below `frm()`. The likelihood
itself agrees: the reference objective evaluated at `frm()`'s pinned
estimates is 745.1140235944, matching `logLik(frm(...))` to 1.82e-12. For
`transform` the same cross-evaluation agrees to 7.22e-09.

The n = 200 runs the audit published still reproduce, which is the check
that the O(d) density changed no number: `transform` -544.1792742616
against a reference -544.1792742715, and `transform2` 161.3844854981
against 161.3844854988, with the pinned rho matching the reference phi to
five decimals in both (0.512781 against 0.512783, and 0.533422 against
0.533424). `dev/tmb-examples-check.R` keeps both sizes as the cases
`transform`/`transform2` and `transform_n200`/`transform2_n200`.

One naming note the harness has to undo. `RTMB::qgamma()` takes a RATE in
its third positional argument, so the coefficient the spelling above
calls `scale` is fitted as 1 / scale. It is the same gamma either way,
which is why the row's log-likelihood matched all along; only the
cross-evaluation into the reference's parameterization has to invert it.

**sdv_multi and sdv_multi_compact, at their published size.** The
multivariate stochastic volatility model of Skaug and Yu (2014). Three
series, 945 time points, a latent AR(1) log-volatility per series, and an
observation density that is multivariate normal with mean zero, an
unstructured correlation matrix, and per-series scales
exp(0.5 (mu_x + h)). Both scripts fit the same model; the `_compact` one
is the same likelihood written with `dautoreg()` and a scaled `dmvnorm()`.

```r
frm(mvbf(bf(x1 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
         bf(x2 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
         bf(x3 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
         rescor = TRUE), data = dd)
```

The mapping is exact rather than approximate. The residual correlation of
a multivariate gaussian fit is the example's `R`, and both sides
parameterize it by the same row-normalized lower-triangular factor, so
`thetar` IS `off_diag_x`. `sigma`'s log link makes its intercept mu_x / 2
and its latent block h / 2, so `ar1()`'s marginal sd absorbs the factor
of 2 exactly. Nothing is approximated.

This row was AWKWARD through v0.42 for cost alone: `ar1()` built a dense
945 x 945 covariance and Cholesky-factorized it on the tape, and the
audit could only reach series lengths of 40 and 80 on simulated data of
the same shape. Since v0.43 the block density is O(d) and the published
size fits in 2.4 s, against a reference that takes 17.4 s.

The two optimizers do not land together, and that is the reference's
side. The example's own single cold start returns `convergence = 1` at
nll 1773.08785102; `frm()` reaches 1771.12753479, which is strictly
better. So the row is established as a `+` row:

* the reference objective evaluated at `frm()`'s estimates is
  1771.12753478 against `frm()`'s own -1771.12753479, agreeing to
  8.33e-09, and
* the reference's `nlminb`, restarted from `frm()`'s point, converges
  there (`convergence = 0`) at 1771.12753478, giving the 8.50e-09 in the
  table.

The 12-parameter box-constrained problem is harder to converge from a
cold start than the likelihood identity is to establish, which is the
same thing the audit saw at series length 80.

**orange_big.** The scaled-up Orange Tree model: 35000 observations,
5000 latent random effects, a three-parameter logistic growth curve with
a random asymptote. The C++ template writes the curve parameters as
`192 + beta[0] + u`, `726 + beta[1]`, `356 + beta[2]`; those constants are
starting-value offsets, so the fitted curve is the same and `frm()`
reaches it from `start = list(beta = c(192, 726, 356))`. Reference 1.9 s,
`frm()` 3.3 s, both including taping.

**socatt.** Cumulative logit with 7 categories, 3 fixed effects and a
random intercept over 264 groups. The ADMB original writes the linear
predictor as `sigma * u` with `u ~ N(0, 1)`; frmtmb's `(1 | g)` is the
same block in its natural parameterization, so the log-likelihoods match
without any pin. The template adds `1e-20` inside the log for safety,
which is below the agreement tolerance.

**lr_test.** The example uses TMB's `map=` to fit three nested models on
a ragged array and run two likelihood-ratio tests. In frmtmb the three
models are three formulas, and the restriction is a formula restriction
rather than a map. Reference and `frm()` agree on all three: 56.9304092428
against 56.9304092414, 59.2697078909 against 59.2697078911, and
64.2211680757 against 64.2211680757. Parameter counts agree too (10, 6, 2),
so `anova()` reproduces the example's chi-square table.

**longlinreg.** One million observations, three parameters, agreeing to
1.86e-09. The point of the example is scale, and the formula front end
costs nothing measurable at this size: reference 22.1 s against `frm()`
14.9 s in the harness run, where the reference pays for the AD Hessian
the example hands to `nlminb` and `frm()` does not. On a quieter machine
the two were 7.9 s and 8.4 s. Read the timings in this document as
indicative rather than as a benchmark; `dev/benchmarks.md` is the place
for measured performance.

## (b) AWKWARD

### matern

`bf(z ~ 0 + mat(pos + 0 | grp))`, `family = gaussian()`.

The `mat()` covariance structure carries the same Matern correlation as
the example's `geoR`-style kernel, but the reference field has unit
marginal variance and no nugget, and the frmtmb spelling has a switch for
neither. The plain formula therefore fits a four-parameter superset:
marginal sd, range, shape, and a residual sigma. Pinning the two extra
parameters through `lower`/`upper` on `theta_1` and `sigma_(Intercept)`
recovers the reference's phi and kappa to four decimals: reference
phi = 1.176572, kappa = 2.944163 against pinned phi = 1.176535,
kappa = 2.944253. Cross-evaluating the reference objective at frmtmb's
pinned estimates gives -100.22169399 against -100.22169404 at its own,
so the kernel itself agrees to 5e-08 and only the optimizer's stopping
point differs. The reported `|diff|` of 7.48e-04 is that stopping-point
difference on a nearly flat ridge, amplified by pinning a residual sigma
at exp(-20) rather than at zero.

Worth new sugar: yes, and it is small. A `nugget = FALSE` argument on the
spatial structures, or a documented `unit = TRUE` that fixes the marginal
variance, would turn this row and both `transform` rows into plain
REPLICATES with no pinning ritual. Three of the 35 examples needed the
same pin, which is the argument for spelling it once.

## (c) OUT OF SCOPE

### Needs a capability frmtmb does not have

| example | source | what blocks it | atomic-function class |
| --- | --- | --- | --- |
| mvrw | RTMB | no latent random-walk term in the grammar | Kalman-filter node |
| mvrw_sparse | adcomp | same model as mvrw | Kalman-filter node |
| sde_linear | adcomp | latent states on a fine time mesh, most unobserved | Kalman-filter node |
| thetalog | adcomp | nonlinear latent state process | none (see below) |
| rickervalidation | RTMB | nonlinear latent state process | none (see below) |
| sam | RTMB | bespoke fisheries state-space model | genuinely never |
| hmm | RTMB | continuous SDE discretized to a CTMC, needs `expm` of a 100 x 100 generator | expm |
| ar1_4D | RTMB | separable AR(1) over four index dimensions | none, needs a Kronecker combinator |
| spde_aniso | adcomp | anisotropic `Q_spde` with an `H` matrix | none, needs an spde option |
| spde_aniso_speedup | adcomp | same | none, needs an spde option |
| fft | adcomp | circulant covariance evaluated in the Fourier domain | FFT node |
| nmix | adcomp | discrete latent abundance summed over its support | genuinely never |
| spa_gauss | RTMB | saddlepoint marginalization (`Tape$laplace(SPA = TRUE)`) | genuinely never |

**thetalog and rickervalidation are not Kalman candidates.** Both have a
Gaussian latent process whose mean is a nonlinear function of the previous
state, so a linear-Gaussian Kalman filter is not an exact marginalization
of either. The Laplace approximation the examples already use is the
better answer. What frmtmb lacks for these is the same thing it lacks for
mvrw: a way to spell a latent state process at all. A Kalman node would
not help them; a state-space term would, with Laplace behind it.

**hmm is the expm case.** The example discretizes a scalar SDE onto a
101-point grid by finite volume, exponentiates the resulting generator
once per fit, and runs a forward filter over the transition matrix.
frmtmb has `hmm()`, but it fits a discrete-state chain whose transition
matrix is parameterized directly by logits, not one obtained as
`expm(A dt)` from a generator built out of model parameters. The missing
piece is a matrix exponential on the tape. That is the same node the
pharmacometrics tier in `dev/feature-gaps.md` wants for linear
compartment models, which is what makes it worth more than one example.

**nmix and spa_gauss are the genuinely-never rows.** `nmix` sums over a
discrete latent abundance N at each site; that is a summation the
formula compiler has no representation for and no plan to grow one.
`spa_gauss` marginalizes by saddlepoint rather than Laplace, which is a
property of the tape's inner problem rather than of the model, and
frmtmb's objective builder does not expose it.

### TMB API demonstrations, not models

These four demonstrate a TMB or RTMB interface feature. The model
underneath each is a row that already replicates, so there is nothing
separate to add.

| example | demonstrates | model underneath |
| --- | --- | --- |
| checkConsistency | `TMB::checkConsistency()` over `OBS()`-marked responses | transform, transform2, sdv_multi, tweedie, compois |
| laplace | a hand-written Laplace approximation using `autodiff::` | spatial (REPLICATES, 1.01e-08) |
| register_atomic | `REGISTER_ATOMIC` to shrink a tape | adaptive_integration (REPLICATES, 1.36e-03) |
| dataeval | `DataEval()` and `Tape$atomic()` for tape reuse | a pooled linear regression (REPLICATES, 1.15e-11) |

`checkConsistency` is the one with a partial frmtmb surface: `frm()` marks
univariate non-matrix responses with `OBS()`, so `TMB::checkConsistency()`
runs on a `frm()` fit's object. It is not wrapped and not tested here.

### OpenMP parallel builds

`linreg_parallel`, `transform_parallel`, `register_atomic_parallel` and
the parallel half of `longlinreg` build the same likelihood with
`parallel_accumulator`. RTMB has no OpenMP parallel accumulation, which
`SPEC.md` records as the one real gap against classic TMB and as not on
the critical path. The serial models replicate; the parallelism does not
and will not. Genuinely never, unless a compiled fast path is added
later, which the fit object's backend-agnostic design leaves open.

## The Kalman prototype

Scoping question: is a Kalman-filter node worth building as an
`RTMB::ADjoint` atomic with a hand-written adjoint?

`mvrw` is the right example to answer it on. It is a three-dimensional
local-level model: a random walk with correlated increments, observed
with independent Gaussian noise, and a flat improper prior on the first
state. It is linear and Gaussian, so the Laplace approximation is exact.
That gives a correctness check nothing else offers: the same model
marginalized two ways, by Laplace over the 300 latent states and by a
Kalman filter in closed form, must agree at every parameter value.

The flat prior on the first state is what makes the two routes agree
constant and all. Integrating `N(y1; u1, R)` over `u1` gives exactly 1, so
the filter starts at `a = y1`, `P = R` and the first observation
contributes nothing. That is exact diffuse initialization, and it is what
the reference's missing initial-state density already encodes.

Three routes, on the example's own simulated data (3 states, 100 time
steps), measured in one R process. Node counts are from RTMB's own tape
printer.

| route | tape nodes | tape time | fit time | logLik |
| --- | --- | --- | --- | --- |
| Laplace, as `mvrw.R` writes it (300 random effects) | 6533 | 0.03 s | 0.05 s | -555.41705224 |
| Kalman filter, plain RTMB operations | 18276 | 0.03 s | 0.02 s | -555.41705224 |
| Kalman filter, `MakeTape(...)$atomic()` | 9 | 0.02 s | 0.02 s | -555.41705224 |

The identity holds:

* `|nll_laplace - nll_kalman|` at four parameter vectors, including the
  starting values, two arbitrary interior points and the optimum:
  2.27e-13, 0.00e+00, 5.68e-13, 1.14e-13.
* `|logLik_laplace - logLik_kalman|` at each route's own optimum:
  1.14e-13.
* Largest difference between the two fitted parameter vectors: 1.37e-12.
* Largest difference between the two `sdreport()` standard errors:
  9.03e-13.

### What the numbers say

**A hand-written adjoint is not needed.** `MakeTape(...)$atomic()`, which
is the same RTMB feature `dataeval.R` demonstrates, collapses the filter
from 18276 outer-tape nodes to 9, with no derivation at all: the adjoint
is the inner tape's own reverse sweep. That is the whole benefit
`REGISTER_ATOMIC` gives in classic TMB. Deriving the Kalman filter's
adjoint by hand is real matrix calculus, mostly in the covariance
recursion `P <- P - K P`, and it would have to be redone for every change
to the parameterization. The measurement says it would buy nothing here,
so the `ADjoint` half of the question is answered without doing it: the
cost is not justified while the automatic route reaches the same tape
size.

**The marginalization is not the bottleneck either.** At the example's own
size the Laplace route fits in 0.05 s with a smaller tape than the plain
filter. Both tapes grow linearly in the series length, and the Laplace
tape grows more slowly: 65 nodes per time step against the filter's 184.
The same holds for the fit, on the same simulated model at four lengths:

| time steps | Laplace fit | Kalman fit | `|nll difference|` |
| --- | --- | --- | --- |
| 100 | 0.20 s | 0.09 s | 2.27e-13 |
| 500 | 0.33 s | 0.36 s | 2.00e-11 |
| 2000 | 1.44 s | 1.77 s | 3.49e-10 |
| 5000 | 4.16 s | 4.47 s | 2.26e-09 |

Both routes are LINEAR in the series length, because the Laplace
Hessian here is block tridiagonal and TMB's sparse factorization
exploits that. Which of the two is the faster constant is
host-dependent: the sweep above had Laplace narrowly ahead past 500
steps, while the pre-merge review's rerun on the same machine had the
filter about twice as fast throughout. Same order either way, and a
constant factor in either direction does not change the conclusion:
marginalizing a linear-Gaussian state space is already a solved
problem in this backend.

The growing `|nll difference|` is the inner Newton solve's tolerance, not
a disagreement: it stays around 1e-13 relative to an objective that grows
with the series length.

**What actually blocks mvrw in frmtmb is the formula grammar.** There is
no way to write a latent random walk. `ar1(rho -> 1)` is not it: an
AR(1) block is stationary with a finite marginal variance, a random walk
is not, and pushing rho toward 1 changes the parameterization of the
process rather than approaching the right one. The prototype confirms
that once the model can be written down, the existing Laplace machinery
marginalizes it exactly at competitive cost. So the roadmap item is a
state-space term in the grammar, and the Kalman node is at most an
optimization behind it, not the thing that unlocks the example.

## The ranked roadmap

Ranked by examples unlocked per unit of work, with the measurements
above as the evidence. The first three items are not atomic functions,
which is the audit's main conclusion: the atomic-function classes in the
original scoping are mostly not what stands between frmtmb and these
examples.

**1. A linear-cost AR(1) covariance structure. DONE in v0.43.**
Unlocked `sdv_multi` and `sdv_multi_compact` at their published size, and
took `transform` and `transform2` from 200 usable time points to their
own 1000. Four examples, and every user with a long repeated-measures
series.

`ar1()` and `hetar1()` no longer build a dense `d x d` covariance. Both
evaluate the AR(1) density in its innovation form,
`z' C^-1 z = z_1^2 + sum (z_i - rho z_(i-1))^2 / (1 - rho^2)` with
`log|C| = (d - 1) log(1 - rho^2)`, which is O(d) per level and vectorized
over levels, with the marginal standard deviations divided out first and
their Jacobian added back. The derivation, including why `RTMB::dautoreg`
assumes unit MARGINAL rather than unit innovation variance, is the
comment block above `ar1_lpdf()` in `R/covstruct.R`.

Three routes were derived. `RTMB::dautoreg` computes the same recursion
but loops over time points in R and takes one vector, so it has to be
called once per level; a `dgmrf` route with the tridiagonal precision
matches the speed but assembles a parameter-dependent sparse matrix and
pays a sparse Cholesky for a log-determinant that is available in closed
form, and `hetar1`'s per-time scaling would have to be folded into that
precision. The closed form is both the cleanest and the fastest. The
other two are kept as gates rather than as code: `tests/testthat/
test-sparsear1.R` holds the density against a `dautoreg` cross-check, and
`dev/sparsear1-bench.R` reports the timings below.

Measured by `Rscript dev/sparsear1-bench.R 5` (R 4.6.1, RTMB 1.9), as
tape construction plus one `fn()` and one `gr()` call on the block
density alone, best of five batches:

| block dimension d | `ar1` dense | `ar1` O(d) | speedup | `hetar1` dense | `hetar1` O(d) | speedup |
| --- | --- | --- | --- | --- | --- | --- |
| 50 | 5.9 ms | 1.5 ms | 3.9x | 6.7 ms | 1.5 ms | 4.3x |
| 200 | 0.035 s | 2.6 ms | 13x | 0.030 s | 2.2 ms | 13x |
| 800 | 0.587 s | 2.0 ms | 294x | 0.767 s | 2.0 ms | 385x |
| 2000 | 22.698 s | 2.1 ms | 11020x | 35.096 s | 2.8 ms | 12561x |

The O(d) column is flat because at these sizes it is tape construction,
not arithmetic. The dense column is the cubic factorization. The dense
values are reconstructed inline by the script, so the comparison needs
nothing but the installed package.

The repeated-measures shape the dense route was written for does not
lose either, which was the thing to check: at 500 levels of d = 4 the two
are 3.7 ms and 2.4 ms, and at 500 levels of d = 10 they are 9.5 ms and
3.4 ms.

The density is IDENTICAL, not merely close. `test-sparsear1.R` keeps the
pre-v0.43 dense computation inline as a reference and holds the new one
against it at 20 random `(b, theta, d)` configurations, drawn wide enough
that |rho| reaches 0.999: values agree to 2.3e-12 relative and gradients
through `MakeTape` to 2.2e-09 relative, against gates of 1e-10 and 1e-08.
A third gate holds the same density against `RTMB::dautoreg` with
`scale` set to the marginal standard deviation, which is the independent
statement that the scale convention is the right one; that agrees to
1.1e-12.
The `vcov`, `start`, `npar`, `sd_idx`, `cor_spec` and `chol_L` accessors
are untouched, so the LKJ prior and the non-centering surface see the
same block they always did.

`ou`, `exp`, `gau` and `mat` stay dense on purpose: they are genuinely
dense kernels over arbitrary positions with no banded inverse to exploit.
`toep` and `homtoep` stay dense too, for the different reason that their
banded parameterization is not positive definite everywhere, so there is
no factor to exploit over the whole parameter space.

**2. A latent state-space term in the grammar.** Unlocks `mvrw`,
`mvrw_sparse` and `sde_linear` outright, and is the only route to
`thetalog` and `rickervalidation`. The prototype settles the hard half
of the question: the marginalization is already solved, exactly and at
competitive speed, by the Laplace machinery frmtmb runs today. What is
missing is a way to say "these coefficients are a random walk" (or a
nonlinear transition, for the last two) in a formula. This is design
work on the parser and one covariance-structure entry, not numerics.

**3. A unit-variance switch on the spatial and correlation structures.**
Turns `matern`, `transform` and `transform2` from pinned fits into plain
ones. Three of the 35 examples needed the identical
`lower = upper = c(theta_1 = 0)` ritual, and `matern` needed a second pin
on the residual sigma. A documented `nugget = FALSE` or `unit = TRUE`
argument would spell it once, and the pins are the only reason those
three rows carry an explanation instead of a formula.

**4. `expm`.** Unlocks `hmm`. This is a real atomic-function item and
the highest-ranked one, because it is not only about this example: the
pharmacometrics tier in `dev/feature-gaps.md` wants linear compartment
models, and those are the same matrix exponential. The example needs
`expm` of a 101 x 101 generator once per objective evaluation.

**5. A separable Kronecker combinator over index dimensions.** Unlocks
`ar1_4D`. frmtmb's existing Kronecker work (`gr_cov`, `gr_prec`, the
`|ID|` merged blocks) is Kronecker over coefficients within a level, not
over several independent index dimensions of one field. One example, and
no other demand recorded, so it ranks below the items above.

**6. An anisotropic option on `spde()`.** Unlocks `spde_aniso` and
`spde_aniso_speedup`. Two examples, and the change is contained: build
`Q_spde` with the `H` matrix of Lindgren et al equation 20 instead of the
isotropic form. Not an atomic function.

**7. Adaptive Gauss-Kronrod in `quadrature = TRUE`.** Not an unlock at
all, a refinement: it is the entire 1.36e-03 in the
`adaptive_integration` row. The spec built in `R/fit.R` sets
`adaptive = FALSE`; turning it on would match `stats::integrate`, which
is the reference's rule.

**8. A Kalman-filter atomic node.** Ranked here, not at the top, because
the measurement says it buys little: both routes are linear and agree
to 1e-13, the constant between them is host-dependent (within about a
factor of two in either direction across two machines), and the node
needs the state-space grammar of item 2 before it can be reached from
a formula at all. If it is ever built,
`MakeTape(...)$atomic()` is the way: 18276 outer nodes down to 9, with no
adjoint to derive.

**9. An FFT node.** Unlocks `fft` and nothing else in this set. A
circulant covariance is a narrow model class and no user has asked for
it.

**10. A root-solve node.** No example in either upstream directory needs
one. It stays on the list only because implicit parameterizations were
part of the original scoping; nothing here argues for it.

**Genuinely never.** `nmix` sums over a discrete latent abundance, which
the formula compiler has no representation for. `spa_gauss` marginalizes
by saddlepoint, a property of the tape's inner problem rather than of the
model. `linreg_parallel`, `transform_parallel` and
`register_atomic_parallel` need OpenMP accumulation, which RTMB does not
have. `sam` is a bespoke fisheries assessment with hard-coded fleet
types and a stock-recruitment switch; it is a program, not a formula.
The four API demonstrations have nothing to unlock.

## Two findings worth carrying into the docs

**The `nl` body resolves functions lexically.** `bf(y ~ qgamma(pnorm(z),
shape, scale), ..., nl = TRUE)` works when the user has run
`library(RTMB)` and fails with "The nonlinear formula body could not be
evaluated: Non-numeric argument to mathematical function" when they have
not, because the body is evaluated in the formula's own environment and
bare `pnorm` finds the `stats` version. Writing `RTMB::pnorm` and
`RTMB::qgamma` in the formula works either way, and gives the identical
fit (both spellings reach logLik = -544.1792742612 on the `transform`
data). Users reaching for a distribution function inside an `nl` body
should be told to prefix it, because the failure message points at the
formula rather than at the search path.

RESOLVED (lane wt-nlenv), the other way round: rather than tell users to
prefix, the body became a scope that sees RTMB's math first, so the bare
spelling is the one that works and the prefixed one still does.
`test-tmb-examples.R` keeps its `RTMB::`-prefixed formulas as the
compatibility case.

**`Vectorize()` over advectors needs RTMB attached.** Building a
per-observation marginalization with `Vectorize()`, as
`adaptive_integration.R` does, fails with "Invalid argument to
'advector' (lost class attribute?)" unless RTMB is on the search path,
because `mapply`'s simplification drops the class.
`do.call("c", lapply(...))` dispatches on the registered method and works
either way. The same trap catches base `matrix()`, which is already
recorded as an RTMB gotcha. This is why `tests/testthat/test-tmb-examples.R`
uses `RTMB::` prefixes throughout rather than attaching the package.

## Showcase nominations

Three REPLICATES rows would carry a vignette section, if one is ever
written. This is a nomination, not vignette content.

1. **orange_big.** Nonlinear growth with a random asymptote, 35000
   observations and 5000 latent variables, fitted in a few seconds and
   within 2.33e-08 of the C++ template. It is the clearest demonstration
   that `nl = TRUE`
   is a real nonlinear mixed-model surface rather than a convenience, and
   the three-parameter logistic curve reads well on a page.
2. **lr_test.** The example exists to show TMB's `map=`. frmtmb needs no
   `map=`: the restriction is a formula, `sigma ~ 0 + g` against
   `sigma ~ 1`, and `anova()` produces the same chi-square table. It is
   the best short argument for the formula-compiler thesis in `SPEC.md`.
3. **transform.** A latent AR(1) field pushed through `pnorm` and a gamma
   quantile, written as one `nl` formula. It shows that the nonlinear
   surface composes with a covariance structure, which nothing else in
   the frequentist ecosystem offers. It would have to carry the pin
   honestly, which is itself an argument for roadmap item 3.

## Reproducing this audit

The vendored upstream copies are gitignored. To refetch:

```sh
gh api repos/kaskr/RTMB/contents/tmb_examples --jq '.[].name'
gh api repos/kaskr/adcomp/contents/tmb_examples --jq '.[].name'
```

then fetch each file's `.content` and base64-decode it into
`dev/tmb-examples/` and `dev/tmb-examples-adcomp/` respectively.

* `Rscript dev/tmb-examples-check.R` runs the RTMB-port rows against
  their own reference scripts, one R subprocess per case. The cases
  `sdv_multi`, `sdv_multi_compact`, `transform` and `transform2` run at
  the examples' published sizes; `transform_n200` and `transform2_n200`
  rerun the shortened series the audit was first written at, so the
  published numbers stay reproducible.
* `Rscript dev/sparsear1-bench.R` reports the AR(1) block density's
  scaling, with the pre-v0.43 dense computation reconstructed inline.
* `Rscript dev/tmbex-adcomp-check.R` runs the adcomp-only rows at their
  published size, against reference likelihoods hand-written from the
  `.cpp` templates (the originals need a compiler, so they cannot be run
  the way the ports can).
* `Rscript dev/tmbex-kalman-prototype.R` runs the Kalman prototype and
  prints the table above.
* `tests/testthat/test-tmb-examples.R` is the regression form, and needs
  no vendored copy at all.

## A replication vignette, deferred on purpose

Maintainer decision (2026-09-03): once the top roadmap items land
(the linear-cost `ar1()`, which did in v0.43; the state-space grammar
term if it happens), a dedicated vignette replicates the relevant TMB examples
side by side, upstream spelling against the one-formula `frm()`
spelling with the logLik agreement shown. Written then rather than
now so the stochastic-volatility and state-space rows appear at their
published sizes instead of pinned or capped, and so the vignette
demonstrates capability rather than apologizing for gaps. The
showcase nominations above are its seed; the tests in
`test-tmb-examples.R` are its validation layer already.
