# The classic TMB examples: what frmtmb replicates today

Audit date: 2026-09-03. Branch `wt-tmbex`, at v0.41.0 (2619816).

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
1e-6 or better. Two rows agree less closely and both are explained below.

The harness is `dev/tmb-examples-check.R` (one R subprocess per case,
because the reference scripts `source()` their data files into the global
environment). The Kalman prototype is `dev/tmbex-kalman-prototype.R`.
The REPLICATES rows are also regression tests in
`tests/testthat/test-tmb-examples.R`, which downloads nothing: the data
generation and the reference negative log-likelihood are both inline.

## Verdict counts

| verdict | count |
| --- | --- |
| (a) REPLICATES | 13 |
| (b) AWKWARD | 3 |
| (c) OUT OF SCOPE | 19 |

Of the 19 out-of-scope rows, 3 are TMB API demonstrations whose
underlying model is another row that replicates, 3 are OpenMP parallel
builds of models that already replicate, and 13 need a capability frmtmb
does not have. (A fourth API demonstration, `dataeval`, is counted under
REPLICATES, because its model fits.)

## (a) REPLICATES

`|diff|` is the absolute difference between the reference maximized
log-likelihood and `logLik(frm(...))`.

| example | source | `frm()` spelling | `|diff|` |
| --- | --- | --- | --- |
| linreg | RTMB | `bf(Y ~ x)`, `family = gaussian()` | 4.19e-13 |
| dataeval | RTMB | `bf(y ~ x)`, `family = gaussian()` | 1.15e-11 |
| tweedie | RTMB | `bf(y ~ 1)`, `family = tweedie()` | 5.46e-12 |
| compois | RTMB | `bf(x ~ 1)`, `family = compois()` | 1.31e-10 |
| spatial | RTMB | `bf(y ~ 1 + x2 + exp(pos + 0 \| grp))`, `family = poisson()` | 1.01e-08 |
| spde | RTMB | `bf(time \| cens(cens) ~ sex + age + wbc + tpi + spde(fem, gr = node))`, `family = weibull()` | 9.35e-09 |
| adaptive_integration | RTMB | `bf(x \| trials(n) ~ 0 + c1 + c2 + c3 + (1 \| obs))`, `family = binomial()`, `quadrature = TRUE` | 1.36e-03 |
| transform | RTMB | `bf(y ~ qgamma(pnorm(z), shape, scale), z ~ 0 + ar1(tim + 0 \| g), shape ~ 1, scale ~ 1, nl = TRUE)` | 3.08e-09 |
| transform2 | RTMB | `bf(y ~ qbeta(pnorm(z), shape1, shape2), z ~ 0 + ar1(tim + 0 \| g), shape1 ~ 1, shape2 ~ 1, nl = TRUE)` | 6.37e-12 |
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
parameter, and it lands 0.135 log-likelihood units above the reference.
Pinning that parameter with `lower = c(theta_1 = 0)` and
`upper = c(theta_1 = 0)` is the `map=` equivalent and recovers the
reference exactly. The fitted `theta_2 = 0.5973` maps to
rho = 0.5973 / sqrt(1 + 0.5973^2) = 0.5128, which is the reference's phi
to four decimals. Both rows were run at n = 200 rather than the example's
n = 1000, for the reason in the AWKWARD section: `ar1()` is a dense
covariance structure and its cost is cubic in the series length.

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

### sdv_multi and sdv_multi_compact

The multivariate stochastic volatility model of Skaug and Yu (2014).
Three series, 945 time points, a latent AR(1) log-volatility per series,
and an observation density that is multivariate normal with mean zero, an
unstructured correlation matrix, and per-series scales
exp(0.5 (mu_x + h)).

This is expressible in frmtmb's grammar, and the mapping is exact rather
than approximate:

```r
frm(mvbf(bf(x1 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
         bf(x2 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
         bf(x3 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
         rescor = TRUE), data = dd)
```

The residual correlation of a multivariate gaussian fit is the example's
`R`; `sigma`'s log link makes its intercept mu_x / 2 and its latent block
h / 2; and h / 2 is AR(1) with the same phi, so `ar1()`'s marginal sd
absorbs the factor of 2 exactly. Nothing is approximated.

Measured on simulated data of the example's own shape, three series with
an AR(1) log-volatility each:

| series length | reference nll | `frm()` logLik | agreement | `frm()` fit time |
| --- | --- | --- | --- | --- |
| 40 | 133.3885404159 | -133.3885404157 | 2.23e-10 | 1.4 s |
| 80 | 259.2112861880 | -259.2038643997 | 7.42e-03 | 4.1 s |

At length 40 the two routes converge to the same point and agree
exactly. At length 80 they do not, and the difference is the reference's
own optimizer, not the model: `frm()` reaches a strictly better objective
than the reference does, and the 259.2113 above is already the best of
twelve reference starts (the example's own single start stops at
260.5668). The 12-parameter box-constrained problem is simply harder to
converge from a cold start than the likelihood identity is to establish.

So what blocks the row is cost, not grammar, and the next section is
about that. The fit-time column already shows it: 1.4 s at length 40,
4.1 s at length 80, against a reference that takes 0.2 s at both.

### The blocker the AWKWARD rows share: `ar1()` is dense

frmtmb's `ar1`, `toep`, `ou`, `exp`, `gau` and `mat` structures all build
a dense `d x d` covariance and hand it to `RTMB::dmvnorm`, where `d` is
the block dimension. That is the right shape for repeated measures, where
`d` is 4 or 10. It is the wrong shape for a time series, where `d` is the
number of time points, because the Cholesky factorization is taped in
full and the cost is cubic.

Measured, as tape construction plus one `fn()` and one `gr()` call:

| block dimension | `ar1()` dense | `RTMB::dautoreg` sparse |
| --- | --- | --- |
| 50 | 0.13 s | 0.02 s |
| 100 | 1.04 s | 0.02 s |
| 200 | 0.22 s | 0.02 s |
| 400 | 1.17 s | 0.05 s |
| 800 | 7.25 s | 0.09 s |
| 1000 | not measured | 0.09 s |

(The dense column is noisy below 400 because taping dominates; the trend
above it is the cubic one.)

At the examples' own sizes this is the difference between a fit and no
fit: `sdv_multi` has d = 945 in three blocks, and `transform` has
d = 1000. It is also why the two `transform` rows were verified at
n = 200.

The fix is not an atomic function. `RTMB::dautoreg` already exists, and
frmtmb already assembles a parameter-dependent sparse precision on the
tape and evaluates it with `RTMB::dgmrf`: that is exactly what the `spde`
covariance structure does, in `R/covstruct.R`. A sparse `ar1` entry in
`covstruct_registry` would follow the same pattern and reuse the same
`dgmrf` call. This is the highest-value item the audit found, and it is
ordinary work rather than new machinery.

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

**1. A sparse AR(1) covariance structure.** Unlocks `sdv_multi` and
`sdv_multi_compact` at their published size, and takes `transform` and
`transform2` from 200 usable time points to 1000. Today `ar1()` builds a
dense `d x d` covariance and Cholesky-factorizes it on the tape; at
d = 800 that is 7.25 s for one tape plus one gradient, against 0.09 s for
`RTMB::dautoreg`. The pattern to copy already exists in the same file:
`covstruct_registry$spde` assembles a parameter-dependent sparse
precision and evaluates it with `RTMB::dgmrf`. A tridiagonal AR(1)
precision is the same shape and simpler. No new machinery, four examples,
and every user with a long repeated-measures series benefits.

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
  their own reference scripts, one R subprocess per case.
* `Rscript dev/tmbex-adcomp-check.R` runs the adcomp-only rows at their
  published size, against reference likelihoods hand-written from the
  `.cpp` templates (the originals need a compiler, so they cannot be run
  the way the ports can).
* `Rscript dev/tmbex-kalman-prototype.R` runs the Kalman prototype and
  prints the table above.
* `tests/testthat/test-tmb-examples.R` is the regression form, and needs
  no vendored copy at all.

