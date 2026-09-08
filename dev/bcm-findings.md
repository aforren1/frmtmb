# Bayesian Cognitive Modeling port: findings

Lane: BCM. Worktree `frmtmb-wt-bcm`, branch `wt-bcm`, off `68d6782`
(frmtmb 0.52.0).

Target: the Stan ports of Lee and Wagenmakers, *Bayesian Cognitive
Modeling: A Practical Course* (2013), at
<https://github.com/stan-dev/example-models/tree/master/Bayesian_Cognitive_Modeling>
(BSD-3). The Stan programs are read as references and adapted; the
models themselves are the book's published equations.

## Inventory

The repository holds 71 model files under `Bayesian_Cognitive_Modeling`
once the `! other` scratch directory is excluded: 65 driver scripts
named `*_Stan.R` with the program in a string, and 6 standalone
`.stan` files for the Binomial chapter. Counting distinct models and
dropping the `GettingStarted` copy of `Rate_1`:

| chapter | models |
|---|---|
| ParameterEstimation/Binomial | Rate_1..5, Survey |
| ParameterEstimation/Gaussian | Gaussian, SevenScientists, IQ |
| ParameterEstimation/DataAnalysis | Correlation_1, Correlation_2, Kappa, ChangeDetection, Planes, ChaSaSoon |
| ParameterEstimation/LatentMixtures | Exams_1, Exams_2, TwoCountryQuiz, TwentyQuestions, Cheating, Malingering_1, Malingering_2 |
| ModelSelection/Means | OneSample, OneSampleOrderRestricted, TwoSample |
| ModelSelection/Rates | Pledgers_1, Pledgers_2, Geurts, GeurtsOrderRestricted, Zeelenberg |
| CaseStudies/PsychophysicalFunctions | PsychophysicalFunction1, PsychophysicalFunction2, PsychometricFunction1_Answers |
| CaseStudies/SignalDetection | SDT_1, SDT_2, SDT_3 |
| CaseStudies/MemoryRetention | Retention_1..3 |
| CaseStudies/Simple | SIMPLE_1, SIMPLE_2 |
| CaseStudies/MPT | MPT_1..5 |
| CaseStudies/GCM | GCM_1, GCM_2, GCM_3, GCM_3_optimized |
| CaseStudies/BART | BART_1, BART_2 |
| CaseStudies/ESP | Ability, Extraversion, OptionalStopping, OptionalStopping_Answer_1 |
| CaseStudies/NumberConcepts | NumberConcept_1..3 |
| CaseStudies/HeuristicDecisionMaking | Search, Stop, SearchStop, TTB |

## Counts

Distinct models in the reference repository, excluding the
`GettingStarted` copy of `Rate_1` and the `! other` scratch directory:
63.

| | count |
|---|---|
| ported and fitted | 52 |
| of those, carrying a Stan identity row of their own | 41 |
| of those, validated by structure or recovery only | 11 |
| refused by name (no free parameter) | 1 (Planes) |
| not ported | 10 |

The identity tier makes **39 `stan_lp_check()` calls over 38 distinct
Stan programs**. Counted from the per-chapter reports with
`options(frmtmb.bcm_report = TRUE)`: binomial 6, gaussian 3,
data-analysis 6, latent-mixtures 5, model-selection 3, psychophysics 2,
signal-detection 3, retention 3, mpt 2, gcm 2, simple 1, bart 1, esp 2.
Two of those numbers do not line up with the model count and both have
a reason: SDT_3 is a second program checked against SDT_2's own fit, so
39 calls run over 38 programs; and one row covers MPT_2 through MPT_5,
which are one model. Netting those, 41 of the 52 fitted models carry an
identity of their own and 11 do not (Cheating, Pledgers_2,
GeurtsOrderRestricted, Zeelenberg, OneSample,
OneSampleOrderRestricted, PsychometricFunction1_Answers, SIMPLE_2,
BART_2, Ability, TwoCountryQuiz).

An earlier draft of this file gave three different numbers for this one
quantity (34 here, "thirty-eight" in the identity section, 39 rows in
the table). The review of 2026-09-06 caught it; the measured value is
the one above.

The 10 not ported are NumberConcept_1 to _3, the four
HeuristicDecisionMaking models (Search, Stop, SearchStop, TTB), GCM_3
and GCM_3_optimized, and ESP's OptionalStopping_Answer_1 (an exercise
answer on the model already ported). None is blocked by a missing seam:
each is a `loglik` slot or a fit of a shape already in the port.

Four of the 52 collapse onto one fit each with a sibling, and the table
counts them separately because the book does: MPT_2 through MPT_5 are
one hierarchical fit, SDT_3 is SDT_2's, Pledgers_2 is Pledgers_1's, and
OneSampleOrderRestricted is OneSample's.

Delivered: 13 `test-bcm-*.R` files with 91 `test_that` blocks over 3469
lines, `helper-stan.R` (283 lines) and `helper-bcm.R` (122), four
family files under `inst/bcm/` (1405 lines), and the vignette (475).

## Early corrections to the desk assessment

- **ChaSaSoon is not a marginalization.** Its Stan program adds
  `nfails * log(binomial_cdf(25, n, theta) - binomial_cdf(14, n, theta))`,
  which is 949 interval-censored binomial observations sharing one
  interval. That is `cens()` plus `weights()` in the grammar, not the
  whole-response slot.
- **Planes and Survey marginalize a discrete PARAMETER, not a latent
  state per row.** Planes sums a hypergeometric over the population
  size; Survey sums a binomial over the number of surveys sent. Both
  are finite sums over one scalar with a discrete uniform prior, and
  both have no regression structure at all.
- **ChangeDetection's Stan port declares `tau` continuous** with an
  `if` on `t[i] - tau`, so its own log density is piecewise constant in
  `tau` and its gradient is zero almost everywhere. The book's model
  puts a discrete uniform on the changepoint. The port here
  marginalizes the discrete changepoint, which is the finite sum the
  loglik slot exists for.

## The identity harness

`tests/testthat/helper-stan.R`. It generalizes the mechanism written
out in place in `helper-rl.R`: a test hands over a Stan program, its
data and a map from frmtmb's estimate to its parameters, and gets back
the assertion that `log_prob` at that point equals frmtmb's log density
plus a stated constant, with a gradient check.

Three things it does that `rl_lp_check()` did not:

- **Names the missing block.** `stan_par_blocks()` reads Stan's own
  `unconstrained_param_names()`, so a map that misses a parameter is
  refused by name instead of surfacing as an rstan length complaint.
- **Indexes the inner block exactly.** `stan_upar_index()` uses the
  same list, which is correct for a simplex or a Cholesky factor where
  the unconstrained length is not the constrained one. The
  bump-and-diff trick `helper-rl.R` uses would break such a constraint.
- **Round-trips the map.** `constrain_pars(unconstrain_pars(pars))`
  must return what went in, which is the unit test on the map alone.

Program compilation is shared with the brms suite through
`brms_stan_model()`, deliberately: that function owns the session cache
that works around a DSO which cannot be re-read from its own RDS in the
session that wrote it, and two caches over one directory would hit it.
`helper-rl.R` already reaches across for the same function, so the
arrangement is precedented.

### The constant rule, as written in the helper

1. The Stan program carries a prior frmtmb does not, and the constant
   is that prior's log density at the estimate. A flat prior on a
   bounded parameter is minus the log of the width, zero only when the
   width is one, and only when the program WRITES it: a bare
   `<lower, upper>` declaration adds nothing to `target`, so on its own
   it contributes nothing.
2. frmtmb carries the same prior through `set_prior()`, and the
   constant is 0. `RTMB::dnorm(..., log = TRUE)` is the full normalized
   density, so a matched `normal()` prior cancels exactly.
3. Neither carries a prior, and the constant is 0.

Every reference program is adapted the same way throughout:
`y ~ dist(...)` becomes `target += dist_lpdf(y | ...)`, because Stan
drops the normalizing constant of a `~` statement and frmtmb's families
do not. Left alone, that would put `lchoose(n, k)` into the constant
and hide a real difference behind it.

### The gradient rule

`grad = "all"` when the two programs hold the same function;
`grad = "inner"` for a fit with random effects, where frmtmb's outer
estimate maximizes the Laplace-approximated marginal and only the inner
block is at a stationary point of the joint; `grad = "none"` with a
stated reason. The tolerance is 1e-3, matching `helper-brms.R`: it is
the scale at which frmtmb's optimizer declares convergence, so a
tighter one tests the optimizer rather than the map. A run with 1e-4
failed on SDT_1 at 1.4e-4, which is the optimizer and not the map.

## The identity, measured

Every row is one `stan_lp_check()` call, `options(frmtmb.bcm_report =
TRUE)`, R 4.6.1 with rstan 2.32.7 and StanHeaders 2.39.1. "residual" is
`measured - stated`; the harness passes it when
`|residual| < 1e-6 max(1, |ours|)`. "grad" is the largest absolute
gradient over the coordinates the row claims, and the tolerance is
1e-3.

| chapter | model | stated const | measured const | residual | grad | ours |
|---|---|---|---|---|---|---|
| binomial | Rate_1 | 0 | -6.2e-15 | -6.2e-15 | 0 | -1.402042718 |
| binomial | Rate_2 | 0 | -1.2e-14 | -1.2e-14 | 5.5e-10 | -2.723193996 |
| binomial | Rate_3 | 0 | -1.3e-14 | -1.3e-14 | 0 | -3.14331251 |
| binomial | Rate_4 | 0 | -4.4e-16 | -4.4e-16 | 2.2e-16 | -0.9659002008 |
| binomial | Rate_5 | 0 | 0 | 0 | 0 | -13.86294361 |
| binomial | Survey | 0 | -7.3e-14 | -7.3e-14 | 3.1e-07 | -15.87642482 |
| gaussian | Gaussian | -log(10) = -2.302585093 | -2.302585093 | 0 | 3.2e-06 | -6.695514875 |
| gaussian | SevenScientists | 0 | -3.6e-15 | -3.6e-15 | 5.7e-06 | -25.02593406 |
| gaussian | IQ | -3 log(300) - log(100) = -21.71651761 | -21.71651761 | 0 | 1.1e-06 | -25.43079502 |
| data-analysis | Correlation_1 | 0 | 0 | 0 | 8.3e-05 | -35.09918872 |
| data-analysis | Correlation_2 | 0 | -7.1e-15 | -7.1e-15 | 2.0e-13 (inner) | -16.31365195 |
| data-analysis | Kappa | 0 | -1.1e-13 | -1.1e-13 | 1.0e-06 | -5.565670845 |
| data-analysis | ChangeDetection | 0 | 5.7e-14 | 5.7e-14 | 3.4e-04 | -164.2869216 |
| data-analysis | Planes | 0 | -3.6e-15 | -3.6e-15 | none (no parameter) | -2.500529574 |
| data-analysis | ChaSaSoon | 0 | -3.7e-13 | -3.7e-13 | 9.1e-08 | -118.1411391 |
| latent-mixtures | Exams_1 | 0 | -1.1e-13 | -1.1e-13 | 4.6e-12 | -44.1104793 |
| latent-mixtures | Exams_2 | 0 | -1.1e-13 | -1.1e-13 | 6.6e-11 (inner) | -27.12901174 |
| latent-mixtures | Malingering_1 | 0 | -4.7e-13 | -4.7e-13 | 6.6e-06 | -76.97499896 |
| latent-mixtures | Malingering_2 | 0 | -2.18e-05 | -2.18e-05 | 2.5e-05 | -58.7146357 |
| latent-mixtures | TwentyQuestions | 0 | 5.7e-14 | 5.7e-14 | none (boundary) | -53.97133923 |
| model-selection | Pledgers_1 | 0 | -1.3e-12 | -1.3e-12 | 1.8e-05 | -8.312756865 |
| model-selection | Geurts | 0 | -1.1e-12 | -1.1e-12 | 5.9e-10 (inner) | -238.0004582 |
| model-selection | TwoSample | 0 | 0 | 0 | 1.5e-05 | -65.38884702 |
| psychophysics | PsychophysicalFunction1 | 0 | 2.3e-13 | 2.3e-13 | 3.8e-10 (inner) | -241.536381 |
| psychophysics | PsychophysicalFunction2 | 0 | 2.3e-13 | 2.3e-13 | 8.8e-09 (inner) | -186.8521235 |
| signal-detection | SDT_1 | 0 | -2.3e-13 | -2.3e-13 | 1.4e-04 | -23.30327956 |
| signal-detection | SDT_2 | 0 | 0 | 0 | 1.7e-10 (inner) | -139.2786419 |
| signal-detection | SDT_3 | 0 | 0 | 0 | 1.7e-10 (inner) | -139.2786419 |
| retention | Retention_1 | 0 | 2.8e-14 | 2.8e-14 | 8.5e-07 | -65.54760139 |
| retention | Retention_2 | 0 | 7.1e-14 | 7.1e-14 | 9.6e-06 | -38.05837414 |
| retention | Retention_3 | 0 | 6.0e-14 | 6.0e-14 | 4.5e-07 (inner) | -11.40692598 |
| mpt | MPT_1 | 0 | -2.8e-13 | -2.8e-13 | 2.4e-07 | -8.291381644 |
| mpt | MPT_2 to MPT_5 | 0 | 5.7e-14 | 5.7e-14 | 2.1e-10 (inner) | -92.94613884 |
| gcm | GCM_1 | 0 | -4.3e-14 | -4.3e-14 | 1.7e-06 | -36.22950444 |
| gcm | GCM_2 | 0 | -1.3e-13 | -1.3e-13 | 5.8e-05 | -122.839191 |
| simple | SIMPLE_1 | 0 | 6.0e-11 | 6.0e-11 | 3.1e-04 | -2083.454577 |
| bart | BART_1 | 0 | 2.8e-14 | 2.8e-14 | 1.2e-05 | -131.083217 |
| esp | OptionalStopping | 0 | 3.6e-15 | 3.6e-15 | 2.2e-07 | -30.37952411 |
| esp | Extraversion | 0 | -2.9e-12 | -2.9e-12 | 1.6e-10 (inner) | -1.747339859 |

Thirty-nine `stan_lp_check()` calls over 38 distinct Stan programs:
SDT_3 reuses SDT_2's fit against a second program, and MPT_2 to MPT_5
share one row because they are one model. Three rows are worth a
sentence.

**SDT_2 and SDT_3 report the same `ours` to ten digits (-139.2786419)
and the same gradient.** That is the claim of the chapter's third
model made numerically: parameter expansion is a reparameterization,
and at `xi = 1` its program is the same function as SDT_2's.

**Malingering_2's residual is -2.18e-05, four orders larger than every
other row, and it is the symptom of a fourth unbounded-under-ML model.**

An earlier draft of this file called it a normalizing-constant
convention between `beta_binomial_lpmf` and RTMBdist. That was wrong,
and the review of 2026-09-06 caught it: a convention would appear in
BOTH components and would not depend on the shapes. Decomposed against
a longhand `lchoose(N, k) + lbeta(k + a, N - k + b) - lbeta(a, b)` and
reproduced independently in `scratchpad/bc-f1-phi2.R`:

- component 1, at `phi1 = 7.52`, agrees with the longhand to 1e-14 and
  its density sums to 1.000000000000 over 0..45;
- component 2, at `phi2 = 7.19e8`, sums to 1.0000022 (the longhand formula gives 0.999999998), and the
  whole residual is there.

It is precision loss in the `lgamma` differences at DIVERGENT shapes,
not a convention. And the reason the shapes are divergent is the
finding. Profiling the exact mixture log-likelihood in `phi2` with
everything else held at the estimate:

```
phi2 = 100         -59.4917902766
phi2 = 1e+04       -58.7234875260
phi2 = 1e+06       -58.7147459055
phi2 = 1e+08       -58.7146583610
phi2 = 7.19e+08    -58.7146576228   <- where the optimizer stopped
phi2 = 1e+10       -58.7146572816
phi2 = Inf         -58.7146574774   (the binomial limit)
```

Monotone up to the binomial limit; the last two differ by 2e-7, which
is floating-point noise on a surface flat to machine precision.
**The maximum likelihood estimate of `phi2` is infinity.** The second
group has no person-level overdispersion at all. It is a plain
binomial, and the beta-binomial's second parameter is buying nothing.
The standard error `summary()` prints for `log(phi2)` is a Hessian on a
flat direction and means nothing.

Confirmed from inside the grammar rather than only by profiling:
`mixture(beta_binomial, binomial)`, the same model with component 2 AT
the limit, reaches -58.71466 with four parameters against -58.71464
with five, and the same `mu = (0.509, 0.993)`. The 2.2e-5 between them
is that same precision loss, with the sign that makes a divergent
beta-binomial look fractionally better than the exact limit it is
approaching.

The test now asserts this instead of asserting that the precisions are
"finite". `fixef()` returns `log(phi)`, so the old assertion passed at
`log(phi2) = 20.39` and could not have failed for the reason its
comment gave.

**SIMPLE_1 needed a tighter optimizer, not a looser tolerance.**
nlminb's default stopping rule left the gradient at 1.2e-3, over the
harness's 1e-3. The log-likelihood is -2083.455 either way and
`rel.tol = 1e-12` brings the gradient to 3.1e-4, so it was the
optimizer's rule and not the map. Eighteen free parameters over 135
rows, each of which reads its whole list, is a flatter surface than the
default rule was chosen for.

Two rows claim the value only, each with a reason in the test:
TwentyQuestions, because person 8 answered nothing and their rate is
driven to a boundary the optimizer stops at rather than zeroes; and
Planes, because the program has no parameters at all.

## Gaps found in the core

Two, both worked around from `inst/bcm/` without touching `R/`.

1. **No probit link.** `frmtmb_links` has identity, log, logit,
   cloglog, inverse, logm1, tan_half and power12; brms has probit,
   probit_approx and cauchit. Signal detection theory is defined on the
   normal scale and so are the multinomial processing trees' latent
   traits, so the port needs one. `get_link()` accepts a LIST with
   `name`, `linkfun`, `linkinv` and `mu_eta`, so `bcm_probit()`
   supplies it in four lines and no seam is missing. Worth adding to
   the core registry all the same; that belongs to whoever owns
   `R/links.R`.
2. **No CDF on the binomial.** `cens()` refuses the core binomial by
   name and its message says how a family supplies one, which is
   exactly what `bcm_binomial_cdf()` does with `RTMB::pbinom`. This is
   the seam working as designed rather than a gap, and the test asserts
   the refusal message as well as the fix.

One seam found and used that `dev/structured-family-protocol.md` does
not document: `dpars[[".eta_<name>"]]` carries the LINEAR PREDICTOR
behind a dpar, which is what the core's own `robust_logit()` reads.
`bcm_contaminant()` uses it to recover `log(phi)` and `log(1 - phi)`
exactly however far the predictor runs; without it a saturating link
turns `log(1 - phi)` into minus infinity where the true value is
finite, and the optimizer walks into it. The protocol document should
say that a custom family may read it.

## Three more gates in the core, found by running into them

None of them is a missing seam. Each is a check keyed on something a
custom family cannot change, and each has a route around it from
`inst/bcm/`.

1. **`se()` is gated on the family NAME.** "se() is supported for
   gaussian and student families only". A custom family that reads
   `aterms[["se"]]` faithfully cannot opt in. ESP's Extraversion needs
   a known measurement standard deviation on a probit gaussian, so the
   port routes it through `vreal(sd)` instead, which is the channel a
   custom family IS given. One word in the formula, one line in the
   family.
2. **`cens()` is gated on family TYPE, not only on a CDF.** The first
   gate, "cens()/trunc() need a family with a CDF", names the seam and
   `bcm_binomial_cdf()` passed it. The second, "cens() is not supported
   for discrete families yet (truncation is)", refuses every discrete
   family whatever it supplies. ChaSaSoon therefore carries its band in
   its own likelihood: the response is the smallest count consistent
   with the row, `vint(hi)` the largest, and
   `log(F(hi) - F(y - 1))` is the exact density when the two agree and
   the band probability otherwise, with no branch. `weights()` carries
   the repeats. The test asserts BOTH refusal messages, so the day
   either gate moves the test says so.
3. **`mixture()` starts a bounded-response mixture at its exchangeable
   point.** Each component's `mu` is initialized from a quantile of the
   response, which for a binomial or a beta-binomial is a COUNT and
   falls outside the logit link's range; the warning says so and both
   components then start at the link's origin, where the two components
   are the same distribution. A binomial mixture escapes;
   `mixture(beta_binomial, beta_binomial)` does not, and comes back
   with the components identical, every latent probability at one half,
   and no error. Measured on Malingering_2: from the default start
   mu1 = mu2 = 0.801 and the log-likelihood is worse than from a
   separated start, which reaches mu = (0.51, 0.99). Every free mixture
   in the port is started from two separated rates.

Item 3 is the one worth acting on in core: a silently degenerate fit is
worse than a refusal, and `init_dpars` for a mixture over a bounded
family could put the components at two quantiles of `y / trials` rather
than of `y`.

### Reproductions

> **SEAM 1 IS CLOSED.** The links lane put probit in the core registry
> (`R/links.R`), so the transcript below is history for that row only.
> `binomial(link = "probit")` now works, and the workaround this file
> recorded has been retired: `bcm_probit()`, `bcm_binomial_probit()`
> and `bcm_gaussian_probit()` are deleted from
> `inst/bcm/binomial-extras.R`, and the models they served are core
> families. `bcm_binomial_probit()` became `binomial(link = "probit")`
> and `bcm_gaussian_probit()` became `gaussian(link = "probit")` with
> the known measurement SD on `se(sd)`.
>
> Closing seam 1 also removed this file's ONE use of the seam 2
> workaround, without touching seam 2. `se()` is still gated on the
> family name and `R/frame.R` is unchanged; the Extraversion model
> simply stopped being a custom family, so the gate no longer applies
> to it. Seams 2, 3 and 4 stand exactly as recorded.
>
> The Stan identities held. Residuals over
> `test-bcm-signal-detection.R`, `test-bcm-binomial.R` and
> `test-bcm-esp.R`, before and after, are in
> `dev/links2-findings.md`.

`scratchpad/bc-seams.R` and `bc-seam2.R` run all four. Transcript, R
4.6.1:

```
SEAM 1  no probit link                          CLOSED, see below
  frm(y | trials(N) ~ x, family = binomial(link = "probit"), data = d)
  ERROR: Unknown link: 'probit'. Available links: identity, log, logit,
         cloglog, inverse, logm1, tan_half, power12
  same call with bcm_binomial_probit():  NO ERROR

SEAM 2  se() gated on the family NAME           R/frame.R:1351
  se_reader <- frmtmb_family("se_reader", dpars = "mu",
    links = list(mu = "identity"), type = "continuous",
    lpdf = function(y, dpars, aterms)
      RTMB::dnorm(y, dpars[["mu"]], aterms[["se"]], log = TRUE),
    init_dpars = list(mu = function(y, aterms) mean(y)))
  frm(v | se(s, sigma = FALSE) ~ 1, family = se_reader, data = d)
  ERROR: se() is supported for gaussian and student families only
  the same term on the core gaussian:    NO ERROR
  workaround, vreal(sd):                 NO ERROR

SEAM 3  cens() refuses a discrete family        R/frame.R:1309 then :1319
  core binomial:      ERROR: cens()/trunc() need a family with a CDF
                             (currently: gaussian, lognormal, poisson,
                             exponential, weibull, inverse.gaussian, cox) ...
  bcm_binomial_band() (which DOES supply lcdf), WITH its vint(hi):
                      ERROR: cens() is not supported for discrete
                             families yet (truncation is)
  workaround, band written into the family:     NO ERROR

  The `vint(hi)` is load-bearing in that second line and the first
  draft of this transcript omitted it. Without it the family's own
  `required_aterms = "vint1"` fires first ("the density needs
  `vint1`") and the discrete gate is never reached, so the transcript
  would have shown the right conclusion by the wrong route. Supplying
  the term isolates the gate; so would a discrete family that carries
  `lcdf` and requires no aterm.

SEAM 4  mixture() starts at the link origin     R/families.R:2917, warned R/fit.R:1914
  frm(bf(k | trials(n) ~ 1),
      family = mixture(beta_binomial, beta_binomial), data = malingering)
  warning: Starting value 30 for mu1 is NaN through its logit link, so it
           was ignored and that intercept starts from zero on the link
           scale. Give init_dpars a value inside the link's range, or
           pass start =
  default start   mu = (0.801336, 0.801336)  logLik -64.25533
                  every latent probability exactly one half
  separated start mu = (0.509,    0.9926)    logLik -58.71464
```

The fourth is the one to act on in core. A silently degenerate fit that
returns identical components and a flat posterior, with only a warning
about a starting value, is worse than a refusal.

### What the trigger actually is

An earlier draft of this file said the trigger was a bounded response.
It is not, and the review of 2026-09-06 supplied the discriminating
experiment. The degeneracy needs EVERY component's `mu` init to fall
outside its link's range, so that all of them fall back to the same
origin and the mixture starts exactly at its exchangeable point, where
the gradient in the separating direction is zero.

A zero-heavy `mixture(poisson, poisson)` hits the same `init` bug and
warns identically, and ESCAPES. Measured in `scratchpad/bc-seam4.R` on
`c(rpois(60, 0.2), rpois(60, 9))`:

```
quantile(y, 1/3) = 0   log -> -Inf     spoiled
quantile(y, 2/3) = 8   log ->  2.079    usable
warning: Starting value 0 for mu1 is -Inf through its log link ...
fit -> mu = (0.2485, 8.847)   logLik -262.60
```

Only `mu1` is spoiled, `mu2`'s quantile is still usable, the components
stay separated, and the fit is fine. So the defect is "all components
collapse to one start", not "the response is bounded".

### What the fix has to be

`R/families.R:2917` takes `stats::quantile(y, k / (K + 1))` as each
component's starting mean. Taking the quantile of `y / size` instead is
necessary and NOT sufficient, and an earlier draft of this file
recommended it unqualified. That recommendation was wrong and would
have broken the very model that motivates the change. Measured on
Malingering_2's data, where 8 of the 22 respondents scored 45 of 45:

```
current    quantile(y, 1/3)      = 30          qlogis -> NaN
current    quantile(y, 2/3)      = 45          qlogis -> NaN
proposed   quantile(y/n, 1/3)    = 0.6666667   qlogis -> 0.6931472
proposed   quantile(y/n, 2/3)    = 1.0         qlogis -> Inf
starting there: ERROR: The optimizer failed on this model
  (mixture(beta_binomial, beta_binomial), ML, nlminb):
  NA/NaN gradient evaluation.
```

The quantile has to be CLAMPED as well, exactly the way each component
family's own `init_dpars$mu` already clamps at
`R/families.R:2150-2155`, `min(max(p, 0.02), 0.98)`. With the clamp it
lands on the separated-start optimum:

```
clamped start (0.6666667, 0.98) -> mu = (0.5091331, 0.992604)
                                   logLik -58.71465
```

So the core change is: take the quantile on the response divided by
`size` where the family has one, and pass it through the same clamp the
component family would have applied to its own init. The edit itself
belongs to whoever owns `R/families.R`; this lane owns only the
diagnosis and `scratchpad/bc-seam4.R`, which reproduces all three
states (degenerate default, unclamped failure, clamped success) plus
the poisson counter-example in one run.

## The capped likelihoods need a feasible start

`min(1, u)` is exact and taped (`bcm_cap1()`), and it turns an
UNBOUNDED problem into a bounded one. It also creates a hard edge: at
theta = 1 a row with fewer successes than trials has a log density of
-Infinity and a gradient of NaN, and the optimizer stops rather than
backing off. So both capped models need a start inside the feasible
region:

- Retention: the default nonlinear start of zero for every parameter
  puts alpha = beta = 1/2, and `exp(-t/2) + 1/2` is above one at the
  shortest lag. Started at the book's neighbourhood instead.
- SIMPLE: a threshold of 0.3 puts a ten-item list over the cap. The
  family's own default is 0.6, and every start at or above it reaches
  the same optimum (-2083.455 for SIMPLE_1); SIMPLE_2 needs its
  intercept at 1 rather than 0.75, because the threshold there is a
  linear function of list length and the intercept has to carry every
  list.

## Models that maximum likelihood cannot fit

Four, and each one is a property of the model rather than of frmtmb.
The fourth was found by the review of 2026-09-06 rather than by this
lane; see the Malingering_2 note under "The identity, measured" for the
profile that establishes it.

- **SevenScientists.** Profiling out the seven standard deviations
  leaves minus the sum of `log|x_i - mu|`, less 7/2, which diverges
  logarithmically as `mu` approaches any observation. There is no
  interior maximum and no local one: between two observations the
  profile has a MINIMUM. The book's `gamma(.001, .001)` on the
  precisions is what makes the Bayesian model proper. Ported as a MAP
  fit with a `normal(0, 1)` on the log standard deviations, carried on
  both sides so the constant stays zero. The test asserts the
  divergence itself rather than describing it.
- **PsychophysicalFunction2's contaminant, left free.** One rate per
  cell with a `beta(1, 1)` prior fits every cell exactly, so the
  mixture likelihood is maximized at a contamination rate of one and
  the psychometric function is not estimated at all. The rate is
  integrated out instead: a binomial rate with a `beta(1, 1)` prior is
  the discrete uniform on 0..n, so the contaminant component
  contributes a flat minus `log(n + 1)`.
- **Retention without its cap.** `exp(-alpha t) + beta` can exceed one,
  and a cell where every item was recalled then has a likelihood of
  `theta^n` with nothing holding `theta` down. Measured: an uncapped
  Retention_2 reached a fitted rate of 1.18 and the identity missed by
  4.3 nats. The cap is carried exactly through
  `min(1, u) = u - (|u - 1| + (u - 1)) / 2`, which is taped because
  `abs()` is, and Retention_2 and Retention_3 then match their programs
  to the last bit.

- **Malingering_2's second precision.** `phi2`'s profile is monotone up
  to the binomial limit, so its maximum likelihood estimate is
  infinity: the bona fide group has no person-level overdispersion and
  the beta-binomial's second parameter is unbounded above. Unlike the
  other three this one does not stop the fit, because the surface is
  flat to machine precision long before the parameter matters, so
  frmtmb returns `phi2 = 7.19e8` and a log-likelihood correct to 2e-5.
  It is still an unbounded direction and the test says so. Cheating,
  which is the same model on different data, is not affected: its two
  precisions are 8.9 and 5.8.

The same closed-form integral that rescues the contaminant model
rescues two more: Malingering_2 and Cheating give every person a
binomial rate with a beta prior, which integrates to a beta-binomial,
so both are `mixture(beta_binomial, beta_binomial)` in frmtmb's own
mean-precision parameterization. That integral is what makes them
estimable at all; it does not make `phi2` bounded.

## The MPT covariance block has no finite standard error

The hierarchical multinomial processing tree is the one fit in the port
whose `diagnose()` is not clean, and the reason is worth writing down
rather than tuning away.

Twenty-one subjects, twenty word pairs each, four categories. Four
categories with a fixed total is THREE free counts per subject, and the
model gives each subject THREE latent traits. An unstructured three by
three covariance over those traits is six more parameters, and 21
subjects is not many to estimate six covariance parameters from.

What is stable:

- the log-likelihood, -102.47794, from five optimizer settings (nlminb
  default, nlminb with `rel.tol = 1e-12`, `optimizer = "optim"`,
  `restarts = 3`, and tightened-with-restarts);
- the estimates: standard deviations 0.517, 0.464, 0.296 and
  correlations -0.529, 0.539, -0.450, which are the case study's own
  latent-trait correlations;
- the identity against MPT's Stan program, constant 5.7e-14 with an
  inner gradient of 2.1e-10.

What is not:

```
Optimizer convergence code: 1 (false convergence (8))
Max |gradient|: 65488633 at theta_3
Hessian positive definite: FALSE
Non-finite standard errors: theta_4, theta_6
min_cov_eigenvalue: -0.2403
```

Every parameter without a finite standard error is a COVARIANCE
parameter; the three group means and the traits themselves are fine.

The fit is NOT saturated, which was the first thing to rule out: three
traits against three free counts per subject would reproduce every
subject exactly, and the largest gap between a fitted and an observed
cell proportion is 0.177. The random effects are pooled and shrinkage
is working. It is the curvature of the covariance block, not the
location of the estimate, that the data do not pin down.

The case study meets the same wall from the Bayesian side. Its own
README records divergent transitions and correlation matrices that are
not positive definite during estimation, and MPT_2's Wishart prior and
MPT_3's LKJ prior are exactly what hold this block down. Maximum
likelihood has no such prior, so it reports the degeneracy instead, and
`tests/testthat/test-bcm-mpt.R` asserts every line of it: both
warnings, `pdHess == FALSE`, a negative minimum eigenvalue, that the
bad standard errors are all `theta`, and that the worst gradient is
`theta_3`. If a later release fixes or changes any of that, the test
says so.

### Warnings in the suite

Before this was written down the port raised eight warnings in the full
suite: two from the MPT fit, two from Retention_3, three from SIMPLE,
and one pre-existing in `test-prior-compat.R`. Each of the seven is a
real property of a hard likelihood, so each is now either asserted (the
MPT pair) or suppressed at its call site with the reason in a comment
(Retention_3's `NA/NaN function evaluation`, which is nlminb probing
past `bcm_cap1()`'s edge and backing off, and SIMPLE's
`false convergence (8)`, which is nlminb's stopping rule on a flat
surface whose gradient is 3.1e-4).

The whole package suite is now FAIL 0, WARN 1, SKIP 217, PASS 4392, and
the one warning is the pre-existing `test-prior-compat.R` one.

## Sampler devices that are not models

Recorded because they collapse the model count.

- **MPT_2 through MPT_5 are one model.** Parameter expansion with a
  Wishart prior, an LKJ prior on a correlation matrix, its Cholesky
  factor, and a non-centered reparameterization. The directory's own
  README says so. One frmtmb fit, `1 + (1 | p | id)` on each of three
  probits.
- **SDT_3 is SDT_2.** Parameter expansion again: `xi` and the raw scale
  are a ridge and only their product is identified, so under maximum
  likelihood SDT_3 has no separate estimate. The identity is asserted
  against SDT_3's own program at `xi = 1`, which is the point of the
  ridge where the two programs are the same function, and it holds.
- **Zeelenberg's "Matt trick" variant** is the same likelihood as its
  centred sibling.
- **`deltaprior` in every model-selection program** is a prior draw
  with no data attached, there so that the sampler produces the
  denominator of a Savage-Dickey ratio.

## The Bayes factors

Every model of chapters 8 and 9 fits, and none of their Bayes factors
is here. Each is a Savage-Dickey density ratio read off the posterior
of one parameter at one point, and a maximum likelihood fit has no
posterior to read.

What `tests/testthat/test-bcm-model-selection.R` does instead is fit
the models and pin the quantity the ratio is taken over: `delta` for
Pledgers, Geurts, Zeelenberg, OneSample and TwoSample. When
`frmtmb.sample` gains the ratio, the evidence test is one call per
model against those fits. The order-restricted variants need nothing
extra: an inequality constraint is not a term of the likelihood, and in
every case here the unrestricted maximum already satisfies it, so the
restricted maximum IS the unrestricted one.

## Refused, by name

One model and one spelling. An earlier draft of this file ran the
refusals and the unported work together in one table, which the review
of 2026-09-06 read as calling unfinished work refused. They are split
here and in the vignette.

| what | why |
|---|---|
| Planes, the model | The model has NO free parameter: its own Stan port runs with `algorithm = "Fixed_param"`. `frm()` estimates parameters, so there is no formula to write. The posterior over the fleet size is computed directly by `bcm_planes_posterior()` and its arithmetic is checked against the reference program with an empty parameter vector. |
| Correlation_2 through `mi(sd)`, the spelling | `mi()` refuses `rescor = TRUE`, and with a correlated random intercept instead there are two variances per response that one observation per subject cannot separate. The MODEL is not refused: `se(sd, sigma = FALSE)` plus `(1 \| p \| id)` is exactly it. |

## Not ported

Ten models, none blocked by a missing seam: each is a `loglik` slot of
a shape `inst/bcm/marginal.R` already uses, or a fit of a shape already
in the port. This is unfinished work, not an obstacle.

| models | count | note |
|---|---|---|
| NumberConcept_1 to _3 | 3 | a finite sum over a per-child knower level, but the level indexes a different response distribution in each of two tasks |
| Search, Stop, SearchStop, TTB | 4 | heuristic decision making |
| GCM_3, GCM_3_optimized | 2 | the case study's contaminant extension of GCM_2: the same family with a mixture over it |
| OptionalStopping_Answer_1 | 1 | the chapter's exercise answers, on the model OptionalStopping already ports |

## The Laplace caveat, and where it can be measured

Rows where a CONTINUOUS random effect sits under something that is not
Gaussian in it: Exams_2, SDT_2, PsychophysicalFunction1 and 2, MPT_2 to
_5, Geurts, Extraversion. The identity is still exact there, because it
compares the JOINT at the conditional modes and both programs hold the
same joint. What is approximate is the MARGINAL frmtmb maximizes, and
the identity cannot see it.

`frm(importance = )` is the measurement, and its scope refuses by name:
one response, one grouping factor, a family with a ROWWISE density, no
nonlinear predictor. So

- admitted: Exams_2, SDT_2, PsychophysicalFunction1 and 2, MPT_2 to _5,
  Geurts;
- refused for having two responses: Extraversion, Correlation_2;
- refused for a nonlinear predictor: Retention_3;
- refused for having no rowwise density: GCM, SIMPLE, Survey,
  ChangeDetection, TwoCountryQuiz. That is the same fact that put them
  on the `loglik` slot, and for them the Laplace approximation is not
  involved at all: they have no random effects, and their sum over a
  discrete latent is exact.

### The Laplace error, measured

Exams_2, fifteen subjects, one binomial mixture with a random intercept
on the free component's rate:

```
laplace logLik    -44.09545
importance logLik -44.10878   (512 draws, mcse 0.0018)
```

So the approximation is off by about 0.013 nats, which is seven times
the Monte Carlo standard error of the number that measures it. Small,
and not zero, and the identity cannot see it at all: the identity
compares the joint at the conditional modes, and both programs hold the
same joint.

## Label switching in TwoCountryQuiz

Worth recording because it looked like a bug and is not. Every question
pattern has a complement under which every person's country flips, and
the two carry the same likelihood, so the marginal posterior over any
one person's country is EXACTLY one half. `latent_probs()` reports
that, because it is the answer.

The quantity that survives the flip is whether two people came from the
same country, and `bcm_two_country_agreement()` computes it from the
same sum: given a question pattern the people are independent, so
within a pattern the pairwise agreement is `p_i1 p_j1 + p_i2 p_j2` and
the patterns are averaged with their own weights. On an eight by eight
quiz with two blocks it recovers the blocks exactly: the smallest
within-block agreement is 1 to eight decimal places and the largest
across-block agreement is 2e-8.

The same fact appears one chapter earlier, in Malingering_1: the two
mixture components are exchangeable under maximum likelihood, so the
book's order restriction is a prior and `sort()` is the honest reading
of the estimate.

## After wt-protocol lands: `accepts_aterms`

The protocol lane adds `frmtmb_family(accepts_aterms =)`, a family-level
declaration of which addition terms a family reads. It defaults to
`NULL` and a family written before it keeps every term, so nothing in
`inst/bcm/` breaks either way. Once both are in, each family here could
declare its own:

| family | would declare |
|---|---|
| `bcm_binomial_band()` | `trials`, `vint1`, `weights` |
| `bcm_gaussian_probit()` | `vreal1` |
| `bcm_contaminant()`, `bcm_binomial_probit()` | `trials` |
| `bcm_mpt_pairs()`, `bcm_kappa()` | `trials` |
| `bcm_gcm()`, `bcm_simple()` | `trials` |

That is an improvement rather than a fix: the declaration produces a
new REFUSAL, not a new acceptance, so it would not open either of the
gates seams 2 and 3 record. `se()` stays gated on the family name and
`cens()` stays refused for discrete families; the protocol lane's own
notes say it hardens around the second of those rather than opening it.

This lane's two cens-gate assertions match the refusal MESSAGES
verbatim (`tests/testthat/test-bcm-data-analysis.R:420-437`), which is
the one place a merge could break something quietly. The protocol lane
states those messages are unchanged, and the assertions were re-run in
a fresh process against this worktree's core after the review; they
should be re-run again after the merge rather than assumed.

## Data

Small published tables are embedded in the test files. Three data sets
are not, and the tests say so where it matters:

- ChangeDetection's thousand-point series: a series of the same shape
  is generated instead, and the test asks the model to find the change
  that was put there.
- GCM_2's forty subjects: ten of them.
- ESP/Ability's hundred pairs of proportions: the spelling is
  Correlation_1's and is validated there and at OptionalStopping.

No third-party data file is added to the package. Everything a test
needs is either a published table typed into the test file or generated
in it.
