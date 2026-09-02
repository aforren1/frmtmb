# Sample a model with NUTS

Runs
[`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html) on
the model's objective and returns the draws with frmtmb parameter names.
Given a [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)
fit it samples the fitted objective, initialized at the ML estimates,
which shortens warmup considerably. Given a formula and `data` it
assembles the same objective without optimizing anything first (see
Sampling from a formula).

## Usage

``` r
frm_sample(
  fit,
  data = NULL,
  family = NULL,
  ...,
  priors = NULL,
  lower = NULL,
  upper = NULL,
  init = NULL,
  init_jitter = 0.25,
  data2 = list(),
  start = NULL,
  control = frmtmb_control(),
  na.action = stats::na.omit,
  REML = FALSE
)
```

## Arguments

- fit:

  A `frmtmb_fit`, or a
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md)/formula to
  assemble and sample directly (then `data` is required).

- data:

  Model data, when `fit` is a formula.

- family:

  Family, when `fit` is a plain formula that does not carry one
  (`frm_sample(bf(y ~ x), data = dd, family = poisson())`; the `+`
  spelling `bf(y ~ x) + poisson()` works too).

- ...:

  Passed to
  [`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html)
  (`chains`, `iter`, `laplace`, `cores`, ...). On Windows more than one
  core falls back to sequential chains with a warning: parallel chains
  run on socket workers, which can evaluate neither the RTMB tape nor
  the objective closure (the known RTMB limitation of tmbstan,
  tmbstan#27). The fallback also covers a core count inherited from
  `options(mc.cores)`, which is what rstan reads when `cores` is not
  given. Fork clusters on unix can, so `cores` works there.

- priors:

  Priors: a
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  specification, or a named list of prior objects (see
  [`prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md))
  whose names are parameter names as in the draws (or whole components:
  `"beta"`, `"theta"`, ...), or the string `"flat"`. The objective is
  re-taped with the prior terms added; a fitted model itself is
  unchanged. On the fit path a parameter without a prior keeps a flat
  improper one. On the formula path the brms default priors apply to
  whatever the specification leaves alone (see Default priors), and
  `priors = "flat"` opts out of them entirely.

- lower, upper:

  Optional named numeric vectors of hard bounds on outer parameters
  (brms `lb`/`ub`), applied on the internal scale through Stan's
  constrained transforms. Chain starting values are clamped strictly
  inside the bounds; a bound that excludes the ML mode itself warns,
  because the chains then no longer start there.

- init:

  Initialization. On a fit the default (`"last.par.best"`) starts chain
  1 exactly at the ML mode and every further chain at the mode plus a
  normal perturbation of sd `init_jitter` on the unconstrained scale.
  The mode anchor keeps warmup short; the jitter keeps the chains
  overdispersed enough for Rhat to retain power against multimodality
  (the standard objection to identical mode starts). `"random"` requests
  Stan's own overdispersed initialization, and is the default from a
  formula, where there is no mode to start at.

- init_jitter:

  Per-chain perturbation sd for the default init; `0` starts every chain
  exactly at the mode. Draws from the R session's RNG, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) makes the inits
  reproducible.

- data2, start, control, na.action, REML:

  As in [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md);
  used only on the formula path.

## Value

An object of class `frmtmb_draws`: list with the `stanfit`, a draws
matrix with named columns
([`as.matrix()`](https://rdrr.io/r/base/matrix.html) method), and the
originating fit.

## Details

**The two routes answer different questions.** `frm_sample(fit)` is a
DIAGNOSTIC: it explores the LIKELIHOOD, with flat improper priors on the
outer parameters, so that
[`check_laplace()`](https://aforren1.github.io/frmtmb/reference/check_laplace.md)
can compare the Laplace and Wald approximations against the shape of the
very same objective the fit maximized. `frm_sample(formula, data)` is a
SAMPLING tool: it samples a POSTERIOR, and defaults to brms's own
weakly-informative priors (see Default priors), because from a formula
there is no estimate for the run to be a diagnostic for.

On the fit path, then, a parameter without a prior keeps a flat improper
one, and the posterior of a variance component with few groups can be
improper. That is the price of measuring the likelihood, and it is why
the fit path is a diagnostic.

## Sampling from a formula

`frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian())`
parses, assembles and tapes exactly as
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) does,
stops before the optimizer, and hands the objective to Stan. Every
pre-optimizer refusal still applies (REML, quadrature, the mixture and
[`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) guards).
There is no mode, so the default `init` is `"random"`: Stan's own
overdispersed initialization on the unconstrained scale, inside any
`lower`/`upper` bounds.

The returned object supports the whole draws surface -
[`summary()`](https://rdrr.io/r/base/summary.html),
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md),
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md),
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md),
[`posterior_epred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md),
[`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md),
[`posterior_linpred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md),
[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md),
[`as_draws()`](https://aforren1.github.io/frmtmb/reference/as_draws.md) -
because those read the model frame and one draw at a time. The methods
that report a maximum-likelihood quantity refuse instead of inventing
one:
[`check_laplace()`](https://aforren1.github.io/frmtmb/reference/check_laplace.md)
(it compares NUTS against a mode that does not exist here), and on the
embedded object reachable as `x$fit`,
[`summary()`](https://rdrr.io/r/base/summary.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md),
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html).

## Default priors

The formula interface defaults to brms 2.23's own weakly-informative
priors, read off
[`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
on matched models. Write `y*` for the response transformed by the `mu`
link and `s = max(2.5, round(mad(y*), 1))`:

|  |  |  |
|----|----|----|
| class | default | scale |
| `b` (slopes) | flat | \- |
| `Intercept` | `student_t(3, round(median(y*), 1), s)` | link |
| `sd` | `student_t(3, 0, s)` | natural sd, log-Jacobian applied |
| `sigma` (intercept only) | `student_t(3, 0, s)` | natural |
| `sigma` (with a predictor) | `student_t(3, 0, 2.5)` | log |

The link is transformed only for `identity`, `log`, `inverse`, `sqrt`
and `1/mu^2` - brms's own list - with a log-scale family
([`lognormal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
and its relatives, whose `mu` link is spelled `identity` but whose
response is logged inside the density) counted as `log`. Under any other
link (`logit`, `cloglog`, ...), for a response that was a factor, and
for a family whose parameters run over categories (ordinal, categorical,
multinomial), the location is 0 and the scale 2.5 - which is what brms
reports for a bernoulli model. Under a `log`, `inverse` or `1/mu^2` link
the ZEROS of the response are replaced by `0.1` before the transform,
element by element, exactly as brms does; the non-zero values are left
alone, so a count response with median 1 keeps `log(1) = 0` as its
location.

The call [`message()`](https://rdrr.io/r/base/message.html)s the priors
it chose, one compact line per class, AND one line per slot it
deliberately left flat, so a model that gets few defaults - or none -
says so rather than looking flat by accident. Wrap it in
[`suppressMessages()`](https://rdrr.io/r/base/message.html) to silence
that.
[`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.md)
on the returned draws reproduces the chosen priors exactly.

*What is deliberately NOT matched.* Each of these is named in the
message whenever the model has one.

- Random-effect CORRELATIONS stay flat. brms puts `lkj(1)` on them,
  which is uniform over correlation matrices; frmtmb parameterizes a
  covariance block by an unconstrained Cholesky `theta` segment whose
  flat prior is NOT uniform on correlations, and no LKJ density is
  implemented on that parameterization, so claiming `lkj(1)` here would
  be false. Set one by hand with `set_prior(class = "theta")` if you
  need it.

- ORDINAL THRESHOLDS stay flat. brms priors them `student_t(3, 0, 2.5)`
  under its `Intercept` class; frmtmb keeps them in the `thres`
  extra-parameter vector, which is not a design column and which
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)'s
  class vocabulary (`b`/`Intercept`/`sd`/`theta`) cannot address. An
  ordinal model still gets its `sd` defaults, and the message names the
  gap.

- The `shape`, `phi` and `nu` dispersion parameters stay flat: brms
  gives them gamma and inverse-gamma defaults, which
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  does not carry.

- MULTIVARIATE models get no defaults at all, because
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  cannot address one response of several.

*Overriding and opting out.* A
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
specification takes over the classes it names and leaves the other
defaults in place, which is brms's partial-override rule; a named list
of prior objects takes over exactly the internal parameters it names.
`priors = "flat"` turns the defaults off entirely and samples the
likelihood, which warns when the model has variance components: their
flat-prior posteriors need not be proper, and neither the chains nor
Rhat can see that.

## Multimodal posteriors

For
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
fits the posterior is multimodal by construction (label switching at
minimum). Mode-centered inits, jittered or not, leave every chain in one
symmetry branch, so Rhat cannot flag the others; use `init = "random"`
there and inspect chains individually.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)
summary(ds)
fixef(ds)
hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)")

# the same model sampled straight from the formula, with no ML fit.
# It reports the brms default priors it chose, and prior_summary()
# gives them back.
ds2 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
                  chains = 1, iter = 500, refresh = 0)
prior_summary(ds2)
fixef(ds2)

# a set_prior() specification takes over the classes it names and
# leaves the rest of the defaults alone
ds3 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
                  chains = 1, iter = 500, refresh = 0,
                  priors = set_prior("exponential(1)", class = "sd"))
prior_summary(ds3)
}
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> frm_sample(): default priors (brms 2.23 defaults; priors = "flat" opts out)
#>   Intercept          student_t(3, 1, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.14, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> frm_sample(): default priors (brms 2.23 defaults; priors = "flat" opts out)
#>   Intercept          student_t(3, 1, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: There were 2 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> student_t(3, 1, 2.5) class=Intercept
#> student_t(3, 0, 2.5) class=Intercept dpar=sigma scale=natural
#> exponential(1) class=sd
# }
```
