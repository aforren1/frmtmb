# Sample a model with NUTS

Runs
[`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html) on
the model's objective and returns the draws with frmtmb parameter names.
Given a
[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
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
  prior = NULL,
  init = NULL,
  init_jitter = 0.25,
  reparameterize = TRUE,
  data2 = list(),
  start = NULL,
  control = frmtmb_control(),
  na.action = stats::na.omit,
  REML = FALSE,
  .diagnostic = FALSE
)
```

## Arguments

- fit:

  A `frmtmb_fit`, or a
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.html)/formula
  to assemble and sample directly (then `data` is required).

- data:

  A data frame of model data, when `fit` is a formula.

- family:

  Family, when `fit` is a plain formula that does not carry one
  (`frm_sample(bf(y ~ x), data = dd, family = poisson())`; the `+`
  spelling `bf(y ~ x) + poisson()` works too).

- ...:

  Passed to
  [`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html)
  (`chains`, `iter`, `laplace`, `cores`, ...). `cores` parallelizes over
  chains on every platform. On Windows the chains run on socket workers,
  each of which rebuilds the tape from the serialized objective closure
  (tmbstan retapes on the worker; the closures are self-contained),
  giving draws identical to a sequential run at the same seed. The
  per-worker startup (a new R process, package load, retape) is a fixed
  several seconds, so short chains gain nothing; long chains approach a
  per-chain speedup. A core count inherited from `options(mc.cores)`,
  which rstan reads when `cores` is not given, behaves the same way.

- prior:

  Priors: a
  [`frmtmb::set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)
  specification, or a named list of prior objects (see
  [`frmtmb::prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.html))
  whose names are parameter names as in the draws (or whole components:
  `"beta"`, `"theta"`, ...), or the string `"flat"`. The objective is
  re-taped with the prior terms added; a fitted model itself is
  unchanged. On both paths the brms default priors apply to whatever
  this argument, and a MAP fit's own prior, leave alone (see Default
  priors), and `prior = "flat"` opts out of them entirely. A `brmsprior`
  object built by brms's own
  [`prior()`](https://aforren1.github.io/frmtmb/reference/prior.html) is
  translated row by row. The argument takes brms's spelling, `prior`;
  the `priors` of releases before 0.43 is gone rather than aliased, and
  because this function's `...` would otherwise swallow it, the old name
  is refused by name.

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

- reparameterize:

  Sample the qualifying random-effect blocks in their non-centered form
  (`b = L(theta) z`), which removes the funnel of the centered joint
  posterior. `TRUE` by default. Not every block qualifies, and the call
  says which did not and why; see Reparameterization, which also
  explains why nothing downstream can tell the two routes apart. The
  mode-anchored `init` still starts at the ML mode: it is mapped through
  `z0 = L(theta_hat)^-1 b_hat`. Seeded draws differ between the two
  routes, so `reparameterize = FALSE` is also the way to reproduce a run
  made before this default existed.

- data2, start, control, na.action, REML:

  As in
  [`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html);
  used only on the formula path.

- .diagnostic:

  Internal, and named with a dot because there is no reason to set it by
  hand. `TRUE` turns the default priors off and silences the disclosure,
  and leaves the fit's own objective alone: a MAP fit's penalty is taped
  into it, so no prior is resolved again and no tape is rebuilt.
  [`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md)
  sets it, because it measures the Laplace and Wald approximations
  against the density the fit maximized and a default prior would change
  what is being measured. Write `prior = "flat"` for the user-facing
  opt-out, which goes further and rebuilds a MAP fit's objective without
  its penalty.

## Value

An object of class `frmtmb_draws`: list with the `stanfit`, a draws
matrix with named columns
([`as.matrix()`](https://rdrr.io/r/base/matrix.html) method), the
originating fit, and, when any block was non-centered, a `reparam` note
saying which. The `stanfit` holds the parameters as Stan sampled them,
so on a non-centered run its random-effect columns are `z` while the
draws matrix holds `b`.

## Details

**Both routes sample a posterior under the same default priors** (see
Default priors), brms 2.23's own weakly-informative ones. What the fit
adds is a starting point, not a different density: the ML mode anchors
the chains and shortens warmup, and the model is the one the fit already
carries.

The two differ only in what else is in the stack. A fit made with
`frm(prior = )` is a MAP fit, and it carries its own prior into
sampling: the priors that apply are then this call's `prior =` first,
the fit's own next, and the defaults under both, most explicit winning
per slot it addresses. `prior = "flat"` opts out of all of it and
samples the bare likelihood, on either route.

[`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md)
is the one caller that stays unpenalized by default. It measures the
Laplace and Wald approximations against the shape of the very objective
the fit maximized, and a default prior would change the thing being
measured, so it asks for the flat density explicitly rather than
inheriting whatever this function defaults to.

## Sampling from a formula

`frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian())`
parses, assembles and tapes exactly as
[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
does, stops before the optimizer, and hands the objective to Stan. Every
pre-optimizer refusal still applies (REML, quadrature, the mixture and
[`frmtmb.latent::hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.html)
guards). There is no mode, so the default `init` is `"random"`: Stan's
own overdispersed initialization on the unconstrained scale, inside any
bounds the prior carries.

The returned object supports the whole draws surface -
[`summary()`](https://rdrr.io/r/base/summary.html),
[`frmtmb::fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.html),
[`frmtmb::VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.html),
[`frmtmb::ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.html),
[`frmtmb::hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html),
[`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
[`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
[`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
[`frmtmb::pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.html),
[`frmtmb::as_draws()`](https://aforren1.github.io/frmtmb/reference/as_draws.html) -
because those read the model frame and one draw at a time. The methods
that report a maximum-likelihood quantity refuse instead of inventing
one:
[`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md)
(it compares NUTS against a mode that does not exist here), and on the
embedded object reachable as `x$fit`,
[`summary()`](https://rdrr.io/r/base/summary.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`frmtmb::fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.html),
[`frmtmb::ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.html),
[`frmtmb::VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html).

## Reparameterization

A random-effect block has a FUNNEL in its centered joint posterior: the
width of the prior on `b` is a standard deviation that is being sampled
at the same time, so the region NUTS must explore narrows as that
standard deviation shrinks, and one step size cannot fit both ends of
it.

`reparameterize = TRUE` (the default) samples `z ~ N(0, I)` instead and
computes `b = L(theta) z` on the tape, with `L` the block's own Cholesky
factor, which is brms's construction. Each draw is mapped back through
ITS OWN `theta`, so the `b[i]` columns of the draws matrix hold the same
quantity in the same order under the same names as
`reparameterize = FALSE` gives, and every method downstream
([`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
[`frmtmb::ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.html),
[`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md),
[`frmtmb::loo()`](https://aforren1.github.io/frmtmb/reference/loo.html),
[`frmtmb::conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html),
[`frmtmb::hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html))
reads them without knowing which route produced them. Only the `stanfit`
inside the object carries `z`.

Priors are unaffected. They are declared on the OUTER parameters
(`beta`, `theta`, `betad`), and `theta` means the same thing on both
routes; nothing prior-able is reparameterized. What priors DO decide is
which blocks are eligible, below.

*Which blocks, and why not all of them.* A block is non-centered when
two things hold, and they are the same thing twice: every parameter it
has is a standard deviation or a correlation with a Cholesky factor
registered for its structure, and every one of those parameters carries
a PRIOR. Both are about not handing the sampler a direction it can run
away in. Non-centering gives NUTS the run of `theta`'s whole range, and
the centered funnel is what was keeping a chain out of the parts of that
range where the posterior is flat or improper. Removing the funnel
without closing those off first trades a slow chain for a wrong one.

In practice that means a default-prior run, on EITHER route, whose
priors cover every standard deviation and every correlation, non-centers
`(1 | g)` and any block written one term at a time,
[`diag()`](https://rdrr.io/r/base/diag.html) and `homdiag()` blocks,
[`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html) smooths, `equalto()`,
`gr(cov = )`, and the CORRELATED blocks `(x | g)`, `cs()`, `ar1()` and
`hetar1()`. A `k =` Hilbert-space `gp()` block stays centered on the
formula route even though its factor is diagonal: its LENGTHSCALES share
the block's `theta`, the default priors cover only standard deviations,
and the gate wants every parameter of a block priored. Prior the whole
block by hand (`set_prior(class = "theta")`) to non-center it. `rr()` is
already non-centered by construction, since its own coefficients are the
standard normal factors. What stays centered is a run whose variance
parameters are flat: `prior = "flat"`, and
[`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md),
which asks for that density on purpose.

The call [`message()`](https://rdrr.io/r/base/message.html)s every block
it left centered, with the reason. A model made only of those samples
centered throughout and says so rather than failing. The reasons:

- A FLAT PRIOR on the block's standard deviation. Send `sd` down with
  `z` where it is: `b = sd z` goes to zero, the likelihood settles on
  the model without that random effect, and the density stops changing:
  a flat tail with nothing to stop a chain in it. Measured on a
  six-group random-intercept model, one chain of 2000, three seeds:
  under flat priors the non-centered chain walks `theta` to -1e15 at a
  bulk-ESS of 1; under the default `student_t(3, 0, s)`, which makes
  that tail integrable, 174 to 284 against 3 to 48 centered.

- A FLAT PRIOR on the block's CORRELATION, which is the same argument:
  flat on frmtmb's unbounded correlation parameter is
  `(1 - rho^2)^-3/2`, improper, with all its mass at `|rho| = 1`. Before
  0.39 that kept every correlated block centered, because no proper
  correlation prior existed to close the tail; `lkj(eta)` now does, and
  it is a default, so `(Days | Subject)`, `cs()`, `ar1()` and `hetar1()`
  non-center like any other block.

- A Student-t latent (`gr(dist = "student")`): a scale mixture, not a
  linear factor.

- `car()`, `spde()` and `gr(prec = )`: sparse precisions whose factor is
  dense.

- The exact `gp()` and the spatial covariances (`ou`, `exp`, `gau`,
  `mat`): a full factorization per gradient evaluation.

- `toep()`: not positive definite everywhere in its parameterization, so
  `b = L z` is not a bijection there.

*What it is worth.* One chain of 2000 iterations, three seeds
(dev/benchmarks.md). Where each group's own data say little (the regime
the funnel lives in) it is decisive: 80 groups of 2 binary observations
run at a min-ESS of 236 against 5 centered, 55 times the effective draws
per second. Where the groups are informative it is a wash: the
`epilepsy` GLMM and an uncorrelated `sleepstudy` are within noise of the
centered chain either way. It is not a blanket speed-up, and the blocks
it declines to touch are the ones where it would have done harm.

`laplace = TRUE` integrates the random effects out, which removes the
funnel by itself, so that route ignores `reparameterize` entirely.

The ML fit is untouched either way. The Laplace approximation integrates
`b` out, and that integral is invariant under a linear change of the
integrated variable, so
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) has
nothing to change and no fitted object ever carries a non-centered tape.

## Default priors

Both routes default to brms 2.23's own weakly-informative priors, read
off
[`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
on matched models. Write `y*` for the response transformed by the `mu`
link and `s = max(2.5, round(mad(y*), 1))`:

|  |  |  |
|----|----|----|
| class | default | scale |
| `b` (slopes) | flat | \- |
| `Intercept` | `student_t(3, round(median(y*), 1), s)` | link |
| `sd` | `student_t(3, 0, s)` | natural sd, log-Jacobian applied |
| `cor` | `lkj(1)` | correlation matrix, Jacobian applied |
| `sigma` (intercept only) | `student_t(3, 0, s)` | natural |
| `sigma` (with a predictor) | `student_t(3, 0, 2.5)` | log |

The link is transformed only for `identity`, `log`, `inverse`, `sqrt`
and `1/mu^2` - brms's own list - with a log-scale family
([`frmtmb::lognormal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.html)
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
[`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.html)
on the returned draws reproduces the chosen priors exactly.

*What is deliberately NOT matched.* Each of these is named in the
message whenever the model has one.

- A `toep()` correlation stays flat. Its banded parameterization is not
  positive definite everywhere, so there is no correlation matrix over
  its whole parameter space for an LKJ density to be about. Every other
  correlated structure gets `lkj(1)`, brms's own default (see
  [`frmtmb::set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)
  for what the density is on frmtmb's parameters).

- ORDINAL THRESHOLDS stay flat. brms priors them `student_t(3, 0, 2.5)`
  under its `Intercept` class; frmtmb keeps them in the `thres`
  extra-parameter vector, which is not a design column and which
  [`frmtmb::set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)'s
  class vocabulary (`b`/`Intercept`/`sd`/`theta`) cannot address. An
  ordinal model still gets its `sd` defaults, and the message names the
  gap.

- The `shape`, `phi` and `nu` dispersion parameters stay flat: brms
  gives them gamma and inverse-gamma defaults, which
  [`frmtmb::set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)
  does not carry.

- NONLINEAR parameters stay flat, which is what brms does: their
  coefficients are class `b`, and the response's median and mad say
  nothing about a rate or a shape sitting inside a nonlinear body. Write
  them with `set_prior(nlpar = )`, as the brms nonlinear vignette does.

- MULTIVARIATE models get no defaults at all: the default location and
  scale are read off ONE response, and frmtmb does not read them per
  response. `set_prior(resp = )` addresses one response, so they can be
  written by hand.

*Overriding and opting out.* A
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)
specification takes over the classes it names and leaves the other
defaults in place, which is brms's partial-override rule; a named list
of prior objects takes over exactly the internal parameters it names.
`prior = "flat"` turns the defaults off entirely and samples the
likelihood, which warns when the model has variance components: their
flat-prior posteriors need not be proper, and neither the chains nor
Rhat can see that. On a MAP fit `"flat"` reaches further than the
defaults: that fit's own prior is taped into its objective, and the bare
likelihood is what `"flat"` names, so the penalty is rebuilt out and the
call says which prior it dropped.

*Three sources, most explicit first.* On a fitted model the stack is
this call's `prior =`, then the fit's own prior (`frm(prior = )` makes a
MAP fit, and that prior is part of the model the fit is), then the
defaults. Each addressed slot is settled by the most explicit source
that names it, one slot at a time: a call-level
`set_prior(class = "sd")` replaces the fit's `sd` prior and the `sd`
default while the `Intercept` default stays.
[`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.html)
on the returned draws prints what the stack came to.

## Multimodal posteriors

For
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
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
                  prior = set_prior("exponential(1)", class = "sd"))
prior_summary(ds3)
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 1, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 1, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: The largest R-hat is 1.07, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 1, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   b                  (flat), as brms leaves slopes
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
