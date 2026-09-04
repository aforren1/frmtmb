# Migrating the posterior workflow from brms

[`vignette("brms-migration", package = "frmtmb")`](https://aforren1.github.io/frmtmb/articles/brms-migration.html)
ports the model GRAMMAR: what
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.html), the
addition terms, the group-level structures and the families do here, and
how a `brm()` call becomes a
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) call.
None of that changes when you sample.
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
takes the same
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.html) object,
assembles it through the same code path, and refuses what
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) refuses.

This page is the other half: what happens to a brms POSTERIOR workflow.
The draws come from Stan’s NUTS sampler, but through tmbstan on an RTMB
objective rather than a generated Stan program, and that one difference
explains every entry below.

## The call

``` r

library(frmtmb)
library(frmtmb.sample)

# brms
brm(y ~ x + (1 | g), data = dd, family = gaussian(),
    chains = 4, iter = 2000, cores = 4)

# here
frm_sample(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           chains = 4, iter = 2000, cores = 4)
```

| brms | here |
|----|----|
| `brm(formula, data, family)` | `frm_sample(formula, data =, family =)`, or `frm_sample(fit)` on a fitted model |
| `chains`, `iter`, `warmup`, `thin`, `cores`, `seed`, `refresh` | the same names, passed straight through to [`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html) |
| `prior = prior(...)` | `prior = set_prior(...)`; a `brmsprior` that brms itself built is translated row by row |
| the default prior set | applied on both routes, read off [`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html); `prior = "flat"` opts out |
| `init =`, `init_r =` | `init =`, `init_jitter =`; from a fit the maximum-likelihood mode anchors the chains |
| `control = list(adapt_delta = )` | not through [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md), whose `control` is [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.html); see below |
| `backend =`, `threads =`, `stan_model_args =` | nothing to configure: there is no Stan program to compile |
| `sample_prior = "only"` | [`frmtmb::frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html), which draws a parameter vector per simulation and returns it with the responses |

The one argument with no route is rstan’s `control` list.
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
already spells `control` as
[`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.html),
the model-assembly options, so `adapt_delta` and `max_treedepth` cannot
travel past it.
[`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
passes everything to tmbstan and takes them:

``` r

sf <- as_tmbstan(fit, chains = 4, control = list(adapt_delta = 0.99))
```

What that costs is the return value.
[`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
gives back a `stanfit`, whose parameters carry Stan’s names and none of
the method surface below.

## Two routes, and which one to take

[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
takes a fit or a formula, and both sample the same posterior under the
same priors. Coming from brms, the formula route is the familiar one: it
is `brm()` with a different function name. The fit route is the one brms
has no counterpart to.

``` r

fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
ds <- frm_sample(fit, chains = 4)
```

Fitting first is worth the seconds it takes. The optimizer either finds
the mode or reports why it could not, `diagnose(fit)` says whether the
model is identified before any chain runs, and the mode then anchors the
chains and shortens warmup. A model that will not converge under
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) is a
model whose posterior geometry is about to be difficult, and finding
that out in a second is better than finding it out in an hour.
[`vignette("sampling")`](https://aforren1.github.io/frmtmb/frmtmb.sample/articles/sampling.md)
has the detail on both routes, the default priors and the non-centered
parameterization.

## The names change

Parameter names drop parentheses on the draws side: `Intercept`, not
`(Intercept)`. That is the vocabulary posterior, bayesplot and
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html)
already speak, so a ported
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html)
string usually needs no edit.

``` r

variables(ds)
hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)")
```

[`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
[`neff_ratio()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
and
[`nuts_params()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
are the exception. They read the `stanfit` directly and report Stan’s
own names.

## The method surface ports

The post-processing generics are brms’s own, on a `frmtmb_draws` object,
so most ported code runs unchanged:
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.html),
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.html),
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.html),
[`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
[`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
[`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
[`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md),
[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html),
[`waic()`](https://aforren1.github.io/frmtmb/reference/loo.html),
[`bayes_R2()`](https://aforren1.github.io/frmtmb/reference/bayes_R2.html),
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html),
[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html),
[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.html),
[`mcmc_plot()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
[`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.html),
[`as_draws()`](https://aforren1.github.io/frmtmb/reference/as_draws.html)
and the `posterior` and `coda` conversions.

[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) here is
the real leave-one-out, not the
[`AIC()`](https://rdrr.io/r/stats/AIC.html) substitute core offers, and
it returns the loo package’s own object.

## What refuses, and what to write instead

Every method below is defined, and every one stops with its reason and
its replacement, so a ported script fails where it went wrong rather
than at “could not find function”.

| brms call | why it stops | write instead |
|----|----|----|
| `stancode(ds)` | there is no Stan program; the R closure is the source | `ds$fit$obj$fn`, `ds$fit$frame` |
| `standata(ds)` | nothing is exported to a Stan data list | `ds$fit$frame`, [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html), `getME()` |
| `expose_functions(ds)` | a custom family is already plain R | call the `lpdf` you passed to [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.html) |
| `plot(ds)` | brms’s default panel is trace and density | `mcmc_plot(ds)`, `mcmc_plot(ds, type = "trace")` |
| `update(ds)` | draws carry no formula to revise | [`update()`](https://rdrr.io/r/stats/update.html) the fit, then [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md) it |
| `restructure(ds)` | no older-object upgrade path exists | re-run [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md) |
| `posterior_samples(ds)` | deprecated brms spelling | `as_draws(ds)`, `as.matrix(ds)` |
| `nsamples(ds)` | deprecated brms spelling | `ndraws(ds)`, `niterations(ds)` |
| `parnames(ds)` | deprecated brms spelling | `variables(ds)` |

Two groups refuse for reasons that are not spelling.

**The refit-based loo remedies.**
[`reloo()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md),
[`kfold()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md),
[`loo_moment_match()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
and
[`loo_subsample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
each need to refit the model on modified data, and that machinery is not
here. When
[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) reports
high Pareto k values, sample under the default priors rather than
`prior = "flat"` and read the number again; if the warning survives,
compare the maximum-likelihood fits by
[`AIC()`](https://rdrr.io/r/stats/AIC.html), or measure the influence
with
[`frmtmb::frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.html).

**The marginal-likelihood family.**
[`bridge_sampler()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md),
[`bayes_factor()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
and
[`post_prob()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
are integrals against the prior. They are undefined under
`prior = "flat"`, and even under the default priors bridge sampling
needs a normalized log-posterior evaluator that the RTMB tape does not
expose. This is the one brms capability with no route here at all.

## When you still want brms

- A Bayes factor or any other marginal-likelihood quantity.
- [`reloo()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
  or
  [`kfold()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
  on a model whose LOO approximation fails.
- A family, a term or an addition term frmtmb does not implement.
  [`frmtmb::frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html)
  answers that per feature, and the answer is the same for
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) and
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).
- Stan code you intend to read, extend or hand to somebody else.

A practical division: screen models with
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html), which
costs seconds, sample the survivors with
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md),
and keep brms for the questions in that list.
