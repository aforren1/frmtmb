# Diagnosing a sampled fit

[`vignette("diagnostics", package = "frmtmb")`](https://aforren1.github.io/frmtmb/articles/diagnostics.html)
asks three questions of a maximum-likelihood fit: did the optimizer
converge, does the family describe the data, and are the reported
uncertainties trustworthy.

A sampled fit reorders that list. There is no optimizer to interrogate,
so the first question becomes whether the sampler explored the posterior
it was given. The third question changes character: the posterior IS the
uncertainty statement, so the interesting comparison is against the
maximum-likelihood fit that the same objective produced.

``` r

library(frmtmb)
library(frmtmb.sample)

set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)

fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
ds <- frm_sample(fit, chains = 4)
```

## Did the sampler explore the posterior?

`summary(ds)` puts `n_eff` and `Rhat` beside each estimate, which is
where to start. Read the pattern and not only the worst entry. One
parameter with a low effective sample size and a clean `Rhat` is a slow
chain. Several at once usually means the geometry is wrong rather than
the run too short, and sampling longer will not fix it.

``` r

summary(ds)
rhat(ds)
neff_ratio(ds)
```

[`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
and
[`neff_ratio()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
read the `stanfit` inside the draws object directly, so they report
Stan’s own parameter names (`par[1]` and its relatives). Everything else
in this package relabels to the frmtmb draws-side names, which
`variables(ds)` lists. That difference is worth knowing before you try
to match two printouts by eye.

The sampler’s own record is in
[`nuts_params()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
with the divergent transitions, the tree depths and the step size per
iteration, and the log posterior is
[`log_posterior()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md).
Both feed bayesplot’s `mcmc_nuts_*()` displays.

``` r

np <- nuts_params(ds)
sum(np$Value[np$Parameter == "divergent__"])   # zero is the target

mcmc_plot(ds, type = "trace")
mcmc_plot(ds, type = "rank_overlay")
pairs(ds, variable = c("sd_g__Intercept", "sigma"))
```

Divergences concentrate in the funnel of a random-effect block, so the
[`pairs()`](https://rdrr.io/r/graphics/pairs.html) display to draw is a
group standard deviation against something it multiplies.

Three moves, in the order worth trying:

1.  **Give the variance parameters a prior.** `reparameterize = TRUE` is
    the default, but a block qualifies for the non-centered form only
    when every one of its variance parameters carries a prior, and the
    call reports each block it left centered and why. Under
    `prior = "flat"` no block qualifies, so a flat run is a centered run
    by construction. `prior = set_prior("exponential(1)", class = "sd")`
    is the usual fix, and it fixes two things at once: it closes the
    flat tail at `sd = 0` and it unlocks the reparameterization that
    removes the funnel.
2.  **Sample longer.** `iter =` and `chains =` pass straight to tmbstan.
    A low `n_eff` with `Rhat` near 1 is answered here.
3.  **Tighten the sampler.** `adapt_delta` does not travel through
    [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md):
    its `control` argument is
    [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.html),
    the assembly options, and rstan’s `control` list has no route past
    it. `as_tmbstan(fit, control = list(adapt_delta = 0.99))` hands
    everything to tmbstan, at the cost of getting a `stanfit` back
    rather than a draws object with the method surface on it.

One check comes before all of these. A tmbstan built against StanHeaders
2.39 samples a standard normal instead of the model, and every
diagnostic on this page looks perfect when it does.
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
refuses such a build before sampling; the `sampling` vignette has the
detail.

## Does the model describe the data?

[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.html)
overlays the observed response with responses drawn from the posterior
predictive distribution. This is the check brms users already run, on
the same bayesplot displays, and it is a stronger statement than the
maximum-likelihood version in core, because the draws carry parameter
uncertainty into the simulated responses.

``` r

pp_check(ds)
pp_check(ds, type = "stat", stat = "sd")
```

[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) is the
real leave-one-out here, not the
[`AIC()`](https://rdrr.io/r/stats/AIC.html) substitute core offers. It
runs Pareto-smoothed importance sampling on the
[`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md)
matrix and returns the loo package’s own object, so
[`print()`](https://rdrr.io/r/base/print.html) and
[`loo::pareto_k_table()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
work on it unchanged.

``` r

loo(ds)

# a comparison is two sampled models, each with its own elpd
ds0 <- frm_sample(bf(y ~ 1 + (1 | g)) + gaussian(), data = dd, chains = 4)
loo_compare(ds, ds0)
```

A high Pareto k is a real finding here, because the refit-based remedies
are absent:
[`reloo()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md),
[`kfold()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md),
[`loo_moment_match()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
and
[`loo_subsample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
all refuse, and say so, since each needs to refit the model on modified
data. What to do instead depends on why k is high. Under
`prior = "flat"` the elpd is unregularized and a model with many
group-level parameters will produce those warnings by construction, so
the first move is to sample under the default priors and read the number
again. If the warnings survive that, the observations loo names are
genuinely influential: compare the models by
[`AIC()`](https://rdrr.io/r/stats/AIC.html) on the maximum-likelihood
fits, or measure the influence directly with
[`frmtmb::frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.html).

Read every one of these numbers as a posterior quantity. Both of
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)’s
routes apply brms’s default priors, so an elpd is regularized the way
brms regularizes it unless the call opted out.

## Does the Laplace approximation hold?

[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
integrates random effects out with a Laplace approximation and reports
Wald standard errors. Both assumptions degrade in the same places:
variance components from few groups, binary data in small clusters. Core
routes around that with `confint(method = "profile")`,
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.html)
and `frm(quadrature = TRUE)`. This package measures it instead.

[`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md)
samples the density the fit maximized, and no other. It hands Stan
`fit$obj` as it stands, which is why it opts out of
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)’s
default priors explicitly rather than inheriting them: a prior nothing
in the fit ever saw would change the very thing being measured. For an
ordinary fit that objective is the likelihood, so the run is flat. For a
MAP fit (`frm(prior = )`) the penalty is taped into the objective, so
the run carries it, which is right, because the mode and the Wald
standard errors it compares against come from the same penalized
Hessian.

``` r

cl <- check_laplace(fit, chains = 2, iter = 1000)
cl
```

One row per outer parameter, with the maximum-likelihood estimate
(`ml`), the posterior mean (`post_mean`), the two uncertainty statements
(`wald_se`, `post_sd`), and the two comparisons worth reading:

- `z_shift` is `(post_mean - ml) / post_sd`, the distance from the mode
  to the posterior mean in posterior standard deviations. It measures
  asymmetry the quadratic approximation cannot see.
- `sd_ratio` is `post_sd / wald_se`. Above one, the Wald interval is too
  narrow.

``` r

# the parameters whose Wald interval to replace with a profile or
# bootstrap one
cl[abs(cl$z_shift) > 0.3 | cl$sd_ratio > 1.3, ]
```

Those thresholds are working values, not a test. A variance component
estimated from eight groups will fail them and should: the log standard
deviation has a skewed likelihood at that sample size, and the skew is
the finding.

Because the defaults are off, this function samples the centered
parameterization on an ordinary fit. That costs it nothing. The
non-centered form is offered only to blocks whose variance parameters
carry priors, and the flat prior that makes this comparison meaningful
is exactly the one that leaves a tail at `sd = 0` for a non-centered
chain to walk into. Give the variance parameters a prior through
`prior =` and the run does non-center, but it is then measuring the
Laplace approximation of a different posterior, which is usually not the
question.

## Which page to read

| question | where |
|----|----|
| did the optimizer converge | [`vignette("diagnostics", package = "frmtmb")`](https://aforren1.github.io/frmtmb/articles/diagnostics.html) |
| does the family describe the data (ML) | same page: [`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.html), `residuals(type = "osa")` |
| did the chains converge | here |
| does the model describe the data (posterior) | here: [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.html), [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) |
| is the Laplace approximation good enough | here: [`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md) |

[`vignette("sampling", package = "frmtmb.sample")`](https://aforren1.github.io/frmtmb/frmtmb.sample/articles/sampling.md)
is the page on getting the draws in the first place: the two routes, the
default priors, and the non-centered parameterization.
