# Sampling a frmtmb model

frmtmb fits by maximum likelihood. This package hands the same objective
to Stan’s NUTS sampler through tmbstan and gives back the posterior,
under brms’s default priors and with frmtmb’s own parameter names.

Nothing about the fit changes. The model is the one
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) already
assembled; what is added is a posterior, and a starting point.

``` r

library(frmtmb)
library(frmtmb.sample)

set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)

fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
ds <- frm_sample(fit, chains = 4)
summary(ds)
```

## Two routes, one posterior

[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
takes either a fit or a formula.

``` r

# from a fit: the ML mode anchors the chains and shortens warmup
ds1 <- frm_sample(fit, chains = 4)

# from a formula: assembled and taped, never optimized, so Stan's own
# overdispersed initialization is the default
ds2 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
                  chains = 4)
```

Both sample the same posterior under the same default priors. What the
fit adds is a starting point, not a different density.

## The default priors

Both routes apply brms 2.23’s weakly-informative defaults, read off
[`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
on matched models. The call reports what it chose, one line per class,
and one line per slot it deliberately left flat, so a model that gets
few defaults says so rather than looking flat by accident.

[`frmtmb::get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.html)
reports them too, once this package is loaded:

``` r

# with frmtmb alone the default column reads "(flat)", which is what
# frm() does. With frmtmb.sample loaded it reads what frm_sample()
# would apply, because this package registers its defaults with
# get_prior()
get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
```

`prior = "flat"` opts out and samples the bare likelihood. That is a
fine diagnostic and fragile inference: under a flat prior on a log
standard deviation the posterior need not be proper, and neither the
chains nor Rhat can see it, so the opt-out warns when the model has
variance components.

A
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)
specification takes over the classes it names and leaves the other
defaults alone, which is brms’s partial-override rule.

``` r

ds <- frm_sample(fit, chains = 4,
                 prior = set_prior("exponential(1)", class = "sd"))
prior_summary(ds)
```

## Non-centering

A random-effect block has a funnel in its centered joint posterior: the
width of the prior on `b` is a standard deviation being sampled at the
same time, so the region NUTS must explore narrows as that standard
deviation shrinks and one step size cannot fit both ends.

`reparameterize = TRUE`, the default, samples `z ~ N(0, I)` and computes
`b = L(theta) z` on the tape. Each draw is mapped back through its own
`theta`, so every method downstream reads the same `b[i]` columns and
cannot tell the two routes apart.

Not every block qualifies, and the call says which did not and why. A
block is non-centered only when every parameter it has is a standard
deviation or a correlation with a registered Cholesky factor, and every
one of those parameters carries a prior. Both conditions are about not
handing the sampler a direction it can run away in: removing the funnel
without closing off the flat tail at `sd = 0` first trades a slow chain
for a wrong one.

See
[`?frm_sample`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
for which structures qualify, what it is worth, and the measurements
behind that.

## Once the draws exist

Whether the chains explored the posterior, whether the model describes
the data, and whether the Laplace approximation of the
maximum-likelihood fit held in the first place are three separate
questions, and
[`vignette("posterior-diagnostics")`](https://aforren1.github.io/frmtmb/frmtmb.sample/articles/posterior-diagnostics.md)
works through them with the tools this package provides for each:
[`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
[`neff_ratio()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
and
[`nuts_params()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
for the sampler,
[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.html)
and [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) for
the model, and
[`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md)
for the approximation.

## The method surface

The draws object supports the post-processing surface a brmsfit does, on
the same generics:

``` r

fixef(ds); ranef(ds); VarCorr(ds)
posterior_epred(ds); posterior_predict(ds); posterior_linpred(ds)
hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)")
conditional_effects(ds)
pp_check(ds)
log_lik(ds); loo(ds); waic(ds); bayes_R2(ds)
as_draws(ds); as.array(ds); as.mcmc(ds)
mcmc_plot(ds, type = "trace")
```

Parameter names drop parentheses on the draws side – `Intercept`, not
`(Intercept)` – because that is the vocabulary posterior, bayesplot and
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html)
already speak. `variables(ds)` lists them.

## Known failure: a tmbstan built against StanHeaders 2.39

stanc 2.39.0 emits two `log_prob_impl` overloads where 2.32 emitted one,
and tmbstan’s install-time code generator patches only the first. HMC
then reads its value and gradient from the unpatched overload, so every
chain samples a standard normal instead of the model, with the
objective, priors and data all silently absent.

[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
and
[`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
check the installed model source once per session and refuse such a
build before sampling, rather than returning plausible garbage after it.
If you hit that refusal, install a binary tmbstan build or reinstall
tmbstan against StanHeaders 2.32.10, and distrust any draws the affected
installation already produced.
