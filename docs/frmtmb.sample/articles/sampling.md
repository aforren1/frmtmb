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
reports them too, when you ask for this route. It has two routes,
because frmtmb has two: `route = "fit"` is the default and reports what
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) applies,
which is flat in every slot. `route = "sample"` reports what
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
applies, and this package is what lets it answer.

``` r

# what frm() applies. Loading this package does not change it
get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# what frm_sample() applies, which is the brms reading of get_prior().
# Without this package loaded, the call is refused rather than
# answered "(flat)"
get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
          route = "sample")
```

The printed table names its route on the first line, so a table copied
out of a session still says which set of defaults it describes.

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

## Conditional effects draw the frame the fit draws

`conditional_effects(ds)` returns what `conditional_effects(fit)`
returns: the same columns in the same order, the same grid, the same
keys, and the same attributes brms’s
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) reads. What
differs is where the band comes from. The curves are posterior
expected-response draws, so `estimate__` is their mean and `lower__` and
`upper__` are their percentiles, and `attr(df, "band")` reads
`"posterior"`.

The default display is the EXPECTED RESPONSE, not the `mu` predictor. On
a zero-inflated or hurdle family, a mixture, or a
[`trunc()`](https://rdrr.io/r/base/Round.html) response those are
different curves, and the one you want is almost always the mean. Name a
`dpar =` to get a predictor instead.

Three arguments mean here what they mean on a fit. `int_conditions =`
labels a moderator, `allow_new_levels =` is accepted in both spellings,
and `re_formula = NULL` conditions on a NEW group: the grouping column
of the returned frame reports `NA`, and each posterior draw carries that
group’s effects drawn from its own covariance parameters. The band is
wider than the population band by what that group’s variance adds, which
is the whole point of asking. `seed =` makes the draw reproducible.

``` r

conditional_effects(ds, effects = "x")                 # population
conditional_effects(ds, effects = "x", re_formula = NULL, seed = 1)
conditional_effects(ds, effects = "x:z",
                    int_conditions = list(z = c(low = -1, high = 1)))
conditional_effects(ds, categorical = FALSE)  # ordinal: E[category]
```

There is no `method =` and no `band =`. The curves ARE the posterior, so
there is no wald/profile/boot choice to make. For a predictive band,
quantile
[`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
over a grid you build yourself.

## Evidence ratios and Bayes factors

[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html)
reports `evid_ratio` and `post_prob` beside the interval. They are
brms’s `Evid.Ratio` and `Post.Prob`, and they answer two different
questions depending on how the hypothesis is written.

A DIRECTIONAL hypothesis is the simple one. Its evidence ratio is the
posterior odds of the claim, counted from the draws. No prior enters it,
so it is always available.

``` r

h <- hypothesis(ds, "x > 0")
h$evid_ratio   # P(x > 0) / P(x < 0)
h$post_prob    # P(x > 0)
```

A POINT hypothesis is the Savage-Dickey density ratio: the posterior
density of the tested quantity at the point, divided by its prior
density there. For a point null nested inside the model, that ratio IS
the Bayes factor. It needs a proper prior on the parameter, and
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
leaves class `"b"` flat by default because brms does, so write one:

``` r

ds <- frm_sample(fit, prior = set_prior("normal(0, 1)", class = "b"))
hypothesis(ds, "x = 0")$evid_ratio
```

Without a proper prior the ratio is `NA` and a warning names the
parameter and says what to write. The other rows of the same call still
report.

The numerator is a kernel density of the drawn quantity, with
[`stats::density()`](https://rdrr.io/r/stats/density.html)‘s default
`nrd0` bandwidth over 4096 grid points, read at the point with a spline.
That is brms’s construction, so the two packages’ numbers are
comparable. The denominator is not brms’s: brms estimates the prior
density from `sample_prior = "yes"` draws, and this package evaluates it
from the specification the model was sampled under. Half the Monte Carlo
error goes away with it.

A density ratio deserves an error bar. The between-chain Monte Carlo
error of each ratio rides on the returned object:

``` r

h <- hypothesis(ds, "x = 0")
attr(h, "evid_ratio_mcse")
```

### What a point ratio refuses

Each refusal names the parameter and its own reason, and leaves that
row’s `evid_ratio` at `NA`:

- the parameter has no proper prior. Write one and resample.
- the hypothesis names a variance component or a correlation. Their
  priors sit on the internal parameter, not on the reported quantity,
  and a point null at zero on a standard deviation is on the boundary,
  where a density ratio answers nothing.
- the hypothesis is not an affine function of the coefficients it names.
  The prior of a nonlinear function of them has no closed form here.
- the hypothesis names the intercept. Its prior is about the intercept
  at the predictor means, which is not the coefficient being tested.

A contrast of several coefficients works when each carries a normal
prior, because the prior of the contrast is then normal too.

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
