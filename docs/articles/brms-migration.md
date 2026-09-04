# frmtmb for brms users

frmtmb reimplements a documented subset of the brms grammar with
identical spelling, so brms model code ports mechanically: change
`brm()` to
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md), keep or
drop the priors, and get a fit in milliseconds-to-seconds instead of an
MCMC run. Estimation is ML (or REML) with the Laplace approximation for
latent effects; inference is frequentist (Wald and profile intervals,
likelihood-ratio tests, AIC).

## What ports directly

| brms | frmtmb | notes |
|----|----|----|
| `bf(y ~ x + (x \| g))` | same | plus `us()`, [`diag()`](https://rdrr.io/r/base/diag.html), `homdiag()`, `cs()`, `ar1()` wrappers |
| `bf(y ~ x, sigma ~ z + (1 \| g))` | same | any dpar, full grammar |
| `bf(y ~ x, sigma = 1)` | same | fixed via TMB `map` |
| `s(x)`, `t2(x, z)` | same | mgcv smooths in any dpar; `te()`/`ti()` unsupported |
| `(1 \| p \| g)` | same | cross-formula RE correlation |
| `mvbind(y1, y2) ~ x`, [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md), [`set_rescor()`](https://aforren1.github.io/frmtmb/reference/mvbf.md) | same | per-response families; `rescor` gaussian-only |
| `y \| trials(n)`, `weights(w)`, `cens(c)`, `trunc(lb=, ub=)` | same | `cens`/`trunc` need a CDF-carrying family |
| `nl = TRUE` | same | a located prior places the start, else provide `start` ([`par_template()`](https://aforren1.github.io/frmtmb/reference/par_template.md) names it); se.fit on the nonlinear mu not yet |
| [`lf()`](https://aforren1.github.io/frmtmb/reference/lf.md), [`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) | same | [`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) on any dpar, bodies chain to any depth (below) |
| [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md) | same idea | the lpdf is plain R over RTMB advectors, not Stan code |
| [`cumulative()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md), [`multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md) | same | multinomial takes `K` explicitly |
| `car(M, gr = g, type =)` | same | all four types |
| `(1 \| mm(g1, g2))`, `mmc(x1, x2)` | same | multi-membership, `weights =` and `scale =` included (below) |
| `ar(t, g, cov = TRUE)`, `ma()`, `arma()`, `cosy()`, `unstr()` | same | gaussian/student only; `cov = TRUE` required (below) |
| `data2 = list(W = W)` | same | also resolves compound expressions (below) |

## What changes

- `brm(...)` becomes `frm(...)`; there are no chains or warmup, and a
  plain call is maximum likelihood rather than a posterior.
  `REML = TRUE` replaces priors as the small-sample correction for
  variance components.

- A `prior = c(prior(...), prior(...))` argument is kept as it stands.
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) spells
  the argument `prior`, as brms does, and
  [`prior()`](https://aforren1.github.io/frmtmb/reference/prior.md),
  [`prior_()`](https://aforren1.github.io/frmtmb/reference/prior.md) and
  [`prior_string()`](https://aforren1.github.io/frmtmb/reference/prior.md)
  build the same specification
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  does; a prior object that brms itself built is translated row by row,
  so the argument works whichever package’s
  [`prior()`](https://aforren1.github.io/frmtmb/reference/prior.md) was
  in scope.

  What the fit MEANS changes. `frm(prior = )` is MAP: the density is a
  PENALTY on the likelihood and the answer is one mode, not a posterior,
  so the reported log likelihood, AIC and
  [`anova()`](https://rdrr.io/r/stats/anova.html) are penalized
  quantities. The estimates land close where the priors are weak: a
  kidney model with brms’s own priors gives `sd(patient)` 0.38 against
  brms’s posterior mean of 0.40. But “close” is a measurement, not a
  guarantee, and the two numbers are not the same quantity. Read a prior
  here as regularization.

  A prior class frmtmb cannot mean the same way brms does (`sigma`,
  brms’s mixture `theta`) is refused by name rather than translated into
  something else.

- A prior with a location DOES place a nonlinear start. brms uses its
  priors to place its sampler, and
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) reads
  them the same way for nonlinear parameters: a `normal()`,
  `student_t()` or `cauchy()` prior on a nonlinear coefficient starts it
  at the prior location, and a message names what was placed. The
  insurance-loss growth curve of the brms nonlinear vignette fits from
  its brms priors alone.

  Without such a prior a nonlinear model still needs `start`, because
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)
  optimizes rather than samples and evaluates the objective AT the
  starting values, where a nonlinear body is rarely defined.
  `par_template(formula, data)` lists the parameter names before there
  is a fit; edit the result and hand it back as `start =`. brms has no
  counterpart to it, and needs none: `brm(init =)` takes Stan program
  names and brms initializes at random in a bounded range instead.

- Posterior summaries become ML estimates:
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md)
  returns point estimates,
  [`confint()`](https://rdrr.io/r/stats/confint.html) gives Wald or
  profile intervals, [`anova()`](https://rdrr.io/r/stats/anova.html)
  gives likelihood-ratio tests, and
  [`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
  replace [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md)
  (or install **frmtmb.sample**, whose
  [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md) is the
  real one, on sampled draws).

- `posterior_predict()` becomes
  [`simulate()`](https://rdrr.io/r/stats/simulate.html);
  `posterior_epred()` over `newdata` becomes
  `predict(newdata, type = "response")`. On a `cens()` response both
  draw the LATENT uncensored value, as brms does;
  `simulate(censored = TRUE)` applies the censoring mechanism instead.

- `sample_prior = "only"` plus `posterior_predict()` becomes
  `frm_simulate(formula, data, prior = set_prior(...))`, which draws a
  parameter vector per simulation and returns it alongside the
  responses.

- `mo()` monotonic effects, `mi()` one-step imputation of continuous
  predictors, `mi(sdx)` measurement error (the `me()` replacement, as in
  current brms), `cs()` category-specific ordinal effects, `gp(x)` /
  `gp(x, k =)` Gaussian processes, and
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  families (including group-level latent classes via `groups = ~g`) all
  work. `mo()`/`mi()` support two-way `:`/`*` interactions with numeric
  terms; discrete predictors cannot be imputed in-model, the same
  restriction Stan has. Multiply imputed data can instead go through
  [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  (Rubin’s rules), the `brm_multiple()` analog. Being frequentist, that
  fit also needs pooled tests, which brms gets for free from the
  posterior:
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  pools any function of the parameters, and
  [`anova()`](https://rdrr.io/r/stats/anova.html) on the result compares
  two nested
  [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  fits with the D1, D2 and D3 rules of `mice` (D3 by default). Mixture
  fits are ML: expect multimodality, compare starts
  ([`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)).

- [`binomial()`](https://rdrr.io/r/stats/family.html) without `trials()`
  is accepted and means Bernoulli, the
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) convention; brms
  rejects it and asks for `trials()` or
  [`bernoulli()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md).
  The divergence is permissive, so brms code ports unchanged, but frmtmb
  code written this way does not port back.

- [`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`, `cosy()`
  and `unstr()` are supported in their covariance form only; see the
  next section.

- glmer’s proportion-response idiom (`weights = size`) becomes
  `y | trials(size)` with either proportions or counts;
  [`weights()`](https://rdrr.io/r/stats/weights.html) in the formula
  stays a frequency weight, as in brms.

- `se()` works as in brms (meta-analysis), including
  `se(x, sigma = TRUE)` and the phylogenetic version with
  `gr(g, cov = A)`.

- `mm()` ports with its `weights =` and `scale =` arguments and with
  `mmc()`; its other arguments have different spellings here (below).

## Nonlinear formulas: `nl = TRUE` and `nlf()`

`nl = TRUE` says the RESPONSE formula is a nonlinear body, so only `mu`
can be nonlinear.
[`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) names the
parameter instead, so any of them can be, and the two spellings meet in
the same place:

``` r

# the same model, both ways
frm(bf(y ~ exp(b * x), b ~ 1, nl = TRUE), data = d)
frm(bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1), data = d)

# and the model nl = TRUE cannot spell: a nonlinear sigma, linear mu
frm(bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1), data = d)
```

A parameter’s link is applied to its body’s value, as in brms: the
`sigma` above is `exp(a + b * z)`. Bodies chain to any depth
(`nlf(a ~ cc * x) + nlf(cc ~ exp(b))`), are computed in dependency
order, and a cycle is refused by name.

Two deliberate differences from brms, both supersets:

- brms requires `nl = TRUE` on the
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) whenever
  the response formula names a nonlinear parameter. Here
  [`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) has
  already declared the name, so `bf(y ~ a) + nlf(a ~ ...)` is enough;
  the flag still works and means the same thing.
- A body may name another distributional parameter of the same response,
  which reads that parameter’s per-row value. brms reads such a name as
  a data column and refuses the model when no column has it. A column
  still wins here, so a ported body keeps its meaning. This is the
  spelling for a variance function of the model’s own mean, which nlme
  writes `varPower(form = ~ fitted(.))`:

``` r

# sd = exp(ls) * |mu|^th, one joint maximum likelihood
frm(bf(y ~ x) + nlf(sigma ~ ls + th * log(abs(mu))) +
      lf(ls ~ 1, th ~ 1), data = d)
```

## Multi-membership

`(1 | mm(g1, g2))` means what it means in brms: the observation belongs
to several levels of one grouping factor, and its design row averages
those levels’ effects. The pooled level set, the `1/J` default weights
and the `scale = TRUE` normalization of a supplied weight matrix are
brms’s, checked against `make_standata()`’s `J_*`, `W_*` and `Z_*`
arrays, so a ported model builds the same design matrix.

`mmc(x1, x2)` also ports: it is one random-slope coefficient of the term
whose covariate value is member specific, member `k` taking argument
`k`.

``` r

# brms and frmtmb, identical spelling
frm(bf(y ~ x + (1 | mm(school1, school2))) + gaussian(), data = d)
frm(bf(y ~ x + (1 | mm(school1, school2,
                       weights = cbind(share1, share2)))) + gaussian(),
    data = d)
frm(bf(y ~ (mmc(hours1, hours2) | mm(school1, school2))) + gaussian(),
    data = d)
```

brms’s remaining `mm()` arguments are spelled the way any other grouping
term spells them here, and the parser names the replacement when it
refuses one:

| brms | frmtmb |
|----|----|
| `mm(g1, g2, cor = FALSE)` | `diag(x \| mm(g1, g2))`, or `(x \|\| mm(g1, g2))` |
| `mm(g1, g2, id = "q")` | the `\|ID\|` key, `(x \| q \| g)` - not yet over `mm()` |
| `mm(g1, g2, cov = A)` | `gr(g, cov = A)` - not yet over `mm()` |
| `mm(g1, g2, by = )`, `pw =`, `dist =` | no equivalent yet |

`?frmtmb-multimembership` is the full page.

## Within-group residual correlation

[`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`, `cosy()`
and `unstr()` are written where brms writes them, in the model formula,
and mean what brms means under `cov = TRUE`: the residuals of one group
become a single correlated draw, `y_g ~ N(mu_g, D R D)` with `D` the
diagonal of that group’s `sigma` values. `?frmtmb-autocor` has the full
page; four things matter when porting.

**`cov = TRUE` is required for
[`ar()`](https://rdrr.io/r/stats/ar.html), `ma()` and `arma()`.** brms
defaults to `cov = FALSE`, which is a different likelihood - a residual
regression that conditions on each group’s first rows - and that form is
not implemented. The call is refused rather than reinterpreted, so a
brms model written with the default has to be changed deliberately.
`cosy()` and `unstr()` have no `cov` argument in brms either and port
unchanged.

``` r

# brms
brm(y ~ x + ar(week, subj, cov = TRUE), data = d)
# frmtmb
frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d)
```

**`sigma` is the marginal residual SD, not the innovation SD.** brms’s
`cholesky_cor_ar1()` divides by `1 - ar^2`, so its `sigma` under
[`ar()`](https://rdrr.io/r/stats/ar.html)/`ma()`/`arma()` is the
innovation scale, while under `cosy()` and `unstr()` it is the marginal
one. Here every structure uses the marginal scale, so
[`sigma()`](https://rdrr.io/r/stats/sigma.html), pearson residuals and a
`sigma ~ x` model mean one thing throughout. The correlation parameters
agree exactly; convert the scales with
`sigma_marginal = sigma_innovation / sqrt(1 - phi^2)` for AR(1).

**Only [`gaussian()`](https://rdrr.io/r/stats/family.html) and
[`student()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md).**
These are the two families brms gives a real residual covariance. brms
accepts the same spelling for poisson, binomial and the rest, but fits a
latent Gaussian AR process added to the linear predictor instead - which
is a random effect, and is spelled as one here:

``` r

# brms, poisson: a latent AR process on the linear predictor
brm(cnt ~ x + ar(week, subj), data = d, family = poisson())
# frmtmb: the same idea as a random effect over the time factor
frm(bf(cnt ~ x + ar1(factor(week) + 0 | subj)) + poisson(), data = d)
```

**Higher orders are allowed, and gaps are honored.** brms limits
`cov = TRUE` to order one; `ar(p = 3)` and `arma(p = 2, q = 1)` work
here. And the lag between two rows is the distance between their
positions in the global set of time levels, so a subject missing week 3
gets `cor(week2, week4) = rho^2`; brms indexes
[`ar()`](https://rdrr.io/r/stats/ar.html)/`ma()`/`arma()`/ `cosy()` by
position within the group and treats the missing row as no gap. On
complete balanced groups the two agree.

The parameters appear in
[`summary()`](https://rdrr.io/r/base/summary.html) under “Within-group
residual correlation” and in
[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
under brms’s names (`ar[1]`, `cosy`, `cortime__1__2`);
[`autocor_matrix()`](https://aforren1.github.io/frmtmb/reference/autocor_matrix.md)
returns the fitted correlation matrix.
[`weights()`](https://rdrr.io/r/stats/weights.html), `cens()`,
[`trunc()`](https://rdrr.io/r/base/Round.html), `se()`, `mi()`,
`rescor = TRUE`, mixtures and `quadrature = TRUE` are refused, because
the likelihood no longer factorizes over rows - brms refuses the same
core set. Random effects alongside the correlated residual are the point
of the feature and are supported.

## Variance functions (nlme `weights = varFunc`)

There are no `varFunc` objects. The `sigma` formula is the variance
function: `sigma` uses a log link, so most of the nlme vocabulary is a
direct translation, and the translation was verified against
[`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) under ML
(log-likelihood agreement 1e-10 or better).

| nlme                        | frmtmb                               |
|-----------------------------|--------------------------------------|
| `varIdent(form = ~ 1 \| g)` | `sigma ~ 0 + g`                      |
| `varExp(form = ~ v)`        | `sigma ~ v`                          |
| `varPower(form = ~ v)`      | `sigma ~ log(abs(v))`                |
| `varFixed(~ v)`             | `sigma ~ offset(0.5 * log(v))`       |
| `varComb(A, B)`             | the two terms in one `sigma` formula |

The coefficient of `log(abs(v))` is nlme’s power `theta`. The
translation is also more general than a `varFunc`: `sigma` takes the
full predictor grammar, so random effects and smooths in the variance
are ordinary spellings.

`varPower(form = ~ fitted(.))`, variance as a function of the model’s
own mean, is spelled with
[`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) (see the
nonlinear-formulas section):
`nlf(sigma ~ ls + th * log(abs(mu))) + lf(ls ~ 1, th ~ 1)` is
`sd = exp(ls) * |mu|^th`, estimated by one joint maximum likelihood.
Expect agreement with `gls` near 1e-2 rather than 1e-10: `gls`
alternates the mean fit against iteratively updated fitted values and
stops without a gradient check, while this is one joint maximization of
the same likelihood. For a positive response the coupling often belongs
in the family instead: [`Gamma()`](https://rdrr.io/r/stats/family.html)
has sd proportional to the mean, and
[`tweedie()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
estimates the power in `Var = phi * mu^p` directly.

## Structural objects: `data2`

Matrices that are not columns of `data` go in `data2`, the brms
spelling: the adjacency matrix of `car()`, the mesh triple of `spde()`,
and the matrices of `gr(prec = )`, `gr(cov = )` and `equalto()`.

``` r

frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
    data = d, data2 = list(W = W))
frm(bf(phen ~ 1 + (1 | gr(species, cov = A))) + gaussian(),
    data = d, data2 = list(A = A))
```

Two differences from brms. First, brms accepts a bare name in `data2`
and nothing else, while frmtmb evaluates the whole expression with
`data2` in front of the data mask, so `gr(g, cov = solve(Q))` finds `Q`
in `data2`. Second, `data2` is optional: a matrix visible from the
formula environment is still found, which is what frmtmb did before
`data2` existed.

Use `data2` anyway. Its objects are stored on the fit, so the fit is
self-contained: [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)
writes the matrices with it, and
[`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md),
[`update()`](https://rdrr.io/r/stats/update.html),
[`drop1()`](https://rdrr.io/r/stats/add1.html),
[`influence()`](https://rdrr.io/r/stats/lm.influence.html) and
[`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
re-assemble in a fresh session where the calling environment that built
the matrix is gone. Without `data2` those refits depend on the name
still resolving where the model was written.

## Method conventions

frmtmb keeps the brms spelling for functions that originate in brms and
the stats/lme4 spelling for standard generics:

- brms-origin, same names and argument spellings:
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  (includes what `conditional_smooths()` covers),
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  (a delta-method Wald test here),
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
  (simulate-based, through bayesplot’s `ppc_*` functions),
  [`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.md),
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
  [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md),
  [`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md),
  [`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.md),
  [`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
  (the usable parameter names, e.g. `sd_Subject__Days`).
- stats-origin generics follow stats/lme4/glmmTMB, not brms:
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns a vector
  (with `se.fit = TRUE`, a list), not a draws matrix with
  `Estimate`/`Q2.5` columns;
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) return
  vectors; [`confint()`](https://rdrr.io/r/stats/confint.html) takes
  `level =`, not `probs =`; `re.form` replaces `re_formula`. The
  exception is an ordinal family, where
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) follows brms
  and returns the `n` by `K` matrix of category probabilities named by
  the response levels, exactly as `predict(type = "response")` does; the
  latent linear predictor is `predict(type = "link")`.
- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  on an ordinal fit draws one probability curve per response category,
  which is what brms draws under `categorical = TRUE`.
- [`coef()`](https://rdrr.io/r/stats/coef.html) is the per-group
  convention shared by all three packages: fixed effects plus
  conditional modes per grouping level.
- `plot(fit)` draws residual diagnostics, not MCMC traces; for
  simulation-based residual checks use
  [`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md).

## Seeing the model code

brms users inspect the generated model with `stancode(fit)`. The frmtmb
analog is the objective function itself: an ordinary R closure that
`RTMB::MakeADFun()` tapes, rebuilt deterministically from the fitted
frame.

``` r

nll <- frmtmb:::build_objective(fit$frame)
nll(fit$estimates)   # joint negative log density at the estimates
```

Two properties matter when reading it:

- The closure takes the named parameter list (the shape of
  `fit$frame$par_template`) and returns the JOINT negative log density:
  random-effect prior terms are inside it. The marginal (Laplace)
  likelihood the optimizer minimized is what `MakeADFun(random = )`
  makes of it, exposed as `fit$obj$fn` over the flat outer-parameter
  vector; at the optimum `fit$obj$fn(fit$opt$par)` equals
  `-logLik(fit)`.
- All data is baked into the closure’s environment (the reason refits
  re-tape instead of recompiling), so `environment(nll)` holds the
  assembled frame, and the closure evaluates on plain numeric parameter
  lists without any tape.

`RTMB::MakeADFun(nll, fit$frame$par_template, random = "b", map = fit$frame$map)`
reproduces `fit$obj`. There is no generated code to print: the closure
IS the model, and [`body()`](https://rdrr.io/r/base/body.html),
[`environment()`](https://rdrr.io/r/base/environment.html), and the
debugger work on it like on any R function.

## When you still want brms

Genuine prior information, full posterior uncertainty for derived
quantities, models with discrete latent structure beyond
observation-level mixtures, or `loo`-based model comparison. A practical
workflow is model screening in frmtmb and a final fit in brms. The
companion package **frmtmb.sample** is the third option:
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.html)
runs NUTS on the frmtmb objective under brms’s default priors and gives
back the posterior method surface, and
[`frmtmb.sample::as_tmbstan()`](https://aforren1.github.io/frmtmb/reference/as_tmbstan.html)
hands the objective to tmbstan directly.
