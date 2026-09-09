# Changelog

## frmtmb.sample 0.4.0

BREAKING. `frm_sample(control =)` is the SAMPLER’s control list now,
brms’s spelling and brms’s meaning, passed to tmbstan unchanged. The
fit-time options move to `fit_control =`.
[`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
already meant `control` that way, so the two entry points agree rather
than offering a third convention, and brms code that tightens
`adapt_delta` ports across unchanged.

- That closes a silent discard. On the fitted-object route
  `control = list(adapt_delta = 0.99)` was IGNORED:
  `stan_args[[1]]$control` came back NULL where 0.97 was asked for, and
  the mean `accept_stat__` was 0.9266 against 0.9836 once the setting
  arrives. `fit_control`, `start`, `data2`, `na.action` and `REML` were
  silently ignored on that route too, and are refused by name. An
  abbreviation is refused as well: `contro =` fell into `...`, reached
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  and partial matched onto ITS `control`, bypassing every check.

- [`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
  no longer returns an empty `stanfit` without saying so. That was worse
  than a missing refusal: it bricked the session, and every later model
  failed with `empty_nested() must be true`. Both doors share one check,
  and a run that produced nothing is worded apart from a model that
  cannot be sampled at all, because those want different fixes.

- [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  reads the compatibility registry before taping and refuses a fit the
  sampler cannot run, in about half a second with rstan never loaded,
  against the roughly three seconds the failure used to take. Where the
  formula position settles it the match is made on position; where it
  cannot, the match is made on names, and
  [`?frm_sample`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  says that is the largest set justifiable that way rather than a sound
  one.

- **`frm_sample(control =)` is the SAMPLER’s control list now, which is
  what brms means by the name.**
  `frm_sample(fit, control = list(adapt_delta = 0.99, max_treedepth = 12))`
  ports straight from a `brm()` call and reaches rstan through tmbstan.
  The fit-time options
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) takes
  are `fit_control =`, and they still apply only on the formula route.
  This is a breaking rename, and the old shape is refused rather than
  reinterpreted: a
  [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.html)
  list arriving in `control` is caught by the field names it carries and
  the refusal names `fit_control`. The two vocabularies share no name,
  which is what makes the two lists decidable, and a test asserts they
  stay disjoint. An option rstan does not have is also refused by name,
  because rstan answers an unknown one by declining to sample and
  returning an empty fit, which used to be reported as a solver failure.
  So is an ABBREVIATION of the name: `control` follows `...` in the
  formals, so only the exact spelling binds to it, and a shorter one
  would land in `...` and then partial-match rstan’s own `control`,
  reaching the sampler unchecked.
  [`vignette("brms-posterior")`](https://aforren1.github.io/frmtmb/frmtmb.sample/articles/brms-posterior.md)
  carried a row saying `adapt_delta` had no route through this function;
  it now says the opposite.

- **[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  refuses an assembly argument on a fitted object instead of accepting
  and discarding it.** `fit_control`, `start`, `data2`, `na.action` and
  `REML` build a model, and a fitted object is already built, so on that
  route they were read by nothing and said nothing:
  `frm_sample(fit, REML = TRUE)` sampled the model it had. `data` and
  `family` were already refused there; these five now are too, and the
  message points at `control` for the caller who meant the sampler.
  `fit_control` matters most, because it is the name the `control`
  refusal sends a confused user to.

- **Breaking, and a return-contract change:
  [`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
  never hands back an empty `stanfit` now.** rstan answers some failures
  by declining to sample, printing its reason to the console and
  returning the shell of a `stanfit` with no error at all.
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  has always refused that;
  [`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
  passed it on, so a failed run looked like a successful one until an
  accessor came back empty. Both now raise the same message, naming the
  function that raised it. The wording distinguishes the two ways to get
  nothing out of this package, because they want different fixes: a
  model that cannot be sampled at all is refused by the compatibility
  pre-flight before any chain starts, and this one is a failure of the
  run in front of you. Nothing in the package or its vignettes relied on
  the empty object.

- **Fixed: a sampler control list handed to a fitted object was
  discarded without a word.** Up to and including 0.3.2,
  `frm_sample(fit, control = list(adapt_delta = 0.99))` on the fit route
  ran the chains with rstan’s own defaults and said nothing: `control`
  was read only where the formula route assembles a model, and the fit
  route does not go there. A brms user porting a call that raised
  `adapt_delta` therefore got the untightened sampler. The argument now
  reaches rstan on both routes, and `ds$stanfit@stan_args[[1]]$control`
  records what was asked for.

- **A pre-flight against the compatibility registry**, on
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  and on
  [`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md).
  Both read the rows that name `frm_sample` with the status `"refused"`
  before anything is taped, before the formula route calls
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html), and
  before the rstan and tmbstan namespaces are loaded, so a model the
  sampler cannot run costs the registry read and nothing else: 0.36 s to
  0.42 s on the first call of a session against the 2.9 s failure it
  replaces, and that failure also left rstan unusable for the rest of
  the session. The refusal repeats the registered note and names the
  row. `frmtmb.ode` registers the first such row.

  What the check refuses on is deliberately narrow, because what
  identifies a feature in a model is WHERE its call sits in the formula.
  It decides three ways. A response FAMILY is matched by name, which is
  exact; it does not look inside
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html).
  An ADDITION TERM is matched by POSITION, since it is written inside
  the bar on the left of the response formula and nowhere else, which
  separates `y | se(v) ~ x` from `y ~ a * se(x)` with a user’s own
  `se()`. Everything else the vocabulary writes as a call is matched by
  NAME, outside addition-term position, and only when the name is not
  also a base R function; that route is justifiable from names but not
  sound, so a refusal registered on `mo()`, `ma()`, `cosy()` or `mm()`
  would also fire on a user’s function of that name. Everything the
  check cannot decide warns and samples rather than refusing a model
  that may be correct. A registered refusal cannot be revoked by a later
  registration: display names are unique, core refuses a second
  registration of one, and a name the vocabulary writes without
  parentheses is declared not to be a call at all.

## frmtmb.sample 0.3.2

Requires frmtmb 0.55.0, for the hazard-container lint that now runs in
this package’s own check.

## frmtmb.sample 0.3.1

- The sampling articles say why they are unevaluated, and the reason
  they used to give was false: four chains on the model in question take
  seconds, not the minutes claimed. They are unevaluated because the
  samplers are optional dependencies and because that page documents a
  Stan build which samples the wrong density silently. The two chunks
  that need no sampler now run, so the page shows output.

- [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) on a
  group-unit matrix says which unit it left out.
  [`loo::loo.matrix()`](https://mc-stan.org/loo/reference/loo.html)
  never sees the attribute the log-likelihood matrix carries, so the
  printed estimate was indistinguishable from a per-observation one.

## frmtmb.sample 0.3.0

[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html)
on draws returns core’s frame, grid and display quantity;
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html)
reports Savage-Dickey evidence ratios validated against brms;
[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) says
which unit it leaves out. Requires frmtmb 0.53.0 for the seam exports.

- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html)
  on a draws object now returns the frame core’s fit method returns,
  which closes the divergence 0.2.0 recorded. The columns and their
  order are core’s
  [`ce_frame()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html):
  the varied predictor(s), then every other model variable at the value
  it is held at, then `cond__`, `cats__`, `effect1__` and `effect2__`,
  then the band. `cond__` is always present as a factor over every
  condition level, so brms’s own
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) has something
  to facet on; a per-category ordinal display is keyed `"x:cats__"`
  rather than `"x"`, so brms’s
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) reads the
  pair back out of the effects attribute; and `int_conditions =`,
  `categorical =` and `seed =` are formals here as they are on the fit
  method. A ported script that indexed the returned frame positionally,
  or that read `ce$x` on an ordinal fit, has to be read again.

- BUG FIX.
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html)
  on draws drew the `mu` predictor where the fit method draws the
  EXPECTED RESPONSE, on every family whose mean is not the inverse link
  of `mu` and on every [`trunc()`](https://rdrr.io/r/base/Round.html)
  response. A zero-inflated fit was plotted `1 / (1 - zi)` times too
  high, a mixture showed one component’s mean instead of the mixture’s,
  and a response truncated below at zero reported a NEGATIVE mean and
  drew it. The column names, the grid and the keys were all correct
  throughout, which is why the first round’s identity checks passed. The
  display quantity now comes from core’s own
  [`ce_pred_dpar()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  so the two surfaces cannot make different choices again. Measured
  after the fix, largest deviation from the fit curve over a five-point
  grid: 1.4% zero-inflated, 1.0% truncated, 0.5% mixture, against 21.7%
  and a sign change before it. An explicitly named `dpar =` was always
  correct and is unchanged.

- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html)
  on draws accepts `allow_new_levels` and lme4’s `allow.new.levels`,
  which the fit method has always accepted and this method warned about
  as unknown. It reads them with core’s own
  [`ce_dots()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  so the accepted set cannot drift again.

- [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html)
  refuses `abs(x) = 0` instead of reporting a Bayes factor twice too
  large. The affine check probed only non-negative points, where
  [`abs()`](https://rdrr.io/r/base/MathFun.html) is the identity; it now
  probes both signs. The prior of `|X|` for `X ~ N(0, 1)` is a
  half-normal with twice the density at zero, so the denominator was
  half what it should have been, silently.

- The `"evid_ratio_mcse"` attribute is measured over at least four
  blocks of draws. With two chains it was the standard deviation of two
  numbers over `sqrt(2)` - one degree of freedom - and understated a
  measured 0.46 gap between two runs of one model as 0.18. Chains are
  split in half when there are fewer than four, which keeps each block
  consecutive and so keeps the autocorrelation the estimator exists to
  respect.

- A point hypothesis on a natural-scale summary says so. `sigma = 0.9`
  used to refuse with “not a population-level coefficient”, which is
  true of the name and false of the model: `sigma` IS a coefficient
  here, reported on a scale its prior is not written on. That is now the
  reason given, for every such quantity.

- BEHAVIOR CHANGE. `conditional_effects(re_formula = NULL)` on draws
  conditions on a NEW group, as it does on a fit, instead of silently
  taking the first observed level. The grid’s grouping column reports
  `NA` and each POSTERIOR DRAW gets that group’s effects drawn from its
  own covariance parameters, which is brms’s construction; the fit
  method draws them per bootstrap replicate instead. The band is wider
  than the population band by the amount the new group’s variance adds
  (1.92 to 3.93 in mean width on the model in
  `dev/sample-ce-findings.md`). `seed =` makes the draw reproducible.

- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html)
  on draws reports an argument it cannot use instead of discarding it in
  silence.

- [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html)
  on a draws object reports `evid_ratio` and `post_prob`, brms’s
  `Evid.Ratio` and `Post.Prob`. A directional hypothesis gets the
  posterior odds of the claim, which needs no prior. A POINT hypothesis
  gets the Savage-Dickey density ratio: the posterior density of the
  tested quantity at the point, from a kernel density with
  [`stats::density()`](https://rdrr.io/r/stats/density.html)’s default
  `nrd0` bandwidth over 4096 grid points read by spline, over the PRIOR
  density there, evaluated from the specification the model was sampled
  under rather than estimated from prior draws. brms estimates both
  halves by kernel density from `sample_prior = "yes"` draws; the two
  agree to the accuracy of that estimate, and the numbers are in
  `dev/sample-ce-findings.md`. The between-chain Monte Carlo error of
  each ratio rides on the returned object as the `"evid_ratio_mcse"`
  attribute, because a density ratio read without one is easy to
  over-believe.

- A point hypothesis whose parameter has no proper prior REFUSES BY
  NAME, with a warning that says which parameter and why, and leaves
  `evid_ratio` `NA`; the rest of the row is unaffected and the other
  hypotheses of the same call still report.
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  leaves class `"b"` flat by default exactly as brms does, so this is
  the commonest answer until a prior is written. A hypothesis on a
  variance component, on a nonlinear function of the coefficients, or on
  the intercept (whose prior is about the intercept at the predictor
  means, not about the coefficient) refuses for its own stated reason.

- [`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md),
  [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) and
  [`waic()`](https://aforren1.github.io/frmtmb/reference/loo.html)
  accept a structured family that declares how its likelihood factorizes
  (`frmtmb::frmtmb_structure(loglik_group = )` or `(loglik_row = )`),
  instead of refusing every family whose likelihood is not rowwise. The
  columns are the pieces the family declares, at the coarsest
  granularity it gives, because a family that groups its likelihood is
  saying its rows are not independently droppable. The matrix carries
  `attr(x, "unit")` naming what a column is when it is not an
  observation. A family that supplies one number for the whole response
  is still refused, in the words it already used.

## frmtmb.sample 0.2.0

- BEHAVIOR CHANGE, following core.
  [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.html) on
  a `frmtmb_draws` object is keyed by the GROUPING FACTOR rather than by
  the block, so `ranef(draws)[["1 | g"]]` returns `NULL` where it used
  to return an array and `ranef(draws)$g` returns it instead. The method
  delegates to core’s
  [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.html)
  once per draw, so it follows core’s re-key exactly; the block label
  rides along in each array’s `"term"` attribute, as it does in core. A
  model with two blocks on ONE factor gives two entries under one name,
  and the method now assembles them BY POSITION: keyed by name,
  `out[[tn]] <- st` wrote the same name twice, so the SECOND block was
  dropped entirely and the list came back with one entry instead of two.

- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html)
  on a draws object inherits four grid changes from core, because it
  builds its grids with core’s
  [`ce_grids_build()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html):
  the two-variable grid now varies the FIRST effect slowest (brms’s row
  order), a numeric moderator is held at the exact `mean +/- sd` rather
  than at `signif(mean +/- sd, 3)`, a `mo()` predictor gets one grid
  point per level instead of a 100-point continuous grid, and a
  `trials()` variable is held at 1 with a message that is new here. What
  it does NOT yet inherit is core’s frame shape: the returned columns,
  the always-present `cond__`, the `"x:cats__"` key for a per-category
  display, `int_conditions =`, `categorical =` and the new-group meaning
  of `re_formula = NULL` are core-only for now, so the two surfaces
  disagree in both directions until the draws method is migrated.

- [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  refuses parallel chains on Windows when the model’s family comes from
  a namespace built by
  [`pkgload::load_all()`](https://pkgload.r-lib.org/reference/load_all.html),
  naming the package and the two remedies. A worker process cannot load
  such a namespace, and R silently substitutes the global environment
  for it, so the family failed at its first internal call and rstan
  reported only that the chains returned no draws.

- The defaults this package registers with
  [`frmtmb::get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.html)
  now answer `route = "sample"` only. Loading this package used to
  change the table
  [`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.html)
  returned for every caller, an author asking about
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
  included; it no longer does. Ask for
  `get_prior(..., route = "sample")` to see what
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  applies, and with this package unloaded that route is refused rather
  than reported as flat.
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  itself is unchanged: the same defaults under the same
  call-over-fit-over-defaults precedence.

- The frmtmb floor rises to 0.51.0. The registration in `.onLoad()` now
  passes `frmtmb_register_compat(expects = )`, which frmtmb gained in
  0.51.0, and an older frmtmb stops the load outright with
  `unused argument (expects = c("hmm", "lca"))`. The old floor of 0.46.0
  permitted an installation that could not load, so this is a hard
  requirement rather than a preference.

- `frmtmb.latent` joins `Suggests`. `hmm()` and `lca()` left frmtmb for
  that package, so the tests here that sample one of them qualify the
  call and skip when it is absent. The two compatibility rules this
  package registers for those pairs are declared in
  `frmtmb_register_compat(expects = )`, which says they name a feature
  another package owns. Before frmtmb checked, those rules were dropped
  without a word in every session that did not also load
  `frmtmb.latent`; now
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html)
  reports them as unresolved until it is loaded.

## frmtmb.sample 0.1.0

First release, extracted from frmtmb 0.46.0.

- NUTS sampling for frmtmb models through tmbstan:
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  on a fit or a formula, with brms 2.23’s default priors on both routes,
  priors stacking most-explicit-first, and non-centered sampling of
  qualifying random-effect blocks. Draws carry frmtmb parameter names.
- The full posterior method surface on core’s generics and on the brms,
  loo, posterior, bayesplot, rstantools, coda and bridgesampling
  generics:
  [`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md),
  [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html),
  [`waic()`](https://aforren1.github.io/frmtmb/reference/loo.html),
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.html),
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.html)
  on draws,
  [`bayes_R2()`](https://aforren1.github.io/frmtmb/reference/bayes_R2.html),
  [`mcmc_plot()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
  and the rest.
- [`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md)
  measures the Laplace and Wald approximations against NUTS on the same
  objective, centered and with no added priors.
- At load, the package registers its sampling default priors with core’s
  [`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.html)
  and its feature rows with the compatibility matrix, so both stay
  truthful whether or not this package is installed.
- Refuses a tmbstan built against StanHeaders 2.39 or later, whose code
  generator leaves one log-density overload unpatched so every chain
  silently samples a standard normal instead of the model.
- Parallel chains work on every platform, Windows included: tmbstan
  rebuilds the tape on each socket worker from the serialized objective
  closure, and frmtmb’s generated closures are self-contained.
