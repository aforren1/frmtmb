# frmtmb.sample 0.3.0

`conditional_effects()` on draws returns core's frame, grid and
display quantity; `hypothesis()` reports Savage-Dickey evidence ratios
validated against brms; `loo()` says which unit it leaves out. Requires
frmtmb 0.53.0 for the seam exports.

* `conditional_effects()` on a draws object now returns the frame
  core's fit method returns, which closes the divergence 0.2.0
  recorded. The columns and their order are core's `ce_frame()`: the
  varied predictor(s), then every other model variable at the value it
  is held at, then `cond__`, `cats__`, `effect1__` and `effect2__`,
  then the band. `cond__` is always present as a factor over every
  condition level, so brms's own `plot()` has something to facet on; a
  per-category ordinal display is keyed `"x:cats__"` rather than `"x"`,
  so brms's `plot()` reads the pair back out of the effects attribute;
  and `int_conditions =`, `categorical =` and `seed =` are formals here
  as they are on the fit method. A ported script that indexed the
  returned frame positionally, or that read `ce$x` on an ordinal fit,
  has to be read again.

* BUG FIX. `conditional_effects()` on draws drew the `mu` predictor
  where the fit method draws the EXPECTED RESPONSE, on every family
  whose mean is not the inverse link of `mu` and on every `trunc()`
  response. A zero-inflated fit was plotted `1 / (1 - zi)` times too
  high, a mixture showed one component's mean instead of the mixture's,
  and a response truncated below at zero reported a NEGATIVE mean and
  drew it. The column names, the grid and the keys were all correct
  throughout, which is why the first round's identity checks passed.
  The display quantity now comes from core's own `ce_pred_dpar()`, so
  the two surfaces cannot make different choices again. Measured after
  the fix, largest deviation from the fit curve over a five-point grid:
  1.4% zero-inflated, 1.0% truncated, 0.5% mixture, against 21.7% and a
  sign change before it. An explicitly named `dpar =` was always
  correct and is unchanged.

* `conditional_effects()` on draws accepts `allow_new_levels` and
  lme4's `allow.new.levels`, which the fit method has always accepted
  and this method warned about as unknown. It reads them with core's
  own `ce_dots()`, so the accepted set cannot drift again.

* `hypothesis()` refuses `abs(x) = 0` instead of reporting a Bayes
  factor twice too large. The affine check probed only non-negative
  points, where `abs()` is the identity; it now probes both signs. The
  prior of `|X|` for `X ~ N(0, 1)` is a half-normal with twice the
  density at zero, so the denominator was half what it should have
  been, silently.

* The `"evid_ratio_mcse"` attribute is measured over at least four
  blocks of draws. With two chains it was the standard deviation of two
  numbers over `sqrt(2)` - one degree of freedom - and understated a
  measured 0.46 gap between two runs of one model as 0.18. Chains are
  split in half when there are fewer than four, which keeps each block
  consecutive and so keeps the autocorrelation the estimator exists to
  respect.

* A point hypothesis on a natural-scale summary says so. `sigma = 0.9`
  used to refuse with "not a population-level coefficient", which is
  true of the name and false of the model: `sigma` IS a coefficient
  here, reported on a scale its prior is not written on. That is now
  the reason given, for every such quantity.

* BEHAVIOR CHANGE. `conditional_effects(re_formula = NULL)` on draws
  conditions on a NEW group, as it does on a fit, instead of silently
  taking the first observed level. The grid's grouping column reports
  `NA` and each POSTERIOR DRAW gets that group's effects drawn from its
  own covariance parameters, which is brms's construction; the fit
  method draws them per bootstrap replicate instead. The band is wider
  than the population band by the amount the new group's variance adds
  (1.92 to 3.93 in mean width on the model in `dev/sample-ce-findings.md`).
  `seed =` makes the draw reproducible.

* `conditional_effects()` on draws reports an argument it cannot use
  instead of discarding it in silence.

* `hypothesis()` on a draws object reports `evid_ratio` and
  `post_prob`, brms's `Evid.Ratio` and `Post.Prob`. A directional
  hypothesis gets the posterior odds of the claim, which needs no
  prior. A POINT hypothesis gets the Savage-Dickey density ratio: the
  posterior density of the tested quantity at the point, from a kernel
  density with `stats::density()`'s default `nrd0` bandwidth over 4096
  grid points read by spline, over the PRIOR density there, evaluated
  from the specification the model was sampled under rather than
  estimated from prior draws. brms estimates both halves by kernel
  density from `sample_prior = "yes"` draws; the two agree to the
  accuracy of that estimate, and the numbers are in
  `dev/sample-ce-findings.md`. The between-chain Monte Carlo error of
  each ratio rides on the returned object as the `"evid_ratio_mcse"`
  attribute, because a density ratio read without one is easy to
  over-believe.

* A point hypothesis whose parameter has no proper prior REFUSES BY
  NAME, with a warning that says which parameter and why, and leaves
  `evid_ratio` `NA`; the rest of the row is unaffected and the other
  hypotheses of the same call still report. `frm_sample()` leaves class
  `"b"` flat by default exactly as brms does, so this is the commonest
  answer until a prior is written. A hypothesis on a variance
  component, on a nonlinear function of the coefficients, or on the
  intercept (whose prior is about the intercept at the predictor means,
  not about the coefficient) refuses for its own stated reason.

* `log_lik()`, `loo()` and `waic()` accept a structured family that
  declares how its likelihood factorizes
  (`frmtmb::frmtmb_structure(loglik_group = )` or `(loglik_row = )`),
  instead of refusing every family whose likelihood is not rowwise.
  The columns are the pieces the family declares, at the coarsest
  granularity it gives, because a family that groups its likelihood is
  saying its rows are not independently droppable. The matrix carries
  `attr(x, "unit")` naming what a column is when it is not an
  observation. A family that supplies one number for the whole
  response is still refused, in the words it already used.

# frmtmb.sample 0.2.0

* BEHAVIOR CHANGE, following core. `ranef()` on a `frmtmb_draws`
  object is keyed by the GROUPING FACTOR rather than by the block, so
  `ranef(draws)[["1 | g"]]` returns `NULL` where it used to return an
  array and `ranef(draws)$g` returns it instead. The method delegates
  to core's `ranef()` once per draw, so it follows core's re-key
  exactly; the block label rides along in each array's `"term"`
  attribute, as it does in core. A model with two blocks on ONE factor
  gives two entries under one name, and the method now assembles them
  BY POSITION: keyed by name, `out[[tn]] <- st` wrote the same name
  twice, so the SECOND block was dropped entirely and the list came
  back with one entry instead of two.

* `conditional_effects()` on a draws object inherits four grid changes
  from core, because it builds its grids with core's
  `ce_grids_build()`: the two-variable grid now varies the FIRST effect
  slowest (brms's row order), a numeric moderator is held at the exact
  `mean +/- sd` rather than at `signif(mean +/- sd, 3)`, a `mo()`
  predictor gets one grid point per level instead of a 100-point
  continuous grid, and a `trials()` variable is held at 1 with a
  message that is new here. What it does NOT yet inherit is core's
  frame shape: the returned columns, the always-present `cond__`, the
  `"x:cats__"` key for a per-category display, `int_conditions =`,
  `categorical =` and the new-group meaning of `re_formula = NULL` are
  core-only for now, so the two surfaces disagree in both directions
  until the draws method is migrated.

* `frm_sample()` refuses parallel chains on Windows when the model's
  family comes from a namespace built by `pkgload::load_all()`, naming the
  package and the two remedies. A worker process cannot load such a
  namespace, and R silently substitutes the global environment for it,
  so the family failed at its first internal call and rstan reported
  only that the chains returned no draws.

* The defaults this package registers with `frmtmb::get_prior()` now
  answer `route = "sample"` only. Loading this package used to change
  the table `get_prior()` returned for every caller, an author asking
  about `frm()` included; it no longer does. Ask for
  `get_prior(..., route = "sample")` to see what `frm_sample()`
  applies, and with this package unloaded that route is refused rather
  than reported as flat. `frm_sample()` itself is unchanged: the same
  defaults under the same call-over-fit-over-defaults precedence.

* The frmtmb floor rises to 0.51.0. The registration in `.onLoad()`
  now passes `frmtmb_register_compat(expects = )`, which frmtmb gained
  in 0.51.0, and an older frmtmb stops the load outright with
  `unused argument (expects = c("hmm", "lca"))`. The old floor of
  0.46.0 permitted an installation that could not load, so this is a
  hard requirement rather than a preference.

* `frmtmb.latent` joins `Suggests`. `hmm()` and `lca()` left frmtmb
  for that package, so the tests here that sample one of them qualify
  the call and skip when it is absent. The two compatibility rules this
  package registers for those pairs are declared in
  `frmtmb_register_compat(expects = )`, which says they name a feature
  another package owns. Before frmtmb checked, those rules were dropped
  without a word in every session that did not also load
  `frmtmb.latent`; now `frm_compat()` reports them as unresolved until
  it is loaded.

# frmtmb.sample 0.1.0

First release, extracted from frmtmb 0.46.0.

* NUTS sampling for frmtmb models through tmbstan: `frm_sample()` on
  a fit or a formula, with brms 2.23's default priors on both routes,
  priors stacking most-explicit-first, and non-centered sampling of
  qualifying random-effect blocks. Draws carry frmtmb parameter
  names.
* The full posterior method surface on core's generics and on the
  brms, loo, posterior, bayesplot, rstantools, coda and
  bridgesampling generics: `log_lik()`, `loo()`, `waic()`,
  `posterior_epred()`, `posterior_predict()`, `hypothesis()`,
  `conditional_effects()` on draws, `bayes_R2()`, `mcmc_plot()`, and
  the rest.
* `check_laplace()` measures the Laplace and Wald approximations
  against NUTS on the same objective, centered and with no added
  priors.
* At load, the package registers its sampling default priors with
  core's `get_prior()` and its feature rows with the compatibility
  matrix, so both stay truthful whether or not this package is
  installed.
* Refuses a tmbstan built against StanHeaders 2.39 or later, whose
  code generator leaves one log-density overload unpatched so every
  chain silently samples a standard normal instead of the model.
* Parallel chains work on every platform, Windows included: tmbstan
  rebuilds the tape on each socket worker from the serialized
  objective closure, and frmtmb's generated closures are
  self-contained.
