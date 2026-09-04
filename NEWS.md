# frmtmb 0.47.0

The monorepo split: the ODE seam and the sampling surface leave core
for companion packages in this repository.

* The sampling surface moved to the companion package
  `frmtmb.sample`: `frm_sample()`, `as_tmbstan()`, `check_laplace()`,
  and every `frmtmb_draws` method, including `log_lik()`, `loo()`,
  `posterior_epred()` and `conditional_effects()` on draws. Install
  it with `remotes::install_github("aforren1/frmtmb", subdir =
  "extensions/frmtmb.sample")`. Absent it, `frm_sample()` is not
  found, and the fit-side methods that need draws name the package
  and the call in their refusal. Core keeps every generic it still
  defines a method for, so one `loo()` serves both packages in a
  session.
* Ordinary differential equation dynamics moved to the companion
  package `frmtmb.ode` (`frm_ode()`, its tests, its vignette). With
  it go the RTMBode suggestion and the `Additional_repositories`
  line: core's dependencies now resolve entirely from CRAN.
* New public interface `?"frmtmb-sampling-api"`: the internals a
  sampling extension is written against, as one documented contract.
  New extension seams: `frmtmb_register_frame_check()` (a frame
  check an extension registers at load), `frmtmb_register_compat()`
  (compatibility-matrix rows), a prior-defaults registry behind
  `get_prior()`, and `frmtmb_ad_overload()` (the tape-safe wrapper
  for user functions, previously internal).
* `get_prior()` reports the default a slot actually has: `(flat)`
  with frmtmb alone, which is what `frm()` does, and the sampling
  defaults once `frmtmb.sample` registers them.
* `?frm` gains "The Laplace approximation, and how to check it",
  naming the in-package remedies and pointing at
  `frmtmb.sample::check_laplace()` for the direct measurement.
* Core's boundary test now asserts zero: no core file names a
  structured family or an ODE symbol outside the two family homes
  that remain until they move out in turn.
* `Suggests`: `coda` and `bridgesampling` dropped; no core code
  registers on them any more.

# frmtmb 0.46.0

The structured-family protocol completes in-package: hmm and lca ride
`fam$structure` alone, and the extension API is exported.

* New exported `latent_probs()` returns the posterior latent-state
  probabilities of any family that declares them. `mixture_probs()`,
  `hmm_probs()` and `lca_probs()` are unchanged and now reach it;
  `hmm_viterbi()` stays its own decoding pass.
* New extension API on one page, `?"frmtmb-extension-api"`:
  `single_response()`, `eval_dpars()`, `fit_extras()`,
  `dpar_linpred()`, `response_mean()` and `as_frmtmb_family()` are
  the read-only accessors a `frmtmb_structure()` slot may use, with a
  stability promise. `uni_resp()` was internal and is now
  `single_response()`.
* `frmtmb_structure()` gains `check_fit`, `unit` and a
  `cluster_robust` capability flag, and `loglik` becomes optional. A
  structure with no `loglik` is a capability declaration: the family
  keeps its rowwise likelihood and carries only the refusals that
  belong to it, and its capability flags default to TRUE so it names
  what it refuses instead of enumerating what already works.
* No core file branches on `hmm()` or `lca()` any more.
  `fit$frame$hmm_g` is gone; a hidden-Markov structure's data is
  `fit$frame$blocks[[response]]`. The compatibility matrix takes both
  families' rows from the families themselves through a registration
  the boundary test polices.
* One wording change: each mixture-type family now states the REML
  and `profile = TRUE` refusals in its own name. The message
  previously listed `mixture()`, `mixture_mvn()` and `lca()`
  whichever one was fitted; everything after the colon is unchanged.
* Fits are bit-identical to 0.45.0, and `frm_compat()` resolves every
  feature pair to the same status.

# frmtmb 0.45.0

The structured-family protocol: the first step of the core and
extensions split.

* New exported `frmtmb_structure()` declares a likelihood that does
  not factorize over the rows of the data. It carries the non-rowwise
  log-likelihood, the frame block that likelihood reads, the frame
  variables and `NA` policy it needs, a conditional mean and
  variance, a structured simulator, and twelve capability flags, each
  with the message its refusal shows. `frmtmb_family()` gains
  `structure =`. A structured family can now be written outside the
  package; `dev/structured-family-protocol.md` is the contract.
* `mixture(groups = )` and `hmm()` reach the core through that one
  slot. The objective's response loop has one branch for both instead
  of one each, the model frame carries one `blocks` slot, and
  `predict()`, `fitted()`, `residuals()`, `simulate()`,
  `conditional_effects()` and the fit-time REML, quadrature and
  profile gates read capability flags instead of naming a family.
  Fits are bit-identical to the previous branches, and no user-facing
  message changed; the only deleted refusals were unreachable behind
  generic guards, which an existing test asserts.
* `fit$frame$mix_g` is gone. A group-level mixture's grouping is
  `fit$frame$blocks[[response]]`, and a structured simulator reads
  `ctx[["block"]]` where it read `ctx[["mix_g"]]`.
* Internal reorganization toward the split: the prior and bound
  machinery the fit route reaches moved from `R/interop.R` to
  `R/priors.R`, and a new boundary test pins where core still names
  the structured families and the ODE seam, as a one-way ratchet the
  remaining protocol steps shrink to empty. No user-visible change.

# frmtmb 0.44.0

Sampling a fit gets the default priors, conditional-effects plots get
real facets, and parameter names become discoverable.

## Default priors on the fit route

* `frm_sample(fit)` now applies the same brms 2.23 default priors the
  formula route applies. Both routes sample a posterior, and a
  prior-free `brm()` call translated onto a fitted model no longer
  pins the chain at whatever boundary mode maximum likelihood found.
  On the kidney model of `dev/brms-vignette-audit.md`, `sd(patient)`
  goes from 3.1e-4 with a 3e-6-wide interval to 0.49 with a
  [0.04, 1.14] interval, and the largest R-hat from 1.57 to 1.10.
  This is a behavior change for existing scripts that sampled a fit:
  `prior = "flat"` restores the bare likelihood.
* Priors on a sampling call now stack, most explicit winning per
  slot: `prior =` on the call, then the fit's own prior, then the
  defaults. Previously an explicit `prior =` on a MAP fit discarded
  that fit's prior entirely.
* Because the defaults cover every standard deviation and
  correlation, the non-centered parameterization is now offered on
  the fit route as well; `reparameterize = FALSE` still asks for the
  centered chain.
* `prior = "flat"` is the opt-out on both routes, and
  `frm_sample(fit, prior = "flat")` is no longer an error. On a MAP
  fit it also removes that fit's own penalty, which is taped into its
  objective, and the call names the prior it dropped. A MAP penalty
  can never be applied twice: every prior-carrying path rebuilds the
  objective from the bare likelihood and adds the merged priors once,
  and a regression asserts it in objective values.
* `check_laplace()` is unchanged in what it measures: it samples the
  density the fit maximized, with no default priors, centered. It now
  asks for that explicitly rather than inheriting `frm_sample()`'s
  default.

## Faceted conditional-effects plots

* `plot()` on a `conditional_effects()` result with several condition
  sets draws one faceted page of small multiples with a shared scale,
  as brms's `facet_wrap("cond__")` does, instead of one full page per
  condition that overwrote its predecessor on a normal device. The
  panels use tinyplot when it is installed; without it the fallback
  is a grid of base-graphics panels, still one page. `ncol =` is
  honored; it was previously absorbed by `...` and dropped.
* `plot(ce, points = TRUE)` shows in each panel only the observations
  belonging to that condition, following brms's `make_point_frame()`:
  a factor condition splits the data per panel, while a numeric
  condition is a reference value that names no observation and keeps
  every point. One deliberate divergence: brms pins factors the
  `conditions` frame does not mention to their first level, silently
  discarding the other levels' observations; frmtmb subsets only on
  the variables the user passed.

## Parameter names, discoverable and usable

* New `par_template()`: the parameter vector's names and starting
  values, for a fitted model or for a formula and data alone. It
  prints as a compact component listing, and the object it returns is
  edited and handed straight back as `frm(start =)` or
  `frm_simulate(newparams =)`. brms has no counterpart: `brm(init =)`
  takes Stan program names and brms initializes at random instead.
* `frm(start =)` matches by name. A named vector inside a component
  overrides those entries and leaves the rest at their defaults, with
  parentheses optional as everywhere a parameter is named. An unnamed
  vector keeps the positional contract exactly. Mixing the two in one
  vector, an unknown name, a name given twice, and a non-numeric
  vector are each refused by name. Note the edge: a full-length named
  vector previously had its names silently ignored; names are now
  interpreted, and a vector whose names match nothing is refused with
  a pointer at `unname()`.
* A `normal()`, `student_t()` or `cauchy()` prior on a nonlinear
  parameter now places its starting value at the prior location,
  which is how brms uses the same information. The insurance-loss
  growth curve of the brms nonlinear vignette fits from its priors
  with no `start` at all. A message names what was placed; `start`
  still wins. Other parameters keep the starts they had.
* Fixed on the way: `autoscale = TRUE` silently recycled a too-short
  start vector into the pre-fit's full-length start. It now goes
  through the same validation every other start does, and refuses.

# frmtmb 0.43.0

The prior interface speaks brms, ar1() evaluates in linear time, the
brms vignettes are measured call by call, and parallel chains work on
Windows.

## The prior interface

* `frm()`, `frm_sample()` and `frm_simulate()` take `prior =`, brms's
  spelling. The `priors =` of earlier releases is renamed, not
  aliased: no ecosystem spells it `priors`, and the package is
  pre-release. A stale `frm(priors =)` fails as an unused argument;
  `frm_sample()`, whose `...` would otherwise pass the name to
  tmbstan and fit with no priors at all, refuses it by name. The fit
  object's field is `fit$prior`, as `brmsfit$prior` is.
* New `prior()`, `prior_()` and `prior_string()`, brms's
  constructors. `prior(normal(5000, 1000), nlpar = "ult")` quotes its
  first argument and reaches the same machinery `set_prior()` does.
  `tag` and `check` are deliberately absent: both are Stan-program
  concerns.
* A `brmsprior` object built by brms itself is translated row by row
  wherever a prior is accepted, so a specification copied out of a
  brms script works even when brms is attached and its `prior()`
  masks frmtmb's. Rows with no faithful frmtmb meaning are refused by
  name rather than turned into a different density: a `tag`, brms's
  mixture `theta` (frmtmb's `theta` is the raw covariance vector),
  and the per-dpar classes such as `sigma`, whose density brms puts
  on the parameter and frmtmb's nearest spelling puts on the link
  scale. brms's own default rows are dropped with a message, because
  frmtmb chooses its own defaults.
* `set_prior()` gains `nlpar` and `resp`, and its argument order now
  follows brms's (`prior, class, coef, group, resp, dpar, nlpar, lb,
  ub`). Under `nlpar`, class `"b"` covers a nonlinear parameter's
  whole coefficient vector, intercept included, because its
  sub-formula is not centered; that is what makes the brms nonlinear
  vignette's insurance-loss priors land. `nlpar` also narrows classes
  `"sd"` and `"cor"`, separating two blocks on one grouping factor
  that `group` alone cannot tell apart. `get_prior()` gains an
  `nlpar` column.
* The formula route of `frm_sample()` no longer invents a location
  default for a nonlinear parameter's intercept. brms leaves those
  flat, and the response's median and mad describe nothing about a
  rate or a shape inside a nonlinear body; the disclosure message
  names the parameters left flat.

## Linear-time ar1()

* `ar1()` and `hetar1()` evaluate their density at linear cost in the
  block dimension instead of building a dense `d x d` covariance and
  factorizing it on the tape. The density, and therefore every fitted
  value, standard error and prior, is unchanged: the new form agrees
  with the old computation to 2.3e-12 in value and 2.2e-9 in
  gradient, gated by a test that keeps the old dense code inline.
  The speedup is orders of magnitude once blocks are large
  (milliseconds where the dense form took seconds at d = 800, and
  tens of seconds at d = 2000); `dev/sparsear1-bench.R` measures it
  on the host at hand.
* A long time series is now an ordinary `ar1()` block. The
  multivariate stochastic-volatility model of Skaug and Yu (2014)
  fits at its published 945 time points in seconds, agreeing with the
  upstream TMB reference to 8.5e-9, with a slightly better optimum
  than the reference's own cold start reaches.

## The brms vignettes, measured

* `dev/brms-vignettes/` translates the brms vignettes call by call,
  each model down both paths: the maximum-likelihood fit with its
  post-processing surface, and `frm_sample()` with the posterior
  workflow. Every one of 567 calls carries a label (clean, spelling,
  behavior, missing, or refusal), and
  `dev/brms-vignette-audit.md` holds the measured scoreboard: about
  7 of 10 model calls port unchanged, about 4 of 10 of all calls.
  README's port claim now states both figures.
* The audit's cheapest finding is applied: the draws and fit surfaces
  refuse with named replacements where they previously failed bare.
  `loo()`, `waic()`, `bayes_R2()`, `LOO()` and `WAIC()` on a
  maximum-likelihood fit now say to sample first or use `AIC()`;
  `loo(a, b)` says comparison belongs to `loo_compare()`;
  `plot()` on draws names `mcmc_plot()`; `update()` on draws says to
  update the fit and re-sample; `rescor_matrix()` on draws names the
  fit; `expose_functions()` on a fit gets the same refusal draws
  already had.
* `conditional_effects(method = "predict")` now works on a nonlinear
  predictor. The guard that refused it asked for a standard error
  that `method = "predict"` never uses: its band is a quantile of
  simulated responses. This is what the brms nonlinear vignette's
  facet-per-year display calls; the same call is also about a third
  faster on every model, because the dead standard-error computation
  is gone.

## Parallel chains on Windows

* `frm_sample(cores = )` parallelizes chains on Windows instead of
  falling back to sequential sampling. The old guard claimed socket
  workers cannot evaluate the RTMB tape; in fact tmbstan rebuilds the
  tape on each worker from the serialized objective closure, and
  frmtmb's generated closures are self-contained, so the rebuild
  works and the draws are identical to a sequential run at the same
  seed. The remaining cost is a fixed several seconds of startup per
  worker, which a message now states; short chains gain nothing,
  long chains approach a per-chain speedup. One caveat: a
  `devtools::load_all()` development namespace cannot be rebuilt on
  a worker, so parallel chains need the installed package.

# frmtmb 0.42.0

The TMB-examples replication audit, and REML's semantics documented
and validated beyond the gaussian case.

## The tmb_examples audit

* Audited all 35 models in the upstream `tmb_examples` directories of
  RTMB and adcomp: 13 replicate through `frm()` at 1e-8 or better,
  including the SPDE Weibull survival model with censoring,
  orange_big's nonlinear growth with 5000 latent effects, socatt's
  ordinal GLMM, lr_test's `map =` restrictions written as plain
  formulas, and a million-row regression; 3 are expressible but
  pinned or cost-capped; 19 name the capability that blocks them.
  `dev/tmb-examples-audit.md` carries the table, and every
  replicating example is a self-contained regression in
  `test-tmb-examples.R` with its reference likelihood written inline.
* The Kalman-filter node the roadmap held open is closed by
  measurement rather than opinion: the filter and the Laplace route
  agree to 1e-13 on the multivariate random walk, both are linear in
  series length, and `MakeTape()$atomic()` collapses the filter to
  nine tape nodes with no adjoint to derive, so the node is free if
  ever wanted. What blocks state-space models is the formula grammar.
  The ranked roadmap that falls out puts a sparse `ar1()` covariance
  first: the dense block is cubic in its dimension, and it is the one
  thing keeping the multivariate stochastic-volatility example from
  fitting at its published size.
* Two RTMB traps recorded for nonlinear formula authors: a bare
  `pnorm()` or `qgamma()` in an `nl` body resolves lexically to
  stats (qualify with `RTMB::`), and `Vectorize()` strips the
  automatic-differentiation class from its arguments.

## REML beyond gaussian

* The `REML` documentation now states what the flag computes: exact
  Patterson-Thompson REML for gaussian models, covariate-dependent
  sigma included, and for every other family the Laplace-approximated
  integrated likelihood under a flat prior on the mean coefficients,
  first-order equivalent to the Cox-Reid adjusted profile likelihood.
  Distributional-parameter coefficients are deliberately not
  integrated, matching nlme's varFunc estimation and the double-GLM
  literature.
* Validated three ways in `test-reml-nongaussian.R`: glmmTMB
  agreement for poisson and bernoulli GLMMs at 1e-5, the Cox-Reid
  identity verified by hand at the pinned REML optimum, and
  `gls(method = "REML")` exactness at 1e-5 for varIdent- and
  varExp-shaped sigma models. The `anova()` design rule is asserted
  for non-gaussian REML fits too.

# frmtmb 0.41.0

The tmbstan wrong-density defense, both argument dialects on the
draws surface, and the pkgcheck container as a local pre-push gate.

## The tmbstan defense

* `frm_sample()` and `as_tmbstan()` refuse, before any sampling, a
  tmbstan installation that silently samples the wrong density.
  tmbstan's install-time generator patches only the first
  `std_normal_lpdf` placeholder in the stanc output; stanc 2.39.0
  (shipped by StanHeaders 2.39.1, on CRAN 2026-09-02) emits two
  log-density overloads, and HMC reads the unpatched one, so every
  chain samples a standard normal instead of the model, likelihood
  and priors alike. The check is one cached read of the installed
  `model.hpp`, and the refusal names the defect and the remedies (a
  prebuilt tmbstan binary, or a rebuild against StanHeaders 2.32.10).
  Affected in the wild: Linux CRAN installs of tmbstan built against
  StanHeaders 2.39, and source installs elsewhere; prebuilt binaries
  are safe. `dev/prior-dropping-investigation.md` carries the proof
  and the upstream issue text.
* Three CI cycles of chain-agreement failures in pkgcheck's container
  were this bug, not chain luck: the container reproduces them
  deterministically, and a one-line patch to tmbstan's generator
  makes it reproduce this machine to every printed digit. The gates
  rationale in the test helpers is corrected accordingly; the ungated
  correctness assertions are what caught it.

## Both dialects on the draws surface

* `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`,
  `predictive_interval()`, `predictive_error()` and `pp_check()` gain
  brms's `re_formula` alongside the `re.form` they have shipped with.
  Either spelling alone is the setting, both together is refused
  naming the pair as one setting, and every default is unchanged.
  brms's own `posterior_epred.brmsfit()` carries the same two
  spellings. One consequence for positional calls: `re.form` now sits
  one position later, so a value landing there by position meets the
  both-set refusal loudly rather than silently changing meaning.
* The random-effect switch now takes effect IN-SAMPLE for
  `posterior_predict()`, `predictive_interval()`,
  `predictive_error()` and `pp_check()` on draws: it used to be
  consulted only under `newdata`, so `re_formula = NA` on training
  rows silently kept every random effect. A structured family (hmm,
  group-level mixtures, residual correlation) refuses the switch by
  name, since the structured draw is the group-level content it
  would remove.
* `pp_check()` on draws previously accepted neither spelling and
  leaked `re.form` into bayesplot as an unknown argument.
* `set_rescor()` takes brms's `rescor`; `rescor_value` stays as an
  alias. New `dev/brms-api-diff.md`: the formals of all 92 exports
  shared with brms 2.23.0, diffed and triaged.

## The pkgcheck container, locally

* `dev/run-pkgcheck-docker.sh` and `dev/pkgcheck-docker.md` run
  rOpenSci's pkgcheck-action container against the source tree on
  Windows docker, mirroring CI line for line, with a gates on/off
  switch for the sampler-agreement tests and a token-scrubbed
  dependency cache. This is the pre-push gate, and it is the
  instrument that found the tmbstan defect.

# frmtmb 0.40.1

* `conditional_effects()` gains `re_formula`, in brms's spelling and
  with brms's default: `NA` draws the population-level curve, `NULL`
  conditions on the grid's reference group levels (pick them with
  `conditions =`), and a one-sided formula keeps the named terms, on
  the fit method and the draws method alike. Previously the setting
  was hard-coded internally, so the lme4 spelling `re.form` died with
  a matched-by-multiple-arguments error on a fit and both spellings
  vanished silently into the dots on draws. A user who reaches for
  `re.form` here is now told which spelling this function takes:
  `conditional_effects()` is brms's function and speaks brms, while
  the fit surface's `predict()` and `simulate()` keep lme4's
  `re.form`. `band = "profile"` exists only for the population-level
  curve and says so.
* The message-uniqueness test failed open on one CI layout that
  offered an existing-but-empty `../../R`; it now requires positive
  identification of the package source tree and skips otherwise.

# frmtmb 0.40.0

Input validation across the export surface after an rOpenSci autotest
sweep, the submission materials themselves, and the retirement of
hand-maintained counts in favor of an enforced uniqueness property.

## Input validation

* Scalar arguments across the export surface now refuse a wrong value
  by name instead of quietly doing something else. `isTRUE()` reads a
  length-2 logical, a string, an integer and `NA` all as `FALSE`, so
  `bf(y ~ x, nl = "yes")` used to build a linear model and
  `frmtmb_control(profile = 1L)` used to turn profiling off. Nineteen
  flags are affected, along with the count, coverage and named-list
  arguments beside them, across `frm()`, `frmtmb_control()`, `bf()`,
  `mvbf()`, `set_rescor()`, `predict()`, `simulate()`, `confint()`,
  `conditional_effects()`, `frm_simulate()`, `frm_bootstrap()`,
  `frm_multiple()`, `diagnose()`, `vcov_cluster()`, `ranef()`,
  `lca()`, `cox()`, `frmtmb_family()`, `prior_normal()`, `prior_t()`
  and `check_custom_family()`. Every refusal is exercised by
  `test-input-validation.R`.
* `confint(level = )` and `conditional_effects(prob = )` refuse a
  length-2 coverage. They used to recycle it: one `confint()` table
  came back with its odd rows at 90 percent and its even rows at 95,
  with nothing in the output recording it.
* A link must be named by a single string:
  `bernoulli(link = c("log", "x"))` returned the string as if it were
  a link, and `link = 1L` selected the identity link by position.
* `frm(start = )` refuses anything that is not a fully named list; an
  unnamed value used to be silently ignored, so the fit ran from the
  defaults and reported them as the user's. `frm(dry_run = )`,
  `frm(data = NULL)`, `frm(control = )`, `frm(na.action = )`,
  `predict(newdata = )` and `predict(re.form = )` all name their
  contract instead of failing downstream.
* One behavior change worth noting: `ranef(fit, condVar = 1)` used to
  work through a truthy `if ()` and now refuses; write
  `condVar = TRUE`.

## Condition messages: property, not count

* Every condition message template in `R/` (`stop()`, `warning()` and
  `message()` alike) is asserted unique by
  `test-message-uniqueness.R`, which parses the sources at test time
  and fails on the first duplicate. The hand-maintained count and its
  ledger are gone from the standards prose: a number there was stale
  the commit after it was written, and the test fails at the moment
  of drift instead. What the walk certifies is template uniqueness (a
  reported message resolves to one line of source); the shared
  validation helpers interpolate the argument name at run time and
  say so.

## rOpenSci submission materials

* An audit of every `@srrstats` tag against the code it sits on
  corrected 30 claims that had drifted over 39 versions, including
  ten internal functions that had lost their roxygen and a
  not-applicable standard that was actually met and tested. The
  Bayesian and Monte Carlo category argument is rewritten around the
  package's documented inferential surface, with every citation
  verified against source.
* README gains the required statement of need and an installation
  section; `codemeta.json` is regenerated; `dev/submission-draft.md`
  fills the submission template, with the r-universe dependency
  precedent verified from brms, bayesplot and bridgesampling's own
  DESCRIPTION files.
* `?frm` documents what `data` accepts: a tibble, a data.table and a
  plain named list all work, a matrix column is a supported model
  variable, and a list column is not one.

# frmtmb 0.39.0

The LKJ correlation prior closes the last default-prior gap against
brms and brings correlated blocks into non-centered sampling;
population-level prediction becomes exact on models with random
smooths; and a mechanical lint pass lands with its refusal policy
written down.

## Breaking: population-level prediction with random smooths

* `predict(re.form = NA)` and `conditional_effects()` now give a true
  population-level prediction on models with random smooths. A smooth
  whose basis is indexed by a grouping factor (`bs = "fs"`, a factor
  `bs = "re"`, a tensor product with an `re` margin, classified from
  the smooth object rather than the `bs` string) is dropped along
  with the `(x | g)` blocks; population smooths, `gp()` and `hsgp()`
  keep their wiggle. The result matches
  `mgcv::predict.gam(exclude = )` on the factor-smooth term to 3e-9.
  Previously all smooth blocks were kept, so `re.form = NA` on such a
  model returned the per-subject curve. This deliberately follows
  mgcv rather than brms, whose `re_formula = NA` keeps factor-smooth
  curves; the divergence and the reason are documented in
  `?predict.frmtmb_fit`.
* `predict(newdata = , re.form = NA)` no longer needs the grouping
  column of a factor-smooth term at all, since the dropped term is
  never rebuilt. When the column is needed and absent, a named error
  lands before mgcv's internal one, and any other missing smooth
  column is named too.
* An unseen level of a factor-smooth term errors by default like any
  other new grouping level; `allow_new_levels = TRUE` predicts it at
  the population level for an `fs` basis (equal to `re.form = NA`
  exactly), and is refused by name for a factor `bs = "re"` basis,
  whose design has no zero row to hand a new level.
* `conditional_effects()` on a fit whose only predictors are matrix
  columns (scalar-on-function regression) names those columns instead
  of claiming nothing is plottable, and no longer offers a
  factor-smooth's grouping factor as an effect.

## The LKJ prior on random-effect correlations

* `prior_lkj(eta)` and `set_prior("lkj(eta)", class = "cor")` put an
  LKJ prior on a random-effect block's correlation. Class `"cor"`
  addresses a whole block, by `group` exactly as class `"sd"` does,
  which is brms's spelling. The density is the LKJ one carried onto
  frmtmb's unconstrained row-normalized Cholesky parameters with the
  exact Jacobian of that map, so on those parameters it means what
  `lkj(eta)` means on the correlation matrix. It covers `us()` and
  `gr(cov = )` blocks of two or more terms, and `cs()`, `ar1()` and
  `hetar1()`, whose single bounded correlation takes the LKJ marginal
  `(1 - rho^2)^(eta - 1)` with that structure's own Jacobian. `toep()`
  is refused by name: its parameterization is not positive definite
  everywhere, so it has no correlation matrix to put a density on.
* The formula route of `frm_sample()` now defaults to `lkj(1)` on every
  correlation, which is brms's own default and the last class where the
  two disagreed. **Seeded draws differ from v0.38 for any model with a
  correlated random-effect block.** `priors = "flat"` opts out, the fit
  route is unchanged (it is a likelihood diagnostic and keeps flat
  priors), `frm()` and `check_laplace()` are untouched, and the
  disclosure message gains a `cor` line.
* Correlated blocks are now sampled NON-CENTERED on the formula route.
  The gate did not move: a block is non-centered when every parameter
  it has carries a prior, and correlations now do. `(x | g)`, `cs()`,
  `ar1()`, `hetar1()` and `gr(cov = )` join `(1 | g)` and the smooths.
* Measured on `sleepstudy (Days | Subject)`, one chain of 2000, three
  seeds: the flat-correlation median of 36 min-ESS/s with a 154
  divergence storm on one seed becomes 120 (centered) and 116
  (non-centered) with no divergences, against ~122 for the matched brms
  model. The prior is what closed the gap; non-centering the block is a
  wash on this model. `dev/benchmarks.md` has the table and the
  distribution-level validation of the density.
* A MAP fit takes the same prior:
  `frm(..., priors = set_prior("lkj(2)", class = "cor"))` penalizes the
  correlation toward zero. Maximum likelihood is unaffected, because
  the default is a sampling default and never enters `frm()`.

## Housekeeping

* A mechanical lint pass: the three `on.exit()` calls gain
  `add = TRUE` (none was clobbering another handler), one
  `any(duplicated())` becomes `anyDuplicated()`, one literal pattern
  gains `fixed = TRUE`. `dev/lint-policy.md` records what was swept
  and what is deliberately not chased, with reasons, ahead of the
  rOpenSci submission.

# frmtmb 0.38.0

Non-centered sampling with the funnel diagnosed as a prior, and four
case studies: hidden Markov, functional regression, Wiener diffusion,
and circular.

## Non-centered sampling

* `frm_sample()` gains `reparameterize` (default `TRUE`): a
  random-effect block whose every parameter is a standard deviation
  carrying a prior is sampled non-centered, `b = L(theta) z`, and
  mapped back per draw. The draws matrix is unchanged in every
  respect (same `b[i]` columns, same names, same order), so
  `log_lik()`, `loo()`, `posterior_epred()`, `ranef()`,
  `conditional_effects()` and `hypothesis()` cannot tell the routes
  apart. Seeded formula-route draws differ from v0.37;
  `reparameterize = FALSE` reproduces the old chains exactly, and the
  fit route, and with it `check_laplace()`, is unchanged.
* Where each group's own data say little (the regime the funnel
  lives in) the gain is decisive: a binary GLMM with 80 groups of 2
  observations goes from a min bulk-ESS of 5 to 236. Where groups
  are informative it is a wash, and it is not offered as a blanket
  speed-up.
* A block is non-centered only when every parameter it has is a
  standard deviation carrying a prior, and `frm_sample()` names every
  block it left centered, with the reason. This is measured, not
  conservatism: under a flat prior the reparameterization hands the
  chain a flat tail the centered geometry was blocking.
* The correlated-slopes sampling deficit against brms is diagnosed,
  and it was not the centering: the flat prior frmtmb puts on a
  correlation's unbounded parameter is `(1 - rho^2)^(-3/2)`,
  improper, and a proper prior on that parameter fixes it in the
  CENTERED parameterization, beating the matched brms throughput.
  brms is ahead there because `lkj(1)` is proper, not because it
  non-centers. An LKJ default on correlations is the queued
  follow-up; `dev/benchmarks.md` carries the measurements.
* `check_laplace()` reports each parameter's bulk effective sample
  size and says so when the chain mixed too poorly to judge the
  approximation, instead of letting a sick chain read as a Laplace
  problem.

## Case studies

* `vignette("case-studies")` gains four sections: a two-state
  gaussian hidden Markov model on an animal track, validated against
  hmmTMB from generic starts (5.9e-10 in the log likelihood) and
  depmixS4 (9.7e-8); function-on-scalar regression through `s(t)`,
  `s(t, by = x)` and a `bs = "fs"` factor smooth, plus
  scalar-on-function regression through a matrix-column linear
  functional term, both validated against `mgcv::gam()` under ML;
  the Wiener first-passage density written as a `custom_family()`
  (RWiener agreement 8.7e-11 over 162 settings), with a data-bounded
  non-decision time through a link object; and von Mises circular
  regression with a cyclic-basis mean and distributional kappa,
  against the closed-form circular MLE to 6 decimals.
* Documented along the way: a smooth model's `logLik()` equals
  mgcv's ML smoothness-selection score `-gam$gcv.ubre` (not
  `logLik.gam`); `frmtmb_family(links = )` accepts a link object,
  which is how a bounded parameter keeps its support without a tape
  branch. RWiener joins Suggests, used only by the vignette and its
  density test.

# frmtmb 0.37.0

Leave-one-out cross-validation and the brmsfit method surface on
draws, the frm_ode() pharmacometrics tier, influence plots and naming
ergonomics, and the root cause of the ODE sampling defect, with
patches for upstream.

## Leave-one-out cross-validation

* `log_lik()` on `frm_sample()` draws gives the ndraws x nobs matrix
  of per-observation log-densities, each row at that draw's own
  parameter vector and conditional on its own group-level values,
  brms's convention (tmbstan samples the random effects too).
  Addition terms enter through the objective's own per-row density
  code, factored out for exactly this purpose, so `cens()`,
  `trunc()`, `weights()` and `trials()` cannot drift from what was
  fitted.
* `loo()` and `waic()` hand that matrix to `loo::loo.matrix()` and
  `loo::waic.matrix()`, with `r_eff` from `loo::relative_eff()` on
  the chain structure. `loo_compare()` ranks the criteria and is
  loo's own comparison when handed criterion objects. `psis()`
  returns the smoothed weights. Validated against brms on a matched
  model: elpd_loo within 4 percent of one standard error, log_lik
  column means correlated at 0.9999.
* `bayes_R2()` implements the residual-based estimator of Gelman,
  Goodrich, Gabry and Vehtari (2019) on expected-prediction draws,
  exactly as brms computes it.
* These are posterior quantities: `?loo` documents that
  `frm_sample(fit)`'s flat improper priors leave them unregularized
  and recommends the formula route or `priors =` for model
  comparison.
* A likelihood with no per-observation column refuses rather than
  inventing one: R-side autocorrelation, `hmm()`, `mixture(groups
  =)`, in-model imputation, and laplace-marginalized draws. Each
  names `AIC()` or `frm_bootstrap()` as the replacement.

## brmsfit method surface on draws

* Draws objects gain `as.array()`, `as.mcmc()`,
  `as_draws_array/df/list/matrix/rvars()`, `posterior_summary()`,
  `posterior_interval()`, `predictive_interval()`,
  `predictive_error()`, `ndraws()`, `nchains()`, `niterations()`,
  `nvariables()`, `nobs()`, `formula()`, `family()`, `getCall()`,
  `ngrps()`, `coef()`, `pp_mixture()`, and the bayesplot delegations
  `mcmc_plot()`, `pairs()`, `nuts_params()`, `log_posterior()`,
  `rhat()`, `neff_ratio()`. 65 of the 96 `brmsfit` methods now
  resolve.
* `conditional_effects()` on draws evaluates the effect grids once per
  posterior draw and bands the curves with their own quantiles: no
  `band =` or `method =` to choose, nonlinear predictors and nominal
  per-category displays work without a delta method, and it runs on
  formula-route draws that have no maximum-likelihood fit behind
  them. `ndraws` thins the draws for a cheaper curve.
* The default `effects` of `conditional_effects()` now include one
  `"a:b"` display per fitted interaction, as in brms; a term of order
  three or more contributes its leading pair, with the remaining
  variables at reference values until `conditions =` pins them.

## Living beside brms

* Attaching brms (or rstantools, nlme, lme4, loo, posterior,
  bridgesampling) after frmtmb no longer breaks dispatch: every frmtmb
  method on a generic those packages also export is additionally
  registered into the foreign generic's own table, so
  `conditional_effects()`, `hypothesis()`, `fixef()`, `log_lik()`,
  `posterior_epred()` and forty-odd others resolve regardless of
  attach order. Verified against the full surface with brms attached.
* `bf()` is a plain function and IS masked by `library(brms)`; a
  brms-built formula reaching `frm()` now says exactly that and names
  the two ways out (`frmtmb::bf()`, or attach brms first).

## Tape-safe scope for user code

* An ODE `dynamics` function and a `custom_family()` density no longer
  need the `"c" <- RTMB::ADoverload("c")` boilerplate: RTMB's
  tape-safe `c()`, `[<-` and `diag<-` are spliced into the function's
  own body (unless it already binds them), and nonlinear formula
  bodies evaluate with the same names in scope. A helper such a
  function CALLS still needs its own bindings, because lexical scope
  does not travel into other functions.
* The brms methods frmtmb does not have fail with the reason and the
  replacement instead of "could not find function":
  `loo_moment_match()`, `loo_subsample()`, `reloo()`, `kfold()`,
  `bridge_sampler()`, `bayes_factor()`, `post_prob()`, `stancode()`,
  `standata()`, `expose_functions()`, `restructure()`, and the
  deprecated spellings `posterior_samples()`, `nsamples()`,
  `parnames()`, `LOO()`, `WAIC()`.
* Suggests gains loo (the estimators are its own, so
  `loo::loo_compare()` and `loo::pareto_k_table()` read the results
  unchanged) and coda (`as.mcmc()` gives `gelman.diag()` one
  component per chain).

## Pharmacometrics in frm_ode()

* Time-varying covariates. `frm_ode(tv = )` takes dynamics inputs
  that change WITHIN a group, as a step function of time: the solve
  splits at each change point and each segment's dynamics see that
  segment's value as an extra parameter. The values may be estimated
  (a covariate-dependent clearance that changes at a known visit);
  the change points must be data, and `tv_break` names the column
  they come from. Last observation carried forward, rxode2's rule
  for covariates, and the state is continuous across a change
  because a covariate moves the derivative, not the state. A
  gradient through an estimated step rate matches central
  differences to 1e-11.
* Steady-state dosing. An `events` row marked `ss = TRUE` with an
  `ii` says its cycle has already settled: every compartment is
  zeroed, as NONMEM and rxode2 read such a record, and the cycle is
  repeated `n_ss` times (default 20) before the record's time. This
  is an approximation, not a closed form, and it says so: off the
  tape the last two cycles are compared and a shortfall past
  `ss_tol` warns. Exact against the analytic one-compartment
  superposition to 2.4e-10.
* Reset events. `method = "reset"` sets every state to `value`,
  which at zero is NONMEM's `EVID = 3` and rxode2's `evid = 3`. A
  reset beside a dose at one instant is `EVID = 4`; the reset goes
  first, so the table's row order does not decide the answer.
* `ii` and `addl` columns on `events` write a repeated schedule
  compactly, expanding to exactly the hand-written rows.
* `output` may be a column instead of a scalar: one state per ROW,
  over one shared solve. A parent and its metabolite in one assay
  column, or two species stacked long, is one call rather than two,
  at half the solves and the identical likelihood.
* Proportional-plus-additive error is documented in `?frm_ode` and
  `vignette("ode")` through `nlf(sigma ~ ...)` over the ODE mean. No
  new machinery, just the spelling.
* `parms` may be empty when every dynamics input is time-varying.

## Influence plots and naming ergonomics

* `plot()` on an `influence()` result draws the Cook's distances
  with the top cases labeled and one dfbetas panel per coefficient,
  each with the conventional 2/sqrt(n) reference band (Belsley, Kuh
  and Welsch 1980). `which` and `ask` follow `plot()` on a fit. Not
  `car::influencePlot()` compatibility, and the man page says why.
* `hypothesis()` says so, once per call, when a fixed-effect
  coefficient shadows a natural-scale name: a covariate literally
  named `sigma` hides the residual SD, and a coefficient spelling
  out `sd_<group>__<term>` or `ar1` hides its summary. The shadowed
  quantity is reachable under a leading dot (`.sigma`,
  `.sd_g__Intercept`), which `variables()` lists exactly when the
  collision exists.
* `vignette("inputs")` gains a naming-collisions section: which
  vocabularies can meet, which meetings are refused, and the actual
  resolution precedence in a nonlinear body (nonlinear parameter,
  then data column, then distributional parameter).
* `frm()` and `frm_sample()` bounds and `confint(parm =)` accept the
  bare name of an intercept-only nonlinear parameter (`la` for
  `la_(Intercept)`). A nonlinear parameter with several coefficients
  is refused, naming the full spellings.

## The ODE sampling defect, root-caused

* The v0.36 finding that `frm_sample()` on an `frm_ode()` fit fails
  at warmup iteration 1 is now fully explained, and it is not
  frmtmb's to fix: RTMBode's solver wrapper lets deSolve failures
  escape as R errors into Stan's C++ (a hard error when a rate
  reaches infinity, and an early return with fewer rows than
  requested). `dev/upstream/` carries the analysis, a three-patch
  series against the RTMB repository (verified: the Lotka-Volterra
  model samples end to end with 0 divergences on the patched
  build), reproduction scripts, and issue text ready to file. The
  same investigation found the more-than-8-states ceiling to be an
  integer overflow in deSolve's lsoda workspace formula.

## Fixed

* `diagnose()` no longer relays "NaNs produced" warnings from the
  pathological fits it exists to inspect; the NaN standard errors
  are the finding, reported through `bad_se`.
* The formula-vs-fit sampler agreement test compares the two taped
  densities pointwise (deterministic) instead of leaning on
  platform-dependent chains.

# frmtmb 0.36.0

Robustness in three senses (numerical, distributional, inferential),
per-parameter nonlinear formulas, one simulator per family, direct
sampling, and t2() newdata prediction.

## Behavior changes

* Draws objects use parenthesis-free parameter names throughout
  (`Intercept`, not `(Intercept)`), matching brms and what
  `variables()` and `hypothesis()` already spoke. This changes
  `summary(ds)`, `fixef(ds)`, `as_draws(ds)` and `check_laplace()`
  output; fits are unchanged, and both spellings still resolve in
  `hypothesis()` and `frm_sample(priors =)`.

## Nonlinear formulas per parameter

* `nlf()` declares a nonlinear formula for one parameter, as in
  brms: `bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1)`, and any
  parameter can carry one - `nlf(sigma ~ a + b * z)` is a nonlinear
  model for the residual SD with a linear mean, which `nl = TRUE`
  cannot spell. Nonlinearity is now a property of one parameter;
  bodies chain to any depth in dependency order, cycles are refused
  by name, and the old `nl = TRUE` models are bit-for-bit unchanged.
* A body may read another distributional parameter's per-row value:
  the spelling for nlme's `varPower(form = ~ fitted(.))`. A data
  column still wins over a parameter name, so ported brms bodies
  keep their meaning. `frm_ode()` works inside `nlf()` bodies.

## Simulation and sampling

* `simulate()`, `posterior_predict()` and `frm_simulate()` now share
  one simulator per family. Three previously silent wrongs are
  fixed: `mixture(groups =)` drew a class per ROW instead of per
  group in `posterior_predict()` and `frm_simulate()`, and residual
  autocorrelation draws were independent there. `hmm()` and
  `mixture_mvn()` gained simulators on all three paths; `cox()`
  explains its refusal at each entry point.
* `frm_sample()` accepts a formula and samples without a maximum-
  likelihood fit first: `frm_sample(bf(y ~ x + (1 | g)), data = dd,
  family = gaussian())` tapes the model and starts from random
  inits. The formula route defaults to brms 2.23's weakly
  informative priors (replicated against `brms::default_prior()`
  exactly, zero-shift rule included), announces every default and
  every deliberate gap (thresholds, shapes, correlations) in one
  message, and `priors = "flat"` opts out. `frm_sample()` on a
  FITTED model is unchanged - flat, likelihood-shaped - which is
  what keeps `check_laplace()` meaningful.
* `frm(dry_run = "objective")` returns the taped model unfitted.

## t2() smooths

* `predict(newdata = )` works for `t2()` smooths, including
  function-on-function terms and `t2(..., by = )`. A t2 basis
  carries a second, prediction-only constraint that the fit does
  not; smooths are now built with `smoothCon(modCon = 3)`, which
  drops it and leaves the design, penalties and log-likelihood
  bit-identical. The refusal is gone.

## Student-t random effects

* `(x | gr(g, dist = "student"))` gives a grouping term a Student-t
  latent instead of a gaussian one, brms's spelling. An outlying group
  then costs the variance component far less: with one group of 30
  displaced by 10 standard deviations, the gaussian latent SD inflates
  to 1.85 against a truth of 1 while a `t(5)` latent holds at 1.25, and
  the intercept's RMSE falls from 0.298 to 0.180. Correlated and
  `diag()` blocks are one multivariate t with a single mixing variable
  per level, which is brms's construction and is verified against
  `mvtnorm::dmvt()` to 1e-10.
* The degrees of freedom are FIXED, at `dist_nu` (default 5), not
  estimated. brms estimates `nu` under a `gamma(2, 0.1)` prior, and the
  prior is carrying it: by maximum likelihood the whole grid from 2.1 to
  500 sits inside the 95% profile interval at 20 groups, and joint ML
  runs to a boundary in 24% to 41% of replicates there. This is the
  frequentist analogue of brms's own
  `prior(constant(3), class = "df")`. Compare a few values with
  `logLik()` instead.
* The reported quantity is the t's SCALE, not its standard deviation -
  which is also what brms's `sd_<group>__<term>` is. `VarCorr()` stores
  the scale matrix, tags it with `nu`, and prints a `Scale` column, a
  converted `Std.Dev.` column (`scale * sqrt(nu/(nu-2))`) and the fixed
  `nu`, so the convention is visible instead of silent.
* `quadrature = TRUE` accepts scalar t blocks and marginalizes them
  EXACTLY, agreeing with adaptive Gauss-Hermite quadrature to 1e-6. It
  is the recommended check: the Laplace approximation over a t latent
  biases the estimated scale UPWARD, by under 2% of a standard error at
  8 observations per group but by a factor of three where a near-null
  variance component meets two-observation groups.
  `dev/tre-feasibility.md` has the measurements.
* `simulate()` and `frm_simulate()` draw a multivariate t;
  `predict(allow_new_levels = TRUE)` inflates the unseen level's
  variance by `nu/(nu-2)`; `REML = TRUE`, `ranef()`, `confint()` and
  `frm_bootstrap()` work unchanged. Refused with their own messages, and
  documented in `?frmtmb-student-re` and `frm_compat()`:
  `gr(cov = )`/`gr(prec = )`, `mm()`, `|ID|` keys, and every covariance
  structure but `us` and `diag`.

## Cluster-robust (sandwich) standard errors

* `vcov_cluster(fit, cluster = ~ g, type = "CR0")` gives the
  cluster-robust covariance of the estimates: the inverse observed
  information of the marginal likelihood sandwiching the outer product
  of the per-cluster scores. `type` takes `"CR0"`, `"CR1"`, `"CR1p"`
  and `"CR1S"`, spelled and defined as in `clubSandwich`. `vcov(fit,
  cluster = ~ g, type = )` is the same thing in the
  `sandwich::vcovCL()` spelling.
* `cluster_scores(fit, ~ g)` returns the per-cluster scores
  themselves, one row per cluster, so any other sandwich can be built
  from them. This is what `sandwich::estfun()` would return if a
  marginalized objective had per-observation contributions; it does
  not, which is why frmtmb still ships no `estfun()` method.
* `confint()`, `hypothesis()` and `summary()` take `vcov = `: a
  covariance over the whole outer parameter vector, or a function of
  the fit returning one. A matrix from `vcov_cluster()` carries `G - 1`
  reference degrees of freedom, so those methods switch from a normal
  to a `t` reference automatically.
* The estimator is refused, with the reason, wherever the marginal
  likelihood does not factor over the clustering factor: a random
  effect whose level spans two clusters (crossed effects, `mm()`
  pooled levels, a global smooth, `gp()`, `car()`, `spde()`), a
  group-level mixture whose groups span clusters, `autocor()`,
  `hmm()`, `rescor = TRUE`, `mi()`/`me()`, `REML = TRUE`,
  `frmtmb_control(profile = TRUE)`, `quadrature = TRUE`, and any fit
  made with priors. `frm_bootstrap()` is the documented fallback.
* `"CR2"`/`"CR3"` are refused rather than approximated: the
  Bell-McCaffrey family is defined through the hat matrix of a linear
  or GLS model, which a Laplace-marginal likelihood with a nonlinear
  link does not have.

## Numerical robustness of the family log-densities

The log-densities now evaluate from the LINEAR-PREDICTOR scale wherever
an exact form exists, instead of undoing the link and recomputing. An
inverse link saturates - `plogis(40)` is exactly 1 in double precision,
and `1 - exp(-exp(4))` is too - so a density written over `1 - mu`
returned `-Inf` with an unusable gradient in a region the linear
predictor describes perfectly well. A fit that converges never visits
it; an optimizer step that overshoots, a separated predictor, a wide
quadrature node and a `frm_sample()` tail draw all do.

* `bernoulli()`, `binomial()` and `zero_inflated_binomial()` use
  `RTMB::dbinom_robust()`, the log-odds parameterization glmmTMB fits
  with. They were wrong in the second decimal at a linear predictor of
  30 and gave `-Inf` at 40; they are now exact over the whole line.
  With the `cloglog` link the old form had no gradient past a
  SINGLE-DIGIT linear predictor.
* `negbinomial()`, `nbinom1()` and `geometric()` use
  `RTMB::dnbinom_robust()`, which takes `log(mu)` and never forms
  `mu / var`. The old form gave `NaN` at a linear predictor of -40.
* Every `zi` and `hu` gate, and `asym_laplace()`'s `quantile`, run in
  log space, so a separated zero-inflation predictor stays
  differentiable.
* `Beta()`, `zero_inflated_beta()` and `beta_binomial()` take both
  shapes off the log-odds, so the second one no longer collapses to
  zero.
* `hurdle_poisson()`'s zero-truncation normalizer `log(1 - exp(-mu))`
  uses `expm1`, which keeps it at a tiny mean.
* `cumulative()`, `sratio()` and `cratio()` with the logit link compute
  their CDF differences in log space. `acat()`, `categorical()` and
  `multinomial()` accumulate their denominators with `logspace_add()`
  instead of summing `exp()`, which used to overflow at a linear
  predictor of 709 / (K - 1).
* Values and gradients at ordinary linear predictors are unchanged to
  1e-15, and every fitted model in the regression set lands on the same
  estimates (largest coefficient change 9.5e-9, largest log-likelihood
  change 4e-12). `cumulative(probit)` keeps the plain CDF difference,
  which has no exact log-space form; it is fragile past a linear
  predictor of 8, so prefer the logit link. New file
  `tests/testthat/test-numerical-robustness.R` pins the behavior.

## New family

* `huber()` fits Huber's least-favorable distribution: gaussian within
  `k` residual standard deviations of `mu` and Laplace outside, so an
  outlier pulls on the fit with bounded influence. It is a proper
  normalized density, so this is ordinary maximum likelihood and
  `logLik()` and `AIC()` mean what they say. Huber's tuning constant is
  a fixed argument, `huber(k = 1.345)`, not an estimated parameter,
  matching how `MASS::rlm()` treats it. Point estimates track
  `MASS::rlm(psi = psi.huber)`; the remaining gap is the scale
  (`rlm()` fixes it at a MAD-type estimate, `huber()` estimates `sigma`
  by ML), and holding `sigma` there reproduces `rlm()` to 1e-5. The
  `asym_laplace()` working-likelihood caveat applies: Wald standard
  errors under a misspecified error distribution are not calibrated, so
  use `frm_bootstrap()`.

# frmtmb 0.35.0

Hidden Markov models, latent class analysis, multi-membership random
effects, three new families (nominal, circular, proportional
hazards), profile and bootstrap effect bands, and a brms-portability
batch measured against the brms vignettes themselves.

## Behavior changes

* `frm()` defaults to `family = gaussian()` when neither the `family`
  argument nor a `+` attachment supplies one, the brms/lme4/glmmTMB
  convention. `frm(y ~ x, data = d)` was an error and is now a linear
  model. The `family` argument overrides a family attached to a
  univariate `bf()` and fills only the empty responses of an
  `mvbf()`.
* Documentation and examples now lead with the separate-family
  spelling `frm(bf(y ~ x), family = gaussian(), data = d)`; the `+`
  attachment stays valid and documented as the alternative.

## Portability (measured against the brms vignettes)

* A scorecard audit ports every model call of the brms vignettes
  through the mechanical `brm` to `frm` transform
  (dev/brms-vignette-port.md): with this release's fixes the
  measured tally is 30 of 42 model calls running mechanically, the
  remaining spelling changes documented in the porting guide.
* `hypothesis()` accepts brms's directional form (`"a > b"`) with
  one-sided p-values and bounds across every method, plus the
  `class=`/`group=` naming shorthand. `update()` speaks
  `formula.`/`newdata` and dotted deltas (`. ~ . + z`). `lf()`
  composes dpar formulas onto `bf()`. One formula may name several
  parameters (`a + b ~ 1`). A bare constructor (`family =
  cumulative`) works.
* One parameter-addressing vocabulary across `hypothesis()`,
  `profile()`, `confint(parm =)` and bounds: parenthesized and
  paren-stripped spellings are interchangeable, and one-to-one
  natural-scale names (`sd_g__x`, `ar1`) reach their internal
  parameter with the scale stated. Ambiguous natural-scale names are
  refused with the alternatives named.

## Multi-membership and effect bands

* Multi-membership random effects: `(x | mm(g1, g2, weights =,
  scale =))` and `mmc()` member-specific covariates, following brms;
  the pooled-level Z construction matches `brms::make_standata()`
  exactly, and every post-fit method reads the block unchanged.
  New-level prediction variance distinguishes a shared new level
  (one draw, weights add) from distinct new levels (independent
  draws).
* `conditional_effects()` gains `band = c("wald", "profile",
  "boot")`: likelihood-ratio bands inverted per grid point, and
  pointwise percentile bands from one shared parametric bootstrap
  (reused across effects and verified against the grid it was run
  over). Effect discovery now works on nonlinear and `mo()`/`mi()`
  fits, and `frm_bootstrap()` works on ordinal fits (previously all
  NA).

## New families

* `categorical()` fits a multinomial logit to an unordered factor,
  brms's spelling of the likelihood `multinomial(K)` already carried
  on a count-matrix response. The first level is the reference and
  each remaining level gets its own linear predictor named
  `mu<Level>`, as in brms; the main formula applies to all of them
  unless a dpar formula overrides one, so `bf(y ~ x, mustout ~ w)`
  gives one category its own predictor. `fitted()` and
  `predict(type = "response")` return the `n x K` matrix of category
  probabilities with the response's own levels as column names, the
  same convention the ordinal families follow, and `simulate()` draws
  factor levels. Validated against `nnet::multinom` (4e-10 in the
  log-likelihood) and against `multinomial(K)` on the one-hot
  response (1e-8). The categories are read off the data by `frm()`
  before the formula is parsed; `categorical(levels =)` or
  `categorical(K =)` states them for the paths that have no data.
  `residuals()` is refused: a nominal response has no scale for one.

* `von_mises(link = "tan_half")` for a circular response in
  `(-pi, pi]`, with the mean direction `mu` through the tan-half link
  and the concentration `kappa` through a log link - brms's
  parameterization. The normalizing constant's `log I0(kappa)` is
  differentiated exactly by RTMB's own `besselI` method, so nothing
  is a series approximation. `kappa ~ x` distributional models work,
  and the simulator is a vectorized Best-Fisher sampler that varies
  both parameters per row. Validated against a hand-rolled RTMB
  likelihood (1e-6) and `circular::mle.vonmises()`.

* `cox()` fits proportional hazards with brms's flexible baseline: an
  M-spline baseline hazard over a simplex of weights, with the
  cumulative baseline hazard the I-spline integral of the same basis
  (`df = 5` cubic splines with an intercept, brms's `bhaz()` default).
  Censoring runs through the ordinary `cens()` addition term - an
  event contributes the density and a censored row the survivor
  function - so right, left, and interval censoring all work, and
  frailty models come free through the Laplace approximation:
  `time | cens(c) ~ x + (1 | g)`. New `cox_baseline()` returns the
  fitted baseline weights. Validated exactly against a hand-rolled
  M-spline PH likelihood (1e-6 in the log-likelihood, 1e-4 in the
  coefficients) and against `survival::coxph()` to 2e-2 on the log
  hazard ratio, which is the honest claim: `coxph()` leaves the
  baseline fully nonparametric. `fitted()`, `predict(type =
  "response")` and `simulate()` are refused - a survival time has no
  mean the censored rows identify, and brms refuses the same
  question. Unpenalized maximum likelihood often puts a baseline
  weight on the simplex boundary, which the optimizer reports as
  singular convergence even at the optimum; `?frmtmb-families`
  explains it and lowering `df` is the remedy.

* New `tan_half` link, for `von_mises()`.

## Internal

* A family may declare `aterm_data(y, aterms)`, family-level data that
  no addition term supplies, built once at frame assembly from the
  validated response. `cox()`'s spline bases use it.

* A family's `lcdf()` may take the family-level extra parameters as a
  fourth argument, which is what lets a survival family's survivor
  function reach its baseline. The three-argument contract is
  unchanged for every other family.

Hidden Markov models as a first-class family.

## New

* `hmm(K, family, time =, group =, init =, trans =)` fits a `K`-state
  hidden Markov model. The response at each time point comes from one
  of `K` state-dependent copies of `family`, and the unobserved state
  follows a first-order Markov chain along `time` within `group`. The
  state sequence is summed out EXACTLY by the forward algorithm,
  evaluated on the same AD tape as everything else, so nothing about
  the Laplace machinery changes: a random effect in a state's linear
  predictor is integrated outside the exact state sum. Each of the
  wrapped family's parameters is copied per state and suffixed
  (`mu1`, `sigma1`, `mu2`, ...) with the full formula grammar,
  random effects included, exactly as `mixture()` does.
  Transition probabilities are a row-wise multinomial logit named
  `tr{i}{j}` with state 1 the reference cell of every row, and
  `trans = ~x` gives every cell a predictor at once. Validated
  against `depmixS4` (logLik to 1e-8 or better and every parameter to
  five decimals, on gaussian, poisson and categorical emissions with
  and without transition covariates) and against `hmmTMB` (1e-12 on
  the stationary fixed-effect model).
* `init = "stationary"` (the default) solves the chain's stationary
  distribution on the tape and costs no parameters; `"estimated"`
  adds `K - 1` free logits, `"uniform"` fixes them. Stationary is
  refused when a transition cell carries a predictor, because a
  time-varying chain has no single stationary distribution.
* `hmm_probs()` returns the smoothed state occupancies
  `P(S_t = k | y)` from a forward-backward pass - the `mixture_probs()`
  analog - and `hmm_viterbi()` the maximum-a-posteriori state path.
  `fitted()`, `predict(type = "response")` and the response and
  pearson residuals all route through `hmm_probs()`, so they report
  the occupancy-weighted mean `sum_k P(S_t = k | y) mu_k(x_t)` rather
  than any single state's. `simulate()` walks the chain forward per
  sequence and then emits, so DHARMa and `pp_check()` see the fitted
  persistence.
* A missing response is a time point the chain passes through without
  emitting: `hmm()` keeps the row and masks its emission instead of
  letting `na.action` drop it, which would shorten the chain and bias
  the transition matrix. `nobs()` counts every row; `fitted()` and
  `residuals()` are `NA` there while `hmm_probs()` is not.

## Refusals

* `REML`, `quadrature = TRUE` and `frmtmb_control(profile = TRUE)` are
  refused on an `hmm()` fit, each with the reason. REML in particular
  used to run and produce a partial restricted likelihood matching no
  standard definition.
* `weights()`, `cens()`, `trunc()`, `se()` and `mi()` on the response,
  multivariate models and `rescor`, `residuals(type = "osa")`,
  `residuals(type = "deviance")`, `predict(se.fit = TRUE)` and
  `predict(newdata =)` on the response scale, and
  `conditional_effects()` are refused, all naming the per-sequence
  likelihood as the reason.
* A grouping in which every sequence has length 1 is refused: the
  transition parameters are then flat directions the reported `df`
  would still count, and the model is a `mixture()`. Holding every
  transition dpar at a constant lifts the refusal, and the fit then
  reproduces `mixture()` to 1e-12.
* A start with every state's location predictor at the same value
  warns: it sits on the label-symmetry axis, where the optimizer
  cannot separate the states. The default quantile-spread starts never
  do this.

## Known limits

* Multimodality is real and no convergence diagnostic reports it. On a
  random-effect model the default cold start has been measured
  converging 8.1 log-likelihood units below the global optimum with
  `convergence == 0`, `max|grad| == 3.5e-4` and a positive-definite
  Hessian. Compare starts before reporting.
* With random effects the Laplace approximation is genuinely
  approximate even for a gaussian response, because the integrand is a
  mixture over state sequences: the measured bias is -0.126 in the
  log-likelihood (8.9e-5 relative) and 4.4e-4 in the parameters
  against adaptive quadrature. `?hmm` says so.
* The tape build grows slightly faster than linearly in the number of
  rows: about 1.9 s at 20 000 rows, against milliseconds to evaluate
  it. Below 5 000 rows nothing is noticeable.

## New

* Latent class analysis, the poLCA measurement model, as the family
  `lca(K)`. The response is a matrix of polytomous item codes, one row
  per subject and one column per item; the items are conditionally
  independent given a subject's latent class, and each class carries
  its own item-response profile. Class membership is the `theta1 ...
  theta{K-1}` dpars with full linear predictors, so a covariate on the
  model formula gives poLCA's latent class regression for free, with
  `fixef()`, `confint()`, `hypothesis()`, `set_prior()` and
  `frm_sample()` on the gating coefficients. Items may have different
  numbers of categories; `ncat` declares them, and by default they are
  inferred as the largest observed code per item, as poLCA does. The
  item profiles are family extra parameters, one vector `pi<j>` per
  item holding its `K * (C_j - 1)` reference-category logits, so they
  appear per item in `summary()` and as `pi<j>_<i>` in `confint()`.
  Validated against poLCA on its own shipped data: the carcinoma
  3-class model agrees to 4.2e-8 in log-likelihood, 2.8e-8 in item
  profiles, 6.2e-9 in class sizes and 9.9e-8 in posterior membership,
  and the election latent class regression agrees to 1.1e-7 in
  log-likelihood, 1.1e-7 in item profiles and 9.9e-7 in gating
  coefficients. A hand-rolled
  `optim()` reference on simulated data agrees to 1.5e-9, and a
  one-item fit reaches the saturated single-categorical likelihood to
  5.6e-11.
* `lca_profiles()` returns the class-conditional item-response
  probability tables (poLCA's `probs`) with the estimated class sizes
  attached, and prints them. `lca_probs()` returns posterior class
  membership per subject (poLCA's `posterior`) with the relative
  entropy of the classification attached; it is `mixture_probs()`
  under an LCA-specific name and check.
* `lca(na.rm = FALSE)` keeps subjects with missing items and masks
  each missing item out of that subject's likelihood, poLCA's
  `na.rm = FALSE` behavior. The default drops incomplete subjects
  through the usual `na.action`, which is poLCA's default.

## Notes

* `lca()` starting values are deterministic: subjects are scored by
  the mean of their item codes rescaled to `[0, 1]`, cut into `K`
  equal-count slices, and each slice's smoothed category proportions
  seed one class. Class 1 is the low-score end, so a data set always
  gets the same labeling. Multimodality is unchanged; `?lca` shows the
  perturbed-`start` loop that replaces poLCA's `nrep`.
* `lca()` v1 refuses random effects, smooths and `gp()` anywhere in
  the model (that is the growth-mixture shape, which
  `mixture(..., groups = ~g)` fits), `REML`, `profile = TRUE`,
  `quadrature`, every addition term, `mvbf()`, and
  `residuals(type = "osa")`. `fitted()`, `residuals()` and
  `predict(type = "response")` are refused because a matrix of nominal
  item codes has no mean; `predict()` returns the gating linear
  predictor.

# frmtmb 0.34.0

Within-group residual correlation (R-side effects), quantile
regression completions, and plot conveniences.

## New

* Within-group residual correlation, brms's R-side autocorrelation,
  as formula terms `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()`
  for `gaussian()` and `student()`. The residuals of a group become
  one multivariate draw, `y_g ~ N(mu_g, D R D)`, with a `sigma ~ x`
  distributional model entering through the diagonal. The ARMA
  autocorrelation function is exact and orders above one are
  supported (brms caps its covariance form at one). Validated
  against `nlme::gls` and `nlme::lme` under ML and REML to 1e-9 or
  better across all five structures, on balanced and ragged data,
  and with random effects alongside the correlated residual. New
  `autocor_matrix()` returns the fitted correlation matrix; the
  parameters appear in `summary()`, `confint()`,
  `confint_varcorr()` and `hypothesis()` under brms's names. See
  `?frmtmb-autocor` for the deliberate divergences: `sigma` is the
  MARGINAL residual SD (brms uses the innovation SD for ar/ma/arma;
  the migration vignette gives the conversion), lags count gaps in
  the global time-level set (nlme semantics), and brms's default
  `cov = FALSE` residual-regression form is refused rather than
  reinterpreted. `simulate()` draws correlated residuals, so
  DHARMa and `pp_check()` see the fitted autocorrelation.
  Combinations that stop the likelihood factorizing over rows
  (weights, cens/trunc, se, mi, rescor, mixtures, quadrature,
  nl, OSA residuals, other families) are refused with the
  alternative named.
* `zero_inflated_asym_laplace()` (brms spelling): a point mass at
  zero mixed with the asymmetric Laplace, for zero-inflated
  quantile regression. A new "Quantile regression inference"
  section on the families page states plainly that ALD-based Wald
  intervals are not calibrated (a property shared with brms) and
  points at `frm_bootstrap()`.
* `plot(conditional_effects(fit), points = TRUE)` overlays the raw
  observations, the brms argument, previously ignored. Points are
  drawn under the bands; displays where raw observations are not
  meaningful (per-category ordinal, non-mean dpars, matrix
  responses) say so instead of failing silently.
* Hidden Markov models: a feasibility study (dev/hmm-feasibility.md)
  establishes that the forward algorithm tapes in RTMB, composes
  with the Laplace approximation over random effects, and is
  expressible today through `custom_family()` with `vint()`
  payloads, validated against depmixS4 and hmmTMB. No user-facing
  grammar yet; the design for a first-class `hmm()` family is
  recorded.

# frmtmb 0.33.0

Two residue fixes, a documentation dark mode, and CI repairs.

## Behavior changes

* `posterior_epred()` on an ordinal fit returns a `draws x
  observations x categories` array with `dimnames` `list(NULL,
  observation names, category levels)`, replacing the flattened
  `draws x (n * K)` matrix of v0.31/v0.32. This is brms's documented
  shape (an S x N x C array for categorical and ordinal models). The
  array is the old matrix reshaped: `matrix(as.vector(ep), nrow =
  ndraws)` recovers the previous value and column order exactly, and
  `ep[, , "high"]` or `ep[k, , ]` replace the old naming recipe.
  Scalar-response families, `posterior_predict()`, and
  `posterior_linpred()` are unchanged. A documentation error is also
  corrected: `multinomial()` response-scale predictions are a
  vector, not a category matrix, and the docs no longer claim
  otherwise.

## New

* A `bf(nl = TRUE)` body can name an object of its formula
  environment that could never be a column of `data` (a data.frame,
  a list, an environment, a formula), so `frm_ode(..., events =
  doses)` takes a dosing table by name instead of needing an inline
  `data.frame()` or a wrapper function. The boundary is exactly what
  `stats::model.frame()` refuses: vectors, factors, and matrices are
  legal model-frame variables and still resolve through the model
  frame, and a column of `data` still wins over a same-named object.
  A body that fails to evaluate reports which names were resolved
  outside the data.
* The documentation site has a light/dark/auto theme toggle.

## CI

* The pkgcheck and dependency workflows resolve RTMBode from
  r-universe: `extra-repositories` on the setup-r action, and an
  `R_PROFILE_USER` profile for the pkgcheck container, whose pak
  cannot read `Additional_repositories` from DESCRIPTION.

# frmtmb 0.32.0

The merged-|ID| Kronecker path, the ordinal prediction surface, ODE
dosing, and vignette figures.

## New

* `frm_ode()` gains `events`, a per-dataset dosing table (`time`,
  `value`, `state`, optional `group`, `method` of add/replace/
  multiply, `duration` for constant-rate infusions) supporting
  repeated doses, and `event_scale`, an estimated multiplier on
  every dose amount, so bioavailability is an ordinary nonlinear
  parameter with covariates and random effects. Dosing splits the
  solve at the event times instead of using deSolve's `events`
  argument: deSolve events do not differentiate correctly through
  RTMBode's augmented system (`replace` and `multiply` measured 42
  and 59 percent relative gradient error; reported upstream). A
  branch on time inside a derivative function is refused by RTMB,
  and an `approxfun()` forcing table inside one is silently wrong;
  both are documented. `vignette("ode")` gains a repeated-dosing
  section.
* Terms sharing an `|ID|` key whose grouping is `gr(g, cov = A)` or
  `gr(g, prec = Q)` now merge into one Kronecker block instead of
  being refused: `mvbf(bf(y1 ~ (1 | q | gr(id, cov = A))), bf(y2 ~
  (1 | q | gr(id, cov = A))))` is the same joint density as the
  long format `(0 + trait | gr(id, cov = A))`, verified to 7e-10 on
  the log likelihood under ML and REML and on both the covariance
  and precision sides. The compatibility rows move from refused to
  conditional.
* Vignettes carry figures where a picture beats a table: a forest
  plot for the meta-analysis, the monotonic step shape, the
  location-scale band with the mgcv overlay, growth-mixture
  trajectories by recovered class, the sleepstudy small multiple,
  and concentration-time curves in the ODE vignette. Figures use
  the suggested tinyplot package and are skipped when it is not
  installed.

## Behavior changes

* `fitted()` on an ordinal family returns the n x K matrix of
  category probabilities (the brms convention), exactly as
  `predict(type = "response")` does; it previously returned the
  latent linear predictor. The invariant `predict(type =
  "response") == fitted()` now holds for every family with no
  documented exception. `predict(type = "link")` keeps the latent
  predictor.
* `residuals(type = "response")` and `type = "pearson"` on an
  ordinal fit score the categories by their codes: `y - sum_k k *
  P(y = k)`, standardized by that distribution's own sd.
  `"response"` previously returned the observed code minus the
  latent predictor without saying so, and `"pearson"` errored.
* One `|ID|` label spread over more than one grouping specification
  is now an error. Previously the merge key included the deparsed
  grouping call, so `(1 | q | g)` next to `(1 | q | gr(g, cov =
  A)))` landed in separate blocks and silently did not correlate at
  all, which is not what the shared label asks. `|ID|`-linked
  `gr()` terms whose matrices resolve to different objects are
  refused at frame assembly, comparing the resolved matrices.

## Corrected behavior

* `emmeans()` on an ordinal fit no longer fails with
  "Non-conformable elements in reference grid"; the basis is built
  by column name (marginal means stay on the latent scale, the
  emmeans convention for `clm`-like models). `pp_check()` and
  `dharma_residuals()` work on ordinal fits (simulated ordered
  factors are compared on their integer codes; DHARMa's rank
  transform runs on simulated categories with `integerResponse =
  TRUE`). `plot()` works on an ordinal fit.
* `conditional_effects()` on an ordinal fit draws one probability
  curve per response category (a `cats__` column), with
  delta-method standard errors over the joint covariance of
  coefficients, thresholds, and `cs()` terms, and bands on the
  logit scale so they cannot leave [0, 1]. `method = "predict"` is
  refused there; `dpar = "mu"` gives the latent display.

# frmtmb 0.31.0

ODE models, an exotic-models case-study vignette, and the defect wave
the vignette work exposed.

## New

* `frm_ode()` fits ordinary differential equation models inside a
  `bf(nl = TRUE)` body. Dynamics parameters and initial states are
  ordinary nonlinear parameters, so they take fixed effects, random
  effects and covariates, and the Laplace approximation is exact
  through the solver's adjoint. It solves one small system per group
  and scatters the solution back into row order; ragged designs,
  unsorted rows, repeated times and an observation at `t0` all work.
  Needs the optional RTMBode package (not on CRAN;
  `Additional_repositories` now names https://kaskr.r-universe.dev),
  and every code path, test, example and vignette chunk degrades
  cleanly without it. New `vignette("ode")` works a population
  pharmacokinetic model on `datasets::Theoph`.
* `frm_ode()` refuses the ways an ODE model goes quietly wrong: a
  dynamics parameter that varies inside a solve group is rejected by
  name at frame assembly; fixed-step integrators such as `rk4` are
  refused because they return a different likelihood, not a noisier
  one; a system near the Laplace ceiling of about eight states warns.
  A failed solve is reported, not absorbed: on numeric evaluation
  paths (`predict()`, `simulate()`, `residuals()`, and a body with no
  estimated dynamics inputs) the group's rows get a penalty value and
  a warning naming the group, and `frm_ode_failures()` reads the
  record back. During fitting most failures cannot be detected at
  all (RTMBode returns an AD object), so they surface as the
  optimizer's NA/NaN gradient; the help page says which is which.
* New `vignette("case-studies")`: eight worked models from the
  showcase literature of brms, MCMCglmm, metafor and mgcv. The
  animal model and its multi-trait form, phylogenetic regression,
  random-effects meta-analysis and meta-regression, monotonic
  ordinal predictors, location-scale regression, growth mixtures,
  measurement error, and sequential ordinal models with
  category-specific effects. Every section cross-checks its fit
  against a reference package or a closed form, the multi-trait
  section proves its pedigree matters with an identity refit, and
  `tests/testthat/test-case-studies.R` pins the agreements.
  Suggests gains ape and metafor for the cross-checks, and
  data.table and tibble for tabular-input tests.

## Breaking

* A `gr(cov = )` or `gr(prec = )` term whose `|ID|` key is shared
  with another term is refused at parse time. The cross-formula
  merge hardcoded an unstructured covariance and silently discarded
  the matrix, so earlier fits of that construct ignored the
  structure they named. The error points at the supported long
  format, `(0 + trait | gr(id, cov = A))`, which takes the verified
  Kronecker path. A lone `gr()` term with an unshared key keeps its
  structure and still fits.

## Corrected behavior

* `hypothesis()` and `variables()` expose `sd_<group>__<term>` and
  `cor_` names for `gr(cov = )`, `gr(prec = )` and `equalto()`
  blocks, so heritability-as-ICC is one line on an animal model.
  Smooth, `gp()`, `car()` and `spde()` blocks stay excluded;
  `?hypothesis` says which and why.
* `vcov(full = TRUE)` works under REML and `profile = TRUE`, taking
  the outer block from the joint precision instead of warning and
  returning only the fixed effects. Its rows are `confint()`'s rows.
* The observation-level-random-effect check no longer fires on
  `gr(cov = )`, `gr(prec = )` or `equalto()` blocks with one row per
  level: a fixed relationship matrix identifies the two variances,
  which is exactly the animal model with one record per individual.
* `predict(type = "response")` on the four ordinal families returns
  an n x K matrix of category probabilities named by the response's
  levels, honoring `cs()` on newdata; `se.fit` there is refused
  rather than faked. marginaleffects gets the categorical `group`
  convention (also fixing `multinomial()` row alignment),
  `posterior_linpred(transform = TRUE)` stays on the latent dpar as
  documented, and `posterior_epred()`'s flattened categorical shape
  is named and documented.
* `anova()` and `drop1()` report NA instead of `< 2.2e-16 ***` when
  the compared models have equal df.
* `confint_varcorr()` no longer reports a zero-width interval at its
  internal clamp for a correlation estimated at plus or minus 1 or
  an sd at zero; those rows keep the estimate, get NA bounds, and
  warn. `VarCorr()` and `ranef()` handle two blocks sharing a term
  label (one printed twice, the other silently dropped a block).
* A `Date`, `POSIXct` or `difftime` predictor is reported at frame
  assembly with its unit and 1970-01-01 origin plus the centering
  fix; the epoch-scale magnitude can stop the optimizer converging.
  Responses and grouping variables are exempt (srr G2.5, G2.9).
* A nonlinear formula body that fails to evaluate reports which of
  its names were resolved to functions instead of found in `data`,
  so a misspelled column sharing a name with a base function (`t`,
  `c`) is named rather than surfacing as a coercion error.
* srr completeness: `srr_stats_pre_submit()` reports no missing
  standards (seven were untagged; six were already satisfied and are
  now tagged with pinning tests, G2.9 is satisfied by the datetime
  reporting above).

# frmtmb 0.30.0

data2, ODE feasibility, and the remaining rOpenSci runway.

## New

* `frm(..., data2 = list())`, the brms spelling for structural
  objects: the matrices of `gr(cov =)`/`gr(prec =)`, `equalto()`,
  `car()` and `spde()` resolve from `data2` before the data and the
  formula environment (compound expressions over `data2` objects
  work, a documented permissive divergence from brms). `data2` is
  stored on the fit, so `saveRDS()` carries the matrices and
  `refit()`/`update()`/`drop1()`/`influence()`/`frm_multiple()`
  re-assemble in a session where the calling environment is gone.
* Population pharmacokinetic models work today through `nl = TRUE`
  bodies calling `RTMBode::ode()` per subject (RTMBode installs
  from kaskr.r-universe.dev): the Laplace approximation is exact
  through the adjoint solver, and a one-compartment population fit
  agrees with nlmixr2's FOCEi to three decimals with no compiler.
  The feasibility study, the sharp edges (never stack subjects
  into one system; adaptive integrators only), and the planned
  `frm_ode()` helper are recorded in dev/ode-feasibility.md.

## rOpenSci runway

* Every exported help page has executing examples (35 pages added;
  the as-cran check runs them all in ~14s).
* Every internal function (288) is documented with its contract
  under `@noRd`; a new "Inputs and preprocessing" vignette covers
  the formula-to-design pipeline, accepted predictor classes,
  argument types, attribute survival, terminology, and reproducible
  runtime claims. srr standards: 105 of 116 tagged in code, docs,
  and tests; the other 11 documented as not applicable.
* Nine standards tests added (degenerate inputs, noise invariance,
  exact-fit behavior, row-name retention, runtime scaling).
* CI gained the rOpenSci `pkgcheck` action (the editors' submission
  checks, srr included) and a coverage workflow.

## Corrected behavior

* Rows dropped by `na.action` are reported with a count (a
  `message()`, so `suppressMessages()` silences it); fitting zero
  rows says so instead of blaming missing values.
* Eleven duplicated error messages across the package were
  disambiguated with call context.
* When a structural expression fails under the `data2` mask (for
  example `solve()` of a singular matrix), the error now reports
  that cause instead of a misleading "not found" from the fallback
  lookup.
* An unordered factor response to an ordinal family warns with the
  level order about to be used, then fits (the brms behavior);
  `mo()` keeps its stricter ordered-factor requirement.

# frmtmb 0.29.0

Review fixes for the spatial wave and rOpenSci review preparation.

## Breaking

* `spde()` grouping values are now explicitly MESH ROW INDICES:
  integer, whole-number character, and integer-level factor
  spellings are accepted and land identically; any other labels
  error stating the contract. Previously a non-mesh-order grouping
  (integer node IDs included, which sort lexicographically) fit a
  silently permuted field with a wrong likelihood and no warning.
  Code that passed non-numeric node labels must now map them to row
  indices.

## Corrected behavior

* `car()`/`spde()` accept call-valued grouping variables
  (`gr = factor(node)`).
* The Windows sequential-chains guard in `frm_sample()` also covers
  `options(mc.cores =)`, which rstan reads as its default.
* Optimizer error wrapping no longer claims a numerical cause for
  non-numerical errors; the underlying message leads.
* The `escar` log-determinant is finite at `rho` near 1 (the
  adjacency eigenspectrum is clamped against floating-point
  overshoot).
* The extreme-parameter heuristics in `diagnose()` and
  `frm_sample()` only test log-SD components and name the parameter;
  a `car`/`bym2` mixing proportion at its boundary no longer trips
  a misleading singularity warning.
* `car()` adjacency validation: NA entries and negative weights
  error informatively; one-sided dimnames work (the rownames-only
  path was broken); duplicate names error.
* Non-finite covariance warnings fire once per fit instead of once
  per `vcov()` call.

## rOpenSci preparation

* Review infrastructure: CONTRIBUTING, code of conduct, CITATION,
  codemeta, `Language` field, return-value documentation, README
  currency fixes.
* srr statistical-software standards: General and Regression
  categories scaffolded with 77 standards tagged at their code
  locations, 11 marked inapplicable with justifications, and the
  28 remaining gaps recorded as the submission roadmap in the
  standards file. Test coverage measured at 87.6% (full suite).

# frmtmb 0.28.0

Spatial GMRF grammar, the mclust covariance taxonomy, and the
close-out of the fuzz tier.

## New

* Spatial conditional autoregressions with the brms spelling:
  `car(M, gr = g, type = "escar"/"esicar"/"icar"/"bym2")` over an
  adjacency matrix, with analytic log-determinants throughout (no
  on-tape factorization) and brms's per-component soft sum-to-zero
  constraint, whose distance to the hard-constrained ML is measured
  and documented (4e-4 logLik at the default, shrinking
  quadratically). Validated against hand-rolled direct ML to 1e-11.
* `spde(fem, gr = node)`: an SPDE-Matern field from fmesher/INLA FEM
  matrices, Q(kappa, tau) assembled on the tape from fixed sparse
  matrices; sd and range reported through the planar identities.
* `gr(g, prec = Q)` supports correlated slopes (sparse Kronecker
  precisions; agrees with the dense `gr(cov = solve(Q))` equivalent
  to 0 ulp).
* `mixture_mvn(K, D, model =)` gains the mclust covariance taxonomy
  (EII, VII, EEI, VEI, EVI, VVI, EEE, VVV) with per-model parameter
  templates; matches `mclust::Mclust` to 1e-12 log-likelihood and
  1e-15 posteriors on intercept-only fits, while keeping
  covariate-dependent means, which mclust cannot fit.
* `anova(refit = TRUE)` refits REML models with ML for comparison
  (warm-started, with a message naming what was refit); `getME()` on
  the lme4 generic with the commonly consumed vocabulary.

## Corrected behavior

* Non-finite covariance matrices warn at `vcov()`/`summary()` and at
  fit time under `se = TRUE`, pointing at `diagnose()` (a
  positive-definite Hessian can still be numerically singular to
  invert; previously silent NaNs).
* Optimizer failures recover automatically where a better start
  provably exists (the plain-Laplace optimum under
  `profile = TRUE`; a cold start when a warm start is itself the
  failure), and every optimizer error carries the model label and
  remedies instead of a bare RTMB message.
* One quadrature configuration is rescued by displaced calibrations;
  the six irrecoverable thin-data configurations refuse with
  messages naming the configuration and the guaranteed
  `quadrature = FALSE` fallback. The fuzz tier runs green.
* `frm_sample(cores > 1)` on Windows falls back to sequential chains
  with a warning instead of failing incomprehensibly (RTMB tapes and
  objective closures cannot reach socket workers; upstream
  kaskr/tmbstan#27). `pp_check()` and `as_draws()` dispatch without
  attaching bayesplot or posterior (frmtmb ships its own generics,
  dual-registered).

## Notes

* dev/benchmarks.md records the RTMBp (parallel RTMB) measurement:
  no gain on Cholesky-dominated mixed models, 1.7-2x on
  accumulation-dominated GLMMs at 4-16 threads; not adopted.
* dev/feature-gaps.md records the ODE roadmap (RTMBode) and the
  Suggests + Additional_repositories packaging plan for non-CRAN
  extensions.

# frmtmb 0.27.0

Fixes from a full review of the v0.22-v0.25 waves.

## Corrected behavior

* `cens()` combined with `trunc()` produced a silently wrong
  likelihood: the censoring contribution was not restricted to the
  truncation window (right-censoring now contributes
  (F(ub) - F(y))/Z), which inflated the dispersion by ~14% on the
  reference problem. The corrected objective matches a hand-rolled
  windowed likelihood to 1e-8, collapses exactly to the old form
  without truncation, and the cens+trunc OSA residuals recalibrate
  to the analytic PIT (5e-14).
* `residuals(type = "deviance")` under `se()` now weights each
  row's unit deviance by its own variance (the glm prior-weight
  form) instead of treating a common dispersion as shared.
* `predict(se.fit = TRUE)` on quadrature fits no longer dies with a
  non-conformable error: it reports mode-conditional standard errors
  with a warning. Singular joint precisions in the predict path now
  degrade like `vcov()` does, and the shared solver uses
  `Matrix::solve` (strictly more robust for ill-conditioned but
  invertible GLMM precisions).
* `frm_sample()` chain inits are clamped inside any bounds passed to
  Stan, with a warning naming the parameters when the ML mode itself
  violates a bound (previously an incomprehensible rstan error).
* `vint()`/`vreal()` variables are required in `newdata`; a custom
  family's prediction previously returned a length-0 vector with no
  message when they were missing.
* Character censoring codes like "0"/"1"/"-1" decode (the error
  message had advertised them).

## Compatibility registry and fuzzer

* The registry's precedence is redesigned (lexicographic comparison
  of sorted side specificities, with a validator that forbids
  file-order ties except documented overrides). 488 of 3750
  resolutions were corrected against probed reality, including 429
  pairs whose covstruct conditions a family-level rule had erased,
  and the newly recorded caveat that spatial structures over a plain
  factor silently use level order as coordinates. `frm_compat()`
  accepts vectors and errors on empty input.
* The fuzz tier's invariants were strengthened (a previously
  unfailable interval check replaced by the Wald identity plus
  parameter-coverage tallies; per-row simulate agreement; grammar
  divergences asserted; convergence demotion keyed to actual
  convergence verdicts after fixing a truncation that discarded
  them). The strengthened tier reproduces only known findings.

# frmtmb 0.26.0

Pooled model comparison across imputations and the diagnostics/UX
backlog.

## New

* `anova()` on `frm_multiple()` fits pools nested-model tests across
  imputations: D1 (multivariate Wald), D2 (chi-square combining),
  and D3 (Meng-Rubin likelihood pooling, with the plug-in leg
  evaluated by re-taping each imputation's objective at the pooled
  parameters - no refits). D1/D2 match `mice` to 1e-7 on an
  exactly-shared reference; D3 is validated against the Meng-Rubin
  formula directly, since `mice::D3`'s `fix.coef` variant is not
  the plug-in statistic. Includes an ARIV clamp and a Reiter-df
  fallback for a small-m case where mice returns NaN.
* `cbind(successes, failures)` binomial responses are accepted (the
  glm/lme4/glmmTMB spelling), rewritten internally to
  `successes | trials(successes + failures)`; bit-identical to the
  `trials()` form.
* `simulate()` returns ordinal draws as ordered factors with the
  original levels and multinomial draws as count matrices (both
  families previously had no simulator), and respects `na.exclude`
  padding.
* `frmtmb_control(check_nlev_1 =, check_olre =)`: lme4-style
  warning/ignore/stop vocabulary for one-level grouping factors and
  gaussian observation-level random effects (the `se()`-based
  meta-analysis idiom is recognized and not flagged).
* `diagnose()` adds a complete-separation heuristic, predictor-scale
  warnings pointing at `autoscale`, and an isSingular-style verdict
  independent of the Hessian; it also no longer errors on fits
  without random effects.

## Corrected behavior

* Models with zero free outer parameters fit degenerately instead of
  dying inside nlminb (fixing latent empty-sdreport and
  empty-gradient bugs found along the way).
* A nonlinear-parameter name colliding with a data column errors
  instead of silently shadowing the column; nonlinear fits that fail
  from default zero starts name `start =` in the error.
* REML `anova()` compares fits whose fixed-effect designs span the
  same column space (term reordering included) instead of refusing
  every REML pair; genuinely different designs and REML/ML mixes
  still refuse with the reason.
* Offsets in distributional-parameter formulas were verified correct
  (to 1e-13) and are now regression-tested against the silent-drop
  failure mode reported upstream (glmmTMB#625).

# frmtmb 0.25.0

Simulation-workflow ergonomics, the last deferred method-surface
items, and CRAN readiness.

## New

* `frm_simulate()` accepts natural-scale parameter names - the same
  vocabulary as `hypothesis()`/`variables()`: coefficients by name,
  `sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`, response-scale
  `sigma` - inverted to the internal parameterization through the
  covariance registry (us/diag/homdiag/smooth/gr structures;
  others refuse by block name). The internal spelling still works
  and is byte-identical.
* Prior-predictive simulation, the `sample_prior = "only"` analog:
  `frm_simulate(..., priors = set_prior(...))` draws a fresh
  parameter vector per replicate on each prior class's documented
  scale (natural SDs for class `sd`, truncation by rejection) and
  attaches the drawn parameters for prior-predictive checks. Every
  coefficient and SD must be pinned by a prior or `newparams`;
  omissions error instead of silently becoming zero effects.
* `simulate(censored = TRUE)` applies the fitted type-I censoring
  mechanism to the draws; the default remains the latent uncensored
  response, which is also brms's `posterior_predict` convention
  (verified against brms source) and is now documented.
* `residuals(type = "deviance")` across the GLM family set
  (gaussian, poisson, binomial/bernoulli, Gamma, exponential,
  inverse.gaussian, nbinom1/2, geometric, beta, tweedie), exact
  against `stats::glm` where glm offers them; families without a
  standard unit deviance refuse by name. `deviance()` itself stays
  `-2 logLik` (lme4 convention).
* `predict(type = "response", se.fit = TRUE)` now works for
  non-identity-mean families (zi/hurdle, lognormal, trials-binomial,
  truncated responses) through a joint delta method across all dpar
  linear predictors, including cross-dpar covariance (agrees with
  glmmTMB's response-scale standard errors to 1.2e-5 on a
  zero-inflated poisson mixed model).

## Corrected behavior

* `conditional_effects(method = "predict")` evaluates addition terms
  on the effect grid, so its bands respect `trials()` and `trunc()`
  (binomial bands are counts in [0, n], not Bernoulli 0/1); aterm
  variables must be pinned in `conditions`, and the error names the
  variable.

## CRAN and infrastructure

* Heavy reference-validation test files are `skip_on_cran()`-gated:
  the CRAN-condition suite drops from ~128s to ~72s on the
  development machine while CI (NOT_CRAN=true) keeps full coverage.
* nlme added to Suggests (a test uses `nlme::Soybean` as reference
  data; CI's `--as-cran` unstated-dependency check halts without the
  declaration).
* Benchmark verdict recorded in dev/benchmarks.md: optimParallel is
  not adopted - with exact AD gradients its concurrency caps at two
  evaluations (measured 1.03-1.22x end to end, slower with cold
  clusters), RTMB tapes cannot ship to PSOCK workers, and 100% of
  InstEval's optimization time is inside the taped objective and the
  inner sparse Cholesky. `frmtmb_control(profile = TRUE)` remains
  the measured lever for many-coefficient models (1.6x there).

# frmtmb 0.24.0

The quadrature and OSA defect clusters surfaced by the fuzzer and
compatibility registry.

## Corrected behavior

* `quadrature = TRUE` is rebuilt around one root cause: TMB's
  Gauss-Kronrod marginalization calibrates each integrand's rescaling
  once, at the parameter values in hand when the tape is built, and
  frmtmb taped at the cold start. A plain Laplace fit now runs first
  and the quadrature tape is built at its optimum. This fixes:
  conditional modes returned as NA (`ranef()`/`fitted()`/`predict()`
  were silently NA for all groups but the first, whose slot held a
  wrong value) - modes now come from the Laplace inner solve and
  match `glmer(nAGQ = 25)`'s `ranef()` to 3e-5; and bare "NA/NaN
  gradient evaluation" crashes for poisson/Gamma/Beta with nested or
  even single scalar blocks - all now fit with gradients < 3e-4.
* `quadrature` combined with `trunc()` produced logLik = +Inf as a
  successful fit; the combination is refused (the CDF difference
  underflows at quadrature nodes; a stable fix needs log-CDF forms).
  `mixture()` with `REML` or `profile = TRUE` is refused: both
  Laplace-expand the mu coefficients around a single mode, and a
  mixture likelihood is permutation-multimodal in exactly those
  coefficients. `mixture()` with quadrature remains supported.
* `residuals(type = "osa")` on censored fits was broken (raw LAPACK
  singularity or NaN for every censored row). Uncensored rows now get
  calibrated one-step residuals by conditioning on the censoring
  window and renormalizing the one-step CDF to it (matches the
  analytic conditional PIT to 7e-9); censored rows are NA (an event
  has no one-step CDF), and row-varying censoring points refuse.
  Note `simulate()` draws the latent uncensored response under
  `cens()`, so DHARMa is not a substitute there (documented).
* `residuals(type = "osa")` on ordinal fits crashed inside the OSA
  machinery; the ordinal log-densities now handle `oneStepPredict`'s
  taped-observation objects (exact Lagrange-basis category
  indicator), and all four ordinal families produce calibrated
  residuals (cumulative matches the analytic randomized-quantile
  residual to 4e-14).
* A singular joint precision under REML or `profile = TRUE` no longer
  throws a raw LAPACK error from `vcov()`/`summary()`/`confint()`/
  `hypothesis()`: it degrades to NaN standard errors with one warning
  pointing at `diagnose()`, matching the ML branch.
* The brms-migration vignette documents that `binomial()` without
  `trials()` means Bernoulli here (glm convention) where brms
  requires `trials()` or `bernoulli()`.

The fuzz tier's open findings drop from 28 to 16 on the identical
plan; most of the remainder are now informative refusals on thin-data
edge cases rather than defects.

# frmtmb 0.23.0

Defect wave driven by an open-issue sweep of brms/lme4/glmmTMB, plus
a feature-compatibility registry.

## Corrected behavior

* Truncation now reaches the whole post-fit surface, not just the
  likelihood: `fitted()`, `predict(type = "response")`, and
  `residuals()` report the truncated mean E[Y | lb <= Y <= ub]
  (closed forms for all six CDF families, validated to 1e-15);
  `simulate()`, `posterior_predict()`, and `frm_simulate()`
  rejection-sample within bounds; newdata re-evaluates variable
  bounds. `residuals(type = "osa")` on truncated (and untruncated
  gaussian) responses was also miscalibrated and now integrates over
  the truncated support (KS uniformity restored from p ~ 1e-15 to
  0.76-0.97). Dpar-scale predictions stay untruncated by design.
  DHARMa/pp_check on truncated models are meaningful again.
* `(f || g)` with a factor now yields independent per-level effects
  (the diag structure the syntax promises) instead of a silently
  fully correlated block; identical to an explicit `diag(f | g)`.
  Numeric double-bars are unchanged. Existing factor-double-bar fits
  change, because the old ones were wrong (lme4 has the same open
  bug, lme4#818).
* `ar1()`/`hetar1()` warn when the ordering factor's integer levels
  have gaps: levels correlate by position (the glmmTMB reading,
  unchanged), so a gap counts as one step; the warning points at
  `ou()` over `num_factor()` for irregular spacing (glmmTMB#1278).
* Predictions at non-estimable points of a rank-deficient design
  return NA with one warning naming the aliased columns, instead of
  silently returning the partial sum (predict.lm semantics;
  lme4#303). Collinear restatements of kept columns stay exact.
* Grouping factors written as calls work: `(1 | factor(x))`,
  `(1 | interaction(a, b))` (lme4#464).
* A random-effect term crossed with `*` or `:` errors instead of
  being silently refit additively (lme4#196); `mo()`/`mi()`
  interaction multipliers must be numeric (brms#1828); `anova()`
  requires equal `nobs` (lme4#622).
* `frm_sample()` chain initialization: chain 1 anchors at the ML
  mode, further chains are jittered on the unconstrained scale
  (`init_jitter`), restoring the overdispersion Rhat needs; a
  boundary-mode warning fires for singular fits; mixture posteriors
  are documented as needing `init = "random"`.

## New

* Feature compatibility registry: `frm_compat()` answers what plays
  with what (works / conditional / refused / broken / untested)
  across 3750 feature pairs, and the new "Feature compatibility"
  article is generated from the registry at build time so it cannot
  drift. Known-broken pairs are listed there; the registry probing
  itself surfaced trunc x quadrature and cens x OSA as broken
  (queued for the next fix wave).
* An audit of 559 currently open brms/lme4/glmmTMB issues is
  recorded in dev/test-backlog.md, including the pathologies frmtmb
  is structurally immune to.
* A pairwise grammar fuzzer (env-gated: `FRMTMB_FUZZ=true`): a
  310-spec covering array over the grammar with metamorphic
  invariants (predict/fitted identity, permutation invariance,
  simulator support membership, vcov sanity) and a brms
  `make_standata` structural oracle. Its first run surfaced the
  quadrature defect cluster now queued for fixing (conditional
  modes left NA, crashes with nested groups and Beta, +Inf logLik
  with trunc()); the tier reports those as failures until they are
  fixed.

# frmtmb 0.22.0

Cross-validation against brms itself, closer brms compatibility, and
fit ergonomics.

## brms agreement suite

* New test tier validating the elaborate grammar structurally against
  `brms::make_standata()`/`brmsterms()` without any Stan compilation
  (designs, RE structures, `|ID|` merging, ordinal thresholds, `cs()`,
  `mo()` codes, gp bases, smooths, addition terms, mvbf, nl, mixture
  naming; most exact to 1e-10 or better), plus an opt-in numeric tier
  (`FRMTMB_BRMS_FIT_TESTS=true`) showing fixed-effects estimates match
  the mode of brms's own generated Stan program to 1e-4 and that our
  `mo()` parameterization is brms's likelihood exactly (via
  `rstan::log_prob`). brms is in Suggests only.
* Found by the suite and fixed: `cs()` predictor variables never
  reached the combined model frame, so `bf(y ~ x + cs(z))` errored
  unless `z` already appeared elsewhere.

## brms compatibility

* Hilbert-space `gp(x, k =)` now uses brms's input convention exactly:
  coordinates are rescaled by the largest pairwise distance (one
  shared factor across dimensions) and centered, so the boundary is
  `L = c` and the same `gp()` call is the same approximation in both
  packages (basis agreement with brms to 1e-16, tied coordinates and
  vector-valued `c =` included). This roughly doubles the effective
  boundary at a given `c`: accuracy against the exact gp improves
  about 25x at the same k (k = 40 now within 2e-4 logLik), and the
  default `c = 1.25` is the right choice in 2-D as well. Reported
  `range(gp)` values stay in data units (exact log-shift
  back-transform); fitted lengthscales differ from pre-0.22 fits by
  the data-dependent scale factor.
* `cens()` accepts brms's character and factor codes
  (`"left"`/`"none"`/`"right"`/`"interval"`, prefix-matched,
  case-insensitively); unknown labels and out-of-range numeric codes
  error informatively instead of coercing to NA.
* `gp()`'s `k`/`c`/`iso`, `rr()`'s `d`, and `se()`'s `sigma` arguments
  now evaluate in the formula environment (brms behavior), so
  variables and expressions work; invalid values error with the
  offending expression named.

## Fit ergonomics

* `frmtmb_control(verbose =)` (and the `frm(verbose =)` shortcut):
  level 1 prints timed stage lines (parse, frame, tape, optimize,
  restarts, sdreport) through `message()` so a slow fit shows where
  it is slow; level 2 adds the optimizer's own iteration trace.
  Bootstrap, influence, and allfit refit loops stay quiet.
* The large-gradient warning now points to `diagnose()` and the new
  "Convergence problems" remedy ladder in the diagnostics vignette
  (scaling, restarts, optimizer comparison, profiles, boundary fits,
  and judging marginal gradients).
* `frm_allfit()`'s nloptr optimizer produced an empty row: nloptr
  rejects callbacks that declare `...` (RTMB's `obj$fn`/`obj$gr` do),
  and the error was swallowed; separately, missing `ftol_rel` made
  NLopt report failure at the converged optimum. Both fixed; the
  agreement test now requires every optimizer to converge and match.
* The brms-migration vignette documents how to recover the model
  function (the `stancode()` analog): the objective closure via
  `build_objective(fit$frame)`, its joint-vs-marginal relationship to
  `fit$obj$fn`, and how to reproduce `fit$obj`.

# frmtmb 0.21.0

Method-surface audit against lme4/glmmTMB/brms, plus fixes from a
full code review.

## Corrected behavior

* `predict(type = "response")` now returns the expected response for
  every family (equal to `fitted()` in-sample), not the primary
  dpar's natural scale; zi/hurdle, lognormal, and trials-binomial
  fits were affected. The glmmTMB aliases `"conditional"`, `"zprob"`,
  `"zlink"`, and `"disp"` are accepted (validated against glmmTMB on
  a zero-inflated Poisson to ~8e-7). Note `"disp"` returns `sigma` on
  its natural scale, where glmmTMB returns the variance-scale
  dispersion for gaussian. `se.fit` for the mean of a
  non-identity-mean family is not yet available and errors with
  guidance.
* `rescor = TRUE` now refuses `cens()`, `trunc()`, and `se()`
  addition terms, which the joint-gaussian likelihood silently
  ignored (wrong likelihood with no warning).
* Bounds (`lower`/`upper`, `set_prior` lb/ub) under
  `frmtmb_control(profile = TRUE)` were positionally misaligned:
  a bound on a fixed coefficient could silently pin a covariance
  parameter instead. Bounding a profiled coefficient now errors;
  covariance-parameter bounds land on the right slot.
* `confint()`, `vcov(full = TRUE)`, and `diagnose()` labels were
  broken or misaligned on `mi()` fits: the latent imputation
  component was counted as an outer parameter.
* `frm_sample(priors = )` failed on fixed-effects-only fits of
  single-dpar families (a `$` partial-match bug), and
  `frm_sample(laplace = TRUE)` was broken twice: the default init
  had the wrong length, and draw columns were mislabeled (a column
  named `b[1]` actually held a covariance parameter).
* `influence()` tables and `cooks.distance()` now align with
  `vcov()` (fits with a constant dpar errored; column names are now
  the coefficient labels).
* `simulate()` follows the stats seed contract: a `"seed"` attribute
  and a restored RNG state instead of clobbering the global stream.
* A covariate literally named `sigma` is no longer shadowed by the
  residual sigma in `hypothesis()`/`variables()`.
* In `predict(newdata = )`, an NA in a variable used only in a
  random-effect design now propagates to the prediction instead of
  being silently zeroed (only genuinely new levels predict at the
  population value).

## New methods

* `drop1()` (AIC/LRT, marginality-aware; matches
  `lme4::drop1.merMod` to 1e-4), `cooks.distance()` directly on a
  fit, `dfbeta()`/`dfbetas()` on `influence()` results (stats sign
  convention; lme4 returns the negation), `na.action()`, and an
  lme4-style "Groups:" line in `summary()`.
* `confint()` accepts lme4's `"Wald"` spelling and gains
  `method = "boot"` (percentile intervals via [frm_bootstrap()]).
* `vcov(full = TRUE)` rows carry per-parameter names (glmmTMB
  convention).
* `predict()` accepts the lme4/glmmTMB `allow.new.levels` dot
  spelling.
* On draws objects: `posterior_linpred()` and `ranef()`
  (brms-shaped arrays).
* `emmeans::recover_data` had been silently missing from NAMESPACE
  since v0.13 (orphaned export tag); emmeans support was broken.

## Internals

* One authoritative outer-parameter map shared by `confint()`,
  bounds resolution, `vcov(full = TRUE)`, and sampling; `graphics`
  and `grDevices` added to Imports; duplicated formula/family
  coercion and gp position-key helpers unified; REML `logLik` df
  counts only outer parameters (documented; lme4 counts the
  integrated fixed effects too, so REML AICs are not comparable
  across packages).

# frmtmb 0.20.0

* `mixture_mvn(K, D)`: multivariate gaussian mixture components for
  model-based clustering of an n x D matrix response (the mclust /
  clustTMB use case). Every class mean is a full linear predictor, so
  cluster means may depend on covariates and random effects; mixing
  weights are multinomial-logit dpars with their own formulas.
  Validated against a hand-rolled ML fit to 1.3e-10 and against
  faithful-data cluster recovery. Limitations (documented in
  `?mixture_mvn`): class covariances are unstructured (`us`) and
  covariate-free - no mclust-style constrained covariance taxonomy
  (EII..VEV); `cens()`/`trunc()`, `mvbf()`, and `simulate()` are not
  supported.
* Multi-dimensional Gaussian processes: `gp(x1, x2, ...)` takes up to
  three variables, with a separate lengthscale per dimension by
  default (the brms convention) or one shared lengthscale via
  `iso = TRUE`, in both the exact and the Hilbert-space (`k =`) form.
  The tensor HSGP basis is capped at `k^D <= 1000` columns.
  `confint_varcorr()` reports one range row per dimension. Validated
  against a direct 2-D ML fit to 3.8e-11.
* Kriging prediction for the exact `gp()`: `predict(newdata = )` at
  positions not seen at fit time now returns the GP conditional mean,
  and `se.fit = TRUE` adds the conditional variance, instead of
  erroring. Kriging mean validated to 1.7e-15 and the standard error
  to the exact decomposition identity. The Hilbert-space form already
  interpolated through its basis.
* `frmtmb_control(sparse_x = TRUE)`: sparse fixed-effect design
  matrices (`Matrix::sparse.model.matrix`), the analog of
  `glmmTMB(sparseX = )`, for models where a many-level fixed factor
  dominates memory. Estimates are identical to the dense path
  (gradient agreement 3e-14; 13.8% lower tape memory in the
  validation model); `model.matrix()` on the fit returns a sparse
  matrix. Prediction falls back to the dense builder for newdata
  containing NA factor rows, which `sparse.model.matrix` would
  silently zero.
* `frmtmb_control(autoscale = TRUE)`: internal standardization of
  badly scaled continuous predictors (the lme4 >= 1.1.37 feature): a
  scaled pre-fit is mapped back exactly and warm-starts the reported
  fit, parameter magnitudes feed `nlminb`'s scale hook, and
  convergence and standard-error Hessian checks run in scaled units.
  Turns 12 scaling warnings into 0 on the validation problem with a
  0.0 log-likelihood difference against a manually standardized
  reference. Reported results are always on the original scale.
* `frm_multiple()` pooling extensions: `$pooled_varcorr` pools
  random-effect SDs and correlations across imputations on
  transformed scales (log for SDs and GP ranges, Fisher z for
  correlations) with Barnard-Rubin degrees of freedom, and
  `hypothesis()` on a `frmtmb_multiple` object pools delta-method
  hypothesis estimates by Rubin's rules with t-based intervals.
  Fixed-effect pooling matches `mice::pool()` to 2.1e-6.
  Likelihood-ratio pooling across imputations (mice's D3) is not
  implemented; `anova()` on pooled fits remains per-fit.

# frmtmb 0.19.0

* Latent classes now combine with continuous random effects
  (growth-mixture models): `mixture(..., groups = ~g)` accepts
  random effects, smooths, and gp() terms in the component formulas.
  The class sum is computed conditional on the latent effects, so
  one Laplace approximation integrates them - the sum-and-integral
  swap makes this exact as an integrand identity, and even crossed
  random effects are structurally fine. Random effects written in a
  component formula are class-specific by construction. Accuracy:
  the Laplace approximation of the class-mixture integrand carries a
  small documented bias (about 0.1 log-likelihood units in the
  validation problem, with parameter estimates matching the exact
  closed-form marginal to 0.01); `quadrature = TRUE` is numerically
  exact when the per-group integrand is univariate (validated to
  5e-3) and approximate when class-specific intercepts couple.
  `simulate()` now draws one class per group and then simulates each
  observation from its group's component; `mixture_probs()` gives
  empirical-Bayes classification conditional on the modes.
* `variables()` (brms spelling): lists every parameter name usable in
  `hypothesis()` expressions (coefficients, `sd_`/`cor_` summaries,
  `sigma`); on `frm_sample()` output it lists the draw columns.
* `get_prior()` (brms spelling): enumerates every slot `set_prior()`
  can target (class/coef/dpar/group rows), from a formula plus data
  or from a fit; the default in every slot is flat.
* README rewritten for the current scope, with related-work pointers
  (glmmTMB, BayesRTMB); getting-started vignette extended to the full
  grammar; SPEC.md status and deviations brought current.

# frmtmb 0.18.0

The last three roadmap features.

* `gp(x)` Gaussian-process terms: exact (a dense squared-exponential
  block over the unique positions, with a standard 1e-6 nugget) or
  the Hilbert-space approximation via `gp(x, k = 30)` (sine basis
  with spectral-density prior SDs), in any linear predictor. Exact
  fits match direct GP marginal ML; the approximation converges to
  the exact answer and predicts at arbitrary positions (the exact
  form predicts at observed positions only). `confint_varcorr()`
  reports `sd(gp)` and `range(gp)`.
* `mo()` and `mi()` interactions: `mo(x):z`, `mo(x)*z`, `mi(x):z`,
  `mi(x)*z` with numeric multipliers; `mo()` interactions share
  their variable's simplex (brms convention). Both validated exact
  against direct ML (the `mi()` interaction stays linear in the
  latent value, so the closed-form marginal still applies).
* Group-level latent-class mixtures: `mixture(..., groups = ~g)`
  sums the class assignment per group (growth-mixture /
  latent-class regression; the tractable nested case of brms#1905).
  Exact against direct ML; `mixture_probs()` returns posterior class
  probabilities per group, or per observation for ordinary mixtures.
  Restrictions: no random effects or smooths alongside, mixing
  predictors evaluated at each group's first row, no `simulate()`
  yet. Crossed-design group mixtures remain out of the Laplace
  class.

# frmtmb 0.17.0

Clearing the roadmap: mixtures, measurement error, and the remaining
grammar gaps.

* `mixture(fam1, fam2, ...)` finite-mixture families: each component
  keeps its own suffixed dpars (`mu1`, `sigma1`, ...), mixing
  proportions are multinomial-logit dpars with their own linear
  predictors (mixture-of-experts comes free), and the likelihood is a
  branch-free logsumexp - so random effects and smooths work in any
  component formula. Exact against direct ML for gaussian and poisson
  mixtures. Component means initialize on spread response quantiles;
  the usual finite-mixture multimodality caveats apply
  (`frm_allfit()`, or order intercepts via bounds).
* `mi(sdx)` measurement error (brms `me()`): known per-observation
  measurement SDs make every true value latent, with the observed
  values entering through a measurement model. Exact against the
  closed-form bivariate-normal marginal; attenuation bias is
  corrected and NAs (missing + mismeasured) combine.
* `cs()` category-specific effects for `sratio`, `cratio`, and `acat`
  (refused for `cumulative`, as in brms). Exact against direct ML.
* `equalto(x + 0 | g, V)` covariance structure: fully fixed known V,
  zero parameters (meta-analytic sampling covariances).
* rr() fits now support `se.fit`: the delta method routes the factor
  columns through the loadings Jacobian and adds loading-parameter
  columns; verified parameterization-invariant against the
  full-rank-vs-`us()` equivalence.
* Registered insight methods (`find_formula`, `find_random`,
  `get_parameters`, `get_varcov`, `find_statistic`, link accessors):
  `insight::is_mixed_model()` and the easystats accessor layer now
  see the random-effect structure.

# frmtmb 0.16.0

The two deferred architecture features.

* `rr(x | g, d = k)` reduced-rank / factor-analytic covariance
  (glmmTMB's GLVM structure): per-level coefficients are loadings
  times iid standard-normal factors, so `b` holds the factors while
  the design matrices span the coefficient space. Matches glmmTMB
  exactly and reduces to `us()` at full rank. `se.fit` on rr models
  is deferred (the loadings enter the coefficients nonlinearly);
  `predict`, `simulate`, `ranef`, and `VarCorr` (the rank-k
  covariance) all work.
* `mi()` one-step imputation for continuous predictors, reversing the
  original exclusion: `bf(y ~ mi(x) + z) + bf(x | mi() ~ z)` turns
  missing `x` values into latent parameters that Laplace integrates
  (exactly, in the linear-gaussian case - validated against the
  closed-form marginal likelihood). Imputation uncertainty
  propagates into the coefficient standard errors automatically.
  Rows with NAs in non-`mi()` variables still drop; imputation
  models must be gaussian or student; discrete predictors remain
  impossible (as in Stan). `me()` measurement error is the natural
  follow-on.

# frmtmb 0.15.0

Tier-2 sweep of dev/feature-gaps.md plus a method surface for sampled
fits.

* `mo()` monotonic effects (brms syntax): scale coefficient times an
  estimated simplex, exact against direct ML; works in prediction,
  `conditional_effects()`, and with ordered factors. No frequentist
  package offers this. Standalone additive terms only for now.
* Sequential ordinal families `sratio`, `cratio`, `acat`, validated
  against direct ML and collapsing to logistic regression at K = 2.
* Covariance structures: `hetar1`, `homcs`, `homtoep` (glmmTMB
  parameterizations, exact where glmmTMB converges) and spatial
  `exp()`, `gau()`, `mat()` over `num_factor(x, y)` coordinates
  (Matern uses bounded internal transforms for stability where both
  we and glmmTMB otherwise diverge).
* Fixed a latent parser bug: `exp(x)` (or any covariance-structure
  name) used as a plain function in a fixed formula was silently
  stripped by the formula splitter - since v0.1. Such calls are now
  protected automatically.
* `influence(fit, groups = )` and `cooks.distance()`: case-deletion
  diagnostics over warm-started refits.
* `frm_multiple()`: fits across multiply-imputed datasets (list or
  `mice::mids`) pooled by Rubin's rules with Barnard-Rubin df.
* `conditional_effects()` gains `method = "predict"` (prediction
  intervals) and data-frame `conditions` (one condition set per row,
  brms style).
* `vint()`/`vreal()` addition terms pass arbitrary data vectors to
  custom families (brms custom-family convention). The test suite
  includes a full Wiener drift-diffusion model written as a
  `custom_family()` in plain R - the workflow brms needs raw Stan
  code for.
* Method surface for `frm_sample()` draws: `summary()` (with
  Rhat/ESS), `fixef()`, `VarCorr()` (natural scale), `hypothesis()`
  (exact posterior version), `posterior_epred()` /
  `posterior_predict()` (each draw runs the full prediction
  machinery), `pp_check()`, and `posterior::as_draws()`.
* Examples added across the reference (run by R CMD check).

# frmtmb 0.14.0

Gap-closing milestone from a sweep of the lme4 / glmmTMB / brms
vignettes (dev/feature-gaps.md).

* `se()` addition term: known per-observation sampling SDs
  (meta-analysis), gaussian and student. `se(x)` fixes the residual SD
  (sigma is mapped out); `se(x, sigma = TRUE)` adds an estimated sigma
  in quadrature. With `(1 | obs)` this is random-effects
  meta-analysis, and with `gr(g, cov = A)` the phylogenetic version.
* Binomial (and beta-binomial) responses may be proportions of
  `trials()`, the glm/glmer idiom, converted to counts internally.
* New families: `bernoulli`, `geometric`, `exponential`, `weibull`
  (mean-parameterized, matches survreg's likelihood),
  `shifted_lognormal`, `hurdle_gamma`, `hurdle_lognormal`,
  `zero_inflated_binomial`, `zero_inflated_beta` (validated against
  glmmTMB), and `asym_laplace` for quantile regression (fixed
  `quantile` dpar reproduces `quantreg::rq` estimates).
* `ranef(condVar = TRUE)`: conditional SDs of the modes (TMB/glmmTMB
  convention: fixed-effect uncertainty propagates), plus tidy
  `as.data.frame()` methods for `ranef` and `VarCorr` output.
* `frm_allfit()`: refit under every available optimizer (nlminb,
  optim, bobyqa, NLopt) and compare, the lme4 `allFit()` analog.
* `frm_simulate()`: de novo simulation from a `bf()` formula and
  supplied parameters with no fitted model (power analysis; the
  glmmTMB `simulate_new()` analog).
* `frmtmb_control(profile = TRUE)` (experimental): profile the
  primary-coefficient vector into the inner Laplace problem, the
  `glmmTMBControl(profile = TRUE)` / `glmer(nAGQ = 0)` speed
  approximation for many-coefficient models; covariance comes from
  the joint precision.
* `anova()` no longer fails on models sharing a primary formula
  (distributional vs plain fits); row labels now include dpar
  formulas.

# frmtmb 0.13.0

Method-surface ("sugar") milestone: the conventional S3 methods that
downstream packages and user muscle memory dispatch on.

* New accessors: `sigma()`, `terms()`, `weights()`, `model.matrix()`,
  `deviance()`, `extractAIC()` (enables `step()` / `drop1()`),
  `ngrps()`, `prior_summary()`, and a `profile()` method wrapping
  `TMB::tmbprofile()`.
* **Breaking**: `coef()` now follows the lme4 / glmmTMB / brms
  convention (per-group coefficients: fixed effects broadcast over
  levels plus the conditional modes). Use `fixef()` for the fixed
  effects alone. Fits without random effects still return the
  coefficient vector.
* `conditional_effects()` with a base-graphics `plot()` method: grid
  predictions per predictor (or `"x:z"` pair) with Wald bands, smooths
  included, matrix covariates held at column means.
* `plot(fit)`: Pearson-residual diagnostics (residuals vs fitted,
  normal QQ).
* `pp_check()` (registered on the bayesplot generic): predictive checks
  from `simulate()` draws through any `ppc_*` function.
* `hypothesis()`: tests of arbitrary expressions of the parameters,
  brms spelling, with three methods: delta-method Wald (default),
  `"profile"` (profile-likelihood intervals via `TMB::tmbroot()`
  lincombs; linear hypotheses, ML fits), and `"boot"` (parametric
  bootstrap percentile intervals; any expression). The expression
  environment includes natural-scale random-effect names
  (`sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`, `sigma`), so ICC
  inference is `hypothesis(fit, "sd_g__Intercept^2 /
  (sd_g__Intercept^2 + sigma^2)", method = "boot")`. The result is a
  `frmtmb_hypothesis` data frame carrying the bootstrap draws or
  profile curves in attributes, with a `plot()` method (bootstrap
  histogram, profile curve, or the implied Wald normal, per
  hypothesis).
* `frm_bootstrap(fit, FUN, nsim)`: parametric bootstrap over
  `refit()`s (the `bootMer` analog), with `confint()` percentile
  intervals.
* `refit(fit, newresp)`: refit to a new response reusing the assembled
  design with a warm start - one re-tape plus optimization; the engine
  for parametric bootstrap over `simulate()` draws.
* The accessor set makes insight's default methods work
  (`find_response()`, `get_data()`, `get_sigma()`, `n_obs()`, ...);
  dedicated insight methods (random-effect formula parts,
  `get_parameters()`) are a roadmap item.

# frmtmb 0.12.0

Pre-release consolidation toward CRAN. Highlights by milestone:

## Model grammar (brms-compatible spelling)

* `bf()` with distributional-parameter formulas carrying the full
  predictor grammar (random effects and smooths in any dpar), constant
  dpars, and addition terms `weights()`, `trials()`, `cens()` (left /
  right / interval), `trunc()` (with the inclusive discrete correction).
* Nonlinear formulas (`nl = TRUE`) with a full linear predictor per
  parameter.
* Multivariate models: `mvbf()` / `bf() + bf()` / `mvbind()`, a family
  per response, gaussian residual correlation (`set_rescor()`), and
  `|ID|` random-effect correlation across formulas.
* mgcv smooths `s()` / `t2()` in any linear predictor, including
  matrix-covariate summation-convention terms: scalar-on-function,
  function-on-scalar, and function-on-function regression validated
  against mgcv.
* Covariance structures: `us`, `diag`, `homdiag`, `cs`, `ar1`, `toep`,
  `ou` (positions via `num_factor()`), and known-structure terms
  `gr(g, cov = A)` (dense, with correlated slopes via a Kronecker
  product) and `gr(g, prec = Q)` (sparse GMRF).

## Families

gaussian, poisson, binomial, Gamma, lognormal, student, negbinomial,
nbinom1, beta, tweedie, compois, beta_binomial, skew_normal,
inverse.gaussian, exgaussian, zero-inflated and hurdle counts,
cumulative ordinal, matrix-response multinomial, and `custom_family()`
as a plain R log-density with `check_custom_family()` AD verification.

## Estimation and inference

* ML and REML by Laplace approximation; `quadrature = TRUE` upgrades
  scalar random effects to adaptive Gauss-Kronrod marginalization
  (matches `glmer(nAGQ = 25)` and GLMMadaptive).
* MAP / regularized ML via brms-style `set_prior()`; hard bounds
  (`lower` / `upper`, and per-prior `lb` / `ub`).
* `sdreport` is deferred until standard errors are first needed
  (roughly a quarter off fit time; `se = TRUE` restores eager mode).
* Pluggable optimizers (`frmtmb_control(optimizer = )`), including
  arbitrary user functions.
* `confint()` (Wald / profile / likelihood-root), natural-scale
  `confint_varcorr()`, `anova()` likelihood-ratio tests, `diagnose()`.

## Diagnostics and ecosystem

* Residuals: response, pearson, and one-step-ahead (`type = "osa"`).
* DHARMa simulation-based residuals (`dharma_residuals()`).
* NUTS on the fitted objective: `frm_sample()` (ML-mode
  initialization, named draws, priors, bounds) and `check_laplace()`;
  bayesplot works on the draws.
* emmeans and marginaleffects support, both verified against glmmTMB.

## Verification

Roughly 500 tests compare fits against glmmTMB, lme4, mgcv, MASS,
survival, nnet, GLMMadaptive, and hand-written RTMB references, plus
an edge-case suite mined from the lme4 / glmmTMB / brms issue
trackers (dev/test-backlog.md).
