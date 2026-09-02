# Changelog

## frmtmb 0.35.0

Hidden Markov models, latent class analysis, multi-membership random
effects, three new families (nominal, circular, proportional hazards),
profile and bootstrap effect bands, and a brms-portability batch
measured against the brms vignettes themselves.

### Behavior changes

- [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) defaults
  to `family = gaussian()` when neither the `family` argument nor a `+`
  attachment supplies one, the brms/lme4/glmmTMB convention.
  `frm(y ~ x, data = d)` was an error and is now a linear model. The
  `family` argument overrides a family attached to a univariate
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) and fills
  only the empty responses of an
  [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md).
- Documentation and examples now lead with the separate-family spelling
  `frm(bf(y ~ x), family = gaussian(), data = d)`; the `+` attachment
  stays valid and documented as the alternative.

### Portability (measured against the brms vignettes)

- A scorecard audit ports every model call of the brms vignettes through
  the mechanical `brm` to `frm` transform (dev/brms-vignette-port.md):
  with this release’s fixes the measured tally is 30 of 42 model calls
  running mechanically, the remaining spelling changes documented in the
  porting guide.
- [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  accepts brms’s directional form (`"a > b"`) with one-sided p-values
  and bounds across every method, plus the `class=`/`group=` naming
  shorthand. [`update()`](https://rdrr.io/r/stats/update.html) speaks
  `formula.`/`newdata` and dotted deltas (`. ~ . + z`).
  [`lf()`](https://aforren1.github.io/frmtmb/reference/lf.md) composes
  dpar formulas onto
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md). One
  formula may name several parameters (`a + b ~ 1`). A bare constructor
  (`family = cumulative`) works.
- One parameter-addressing vocabulary across
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md),
  [`profile()`](https://rdrr.io/r/stats/profile.html), `confint(parm =)`
  and bounds: parenthesized and paren-stripped spellings are
  interchangeable, and one-to-one natural-scale names (`sd_g__x`, `ar1`)
  reach their internal parameter with the scale stated. Ambiguous
  natural-scale names are refused with the alternatives named.

### Multi-membership and effect bands

- Multi-membership random effects:
  `(x | mm(g1, g2, weights =, scale =))` and `mmc()` member-specific
  covariates, following brms; the pooled-level Z construction matches
  [`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html)
  exactly, and every post-fit method reads the block unchanged.
  New-level prediction variance distinguishes a shared new level (one
  draw, weights add) from distinct new levels (independent draws).
- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  gains `band = c("wald", "profile", "boot")`: likelihood-ratio bands
  inverted per grid point, and pointwise percentile bands from one
  shared parametric bootstrap (reused across effects and verified
  against the grid it was run over). Effect discovery now works on
  nonlinear and `mo()`/`mi()` fits, and
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  works on ordinal fits (previously all NA).

### New families

- [`categorical()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  fits a multinomial logit to an unordered factor, brms’s spelling of
  the likelihood `multinomial(K)` already carried on a count-matrix
  response. The first level is the reference and each remaining level
  gets its own linear predictor named `mu<Level>`, as in brms; the main
  formula applies to all of them unless a dpar formula overrides one, so
  `bf(y ~ x, mustout ~ w)` gives one category its own predictor.
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  `predict(type = "response")` return the `n x K` matrix of category
  probabilities with the response’s own levels as column names, the same
  convention the ordinal families follow, and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws factor
  levels. Validated against
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) (4e-10
  in the log-likelihood) and against `multinomial(K)` on the one-hot
  response (1e-8). The categories are read off the data by
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) before
  the formula is parsed; `categorical(levels =)` or `categorical(K =)`
  states them for the paths that have no data.
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) is refused: a
  nominal response has no scale for one.

- `von_mises(link = "tan_half")` for a circular response in `(-pi, pi]`,
  with the mean direction `mu` through the tan-half link and the
  concentration `kappa` through a log link - brms’s parameterization.
  The normalizing constant’s `log I0(kappa)` is differentiated exactly
  by RTMB’s own `besselI` method, so nothing is a series approximation.
  `kappa ~ x` distributional models work, and the simulator is a
  vectorized Best-Fisher sampler that varies both parameters per row.
  Validated against a hand-rolled RTMB likelihood (1e-6) and
  [`circular::mle.vonmises()`](https://rdrr.io/pkg/circular/man/mle.vonmises.html).

- [`cox()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  fits proportional hazards with brms’s flexible baseline: an M-spline
  baseline hazard over a simplex of weights, with the cumulative
  baseline hazard the I-spline integral of the same basis (`df = 5`
  cubic splines with an intercept, brms’s `bhaz()` default). Censoring
  runs through the ordinary `cens()` addition term - an event
  contributes the density and a censored row the survivor function - so
  right, left, and interval censoring all work, and frailty models come
  free through the Laplace approximation:
  `time | cens(c) ~ x + (1 | g)`. New
  [`cox_baseline()`](https://aforren1.github.io/frmtmb/reference/cox_baseline.md)
  returns the fitted baseline weights. Validated exactly against a
  hand-rolled M-spline PH likelihood (1e-6 in the log-likelihood, 1e-4
  in the coefficients) and against
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) to
  2e-2 on the log hazard ratio, which is the honest claim: `coxph()`
  leaves the baseline fully nonparametric.
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  `predict(type = "response")` and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) are refused - a
  survival time has no mean the censored rows identify, and brms refuses
  the same question. Unpenalized maximum likelihood often puts a
  baseline weight on the simplex boundary, which the optimizer reports
  as singular convergence even at the optimum; `?frmtmb-families`
  explains it and lowering `df` is the remedy.

- New `tan_half` link, for
  [`von_mises()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md).

### Internal

- A family may declare `aterm_data(y, aterms)`, family-level data that
  no addition term supplies, built once at frame assembly from the
  validated response.
  [`cox()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)’s
  spline bases use it.

- A family’s `lcdf()` may take the family-level extra parameters as a
  fourth argument, which is what lets a survival family’s survivor
  function reach its baseline. The three-argument contract is unchanged
  for every other family.

Hidden Markov models as a first-class family.

### New

- `hmm(K, family, time =, group =, init =, trans =)` fits a `K`-state
  hidden Markov model. The response at each time point comes from one of
  `K` state-dependent copies of `family`, and the unobserved state
  follows a first-order Markov chain along `time` within `group`. The
  state sequence is summed out EXACTLY by the forward algorithm,
  evaluated on the same AD tape as everything else, so nothing about the
  Laplace machinery changes: a random effect in a state’s linear
  predictor is integrated outside the exact state sum. Each of the
  wrapped family’s parameters is copied per state and suffixed (`mu1`,
  `sigma1`, `mu2`, …) with the full formula grammar, random effects
  included, exactly as
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  does. Transition probabilities are a row-wise multinomial logit named
  `tr{i}{j}` with state 1 the reference cell of every row, and
  `trans = ~x` gives every cell a predictor at once. Validated against
  `depmixS4` (logLik to 1e-8 or better and every parameter to five
  decimals, on gaussian, poisson and categorical emissions with and
  without transition covariates) and against `hmmTMB` (1e-12 on the
  stationary fixed-effect model).
- `init = "stationary"` (the default) solves the chain’s stationary
  distribution on the tape and costs no parameters; `"estimated"` adds
  `K - 1` free logits, `"uniform"` fixes them. Stationary is refused
  when a transition cell carries a predictor, because a time-varying
  chain has no single stationary distribution.
- [`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md)
  returns the smoothed state occupancies `P(S_t = k | y)` from a
  forward-backward pass - the
  [`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
  analog - and
  [`hmm_viterbi()`](https://aforren1.github.io/frmtmb/reference/hmm_viterbi.md)
  the maximum-a-posteriori state path.
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  `predict(type = "response")` and the response and pearson residuals
  all route through
  [`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md),
  so they report the occupancy-weighted mean
  `sum_k P(S_t = k | y) mu_k(x_t)` rather than any single state’s.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) walks the chain
  forward per sequence and then emits, so DHARMa and
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
  see the fitted persistence.
- A missing response is a time point the chain passes through without
  emitting:
  [`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) keeps
  the row and masks its emission instead of letting `na.action` drop it,
  which would shorten the chain and bias the transition matrix.
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) counts every row;
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) are `NA` there
  while
  [`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md)
  is not.

### Refusals

- `REML`, `quadrature = TRUE` and `frmtmb_control(profile = TRUE)` are
  refused on an
  [`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) fit,
  each with the reason. REML in particular used to run and produce a
  partial restricted likelihood matching no standard definition.
- [`weights()`](https://rdrr.io/r/stats/weights.html), `cens()`,
  [`trunc()`](https://rdrr.io/r/base/Round.html), `se()` and `mi()` on
  the response, multivariate models and `rescor`,
  `residuals(type = "osa")`, `residuals(type = "deviance")`,
  `predict(se.fit = TRUE)` and `predict(newdata =)` on the response
  scale, and
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  are refused, all naming the per-sequence likelihood as the reason.
- A grouping in which every sequence has length 1 is refused: the
  transition parameters are then flat directions the reported `df` would
  still count, and the model is a
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md).
  Holding every transition dpar at a constant lifts the refusal, and the
  fit then reproduces
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  to 1e-12.
- A start with every state’s location predictor at the same value warns:
  it sits on the label-symmetry axis, where the optimizer cannot
  separate the states. The default quantile-spread starts never do this.

### Known limits

- Multimodality is real and no convergence diagnostic reports it. On a
  random-effect model the default cold start has been measured
  converging 8.1 log-likelihood units below the global optimum with
  `convergence == 0`, `max|grad| == 3.5e-4` and a positive-definite
  Hessian. Compare starts before reporting.
- With random effects the Laplace approximation is genuinely approximate
  even for a gaussian response, because the integrand is a mixture over
  state sequences: the measured bias is -0.126 in the log-likelihood
  (8.9e-5 relative) and 4.4e-4 in the parameters against adaptive
  quadrature.
  [`?hmm`](https://aforren1.github.io/frmtmb/reference/hmm.md) says so.
- The tape build grows slightly faster than linearly in the number of
  rows: about 1.9 s at 20 000 rows, against milliseconds to evaluate it.
  Below 5 000 rows nothing is noticeable.

### New

- Latent class analysis, the poLCA measurement model, as the family
  `lca(K)`. The response is a matrix of polytomous item codes, one row
  per subject and one column per item; the items are conditionally
  independent given a subject’s latent class, and each class carries its
  own item-response profile. Class membership is the
  `theta1 ... theta{K-1}` dpars with full linear predictors, so a
  covariate on the model formula gives poLCA’s latent class regression
  for free, with
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md),
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  and
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  on the gating coefficients. Items may have different numbers of
  categories; `ncat` declares them, and by default they are inferred as
  the largest observed code per item, as poLCA does. The item profiles
  are family extra parameters, one vector `pi<j>` per item holding its
  `K * (C_j - 1)` reference-category logits, so they appear per item in
  [`summary()`](https://rdrr.io/r/base/summary.html) and as `pi<j>_<i>`
  in [`confint()`](https://rdrr.io/r/stats/confint.html). Validated
  against poLCA on its own shipped data: the carcinoma 3-class model
  agrees to 4.2e-8 in log-likelihood, 2.8e-8 in item profiles, 6.2e-9 in
  class sizes and 9.9e-8 in posterior membership, and the election
  latent class regression agrees to 1.1e-7 in log-likelihood, 1.1e-7 in
  item profiles and 9.9e-7 in gating coefficients. A hand-rolled
  [`optim()`](https://rdrr.io/r/stats/optim.html) reference on simulated
  data agrees to 1.5e-9, and a one-item fit reaches the saturated
  single-categorical likelihood to 5.6e-11.
- [`lca_profiles()`](https://aforren1.github.io/frmtmb/reference/lca_profiles.md)
  returns the class-conditional item-response probability tables
  (poLCA’s `probs`) with the estimated class sizes attached, and prints
  them.
  [`lca_probs()`](https://aforren1.github.io/frmtmb/reference/lca_probs.md)
  returns posterior class membership per subject (poLCA’s `posterior`)
  with the relative entropy of the classification attached; it is
  [`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
  under an LCA-specific name and check.
- `lca(na.rm = FALSE)` keeps subjects with missing items and masks each
  missing item out of that subject’s likelihood, poLCA’s `na.rm = FALSE`
  behavior. The default drops incomplete subjects through the usual
  `na.action`, which is poLCA’s default.

### Notes

- [`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md) starting
  values are deterministic: subjects are scored by the mean of their
  item codes rescaled to `[0, 1]`, cut into `K` equal-count slices, and
  each slice’s smoothed category proportions seed one class. Class 1 is
  the low-score end, so a data set always gets the same labeling.
  Multimodality is unchanged;
  [`?lca`](https://aforren1.github.io/frmtmb/reference/lca.md) shows the
  perturbed-`start` loop that replaces poLCA’s `nrep`.
- [`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md) v1
  refuses random effects, smooths and `gp()` anywhere in the model (that
  is the growth-mixture shape, which `mixture(..., groups = ~g)` fits),
  `REML`, `profile = TRUE`, `quadrature`, every addition term,
  [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md), and
  `residuals(type = "osa")`.
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and
  `predict(type = "response")` are refused because a matrix of nominal
  item codes has no mean;
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns the gating
  linear predictor.

## frmtmb 0.34.0

Within-group residual correlation (R-side effects), quantile regression
completions, and plot conveniences.

### New

- Within-group residual correlation, brms’s R-side autocorrelation, as
  formula terms [`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`,
  `arma()`, `cosy()` and `unstr()` for
  [`gaussian()`](https://rdrr.io/r/stats/family.html) and
  [`student()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md).
  The residuals of a group become one multivariate draw,
  `y_g ~ N(mu_g, D R D)`, with a `sigma ~ x` distributional model
  entering through the diagonal. The ARMA autocorrelation function is
  exact and orders above one are supported (brms caps its covariance
  form at one). Validated against
  [`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html) and
  [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) under ML and REML
  to 1e-9 or better across all five structures, on balanced and ragged
  data, and with random effects alongside the correlated residual. New
  [`autocor_matrix()`](https://aforren1.github.io/frmtmb/reference/autocor_matrix.md)
  returns the fitted correlation matrix; the parameters appear in
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
  and
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  under brms’s names. See `?frmtmb-autocor` for the deliberate
  divergences: `sigma` is the MARGINAL residual SD (brms uses the
  innovation SD for ar/ma/arma; the migration vignette gives the
  conversion), lags count gaps in the global time-level set (nlme
  semantics), and brms’s default `cov = FALSE` residual-regression form
  is refused rather than reinterpreted.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws correlated
  residuals, so DHARMa and
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
  see the fitted autocorrelation. Combinations that stop the likelihood
  factorizing over rows (weights, cens/trunc, se, mi, rescor, mixtures,
  quadrature, nl, OSA residuals, other families) are refused with the
  alternative named.
- [`zero_inflated_asym_laplace()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  (brms spelling): a point mass at zero mixed with the asymmetric
  Laplace, for zero-inflated quantile regression. A new “Quantile
  regression inference” section on the families page states plainly that
  ALD-based Wald intervals are not calibrated (a property shared with
  brms) and points at
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md).
- `plot(conditional_effects(fit), points = TRUE)` overlays the raw
  observations, the brms argument, previously ignored. Points are drawn
  under the bands; displays where raw observations are not meaningful
  (per-category ordinal, non-mean dpars, matrix responses) say so
  instead of failing silently.
- Hidden Markov models: a feasibility study (dev/hmm-feasibility.md)
  establishes that the forward algorithm tapes in RTMB, composes with
  the Laplace approximation over random effects, and is expressible
  today through
  [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  with `vint()` payloads, validated against depmixS4 and hmmTMB. No
  user-facing grammar yet; the design for a first-class
  [`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) family
  is recorded.

## frmtmb 0.33.0

Two residue fixes, a documentation dark mode, and CI repairs.

### Behavior changes

- [`posterior_epred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  on an ordinal fit returns a `draws x observations x categories` array
  with `dimnames` `list(NULL, observation names, category levels)`,
  replacing the flattened `draws x (n * K)` matrix of v0.31/v0.32. This
  is brms’s documented shape (an S x N x C array for categorical and
  ordinal models). The array is the old matrix reshaped:
  `matrix(as.vector(ep), nrow = ndraws)` recovers the previous value and
  column order exactly, and `ep[, , "high"]` or `ep[k, , ]` replace the
  old naming recipe. Scalar-response families,
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md),
  and
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  are unchanged. A documentation error is also corrected:
  [`multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  response-scale predictions are a vector, not a category matrix, and
  the docs no longer claim otherwise.

### New

- A `bf(nl = TRUE)` body can name an object of its formula environment
  that could never be a column of `data` (a data.frame, a list, an
  environment, a formula), so `frm_ode(..., events = doses)` takes a
  dosing table by name instead of needing an inline
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html) or a wrapper
  function. The boundary is exactly what
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  refuses: vectors, factors, and matrices are legal model-frame
  variables and still resolve through the model frame, and a column of
  `data` still wins over a same-named object. A body that fails to
  evaluate reports which names were resolved outside the data.
- The documentation site has a light/dark/auto theme toggle.

### CI

- The pkgcheck and dependency workflows resolve RTMBode from r-universe:
  `extra-repositories` on the setup-r action, and an `R_PROFILE_USER`
  profile for the pkgcheck container, whose pak cannot read
  `Additional_repositories` from DESCRIPTION.

## frmtmb 0.32.0

The merged-\|ID\| Kronecker path, the ordinal prediction surface, ODE
dosing, and vignette figures.

### New

- [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
  gains `events`, a per-dataset dosing table (`time`, `value`, `state`,
  optional `group`, `method` of add/replace/ multiply, `duration` for
  constant-rate infusions) supporting repeated doses, and `event_scale`,
  an estimated multiplier on every dose amount, so bioavailability is an
  ordinary nonlinear parameter with covariates and random effects.
  Dosing splits the solve at the event times instead of using deSolve’s
  `events` argument: deSolve events do not differentiate correctly
  through RTMBode’s augmented system (`replace` and `multiply` measured
  42 and 59 percent relative gradient error; reported upstream). A
  branch on time inside a derivative function is refused by RTMB, and an
  [`approxfun()`](https://rdrr.io/r/stats/approxfun.html) forcing table
  inside one is silently wrong; both are documented.
  [`vignette("ode")`](https://aforren1.github.io/frmtmb/articles/ode.md)
  gains a repeated-dosing section.
- Terms sharing an `|ID|` key whose grouping is `gr(g, cov = A)` or
  `gr(g, prec = Q)` now merge into one Kronecker block instead of being
  refused:
  `mvbf(bf(y1 ~ (1 | q | gr(id, cov = A))), bf(y2 ~ (1 | q | gr(id, cov = A))))`
  is the same joint density as the long format
  `(0 + trait | gr(id, cov = A))`, verified to 7e-10 on the log
  likelihood under ML and REML and on both the covariance and precision
  sides. The compatibility rows move from refused to conditional.
- Vignettes carry figures where a picture beats a table: a forest plot
  for the meta-analysis, the monotonic step shape, the location-scale
  band with the mgcv overlay, growth-mixture trajectories by recovered
  class, the sleepstudy small multiple, and concentration-time curves in
  the ODE vignette. Figures use the suggested tinyplot package and are
  skipped when it is not installed.

### Behavior changes

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on an ordinal
  family returns the n x K matrix of category probabilities (the brms
  convention), exactly as `predict(type = "response")` does; it
  previously returned the latent linear predictor. The invariant
  `predict(type = "response") == fitted()` now holds for every family
  with no documented exception. `predict(type = "link")` keeps the
  latent predictor.
- `residuals(type = "response")` and `type = "pearson"` on an ordinal
  fit score the categories by their codes: `y - sum_k k * P(y = k)`,
  standardized by that distribution’s own sd. `"response"` previously
  returned the observed code minus the latent predictor without saying
  so, and `"pearson"` errored.
- One `|ID|` label spread over more than one grouping specification is
  now an error. Previously the merge key included the deparsed grouping
  call, so `(1 | q | g)` next to `(1 | q | gr(g, cov = A)))` landed in
  separate blocks and silently did not correlate at all, which is not
  what the shared label asks. `|ID|`-linked `gr()` terms whose matrices
  resolve to different objects are refused at frame assembly, comparing
  the resolved matrices.

### Corrected behavior

- `emmeans()` on an ordinal fit no longer fails with “Non-conformable
  elements in reference grid”; the basis is built by column name
  (marginal means stay on the latent scale, the emmeans convention for
  `clm`-like models).
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
  and
  [`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
  work on ordinal fits (simulated ordered factors are compared on their
  integer codes; DHARMa’s rank transform runs on simulated categories
  with `integerResponse = TRUE`).
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) works on an
  ordinal fit.
- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  on an ordinal fit draws one probability curve per response category (a
  `cats__` column), with delta-method standard errors over the joint
  covariance of coefficients, thresholds, and `cs()` terms, and bands on
  the logit scale so they cannot leave \[0, 1\]. `method = "predict"` is
  refused there; `dpar = "mu"` gives the latent display.

## frmtmb 0.31.0

ODE models, an exotic-models case-study vignette, and the defect wave
the vignette work exposed.

### New

- [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
  fits ordinary differential equation models inside a `bf(nl = TRUE)`
  body. Dynamics parameters and initial states are ordinary nonlinear
  parameters, so they take fixed effects, random effects and covariates,
  and the Laplace approximation is exact through the solver’s adjoint.
  It solves one small system per group and scatters the solution back
  into row order; ragged designs, unsorted rows, repeated times and an
  observation at `t0` all work. Needs the optional RTMBode package (not
  on CRAN; `Additional_repositories` now names
  <https://kaskr.r-universe.dev>), and every code path, test, example
  and vignette chunk degrades cleanly without it. New
  [`vignette("ode")`](https://aforren1.github.io/frmtmb/articles/ode.md)
  works a population pharmacokinetic model on
  [`datasets::Theoph`](https://rdrr.io/r/datasets/Theoph.html).
- [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
  refuses the ways an ODE model goes quietly wrong: a dynamics parameter
  that varies inside a solve group is rejected by name at frame
  assembly; fixed-step integrators such as `rk4` are refused because
  they return a different likelihood, not a noisier one; a system near
  the Laplace ceiling of about eight states warns. A failed solve is
  reported, not absorbed: on numeric evaluation paths
  ([`predict()`](https://rdrr.io/r/stats/predict.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), and a body
  with no estimated dynamics inputs) the group’s rows get a penalty
  value and a warning naming the group, and
  [`frm_ode_failures()`](https://aforren1.github.io/frmtmb/reference/frm_ode_failures.md)
  reads the record back. During fitting most failures cannot be detected
  at all (RTMBode returns an AD object), so they surface as the
  optimizer’s NA/NaN gradient; the help page says which is which.
- New
  [`vignette("case-studies")`](https://aforren1.github.io/frmtmb/articles/case-studies.md):
  eight worked models from the showcase literature of brms, MCMCglmm,
  metafor and mgcv. The animal model and its multi-trait form,
  phylogenetic regression, random-effects meta-analysis and
  meta-regression, monotonic ordinal predictors, location-scale
  regression, growth mixtures, measurement error, and sequential ordinal
  models with category-specific effects. Every section cross-checks its
  fit against a reference package or a closed form, the multi-trait
  section proves its pedigree matters with an identity refit, and
  `tests/testthat/test-case-studies.R` pins the agreements. Suggests
  gains ape and metafor for the cross-checks, and data.table and tibble
  for tabular-input tests.

### Breaking

- A `gr(cov = )` or `gr(prec = )` term whose `|ID|` key is shared with
  another term is refused at parse time. The cross-formula merge
  hardcoded an unstructured covariance and silently discarded the
  matrix, so earlier fits of that construct ignored the structure they
  named. The error points at the supported long format,
  `(0 + trait | gr(id, cov = A))`, which takes the verified Kronecker
  path. A lone `gr()` term with an unshared key keeps its structure and
  still fits.

### Corrected behavior

- [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  and
  [`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
  expose `sd_<group>__<term>` and `cor_` names for `gr(cov = )`,
  `gr(prec = )` and `equalto()` blocks, so heritability-as-ICC is one
  line on an animal model. Smooth, `gp()`, `car()` and `spde()` blocks
  stay excluded;
  [`?hypothesis`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  says which and why.
- `vcov(full = TRUE)` works under REML and `profile = TRUE`, taking the
  outer block from the joint precision instead of warning and returning
  only the fixed effects. Its rows are
  [`confint()`](https://rdrr.io/r/stats/confint.html)’s rows.
- The observation-level-random-effect check no longer fires on
  `gr(cov = )`, `gr(prec = )` or `equalto()` blocks with one row per
  level: a fixed relationship matrix identifies the two variances, which
  is exactly the animal model with one record per individual.
- `predict(type = "response")` on the four ordinal families returns an n
  x K matrix of category probabilities named by the response’s levels,
  honoring `cs()` on newdata; `se.fit` there is refused rather than
  faked. marginaleffects gets the categorical `group` convention (also
  fixing
  [`multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  row alignment), `posterior_linpred(transform = TRUE)` stays on the
  latent dpar as documented, and
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)’s
  flattened categorical shape is named and documented.
- [`anova()`](https://rdrr.io/r/stats/anova.html) and
  [`drop1()`](https://rdrr.io/r/stats/add1.html) report NA instead of
  `< 2.2e-16 ***` when the compared models have equal df.
- [`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
  no longer reports a zero-width interval at its internal clamp for a
  correlation estimated at plus or minus 1 or an sd at zero; those rows
  keep the estimate, get NA bounds, and warn.
  [`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
  and [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md)
  handle two blocks sharing a term label (one printed twice, the other
  silently dropped a block).
- A `Date`, `POSIXct` or `difftime` predictor is reported at frame
  assembly with its unit and 1970-01-01 origin plus the centering fix;
  the epoch-scale magnitude can stop the optimizer converging. Responses
  and grouping variables are exempt (srr G2.5, G2.9).
- A nonlinear formula body that fails to evaluate reports which of its
  names were resolved to functions instead of found in `data`, so a
  misspelled column sharing a name with a base function (`t`, `c`) is
  named rather than surfacing as a coercion error.
- srr completeness: `srr_stats_pre_submit()` reports no missing
  standards (seven were untagged; six were already satisfied and are now
  tagged with pinning tests, G2.9 is satisfied by the datetime reporting
  above).

## frmtmb 0.30.0

data2, ODE feasibility, and the remaining rOpenSci runway.

### New

- `frm(..., data2 = list())`, the brms spelling for structural objects:
  the matrices of `gr(cov =)`/`gr(prec =)`, `equalto()`, `car()` and
  `spde()` resolve from `data2` before the data and the formula
  environment (compound expressions over `data2` objects work, a
  documented permissive divergence from brms). `data2` is stored on the
  fit, so [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) carries the
  matrices and
  [`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md)/[`update()`](https://rdrr.io/r/stats/update.html)/[`drop1()`](https://rdrr.io/r/stats/add1.html)/[`influence()`](https://rdrr.io/r/stats/lm.influence.html)/[`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  re-assemble in a session where the calling environment is gone.
- Population pharmacokinetic models work today through `nl = TRUE`
  bodies calling
  [`RTMBode::ode()`](https://rdrr.io/pkg/RTMBode/man/ode.html) per
  subject (RTMBode installs from kaskr.r-universe.dev): the Laplace
  approximation is exact through the adjoint solver, and a
  one-compartment population fit agrees with nlmixr2’s FOCEi to three
  decimals with no compiler. The feasibility study, the sharp edges
  (never stack subjects into one system; adaptive integrators only), and
  the planned
  [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
  helper are recorded in dev/ode-feasibility.md.

### rOpenSci runway

- Every exported help page has executing examples (35 pages added; the
  as-cran check runs them all in ~14s).
- Every internal function (288) is documented with its contract under
  `@noRd`; a new “Inputs and preprocessing” vignette covers the
  formula-to-design pipeline, accepted predictor classes, argument
  types, attribute survival, terminology, and reproducible runtime
  claims. srr standards: 105 of 116 tagged in code, docs, and tests; the
  other 11 documented as not applicable.
- Nine standards tests added (degenerate inputs, noise invariance,
  exact-fit behavior, row-name retention, runtime scaling).
- CI gained the rOpenSci `pkgcheck` action (the editors’ submission
  checks, srr included) and a coverage workflow.

### Corrected behavior

- Rows dropped by `na.action` are reported with a count (a
  [`message()`](https://rdrr.io/r/base/message.html), so
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) silences
  it); fitting zero rows says so instead of blaming missing values.
- Eleven duplicated error messages across the package were disambiguated
  with call context.
- When a structural expression fails under the `data2` mask (for example
  [`solve()`](https://rdrr.io/r/base/solve.html) of a singular matrix),
  the error now reports that cause instead of a misleading “not found”
  from the fallback lookup.
- An unordered factor response to an ordinal family warns with the level
  order about to be used, then fits (the brms behavior); `mo()` keeps
  its stricter ordered-factor requirement.

## frmtmb 0.29.0

Review fixes for the spatial wave and rOpenSci review preparation.

### Breaking

- `spde()` grouping values are now explicitly MESH ROW INDICES: integer,
  whole-number character, and integer-level factor spellings are
  accepted and land identically; any other labels error stating the
  contract. Previously a non-mesh-order grouping (integer node IDs
  included, which sort lexicographically) fit a silently permuted field
  with a wrong likelihood and no warning. Code that passed non-numeric
  node labels must now map them to row indices.

### Corrected behavior

- `car()`/`spde()` accept call-valued grouping variables
  (`gr = factor(node)`).
- The Windows sequential-chains guard in
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  also covers `options(mc.cores =)`, which rstan reads as its default.
- Optimizer error wrapping no longer claims a numerical cause for
  non-numerical errors; the underlying message leads.
- The `escar` log-determinant is finite at `rho` near 1 (the adjacency
  eigenspectrum is clamped against floating-point overshoot).
- The extreme-parameter heuristics in
  [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
  and
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  only test log-SD components and name the parameter; a `car`/`bym2`
  mixing proportion at its boundary no longer trips a misleading
  singularity warning.
- `car()` adjacency validation: NA entries and negative weights error
  informatively; one-sided dimnames work (the rownames-only path was
  broken); duplicate names error.
- Non-finite covariance warnings fire once per fit instead of once per
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) call.

### rOpenSci preparation

- Review infrastructure: CONTRIBUTING, code of conduct, CITATION,
  codemeta, `Language` field, return-value documentation, README
  currency fixes.
- srr statistical-software standards: General and Regression categories
  scaffolded with 77 standards tagged at their code locations, 11 marked
  inapplicable with justifications, and the 28 remaining gaps recorded
  as the submission roadmap in the standards file. Test coverage
  measured at 87.6% (full suite).

## frmtmb 0.28.0

Spatial GMRF grammar, the mclust covariance taxonomy, and the close-out
of the fuzz tier.

### New

- Spatial conditional autoregressions with the brms spelling:
  `car(M, gr = g, type = "escar"/"esicar"/"icar"/"bym2")` over an
  adjacency matrix, with analytic log-determinants throughout (no
  on-tape factorization) and brms’s per-component soft sum-to-zero
  constraint, whose distance to the hard-constrained ML is measured and
  documented (4e-4 logLik at the default, shrinking quadratically).
  Validated against hand-rolled direct ML to 1e-11.
- `spde(fem, gr = node)`: an SPDE-Matern field from fmesher/INLA FEM
  matrices, Q(kappa, tau) assembled on the tape from fixed sparse
  matrices; sd and range reported through the planar identities.
- `gr(g, prec = Q)` supports correlated slopes (sparse Kronecker
  precisions; agrees with the dense `gr(cov = solve(Q))` equivalent to 0
  ulp).
- `mixture_mvn(K, D, model =)` gains the mclust covariance taxonomy
  (EII, VII, EEI, VEI, EVI, VVI, EEE, VVV) with per-model parameter
  templates; matches
  [`mclust::Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html)
  to 1e-12 log-likelihood and 1e-15 posteriors on intercept-only fits,
  while keeping covariate-dependent means, which mclust cannot fit.
- `anova(refit = TRUE)` refits REML models with ML for comparison
  (warm-started, with a message naming what was refit); `getME()` on the
  lme4 generic with the commonly consumed vocabulary.

### Corrected behavior

- Non-finite covariance matrices warn at
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  and at fit time under `se = TRUE`, pointing at
  [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
  (a positive-definite Hessian can still be numerically singular to
  invert; previously silent NaNs).
- Optimizer failures recover automatically where a better start provably
  exists (the plain-Laplace optimum under `profile = TRUE`; a cold start
  when a warm start is itself the failure), and every optimizer error
  carries the model label and remedies instead of a bare RTMB message.
- One quadrature configuration is rescued by displaced calibrations; the
  six irrecoverable thin-data configurations refuse with messages naming
  the configuration and the guaranteed `quadrature = FALSE` fallback.
  The fuzz tier runs green.
- `frm_sample(cores > 1)` on Windows falls back to sequential chains
  with a warning instead of failing incomprehensibly (RTMB tapes and
  objective closures cannot reach socket workers; upstream
  kaskr/tmbstan#27).
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
  and
  [`as_draws()`](https://aforren1.github.io/frmtmb/reference/as_draws.md)
  dispatch without attaching bayesplot or posterior (frmtmb ships its
  own generics, dual-registered).

### Notes

- dev/benchmarks.md records the RTMBp (parallel RTMB) measurement: no
  gain on Cholesky-dominated mixed models, 1.7-2x on
  accumulation-dominated GLMMs at 4-16 threads; not adopted.
- dev/feature-gaps.md records the ODE roadmap (RTMBode) and the
  Suggests + Additional_repositories packaging plan for non-CRAN
  extensions.

## frmtmb 0.27.0

Fixes from a full review of the v0.22-v0.25 waves.

### Corrected behavior

- `cens()` combined with [`trunc()`](https://rdrr.io/r/base/Round.html)
  produced a silently wrong likelihood: the censoring contribution was
  not restricted to the truncation window (right-censoring now
  contributes (F(ub) - F(y))/Z), which inflated the dispersion by ~14%
  on the reference problem. The corrected objective matches a
  hand-rolled windowed likelihood to 1e-8, collapses exactly to the old
  form without truncation, and the cens+trunc OSA residuals recalibrate
  to the analytic PIT (5e-14).
- `residuals(type = "deviance")` under `se()` now weights each row’s
  unit deviance by its own variance (the glm prior-weight form) instead
  of treating a common dispersion as shared.
- `predict(se.fit = TRUE)` on quadrature fits no longer dies with a
  non-conformable error: it reports mode-conditional standard errors
  with a warning. Singular joint precisions in the predict path now
  degrade like [`vcov()`](https://rdrr.io/r/stats/vcov.html) does, and
  the shared solver uses
  [`Matrix::solve`](https://rdrr.io/pkg/Matrix/man/solve-methods.html)
  (strictly more robust for ill-conditioned but invertible GLMM
  precisions).
- [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  chain inits are clamped inside any bounds passed to Stan, with a
  warning naming the parameters when the ML mode itself violates a bound
  (previously an incomprehensible rstan error).
- `vint()`/`vreal()` variables are required in `newdata`; a custom
  family’s prediction previously returned a length-0 vector with no
  message when they were missing.
- Character censoring codes like “0”/“1”/“-1” decode (the error message
  had advertised them).

### Compatibility registry and fuzzer

- The registry’s precedence is redesigned (lexicographic comparison of
  sorted side specificities, with a validator that forbids file-order
  ties except documented overrides). 488 of 3750 resolutions were
  corrected against probed reality, including 429 pairs whose covstruct
  conditions a family-level rule had erased, and the newly recorded
  caveat that spatial structures over a plain factor silently use level
  order as coordinates.
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
  accepts vectors and errors on empty input.
- The fuzz tier’s invariants were strengthened (a previously unfailable
  interval check replaced by the Wald identity plus parameter-coverage
  tallies; per-row simulate agreement; grammar divergences asserted;
  convergence demotion keyed to actual convergence verdicts after fixing
  a truncation that discarded them). The strengthened tier reproduces
  only known findings.

## frmtmb 0.26.0

Pooled model comparison across imputations and the diagnostics/UX
backlog.

### New

- [`anova()`](https://rdrr.io/r/stats/anova.html) on
  [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  fits pools nested-model tests across imputations: D1 (multivariate
  Wald), D2 (chi-square combining), and D3 (Meng-Rubin likelihood
  pooling, with the plug-in leg evaluated by re-taping each imputation’s
  objective at the pooled parameters - no refits). D1/D2 match `mice` to
  1e-7 on an exactly-shared reference; D3 is validated against the
  Meng-Rubin formula directly, since
  [`mice::D3`](https://amices.org/mice/reference/D3.html)’s `fix.coef`
  variant is not the plug-in statistic. Includes an ARIV clamp and a
  Reiter-df fallback for a small-m case where mice returns NaN.
- `cbind(successes, failures)` binomial responses are accepted (the
  glm/lme4/glmmTMB spelling), rewritten internally to
  `successes | trials(successes + failures)`; bit-identical to the
  `trials()` form.
- [`simulate()`](https://rdrr.io/r/stats/simulate.html) returns ordinal
  draws as ordered factors with the original levels and multinomial
  draws as count matrices (both families previously had no simulator),
  and respects `na.exclude` padding.
- `frmtmb_control(check_nlev_1 =, check_olre =)`: lme4-style
  warning/ignore/stop vocabulary for one-level grouping factors and
  gaussian observation-level random effects (the `se()`-based
  meta-analysis idiom is recognized and not flagged).
- [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
  adds a complete-separation heuristic, predictor-scale warnings
  pointing at `autoscale`, and an isSingular-style verdict independent
  of the Hessian; it also no longer errors on fits without random
  effects.

### Corrected behavior

- Models with zero free outer parameters fit degenerately instead of
  dying inside nlminb (fixing latent empty-sdreport and empty-gradient
  bugs found along the way).
- A nonlinear-parameter name colliding with a data column errors instead
  of silently shadowing the column; nonlinear fits that fail from
  default zero starts name `start =` in the error.
- REML [`anova()`](https://rdrr.io/r/stats/anova.html) compares fits
  whose fixed-effect designs span the same column space (term reordering
  included) instead of refusing every REML pair; genuinely different
  designs and REML/ML mixes still refuse with the reason.
- Offsets in distributional-parameter formulas were verified correct (to
  1e-13) and are now regression-tested against the silent-drop failure
  mode reported upstream (glmmTMB#625).

## frmtmb 0.25.0

Simulation-workflow ergonomics, the last deferred method-surface items,
and CRAN readiness.

### New

- [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
  accepts natural-scale parameter names - the same vocabulary as
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)/[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md):
  coefficients by name, `sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`,
  response-scale `sigma` - inverted to the internal parameterization
  through the covariance registry (us/diag/homdiag/smooth/gr structures;
  others refuse by block name). The internal spelling still works and is
  byte-identical.
- Prior-predictive simulation, the `sample_prior = "only"` analog:
  `frm_simulate(..., priors = set_prior(...))` draws a fresh parameter
  vector per replicate on each prior class’s documented scale (natural
  SDs for class `sd`, truncation by rejection) and attaches the drawn
  parameters for prior-predictive checks. Every coefficient and SD must
  be pinned by a prior or `newparams`; omissions error instead of
  silently becoming zero effects.
- `simulate(censored = TRUE)` applies the fitted type-I censoring
  mechanism to the draws; the default remains the latent uncensored
  response, which is also brms’s `posterior_predict` convention
  (verified against brms source) and is now documented.
- `residuals(type = "deviance")` across the GLM family set (gaussian,
  poisson, binomial/bernoulli, Gamma, exponential, inverse.gaussian,
  nbinom1/2, geometric, beta, tweedie), exact against
  [`stats::glm`](https://rdrr.io/r/stats/glm.html) where glm offers
  them; families without a standard unit deviance refuse by name.
  [`deviance()`](https://rdrr.io/r/stats/deviance.html) itself stays
  `-2 logLik` (lme4 convention).
- `predict(type = "response", se.fit = TRUE)` now works for
  non-identity-mean families (zi/hurdle, lognormal, trials-binomial,
  truncated responses) through a joint delta method across all dpar
  linear predictors, including cross-dpar covariance (agrees with
  glmmTMB’s response-scale standard errors to 1.2e-5 on a zero-inflated
  poisson mixed model).

### Corrected behavior

- `conditional_effects(method = "predict")` evaluates addition terms on
  the effect grid, so its bands respect `trials()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) (binomial bands are
  counts in \[0, n\], not Bernoulli 0/1); aterm variables must be pinned
  in `conditions`, and the error names the variable.

### CRAN and infrastructure

- Heavy reference-validation test files are `skip_on_cran()`-gated: the
  CRAN-condition suite drops from ~128s to ~72s on the development
  machine while CI (NOT_CRAN=true) keeps full coverage.
- nlme added to Suggests (a test uses
  [`nlme::Soybean`](https://rdrr.io/pkg/nlme/man/Soybean.html) as
  reference data; CI’s `--as-cran` unstated-dependency check halts
  without the declaration).
- Benchmark verdict recorded in dev/benchmarks.md: optimParallel is not
  adopted - with exact AD gradients its concurrency caps at two
  evaluations (measured 1.03-1.22x end to end, slower with cold
  clusters), RTMB tapes cannot ship to PSOCK workers, and 100% of
  InstEval’s optimization time is inside the taped objective and the
  inner sparse Cholesky. `frmtmb_control(profile = TRUE)` remains the
  measured lever for many-coefficient models (1.6x there).

## frmtmb 0.24.0

The quadrature and OSA defect clusters surfaced by the fuzzer and
compatibility registry.

### Corrected behavior

- `quadrature = TRUE` is rebuilt around one root cause: TMB’s
  Gauss-Kronrod marginalization calibrates each integrand’s rescaling
  once, at the parameter values in hand when the tape is built, and
  frmtmb taped at the cold start. A plain Laplace fit now runs first and
  the quadrature tape is built at its optimum. This fixes: conditional
  modes returned as NA
  ([`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md)/[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`predict()`](https://rdrr.io/r/stats/predict.html)
  were silently NA for all groups but the first, whose slot held a wrong
  value) - modes now come from the Laplace inner solve and match
  `glmer(nAGQ = 25)`’s
  [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md) to
  3e-5; and bare “NA/NaN gradient evaluation” crashes for
  poisson/Gamma/Beta with nested or even single scalar blocks - all now
  fit with gradients \< 3e-4.
- `quadrature` combined with
  [`trunc()`](https://rdrr.io/r/base/Round.html) produced logLik = +Inf
  as a successful fit; the combination is refused (the CDF difference
  underflows at quadrature nodes; a stable fix needs log-CDF forms).
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  with `REML` or `profile = TRUE` is refused: both Laplace-expand the mu
  coefficients around a single mode, and a mixture likelihood is
  permutation-multimodal in exactly those coefficients.
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  with quadrature remains supported.
- `residuals(type = "osa")` on censored fits was broken (raw LAPACK
  singularity or NaN for every censored row). Uncensored rows now get
  calibrated one-step residuals by conditioning on the censoring window
  and renormalizing the one-step CDF to it (matches the analytic
  conditional PIT to 7e-9); censored rows are NA (an event has no
  one-step CDF), and row-varying censoring points refuse. Note
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws the latent
  uncensored response under `cens()`, so DHARMa is not a substitute
  there (documented).
- `residuals(type = "osa")` on ordinal fits crashed inside the OSA
  machinery; the ordinal log-densities now handle `oneStepPredict`’s
  taped-observation objects (exact Lagrange-basis category indicator),
  and all four ordinal families produce calibrated residuals (cumulative
  matches the analytic randomized-quantile residual to 4e-14).
- A singular joint precision under REML or `profile = TRUE` no longer
  throws a raw LAPACK error from
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`confint()`](https://rdrr.io/r/stats/confint.html)/
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md):
  it degrades to NaN standard errors with one warning pointing at
  [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md),
  matching the ML branch.
- The brms-migration vignette documents that
  [`binomial()`](https://rdrr.io/r/stats/family.html) without `trials()`
  means Bernoulli here (glm convention) where brms requires `trials()`
  or
  [`bernoulli()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md).

The fuzz tier’s open findings drop from 28 to 16 on the identical plan;
most of the remainder are now informative refusals on thin-data edge
cases rather than defects.

## frmtmb 0.23.0

Defect wave driven by an open-issue sweep of brms/lme4/glmmTMB, plus a
feature-compatibility registry.

### Corrected behavior

- Truncation now reaches the whole post-fit surface, not just the
  likelihood: [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  `predict(type = "response")`, and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) report the
  truncated mean E\[Y \| lb \<= Y \<= ub\] (closed forms for all six CDF
  families, validated to 1e-15);
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md),
  and
  [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
  rejection-sample within bounds; newdata re-evaluates variable bounds.
  `residuals(type = "osa")` on truncated (and untruncated gaussian)
  responses was also miscalibrated and now integrates over the truncated
  support (KS uniformity restored from p ~ 1e-15 to 0.76-0.97).
  Dpar-scale predictions stay untruncated by design. DHARMa/pp_check on
  truncated models are meaningful again.
- `(f || g)` with a factor now yields independent per-level effects (the
  diag structure the syntax promises) instead of a silently fully
  correlated block; identical to an explicit `diag(f | g)`. Numeric
  double-bars are unchanged. Existing factor-double-bar fits change,
  because the old ones were wrong (lme4 has the same open bug,
  lme4#818).
- `ar1()`/`hetar1()` warn when the ordering factor’s integer levels have
  gaps: levels correlate by position (the glmmTMB reading, unchanged),
  so a gap counts as one step; the warning points at `ou()` over
  [`num_factor()`](https://aforren1.github.io/frmtmb/reference/num_factor.md)
  for irregular spacing (glmmTMB#1278).
- Predictions at non-estimable points of a rank-deficient design return
  NA with one warning naming the aliased columns, instead of silently
  returning the partial sum (predict.lm semantics; lme4#303). Collinear
  restatements of kept columns stay exact.
- Grouping factors written as calls work: `(1 | factor(x))`,
  `(1 | interaction(a, b))` (lme4#464).
- A random-effect term crossed with `*` or `:` errors instead of being
  silently refit additively (lme4#196); `mo()`/`mi()` interaction
  multipliers must be numeric (brms#1828);
  [`anova()`](https://rdrr.io/r/stats/anova.html) requires equal `nobs`
  (lme4#622).
- [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  chain initialization: chain 1 anchors at the ML mode, further chains
  are jittered on the unconstrained scale (`init_jitter`), restoring the
  overdispersion Rhat needs; a boundary-mode warning fires for singular
  fits; mixture posteriors are documented as needing `init = "random"`.

### New

- Feature compatibility registry:
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
  answers what plays with what (works / conditional / refused / broken /
  untested) across 3750 feature pairs, and the new “Feature
  compatibility” article is generated from the registry at build time so
  it cannot drift. Known-broken pairs are listed there; the registry
  probing itself surfaced trunc x quadrature and cens x OSA as broken
  (queued for the next fix wave).
- An audit of 559 currently open brms/lme4/glmmTMB issues is recorded in
  dev/test-backlog.md, including the pathologies frmtmb is structurally
  immune to.
- A pairwise grammar fuzzer (env-gated: `FRMTMB_FUZZ=true`): a 310-spec
  covering array over the grammar with metamorphic invariants
  (predict/fitted identity, permutation invariance, simulator support
  membership, vcov sanity) and a brms `make_standata` structural oracle.
  Its first run surfaced the quadrature defect cluster now queued for
  fixing (conditional modes left NA, crashes with nested groups and
  Beta, +Inf logLik with trunc()); the tier reports those as failures
  until they are fixed.

## frmtmb 0.22.0

Cross-validation against brms itself, closer brms compatibility, and fit
ergonomics.

### brms agreement suite

- New test tier validating the elaborate grammar structurally against
  [`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html)/`brmsterms()`
  without any Stan compilation (designs, RE structures, `|ID|` merging,
  ordinal thresholds, `cs()`, `mo()` codes, gp bases, smooths, addition
  terms, mvbf, nl, mixture naming; most exact to 1e-10 or better), plus
  an opt-in numeric tier (`FRMTMB_BRMS_FIT_TESTS=true`) showing
  fixed-effects estimates match the mode of brms’s own generated Stan
  program to 1e-4 and that our `mo()` parameterization is brms’s
  likelihood exactly (via
  [`rstan::log_prob`](https://mc-stan.org/rstan/reference/stanfit-method-logprob.html)).
  brms is in Suggests only.
- Found by the suite and fixed: `cs()` predictor variables never reached
  the combined model frame, so `bf(y ~ x + cs(z))` errored unless `z`
  already appeared elsewhere.

### brms compatibility

- Hilbert-space `gp(x, k =)` now uses brms’s input convention exactly:
  coordinates are rescaled by the largest pairwise distance (one shared
  factor across dimensions) and centered, so the boundary is `L = c` and
  the same `gp()` call is the same approximation in both packages (basis
  agreement with brms to 1e-16, tied coordinates and vector-valued `c =`
  included). This roughly doubles the effective boundary at a given `c`:
  accuracy against the exact gp improves about 25x at the same k (k = 40
  now within 2e-4 logLik), and the default `c = 1.25` is the right
  choice in 2-D as well. Reported `range(gp)` values stay in data units
  (exact log-shift back-transform); fitted lengthscales differ from
  pre-0.22 fits by the data-dependent scale factor.
- `cens()` accepts brms’s character and factor codes
  (`"left"`/`"none"`/`"right"`/`"interval"`, prefix-matched,
  case-insensitively); unknown labels and out-of-range numeric codes
  error informatively instead of coercing to NA.
- `gp()`’s `k`/`c`/`iso`, `rr()`’s `d`, and `se()`’s `sigma` arguments
  now evaluate in the formula environment (brms behavior), so variables
  and expressions work; invalid values error with the offending
  expression named.

### Fit ergonomics

- `frmtmb_control(verbose =)` (and the `frm(verbose =)` shortcut): level
  1 prints timed stage lines (parse, frame, tape, optimize, restarts,
  sdreport) through [`message()`](https://rdrr.io/r/base/message.html)
  so a slow fit shows where it is slow; level 2 adds the optimizer’s own
  iteration trace. Bootstrap, influence, and allfit refit loops stay
  quiet.
- The large-gradient warning now points to
  [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
  and the new “Convergence problems” remedy ladder in the diagnostics
  vignette (scaling, restarts, optimizer comparison, profiles, boundary
  fits, and judging marginal gradients).
- [`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)’s
  nloptr optimizer produced an empty row: nloptr rejects callbacks that
  declare `...` (RTMB’s `obj$fn`/`obj$gr` do), and the error was
  swallowed; separately, missing `ftol_rel` made NLopt report failure at
  the converged optimum. Both fixed; the agreement test now requires
  every optimizer to converge and match.
- The brms-migration vignette documents how to recover the model
  function (the `stancode()` analog): the objective closure via
  `build_objective(fit$frame)`, its joint-vs-marginal relationship to
  `fit$obj$fn`, and how to reproduce `fit$obj`.

## frmtmb 0.21.0

Method-surface audit against lme4/glmmTMB/brms, plus fixes from a full
code review.

### Corrected behavior

- `predict(type = "response")` now returns the expected response for
  every family (equal to
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) in-sample),
  not the primary dpar’s natural scale; zi/hurdle, lognormal, and
  trials-binomial fits were affected. The glmmTMB aliases
  `"conditional"`, `"zprob"`, `"zlink"`, and `"disp"` are accepted
  (validated against glmmTMB on a zero-inflated Poisson to ~8e-7). Note
  `"disp"` returns `sigma` on its natural scale, where glmmTMB returns
  the variance-scale dispersion for gaussian. `se.fit` for the mean of a
  non-identity-mean family is not yet available and errors with
  guidance.
- `rescor = TRUE` now refuses `cens()`,
  [`trunc()`](https://rdrr.io/r/base/Round.html), and `se()` addition
  terms, which the joint-gaussian likelihood silently ignored (wrong
  likelihood with no warning).
- Bounds (`lower`/`upper`, `set_prior` lb/ub) under
  `frmtmb_control(profile = TRUE)` were positionally misaligned: a bound
  on a fixed coefficient could silently pin a covariance parameter
  instead. Bounding a profiled coefficient now errors;
  covariance-parameter bounds land on the right slot.
- [`confint()`](https://rdrr.io/r/stats/confint.html),
  `vcov(full = TRUE)`, and
  [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
  labels were broken or misaligned on `mi()` fits: the latent imputation
  component was counted as an outer parameter.
- `frm_sample(priors = )` failed on fixed-effects-only fits of
  single-dpar families (a `$` partial-match bug), and
  `frm_sample(laplace = TRUE)` was broken twice: the default init had
  the wrong length, and draw columns were mislabeled (a column named
  `b[1]` actually held a covariance parameter).
- [`influence()`](https://rdrr.io/r/stats/lm.influence.html) tables and
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
  now align with [`vcov()`](https://rdrr.io/r/stats/vcov.html) (fits
  with a constant dpar errored; column names are now the coefficient
  labels).
- [`simulate()`](https://rdrr.io/r/stats/simulate.html) follows the
  stats seed contract: a `"seed"` attribute and a restored RNG state
  instead of clobbering the global stream.
- A covariate literally named `sigma` is no longer shadowed by the
  residual sigma in
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)/[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md).
- In `predict(newdata = )`, an NA in a variable used only in a
  random-effect design now propagates to the prediction instead of being
  silently zeroed (only genuinely new levels predict at the population
  value).

### New methods

- [`drop1()`](https://rdrr.io/r/stats/add1.html) (AIC/LRT,
  marginality-aware; matches
  [`lme4::drop1.merMod`](https://rdrr.io/pkg/lme4/man/drop1.merMod.html)
  to 1e-4),
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
  directly on a fit,
  [`dfbeta()`](https://rdrr.io/r/stats/influence.measures.html)/[`dfbetas()`](https://rdrr.io/r/stats/influence.measures.html)
  on [`influence()`](https://rdrr.io/r/stats/lm.influence.html) results
  (stats sign convention; lme4 returns the negation),
  [`na.action()`](https://rdrr.io/r/stats/na.action.html), and an
  lme4-style “Groups:” line in
  [`summary()`](https://rdrr.io/r/base/summary.html).
- [`confint()`](https://rdrr.io/r/stats/confint.html) accepts lme4’s
  `"Wald"` spelling and gains `method = "boot"` (percentile intervals
  via \[frm_bootstrap()\]).
- `vcov(full = TRUE)` rows carry per-parameter names (glmmTMB
  convention).
- [`predict()`](https://rdrr.io/r/stats/predict.html) accepts the
  lme4/glmmTMB `allow.new.levels` dot spelling.
- On draws objects:
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  and [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md)
  (brms-shaped arrays).
- [`emmeans::recover_data`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.html)
  had been silently missing from NAMESPACE since v0.13 (orphaned export
  tag); emmeans support was broken.

### Internals

- One authoritative outer-parameter map shared by
  [`confint()`](https://rdrr.io/r/stats/confint.html), bounds
  resolution, `vcov(full = TRUE)`, and sampling; `graphics` and
  `grDevices` added to Imports; duplicated formula/family coercion and
  gp position-key helpers unified; REML `logLik` df counts only outer
  parameters (documented; lme4 counts the integrated fixed effects too,
  so REML AICs are not comparable across packages).

## frmtmb 0.20.0

- `mixture_mvn(K, D)`: multivariate gaussian mixture components for
  model-based clustering of an n x D matrix response (the mclust /
  clustTMB use case). Every class mean is a full linear predictor, so
  cluster means may depend on covariates and random effects; mixing
  weights are multinomial-logit dpars with their own formulas. Validated
  against a hand-rolled ML fit to 1.3e-10 and against faithful-data
  cluster recovery. Limitations (documented in
  [`?mixture_mvn`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)):
  class covariances are unstructured (`us`) and covariate-free - no
  mclust-style constrained covariance taxonomy (EII..VEV);
  `cens()`/[`trunc()`](https://rdrr.io/r/base/Round.html),
  [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md), and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) are not
  supported.
- Multi-dimensional Gaussian processes: `gp(x1, x2, ...)` takes up to
  three variables, with a separate lengthscale per dimension by default
  (the brms convention) or one shared lengthscale via `iso = TRUE`, in
  both the exact and the Hilbert-space (`k =`) form. The tensor HSGP
  basis is capped at `k^D <= 1000` columns.
  [`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
  reports one range row per dimension. Validated against a direct 2-D ML
  fit to 3.8e-11.
- Kriging prediction for the exact `gp()`: `predict(newdata = )` at
  positions not seen at fit time now returns the GP conditional mean,
  and `se.fit = TRUE` adds the conditional variance, instead of
  erroring. Kriging mean validated to 1.7e-15 and the standard error to
  the exact decomposition identity. The Hilbert-space form already
  interpolated through its basis.
- `frmtmb_control(sparse_x = TRUE)`: sparse fixed-effect design matrices
  ([`Matrix::sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)),
  the analog of `glmmTMB(sparseX = )`, for models where a many-level
  fixed factor dominates memory. Estimates are identical to the dense
  path (gradient agreement 3e-14; 13.8% lower tape memory in the
  validation model);
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) on the
  fit returns a sparse matrix. Prediction falls back to the dense
  builder for newdata containing NA factor rows, which
  `sparse.model.matrix` would silently zero.
- `frmtmb_control(autoscale = TRUE)`: internal standardization of badly
  scaled continuous predictors (the lme4 \>= 1.1.37 feature): a scaled
  pre-fit is mapped back exactly and warm-starts the reported fit,
  parameter magnitudes feed `nlminb`’s scale hook, and convergence and
  standard-error Hessian checks run in scaled units. Turns 12 scaling
  warnings into 0 on the validation problem with a 0.0 log-likelihood
  difference against a manually standardized reference. Reported results
  are always on the original scale.
- [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  pooling extensions: `$pooled_varcorr` pools random-effect SDs and
  correlations across imputations on transformed scales (log for SDs and
  GP ranges, Fisher z for correlations) with Barnard-Rubin degrees of
  freedom, and
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  on a `frmtmb_multiple` object pools delta-method hypothesis estimates
  by Rubin’s rules with t-based intervals. Fixed-effect pooling matches
  [`mice::pool()`](https://amices.org/mice/reference/pool.html) to
  2.1e-6. Likelihood-ratio pooling across imputations (mice’s D3) is not
  implemented; [`anova()`](https://rdrr.io/r/stats/anova.html) on pooled
  fits remains per-fit.

## frmtmb 0.19.0

- Latent classes now combine with continuous random effects
  (growth-mixture models): `mixture(..., groups = ~g)` accepts random
  effects, smooths, and gp() terms in the component formulas. The class
  sum is computed conditional on the latent effects, so one Laplace
  approximation integrates them - the sum-and-integral swap makes this
  exact as an integrand identity, and even crossed random effects are
  structurally fine. Random effects written in a component formula are
  class-specific by construction. Accuracy: the Laplace approximation of
  the class-mixture integrand carries a small documented bias (about 0.1
  log-likelihood units in the validation problem, with parameter
  estimates matching the exact closed-form marginal to 0.01);
  `quadrature = TRUE` is numerically exact when the per-group integrand
  is univariate (validated to 5e-3) and approximate when class-specific
  intercepts couple.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) now draws one
  class per group and then simulates each observation from its group’s
  component;
  [`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
  gives empirical-Bayes classification conditional on the modes.
- [`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
  (brms spelling): lists every parameter name usable in
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  expressions (coefficients, `sd_`/`cor_` summaries, `sigma`); on
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  output it lists the draw columns.
- [`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md)
  (brms spelling): enumerates every slot
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  can target (class/coef/dpar/group rows), from a formula plus data or
  from a fit; the default in every slot is flat.
- README rewritten for the current scope, with related-work pointers
  (glmmTMB, BayesRTMB); getting-started vignette extended to the full
  grammar; SPEC.md status and deviations brought current.

## frmtmb 0.18.0

The last three roadmap features.

- `gp(x)` Gaussian-process terms: exact (a dense squared-exponential
  block over the unique positions, with a standard 1e-6 nugget) or the
  Hilbert-space approximation via `gp(x, k = 30)` (sine basis with
  spectral-density prior SDs), in any linear predictor. Exact fits match
  direct GP marginal ML; the approximation converges to the exact answer
  and predicts at arbitrary positions (the exact form predicts at
  observed positions only).
  [`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
  reports `sd(gp)` and `range(gp)`.
- `mo()` and `mi()` interactions: `mo(x):z`, `mo(x)*z`, `mi(x):z`,
  `mi(x)*z` with numeric multipliers; `mo()` interactions share their
  variable’s simplex (brms convention). Both validated exact against
  direct ML (the `mi()` interaction stays linear in the latent value, so
  the closed-form marginal still applies).
- Group-level latent-class mixtures: `mixture(..., groups = ~g)` sums
  the class assignment per group (growth-mixture / latent-class
  regression; the tractable nested case of brms#1905). Exact against
  direct ML;
  [`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
  returns posterior class probabilities per group, or per observation
  for ordinary mixtures. Restrictions: no random effects or smooths
  alongside, mixing predictors evaluated at each group’s first row, no
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) yet.
  Crossed-design group mixtures remain out of the Laplace class.

## frmtmb 0.17.0

Clearing the roadmap: mixtures, measurement error, and the remaining
grammar gaps.

- `mixture(fam1, fam2, ...)` finite-mixture families: each component
  keeps its own suffixed dpars (`mu1`, `sigma1`, …), mixing proportions
  are multinomial-logit dpars with their own linear predictors
  (mixture-of-experts comes free), and the likelihood is a branch-free
  logsumexp - so random effects and smooths work in any component
  formula. Exact against direct ML for gaussian and poisson mixtures.
  Component means initialize on spread response quantiles; the usual
  finite-mixture multimodality caveats apply
  ([`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md),
  or order intercepts via bounds).
- `mi(sdx)` measurement error (brms `me()`): known per-observation
  measurement SDs make every true value latent, with the observed values
  entering through a measurement model. Exact against the closed-form
  bivariate-normal marginal; attenuation bias is corrected and NAs
  (missing + mismeasured) combine.
- `cs()` category-specific effects for `sratio`, `cratio`, and `acat`
  (refused for `cumulative`, as in brms). Exact against direct ML.
- `equalto(x + 0 | g, V)` covariance structure: fully fixed known V,
  zero parameters (meta-analytic sampling covariances).
- rr() fits now support `se.fit`: the delta method routes the factor
  columns through the loadings Jacobian and adds loading-parameter
  columns; verified parameterization-invariant against the
  full-rank-vs-`us()` equivalence.
- Registered insight methods (`find_formula`, `find_random`,
  `get_parameters`, `get_varcov`, `find_statistic`, link accessors):
  [`insight::is_mixed_model()`](https://easystats.github.io/insight/reference/is_mixed_model.html)
  and the easystats accessor layer now see the random-effect structure.

## frmtmb 0.16.0

The two deferred architecture features.

- `rr(x | g, d = k)` reduced-rank / factor-analytic covariance
  (glmmTMB’s GLVM structure): per-level coefficients are loadings times
  iid standard-normal factors, so `b` holds the factors while the design
  matrices span the coefficient space. Matches glmmTMB exactly and
  reduces to `us()` at full rank. `se.fit` on rr models is deferred (the
  loadings enter the coefficients nonlinearly); `predict`, `simulate`,
  `ranef`, and `VarCorr` (the rank-k covariance) all work.
- `mi()` one-step imputation for continuous predictors, reversing the
  original exclusion: `bf(y ~ mi(x) + z) + bf(x | mi() ~ z)` turns
  missing `x` values into latent parameters that Laplace integrates
  (exactly, in the linear-gaussian case - validated against the
  closed-form marginal likelihood). Imputation uncertainty propagates
  into the coefficient standard errors automatically. Rows with NAs in
  non-`mi()` variables still drop; imputation models must be gaussian or
  student; discrete predictors remain impossible (as in Stan). `me()`
  measurement error is the natural follow-on.

## frmtmb 0.15.0

Tier-2 sweep of dev/feature-gaps.md plus a method surface for sampled
fits.

- `mo()` monotonic effects (brms syntax): scale coefficient times an
  estimated simplex, exact against direct ML; works in prediction,
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md),
  and with ordered factors. No frequentist package offers this.
  Standalone additive terms only for now.
- Sequential ordinal families `sratio`, `cratio`, `acat`, validated
  against direct ML and collapsing to logistic regression at K = 2.
- Covariance structures: `hetar1`, `homcs`, `homtoep` (glmmTMB
  parameterizations, exact where glmmTMB converges) and spatial
  [`exp()`](https://rdrr.io/r/base/Log.html), `gau()`, `mat()` over
  `num_factor(x, y)` coordinates (Matern uses bounded internal
  transforms for stability where both we and glmmTMB otherwise diverge).
- Fixed a latent parser bug: `exp(x)` (or any covariance-structure name)
  used as a plain function in a fixed formula was silently stripped by
  the formula splitter - since v0.1. Such calls are now protected
  automatically.
- `influence(fit, groups = )` and
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html):
  case-deletion diagnostics over warm-started refits.
- [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md):
  fits across multiply-imputed datasets (list or
  [`mice::mids`](https://amices.org/mice/reference/mids.html)) pooled by
  Rubin’s rules with Barnard-Rubin df.
- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  gains `method = "predict"` (prediction intervals) and data-frame
  `conditions` (one condition set per row, brms style).
- `vint()`/`vreal()` addition terms pass arbitrary data vectors to
  custom families (brms custom-family convention). The test suite
  includes a full Wiener drift-diffusion model written as a
  [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  in plain R - the workflow brms needs raw Stan code for.
- Method surface for
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  draws: [`summary()`](https://rdrr.io/r/base/summary.html) (with
  Rhat/ESS),
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
  [`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
  (natural scale),
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  (exact posterior version),
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  /
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  (each draw runs the full prediction machinery),
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md),
  and
  [`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html).
- Examples added across the reference (run by R CMD check).

## frmtmb 0.14.0

Gap-closing milestone from a sweep of the lme4 / glmmTMB / brms
vignettes (dev/feature-gaps.md).

- `se()` addition term: known per-observation sampling SDs
  (meta-analysis), gaussian and student. `se(x)` fixes the residual SD
  (sigma is mapped out); `se(x, sigma = TRUE)` adds an estimated sigma
  in quadrature. With `(1 | obs)` this is random-effects meta-analysis,
  and with `gr(g, cov = A)` the phylogenetic version.
- Binomial (and beta-binomial) responses may be proportions of
  `trials()`, the glm/glmer idiom, converted to counts internally.
- New families: `bernoulli`, `geometric`, `exponential`, `weibull`
  (mean-parameterized, matches survreg’s likelihood),
  `shifted_lognormal`, `hurdle_gamma`, `hurdle_lognormal`,
  `zero_inflated_binomial`, `zero_inflated_beta` (validated against
  glmmTMB), and `asym_laplace` for quantile regression (fixed `quantile`
  dpar reproduces
  [`quantreg::rq`](https://rdrr.io/pkg/quantreg/man/rq.html) estimates).
- `ranef(condVar = TRUE)`: conditional SDs of the modes (TMB/glmmTMB
  convention: fixed-effect uncertainty propagates), plus tidy
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  for `ranef` and `VarCorr` output.
- [`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md):
  refit under every available optimizer (nlminb, optim, bobyqa, NLopt)
  and compare, the lme4 `allFit()` analog.
- [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md):
  de novo simulation from a
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formula
  and supplied parameters with no fitted model (power analysis; the
  glmmTMB `simulate_new()` analog).
- `frmtmb_control(profile = TRUE)` (experimental): profile the
  primary-coefficient vector into the inner Laplace problem, the
  `glmmTMBControl(profile = TRUE)` / `glmer(nAGQ = 0)` speed
  approximation for many-coefficient models; covariance comes from the
  joint precision.
- [`anova()`](https://rdrr.io/r/stats/anova.html) no longer fails on
  models sharing a primary formula (distributional vs plain fits); row
  labels now include dpar formulas.

## frmtmb 0.13.0

Method-surface (“sugar”) milestone: the conventional S3 methods that
downstream packages and user muscle memory dispatch on.

- New accessors: [`sigma()`](https://rdrr.io/r/stats/sigma.html),
  [`terms()`](https://rdrr.io/r/stats/terms.html),
  [`weights()`](https://rdrr.io/r/stats/weights.html),
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
  [`deviance()`](https://rdrr.io/r/stats/deviance.html),
  [`extractAIC()`](https://rdrr.io/r/stats/extractAIC.html) (enables
  [`step()`](https://rdrr.io/r/stats/step.html) /
  [`drop1()`](https://rdrr.io/r/stats/add1.html)),
  [`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.md),
  [`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.md),
  and a [`profile()`](https://rdrr.io/r/stats/profile.html) method
  wrapping
  [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html).
- **Breaking**: [`coef()`](https://rdrr.io/r/stats/coef.html) now
  follows the lme4 / glmmTMB / brms convention (per-group coefficients:
  fixed effects broadcast over levels plus the conditional modes). Use
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) for
  the fixed effects alone. Fits without random effects still return the
  coefficient vector.
- [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  with a base-graphics
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method: grid
  predictions per predictor (or `"x:z"` pair) with Wald bands, smooths
  included, matrix covariates held at column means.
- `plot(fit)`: Pearson-residual diagnostics (residuals vs fitted, normal
  QQ).
- [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
  (registered on the bayesplot generic): predictive checks from
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws through
  any `ppc_*` function.
- [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md):
  tests of arbitrary expressions of the parameters, brms spelling, with
  three methods: delta-method Wald (default), `"profile"`
  (profile-likelihood intervals via
  [`TMB::tmbroot()`](https://rdrr.io/pkg/TMB/man/tmbroot.html) lincombs;
  linear hypotheses, ML fits), and `"boot"` (parametric bootstrap
  percentile intervals; any expression). The expression environment
  includes natural-scale random-effect names (`sd_<group>__<term>`,
  `cor_<group>__<t1>__<t2>`, `sigma`), so ICC inference is
  `hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)", method = "boot")`.
  The result is a `frmtmb_hypothesis` data frame carrying the bootstrap
  draws or profile curves in attributes, with a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  (bootstrap histogram, profile curve, or the implied Wald normal, per
  hypothesis).
- `frm_bootstrap(fit, FUN, nsim)`: parametric bootstrap over
  [`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md)s
  (the `bootMer` analog), with
  [`confint()`](https://rdrr.io/r/stats/confint.html) percentile
  intervals.
- `refit(fit, newresp)`: refit to a new response reusing the assembled
  design with a warm start - one re-tape plus optimization; the engine
  for parametric bootstrap over
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws.
- The accessor set makes insight’s default methods work
  (`find_response()`, `get_data()`, `get_sigma()`, `n_obs()`, …);
  dedicated insight methods (random-effect formula parts,
  `get_parameters()`) are a roadmap item.

## frmtmb 0.12.0

Pre-release consolidation toward CRAN. Highlights by milestone:

### Model grammar (brms-compatible spelling)

- [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) with
  distributional-parameter formulas carrying the full predictor grammar
  (random effects and smooths in any dpar), constant dpars, and addition
  terms [`weights()`](https://rdrr.io/r/stats/weights.html), `trials()`,
  `cens()` (left / right / interval),
  [`trunc()`](https://rdrr.io/r/base/Round.html) (with the inclusive
  discrete correction).
- Nonlinear formulas (`nl = TRUE`) with a full linear predictor per
  parameter.
- Multivariate models:
  [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md) /
  `bf() + bf()` / `mvbind()`, a family per response, gaussian residual
  correlation
  ([`set_rescor()`](https://aforren1.github.io/frmtmb/reference/mvbf.md)),
  and `|ID|` random-effect correlation across formulas.
- mgcv smooths `s()` / `t2()` in any linear predictor, including
  matrix-covariate summation-convention terms: scalar-on-function,
  function-on-scalar, and function-on-function regression validated
  against mgcv.
- Covariance structures: `us`, `diag`, `homdiag`, `cs`, `ar1`, `toep`,
  `ou` (positions via
  [`num_factor()`](https://aforren1.github.io/frmtmb/reference/num_factor.md)),
  and known-structure terms `gr(g, cov = A)` (dense, with correlated
  slopes via a Kronecker product) and `gr(g, prec = Q)` (sparse GMRF).

### Families

gaussian, poisson, binomial, Gamma, lognormal, student, negbinomial,
nbinom1, beta, tweedie, compois, beta_binomial, skew_normal,
inverse.gaussian, exgaussian, zero-inflated and hurdle counts,
cumulative ordinal, matrix-response multinomial, and
[`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
as a plain R log-density with
[`check_custom_family()`](https://aforren1.github.io/frmtmb/reference/check_custom_family.md)
AD verification.

### Estimation and inference

- ML and REML by Laplace approximation; `quadrature = TRUE` upgrades
  scalar random effects to adaptive Gauss-Kronrod marginalization
  (matches `glmer(nAGQ = 25)` and GLMMadaptive).
- MAP / regularized ML via brms-style
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md);
  hard bounds (`lower` / `upper`, and per-prior `lb` / `ub`).
- `sdreport` is deferred until standard errors are first needed (roughly
  a quarter off fit time; `se = TRUE` restores eager mode).
- Pluggable optimizers (`frmtmb_control(optimizer = )`), including
  arbitrary user functions.
- [`confint()`](https://rdrr.io/r/stats/confint.html) (Wald / profile /
  likelihood-root), natural-scale
  [`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md),
  [`anova()`](https://rdrr.io/r/stats/anova.html) likelihood-ratio
  tests,
  [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md).

### Diagnostics and ecosystem

- Residuals: response, pearson, and one-step-ahead (`type = "osa"`).
- DHARMa simulation-based residuals
  ([`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)).
- NUTS on the fitted objective:
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  (ML-mode initialization, named draws, priors, bounds) and
  [`check_laplace()`](https://aforren1.github.io/frmtmb/reference/check_laplace.md);
  bayesplot works on the draws.
- emmeans and marginaleffects support, both verified against glmmTMB.

### Verification

Roughly 500 tests compare fits against glmmTMB, lme4, mgcv, MASS,
survival, nnet, GLMMadaptive, and hand-written RTMB references, plus an
edge-case suite mined from the lme4 / glmmTMB / brms issue trackers
(dev/test-backlog.md).
