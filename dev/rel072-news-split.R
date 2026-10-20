# NEWS sentences over 25 words (dev/rel072-news-long.R), split after the
# NEWS merge (dev/rel072-news-merge.R sources this), so the lanes' NEWS
# stays as each lane wrote it. Each edit is a fixed string that must
# match once in the file it names. Adapted from dev/rel071-news-split.R.
# Raw strings (r"--(...)--"), because the text holds both quote marks.
list(
  "NEWS.md" = list(
    c(r"--(  `lognormal()` each add a guard that is `NaN` outside the support and
  exactly 0 on it, its edges included as in Stan. `beta()` and `gamma()`
  relied on a `log()` that the tape drops where it is multiplied by a
  constant 0, so `beta(1, 1)` on a coefficient was flat over the whole
  line and `gamma(1, r)` finite below 0; `exponential()` had no guard at
  all and was finite below 0 on doubles too.)--",
      r"--(  `lognormal()` each add a guard. The guard is `NaN` outside the
  support and exactly 0 on it, its edges included as in Stan. `beta()`
  and `gamma()` relied on a `log()` that the tape drops where it is
  multiplied by a constant 0. So `beta(1, 1)` on a coefficient was flat
  over the whole line, and `gamma(1, r)` was finite below 0.
  `exponential()` had no guard at all and was finite below 0 on doubles
  too.)--"),
    c(r"--(* **A bounded density on a parameter with no such bound of its own
  warns, in brms's words** (`uniform()` or `beta()` on a coefficient,
  `gamma()` on an intercept), when neither its row nor the class row
  above it gives `lb` or `ub`, as brms 2.23.0 does.)--",
      r"--(* **A bounded density on a parameter with no such bound of its own
  warns, in brms's words.** Examples are `uniform()` or `beta()` on a
  coefficient and `gamma()` on an intercept. It warns when neither its
  row nor the class row above it gives `lb` or `ub`, as brms 2.23.0
  does.)--"),
    c(r"--(  Before, newdata was coded by its own levels: a relevelled or
  character `dec()` swapped the two boundaries in silence, and one value
  alone (a one-row newdata, or a logical `dec()` given `TRUE`) was
  refused for holding one level. A value the fit did not see is now
  refused by name, where it used to be coded as if it were one of the
  fit's.)--",
      r"--(  Before, newdata was coded by its own levels. A relevelled or
  character `dec()` swapped the two boundaries in silence. One value
  alone (a one-row newdata, or a logical `dec()` given `TRUE`) was
  refused for holding one level. A value the fit did not see is now
  refused by name. It used to be coded as if it were one of the fit's.)--"),
    c(r"--(  standard errors 6.2 percent off. Over ten designs with covariates of
  very different scales, poisson and bernoulli, the error was up to 21
  percent, against an exact covariance formed in 256 bits; it is now
  below 1e-11 on 18 of the 20 fits and 2.6e-8 on the other two, a
  covariate on [1000, 1001].)--",
      r"--(  standard errors 6.2 percent off. Take ten poisson and bernoulli
  designs with covariates of very different scales, against an exact
  covariance formed in 256 bits. The error was up to 21 percent. It is
  now below 1e-11 on 18 of the 20 fits, and 2.6e-8 on the other two (a
  covariate on [1000, 1001]).)--"),
    c(r"--(  Every fit without random effects moves: on 27 ordinary, well-scaled
  fits the standard errors moved by at most 5.2e-7 relative and the
  correlations by at most 1.1e-6, all toward the exact values.)--",
      r"--(  Every fit without random effects moves. On 27 ordinary, well-scaled
  fits the standard errors moved by at most 5.2e-7 relative, and the
  correlations by at most 1.1e-6. All moved toward the exact values.)--"),
    c(r"--(  the first `vcov()` builds it: faster than the differences, but on a
  bernoulli GLM with 1e4 rows the peak memory rose by 350 MB at 200
  coefficients and by 870 MB at 500.)--",
      r"--(  the first `vcov()` builds it. That is faster than the differences,
  but it takes more memory. On a bernoulli GLM with 1e4 rows the peak
  memory rose by 350 MB at 200 coefficients and by 870 MB at 500.)--"),
    c(r"--(  the largest. A poisson `y ~ x + (1 | g)` fit by importance sampling
  with x ~ N(500, 100) had an intercept standard error of 0.498 where
  the same model on x / 100, an exact reparameterization, gives 0.270:
  85 percent off.)--",
      r"--(  the largest. Take a poisson `y ~ x + (1 | g)` fit by importance
  sampling with x ~ N(500, 100). Its intercept standard error was 0.498.
  The same model on x / 100, an exact reparameterization, gives 0.270,
  so the first was 85 percent off.)--"),
    c(r"--(  A coefficient's step is 1e-3 divided by the root mean square of its
  design column, where that is above one, so that a step moves the
  linear predictor by about 1e-3 in root mean square over the rows. A
  row far out in a heavy tail moves by more: on a lognormal column with
  a root mean square of 50 and a largest value of 836, that row moves
  by 1.7e-2, and the standard errors are within 7e-6 of a Richardson
  reference, where 0.71.0 was 1.6e-2 off. A poisson GLMM with
  x ~ N(500, 100) and `(1 | g)` had standard errors 0.37 percent off
  and a correlation 0.10 off, and its group standard deviation lost its
  standard error as "flat", with a warning.)--",
      r"--(  Where the root mean square of a coefficient's design column is above
  one, the step is 1e-3 divided by it. So a step moves the linear
  predictor by about 1e-3 in root mean square over the rows. A row far
  out in a heavy tail moves by more. Take a lognormal column with a root
  mean square of 50 and a largest value of 836: that row moves by
  1.7e-2. The standard errors are within 7e-6 of a Richardson reference,
  where 0.71.0 was 1.6e-2 off. A poisson GLMM with x ~ N(500, 100) and
  `(1 | g)` had standard errors 0.37 percent off and a correlation 0.10
  off. Its group standard deviation lost its standard error as "flat",
  with a warning.)--"),
    c(r"--(  in its location, so its exact Hessian has no curvature there, and
  inverting it gave standard errors of exactly 0 and `confint()`
  intervals of zero width.)--",
      r"--(  in its location, so its exact Hessian has no curvature there.
  Inverting it gave standard errors of exactly 0 and `confint()`
  intervals of zero width.)--"),
    c(r"--(  start and the flags. Since `update()` was added (0.6.0) it evaluated
  the symbols of the stored call again, so the refit used the objects
  of that name at that time. With `x2` changed in the workspace,
  `update(fit)` gave a log-likelihood of -3597.28 for a fit of
  -3536.91, and with 10 rows deleted it had 1990 observations for the
  fit's 2000. With the data removed, after `saveRDS()` and `readRDS()`
  in a new session, or for a fit made inside a function, it stopped
  with "object 'd' not found" (or `'form'`, `'fam'`). A global `form`
  of the same name as the call's argument became the model without a
  word. The data the fit holds are the model frame and the variables
  the formula reads through a transform, on the rows the fit used,
  after `na.action`. The fit's bases go with them, so `update(fit)`
  with nothing changed gives the fit's log-likelihood, also when
  `na.action` removed rows and the formula has `splines::ns()`,
  `poly()` or `scale()`.)--",
      r"--(  start and the flags. Since `update()` was added (0.6.0) it evaluated
  the symbols of the stored call again. So the refit used the objects
  of that name at that time. With `x2` changed in the workspace,
  `update(fit)` gave a log-likelihood of -3597.28 for a fit of
  -3536.91. With 10 rows deleted it had 1990 observations for the
  fit's 2000. It stopped with "object 'd' not found" (or `'form'`,
  `'fam'`) in three cases. The data were removed, the fit was read
  back with `readRDS()` in a new session, or it was made inside a
  function. A global `form` of the same name as the call's argument
  became the model without a word. The data the fit holds are the
  model frame and the variables the formula reads through a transform.
  They are on the rows the fit used, after `na.action`. The fit's bases
  go with them. So `update(fit)` with nothing changed gives the fit's
  log-likelihood, also when `na.action` removed rows and the formula
  has `splines::ns()`, `poly()` or `scale()`.)--"),
    c(r"--(* **`drop1()` keeps the fit's bases.** A submodel keeps the knots of
  `ns()`, the coefficients of `poly()` and the center and scale of
  `scale()` that the fit computed from every row before `na.action`,
  as `stats::drop1()` keeps the fit's design columns, and their
  coefficient names. With three missing values in another variable,
  the likelihood-ratio statistics now equal `stats::drop1()`'s on
  `lm()` and `glm()`, and the AIC of a submodel equals that of
  `update(fit, ~ . - x3)`.)--",
      r"--(* **`drop1()` keeps the fit's bases.** A submodel keeps the bases
  that the fit computed from every row before `na.action`: the knots of
  `ns()`, the coefficients of `poly()`, and the center and scale of
  `scale()`. It keeps their coefficient names too, as `stats::drop1()`
  keeps the fit's design columns. With three missing values in another
  variable, the likelihood-ratio statistics now equal
  `stats::drop1()`'s on `lm()` and `glm()`. The AIC of a submodel
  equals that of `update(fit, ~ . - x3)`.)--"),
    c(r"--(  largest displacement, because each deletion computed the basis again
  from its own rows (since 0.15.0); they now agree to 4.6e-6 of it.)--",
      r"--(  largest displacement. Each deletion computed the basis again from its
  own rows (since 0.15.0). They now agree to 4.6e-6 of it.)--"),
    c(r"--(  `refcat = NA` a maximum-likelihood fit determines only the category
  probabilities, and its standard-error warning names the common shift
  of every category's predictor; a prior on each predictor determines
  the coefficients. `logLik()`, `AIC()`, `BIC()`, `anova()` and
  `drop1()` count the degrees of freedom of the same model with a
  reference category, which has the same likelihood; with a category
  formula of its own the count is not defined, the df is `NA`, and
  `anova()` and `drop1()` say why.)--",
      r"--(  `refcat = NA` a maximum-likelihood fit determines only the category
  probabilities. Its standard-error warning names the common shift of
  every category's predictor. A prior on each predictor determines the
  coefficients. `logLik()`, `AIC()`, `BIC()`, `anova()` and `drop1()`
  count the degrees of freedom of the same model with a reference
  category, which has the same likelihood. With a category formula of
  its own the count is not defined: the df is `NA`, and `anova()` and
  `drop1()` say why.)--"),
    c(r"--(* `logistic_normal()` fits inside `mvbf()`. Its correlations are a
  block of its response, reported as `lncor_<resp>(b,c)` in the summary
  and `lncor_<resp>__b__c` in `variables()` and `hypothesis()`, with
  `<resp>` brms's Stan name of the response (`Ya` for `Y_a`), and set
  with `set_prior(class = "lncor", resp = )`.)--",
      r"--(* `logistic_normal()` fits inside `mvbf()`. Its correlations are a
  block of its response, set with
  `set_prior(class = "lncor", resp = )`. The summary reports them as
  `lncor_<resp>(b,c)`, and `variables()` and `hypothesis()` as
  `lncor_<resp>__b__c`. `<resp>` is brms's Stan name of the response
  (`Ya` for `Y_a`).)--"),
    c(r"--(  `NaN` on doubles for an end category whose neighboring thresholds are
  tied, under every link with a log-odds form (logit, probit,
  probit_approx, cloglog), where the category has a positive
  probability.)--",
      r"--(  `NaN` on doubles for an end category whose neighboring thresholds are
  tied. The category has a positive probability there. This held under
  every link with a log-odds form (logit, probit, probit_approx,
  cloglog).)--"),
    c(r"--(  simplex families; `dharma_residuals()` on it is refused by name.)--",
      r"--(  simplex families. `dharma_residuals()` on it is refused by name.)--"),
    c(r"--(  terms.** The simulator gave each component the whole term list, so
  the j-th row of a class read the `trials()` of the j-th row of the
  data. In `simulate()`, `predict()` and frmtmb.sample's
  `posterior_predict()` of a binomial mixture with `N` of 2 and 50,
  25% of the draws exceeded `N`, and the `N = 2` rows averaged 11.8
  where `fitted()` says 0.98.)--",
      r"--(  terms.** The simulator gave each component the whole term list. So
  the j-th row of a class read the `trials()` of the j-th row of the
  data. Take a binomial mixture with `N` of 2 and 50. In `simulate()`,
  `predict()` and frmtmb.sample's `posterior_predict()`, 25% of the
  draws exceeded `N`, and the `N = 2` rows averaged 11.8 where
  `fitted()` says 0.98.)--"),
    c(r"--(* `drop1()` refits a term that reads a raw variable. `z ~ x + log(x2)`,
  `scale(x)`, `I(x^2)`, `splines::ns(x2, 3)`, `poly(x, 2)` beside
  another term, `x * log(x2)`, a transformed response, and
  `drop1(scope = )` stopped with "The model uses `x2`, which is not a
  column of `data`" (since 0.49.0), because the refits read the model
  frame, which holds `log(x2)` and not `x2`.)--",
      r"--(* `drop1()` refits a term that reads a raw variable. These stopped
  with "The model uses `x2`, which is not a column of `data`" (since
  0.49.0): `z ~ x + log(x2)`, `scale(x)`, `I(x^2)`,
  `splines::ns(x2, 3)`, `poly(x, 2)` beside another term,
  `x * log(x2)`, a transformed response, and `drop1(scope = )`. The
  refits read the model frame, which holds `log(x2)` and not `x2`.)--"),
    c(r"--(  did not read, so every step that removed a term stopped with "object
  'change' not found" (since 0.21.0, when drop1() was added). The call
  it evaluates now holds the fit's data and formula by value, so
  `step()` gives the same model with the data object changed,
  removed, or in a new session after `readRDS()`.)--",
      r"--(  did not read. So every step that removed a term stopped with "object
  'change' not found" (since 0.21.0, when drop1() was added). The call
  it evaluates now holds the fit's data and formula by value. So
  `step()` gives the same model when the data object is changed or
  removed, or in a new session after `readRDS()`.)--"),
    c(r"--(  was refused as matching no group-level term (since 0.61.0), and
  `make_conditions(fit, "x2")` refused `x2` of `log(x2)` (since
  0.65.0); both read the model frame.)--",
      r"--(  was refused as matching no group-level term (since 0.61.0).
  `make_conditions(fit, "x2")` refused `x2` of `log(x2)` (since
  0.65.0). Both read the model frame.)--"),
    c(r"--(  row names on a tibble is deprecated", and a fit on a named list keeps
  the raw variables a refit reads (both since 0.66.0).)--",
      r"--(  row names on a tibble is deprecated". A fit on a named list keeps
  the raw variables a refit reads. Both were defects since 0.66.0.)--"),
    c(r"--(  `multinomial()`, shims brms's Frechet, generalized extreme value
  and categorical functions onto frmtmb's families, and reaches
  frmtmb.eam's `wiener()`. Where brms checks only the shape of a log
  density or a CDF, a value check now compares it with brms's own
  function on the same arguments.)--",
      r"--(  `multinomial()`. It shims brms's Frechet, generalized extreme value
  and categorical functions onto frmtmb's families, and reaches
  frmtmb.eam's `wiener()`. Some of brms's checks read only the shape of
  a log density or a CDF. A value check now compares such a value with
  brms's own function on the same arguments.)--")),
  "extensions/frmtmb.sample/NEWS.md" = list(
    c(r"--(  link scale: `shape` (`gamma(0.01, 0.01)`, and `inv_gamma(0.4, 0.3)`
  for the negative binomials), `phi` (`gamma(0.01, 0.01)`), `nu`
  (`gamma(2, 0.1)`), `kappa`, the exgaussian `beta`, `alpha`, `zi`,
  `hu`, `zoi`, `coi`, `quantile` (`beta(1, 1)`), and frmtmb.eam's
  `bs`, `bias` and `ndt`.)--",
      r"--(  link scale: `shape` (`gamma(0.01, 0.01)`, and `inv_gamma(0.4, 0.3)`
  for the negative binomials), `phi` (`gamma(0.01, 0.01)`), `nu`
  (`gamma(2, 0.1)`) and `kappa`. So were the exgaussian `beta`,
  `alpha`, `zi`, `hu`, `zoi`, `coi` and `quantile` (`beta(1, 1)`), and
  frmtmb.eam's `bs`, `bias` and `ndt`.)--"),
    c(r"--(  it. The default-prior message now says when a specification applies
  to no parameter for this reason, and names the coef to write.)--",
      r"--(  it. The default-prior message now names a specification that
  applies to no parameter. That happens for this reason, or when a more
  specific specification of your own covers its parameters. The message
  names the coef to write.)--"),
    c(r"--(  move: a beta `phi` near 300 on 10 rows moved by 30 Monte Carlo
  standard errors, and a negative binomial `shape` and a Student-t `nu`
  had no proper posterior before.)--",
      r"--(  move. A beta `phi` near 300 on 10 rows moved by 30 Monte Carlo
  standard errors. A negative binomial `shape` and a Student-t `nu` had
  no proper posterior before.)--"),
    c(r"--(  natural-scale default that addressed a slot the model does not have.
  `sigma` without a formula now reads the response on the identity
  scale, as brms does, where a log `mu` link read it on the log scale.)--",
      r"--(  natural-scale default that addressed a slot the model does not have.
  `sigma` without a formula now reads the response on the identity
  scale, as brms does. Under a log `mu` link it read it on the log
  scale.)--"),
    c(r"--(  frmtmb's link does not hold the parameter inside it: a log-linked
  `ndt` is bounded above by the smallest response and an identity- or
  log-linked `nu` below by 1, as brms declares them, so Stan transforms
  the draws inside rather than meeting the density's edge.)--",
      r"--(  frmtmb's link does not hold the parameter inside it. A log-linked
  `ndt` is bounded above by the smallest response. An identity- or
  log-linked `nu` is bounded below by 1. brms declares them so, and
  Stan then transforms the draws inside rather than meeting the
  density's edge.)--"),
    c(r"--(  on. A binomial with `y | trials(N)` and newdata `N = 1000` drew 0 or
  1 (column mean 0.595 where `posterior_epred()` gave 576.3), and a
  gaussian with `se(s)` drew at sd 0.033 where newdata had `s = 5`.)--",
      r"--(  on. A binomial with `y | trials(N)` and newdata `N = 1000` drew 0 or
  1 (column mean 0.595 where `posterior_epred()` gave 576.3). A
  gaussian with `se(s)` drew at sd 0.033 where newdata had `s = 5`.)--"),
    c(r"--(* Draws of a `logistic_normal()` response inside `mvbf()` report its
  correlations as `lncor_<resp>__b__c`, with `<resp>` brms's Stan name
  of the response, and the sampler's default puts
  brms's `lkj(1)` on each such response.)--",
      r"--(* Draws of a `logistic_normal()` response inside `mvbf()` report its
  correlations as `lncor_<resp>__b__c`, with `<resp>` brms's Stan name
  of the response. The sampler's default puts brms's `lkj(1)` on each
  such response.)--"),
    c(r"--(  `cor_<group>__<coef1>__<coef2>`**, and a smooth's `sds_`, `rescor__`,
  `ar[1]`, `cosy` and `cortime__`, as columns of the draws matrix, so
  `as_draws_df(draws)$sd_g__Intercept` works as on a brms fit.)--",
      r"--(  `cor_<group>__<coef1>__<coef2>`**, and a smooth's `sds_`, `rescor__`,
  `ar[1]`, `cosy` and `cortime__`, as columns of the draws matrix. So
  `as_draws_df(draws)$sd_g__Intercept` works as on a brms fit.)--"),
    c(r"--(  scale. They handed it the reporting scale, which for a mixture's
  `theta` is already the softmax probability, and the simulator
  applied the softmax again, so the mixing weights were wrong.)--",
      r"--(  scale. They handed it the reporting scale, which for a mixture's
  `theta` is already the softmax probability. The simulator applied
  the softmax again, so the mixing weights were wrong.)--")),
  "extensions/frmtmb.eam/NEWS.md" = list(
    c(r"--(  refusal: frmtmb holds a numeric `dec()` variable at
  `round(median())`, a boundary, where it used to hold it at its mean
  (0.573 on the test data), which was drawn at the lower boundary in
  silence.)--",
      r"--(  refusal. frmtmb holds a numeric `dec()` variable at
  `round(median())`, a boundary. It used to hold it at its mean (0.573
  on the test data), which was drawn at the lower boundary in silence.)--"),
    c(r"--(  fix of the same release. A relevelled factor or a character vector
  was coded by its own levels, so the two boundaries swapped in
  silence; one value alone, and a one-row logical `dec()`, were refused
  for holding one level; a value the fit did not see was coded as one
  of the two boundaries.)--",
      r"--(  fix of the same release. A relevelled factor or a character vector
  was coded by its own levels, so the two boundaries swapped in
  silence. One value alone, and a one-row logical `dec()`, were refused
  for holding one level. A value the fit did not see was coded as one
  of the two boundaries.)--")),
  "extensions/frmtmb.latent/NEWS.md" = list(
    c(r"--(  simulator gave the state-dependent family the whole term list, so
  the j-th row of a state read the `trials()` of the j-th row of the
  data. In `simulate()` of a binomial `hmm(K = 2)` with `N` of 2 and
  50, 26% of the draws exceeded `N`, and the `N = 2` rows averaged 12.2
  where `fitted()` says 0.91.)--",
      r"--(  simulator gave the state-dependent family the whole term list. So
  the j-th row of a state read the `trials()` of the j-th row of the
  data. Take a binomial `hmm(K = 2)` with `N` of 2 and 50. In its
  `simulate()`, 26% of the draws exceeded `N`, and the `N = 2` rows
  averaged 12.2 where `fitted()` says 0.91.)--"),
    c(r"--(  formula, the family and the control from the fit. A fit whose data
  object is gone from the workspace, that was saved and read back in a
  new session, or that was made inside a function, now refits; since
  `hmm_starts()` was added (0.55.2 of frmtmb) it stopped with "the
  model's own call could not be re-evaluated from here". A symbol of
  the call that `update()` does not take from the fit, such as
  `na.action`, still has to be visible from where `hmm_starts()` is
  called.)--",
      r"--(  formula, the family and the control from the fit. Now a fit refits
  in three cases: its data object is gone from the workspace, it was
  saved and read back in a new session, or it was made inside a
  function. Since `hmm_starts()` was added (0.55.2 of frmtmb) such a fit
  stopped with "the model's own call could not be re-evaluated from
  here". Some symbols of the call are not taken from the fit, such as
  `na.action`. They still have to be visible from where `hmm_starts()`
  is called.)--")),
  "extensions/frmtmb.coupling/NEWS.md" = list(
    c(r"--(  interval is documented to stay inside (0, 1). frmtmb now takes the
  exact Hessian on a fit without random effects, so a coherence run off
  to 1 has a finite, very large standard error rather than a NaN one,
  and its interval came back with an upper end of exactly 1. The bound
  is one number on the logit scale, `qlogis(1 - eps)` (about 36.04),
  and it is the same at both ends, so a coherence run off to 0 is
  refused in the same way.)--",
      r"--(  interval is documented to stay inside (0, 1). frmtmb now takes the
  exact Hessian on a fit without random effects. So a coherence run off
  to 1 has a finite, very large standard error rather than a NaN one.
  Its interval came back with an upper end of exactly 1. The bound is
  one number on the logit scale, `qlogis(1 - eps)` (about 36.04). It is
  the same at both ends, so a coherence run off to 0 is refused in the
  same way.)--"))
)
