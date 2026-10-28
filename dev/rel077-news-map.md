| package | section (as filed) | source | record | bullet |
|---|---|---|---|---|
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | Under the covariance form of `ar()`, `ma()` and `arma()`, |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | So `sigma()`, `VarCorr()`, `summary()`, `confint()`, `vcov() |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | A prior on class `sigma` is a density on the innovation sd i |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | The standard error of `log(sigma)` is the one of the innovat |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | `ar(p = 2, cov = TRUE)` and the other higher orders, which b |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | `conditional_effects(method = "posterior_predict")` draws it |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | `cooks.distance()`, `dfbeta()`, `dfbetas()` and `influence() |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | A fit saved by frmtmb 0.76.0 or earlier with one of these te |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | Such a fit still answers where only the mean is read, becaus |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | `insight::get_sigma()` reads `sigma()` for a frmtmb fit. Bef |
| frmtmb | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | The sampling-API function `autocor_innov_sd()` is removed: n |
| frmtmb | Breaking changes | mointerop | `dev/mointerop-findings.md` | `mo()` takes a variable with two values, or with one, as brm |
| frmtmb | Breaking changes | mointerop | `dev/mointerop-findings.md` | marginaleffects' `vcov =` argument reaches the fit. Through  |
| frmtmb | Breaking changes | mointerop | `dev/mointerop-findings.md` | On a fit whose `mo()` term has a category that no row is in, |
| frmtmb | Breaking changes | mointerop | `dev/mointerop-findings.md` | marginaleffects' prediction arguments that frmtmb does not t |
| frmtmb | Bug fixes | mointerop | `dev/mointerop-findings.md` | `predictions()`, `comparisons()` and `slopes()` of marginale |
| frmtmb | Bug fixes | mointerop | `dev/mointerop-findings.md` | `predictions()` on new data with a factor that has levels no |
| frmtmb | Bug fixes | mointerop | `dev/mointerop-findings.md` | `vcov_cluster()` on a fit that lost a parameter (a `mo()` ca |
| frmtmb | Documentation | mointerop | `dev/mointerop-findings.md` | New help topic `?frmtmb-marginaleffects`: the parameters |
| frmtmb | Breaking changes | reslot | `dev/reslot-findings.md` | Each random-effect parameter, and the residual correlation o |
| frmtmb | Breaking changes | reslot | `dev/reslot-findings.md` | A centered intercept has one prior slot, class `Intercept`. |
| frmtmb | Breaking changes | reslot | `dev/reslot-findings.md` | `default_prior(route = "sample")` prints brms's `(flat)` wit |
| frmtmb | Breaking changes | reslot | `dev/reslot-findings.md` | The named list of prior objects, `prior = list(theta_1 = |
| frmtmb | Breaking changes | reslot | `dev/reslot-findings.md` | `set_prior("", class = ...)` with no bound is brms's `(flat) |
| frmtmb | Breaking changes | reslot | `dev/reslot-findings.md` | `c()` of a prior specification with something that is not on |
| frmtmb | For extension authors | reslot | `dev/reslot-findings.md` | A prior density can be brms's `(flat)`, which the prior reso |
| frmtmb | New features | priordraws | `dev/priordraws-findings.md` | `build_objective()` reads a third field on a copy of the fra |
| frmtmb | New features | priordraws | `dev/priordraws-findings.md` | `plot()` of a hypothesis draws the density of its prior draw |
| frmtmb | Bug fixes | consolidation | `dev/round-20261028.md` | `prior =` given a list of `set_prior()` specifications, such |
| frmtmb.sample | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | Under the covariance form of `ar()`, `ma()` and `arma()`, |
| frmtmb.sample | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | The default `student_t(3, 0, s)` on `sigma` and a prior of y |
| frmtmb.sample | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | We compared the draws with brms's own (4 chains of 4000 on e |
| frmtmb.sample | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | `conditional_effects(method = "posterior_predict")` on draws |
| frmtmb.sample | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | Draws saved with frmtmb.sample 0.23.0 for a model with one o |
| frmtmb.sample | Breaking changes | sigmainn | `dev/sigmainn-findings.md` | The message of the default priors no longer marks `sigma` wi |
| frmtmb.sample | New features | mointerop | `dev/mointerop-findings.md` | `frm_sample()` samples a `mo()` term on two values, as brms  |
| frmtmb.sample | New features | mointerop | `dev/mointerop-findings.md` | `frm_sample()` refuses a `mo()` term on one value, which is  |
| frmtmb.sample | Breaking changes | reslot | `dev/reslot-findings.md` | `prior = "flat"` is brms's `(flat)` in every slot: no densit |
| frmtmb.sample | Breaking changes | reslot | `dev/reslot-findings.md` | A bounds-only specification that replaces a default, |
| frmtmb.sample | Breaking changes | reslot | `dev/reslot-findings.md` | The named list of prior objects (`prior = list(theta_1 = |
| frmtmb.sample | Breaking changes | reslot | `dev/reslot-findings.md` | frmtmb.spatial's `car()` dependence parameter is the one par |
| frmtmb.sample | Breaking changes | reslot | `dev/reslot-findings.md` | `prior =` refuses a `default_prior()` or `validate_prior()`  |
| frmtmb.sample | Breaking changes | reslot | `dev/reslot-findings.md` | `delta`, the distance between equidistant ordinal thresholds |
| frmtmb.sample | Breaking changes | reslot | `dev/reslot-findings.md` | `prior_summary()` and the default-prior message print brms's |
| frmtmb.sample | Breaking changes | reslot | `dev/reslot-findings.md` | `set_prior(class = "theta")` on a standard deviation, a corr |
| frmtmb.sample | New features | priordraws | `dev/priordraws-findings.md` | `frm_sample(sample_prior = )` takes brms's `"no"` (the defau |
| frmtmb.sample | New features | priordraws | `dev/priordraws-findings.md` | New `prior_draws()`, brms's, with brms's `variable` and its  |
| frmtmb.sample | New features | priordraws | `dev/priordraws-findings.md` | The warmup iterations rstan saves (`save_warmup = TRUE`, the |
| frmtmb.sample | New features | priordraws | `dev/priordraws-findings.md` | `update()` of draws samples the same model again with change |
| frmtmb.sample | Breaking changes | priordraws | `dev/priordraws-findings.md` | The evidence ratio of a point hypothesis is brms's. `hypothe |
| frmtmb.sample | Breaking changes | priordraws | `dev/priordraws-findings.md` | `as_draws(inc_warmup = TRUE)` and `as.mcmc(inc_warmup = TRUE |
