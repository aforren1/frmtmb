# Translation of brms's "Estimating Distributional Models with brms"
# (brms 2.23.0, doc/brms_distreg.Rmd) onto the frmtmb surface.
#
# Read the comments, not just the code: every `<EDGE>` marker on a bv()
# call is a place where the brms line did not carry over verbatim. The
# reference for a BEHAVIOR edge is doc/brms_distreg.html, which holds
# brms's own printed output.
#
# The vignette is translated TWICE, once onto each inference path.
# PATH 1 is frm(), maximum likelihood with a Laplace approximation for
# the random effects. PATH 2 is frm_sample(), which hands the same TMB
# objective to tmbstan and gives back draws. The brms vignette shows a
# posterior workflow, so PATH 2 is the closer analogue and its edges
# count the same as PATH 1's.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_distreg")

options(mc.cores = 1)
set.seed(1234)

## ---- shared data ----------------------------------------------------

group <- rep(c("treat", "placebo"), each = 30)
symptom_post <- c(rnorm(30, mean = 1, sd = 2), rnorm(30, mean = 0, sd = 1))
dat1 <- data.frame(group, symptom_post)

bv("data", "head(dat1)", print(utils::head(dat1)), NA_character_, "")

# The vignette's own URL. brms ships no copy of these 250 rows, so the
# read stays where the vignette puts it.
zinb <- utils::read.csv("https://paul-buerkner.github.io/data/fish.csv")
bv("data", "head(zinb)", print(utils::head(zinb)), NA_character_, "")

# The vignette's own size (n = 200), which is small enough to keep.
dat_smooth <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
bv("data", "head(dat_smooth[, 1:6])",
   print(utils::head(dat_smooth[, 1:6])), NA_character_, "")

## ============ PATH 1: ML / Laplace (frm) ============

## ---- A simple distributional model ---------------------------------

# The vignette's line already carries family = gaussian(), so nothing
# changes but brm -> frm.
fit1 <- bv("model", "ML: fit1", {
  frm(bf(symptom_post ~ group, sigma ~ group),
      data = dat1, family = gaussian())
}, NA_character_, "")

# brms prints the sigma formula on a second "Formula:" line and gives
# every coefficient a posterior interval with Rhat and ESS. frmtmb keeps
# only the mu formula in the header and reports Wald z and p instead.
bv("post", "ML: summary(fit1)", summary(fit1), "BEHAVIOR",
   "the header drops the sigma ~ group formula that brms prints on its second Formula line; the coefficient blocks themselves match brms position for position")

# brms draws two trace and density pages here. frmtmb has no chain on
# this path, so plot.frmtmb_fit draws the regression diagnostics.
bv("post", "ML: plot(fit1, N = 2, ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit1, N = 2, ask = FALSE)
}, "BEHAVIOR", "brms draws trace and density pages; frmtmb draws residual and QQ diagnostics, and N = is absorbed by ...")

bv("post", "ML: plot(conditional_effects(fit1), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit1), points = TRUE)
}, NA_character_, "")

# The sigma predictor is reachable the same way it is in brms, which is
# the point of the section: the log-scale contrast is a real parameter.
bv("post", "ML: conditional_effects(fit1, dpar = 'sigma')", {
  conditional_effects(fit1, dpar = "sigma")
}, NA_character_, "")

# The two point hypotheses port verbatim. brms answers with an evidence
# ratio and a posterior probability, both of which need a prior and
# draws; frmtmb answers with a delta-method interval, a z and a p.
hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_grouptreat) = 0")
bv("post", "ML: hypothesis(fit1, hyp) [point]", {
  hypothesis(fit1, hyp)
}, "BEHAVIOR", "Evid.Ratio, Post.Prob and Star become a delta-method se, a z and a p; the estimates match brms position for position")

# The directional form is accepted. The v0.34 audit recorded this as
# FN-5, a refusal that did not name a replacement; it is fixed, and the
# print method now says which rows are one-sided.
hyp2 <- bv("post", "ML: hypothesis(fit1, hyp) [directional >]", {
  hypothesis(fit1,
             "exp(sigma_Intercept + sigma_grouptreat) > exp(sigma_Intercept)")
}, "BEHAVIOR", "brms reports Evid.Ratio 3999 and Post.Prob 1; frmtmb reports a one-sided p and a one-sided interval, so the direction is answered but the Bayes factor is not")

bv("post", "ML: plot(hyp, chars = NULL)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(hyp2, chars = NULL)
}, "BEHAVIOR", "brms draws the posterior density of the contrast; frmtmb draws the estimate with its confidence interval, because there is no posterior to shade")

## ---- Zero-inflated models ------------------------------------------

fit_zinb1 <- bv("model", "ML: fit_zinb1", {
  frm(bf(count ~ persons + child + camper),
      data = zinb, family = zero_inflated_poisson())
}, NA_character_, "")

# The four mu coefficients land on brms's posterior means. The zi row
# does not, and the difference is a scale, not an estimate: brms prints
# a constant zi under "Further Distributional Parameters" on the
# response scale (0.41), frmtmb prints the same quantity as a zi
# intercept on the logit scale (-0.37, whose inverse logit is 0.409).
bv("post", "ML: summary(fit_zinb1)", summary(fit_zinb1), "BEHAVIOR",
   "a constant zi is printed as a logit-scale intercept (-0.37) where brms prints the probability (0.41); plogis() converts")

bv("post", "ML: plot(conditional_effects(fit_zinb1), ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit_zinb1), ask = FALSE)
}, NA_character_, "")

fit_zinb2 <- bv("model", "ML: fit_zinb2", {
  frm(bf(count ~ persons + child + camper, zi ~ child),
      data = zinb, family = zero_inflated_poisson())
}, NA_character_, "")

# With zi predicted, brms is on the logit scale too, so the two
# summaries agree on every row.
bv("post", "ML: summary(fit_zinb2)", summary(fit_zinb2), NA_character_, "")

bv("post", "ML: plot(conditional_effects(fit_zinb2), ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit_zinb2), ask = FALSE)
}, NA_character_, "")

## ---- Additive distributional models --------------------------------

# brms adds chains = 2 and control = list(adapt_delta = 0.95). Both are
# MCMC-only and are dropped.
fit_smooth1 <- bv("model", "ML: fit_smooth1", {
  frm(bf(y ~ s(x1) + s(x2) + (1 | fac), sigma ~ s(x0) + (1 | fac)),
      data = dat_smooth, family = gaussian())
}, NA_character_, "")

# Three divergences in one printout, none of them an error. The header
# again drops the sigma formula. brms reports a smoothing hyperparameter
# sds() per smooth; frmtmb reports the wiggliness standard deviation in
# the random-effect block plus an effective degrees of freedom line. That
# edf line can print a negative value, which no effective degrees of
# freedom can be.
bv("post", "ML: summary(fit_smooth1)", summary(fit_smooth1), "BEHAVIOR",
   "the sigma formula is missing from the header, sds() becomes a wiggliness SD in the random-effect block, and the printed edf of s(x1) comes out negative")

bv("post", "ML: plot(conditional_effects(fit_smooth1), points = TRUE, ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit_smooth1), points = TRUE, ask = FALSE)
}, NA_character_, "")

## ---- ML path: named by the audit brief, absent from this vignette ---
#
# brms 2.23.0's brms_distreg.Rmd contains no conditional_smooths() call
# and no mixture model. Both are exercised here because the audit brief
# names them, and both are labeled so the vignette tally stays honest.

# conditional_smooths() is the smooth-term companion of
# conditional_effects() and the natural next call after a two-smooth
# distributional fit. It is still absent (v0.34 audit FN-14).
bv("post", "ML: appendix conditional_smooths(fit_smooth1)", {
  conditional_smooths(fit_smooth1)
}, "MISSING", "no frmtmb path; the smooth curves are reachable only through conditional_effects(), which gives the fitted mean and not the centered smooth")

# Mixture models: the brms family constructor is accepted verbatim.
fit_mix <- bv("model", "ML: appendix mixture(gaussian(), gaussian())", {
  frm(bf(y ~ x1), data = dat_smooth,
      family = mixture(gaussian(), gaussian()))
}, NA_character_, "")

# brms writes per-component predictors as bf(y ~ 1, mu1 ~ x1, mu2 ~ x1).
# frmtmb reads the main formula AS mu1, so naming mu1 a second time is a
# collision rather than a missing feature.
bv("model", "ML: appendix bf(y ~ 1, mu1 ~ x1, mu2 ~ x1)", {
  frm(bf(y ~ 1, mu1 ~ x1, mu2 ~ x1), data = dat_smooth,
      family = mixture(gaussian(), gaussian()))
}, "REFUSAL", "the message lists the free dpars but does not say that the main formula already occupies mu1, so it reads as if mu1 were unsupported")

fit_mix2 <- bv("model", "ML: appendix per-component predictors, frmtmb spelling", {
  frm(bf(y ~ x1, mu2 ~ x1), data = dat_smooth,
      family = mixture(gaussian(), gaussian()))
}, "SPELLING", "mu1 ~ x1 becomes the main formula y ~ x1, and only mu2 is named")

bv("post", "ML: appendix pp_mixture(fit_mix2)", {
  utils::head(pp_mixture(fit_mix2))
}, "MISSING", "pp_mixture() is a frmtmb_draws method; the point-fit route is mixture_probs()")

bv("post", "ML: appendix mixture_probs(fit_mix2) [the substitute]", {
  utils::head(mixture_probs(fit_mix2))
}, "SPELLING", "pp_mixture() becomes mixture_probs(); the columns are the same responsibilities without a posterior interval")

bv("post", "ML: appendix pp_check(fit_zinb1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(fit_zinb1))
}, NA_character_, "")

## ============ PATH 2: sampling (frm_sample) ============
#
# The FIT route is used throughout: this vignette sets no prior on any
# model, so there is nothing for the formula route's brms-like defaults
# to carry, and starting from the frm() fit reuses the objective that
# PATH 1 already built.
#
# Chain settings are fixed for every call below: one chain, 400
# iterations with 200 warmup, seed 1. brms runs 4 chains of 2000. A
# 200-draw posterior is enough to show which parts of the workflow
# EXIST; where it is not enough to reproduce the vignette's numbers,
# that is recorded as a BEHAVIOR row rather than fixed by lengthening
# the chain.

SAMP <- function(fit) {
  frm_sample(fit, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}

## ---- A simple distributional model ---------------------------------

s_fit1 <- bv("model", "SAMPLE: fit1", SAMP(fit1), "SPELLING",
   "brm(...) becomes frm() then frm_sample(fit, chains =, iter =, warmup =, seed =); the MCMC arguments come back on the second call, not the first")

# 200 draws is short enough that tmbstan raises a low-ESS warning on
# sigma_grouptreat. brms's own run does not, because it takes 4000.
bv("post", "SAMPLE: summary(fit1)", summary(s_fit1), "BEHAVIOR",
   "a bare mean/sd/quantile/n_eff/Rhat matrix, without brms's Family, Formula, Data and Draws header, and with a low-ESS warning that a 200-draw chain earns and brms's 4000-draw run does not")

# There is no plot method for a draws object, so plot() reaches
# plot.default and fails inside base graphics. The message is about
# list components and names nothing useful.
bv("post", "SAMPLE: plot(fit1, N = 2, ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(s_fit1, N = 2, ask = FALSE)
}, "MISSING", "no plot method for frmtmb_draws; the base-graphics message does not name mcmc_plot(), which is the replacement")

bv("post", "SAMPLE: mcmc_plot(fit1) [the substitute]", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(mcmc_plot(s_fit1))
}, "SPELLING", "plot(fit, N =) becomes mcmc_plot(); this is the call that draws what brms's plot() draws")

bv("post", "SAMPLE: plot(conditional_effects(fit1), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s_fit1), points = TRUE)
}, NA_character_, "")

# On draws the hypothesis method switches to method = "posterior", so
# the interval is a quantile interval and not a delta-method one. The
# evidence ratio and posterior probability are still absent.
bv("post", "SAMPLE: hypothesis(fit1, hyp) [point]", {
  hypothesis(s_fit1, hyp)
}, "BEHAVIOR", "method = posterior gives quantile intervals close to brms (0.97 and 1.83), but Evid.Ratio, Post.Prob and Star are still not reported")

s_hyp2 <- bv("post", "SAMPLE: hypothesis(fit1, hyp) [directional >]", {
  hypothesis(s_fit1,
             "exp(sigma_Intercept + sigma_grouptreat) > exp(sigma_Intercept)")
}, "BEHAVIOR", "the one-sided p is bounded below by 1/ndraws, so a 200-draw chain reports p = 0.00995 where brms reports Evid.Ratio 3999")

bv("post", "SAMPLE: plot(hyp, chars = NULL)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(s_hyp2, chars = NULL)
}, "BEHAVIOR", "still the estimate with its interval, not brms's shaded posterior density, even though the draws to build one are now present")

## ---- Zero-inflated models ------------------------------------------

s_zinb1 <- bv("model", "SAMPLE: fit_zinb1", SAMP(fit_zinb1), "SPELLING",
   "the two-call spelling again")

bv("post", "SAMPLE: summary(fit_zinb1)", summary(s_zinb1), "BEHAVIOR",
   "zi stays on the logit scale here too, so the row to compare with brms's zi = 0.41 is plogis(zi_Intercept)")

bv("post", "SAMPLE: plot(conditional_effects(fit_zinb1), ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s_zinb1), ask = FALSE)
}, NA_character_, "")

s_zinb2 <- bv("model", "SAMPLE: fit_zinb2", SAMP(fit_zinb2), "SPELLING",
   "the two-call spelling again")

bv("post", "SAMPLE: summary(fit_zinb2)", summary(s_zinb2), NA_character_, "")

bv("post", "SAMPLE: plot(conditional_effects(fit_zinb2), ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s_zinb2), ask = FALSE)
}, NA_character_, "")

# The vignette itself never calls loo(), but the ML path's only MISSING
# post-processing gap is loo(), so it is checked where it now works.
bv("post", "SAMPLE: loo(fit_zinb1) [recovered from the ML path]", {
  loo(s_zinb1)
}, "BEHAVIOR", "the elpd is reachable, but 200 draws leave 7 of 250 Pareto k above 0.57 and p_loo at 54 against 250 observations, so the number is not the one brms's 4000 draws would print")

bv("post", "SAMPLE: loo_compare(zinb1, zinb2)", {
  loo_compare(loo(s_zinb1), loo(s_zinb2))
}, "BEHAVIOR", "the comparison runs and separates the two models, but the printed table carries frmtmb's own p_worse and diag_ columns and not brms's elpd_diff/se_diff pair alone")

## ---- Additive distributional models --------------------------------

# The one model in this vignette that needs the FORMULA route. Its ML
# mode is singular in two variance components, and the fit route starts
# the chain there with a flat prior on every one of them: measured, that
# run takes 142 s, hits maximum treedepth on all 200 transitions and
# ends at Rhat 2.07. The formula route supplies the brms 2.23 default
# student_t priors on the intercepts and the standard deviations, and
# the same chain then takes 13 s and reaches Rhat 1.05.
s_smooth1 <- bv("model", "SAMPLE: fit_smooth1", {
  frm_sample(bf(y ~ s(x1) + s(x2) + (1 | fac), sigma ~ s(x0) + (1 | fac)),
             data = dat_smooth, family = gaussian(),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "SPELLING",
   "the formula route is required here: frm_sample(fit) inherits flat variance priors from the ML fit and does not mix, where brms's own defaults do")

bv("post", "SAMPLE: summary(fit_smooth1)", summary(s_smooth1), "BEHAVIOR",
   "the five variance components print as unnamed theta_1 to theta_5 rows on their internal scale, where brms names them sds(sx1_1) and sd(Intercept) under two labeled blocks; VarCorr() below is the only route to the named, natural-scale values")

bv("post", "SAMPLE: ranef(fit_smooth1)", ranef(s_smooth1), NA_character_, "")

bv("post", "SAMPLE: VarCorr(fit_smooth1)", VarCorr(s_smooth1), NA_character_, "")

bv("post", "SAMPLE: plot(conditional_effects(fit_smooth1), points = TRUE, ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s_smooth1), points = TRUE, ask = FALSE)
}, NA_character_, "")

# conditional_effects() on draws has neither `method =` nor `band =`.
# The band is built from the draws themselves, so there is nothing to
# choose, and the refusal says where the prediction interval lives.
bv("post", "SAMPLE: conditional_effects(fit1, method = 'predict')", {
  conditional_effects(s_fit1, method = "predict")
}, "REFUSAL", "the draws method takes no method = argument; the message names posterior_predict(), which is the right replacement")

## ---- SAMPLE path: what the ML path could not reach ------------------

bv("post", "SAMPLE: waic(fit1)", waic(s_fit1), NA_character_, "")

bv("post", "SAMPLE: bayes_R2(fit1)", bayes_R2(s_fit1), NA_character_, "")

bv("post", "SAMPLE: pp_check(fit_zinb1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(s_zinb1))
}, NA_character_, "")

bv("post", "SAMPLE: posterior_epred / posterior_predict (fit1)", {
  ep <- posterior_epred(s_fit1)
  pr <- posterior_predict(s_fit1)
  c(epred = dim(ep), predict = dim(pr))
}, NA_character_, "")

bv("post", "SAMPLE: as_draws_df(fit1)", {
  utils::head(as_draws_df(s_fit1), 3)
}, NA_character_, "")

bv_done()
