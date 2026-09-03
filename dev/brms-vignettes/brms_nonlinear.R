# Translation of brms's "Estimating Non-Linear Models with brms"
# (brms 2.23.0, doc/brms_nonlinear.Rmd) onto the frmtmb surface.
#
# Read the comments, not just the code: every `<EDGE>` marker on a bv()
# call is a place where the brms line did not carry over verbatim.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_nonlinear")

set.seed(1234)

## ---- A simple non-linear model -------------------------------------

b <- c(2, 0.75)
x <- rnorm(100)
y <- rnorm(100, mean = b[1] * exp(b[2] * x))
dat1 <- data.frame(x, y)

# brms: prior(normal(1, 2), nlpar = "b1") + prior(normal(0, 2), nlpar = "b2").
# Those priors exist to put the sampler in the right region, so they are
# MCMC-only and are dropped. What replaces them is `start`: frmtmb begins
# at zero, where the gradient of b1 * exp(b2 * x) is NaN.
#
# brms's `b1 + b2 ~ 1`, one formula naming two parameters, is accepted
# verbatim. (The v0.34 audit recorded this as FN-8, a refusal; it is
# fixed.)
fit1 <- bv("model", "ML: fit1", {
  frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
      data = dat1, family = gaussian(),
      start = list(beta = c(1, 0)))
}, "SPELLING",
"start = replaces the identifying priors, because frmtmb begins at zero where the gradient of this body is NaN. family = gaussian() is house style only: frm() DOES default to gaussian, so the v0.34 audit's FN-1 is fixed")

# The split form is the same model, checked rather than assumed.
bv("post", "ML: b1 + b2 ~ 1 equals b1 ~ 1, b2 ~ 1", {
  alt <- frm(bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1, nl = TRUE),
             data = dat1, family = gaussian(), start = list(beta = c(1, 0)))
  stopifnot(all.equal(unlist(fixef(alt)), unlist(fixef(fit1)), tolerance = 1e-6))
  "identical"
}, NA_character_, "")

bv("post", "ML: summary(fit1)", summary(fit1), NA_character_, "")

# brms's plot(fit1) draws trace and density plots per parameter. frmtmb's
# plot.frmtmb_fit draws the two regression diagnostics instead: there is
# no chain to trace.
bv("post", "ML: plot(fit1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit1)
}, "BEHAVIOR", "brms draws trace/density per parameter; frmtmb draws residual and QQ diagnostics")

# A nonlinear predictor has no delta-method standard error, so the
# default Wald band is refused. band = "boot" refits.
bv("post", "ML: conditional_effects(fit1) [default band]", {
  conditional_effects(fit1)
}, "REFUSAL", "wald band on a nonlinear predictor; the message names band = 'boot'")

ce1 <- bv("post", "ML: conditional_effects(fit1, band = 'boot')", {
  conditional_effects(fit1, band = "boot", boot = 25, seed = 1)
}, "SPELLING", "needs band = 'boot' (25 refits here, 200 by default)")

bv("post", "ML: plot(ce1, points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(ce1, points = TRUE)
}, NA_character_, "")

## ---- The linear comparison -----------------------------------------

# The vignette writes brm(y ~ x, data = dat1), with no family at all.
# That ports verbatim. The bare call is run here rather than assumed, so
# the "frm() has no family default" claim in the v0.34 audit is measured
# and retired.
fit2 <- bv("model", "ML: fit2", {
  f <- frm(bf(y ~ x), data = dat1)
  stopifnot(identical(family(f)$family, "gaussian"))
  f
}, NA_character_, "")

bv("post", "ML: summary(fit2)", summary(fit2), NA_character_, "")

bv("post", "ML: pp_check(fit1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(fit1))
}, NA_character_, "")

bv("post", "ML: pp_check(fit2)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(fit2))
}, NA_character_, "")

# loo() needs draws. On a frmtmb_fit there are none.
bv("post", "ML: loo(fit1, fit2)", loo(fit1, fit2),
   "MISSING", "loo() is a frmtmb_draws method; the frequentist route is AIC()/anova()")

bv("post", "ML: AIC(fit1, fit2) [the substitute]", {
  print(AIC(fit1, fit2))
}, "SPELLING", "loo comparison becomes AIC/BIC; no elpd, no Pareto k")

## ---- The insurance loss model --------------------------------------

loss <- brms::loss
bv("data", "head(loss)", print(utils::head(loss)), NA_character_, "")

# brms: three nlpar priors plus control = list(adapt_delta = 0.9). The
# priors' means become `start`; adapt_delta is MCMC-only.
fit_loss <- bv("model", "ML: fit_loss", {
  frm(bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
         ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE),
      data = loss, family = gaussian(),
      start = list(beta = c(5000, 1, 45)))
}, "SPELLING", "the three nlpar prior means become start = list(beta = ...)")

bv("post", "ML: summary(fit_loss)", summary(fit_loss), "BEHAVIOR",
   "two divergences from brms's printed summary: an empty 'Coefficients (mu):' block appears on a nonlinear fit, and sigma is printed on its LINK scale (4.89) where brms reports the response-scale residual SD (exp(4.89) = 133)")

# brms: plot(fit_loss, N = 3, ask = FALSE). `N` selects how many
# parameters per trace page and has no meaning without chains.
bv("post", "ML: plot(fit_loss, N = 3, ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit_loss, ask = FALSE)
}, "BEHAVIOR", "N = is silently absorbed by ...; the panels are diagnostics, not traces")

bv("post", "ML: conditional_effects(fit_loss)", {
  conditional_effects(fit_loss)
}, "REFUSAL", "same nonlinear wald refusal as fit1")

## ---- THE DEEP DIVE: one facet per accident year --------------------
#
# brms's vignette line, verbatim:
#
#   conditions <- data.frame(AY = unique(loss$AY))
#   rownames(conditions) <- unique(loss$AY)
#   me_loss <- conditional_effects(fit_loss, conditions = conditions,
#                                  re_formula = NULL, method = "predict")
#   plot(me_loss, ncol = 5, points = TRUE)
#
# See dev/brms-vignette-audit.md, section "The conditional_effects
# faceting diagnosis". In short: the DATA layer is right and the
# refusal is spurious. `method = "predict"` never uses a Wald band -
# it overwrites estimate__/lower__/upper__ from simulated responses -
# but the nonlinear guard fires before the method is consulted.

conditions <- data.frame(AY = unique(loss$AY))
rownames(conditions) <- unique(loss$AY)

bv("post", "ML: conditional_effects(conditions, re_formula = NULL, method = 'predict')", {
  conditional_effects(fit_loss, conditions = conditions,
                      re_formula = NULL, method = "predict")
}, "REFUSAL",
"the vignette's exact call; refused by the nonlinear-wald guard even though method = 'predict' never builds a wald band")

# The workaround that gets the per-year curves out today.
me_loss <- bv("post", "ML: conditional_effects(conditions, re_formula = NULL, band = 'boot')", {
  conditional_effects(fit_loss, conditions = conditions,
                      re_formula = NULL, band = "boot", boot = 25, seed = 1)
}, "SPELLING", "method = 'predict' becomes band = 'boot'; the band is then a confidence band, not a prediction interval")

# Proof that the grid and the random-effect conditioning are correct:
# ten condition sets, ten different curves, labeled by the row names.
bv("post", "ML: per-year curves present in the returned data", {
  d <- me_loss[[1]]
  stopifnot(!is.null(d$cond__), length(unique(d$cond__)) == 10L)
  print(round(tapply(d$estimate__, d$cond__, max), 1))
}, NA_character_, "")

# And the contrast: re_formula = NA (the default) collapses all ten to
# the population curve, exactly as brms documents.
bv("post", "ML: re_formula = NA collapses the ten conditions", {
  d <- conditional_effects(fit_loss, conditions = conditions,
                           band = "boot", boot = 25, seed = 1)[[1]]
  print(round(tapply(d$estimate__, d$cond__, max), 1))
}, NA_character_, "")

# plot(): frmtmb loops over cond__ and draws one BASE-GRAPHICS panel per
# condition, so ten plots go to the device in sequence. `ncol = 5` is
# absorbed by ... and ignored; there is no facet grid, and in a script
# only the last panel survives on a single-page device. That is the
# "group-level single plot" report.
bv("post", "ML: plot(me_loss, ncol = 5, points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(me_loss, ncol = 5, points = TRUE, ask = FALSE)
}, "BEHAVIOR", "ten sequential base-graphics panels, not a 2x5 facet grid; ncol is silently ignored")

# The workaround for the grid.
bv("post", "ML: par(mfrow) workaround for the facet grid", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  op <- graphics::par(mfrow = c(2, 5)); on.exit(graphics::par(op), add = TRUE)
  plot(me_loss, points = TRUE, ask = FALSE)
}, "SPELLING", "the caller sets par(mfrow = c(2, 5)) instead of passing ncol")

## ---- Item response theory ------------------------------------------

inv_logit <- function(x) 1 / (1 + exp(-x))
ability <- rnorm(300)
p <- 0.33 + 0.67 * inv_logit(ability)
answer <- ifelse(runif(300, 0, 1) < p, 1, 0)
dat_ir <- data.frame(ability, answer)

fit_ir1 <- bv("model", "ML: fit_ir1", {
  frm(bf(answer ~ ability), data = dat_ir, family = bernoulli())
}, NA_character_, "")

bv("post", "ML: summary(fit_ir1)", summary(fit_ir1), NA_character_, "")

bv("post", "ML: plot(conditional_effects(fit_ir1), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit_ir1), points = TRUE)
}, NA_character_, "")

# bernoulli("identity"): brms's guessing-floor model needs the identity
# link so that 0.33 + 0.67 * inv_logit(eta) IS the probability.
fit_ir2 <- bv("model", "ML: fit_ir2", {
  frm(bf(answer ~ 0.33 + 0.67 * inv_logit(eta), eta ~ ability, nl = TRUE),
      data = dat_ir, family = bernoulli("identity"),
      start = list(beta = c(0, 0)))
}, "SPELLING", "prior(normal(0, 5), nlpar = 'eta') is dropped; start is at the same zero")

bv("post", "ML: summary(fit_ir2)", summary(fit_ir2), NA_character_, "")

bv("post", "ML: plot(conditional_effects(fit_ir2), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit_ir2, band = "boot", boot = 25, seed = 1),
       points = TRUE)
}, "SPELLING", "band = 'boot' again")

bv("post", "ML: loo(fit_ir1, fit_ir2)", loo(fit_ir1, fit_ir2),
   "MISSING", "as above")

# The 3PL model. brms bounds `guess` in [0, 1] with a beta(1, 1) prior
# carrying lb/ub; frmtmb spells the same bounds as lower =/upper =, which
# is a genuine identification constraint and not an MCMC convenience.
fit_ir3 <- bv("model", "ML: fit_ir3", {
  frm(bf(answer ~ guess + (1 - guess) * inv_logit(eta),
         eta ~ 0 + ability, guess ~ 1, nl = TRUE),
      data = dat_ir, family = bernoulli("identity"),
      start = list(beta = c(0, 0.3)),
      lower = c(guess = 0), upper = c(guess = 1))
}, "SPELLING", "prior(beta(1, 1), nlpar = 'guess', lb = 0, ub = 1) becomes lower =/upper =")

bv("post", "ML: summary(fit_ir3)", summary(fit_ir3), NA_character_, "")

bv("post", "ML: plot(fit_ir3)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit_ir3)
}, "BEHAVIOR", "diagnostics, not traces")

bv("post", "ML: plot(conditional_effects(fit_ir3), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit_ir3, band = "boot", boot = 25, seed = 1),
       points = TRUE)
}, "SPELLING", "band = 'boot'")


## ============ PATH 2: sampling (frm_sample) ==========================
#
# The brms vignette shows posterior output, so this path is the closer
# analogue. Route choice: the FIT route, frm_sample(fit, ...), for every
# model here. The brms priors on the nonlinear parameters were reproduced
# as `start` on the ML path, and the fit route carries the same
# parameterization and the same starting point into the chain. The
# formula route is used once, for fit_loss, to show what the default
# priors change.
#
# Sizes: chains = 1, iter = 400, warmup = 200. brms's vignette runs
# 4 chains of 2000. Short chains are recorded where they cost the
# vignette's illustrated output, not hidden.

options(mc.cores = 1)

s1 <- bv("model", "SAMPLE: fit1", {
  frm_sample(fit1, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "SPELLING", "brm(...) becomes frm() then frm_sample(); the priors stay dropped, so this chain runs under frmtmb's defaults and not brms's normal(1, 2) / normal(0, 2)")

bv("post", "SAMPLE: summary(fit1)", summary(s1), "BEHAVIOR",
   "the draws summary prints mean/sd/2.5%/97.5%/n_eff/Rhat, not brms's Estimate/Est.Error/l-95% CI/u-95% CI/Rhat/Bulk_ESS/Tail_ESS; at 200 post-warmup draws the ESS columns are far below what the vignette shows")

# brms's plot(fit1) is the trace and density display. On draws frmtmb
# routes that through bayesplot, as brms does.
bv("post", "SAMPLE: plot(fit1) -> mcmc_plot(s1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(mcmc_plot(s1, type = "trace"))
}, "SPELLING", "plot(fit) becomes mcmc_plot(draws, type = 'trace'); the default plot() method is still the residual diagnostics")

bv("post", "SAMPLE: conditional_effects(fit1)", {
  conditional_effects(s1)
}, NA_character_,
"the nonlinear wald refusal is gone: the draws ARE the band, so no delta method is needed")

s2 <- bv("model", "SAMPLE: fit2", {
  frm_sample(fit2, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "SPELLING", "as above")

bv("post", "SAMPLE: pp_check(fit1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(s1))
}, NA_character_, "")

# The call the ML path had to give up. It works here.
bv("post", "SAMPLE: loo(fit1, fit2)", {
  l1 <- loo(s1); l2 <- loo(s2)
  print(loo_compare(l1, l2))
}, "SPELLING",
"brms writes loo(fit1, fit2) with two fits in one call; frmtmb takes one fit per loo() and compares with loo_compare()")

fl_s <- bv("model", "SAMPLE: fit_loss (fit route)", {
  frm_sample(fit_loss, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "BEHAVIOR",
"frm_sample() reports that the AY block stays CENTERED because its sd has a flat prior on the fit route; brms always has a prior there, so brms's chain is non-centered and mixes differently")

# The formula route is where the default priors live. This is the closer
# match to brms, which puts a half-t on every sd.
bv("model", "SAMPLE: fit_loss (formula route)", {
  frm_sample(bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
                ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE),
             data = loss, family = gaussian(),
             start = list(beta = c(5000, 1, 45)),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "SPELLING",
"the formula route supplies brms-like default priors, which is what makes the non-centered parameterization available; brms needs no such choice")

bv("post", "SAMPLE: summary(fit_loss)", summary(fl_s), "BEHAVIOR",
   "the AY standard deviation appears as the raw internal parameter 'theta_1' on the log scale (6.53, so exp(6.53) = 686 against brms's 747.55) where brms prints 'sd(ult_Intercept)' on the natural scale; sigma_Intercept is on the link scale for the same reason")

# VarCorr() is the natural-scale reading of that same quantity.
bv("post", "SAMPLE: VarCorr(fit_loss) [the natural-scale reading]",
   VarCorr(fl_s), "SPELLING",
   "brms's summary shows the sd inline; on frmtmb draws the natural-scale value comes from VarCorr(), not summary()")

## ---- THE DEEP DIVE on the sampling path ----------------------------
#
# The draws method breaks in a DIFFERENT place from the fit method.
# `method =` does not exist here at all, and the refusal names the
# replacement. Drop it and the vignette's display comes out right.

bv("post", "SAMPLE: the vignette's exact call (method = 'predict')", {
  conditional_effects(fl_s, conditions = conditions, re_formula = NULL,
                      method = "predict")
}, "REFUSAL",
"conditional_effects() on draws has no method = at all; the message names posterior_predict() over your own grid, which is the right redirect")

me_s <- bv("post", "SAMPLE: conditions + re_formula = NULL", {
  conditional_effects(fl_s, conditions = conditions, re_formula = NULL)
}, "SPELLING",
"only method = 'predict' has to go; the per-year curves and cond__ come out with no other change, and the band is the posterior interval")

bv("post", "SAMPLE: per-year curves present in the returned data", {
  d <- me_s[[1]]
  stopifnot(!is.null(d$cond__), length(unique(d$cond__)) == 10L)
  print(round(tapply(d$estimate__, d$cond__, max), 1))
}, NA_character_, "")

bv("post", "SAMPLE: plot(me_loss, ncol = 5, points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(me_s, ncol = 5, points = TRUE, ask = FALSE)
}, "BEHAVIOR",
"identical to the fit path: ten sequential base-graphics panels, ncol ignored, no facet grid. The plot method is shared, so the faceting gap is not path-specific")

## ---- Item response theory, sampled ---------------------------------

s_ir1 <- bv("model", "SAMPLE: fit_ir1", {
  frm_sample(fit_ir1, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, NA_character_, "")

s_ir2 <- bv("model", "SAMPLE: fit_ir2", {
  frm_sample(fit_ir2, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, NA_character_, "")

bv("post", "SAMPLE: loo(fit_ir1, fit_ir2)", {
  print(loo_compare(loo(s_ir1), loo(s_ir2)))
}, "SPELLING", "one loo() per fit, then loo_compare()")

bv("post", "SAMPLE: plot(conditional_effects(fit_ir2), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s_ir2), points = TRUE)
}, NA_character_,
"no band = argument is needed on draws, so this line ports unchanged where the ML path needed band = 'boot'")

bv("post", "SAMPLE: hypothesis(s_ir1, 'ability > 0')", {
  hypothesis(s_ir1, "ability > 0")
}, NA_character_,
"brms's one-sided form ports verbatim, and reports method = posterior here")

# The same line on the fit object, checked rather than assumed. The
# v0.34 audit recorded it as refused (FN-5); it is not.
bv("post", "ML: hypothesis(fit_ir1, 'ability > 0')", {
  hypothesis(fit_ir1, "ability > 0")
}, NA_character_,
"accepted on a point fit too, reporting method = wald with a one-sided p and a half-open interval; only the reported quantity differs from brms's Evid.Ratio")

bv_done()
