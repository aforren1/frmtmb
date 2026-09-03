# Translation of brms's "Estimating Multivariate Models with brms"
# (brms 2.23.0, doc/brms_multivariate.Rmd) onto the frmtmb surface.
#
# Read the comments, not just the code: every edge label on a bv() call
# marks a place where the brms line did not carry over verbatim.
#
# The script covers both inference paths. PATH 1 is `frm()` and the
# frequentist post-processing surface; PATH 2 is `frm_sample()` and the
# draws surface, which is the closer analogue of what the vignette
# literally prints. Every label carries an "ML: " or "SAMPLE: " prefix.
#
# Sizes. PATH 1 uses the whole of BTdata, the vignette's own 828 rows
# over 106 dams and 104 fosternests, because each fit takes under 2 s.
# PATH 2 keeps 30 dams, about 220 rows, for the reason measured in a
# BEHAVIOR row at the head of that section.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_multivariate")

options(mc.cores = 1)

## ---- The data --------------------------------------------------------
#
# The vignette loads BTdata from MCMCglmm. brms ships no copy of it, so
# a missing MCMCglmm falls back to a simulated bivariate set with the
# same crossed structure. The edge labels stand either way; only the
# printed numbers change.

bv_sim <- !requireNamespace("MCMCglmm", quietly = TRUE)
if (!bv_sim) {
  utils::data("BTdata", package = "MCMCglmm", envir = environment())
  bv_sim <- !exists("BTdata", inherits = FALSE)
}

if (bv_sim) {
  # SUBSTITUTED DATA. Neither MCMCglmm nor brms supplies BTdata here.
  message("NOTE: MCMCglmm absent; simulated bivariate data substituted")
  set.seed(11)
  bv_n <- 500
  dam <- factor(sample(paste0("d", 1:60), bv_n, TRUE))
  fosternest <- factor(sample(paste0("f", 1:60), bv_n, TRUE))
  u_d <- matrix(rnorm(120), 60, 2) %*%
    chol(matrix(c(0.22, -0.06, -0.06, 0.07), 2))
  u_f <- matrix(rnorm(120), 60, 2) %*%
    chol(matrix(c(0.07, 0.06, 0.06, 0.11), 2))
  sex <- factor(sample(c("Fem", "Male", "UNK"), bv_n, TRUE,
                       c(0.47, 0.47, 0.06)))
  hatchdate <- rnorm(bv_n)
  BTdata <- data.frame(
    tarsus = -0.4 + 0.77 * (sex == "Male") - 0.04 * hatchdate +
      u_d[dam, 1] + u_f[fosternest, 1] + rnorm(bv_n, 0, 0.76),
    back = -0.01 + 0.01 * (sex == "Male") - 0.09 * hatchdate +
      u_d[dam, 2] + u_f[fosternest, 2] + rnorm(bv_n, 0, 0.9),
    dam = dam, fosternest = fosternest,
    hatchdate = hatchdate, sex = sex)
}

bv("data", "head(BTdata)", {
  cat("rows:", nrow(BTdata), " simulated:", bv_sim, "\n")
  print(utils::head(BTdata))
}, NA_character_, "")


## ============ PATH 1: ML / Laplace (frm) ============

## ---- Basic multivariate models --------------------------------------
#
# `mvbind()`, the `|p|` and `|q|` cross-response correlation tags and
# `set_rescor(TRUE)` are all accepted verbatim. brms omits `family`, and
# so may frmtmb: the missing gaussian default was the v0.34 audit's FN-1
# and it is fixed on the multivariate path too. The house style states
# the family, so it is stated.

bform1 <- bf(mvbind(tarsus, back) ~ sex + hatchdate +
               (1 | p | fosternest) + (1 | q | dam)) + set_rescor(TRUE)

fit1 <- bv("model", "ML: fit1", {
  frm(bform1, data = BTdata, family = gaussian())
}, NA_character_, "")

# brms attaches a criterion to the fit and reuses it later. frmtmb has
# no such function, on either path.
bv("post", "ML: fit1 <- add_criterion(fit1, 'loo')", {
  add_criterion(fit1, "loo")
}, "MISSING", "add_criterion() does not exist; AIC()/BIC() are the frequentist payload but there is no way to attach one to the fit")

# The estimates agree with the vignette. Fosternest correlation 0.83
# against brms's 0.70, dam correlation -0.55 against -0.52, residual
# correlation -0.05 against -0.05, and every regression coefficient to
# two decimals. What differs is the report.
bv("post", "ML: summary(fit1)", summary(fit1), "BEHAVIOR",
   "the header prints an empty 'Family:' and only the first response's formula, where brms prints MV(gaussian, gaussian) and both formulas; each sigma is its own coefficient block on the log link scale rather than brms's response-scale sigma_tarsus and sigma_back")

bv("post", "ML: pp_check(fit1, resp = 'tarsus')", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(fit1, resp = "tarsus"))
}, "REFUSAL", "pp_check() refuses every multivariate fit and names no alternative, even though resp = would make the check well defined")

bv("post", "ML: pp_check(fit1, resp = 'back')", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(fit1, resp = "back"))
}, "REFUSAL", "the same refusal")

bv("post", "ML: bayes_R2(fit1)", bayes_R2(fit1), "MISSING",
   "bayes_R2() is a frmtmb_draws method; on a point fit there is no R2 of any kind and nothing is named")

## ---- The bf() + bf() spelling ---------------------------------------
#
# Adding two bf() objects, each with its own predictors, ports verbatim.

bf_tarsus <- bf(tarsus ~ sex + (1 | p | fosternest) + (1 | q | dam))
bf_back <- bf(back ~ hatchdate + (1 | p | fosternest) + (1 | q | dam))

fit2 <- bv("model", "ML: fit2", {
  frm(bf_tarsus + bf_back + set_rescor(TRUE),
      data = BTdata, family = gaussian())
}, NA_character_, "")

bv("post", "ML: fit2 <- add_criterion(fit2, 'loo')", {
  add_criterion(fit2, "loo")
}, "MISSING", "as above")

bv("post", "ML: summary(fit2)", summary(fit2), "BEHAVIOR",
   "the same header and link-scale divergences as summary(fit1)")

bv("post", "ML: loo(fit1, fit2)", loo(fit1, fit2), "MISSING",
   "loo() is a frmtmb_draws method; on the ML path the comparison is AIC() or anova(), and PATH 2 recovers the real thing")

# The vignette's reading, that the two models fit about equally well,
# does survive the rewrite: the likelihood-ratio test does not reject.
bv("post", "ML: AIC(fit1, fit2) and anova() [the substitute]", {
  print(AIC(fit1, fit2))
  anova(fit2, fit1)
}, "SPELLING", "the loo comparison becomes AIC plus a likelihood-ratio test; no elpd, no se_diff, no Pareto k, and anova() labels its rows with the whole multivariate formula")

## ---- Per-response families, lf() and a smooth ------------------------
#
# `lf(sigma ~ 0 + sex)` added to a bf() is the brms spelling for a
# distributional formula on one response. The v0.34 audit recorded lf()
# as missing (FN-10); it is exported now and the whole block ports
# verbatim, per-response families included.

bf_tarsus3 <- bf(tarsus ~ sex + (1 | p | fosternest) + (1 | q | dam)) +
  lf(sigma ~ 0 + sex) + skew_normal()
bf_back3 <- bf(back ~ s(hatchdate) + (1 | p | fosternest) + (1 | q | dam)) +
  gaussian()

fit3 <- bv("model", "ML: fit3", {
  frm(bf_tarsus3 + bf_back3 + set_rescor(FALSE), data = BTdata)
}, NA_character_, "")

bv("post", "ML: fit3 <- add_criterion(fit3, 'loo')", {
  add_criterion(fit3, "loo")
}, "MISSING", "as above")

# Both readings the vignette takes from this summary port: the log sigma
# of tarsus is largest for unknown sex (-0.41 against brms's -0.39), and
# alpha is negative (-1.30 against -1.24).
bv("post", "ML: summary(fit3)", summary(fit3), "BEHAVIOR",
   "the smooth is reported as an sd(wiggle) in the random-effect block plus an edf line, where brms prints sds(back_shatchdate_1) under a Smoothing Spline Hyperparameters heading")

bv("post", "ML: conditional_effects(fit3, 'hatchdate', resp = 'back')", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit3, "hatchdate", resp = "back"))
}, NA_character_, "")

## ---- The frequentist surface on a multivariate fit -------------------
#
# Not in the vignette, which is a posterior document throughout. These
# are the calls a frmtmb user reaches for instead, measured here because
# a multivariate fit is where they are most likely to run out.

bv("post", "ML: fixef(fit1)", fixef(fit1), NA_character_, "")

bv("post", "ML: confint(fit1)", confint(fit1), "BEHAVIOR",
   "the covariance parameters are named theta_1, theta_2, ... on their internal scale, so a reader cannot tell an sd row from a cor row")

bv("post", "ML: VarCorr(fit1)", VarCorr(fit1), NA_character_, "")

bv("post", "ML: rescor_matrix(fit1)", rescor_matrix(fit1), "SPELLING",
   "the summary does print the residual correlation, but the numeric matrix brms exposes through the summary rows needs its own accessor here")

bv("post", "ML: ranef(fit1)", {
  lapply(ranef(fit1), utils::head, 3)
}, "BEHAVIOR", "the list is keyed by frmtmb block name, 'tarsus 1 | dam + back 1 | dam [ID]', where brms keys it by grouping factor and names the columns tarsus_Intercept and back_Intercept")

bv("post", "ML: coef(fit1)", {
  lapply(coef(fit1), function(z) lapply(z, utils::head, 2))
}, NA_character_, "")

bv("post", "ML: predict(fit1, resp = 'back')", {
  utils::head(predict(fit1, resp = "back"))
}, NA_character_, "")

bv("post", "ML: residuals(fit1)", utils::head(residuals(fit1)),
   "REFUSAL", "residuals() refuses multivariate fits and names no substitute, although predict(resp = ) against the observed column does the job")

bv("post", "ML: fitted(fit3)", utils::head(stats::fitted(fit3), 3),
   "REFUSAL", "fitted() refuses multivariate fits as well, so the two accessors that would give per-response predictions are split: predict() works, fitted() does not")


## ============ PATH 2: sampling (frm_sample) ============
#
# The vignette prints posterior summaries, a loo comparison and a
# Bayesian R2, so this path is the closer analogue and it is where two
# of the MISSING items above come back.
#
# Size. Sampling time on this model follows the geometry of the two
# crossed correlated blocks, not the row count: a random 30-dam subset
# of about 220 rows samples in 8 s per model, while the first 30 dam
# levels, 201 rows, took 30 to 60 s, and the full 828 rows do not finish
# in the budget. PATH 2 keeps the random 30 dams.
#
# Route. brms sets no prior on any of the three models, so the FIT route
# is the literal translation. It is measured once, on fit1, and it is
# the wrong choice here: a fit carries flat priors, flat on a
# cross-response correlation is the improper (1 - rho^2)^-3/2, and the
# chain drifts to correlations at plus or minus one. frm_sample() says
# exactly that and names set_prior(class = "cor"). The rest of the
# section therefore uses the FORMULA route, which supplies brms 2.23's
# own defaults.

set.seed(1234)
bv_dams <- sample(levels(factor(BTdata$dam)), 30)
BTsmall <- droplevels(BTdata[BTdata$dam %in% bv_dams, ])

bv("data", "SAMPLE: the shrunken data", {
  cat("rows:", nrow(BTsmall), " dams:", nlevels(BTsmall$dam),
      " fosternests:", nlevels(BTsmall$fosternest), "\n")
  nrow(BTsmall)
}, "BEHAVIOR", "the sample path needs about 220 of the vignette's 828 rows; the full data does not finish a chain in the time budget")

bv_draws <- function(...) {
  frm_sample(..., chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}

# The fit route, at half the iterations to keep it inside the budget.
# The note it prints, and the boundary correlations that follow, are the
# finding.
s1_fitroute <- bv("model", "SAMPLE: fit1 [fit route]", {
  frm_sample(frm(bform1, data = BTsmall, family = gaussian()),
             chains = 1, iter = 200, warmup = 100, seed = 1,
             cores = 1, refresh = 0)
}, "BEHAVIOR", "the fit route leaves both correlated blocks with flat priors, improper on a correlation; frm_sample() prints the reason and names set_prior(class = 'sd') and set_prior(class = 'cor')")

# The evidence for that row: both cross-response correlations land on
# the boundary, where brms reports 0.70 and -0.52.
bv("post", "SAMPLE: VarCorr(fit1) [fit route]", VarCorr(s1_fitroute),
   "BEHAVIOR", "the flat correlation prior pushes both cross-response correlations to the boundary at plus or minus one, so the vignette's 0.70 and -0.52 are not recoverable by this route")

s1 <- bv("model", "SAMPLE: fit1 [formula route]", {
  bv_draws(bform1, data = BTsmall, family = gaussian())
}, "SPELLING", "brm(bform1, chains = 2, cores = 2) becomes frm_sample(bform1, ...); the formula route is used because brms's implicit defaults, not flat priors, are what keep the correlations proper")

bv("post", "SAMPLE: fit1 <- add_criterion(fit1, 'loo')", {
  add_criterion(s1, "loo")
}, "MISSING", "add_criterion() does not exist on the draws path either; loo() computes the same criterion but nothing attaches it to the object")

bv("post", "SAMPLE: loo(fit1) [what add_criterion would have stored]", {
  loo(s1)
}, "SPELLING", "the criterion is printed rather than attached; a 200-draw chain flags many Pareto k values that brms's 2000 draws do not")

bv("post", "SAMPLE: summary(fit1)", summary(s1), "BEHAVIOR",
   "one flat mean/sd/quantile/n_eff/Rhat matrix where brms prints Multilevel Hyperparameters, per-response Regression Coefficients, Further Distributional Parameters and Residual Correlations as named blocks; the six covariance parameters arrive as theta_1 to theta_6 and the residual correlation as thetar_1")

# The named, natural-scale reading of the same parameters. This is the
# closest thing on the draws path to brms's Multilevel Hyperparameters
# block, and the vignette's two correlations are legible in it.
bv("post", "SAMPLE: VarCorr(fit1) [the natural-scale reading]",
   VarCorr(s1), "SPELLING",
   "brms shows the SDs and correlations inline in summary(); on frmtmb draws they need VarCorr()")

bv("post", "SAMPLE: rescor_matrix(fit1)", {
  rescor_matrix(s1)
}, "MISSING", "rescor_matrix() has no frmtmb_draws method and returns NULL in silence rather than refusing, so the residual correlation the vignette reports has no posterior summary here")

bv("post", "SAMPLE: pp_check(fit1, resp = 'tarsus')", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(s1, resp = "tarsus"))
}, "REFUSAL", "the same multivariate refusal as the ML path; having draws does not help")

bv("post", "SAMPLE: bayes_R2(fit1)", bayes_R2(s1), "BEHAVIOR",
   "the call the ML path had to give up runs here, but it returns R2tarsus alone where brms returns one row per response")

s2 <- bv("model", "SAMPLE: fit2", {
  bv_draws(bf_tarsus + bf_back + set_rescor(TRUE), data = BTsmall,
           family = gaussian())
}, "SPELLING", "the formula route again, for the same reason")

bv("post", "SAMPLE: summary(fit2)", summary(s2), "BEHAVIOR",
   "flat matrix and internal theta names, as above")

# The vignette's own line, two fits inside one loo() call.
bv("post", "SAMPLE: loo(fit1, fit2)", loo(s1, s2), "MISSING",
   "the two-model form errors with \"'list' object cannot be coerced to type 'integer'\" because the second fit lands in the ndraws argument; the message neither refuses nor names loo_compare()")

# The vignette's conclusion, no noteworthy difference, does survive.
bv("post", "SAMPLE: loo_compare(loo(fit1), loo(fit2)) [the substitute]", {
  loo_compare(loo(s1), loo(s2))
}, "SPELLING", "one loo() per fit and then loo_compare(); the elpd difference stays well inside its standard error, which is the vignette's reading")

s3 <- bv("model", "SAMPLE: fit3", {
  bv_draws(bf_tarsus3 + bf_back3 + set_rescor(FALSE), data = BTsmall)
}, "SPELLING", "the formula route on the skew-normal, lf() and smooth model")

bv("post", "SAMPLE: summary(fit3)", summary(s3), "BEHAVIOR",
   "flat matrix again; the smoothing SD and the covariance blocks all arrive as theta_ entries, and a 200-draw chain leaves Rhat well above 1.01, so the numbers are not the vignette's")

bv("post", "SAMPLE: conditional_effects(fit3, 'hatchdate', resp = 'back')", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s3, "hatchdate", resp = "back"))
}, NA_character_, "")

# brms's conditional_effects() takes method = "predict" for a predictive
# band. The draws method has neither method = nor band =.
bv("post", "SAMPLE: conditional_effects(..., method = 'predict')", {
  conditional_effects(s3, "hatchdate", resp = "back", method = "predict")
}, "REFUSAL", "conditional_effects() on draws has no method =; the message says the curves already are posterior draws and names posterior_predict() for a predictive band, which points right")

bv("post", "SAMPLE: posterior_epred(fit1)", {
  dim(posterior_epred(s1))
}, "BEHAVIOR", "the matrix is draws by observations for one response only; brms returns a draws by observations by response array, so the second response has no way out through this entry point")

bv("post", "SAMPLE: posterior_predict(fit1)", {
  dim(posterior_predict(s1))
}, "BEHAVIOR", "the same shape and the same missing response dimension")

bv("post", "SAMPLE: as_draws_df(fit1)", {
  utils::head(as_draws_df(s1)[, 1:5], 3)
}, NA_character_, "")

bv_done()
