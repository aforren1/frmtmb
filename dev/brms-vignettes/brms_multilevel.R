# Translation of brms's "Advanced Bayesian Multilevel Modeling with the
# R Package brms" (brms 2.23.0, doc/brms_multilevel.ltx) onto the
# frmtmb surface.
#
# The vignette ships as pre-rendered JSS LaTeX. The code comes out of
# the Sinput environments; in this file the printed output shares those
# environments with the code, so the reference numbers quoted in the
# edge comments below come from the Sinput blocks that parse as output
# rather than as R.
#
# The script covers both inference paths. Section 1 is frm(), maximum
# likelihood through the Laplace approximation. Section 2 samples the
# same models with frm_sample() and runs the vignette's own posterior
# workflow, which is the closer analogue of what the paper shows.
#
# Read the comments, not just the code: every `<EDGE>` marker on a bv()
# call is a place where the brms line did not carry over verbatim.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_multilevel")

options(mc.cores = 1)
set.seed(4321)

## ---- Data setup ----------------------------------------------------
#
# Example 1. brms reads the UCLA fish data:
#   zinb <- read.csv("http://stats.idre.ucla.edu/stat/data/fish.csv")
# That host no longer serves the file, so the data are SIMULATED here
# from the coefficients the vignette itself reports (Intercept -1.01,
# persons 0.87, child -1.36, camper 0.80, zi 0.41) at the vignette's
# own 250 rows. The columns the vignette prints but never models
# (nofish, livebait, xb, zg) are carried so head() looks the same.
n_zinb <- 250L
zinb <- data.frame(
  nofish = rbinom(n_zinb, 1, 0.2),
  livebait = rbinom(n_zinb, 1, 0.85),
  camper = sample(0:1, n_zinb, TRUE),
  persons = sample(1:4, n_zinb, TRUE),
  child = sample(0:3, n_zinb, TRUE),
  xb = rnorm(n_zinb), zg = rnorm(n_zinb)
)
zinb$camper <- factor(zinb$camper, labels = c("no", "yes"))
zinb$count <- ifelse(
  runif(n_zinb) < 0.41, 0L,
  rpois(n_zinb, exp(-1.01 + 0.87 * zinb$persons - 1.36 * zinb$child +
                      0.80 * (zinb$camper == "yes"))))

bv("data", "head(zinb)", print(utils::head(zinb)), NA_character_, "")

# Example 2. brms: data("rent99", package = "gamlss.data"). The package
# is installed, so the data are the vignette's own. The vignette fits
# all 3082 rows; a 400-row subsample is used here because the sampled
# distributional model on the full data takes minutes.
data("rent99", package = "gamlss.data", envir = environment())
rent <- rent99[sample(nrow(rent99), 400L), ]

bv("data", "head(rent99)", print(utils::head(rent99)), NA_character_, "")

# Example 3. brms downloads ClarkTriangle.csv from GitHub. brms ships
# the same table as `loss`, so the download is replaced by the packaged
# copy. It carries one extra column, `premium`, which no model uses.
loss <- brms::loss

bv("data", "head(loss)", print(utils::head(loss)), NA_character_, "")

# Example 4. brms: data_mm <- sim_multi_mem(nschools = 10,
# nstudents = 1000, change = 0.1). `sim_multi_mem` lives in the paper's
# online supplement and is not exported by brms, so the simulator is
# written out here. The vignette's 1000 students become 300 to keep the
# fit small; 10 schools and a 10 percent change rate are unchanged.
n_school <- 10L
n_stud <- 300L
school_eff <- rnorm(n_school, 0, 3)
s1 <- sample(n_school, n_stud, TRUE)
s2 <- s1
moved <- seq_len(round(0.1 * n_stud))
s2[moved] <- vapply(s1[moved],
                    function(k) sample(setdiff(seq_len(n_school), k), 1L),
                    integer(1))
data_mm <- data.frame(s1 = s1, s2 = s2, w1 = 0.5, w2 = 0.5)
data_mm$y <- 19 + 0.5 * (school_eff[s1] + school_eff[s2]) +
  rnorm(n_stud, 0, 3.5)

bv("data", "head(data_mm)", print(utils::head(data_mm)), NA_character_, "")

bv("data", "data_mm[101:106, ]", print(data_mm[101:106, ]), NA_character_, "")

## ============ PATH 1: ML / Laplace (frm) ============

## ---- Example 1: Catching fish --------------------------------------
#
# brms: fit_zinb1 <- brm(count ~ persons + child + camper, data = zinb,
#                        family = zero_inflated_poisson("log"))
fit_zinb1 <- bv("model", "ML: fit_zinb1", {
  frm(bf(count ~ persons + child + camper), data = zinb,
      family = zero_inflated_poisson("log"))
}, NA_character_, "")

bv("post", "ML: summary(fit_zinb1)", summary(fit_zinb1), "BEHAVIOR",
   "brms prints zi as one 'Family Specific Parameters' row on the probability scale (0.41); frmtmb prints a 'Coefficients (zi):' block whose intercept is on the logit link")

# Figure me_zinb1. All three predictors come back.
bv("post", "ML: conditional_effects(fit_zinb1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(fit_zinb1)
  print(names(ce))
  plot(ce, ask = FALSE)
  names(ce)
}, NA_character_, "")

# brms: fit_zinb2 <- brm(bf(count ~ persons + child + camper,
#                           zi ~ child),
#                        data = zinb, family = zero_inflated_poisson())
# The two-part bf() ports verbatim.
fit_zinb2 <- bv("model", "ML: fit_zinb2", {
  frm(bf(count ~ persons + child + camper, zi ~ child), data = zinb,
      family = zero_inflated_poisson())
}, NA_character_, "")

bv("post", "ML: summary(fit_zinb2)", summary(fit_zinb2), "BEHAVIOR",
   "brms folds zi_Intercept and zi_child into the one Population-Level table; frmtmb prints a separate 'Coefficients (zi):' block")

# brms: LOO(fit_zinb1, fit_zinb2), printing 1639.52 and 1621.35.
bv("post", "ML: LOO(fit_zinb1, fit_zinb2)", LOO(fit_zinb1, fit_zinb2),
   "MISSING",
   "LOO() is a frmtmb_draws method; on a point estimate the frequentist route is AIC()/BIC()/anova()")

bv("post", "ML: AIC(fit_zinb1, fit_zinb2) [the substitute]", {
  print(stats::AIC(fit_zinb1, fit_zinb2))
}, "SPELLING", "the LOOIC table becomes an AIC table with no standard error on the difference")

## ---- Example 2: Housing rents --------------------------------------
#
# brms: fit_rent1 <- brm(rentsqm ~ t2(area, yearc) + (1|district),
#                        data = rent99, chains = 2, cores = 2)
# chains and cores go; the family default is gaussian in brms and now
# in frm() as well, so it is stated only for house style. (The v0.34
# audit recorded the missing default as FN-1; it is fixed.)
fit_rent1 <- bv("model", "ML: fit_rent1", {
  frm(bf(rentsqm ~ t2(area, yearc) + (1 | district)), data = rent,
      family = gaussian())
}, NA_character_, "")

bv("post", "ML: summary(fit_rent1)", summary(fit_rent1), "BEHAVIOR",
   "brms prints one 'Smooth Terms' block of sds(t2areayearc_1..3); frmtmb prints each wiggle SD as its own random-effect block plus an 'edf of the penalized part' line that brms has no counterpart for")

# brms: conditional_effects(fit_rent1, surface = TRUE), Figure me_rent2.
# (The v0.34 audit recorded this as FN-4, a numeric-matrix error; it is
# now a refusal that names the replacement.)
bv("post", "ML: conditional_effects(fit_rent1, surface = TRUE)", {
  conditional_effects(fit_rent1, surface = TRUE)
}, "REFUSAL",
"no surface display exists; the message names effects = 'x1:x2' as the substitute, which draws curves at three values of the second variable instead of a heat map")

bv("post", "ML: conditional_effects(fit_rent1, 'area:yearc') [the substitute]", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(fit_rent1, "area:yearc")
  plot(ce, ask = FALSE)
  names(ce)
}, "SPELLING", "surface = TRUE becomes a named two-variable effect")

# brms:
#   bform <- bf(rentsqm ~ t2(area, yearc) + (1|ID1|district),
#               sigma ~ t2(area, yearc) + (1|ID1|district))
#   fit_rent2 <- brm(bform, data = rent99, chains = 2, cores = 2)
# The |ID1| cross-formula correlation syntax ports verbatim.
fit_rent2 <- bv("model", "ML: fit_rent2", {
  bform <- bf(rentsqm ~ t2(area, yearc) + (1 | ID1 | district),
              sigma ~ t2(area, yearc) + (1 | ID1 | district))
  frm(bform, data = rent, family = gaussian())
}, NA_character_, "")

# brms's reference: sd(Intercept) 0.60, sd(sigma_Intercept) 0.11,
# cor 0.72.
bv("post", "ML: VarCorr(fit_rent2) [the group-level block]", VarCorr(fit_rent2),
   "BEHAVIOR",
   "the |ID1| block is found and the correlation between the two intercepts is positive as brms reports, but the ML estimate is pinned at the 1.00 boundary against brms's posterior mean of 0.72, the two SDs read 0.34 and 0.057 against 0.60 and 0.11, and three of the six smooth wiggle SDs collapse to zero")

# brms: conditional_smooths(fit_rent2), Figure me_rent3.
bv("post", "ML: conditional_smooths(fit_rent2)", conditional_smooths(fit_rent2),
   "MISSING",
   "no such function; ?conditional_effects says smooths are folded into it, but nothing isolates the spline part and nothing points a caller there from the missing name")

## ---- Example 3: Insurance loss payments ----------------------------
#
# brms:
#   nlform <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
#                ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1, nl = TRUE)
#   nlprior <- c(prior(normal(5000, 1000), nlpar = "ult"),
#                prior(normal(1, 2), nlpar = "omega"),
#                prior(normal(45, 10), nlpar = "theta"))
#   fit_loss1 <- brm(formula = nlform, data = loss, family = gaussian(),
#                    prior = nlprior, control = list(adapt_delta = 0.9))
#
# The vignette calls these priors mandatory for identifiability, so
# they are not MCMC-only and cannot simply be deleted. On the ML path
# they have to become `start`: a prior is added to the objective AFTER
# it is built, and the objective's gradient at the default zero start
# is already NaN, so priors alone do not rescue the fit. The next row
# measures that.
fit_loss1 <- bv("model", "ML: fit_loss1", {
  nlform <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
               ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE)
  frm(nlform, data = loss, family = gaussian(),
      start = list(beta = c(5000, 1, 45)))
}, "SPELLING",
"the three nlpar prior means become start = list(beta = ...); brms's `nlpar =` argument of prior() has no frmtmb spelling and the nearest target is dpar =")

bv("model", "ML: fit_loss1 with the priors instead of start", {
  frm(bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
         ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE),
      data = loss, family = gaussian(),
      prior = set_prior("normal(5000,1000)", class = "Intercept", dpar = "ult") +
        set_prior("normal(1,2)", class = "Intercept", dpar = "omega") +
        set_prior("normal(45,10)", class = "Intercept", dpar = "theta"))
}, "BEHAVIOR",
"the priors are accepted with dpar in place of brms's nlpar but they do not identify the fit: the objective is still evaluated at zero first and dies on a NaN gradient, and the message names start = and not the priors")

bv("post", "ML: summary(fit_loss1)", summary(fit_loss1), "BEHAVIOR",
   "ult 5292 against brms's 5273.70 and theta 45.90 against 46.07, but sd(ult_Intercept) reads 598 against brms's 745.74, sigma is on its LOG link scale (4.89, so 133 against brms's 139.93), and an empty 'Coefficients (mu):' block is printed on a nonlinear fit")

# brms: conditional_effects(fit_loss1), Figure me_loss1.
bv("post", "ML: conditional_effects(fit_loss1)", conditional_effects(fit_loss1),
   "REFUSAL",
   "a nonlinear predictor has no delta-method standard error, so the default Wald band is refused; the message names band = 'boot' and dpar =, both of which work")

me_loss <- bv("post", "ML: conditional_effects(fit_loss1, band = 'boot')", {
  conditional_effects(fit_loss1, band = "boot", boot = 25, seed = 1)
}, "SPELLING", "needs band = 'boot' (25 refits here, 200 by default)")

# brms:
#   conditions <- data.frame(AY = unique(loss$AY))
#   rownames(conditions) <- unique(loss$AY)
#   me_year <- conditional_effects(fit_loss1, conditions = conditions,
#                                  re_formula = NULL, method = "predict")
#   plot(me_year, ncol = 5, points = TRUE)
conditions <- data.frame(AY = unique(loss$AY))
rownames(conditions) <- unique(loss$AY)

bv("post", "ML: conditional_effects(conditions, re_formula = NULL, method = 'predict')", {
  conditional_effects(fit_loss1, conditions = conditions,
                      re_formula = NULL, method = "predict")
}, "REFUSAL",
"the vignette's exact call; the nonlinear Wald guard fires before method = 'predict' is consulted, even though that method overwrites the band from simulated responses and never builds a Wald one")

me_year <- bv("post", "ML: conditional_effects(conditions, re_formula = NULL, band = 'boot')", {
  conditional_effects(fit_loss1, conditions = conditions,
                      re_formula = NULL, band = "boot", boot = 25, seed = 1)
}, "SPELLING",
"method = 'predict' becomes band = 'boot', so the band is a confidence band and not the prediction interval the vignette's figure shows")

bv("post", "ML: plot(me_year, ncol = 5, points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(me_year, ncol = 5, points = TRUE, ask = FALSE)
}, "BEHAVIOR",
"ten sequential base-graphics panels, not brms's 2x5 facet grid; ncol is absorbed by ... and ignored, and par(mfrow) is the workaround")

# brms:
#   nlform2 <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
#                 ult ~ 1 + (1|ID1|AY), omega ~ 1 + (1|ID1|AY),
#                 theta ~ 1 + (1|ID1|AY), nl = TRUE)
#   fit_loss2 <- update(fit_loss1, formula = nlform2,
#                       control = list(adapt_delta = 0.90))
fit_loss2 <- bv("model", "ML: fit_loss2", {
  nlform2 <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
                ult ~ 1 + (1 | ID1 | AY), omega ~ 1 + (1 | ID1 | AY),
                theta ~ 1 + (1 | ID1 | AY), nl = TRUE)
  update(fit_loss1, formula = nlform2)
}, NA_character_, "")

bv("post", "ML: VarCorr(fit_loss2)", VarCorr(fit_loss2), NA_character_, "")

# The one-formula convenience the vignette names next:
#   ult + omega + theta ~ 1 + (1|ID1|AY)
# (The v0.34 audit recorded a multi-parameter nlpar formula as FN-8, a
# rejection; it is fixed.)
bv("model", "ML: ult + omega + theta ~ 1 + (1|ID1|AY)", {
  frm(bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
         ult + omega + theta ~ 1 + (1 | ID1 | AY), nl = TRUE),
      data = loss, family = gaussian(), start = list(beta = c(5000, 1, 45)))
}, NA_character_, "")

# brms: LOO(fit_loss1, fit_loss2), printing 715.44 and 720.60.
bv("post", "ML: LOO(fit_loss1, fit_loss2)", LOO(fit_loss1, fit_loss2),
   "MISSING", "as above")

bv("post", "ML: AIC(fit_loss1, fit_loss2) [the substitute]", {
  print(stats::AIC(fit_loss1, fit_loss2))
}, "SPELLING",
"AIC 733.4 against 737.0 keeps brms's ordering (the simpler model wins) but is not on the LOOIC scale")

## ---- Example 4: Performance of school children ---------------------
#
# brms: fit_mm <- brm(y ~ 1 + (1 | mm(s1, s2)), data = data_mm)
# (The v0.34 audit recorded mm() as FN-2, a missing function; it is
# fixed, and both the plain and the weighted form now run.)
fit_mm <- bv("model", "ML: fit_mm", {
  frm(bf(y ~ 1 + (1 | mm(s1, s2))), data = data_mm, family = gaussian())
}, NA_character_, "")

bv("post", "ML: summary(fit_mm)", summary(fit_mm), "BEHAVIOR",
   "the grouping factor is labeled 'mm(s1, s2)' where brms labels it '~mms1s2', and sigma is on its LOG link scale (1.25, so 3.5) where brms reports 3.58")

bv("post", "ML: ngrps(fit_mm)", ngrps(fit_mm), NA_character_, "")

# brms amends the weights and refits:
#   data_mm[1:100, "w1"] <- runif(100, 0, 1)
#   data_mm[1:100, "w2"] <- 1 - data_mm[1:100, "w1"]
#   fit_mm2 <- brm(y ~ 1 + (1 | mm(s1, s2, weights = cbind(w1, w2))),
#                  data = data_mm)
data_mm[moved, "w1"] <- runif(length(moved), 0, 1)
data_mm[moved, "w2"] <- 1 - data_mm[moved, "w1"]

bv("data", "head(data_mm) after the weights are amended",
   print(utils::head(data_mm)), NA_character_, "")

fit_mm2 <- bv("model", "ML: fit_mm2", {
  frm(bf(y ~ 1 + (1 | mm(s1, s2, weights = cbind(w1, w2)))),
      data = data_mm, family = gaussian())
}, NA_character_, "")

bv("post", "ML: summary(fit_mm2)", summary(fit_mm2), NA_character_, "")

## ============ PATH 2: sampling (frm_sample) ============
#
# Every chain here is one chain of 400 iterations with 200 warmup,
# against the vignette's 4 x 2000 (2 x 2000 for the rent models) with
# 1000 warmup. That is deliberate and it costs accuracy: expect R-hat
# above 1.05 and low ESS on most of these fits. Where the short chain
# makes the vignette's own output unreachable, the row says so instead
# of claiming agreement.
#
# Route choice: the FIT route is used for the fits whose chains move
# freely from the ML mode, because it reuses the assembled objective.
# The FORMULA route is used for fit_loss1, where the vignette's priors
# are the point of the model, and for the two rent models, where a flat
# prior on a smooth's wiggle SD pins the chain at a boundary mode. The
# cost of the wrong choice is measured below.

s_zinb1 <- bv("model", "SAMPLE: fit_zinb1", {
  frm_sample(fit_zinb1, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, NA_character_, "")

bv("post", "SAMPLE: summary(fit_zinb1)", summary(s_zinb1), "BEHAVIOR",
   "the population-level estimates track brms but zi stays on the logit link as zi_Intercept, where brms reports zi = 0.41 on the probability scale")

bv("post", "SAMPLE: conditional_effects(fit_zinb1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(s_zinb1)
  print(names(ce))
  names(ce)
}, NA_character_, "")

bv("post", "SAMPLE: pp_check(fit_zinb1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(s_zinb1))
}, NA_character_, "")

s_zinb2 <- bv("model", "SAMPLE: fit_zinb2", {
  frm_sample(fit_zinb2, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, NA_character_, "")

# brms: LOO(fit_zinb1, fit_zinb2), printing 1639.52, 1621.35,
# difference 18.16 (15.71), read as "better fit but modest".
bv("post", "SAMPLE: LOO(fit_zinb1, fit_zinb2)", LOO(s_zinb1, s_zinb2),
   "REFUSAL",
   "LOO() is refused even on draws as a deprecated brms spelling; the message names loo() and loo::loo_compare(), which is exactly the replacement")

bv("post", "SAMPLE: loo_compare(loo(fit_zinb1), loo(fit_zinb2))", {
  loo_compare(loo(s_zinb1), loo(s_zinb2))
}, "BEHAVIOR",
"the two models cannot be separated here: elpd_diff is 0.8 with se 0.4 and two Pareto k values are flagged, where brms's 4000 draws gave 18.16 (15.71). The layout is loo_compare's, not brms's three-row LOOIC block")

s_rent1 <- bv("model", "SAMPLE: fit_rent1", {
  frm_sample(bf(rentsqm ~ t2(area, yearc) + (1 | district)), data = rent,
             family = gaussian(),
             prior = set_prior("cauchy(0,2)", class = "sd"),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "BEHAVIOR",
"frm_sample() reports the brms 2.23 default priors it supplied for the intercepts and says the slopes stay flat, which brms does not print at fit time")

# Why the formula route. brms's call carried no prior, so the fit route
# is the literal translation, and it takes 98 s here against 14 s for
# the formula route, ending with R-hat NA and every post-warmup
# transition at maximum treedepth. A smooth's wiggle SD has a flat
# prior on that route, so the chain sits on the boundary mode the ML
# fit handed it. The note frm_sample() prints names set_prior(class =
# "sd"), which is the fix. A shorter chain shows the same shape.
bv("model", "SAMPLE: fit_rent1 through the fit route, no priors", {
  frm_sample(fit_rent1, chains = 1, iter = 100, warmup = 50, seed = 1,
             cores = 1, refresh = 0)
}, "BEHAVIOR",
"the literal translation of brms's prior-free call pins the chain at a boundary mode: the largest R-hat is 2.1 on a quarter of the draws, and the full-length version of this run took 98 s against 14 s for the same model with a wiggle-SD prior")

bv("post", "SAMPLE: summary(fit_rent1)", summary(s_rent1), "BEHAVIOR",
   "the smoothing SDs appear as raw theta_1..n rows, not as brms's named sds(t2areayearc_1..3) block")

bv("post", "SAMPLE: conditional_effects(fit_rent1, surface = TRUE)", {
  conditional_effects(s_rent1, surface = TRUE)
}, "REFUSAL", "the draws method refuses the surface display for the same reason the ML method does")

bv("post", "SAMPLE: bayes_R2(fit_rent1)", bayes_R2(s_rent1), NA_character_, "")

s_rent2 <- bv("model", "SAMPLE: fit_rent2", {
  frm_sample(bf(rentsqm ~ t2(area, yearc) + (1 | ID1 | district),
                sigma ~ t2(area, yearc) + (1 | ID1 | district)),
             data = rent, family = gaussian(),
             prior = set_prior("cauchy(0,2)", class = "sd") +
               set_prior("lkj(2)", class = "cor"),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "SPELLING",
"the fit route is unusable here (283 s and a chain that never leaves the boundary mode), so the model is respecified through the formula route with a wiggle-SD prior brms's call did not carry")

# brms's reference for this block: sd(Intercept) 0.60,
# sd(sigma_Intercept) 0.11, cor(Intercept, sigma_Intercept) 0.72.
bv("post", "SAMPLE: VarCorr(fit_rent2) [the group-level block]", VarCorr(s_rent2),
   "BEHAVIOR",
   "the two SDs read 0.38 and 0.085 against brms's 0.60 and 0.11, but the 200-draw chain cannot settle the cross-formula correlation: 0.07 with an interval from -0.71 to 0.83, where brms reports 0.72 (0.35 to 0.98). The row layout is a data frame of grp/var1/var2 rather than brms's named sd() and cor() lines")

bv("post", "SAMPLE: conditional_smooths(fit_rent2)", conditional_smooths(s_rent2),
   "MISSING", "missing on both paths")

# The nonlinear model, through the FORMULA route because the vignette's
# priors are what identifies it. `start` is still needed: the priors do
# not move the starting point.
s_loss1 <- bv("model", "SAMPLE: fit_loss1", {
  frm_sample(bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
                ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE),
             data = loss, family = gaussian(),
             start = list(beta = c(5000, 1, 45)),
             prior = set_prior("normal(5000,1000)", class = "Intercept", dpar = "ult") +
               set_prior("normal(1,2)", class = "Intercept", dpar = "omega") +
               set_prior("normal(45,10)", class = "Intercept", dpar = "theta") +
               set_prior("cauchy(0,2)", class = "sd"),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "SPELLING",
"prior(..., nlpar = 'ult') becomes set_prior(..., class = 'Intercept', dpar = 'ult'), and start = is still required because the priors do not choose the starting point")

bv("post", "SAMPLE: summary(fit_loss1)", summary(s_loss1), "BEHAVIOR",
   "ult 5154 against brms's 5273.70 and omega 1.33 against 1.34, but the sampler reports 58 divergent transitions and a largest R-hat of 1.34, and the group-level SD brms reports as 745.74 appears only as the raw theta_1")

# The band the ML path had to bootstrap comes free from the draws.
bv("post", "SAMPLE: conditional_effects(fit_loss1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(s_loss1)
  print(names(ce))
  names(ce)
}, NA_character_, "")

bv("post", "SAMPLE: conditional_effects(conditions, re_formula = NULL, method = 'predict')", {
  conditional_effects(s_loss1, conditions = conditions,
                      re_formula = NULL, method = "predict")
}, "REFUSAL",
"the draws method has no method = at all; the message says to quantile posterior_predict() over your own grid, which is right but is not a one-line replacement")

bv("post", "SAMPLE: posterior_predict(fit_loss1) [the named substitute]", {
  pp <- posterior_predict(s_loss1)
  print(dim(pp))
  round(apply(pp, 2, stats::quantile, c(0.025, 0.5, 0.975))[, 1:4], 1)
}, "SPELLING", "the prediction interval has to be quantiled by hand")

s_mm <- bv("model", "SAMPLE: fit_mm", {
  frm_sample(fit_mm, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, NA_character_, "")

# brms's reference: sd(Intercept) 2.76, Intercept 19, sigma 3.58.
bv("post", "SAMPLE: summary(fit_mm)", summary(s_mm), "BEHAVIOR",
   "the Intercept reads 18.1 against brms's 19, but the multi-membership SD is reported as the raw theta_1 rather than as brms's sd(Intercept) row; VarCorr() is the place that names it")

bv("post", "SAMPLE: VarCorr(fit_mm)", VarCorr(s_mm), NA_character_, "")

bv("post", "SAMPLE: ngrps(fit_mm)", ngrps(s_mm), NA_character_, "")

s_mm2 <- bv("model", "SAMPLE: fit_mm2", {
  frm_sample(fit_mm2, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, NA_character_, "")

bv("post", "SAMPLE: summary(fit_mm2)", summary(s_mm2), NA_character_, "")

bv_done()
