# Translation of brms's "brms: An R Package for Bayesian Multilevel
# Models Using Stan" (brms 2.23.0, doc/brms_overview.ltx) onto the
# frmtmb surface.
#
# The vignette ships as pre-rendered JSS LaTeX. The code comes out of
# the Sinput environments and the printed reference output out of the
# Soutput environments, so every BEHAVIOR edge below names the number
# brms printed next to the number frmtmb prints.
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
bv_init("brms_overview")

options(mc.cores = 1)

## ============ PATH 1: ML / Laplace (frm) ============

## ---- A worked example ----------------------------------------------
#
# brms: library("brms"); data("kidney"); head(kidney, n = 3).
# The harness does not attach brms on purpose, so the data set is
# reached with the `brms::` prefix. Full size, 76 rows: no shrinking is
# needed because the fit takes under a second.

kidney <- brms::kidney
inhaler <- brms::inhaler

bv("data", "head(kidney, n = 3)", print(utils::head(kidney, n = 3)),
   NA_character_, "")

## ---- fit1 ----------------------------------------------------------
#
# brms:
#   fit1 <- brm(formula = time | cens(censored) ~ age * sex + disease
#               + (1 + age|patient),
#               data = kidney, family = lognormal(),
#               prior = c(set_prior("normal(0,5)", class = "b"),
#                         set_prior("cauchy(0,2)", class = "sd"),
#                         set_prior("lkj(2)", class = "cor")),
#               warmup = 1000, iter = 2000, chains = 4,
#               control = list(adapt_delta = 0.95))
#
# The three priors are weakly informative regularizers that keep the
# sampler out of the tails; the model is identified without them. They
# go, together with warmup/iter/chains/control. What is left is the
# mechanical transform: brm -> frm plus argument removal. The formula
# is wrapped in bf() for house style only; the bare formula is also
# accepted. (The v0.34 audit recorded a missing gaussian default as
# FN-1; frm() now defaults, and this model names its family anyway.)
fit1 <- bv("model", "ML: fit1", {
  frm(bf(time | cens(censored) ~ age * sex + disease + (1 + age | patient)),
      data = kidney, family = lognormal())
}, NA_character_, "")

# THE PRIOR INTERFACE, measured separately. The same three priors,
# spelled on frmtmb's current surface. They are accepted, but they are
# penalties on the likelihood and not densities in a posterior, so the
# fit is a MAP estimate. This is the one place where the brms number
# and the frmtmb number agree BECAUSE of the prior: with no penalty
# sd(patient) collapses to 3e-4 (a boundary maximum), and with
# cauchy(0,2) it comes back to 0.38, against brms's posterior mean of
# 0.40.
bv("model", "ML: fit1 with the vignette's priors kept", {
  pr <- set_prior("normal(0,5)", class = "b") +
    set_prior("cauchy(0,2)", class = "sd") +
    set_prior("lkj(2)", class = "cor")
  fp <- frm(bf(time | cens(censored) ~ age * sex + disease +
                 (1 + age | patient)),
            data = kidney, family = lognormal(), priors = pr)
  print(VarCorr(fp))
  fp
}, "BEHAVIOR",
"set_prior() takes brms's own strings but adds a penalty, so the fit is MAP: sd(patient) reads 0.38 against brms's posterior mean 0.40, and 3e-4 without the penalty")

## ---- prior: the standalone specification ---------------------------
#
# brms:
#   prior <- c(set_prior("normal(0,10)", class = "b", coef = "age"),
#              set_prior("cauchy(1,2)", class = "b", coef = "sexfemale"))
#
# The line ports character for character. The object it builds means a
# penalty rather than a prior, which is why the model call above is
# labeled and this one is not: the spelling itself is clean.
bv("post", "ML: prior <- c(set_prior(...), set_prior(...))", {
  c(set_prior("normal(0,10)", class = "b", coef = "age"),
    set_prior("cauchy(1,2)", class = "b", coef = "sexfemale"))
}, NA_character_, "")

# The vectorized shortcut named in the same paragraph.
bv("post", "ML: set_prior('normal(0,10)', class = 'b')", {
  set_prior("normal(0,10)", class = "b")
}, NA_character_, "")

# brms's shrinkage prior for sparse population-level effects.
bv("post", "ML: set_prior('horseshoe(1)')", {
  set_prior("horseshoe(1)")
}, "MISSING",
"only normal, student_t, cauchy, exponential and lkj parse; the message lists them but names no shrinkage substitute")

# The footnote's get_prior() call, verbatim.
bv("post", "ML: get_prior(time | cens(censored) ~ ..., family = lognormal())", {
  get_prior(bf(time | cens(censored) ~ age * sex + disease +
                 (1 + age | patient)),
            data = kidney, family = lognormal())
}, "BEHAVIOR",
"the row set matches brms but every prior column reads '(flat)' and there is no source column, because frmtmb has no default priors to report on this path")

## ---- Analyzing the results -----------------------------------------

# brms's stancode(fit1) and standata(fit1), named in the same paragraph
# as summary().
bv("post", "ML: stancode(fit1)", stancode(fit1), "MISSING",
   "stancode() is a frmtmb_draws method; a frmtmb_fit holds a TMB objective and no Stan program")

bv("post", "ML: standata(fit1)", standata(fit1), "MISSING",
   "same as stancode(); frm_sample() produces the draws object these methods want")

# brms prints the WAIC on the header line when waic = TRUE (673.51 in
# the vignette). frmtmb absorbs the argument into ... and prints AIC and
# BIC instead.
bv("post", "ML: summary(fit1, waic = TRUE)", summary(fit1, waic = TRUE),
   "BEHAVIOR",
   "waic = is silently absorbed; the header carries AIC 679.13 and BIC 704.77 where brms printed WAIC 673.51, and sigma is on its LOG link scale (0.114) where brms reports the response-scale 1.15")

# brms draws trace and density plots per parameter (Figure kidney_plot).
# There is no chain to trace, so plot.frmtmb_fit draws the two
# regression diagnostics.
bv("post", "ML: plot(fit1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit1)
}, "BEHAVIOR", "brms draws trace/density per parameter; frmtmb draws residual and QQ diagnostics")

# Figure kidney_conditional_effects. The four effects brms plots
# (age, sex, disease, age:sex) all come back.
bv("post", "ML: conditional_effects(fit1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(fit1)
  print(names(ce))
  plot(ce, ask = FALSE)
  names(ce)
}, NA_character_, "")

# brms's older name for the same method, which this vignette's era
# still used.
bv("post", "ML: marginal_effects(fit1)", marginal_effects(fit1), "MISSING",
   "brms's pre-2.10 alias is not re-exported; conditional_effects() is the replacement but nothing says so")

# "An even more detailed investigation can be achieved by applying the
# shinystan package through method launch_shiny."
bv("post", "ML: launch_shiny(fit1)", launch_shiny(fit1), "MISSING",
   "shinystan inspects chains; a point estimate has none, and no alternative is named")

## ---- hypothesis ----------------------------------------------------
#
# brms: hypothesis(fit1, "Intercept - age > 0", class = "sd",
#                  group = "patient")
# printing Estimate 0.39, l-95% CI 0.03, Evid.Ratio 67.97.
#
# The line runs verbatim, including the one-sided operator and the
# class/group narrowing. (The v0.34 audit recorded the one-sided form
# as FN-5, a rejection; it is fixed.) What comes back is a frequentist
# test.
bv("post", "ML: hypothesis(fit1, 'Intercept - age > 0', class = 'sd', group = 'patient')", {
  hypothesis(fit1, "Intercept - age > 0", class = "sd", group = "patient")
}, "BEHAVIOR",
"the last column is a one-sided p value, not an Evid.Ratio; and because both group-level SDs sit at the boundary under ML the estimate is 3e-4 where brms reports 0.39")

## ---- fit2 ----------------------------------------------------------
#
# brms: fit2 <- update(fit1, formula. = ~ . - (1 + age|patient)
#                      + (1|patient))
# The delta formula and the `formula.` argument name are both accepted.
# (The v0.34 audit recorded both as FN-9; they are fixed.)
fit2 <- bv("model", "ML: fit2", {
  update(fit1, formula. = ~ . - (1 + age | patient) + (1 | patient))
}, NA_character_, "")

bv("post", "ML: formula(fit2)", print(stats::formula(fit2)), NA_character_, "")

# brms: LOO(fit1, fit2), printing LOOIC 675.45 and 674.17.
bv("post", "ML: LOO(fit1, fit2)", LOO(fit1, fit2), "MISSING",
   "LOO() is a frmtmb_draws method; on a point estimate the frequentist route is AIC()/BIC()/anova()")

bv("post", "ML: AIC(fit1, fit2) [the substitute]", {
  print(stats::AIC(fit1, fit2))
}, "SPELLING",
"the LOOIC table becomes an AIC table: 679.13 against 675.13, the same ordering brms found, but with no standard error on the difference")

## ---- Modeling ordinal data -----------------------------------------

bv("data", "head(inhaler, n = 1)", print(utils::head(inhaler, n = 1)),
   NA_character_, "")

# brms: fit3 <- brm(formula = rating ~ treat + period + carry
#                   + (1|subject), data = inhaler, family = cumulative)
# The bare family constructor, without parentheses, is accepted.
# (The v0.34 audit recorded this as FN-6, a rejection; it is fixed.)
fit3 <- bv("model", "ML: fit3", {
  frm(bf(rating ~ treat + period + carry + (1 | subject)),
      data = inhaler, family = cumulative)
}, NA_character_, "")

bv("post", "ML: summary(fit3)", summary(fit3), "BEHAVIOR",
   "the K-1 thresholds print as 'Family parameters (tau_raw, internal scale)' in a (first threshold, log increments) parameterization, not as brms's Intercept[1..3] on the latent scale")

## ---- fit4 ----------------------------------------------------------
#
# brms:
#   fit4 <- brm(formula = rating ~ period + carry + cs(treat)
#               + (1|subject),
#               data = inhaler, family = sratio,
#               threshold = "equidistant",
#               prior = set_prior("normal(-1,2)", coef = "treat"))
#
# The vignette's exact call. `threshold = "equidistant"` has no frmtmb
# argument at all, so it lands as an unused argument rather than as a
# refusal that names the supported parameterization.
bv("model", "ML: fit4 (the vignette's exact call)", {
  frm(bf(rating ~ period + carry + cs(treat) + (1 | subject)),
      data = inhaler, family = sratio, threshold = "equidistant",
      priors = set_prior("normal(-1,2)", coef = "treat"))
}, "MISSING",
"threshold = 'equidistant' has no frmtmb spelling; the error is R's generic 'unused argument' and points nowhere")

# The prior alone, with the threshold restriction removed, so the two
# gaps are measured apart. A cs() coefficient is not addressable by
# class "b" and coefficient name: it lives in its own bcs2 block.
bv("model", "ML: fit4 prior on the cs(treat) coefficient", {
  frm(bf(rating ~ period + carry + cs(treat) + (1 | subject)),
      data = inhaler, family = sratio(),
      priors = set_prior("normal(-1,2)", coef = "treat"))
}, "MISSING",
"a category-specific coefficient has no prior target: the message says 'Prior target not found (class=b, coef=treat)' and names no class that would reach it")

# What is left after both gaps: the model with free thresholds and no
# prior. This is the only form that runs today.
fit4 <- bv("model", "ML: fit4 (threshold and prior dropped)", {
  frm(bf(rating ~ period + carry + cs(treat) + (1 | subject)),
      data = inhaler, family = sratio())
}, "SPELLING",
"threshold = and the informative prior are both deleted, so the thresholds are free and brms's delta parameter has nothing to correspond to")

bv("post", "ML: summary(fit4, waic = TRUE)", summary(fit4, waic = TRUE),
   "BEHAVIOR",
   "brms prints Intercept[1..3], treat[1..3] and delta; frmtmb prints period and carry, then tau_raw_1..3 and bcs2_1..3 on their internal scales, and the third threshold runs off to -18 because nothing holds the thresholds equidistant")

# brms: plot(fit4), Figure inhaler_plot.
bv("post", "ML: plot(fit4)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit4)
}, "BEHAVIOR", "diagnostics, not traces")

## ============ PATH 2: sampling (frm_sample) ============
#
# Every chain here is one chain of 400 iterations with 200 warmup,
# against the vignette's 4 x 2000 with 1000 warmup. That is deliberate
# and it costs accuracy: expect R-hat above 1.05 and low ESS on most of
# these fits. Where the short chain makes the vignette's own output
# unreachable, the row says so instead of claiming agreement.
#
# Route choice: fit1 goes through the FORMULA route, because the
# vignette's priors are the whole point of the sampled fit and only the
# formula route takes `priors =`. The next row shows what the fit route
# does without them.

s_fit1 <- bv("model", "SAMPLE: fit1", {
  frm_sample(bf(time | cens(censored) ~ age * sex + disease +
                  (1 + age | patient)),
             data = kidney, family = lognormal(),
             priors = set_prior("normal(0,5)", class = "b") +
               set_prior("cauchy(0,2)", class = "sd") +
               set_prior("lkj(2)", class = "cor"),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, NA_character_, "")

# The FIT route on the same model, to show why the formula route was
# chosen. frm_sample() initializes at the ML mode, and the ML mode of
# this model is singular, so with flat variance priors the chain never
# leaves the boundary. frmtmb warns about exactly this and names
# set_prior(class = "sd") and set_prior(class = "cor").
bv("model", "SAMPLE: fit1 through the fit route, no priors", {
  sf <- frm_sample(fit1, chains = 1, iter = 200, warmup = 100, seed = 1,
                   cores = 1, refresh = 0)
  print(VarCorr(sf))
  sf
}, "BEHAVIOR",
"without priors the chain stays pinned at the singular ML mode: sd(patient) reads 3.1e-4 with an interval 3e-7 wide and R-hat reaches 1.87; the note names the set_prior() classes that fix it")

# brms: summary(fit1, waic = TRUE). The reference numbers are
# sd(Intercept) 0.40, sd(age) 0.01, cor -0.13, sexfemale 2.42,
# sigma 1.15.
bv("post", "SAMPLE: summary(fit1, waic = TRUE)", summary(s_fit1, waic = TRUE),
   "BEHAVIOR",
   "the estimates land on brms (sexfemale 2.29 against 2.42) but waic = is absorbed and no WAIC line is printed, the group-level block appears as raw theta_1..3 instead of brms's sd()/cor() rows (VarCorr() carries those), and the poor n_eff comes from the 200-draw chain and not from a porting gap")

bv("post", "SAMPLE: VarCorr(fit1)", VarCorr(s_fit1), NA_character_, "")

# brms's Figure kidney_plot is trace and density per parameter.
bv("post", "SAMPLE: plot(fit1) [trace and density]", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(mcmc_plot(s_fit1))
}, "SPELLING",
"plot() on a draws object is not the brms trace panel; mcmc_plot() is the bayesplot entry point that draws it")

bv("post", "SAMPLE: conditional_effects(fit1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(s_fit1)
  print(names(ce))
  names(ce)
}, NA_character_, "")

# The draws method has no method = argument at all. brms's
# method = "predict" for a predictive band has to be built by hand.
bv("post", "SAMPLE: conditional_effects(fit1, method = 'predict')", {
  conditional_effects(s_fit1, method = "predict")
}, "REFUSAL",
"the draws method has no method = or band =; the message says to quantile posterior_predict() over your own grid, which is right but is not a one-line replacement")

bv("post", "SAMPLE: pp_check(fit1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(s_fit1))
}, NA_character_, "")

# The one-sided hypothesis, now against a posterior. brms reports
# Estimate 0.39, l-95% CI 0.03, Evid.Ratio 67.97.
bv("post", "SAMPLE: hypothesis(fit1, 'Intercept - age > 0', class = 'sd', group = 'patient')", {
  hypothesis(s_fit1, "Intercept - age > 0", class = "sd", group = "patient")
}, "BEHAVIOR",
"the estimate and lower bound now match brms (0.44 and 0.07 against 0.39 and 0.03) but the reported quantity is a posterior tail probability, not brms's Evid.Ratio")

bv("post", "SAMPLE: bayes_R2(fit1)", bayes_R2(s_fit1), NA_character_, "")

bv("post", "SAMPLE: posterior_summary(fit1)", posterior_summary(s_fit1),
   NA_character_, "")

# fit2 as its own sampled model. The fit route would be the natural
# choice here (no priors in brms's update() call), but the reduced
# model still has a boundary variance, so the formula route with the
# vignette's own sd prior is used and the divergence is recorded above.
s_fit2 <- bv("model", "SAMPLE: fit2", {
  frm_sample(bf(time | cens(censored) ~ age * sex + disease +
                  (1 | patient)),
             data = kidney, family = lognormal(),
             priors = set_prior("normal(0,5)", class = "b") +
               set_prior("cauchy(0,2)", class = "sd"),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, NA_character_, "")

# brms: LOO(fit1, fit2), printing 675.45, 674.17, difference 1.28 (0.99).
bv("post", "SAMPLE: LOO(fit1, fit2)", LOO(s_fit1, s_fit2), "REFUSAL",
   "LOO() is refused even on draws as a deprecated brms spelling; the message names loo() and loo::loo_compare(), which is exactly the replacement")

bv("post", "SAMPLE: loo(fit1); loo(fit2) [the LOOIC values]", {
  print(loo(s_fit1))
  loo(s_fit2)
}, "SPELLING",
"looic 674.7 against brms's 675.45, but a third of the Pareto k values are flagged bad on 200 draws, which brms's 4000 draws did not show")

bv("post", "SAMPLE: loo_compare(loo(fit1), loo(fit2))", {
  loo_compare(loo(s_fit1), loo(s_fit2))
}, "BEHAVIOR",
"elpd_diff 1.2 (se 0.9) reproduces brms's 1.28 (0.99), but the table is loo_compare's layout with diagnostic flags, not brms's three-row LOOIC block")

bv("post", "SAMPLE: waic(fit1)", waic(s_fit1), "BEHAVIOR",
   "waic 674.1 against brms's summary header of 673.51, with a p_waic warning that the 200-draw chain causes")

# The ordinal models.
s_fit3 <- bv("model", "SAMPLE: fit3", {
  frm_sample(bf(rating ~ treat + period + carry + (1 | subject)),
             data = inhaler, family = cumulative(),
             priors = set_prior("normal(0,5)", class = "b") +
               set_prior("cauchy(0,2)", class = "sd"),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "BEHAVIOR",
"frm_sample() reports that it has no default prior for this family's thresholds, where brms priors them under its Intercept class")

bv("post", "SAMPLE: summary(fit3)", summary(s_fit3), "BEHAVIOR",
   "the thresholds appear as theta_1 and tau_raw_1..3 on internal scales, not as brms's Intercept[1..3]")

# brms's fit4 keeps the informative prior on cs(treat). On the sampling
# path the prior has to be widened to the whole b class, because a
# category-specific coefficient still has no target of its own.
s_fit4 <- bv("model", "SAMPLE: fit4", {
  frm_sample(bf(rating ~ period + carry + cs(treat) + (1 | subject)),
             data = inhaler, family = sratio(),
             priors = set_prior("normal(-1,2)", class = "b") +
               set_prior("cauchy(0,2)", class = "sd"),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "SPELLING",
"prior(normal(-1,2), coef = 'treat') has to become a class-wide b prior, and threshold = 'equidistant' is still dropped")

bv("post", "SAMPLE: summary(fit4, waic = TRUE)", summary(s_fit4, waic = TRUE),
   "BEHAVIOR",
   "the cs(treat) coefficients bcs2_1..3 read -0.93, -1.04, -0.68 against brms's treat[1..3] of -0.96, -0.65, -2.65, but tau_raw_3 runs to -70 with R-hat 1.5 because nothing holds the thresholds equidistant; brms's delta has no counterpart")

bv_done()
