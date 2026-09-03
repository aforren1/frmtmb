# Translation of brms's "Estimating Monotonic Effects with brms"
# (brms 2.23.0, doc/brms_monotonic.Rmd) onto the frmtmb surface.
#
# Read the comments, not just the code: every `<EDGE>` marker on a bv()
# call is a place where the brms line did not carry over verbatim. The
# reference for a BEHAVIOR edge is doc/brms_monotonic.html, which holds
# brms's own printed output.
#
# The vignette is translated TWICE, once onto each inference path.
# PATH 1 is frm(), maximum likelihood with a Laplace approximation for
# the random effects. PATH 2 is frm_sample(), which hands the same TMB
# objective to tmbstan and gives back draws. This vignette is the harder
# case for PATH 2: mo()'s simplex has no prior on the frmtmb side, and
# the section below shows what that costs.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_monotonic")

options(mc.cores = 1)
set.seed(1234)

## ---- shared data ----------------------------------------------------
#
# The vignette's own sizes (100 people, 10 cities) are small enough to
# keep. The data is built once here because the vignette mutates `ls`
# between fit5 and fit6, and both inference paths must fit fit1 to fit5
# on the pre-mutation response and fit6 on the post-mutation one.

income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(income_options, 100, TRUE),
                 levels = income_options, ordered = TRUE)
mean_ls <- c(30, 60, 70, 75)
ls <- mean_ls[income] + rnorm(100, sd = 7)
dat <- data.frame(income, ls)
dat$income_num <- as.numeric(dat$income)
dat$age <- rnorm(100, mean = 40, sd = 10)

# The vignette's "assume income is an unordered factor" step. It is a
# contrast assignment on the ordered factor, so it also changes what
# mo() sees; a separate frame keeps the monotonic fits away from it.
dat_unord <- dat
contrasts(dat_unord$income) <- contr.treatment(4)

dat_city <- dat
dat_city$city <- rep(1:10, each = 10)
var_city <- rnorm(10, sd = 10)
dat_city$ls <- dat_city$ls + var_city[dat_city$city]

## ============ PATH 1: ML / Laplace (frm) ============

## ---- A simple monotonic model --------------------------------------

fit1 <- bv("model", "ML: fit1", {
  frm(bf(ls ~ mo(income)), data = dat, family = gaussian())
}, NA_character_, "")

# brms writes `brm(ls ~ mo(income), data = dat)` with no family at all.
# The v0.34 audit recorded that as FN-1, the single highest-leverage
# failure in the whole port. It is fixed: the bare line now runs, and
# the family argument above is house style rather than a requirement.
bv("post", "ML: the brms line with no family at all", {
  alt <- frm(ls ~ mo(income), data = dat)
  stopifnot(all.equal(unlist(fixef(alt)), unlist(fixef(fit1)),
                      tolerance = 1e-8))
  "identical"
}, NA_character_, "")

# brms prints a "Monotonic Simplex Parameters" block with the three
# elements of moincome1 on the simplex itself (0.65, 0.26, 0.09). frmtmb
# prints two unconstrained zeta values instead, so the reader cannot see
# the quantity the whole vignette is about without transforming them.
bv("post", "ML: summary(fit1)", summary(fit1), "BEHAVIOR",
   "brms reports the 3-element simplex moincome1; frmtmb reports 2 zeta parameters on their internal unconstrained scale, under the heading 'Family parameters', and never prints the simplex")

# brms draws the posterior densities of the simplex here. There is no
# chain on this path, and `variable` and `regex` are absorbed.
bv("post", "ML: plot(fit1, variable = 'simo', regex = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit1, variable = "simo", regex = TRUE)
}, "BEHAVIOR", "variable = and regex = are absorbed by ...; the panels are residual and QQ diagnostics, not the simplex densities the vignette's paragraph refers to")

# The v0.34 audit recorded this as FN-3: conditional_effects() found no
# plottable predictor on a fit whose only fixed term is mo(). Fixed.
bv("post", "ML: plot(conditional_effects(fit1))", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(fit1))
}, NA_character_, "")

## ---- The two comparison models -------------------------------------

fit2 <- bv("model", "ML: fit2", {
  frm(bf(ls ~ income_num), data = dat, family = gaussian())
}, NA_character_, "")

bv("post", "ML: summary(fit2)", summary(fit2), NA_character_, "")

fit3 <- bv("model", "ML: fit3", {
  frm(bf(ls ~ income), data = dat_unord, family = gaussian())
}, NA_character_, "")

bv("post", "ML: summary(fit3)", summary(fit3), NA_character_, "")

bv("post", "ML: loo(fit1, fit2, fit3)", loo(fit1, fit2, fit3),
   "MISSING", "loo() is a frmtmb_draws method; on a point fit the frequentist route is AIC()/BIC(), and PATH 2 below recovers the real thing")

bv("post", "ML: AIC(fit1, fit2, fit3) [the substitute]", {
  print(AIC(fit1, fit2, fit3))
}, "SPELLING", "the elpd comparison becomes AIC; the ordering agrees with brms (fit1 and fit3 tied, fit2 far worse) but there is no elpd_diff, no se_diff and no Pareto k")

## ---- Setting prior distributions -----------------------------------
#
# This is the one section whose priors are NOT an MCMC convenience: the
# dirichlet on the simplex IS the subject, so it is translated rather
# than dropped.

bv("post", "ML: prior4 = prior(dirichlet(c(2, 1, 1)), class = 'simo', coef = 'moincome1')", {
  set_prior("dirichlet(c(2, 1, 1))", class = "simo", coef = "moincome1")
}, "MISSING", "set_prior() has no dirichlet; the parse error lists the five supported densities and none of them is a density over a simplex")

# The same call with a density set_prior() does know, to show that the
# class is refused independently of the distribution.
bv("post", "ML: set_prior(class = 'simo') with a supported density", {
  set_prior("normal(0, 1)", class = "simo", coef = "moincome1")
}, "MISSING", "class = 'simo' is not among the five accepted classes; the match.arg message lists them but does not say that monotonic simplexes cannot be given a prior at all")

bv("model", "ML: fit4 [the vignette's exact translation]", {
  frm(bf(ls ~ mo(income)), data = dat, family = gaussian(),
      priors = set_prior("dirichlet(c(2, 1, 1))", class = "simo",
                         coef = "moincome1"))
}, "MISSING", "the model cannot be built at all, because its prior cannot be constructed; sample_prior = TRUE is MCMC-only and is dropped separately")

# The section still needs a fit to talk about, so fit4 is the flat one.
# It is the same model as fit1, which is itself the finding: on this
# path there is no way to express the section's assumption.
fit4 <- bv("model", "ML: fit4 [without the prior]", {
  frm(bf(ls ~ mo(income)), data = dat, family = gaussian())
}, "MISSING", "the dirichlet(c(2, 1, 1)) that the section exists to demonstrate has no spelling, so fit4 is fit1 again")

bv("post", "ML: summary(fit4)", summary(fit4), "BEHAVIOR",
   "identical to summary(fit1) because the prior was dropped; brms's fit4 differs from its fit1 in exactly the simplex block frmtmb does not print")

bv("post", "ML: plot(fit4, variable = 'prior_simo', regex = TRUE, N = 3)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit4, variable = "prior_simo", regex = TRUE, N = 3)
}, "MISSING", "there is no prior to draw from and no draws to plot; the call is silently absorbed and draws diagnostics instead of refusing")

## ---- Interactions of monotonic variables ---------------------------

fit5 <- bv("model", "ML: fit5", {
  frm(bf(ls ~ mo(income) * age), data = dat, family = gaussian())
}, NA_character_, "")

# brms estimates a SECOND simplex for the interaction (moincome:age1),
# so its summary carries six simplex rows. frmtmb shares one simplex
# between the main effect and the interaction, which is a different
# model, not a different printout.
bv("post", "ML: summary(fit5)", summary(fit5), "BEHAVIOR",
   "brms fits two simplexes (moincome1 and moincome:age1); frmtmb makes the interaction share the main effect's simplex, so the shape of the monotonic effect cannot vary with age")

bv("post", "ML: conditional_effects(fit5, 'income:age')", {
  conditional_effects(fit5, "income:age")
}, NA_character_, "")

## ---- Monotonic group-level effects ---------------------------------

bv("model", "ML: fit6", {
  frm(bf(ls ~ mo(income) * age + (mo(income) | city)), data = dat_city,
      family = gaussian())
}, "REFUSAL", "group-level mo() is refused; the message names the supported forms (standalone term, two-way interaction) and so points at the workaround below")

fit6 <- bv("model", "ML: fit6 [the varying-intercept workaround]", {
  frm(bf(ls ~ mo(income) * age + (1 | city)), data = dat_city,
      family = gaussian())
}, "SPELLING", "(mo(income) | city) becomes (1 | city); the city intercepts are kept and the varying monotonic effect that the section is about is lost")

bv("post", "ML: summary(fit6)", summary(fit6), "BEHAVIOR",
   "no sd(moincome) and no cor(Intercept, moincome) rows, because the varying slope is not in the model")

## ---- ML path: named by the audit brief, absent from this vignette ---
#
# brms 2.23.0's brms_monotonic.Rmd fits nothing ordinal, so the
# threshold = argument the audit brief names is exercised here on the
# same data. These rows are labeled so the vignette tally stays honest.

dat$ls_cat <- as.integer(cut(dat$ls, 3))

bv("model", "ML: appendix cumulative() with threshold = 'equidistant'", {
  frm(bf(ls_cat ~ mo(income)), data = dat, family = cumulative(),
      threshold = "equidistant")
}, "MISSING", "frm() has no threshold argument, so the ordinal threshold parameterizations have no spelling; the message is R's generic unused-argument error and names no alternative")

bv("model", "ML: appendix family = cumulative [bare constructor]", {
  frm(bf(ls_cat ~ mo(income)), data = dat, family = cumulative)
}, NA_character_, "")

## ============ PATH 2: sampling (frm_sample) ============
#
# The FIT route is used for every model except fit4. The vignette sets
# no prior on fit1, fit2, fit3, fit5 or fit6, so there is nothing for
# the formula route's brms-like defaults to carry, and starting from the
# frm() fit reuses the objective PATH 1 already built. fit4 IS a prior
# model, so it takes the formula route with a priors = argument, which
# is where the dirichlet gap shows up a second time.
#
# Chain settings are fixed for every call below: one chain, 400
# iterations with 200 warmup, seed 1. brms runs 4 chains of 2000.

SAMP <- function(fit) {
  frm_sample(fit, chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}

## ---- A simple monotonic model --------------------------------------

# The fit route starts the chain at the ML mode, which is inside the
# well-behaved part of the zeta surface. That is what keeps this run
# tame; the formula route below starts elsewhere and shows what the
# missing simplex prior actually costs.
s_fit1 <- bv("model", "SAMPLE: fit1", SAMP(fit1), "SPELLING",
   "brm(...) becomes frm() then frm_sample(fit, chains =, iter =, warmup =, seed =); the MCMC arguments come back on the second call, not the first")

bv("post", "SAMPLE: summary(fit1)", summary(s_fit1), "BEHAVIOR",
   "the regression coefficients agree with brms (Intercept 30.2 against 30.01, moincome 15.3 against 15.73), but the simplex rows are zeta1_1 and zeta1_2 on their internal unconstrained scale and not brms's moincome1[1], [2], [3] on the simplex, so the quantity the vignette interprets is still not printed")

# There is no plot method for a draws object, so plot() falls through to
# plot.default and dies inside base graphics.
bv("post", "SAMPLE: plot(fit1, variable = 'simo', regex = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(s_fit1, variable = "simo", regex = TRUE)
}, "MISSING", "no plot method for frmtmb_draws; the base-graphics message about list components names nothing, and mcmc_plot() is the replacement")

bv("post", "SAMPLE: mcmc_plot(fit1) [the substitute]", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(mcmc_plot(s_fit1))
}, "SPELLING", "plot(fit, variable =, regex =) becomes mcmc_plot(); this is the call that draws what brms's plot() draws")

bv("post", "SAMPLE: plot(conditional_effects(fit1))", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s_fit1))
}, NA_character_, "")

## ---- The two comparison models -------------------------------------

s_fit2 <- bv("model", "SAMPLE: fit2", SAMP(fit2), "SPELLING",
   "brm(...) becomes frm() then frm_sample(fit, chains =, iter =, warmup =, seed =)")

s_fit3 <- bv("model", "SAMPLE: fit3", SAMP(fit3), "SPELLING",
   "the two-call spelling again")

# The ML path's one MISSING post-processing gap, recovered.
bv("post", "SAMPLE: loo(fit1, fit2, fit3)", {
  loo_compare(loo(s_fit1), loo(s_fit2), loo(s_fit3))
}, "SPELLING", "brms's loo(fit1, fit2, fit3) becomes loo_compare(loo(a), loo(b), loo(c)); the elpd_diff ordering reproduces brms's, and the printed table adds frmtmb's p_worse and diag_ columns")

## ---- Setting prior distributions -----------------------------------
#
# The formula route, because this is the one model in the vignette that
# exists to carry a prior.

bv("model", "SAMPLE: fit4 [the vignette's exact translation]", {
  frm_sample(bf(ls ~ mo(income)), data = dat, family = gaussian(),
             priors = set_prior("dirichlet(c(2, 1, 1))", class = "simo",
                                coef = "moincome1"),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "MISSING", "the same dirichlet and simo gap as on the ML path, and it costs more here, as the next call shows")

# THE FINDING OF THIS SECTION. The formula route announces the brms 2.23
# defaults it supplies, and the monotonic simplex is not among them.
# brms's implicit dirichlet(1) is what makes that simplex a proper
# posterior; without it the zeta direction is flat and the chain walks
# it. Measured on this run: 14 divergent transitions, largest Rhat 1.65.
s_fit4 <- bv("model", "SAMPLE: fit4 [formula route, default priors]", {
  frm_sample(bf(ls ~ mo(income)), data = dat, family = gaussian(),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "BEHAVIOR", "the announced defaults cover the intercepts and the standard deviations but not the monotonic simplex, and the chain then reports 14 divergent transitions and Rhat 1.65 where brms's implicit dirichlet(1) would hold the simplex proper")

bv("post", "SAMPLE: summary(fit4)", summary(s_fit4), "BEHAVIOR",
   "zeta1_2 has a posterior mean near -2e6 with n_eff 3 and Rhat 1.80, which is the flat zeta direction and not a posterior; brms's own fit4 differs from its fit1 in exactly this block, and here neither model can express the difference")

## ---- Interactions of monotonic variables ---------------------------

s_fit5 <- bv("model", "SAMPLE: fit5", SAMP(fit5), "SPELLING",
   "the two-call spelling again")

bv("post", "SAMPLE: summary(fit5)", summary(s_fit5), "BEHAVIOR",
   "two zeta rows where brms prints six simplex rows, because the interaction shares the main effect's simplex instead of getting its own")

bv("post", "SAMPLE: conditional_effects(fit5, 'income:age')", {
  conditional_effects(s_fit5, "income:age")
}, NA_character_, "")

## ---- Monotonic group-level effects ---------------------------------

# The refusal happens in frm(), before anything is sampled, so PATH 2
# meets it in the same place PATH 1 did.
bv("model", "SAMPLE: fit6", {
  frm_sample(bf(ls ~ mo(income) * age + (mo(income) | city)),
             data = dat_city, family = gaussian(),
             chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}, "REFUSAL", "the same group-level mo() refusal, raised while the formula is parsed and before any sampling starts")

## ---- SAMPLE path: what the ML path could not reach ------------------

bv("post", "SAMPLE: waic(fit1)", waic(s_fit1), NA_character_, "")

bv("post", "SAMPLE: bayes_R2(fit1)", bayes_R2(s_fit1), NA_character_, "")

bv("post", "SAMPLE: as_draws_df(fit1)", {
  utils::head(as_draws_df(s_fit1), 3)
}, NA_character_, "")

bv_done()
