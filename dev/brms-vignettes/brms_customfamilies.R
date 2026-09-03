# Translation of brms's "Define Custom Response Distributions with brms"
# (brms 2.23.0, doc/brms_customfamilies.Rmd) onto the frmtmb surface.
#
# The whole vignette is one worked example: define a beta-binomial by
# hand and then make the post-processing methods work for it. brms
# builds it out of Stan source; frmtmb builds it out of an R
# log-density. That single difference reshapes almost every line, so
# read the comments for what each brms step turns into.
#
# Read the comments, not just the code: every `<EDGE>` marker on a bv()
# call is a place where the brms line did not carry over verbatim.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_customfamilies")

options(mc.cores = 1)
set.seed(1234)

## ---- A case study ---------------------------------------------------

bv("data", "data(cbpp, package = 'lme4'); head(cbpp)", {
  data("cbpp", package = "lme4", envir = globalenv())
  print(utils::head(cbpp))
}, NA_character_, "")

## ============ PATH 1: ML / Laplace (frm) ============

fit1 <- bv("model", "ML: fit1", {
  frm(bf(incidence | trials(size) ~ period + (1 | herd)),
      data = cbpp, family = binomial())
}, NA_character_, "")

# brms reports posterior means -1.40, -1.00, -1.14, -1.62 with
# sd(herd) = 0.76. The fixed effects land on top of those. sd(herd)
# comes back at 0.64: brms's half-t prior holds the variance component
# up and maximum likelihood does not.
bv("post", "ML: summary(fit1)", summary(fit1), "BEHAVIOR",
   "fixed effects match brms's posterior means; sd(herd) is 0.64 against brms's 0.76 because ML has no prior holding the variance component up")

## ---- Fitting custom family models -----------------------------------

# brms's call, verbatim:
#
#   beta_binomial2 <- custom_family(
#     "beta_binomial2", dpars = c("mu", "phi"),
#     links = c("logit", "log"),
#     lb = c(0, 0), ub = c(1, NA),
#     type = "int", vars = "vint1[n]"
#   )
#
# What carries over: the name, the dpars, and the links (as a named
# list, not a bare vector). What does not:
#   lb / ub    the bounds are implied by the links here, so there is
#              nothing to declare. frmtmb's own lower =/upper = are
#              frm() arguments over coefficients, not family metadata.
#   type       "int" becomes "discrete". frmtmb's `type` has four
#              values and also selects what fitted() returns.
#   vars       there is no declaration of extra variables. The lpdf
#              reads `aterms$vint1` directly, so vint(size) in the
#              formula is the whole mechanism.
# What is added: `lpdf`, which is the point. brms has no argument for
# it because the density lives in Stan.
beta_binomial2 <- bv("post", "ML: beta_binomial2 <- custom_family(...)", {
  custom_family(
    "beta_binomial2", dpars = c("mu", "phi"),
    links = list(mu = "logit", phi = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dbetabinom(y, aterms$vint1,
                           dpars$mu * dpars$phi,
                           (1 - dpars$mu) * dpars$phi, log = TRUE)
    },
    type = "discrete",
    # brms gets these three from the exposed Stan functions further
    # down the vignette. In frmtmb they are family metadata, declared
    # once here, so the three post-processing sections below have
    # nothing left to define.
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu * aterms$vint1,
      var_fn = function(dpars, aterms) {
        tt <- aterms$vint1
        tt * dpars$mu * (1 - dpars$mu) * (tt + dpars$phi) / (1 + dpars$phi)
      }
    ),
    sim = function(dpars, aterms, n) {
      p <- rbeta(n, dpars$mu * dpars$phi, (1 - dpars$mu) * dpars$phi)
      rbinom(n, aterms$vint1, p)
    }
  )
}, "SPELLING",
"lb/ub and vars have no frmtmb counterpart, type = 'int' becomes 'discrete', links is a named list, and an R lpdf replaces the Stan block")

# frmtmb has no brms counterpart for this: it checks the hand-written
# density against a numerical gradient before any fit is attempted.
bv("post", "ML: check_custom_family(beta_binomial2)", {
  check_custom_family(beta_binomial2, y = cbpp$incidence,
                      dpars = list(mu = rep(0.2, nrow(cbpp)),
                                   phi = rep(5, nrow(cbpp))),
                      aterms = list(vint1 = cbpp$size))
}, NA_character_, "")

# brms: stan_funs <- "real beta_binomial2_lpmf(...) {...}". There is no
# Stan program, so there is nothing for this string to be.
bv("post", "ML: stan_funs <- '<Stan source>'", {
  stop("no Stan source is involved: the density is the R `lpdf` above")
}, "MISSING", "the Stan functions block has no frmtmb analog; the lpdf argument of custom_family() carries the same content in R")

# brms: stanvars <- stanvar(scode = stan_funs, block = "functions").
bv("post", "ML: stanvars <- stanvar(scode = stan_funs, block = 'functions')", {
  stanvar(scode = "", block = "functions")
}, "MISSING", "stanvar() does not exist in frmtmb and nothing is named as a replacement; the lpdf argument covers the functions block, and there is no route for the other Stan blocks")

# brms: brm(..., family = beta_binomial2, stanvars = stanvars).
# `stanvars` is dropped with the Stan source. The family attaches with
# `+` because a custom family object is not a name frm(family =) can
# look up; the `+` form is the documented spelling for it.
# vint(size) carries the trials, exactly as in brms.
fit2 <- bv("model", "ML: fit2 (the custom family)", {
  frm(bf(incidence | vint(size) ~ period + (1 | herd)) + beta_binomial2,
      data = cbpp)
}, "SPELLING", "stanvars = is dropped, and the family attaches with bf(...) + fam rather than family =")

# brms reports Intercept -1.35 with sd(herd) = 0.38. The fixed effects
# match. sd(herd) collapses to the boundary here: once phi absorbs the
# overdispersion, maximum likelihood puts nothing left on the herd
# variance, and brms's prior is what keeps its posterior at 0.38.
bv("post", "ML: summary(fit2)", summary(fit2), "BEHAVIOR",
   "Intercept matches brms's -1.35, but sd(herd) goes to the zero boundary where brms reports 0.38: phi absorbs the overdispersion and ML has no prior holding the variance up")

# The vignette's own aside is that beta-binomial is native in brms now.
# It is native in frmtmb too, and it agrees with the hand-written
# family to the last digit, which is the check that the translation is
# faithful.
bv("post", "ML: the native beta_binomial() agrees with the hand-written one", {
  native <- frm(bf(incidence | trials(size) ~ period + (1 | herd)),
                data = cbpp, family = beta_binomial())
  stopifnot(all.equal(as.numeric(logLik(native)),
                      as.numeric(logLik(fit2)), tolerance = 1e-6))
  "identical"
}, NA_character_, "")

## ---- Post-processing custom family models ---------------------------

# brms: expose_functions(fit2, vectorize = TRUE), which compiles the
# Stan functions so the three R methods below can call them.
# The informative refusal is registered on frmtmb_draws only, so on the
# ML path the user meets a bare dispatch failure. The draws path further
# down gets the good message.
bv("post", "ML: expose_functions(fit2, vectorize = TRUE)", {
  expose_functions(fit2)
}, "MISSING", "on a frmtmb_fit this is a bare 'no applicable method' with no guidance; the informative refusal exists only on the draws object")

# brms: log_lik_beta_binomial2 <- function(i, prep) {...}, a method
# looked up by name. frmtmb has no such lookup: the lpdf given to
# custom_family() IS the log-likelihood, so this step is already done.
bv("post", "ML: log_lik_beta_binomial2 <- function(i, prep)", {
  ll <- beta_binomial2$lpdf(cbpp$incidence,
                            list(mu = rep(0.2, nrow(cbpp)),
                                 phi = rep(5, nrow(cbpp))),
                            list(vint1 = cbpp$size))
  print(round(sum(ll), 4))
}, "SPELLING", "no naming convention and no prep object: the lpdf passed to custom_family() is the log-likelihood already")

# loo() needs draws, so on a frmtmb_fit there are none. Section 2 runs
# the vignette's comparison properly.
bv("post", "ML: loo(fit1, fit2)", loo(fit1, fit2), "MISSING",
   "loo() is a frmtmb_draws method; on the ML path the comparison is AIC()/BIC(), and the sampling path below recovers the real thing")

# brms reports elpd_diff = -4.7 (se 4.1) in favor of fit2. AIC agrees on
# the direction and the rough size, with no elpd and no Pareto k.
bv("post", "ML: AIC(fit1, fit2) [the substitute]", {
  print(AIC(fit1, fit2))
}, "SPELLING", "the loo comparison becomes AIC/BIC; the direction matches brms's elpd_diff of -4.7 in favor of fit2")

# brms: posterior_predict_beta_binomial2 <- function(i, prep, ...).
# Already supplied, as the `sim` argument of custom_family().
bv("post", "ML: posterior_predict_beta_binomial2 <- function(i, prep, ...)", {
  print(utils::head(simulate(fit2, nsim = 1)[[1]]))
}, "SPELLING", "the random-draw function is the `sim` argument of custom_family(), given up front rather than defined later by name")

bv("post", "ML: pp_check(fit2)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(fit2))
}, NA_character_, "")

# The failure mode is worth pinning: a custom family without `sim`
# refuses pp_check() cleanly rather than producing wrong draws.
bv("post", "ML: pp_check() on a custom family with no sim", {
  bare <- custom_family(
    "beta_binomial2_bare", dpars = c("mu", "phi"),
    links = list(mu = "logit", phi = "log"),
    lpdf = beta_binomial2$lpdf, type = "discrete")
  fb <- frm(bf(incidence | vint(size) ~ period + (1 | herd)) + bare,
            data = cbpp)
  pp_check(fb)
}, "REFUSAL", "the refusal names the family and says it has no simulator; nothing names the `sim` argument that would fix it")

# brms: posterior_epred_beta_binomial2 <- function(prep). Already
# supplied, as post$mean_fn.
bv("post", "ML: posterior_epred_beta_binomial2 <- function(prep)", {
  print(round(utils::head(fitted(fit2)), 3))
}, "SPELLING", "the mean function is post$mean_fn inside custom_family(), not a separately named method")

# The vignette sets size = 1 so the y axis reads as a probability. That
# works here because post$mean_fn multiplies mu by vint1.
bv("post", "ML: conditional_effects(fit2, conditions = data.frame(size = 1))", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(fit2, conditions = data.frame(size = 1))
  print(ce[[1]][, c("period", "estimate__", "lower__", "upper__")])
  plot(ce)
  "plotted"
}, NA_character_, "")

## ============ PATH 2: sampling (frm_sample) ============
#
# The FIT route is used for both models: the vignette sets no priors, so
# there is nothing for the formula route's brms-like defaults to carry,
# and reusing the fitted objects keeps the two paths on one tape each.
# 400 iterations with 200 warmup is short. It is enough to show which
# entry points exist, and not enough for the elpd numbers to be stable,
# which is recorded rather than papered over.

s1 <- bv("model", "SAMPLE: fit1", {
  frm_sample(fit1, chains = 1, iter = 400, warmup = 200,
             seed = 1, cores = 1, refresh = 0)
}, NA_character_, "")

# The custom family samples with no extra declaration: the same R lpdf
# is the Stan model here. This is the single largest thing the sampling
# path adds to this vignette.
s2 <- bv("model", "SAMPLE: fit2 (the custom family)", {
  frm_sample(fit2, chains = 1, iter = 400, warmup = 200,
             seed = 1, cores = 1, refresh = 0)
}, NA_character_, "")

bv("post", "SAMPLE: summary(fit2 draws)", summary(s2), "BEHAVIOR",
   "the herd variance that ML pinned at zero now has a posterior, but 200 post-warmup draws leave its Rhat far above 1.01, so the number is not comparable to brms's 0.38")

bv("post", "SAMPLE: log_lik(fit2 draws)", {
  print(dim(log_lik(s2)))
}, NA_character_, "loo() and everything built on it work for a custom family with no extra user input")

bv("post", "SAMPLE: loo(fit1, fit2) [the vignette's comparison]", {
  print(loo_compare(loo(s1), loo(s2)))
}, "SPELLING",
"loo() takes one draws object, so the two-model form becomes loo_compare(loo(s1), loo(s2)); the short chains put the elpd_diff nowhere near brms's -4.7")

bv("post", "SAMPLE: posterior_predict(fit2 draws)", {
  print(dim(posterior_predict(s2)))
}, NA_character_, "")

bv("post", "SAMPLE: posterior_epred(fit2 draws)", {
  print(dim(posterior_epred(s2)))
}, NA_character_, "")

bv("post", "SAMPLE: pp_check(fit2 draws)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(s2))
}, NA_character_, "")

bv("post", "SAMPLE: expose_functions(fit2 draws)", expose_functions(s2),
   "REFUSAL", "same refusal on the draws object as on the fit, and the same right answer")

bv("post", "SAMPLE: conditional_effects(s2, conditions = data.frame(size = 1))", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s2, conditions = data.frame(size = 1)))
  "plotted"
}, NA_character_, "")

bv("post", "SAMPLE: stancode(s2)", stancode(s2), "MISSING",
   "there is no Stan program to print; the message points at the RTMB closure and the assembled frame instead")

bv_done()
