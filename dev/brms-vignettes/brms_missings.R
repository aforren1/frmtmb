# Translation of brms's "Handle Missing Values with brms"
# (brms 2.23.0, doc/brms_missings.Rmd) onto the frmtmb surface.
#
# Read the comments, not just the code: every `<EDGE>` marker on a bv()
# call is a place where the brms line did not carry over verbatim.
#
# Two inference paths are covered. Section 1 is the maximum-likelihood
# path, which is the one frmtmb is built for. Section 2 samples the same
# models with frm_sample(), because a posterior is what the brms
# vignette actually shows and the draws surface answers a different set
# of the vignette's questions.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_missings")

options(mc.cores = 1)
set.seed(1234)

## ---- data ----------------------------------------------------------

bv("data", "data(nhanes, package = 'mice'); head(nhanes)", {
  data("nhanes", package = "mice", envir = globalenv())
  print(utils::head(nhanes))
}, NA_character_, "")

## ============ PATH 1: ML / Laplace (frm) ============

## ---- Imputation before model fitting --------------------------------

# The vignette uses m = 5. Three imputations keep the whole script
# inside the time budget and change nothing about which frmtmb surface
# is exercised.
imp <- bv("data", "mice(nhanes, m = 3, print = FALSE)", {
  mice::mice(nhanes, m = 3, print = FALSE)
}, NA_character_, "")

# brms: brm_multiple(bmi ~ age*chl, data = imp, chains = 2). `chains` is
# MCMC-only. A mids object is accepted directly, as in brms, and the
# gaussian family is defaulted, so nothing else changes.
# (The v0.34 audit recorded the missing family default as FN-1; it is
# fixed. frm() and frm_multiple() both default to gaussian() now.)
fit_imp1 <- bv("model", "ML: fit_imp1", {
  frm_multiple(bmi ~ age * chl, data = imp)
}, NA_character_, "")

# brms prints one pooled summary. frmtmb has no summary method for a
# frm_multiple result, so summary() falls through to summary.default and
# reports the object's storage layout instead of the pooled table.
bv("post", "ML: summary(fit_imp1)", summary(fit_imp1), "BEHAVIOR",
   "no summary method for frmtmb_multiple: summary.default reports the list layout, not the pooled coefficient table")

# print() is where the pooled table lives. It carries Rubin's df, p and
# fmi, which brms's Rhat/ESS block has no counterpart for, and it has no
# Rhat, so the vignette's whole "Rhat is a false positive here"
# discussion has nothing to attach to.
bv("post", "ML: print(fit_imp1) [the pooled table]", {
  print(fit_imp1)
  "printed"
}, "SPELLING", "the pooled table is on print(), not summary(); it reports df/p/fmi where brms reports Rhat/ESS")

# brms: plot(fit_imp1, variable = "^b", regex = TRUE) draws the per-chain
# traces the vignette's false-positive-Rhat argument rests on. There are
# no chains here, so this is refused.
# (The v0.34 audit recorded FN-11's plot() as a raw base-graphics error,
# "'x' is a list, but does not have components 'x' and 'y'". It is now a
# refusal that names the workaround.)
bv("post", "ML: plot(fit_imp1, variable = '^b', regex = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fit_imp1, variable = "^b", regex = TRUE)
}, "REFUSAL", "no pooled display; the message names plot(x$fits[[1]]) and it is the right replacement")

# The vignette's convergence block:
#   draws <- as_draws_array(fit_imp1)
#   nc <- nchains(fit_imp1) / m
#   lapply(split, summarise_draws, default_convergence_measures())
# All four entry points need draws.
bv("post", "ML: as_draws_array(fit_imp1)", posterior::as_draws_array(fit_imp1),
   "REFUSAL", "the message names x$pooled, hypothesis() and frm_sample(x$fits[[1]]) and all three are right")

# nchains() is the one inconsistency in that refusal set. frmtmb exports
# its own nchains generic and registers the frmtmb_multiple method on
# posterior's, so an unqualified call finds no method while the
# qualified one gives the good refusal.
bv("post", "ML: nchains(fit_imp1)", nchains(fit_imp1), "BEHAVIOR",
   "frmtmb's own nchains() generic has no frmtmb_multiple method; posterior::nchains() has the informative refusal, so the answer depends on which nchains you call")

bv("post", "ML: posterior::nchains(fit_imp1)", posterior::nchains(fit_imp1),
   "REFUSAL", "the qualified call gives the refusal with the workaround")

# The convergence question itself still has a frequentist answer: the
# per-imputation fits are stored, so the spread of their estimates is
# readable directly.
bv("post", "ML: per-submodel check (the substitute for summarise_draws)", {
  est <- vapply(fit_imp1$fits, function(f) unlist(fixef(f)), numeric(5))
  print(round(t(est), 3))
}, "SPELLING", "default_convergence_measures() becomes reading the m stored fits; there is no Rhat to compute")

# brms: conditional_effects(fit_imp1, "age:chl").
# (The v0.34 audit recorded FN-11's conditional_effects as "no
# applicable method"; a frmtmb_multiple method is now registered and
# refuses on purpose.)
bv("post", "ML: conditional_effects(fit_imp1, 'age:chl')", {
  conditional_effects(fit_imp1, "age:chl")
}, "REFUSAL", "pooling an effect curve across imputations is not implemented; the message names conditional_effects(x$fits[[1]])")

bv("post", "ML: conditional_effects(fit_imp1$fits[[1]], 'age:chl') [workaround]", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(fit_imp1$fits[[1]], "age:chl")
  plot(ce)
  "plotted"
}, "SPELLING", "the curve comes from one imputation, so it carries no between-imputation uncertainty")

# hypothesis() is the one inferential post-processing entry point that
# does read the pooled table, so a coefficient test survives pooling.
bv("post", "ML: hypothesis(fit_imp1, 'age')", {
  hypothesis(fit_imp1, "age")
}, NA_character_, "")

## ---- Compatibility with other multiple imputation packages ----------

# The vignette says brm_multiple also takes a plain list of data frames.
bv("data", "mice::complete(imp, action = 'all')", {
  imp_list <- mice::complete(imp, action = "all")
  assign("imp_list", imp_list, envir = globalenv())
  length(imp_list)
}, NA_character_, "")

bv("model", "ML: fit_imp1 from a list of data frames", {
  frm_multiple(bmi ~ age * chl, data = imp_list)
}, NA_character_, "")

## ---- Imputation during model fitting --------------------------------

# brms: brm(bform, data = nhanes). The bf() + bf() + set_rescor(FALSE)
# multivariate spelling carries over verbatim, and so do both mi()
# forms.
bform <- bf(bmi | mi() ~ age * mi(chl)) +
  bf(chl | mi() ~ age) + set_rescor(FALSE)

fit_imp2 <- bv("model", "ML: fit_imp2", {
  frm(bform, data = nhanes)
}, NA_character_, "")

# One divergence from brms's printed summary: the "Family:" line is
# empty on a multivariate fit, where brms names the per-response
# families. The coefficient blocks themselves are right, and the
# mi(chl) coefficient is spelled `michl` where brms spells it
# `bmi_michl`.
bv("post", "ML: summary(fit_imp2)", summary(fit_imp2), "BEHAVIOR",
   "the Family: line prints empty on a multivariate fit, and the mi() coefficient is named michl where brms names it bmi_michl")

# (The v0.34 audit recorded FN-13: conditional_effects on an mi() fit
# could not build its grid, "mi(chl): newdata must supply complete
# values". It is fixed.)
bv("post", "ML: conditional_effects(fit_imp2, 'age:chl', resp = 'bmi')", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  ce <- conditional_effects(fit_imp2, "age:chl", resp = "bmi")
  plot(ce)
  "plotted"
}, NA_character_, "")

## ---- Combining measurement error and missing values -----------------

bv("data", "nhanes$se <- rexp(nrow(nhanes), 2)", {
  nhanes <- get("nhanes", envir = globalenv())
  nhanes$se <- rexp(nrow(nhanes), 2)
  assign("nhanes", nhanes, envir = globalenv())
  summary(nhanes$se)
}, NA_character_, "")

# The vignette's fit_imp3 chunk is eval = FALSE, but it is the code a
# reader copies, so it is translated and run.
fit_imp3 <- bv("model", "ML: fit_imp3 (mi(se) measurement error)", {
  bform3 <- bf(bmi | mi() ~ age * mi(chl)) +
    bf(chl | mi(se) ~ age) + set_rescor(FALSE)
  frm(bform3, data = nhanes)
}, NA_character_, "")

bv("post", "ML: summary(fit_imp3)", summary(fit_imp3), "BEHAVIOR",
   "same empty Family: line as fit_imp2")

## ============ PATH 2: sampling (frm_sample) ============
#
# The vignette sets no priors, so every model here takes the FIT route:
# frm_sample(fit, ...) reuses the fitted object's own tape and starting
# values. 400 iterations with 200 warmup is far short of what the
# vignette's 2 chains x 2000 give, so the chains are noisy on purpose;
# what is being measured is which entry points exist, not agreement.

# frm_multiple() has no sampling analogue: Rubin's rules pool point
# estimates, and there is nothing in frmtmb that pools m posteriors the
# way brm_multiple() concatenates m chains.
bv("model", "SAMPLE: fit_imp1 (the pooled object)", {
  frm_sample(fit_imp1, chains = 1, iter = 400, warmup = 200,
             seed = 1, cores = 1, refresh = 0)
}, "MISSING",
"frm_sample() takes a fit or a formula, not a frmtmb_multiple; brm_multiple()'s pooled posterior has no frmtmb analog. The workaround is one imputation at a time")

s_imp1 <- bv("model", "SAMPLE: fit_imp1$fits[[1]] (one imputation)", {
  frm_sample(fit_imp1$fits[[1]], chains = 1, iter = 400, warmup = 200,
             seed = 1, cores = 1, refresh = 0)
}, "SPELLING", "one imputation's posterior stands in for the m-fold pooled posterior brm_multiple() returns")

bv("post", "SAMPLE: summary(fit_imp1 draws)", summary(s_imp1), NA_character_, "")

# This is the vignette's own point, recovered: nchains/ndraws exist on
# the draws object, so the convergence block has somewhere to run. It
# runs on one imputation, not on m of them, so the non-overlaying chains
# the vignette is about cannot appear.
bv("post", "SAMPLE: nchains/ndraws/as_draws_df", {
  print(c(nchains = nchains(s_imp1), ndraws = ndraws(s_imp1)))
  print(dim(posterior::as_draws_df(s_imp1)))
}, "BEHAVIOR",
"the entry points work, but with one imputation per draws object the vignette's between-imputation Rhat inflation cannot be reproduced")

bv("post", "SAMPLE: conditional_effects(s_imp1, 'age:chl')", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s_imp1, "age:chl"))
  "plotted"
}, NA_character_, "")

bv("post", "SAMPLE: pp_check(s_imp1)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(s_imp1))
}, NA_character_, "")

bv("post", "SAMPLE: loo(s_imp1)", loo(s_imp1), NA_character_,
   "loo() is a draws method, so the sampling path recovers what the ML path has to answer with AIC()")

# The one-sided form the ML path used to refuse.
bv("post", "SAMPLE: hypothesis(s_imp1, 'age > 0')", {
  hypothesis(s_imp1, "age > 0")
}, NA_character_, "")

s_imp2 <- bv("model", "SAMPLE: fit_imp2", {
  frm_sample(fit_imp2, chains = 1, iter = 400, warmup = 200,
             seed = 1, cores = 1, refresh = 0)
}, NA_character_, "")

bv("post", "SAMPLE: summary(fit_imp2 draws)", summary(s_imp2), NA_character_,
   "")

# The imputed values are parameters here, so a posterior over them
# exists. This is the closest thing in frmtmb to what the vignette calls
# imputation during model fitting.
bv("post", "SAMPLE: conditional_effects(s_imp2, 'age:chl', resp = 'bmi')", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(s_imp2, "age:chl", resp = "bmi"))
  "plotted"
}, NA_character_, "")

# conditional_effects() on draws has no `method =` and no `band =`; the
# band comes from the draws themselves.
bv("post", "SAMPLE: conditional_effects(s_imp2, method = 'predict')", {
  conditional_effects(s_imp2, "age:chl", resp = "bmi", method = "predict")
}, "REFUSAL", "the draws method has no method = argument; the message names posterior_predict()")

bv("post", "SAMPLE: pp_check(s_imp2)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(s_imp2))
}, "REFUSAL", "pp_check() refuses multivariate fits on the draws path as well; nothing is named as a replacement")

bv("post", "SAMPLE: loo(s_imp2)", loo(s_imp2), "REFUSAL",
   "log_lik() refuses an mi() model because a latent row's column is not the observation's own likelihood; the message names frm_multiple() plus AIC() and frm_bootstrap()")

bv("post", "SAMPLE: posterior_predict(s_imp2)", {
  pp <- posterior_predict(s_imp2)
  print(dim(pp))
}, NA_character_, "")

bv_done()
