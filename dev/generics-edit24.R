source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")

# B1. The NEWS bullet.  The version NUMBER is deferred under the lane
# rules; the bullet is not a number and is this lane's to write.
news <- paste0(
"# frmtmb (development version)\n",
"\n",
"* **frmtmb no longer breaks brms, lme4, posterior, loo, rstantools or\n",
"  bayesplot for objects that already exist.** Attaching frmtmb after\n",
"  one of them used to cost that package every S3 method it had on a\n",
"  name the two share, because frmtmb defined and exported its own\n",
"  generic and `UseMethod()` reads the method table of the namespace\n",
"  where the generic it reached was defined. Measured on the previous\n",
"  release: `library(brms); library(frmtmb)` then `loo(fit)` on an\n",
"  existing `brmsfit` lost all 27 of the generics the two share, two of\n",
"  them SILENTLY by falling into frmtmb's own `.default`; it is now 0\n",
"  of 27, in every load order tested. `library(lme4); library(frmtmb)`\n",
"  lost all 5 of the names lme4 and frmtmb share; it is now 0 of 5.\n",
"\n",
"  frmtmb now shares the owner's generic instead of shadowing it.\n",
"  `fixef()`, `ranef()` and `VarCorr()` are imported from nlme, which\n",
"  is Recommended and is where lme4, glmmTMB and brms get them too, and\n",
"  `refit()` from generics. The rest resolve to the owner's generic at\n",
"  run time when the owner is loaded, and to frmtmb's own when it is\n",
"  not, so nothing new appears in Imports.\n",
"\n",
"* **BREAKING, one call.** `posterior_summary()`'s first argument is\n",
"  now `x`, not `object`, matching brms's generic, so\n",
"  `posterior_summary(object = m)` fails and `posterior_summary(x = m)`\n",
"  or `posterior_summary(m)` is the spelling. That is the only\n",
"  caller-visible break: `nvariables()` and `VarCorr()` change too, but\n",
"  both only GAIN an argument (`...` from posterior, `sigma` from\n",
"  nlme), so no call that named their first argument breaks.\n",
"\n",
"  `refit()` is a fourth changed signature and a subtler one: it is\n",
"  `(object, ...)` from generics normally and `(object, newresp, ...)`\n",
"  when lme4 is loaded, because lme4 defines its own `refit` generic\n",
"  rather than importing one, so `args(refit)` now depends on what else\n",
"  is loaded. No caller breaks on it.\n",
"\n",
"  **If you write an S3 method on any of these names, match the\n",
"  OWNER's formals, not frmtmb's.** `frmtmb.sample` is updated here for\n",
"  that reason and its `frmtmb (>= ...)` floor moves with this release.\n",
"\n",
"* `as_draws()`, `as_draws_array()`, `as_draws_df()`,\n",
"  `as_draws_list()`, `as_draws_matrix()`, `as_draws_rvars()`,\n",
"  `ndraws()`, `nchains()`, `niterations()` and `nvariables()` now\n",
"  refuse a `frmtmb_fit` by name and point at\n",
"  `frmtmb.sample::frm_sample()`. They had no method for the class at\n",
"  all, which was R's own \"no applicable method\" before and would have\n",
"  become posterior's \"All list elements must be lists themselves\"\n",
"  after, since a fit is a bare list.\n",
"\n",
"* `?frmtmb-scales` is new and says, per method, whether a number is on\n",
"  the link scale, the response scale or neither. Two places where\n",
"  frmtmb and brms disagree are written down there rather than left to\n",
"  be discovered: `predict()` returns the linear predictor where brms\n",
"  returns the response scale, and `summary()` prints `sigma` on its\n",
"  log link where `sigma()` back-transforms.\n",
"\n")
p <- file.path("C:/Users/adf44/source/r/frmtmb-wt-generics", "NEWS.md")
cur <- rawToChar(readBin(p, "raw", file.info(p)$size))
if (grepl("frmtmb (development version)", cur, fixed = TRUE)) {
  stop("a development section already exists; edit it instead")
}
writeBin(charToRaw(paste0(news, cur)), p)
cat("prepended the NEWS section\n")

# The frmtmb.sample floor must move with the bump, and the NUMBER is
# the user's, so nothing is changed here: it is named in NEWS and in
# dev/generics-findings.md instead.
cat("DONE\n")
