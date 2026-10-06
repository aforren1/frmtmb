# Copied unchanged from dev/ordmix-rev2-b1.R at the 0.68.0 release, except that it loads
# the release library rellib-r6 (dev/round-20261005.md).
# Reviewer of lane ordmix, re-check of B1: the edges of the fixed
# flat-direction warnings. Data: gen() of dev/ordmix-rev-falsealarm.R,
# seeds 1..R, n = 500. Per fit: whether a flat-direction warning fired,
# and the largest standard error (a flat direction shows as NaN or a
# huge value).
# Usage: Rscript dev/ordmix-rev2-b1.R <case> [R]
args <- commandArgs(TRUE)
case <- args[1]
R <- if (length(args) >= 2) as.integer(args[2]) else 10L
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ex <- parse("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel068-falsealarm.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("cut4", "gen")) eval(e)
}
hm <- function() mixture(hurdle_cumulative(), hurdle_cumulative())
spec <- switch(case,
  # K = 2 hurdle: one prior on hu1 alone identifies (3 parameters, 2
  # numbers seen, 1 flat direction)
  hu1_only = list(f = bf(yh ~ x), fam = hm(),
                  prior = set_prior("beta(4, 16)", class = "hu1")),
  # a prior on the thresholds alone leaves hu and theta flat
  thres_only = list(f = bf(yh ~ x), fam = hm(),
                    prior = set_prior("normal(0, 2)", class = "Intercept",
                                      dpar = "mu1")),
  # the prior the warning names, brms's beta(1, 1), on both
  beta11 = list(f = bf(yh ~ x), fam = hm(),
                prior = set_prior("beta(1, 1)", class = "hu1") +
                  set_prior("beta(1, 1)", class = "hu2")),
  # no prior: the warning must fire
  none = list(f = bf(yh ~ x), fam = hm()),
  # no predictor anywhere, a prior on component 1's thresholds only
  nopred_mu1 = list(f = bf(y ~ 1), fam = mixture(cumulative(), cumulative()),
                    prior = set_prior("normal(0, 2)", class = "Intercept",
                                      dpar = "mu1")),
  # shared thresholds, no mu predictor, identical except a held disc
  disc_differs = list(f = bf(ysh ~ 1, theta1 ~ z, disc2 = 2),
                      fam = mixture(cumulative(), cumulative(),
                                    order = "mu")),
  disc_same = list(f = bf(ysh ~ 1, theta1 ~ z, disc1 = 2, disc2 = 2),
                   fam = mixture(cumulative(), cumulative(), order = "mu"))
)
flat_pat <- paste("no distributional parameter has a predictor",
                  "the likelihood sees only P(Y = 0)",
                  "the components are one distribution", sep = "|")
for (s in seq_len(R)) {
  d <- gen(s)
  w <- character(0)
  fit <- tryCatch(withCallingHandlers(
    frm(spec$f, family = spec$fam, data = d, prior = spec$prior),
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
  if (is.character(fit)) {
    cat(sprintf("REP case=%s seed=%d ERROR %s\n", case, s,
                substr(fit, 1, 200)))
    next
  }
  se <- suppressWarnings(fixef(fit)[, "Est.Error"])
  hu_se <- tryCatch(suppressWarnings({
    v <- sqrt(diag(vcov(fit, full = TRUE)))
    max(v)
  }), error = function(e) NA_real_)
  cat(sprintf("REP case=%s seed=%d flatwarn=%s maxse_fixef=%.3g maxse_full=%.3g nwarn=%d\n",
              case, s, any(grepl(flat_pat, w)), max(se), hu_se, length(w)))
  for (m in unique(w)) cat("   WARN", substr(m, 1, 150), "\n")
}
