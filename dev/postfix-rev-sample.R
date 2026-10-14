# Reviewer of lane postfix: frmtmb.sample's new surfaces.
#  A. conditional_effects(method = "predict") on draws against brms's
#     (dev/postfix-rev-brms.R writes dev/postfix-rev-log/brms-ce.rds).
#  B. loo(d1, d2)'s ic_diffs__ and diffs against brms's own compare_ic()
#     and loo_compare() run on the same loo objects.
#  C. frm_sample() of a frm_multiple() result: R-hat, summary(), loo's
#     r_eff (brms sets r_eff = 1 for brm_multiple draws).
# Usage: Rscript dev/postfix-rev-sample.R
.libPaths(c("C:/Users/adf44/source/r/wt-postfix-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
source("dev/postfix-rev-data.R")
q <- function(e) suppressWarnings(suppressMessages(e))
ic <- list(x = c(-1, -0.5, 0, 0.5, 1))
cols <- c("x", "estimate__", "se__", "lower__", "upper__")

## A
sp <- q(frm_sample(frm(bf(y ~ x), family = poisson(), data = dp),
                   chains = 4, iter = 2000, seed = 1, refresh = 0))
sb <- q(frm_sample(frm(bf(y ~ x), family = bernoulli(), data = db),
                   chains = 4, iter = 2000, seed = 1, refresh = 0))
fr <- list(
  pois = conditional_effects(sp, "x", method = "predict",
                             int_conditions = ic, seed = 2)$x,
  pois_pp = conditional_effects(sp, "x", method = "posterior_predict",
                                int_conditions = ic, seed = 2)$x,
  pois_sqrt = conditional_effects(sp, "x", method = "predict",
                                  int_conditions = ic, transform = sqrt,
                                  seed = 2)$x,
  pois_epred = conditional_effects(sp, "x", int_conditions = ic)$x,
  bern = conditional_effects(sb, "x", method = "predict",
                             int_conditions = ic, seed = 2)$x)
br <- readRDS("dev/postfix-rev-log/brms-ce.rds")
cat("fixef poisson frmtmb.sample vs brms\n")
print(cbind(fixef(sp)[, 1:2], br$fixef_pois[, 1:2]))
for (k in c("pois", "pois_sqrt", "pois_epred", "bern")) {
  cat("--", k, ": frmtmb.sample then brms\n")
  print(fr[[k]][, cols]); print(br[[k]][, cols])
}
cat("predict and posterior_predict identical:",
    identical(fr$pois[, cols], fr$pois_pp[, cols]), "\n")
## the convention, read off posterior_predict() of the same draws
pp <- posterior_predict(sp, newdata = data.frame(x = ic$x))
cat("pp median", apply(pp, 2, median), "\n")
cat("pp mad   ", apply(pp, 2, mad), "\n")
cat("pp q2.5  ", apply(pp, 2, quantile, 0.025), "\n")
cat("pp q97.5 ", apply(pp, 2, quantile, 0.975), "\n")
cat("ce robust = FALSE (mean, sd):\n")
print(conditional_effects(sp, "x", method = "predict", int_conditions = ic,
                          robust = FALSE, seed = 2)$x[, cols])
cat("seed reproducible:", identical(
  conditional_effects(sp, "x", method = "predict", int_conditions = ic,
                      seed = 3)$x,
  conditional_effects(sp, "x", method = "predict", int_conditions = ic,
                      seed = 3)$x), "\n")
## a transform that is not monotone: brms applies it to the draws, so
## the median of |y - 3| is not |median - 3|
tr <- function(y) abs(y - 3)
ct <- conditional_effects(sp, "x", method = "predict", int_conditions = ic,
                          transform = tr, seed = 2)$x
cat("transform abs(y - 3): estimate__", ct$estimate__, "\n")
cat("  median of tr(pp)   ", apply(tr(pp), 2, median), "\n")

## B
s0 <- q(frm_sample(frm(bf(y ~ 1), family = poisson(), data = dp),
                   chains = 4, iter = 2000, seed = 1, refresh = 0))
two <- q(loo(sp, s0))
bc <- q(brms:::compare_ic(x = two$loos))
cat("ic_diffs__ lane:\n"); print(two$ic_diffs__)
cat("ic_diffs__ brms compare_ic on the same loos:\n"); print(bc$ic_diffs__)
cat("max abs diff", max(abs(two$ic_diffs__ - bc$ic_diffs__)), "\n")
cat("diffs vs loo::loo_compare:",
    isTRUE(all.equal(two$diffs, loo::loo_compare(two$loos))), "\n")
cat("names:", names(two), " loos names:", names(two$loos), "\n")
## models on different responses, same n: brms warns (yhash)
dp2 <- dp; dp2$y <- rev(dp2$y)
s2 <- q(frm_sample(frm(bf(y ~ x), family = poisson(), data = dp2),
                   chains = 2, iter = 1000, seed = 1, refresh = 0))
w <- character()
invisible(withCallingHandlers(loo(sp, s2), warning = function(x) {
  w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
}))
cat("loo(sp, s2) on different responses, warnings:\n")
cat(paste(" ", unique(w)), sep = "\n")
cat("\n")

## C
set.seed(5)
n <- 60
x <- rnorm(n); y <- 1 + 0.8 * x + rnorm(n, 0, 0.5)
x[1:25] <- NA
imps <- lapply(1:3, function(j) {
  xi <- x
  # deliberately different imputations, so the posteriors differ
  xi[is.na(xi)] <- rnorm(25, c(-2, 0, 2)[j], 0.3)
  data.frame(x = xi, y = y)
})
fm <- q(frm_multiple(bf(y ~ x), family = gaussian(), data = imps))
sm <- q(frm_sample(fm, chains = 2, iter = 1000, seed = 7, refresh = 0))
cat("nchains", rstan::nchains(sm$stanfit), " draws", nrow(sm$draws), "\n")
w <- character()
sumtab <- withCallingHandlers(summary(sm), warning = function(x) {
  w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
})
print(sumtab)
cat("summary warnings:", length(w), paste(unique(w), collapse = " | "), "\n")
cat("rhat():\n"); print(rhat(sm))
w <- character()
pr <- capture.output(withCallingHandlers(print(sm), warning = function(x) {
  w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
}))
cat(head(pr, 40), sep = "\n")
cat("print warnings:", length(w), paste(unique(w), collapse = " | "), "\n")
## loo r_eff across imputations vs brms's 1
ll <- log_lik(sm)
re <- loo::relative_eff(exp(ll), chain_id = attr(ll, "chain_id"))
cat("r_eff from chains across imputations: range", range(re), "\n")
l_lane <- q(loo(sm))
l_one <- q(loo::loo(ll, r_eff = rep(1, ncol(ll))))
print(l_lane$estimates); print(l_one$estimates)
cat("in-sample rows of the combined draws read imputation 1:",
    isTRUE(all.equal(sm$fit$data$x, imps[[1]]$x)), "\n")
