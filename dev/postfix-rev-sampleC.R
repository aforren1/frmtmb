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
cat("nchains", sm$stanfit@sim$chains, " draws", nrow(sm$draws), "\n")
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
