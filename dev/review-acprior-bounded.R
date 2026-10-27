# Reviewer: frm_sample()'s bounded natural-scale defaults (sigma under
# an identity link, nu, ndt) under a user prior of the same slot. brms
# keeps the declared bound (dev/review-acprior-log/brms.txt).
# Usage: Rscript dev/review-acprior-bounded.R before|after
# Log: dev/review-acprior-log/bounded-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(9)
dd <- data.frame(x = rnorm(60))
dd$y <- 1 + dd$x + rnorm(60)
dd$yp <- 0.3 + exp(rnorm(60, 0, 0.5))
fmt <- function(v) {
  if (!length(v)) return("(none)")
  paste(names(v), signif(unlist(v), 6), collapse = "; ")
}
showb <- function(lab, fit, pr) {
  cat("\n---", lab, "\n")
  tryCatch({
    ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pr))$ri
    cat("  lower:", fmt(ri$lower), "\n  upper:", fmt(ri$upper), "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
fs <- frm(y ~ x, data = dd,
          family = brmsfamily("gaussian", link_sigma = "identity"),
          dry_run = "objective")
showb("sigma identity: default", fs, NULL)
showb("sigma identity: normal(1, 1)", fs,
      set_prior("normal(1, 1)", class = "sigma"))
ft <- frm(y ~ x, data = dd, family = brmsfamily("student", link_nu = "identity"),
          dry_run = "objective")
showb("nu identity: default", ft, NULL)
showb("nu identity: normal(5, 2)", ft, set_prior("normal(5, 2)", class = "nu"))
fn <- frm(yp ~ x, data = dd, family = brmsfamily("shifted_lognormal"),
          dry_run = "objective")
showb("ndt log link: default", fn, NULL)
showb("ndt log link: normal(0.2, 0.1)", fn,
      set_prior("normal(0.2, 0.1)", class = "ndt"))
# the cosy analogue, for the record
set.seed(5)
d <- data.frame(t = rep(1:5, 10), g = rep(1:10, each = 5), x = rnorm(50))
d$y <- d$x + rnorm(50)
fc <- frm(y ~ x + cosy(t, g), data = d, dry_run = "objective")
showb("cosy: default", fc, NULL)
showb("cosy: normal(0, 0.1)", fc, set_prior("normal(0, 0.1)", class = "cosy"))
# fix round 1: a user bound on one side is kept, the other side lent
showb("nu identity: normal(5, 2), lb = 1.5", ft,
      set_prior("normal(5, 2)", class = "nu", lb = 1.5))
showb("ndt log link: normal(0.2, 0.1), lb = 0.1", fn,
      set_prior("normal(0.2, 0.1)", class = "ndt", lb = 0.1))
showb("cosy: normal(0, 0.1), ub = 0.5", fc,
      set_prior("normal(0, 0.1)", class = "cosy", ub = 0.5))
showb("cosy: bounds only lb = -0.1", fc,
      set_prior("", class = "cosy", lb = -0.1))
showb("sigma identity: bounds only lb = 0.5", fs,
      set_prior("", class = "sigma", lb = 0.5))
