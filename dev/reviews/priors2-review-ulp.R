# Review of lane priors2, item 5: does bound_inside() land inside on
# every end, or can its 64-step cap leave a bound outside? Where the
# inverse link is flat in ulps of the internal value (log near 1, logit
# near 1, probit and cloglog tails), one ulp of v takes many ulps of b.
# Rscript priors2-review-ulp.R
.libPaths(c("C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
bi <- frmtmb:::internal_bound
miss <- function(lab, link, ends, sd = FALSE) {
  pm <- if (sd) list(scale = "sd") else
    list(scale = "natural", link = frmtmb:::get_link(link))
  inv <- if (sd) exp else pm$link$linkinv
  lo <- vapply(ends, bi, 0, pm = pm, which = "lb")
  hi <- vapply(ends, bi, 0, pm = pm, which = "ub")
  bad_lo <- ends[inv(lo) < ends]
  bad_hi <- ends[inv(hi) > ends]
  cat(sprintf("%-28s %4d ends: lb outside %3d, ub outside %3d", lab,
              length(ends), length(bad_lo), length(bad_hi)), "\n")
  if (length(bad_lo)) cat("   lb e.g.", format(head(bad_lo, 5), digits = 17),
                          "\n")
  if (length(bad_hi)) cat("   ub e.g.", format(head(bad_hi, 5), digits = 17),
                          "\n")
  t0 <- proc.time()[["elapsed"]]
  for (i in 1:20) vapply(ends, bi, 0, pm = pm, which = "ub")
  cat(sprintf("   %.2f ms per bound", 1000 * (proc.time()[["elapsed"]] - t0) /
                (20 * length(ends))), "\n")
}
miss("log, 0.90 to 1.10", "log", seq(0.9, 1.1, by = 0.001))
miss("log, 0.05 to 5", "log", seq(0.05, 5, by = 0.05))
miss("sd scale, 0.90 to 1.10", NULL, seq(0.9, 1.1, by = 0.001), sd = TRUE)
miss("logit, 0.01 to 0.99", "logit", seq(0.01, 0.99, by = 0.01))
miss("logit, 0.99 to 0.99999", "logit", 1 - 10^-seq(2, 5, by = 0.01))
miss("probit, 0.01 to 0.99", "probit", seq(0.01, 0.99, by = 0.01))
miss("cloglog, 0.01 to 0.99", "cloglog", seq(0.01, 0.99, by = 0.01))
# a fit at such a bound: sigma's ub just above 1 where the data's sigma
# is above it, 5027's NA/NaN shape
set.seed(2)
ds <- data.frame(z = rnorm(120))
ds$y <- 1 + 0.3 * ds$z + rnorm(120, sd = 1.4)
pm <- list(scale = "natural", link = frmtmb:::get_link("log"))
cand <- seq(1.001, 1.1, by = 0.001)
hi <- vapply(cand, bi, 0, pm = pm, which = "ub")
badv <- cand[exp(hi) > cand]
cat("ub values on sigma that stay outside:", length(badv), "\n")
for (ub in head(badv, 2)) {
  pr <- set_prior(sprintf("uniform(0.1, %.17g)", ub), class = "sigma",
                  lb = 0.1, ub = ub)
  r <- tryCatch(withCallingHandlers(frm(y ~ z, data = ds, prior = pr),
                                    warning = function(w) {
                                      cat("  [warning]",
                                          substr(conditionMessage(w), 1, 120),
                                          "\n")
                                      invokeRestart("muffleWarning")
                                    }),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat("ub =", format(ub, digits = 17), ": ERROR",
        substr(conditionMessage(r), 1, 160), "\n")
  } else {
    cat("ub =", format(ub, digits = 17), ": code", r$opt$convergence,
        "objective", format(r$opt$objective, digits = 12), "sigma",
        format(exp(r$estimates$betad[1]), digits = 17), "\n")
  }
}
