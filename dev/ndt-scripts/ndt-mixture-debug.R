# ndt lane: the mixture arm of test-defects.R, on either build.
#   Rscript --vanilla dev/ndt-scripts/ndt-mixture-debug.R ref|new
# Seed 3, the test's own.
a <- commandArgs(trailingOnly = TRUE)
arm <- if (length(a)) a[1L] else "new"
.libPaths(if (identical(arm, "ref")) {
  c("C:/Users/adf44/source/r/reflib-r2",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
} else {
  c("C:/Users/adf44/source/r/ndt-lib",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
})
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("arm:", arm, "\n")

set.seed(3)
dat <- ddm_simulate(300, mu = 0.9, bs = 1.4, ndt = 0.30)
k <- sample(300, 18)
dat$rt[k] <- stats::runif(18, 0.12, 2.5)
cat("min(rt):", format(min(dat$rt), digits = 8), "\n")

fam <- wiener(max_ndt = 0.4, allow_unreachable = TRUE)
cat("component ndt link name:", fam$links$ndt$name, "\n")
cat("component init ndt     :", fam$init_dpars$ndt(dat$rt, list()), "\n")
mx <- mixture(fam, lognormal())
cat("mixture ndt1 link name :", mx$links$ndt1$name, "\n")
cat("mixture ndt1 init      :",
    tryCatch(mx$init_dpars$ndt1(dat$rt, list()), error = function(e) NA),
    "\n")

fit <- frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5), family = mx, data = dat)
e <- unlist(fixef(fit))
cat("logLik:", format(as.numeric(logLik(fit)), digits = 10),
    " conv:", fit$opt$convergence, "\n")
print(round(e, 6))
cat("ndt (0.4 * plogis(eta)):",
    format(0.4 / (1 + exp(-e[["ndt1.(Intercept)"]])), digits = 8), "\n")
cat("finite residuals:",
    all(is.finite(stats::residuals(fit, type = "response"))), "\n")

# the log density the component gives at a few ndt values, holding the
# rest at the fit's own estimates: is the optimum at the wall or is the
# surface different from the reference build's?
lp <- function(t0) {
  dp <- list(mu = e[["mu1.(Intercept)"]], bs = exp(e[["bs1.(Intercept)"]]),
             ndt = t0, bias = 0.5)
  f <- if (is.null(fam[["ndt_bound"]])) fam else fam
  sum(as.numeric(f[["lpdf"]](dat$rt,
      lapply(if (is.null(fam[["ndt_bound"]])) dp else
             within(dp, ndt <- t0 / 0.4), function(v) rep(v, nrow(dat))),
      list(dec = as.numeric(dat$upper)))))
}
for (t0 in c(0.05, 0.1, 0.15, 0.18, 0.25, 0.3, 0.35)) {
  cat(sprintf("  component sum lpdf at ndt=%.2f : %12.4f\n", t0, lp(t0)))
}
