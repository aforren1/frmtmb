# Reviewer: frmtmb.spline test-difference.R's by-factor smooth beside an
# exact gp(), which the trial merge warns on ("curves downward" on 6 of
# 9 parameters). Is the fit a saddle, or is the finite-difference
# Hessian of the Laplace objective wrong? Compares sdreport() on the
# release build, optimHess at two steps, and a refit from perturbed
# starts.
#   Rscript dev/nanse-rev-spline.R release|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  release = "C:/Users/adf44/source/r/rellib-r6",
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
set.seed(8)
n <- 120
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
w <- character()
fit <- withCallingHandlers(
  frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), family = gaussian(),
      data = d),
  warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
nm <- frmtmb:::outer_par_names(fit)
p <- fit$opt$par
cat("code", fit$opt$convergence, "logLik", format(logLik(fit), digits = 12),
    "\n")
for (x in w) cat("warn:", substr(x, 1, 250), "\n")
cat("par:", paste(sprintf("%s=%.4g", nm, p), collapse = " "), "\n")
cat("grad:", paste(sprintf("%.2g", fit$obj$gr(p)), collapse = " "), "\n")
sdr <- RTMB::sdreport(fit$obj)
cat("plain sdreport: pdHess", sdr$pdHess, " SE",
    paste(signif(sqrt(pmax(diag(sdr$cov.fixed), NA)), 3), collapse = " "),
    "\n")
for (h in c(1e-3, 1e-4)) {
  H <- stats::optimHess(p, fit$obj$fn, fit$obj$gr,
                        control = list(ndeps = rep(h, length(p))))
  H <- (H + t(H)) / 2
  D <- sqrt(abs(diag(H)))
  ev <- eigen(H / outer(D, D), symmetric = TRUE)
  cat(sprintf("optimHess ndeps %g: unit-diag eigenvalues %s\n", h,
              paste(signif(ev$values, 3), collapse = " ")))
  k <- which.min(ev$values)
  cat("   min direction:", paste(sprintf("%s=%.2f", nm,
                                          ev$vectors[, k])[
                                            abs(ev$vectors[, k]) > 0.1],
                                  collapse = " "), "\n")
}
# step along the most negative direction: does the likelihood rise?
H <- stats::optimHess(p, fit$obj$fn, fit$obj$gr)
H <- (H + t(H)) / 2
e <- eigen(H, symmetric = TRUE)
v <- e$vectors[, which.min(e$values)]
f0 <- fit$obj$fn(p)
for (t in c(-1, -0.3, -0.1, 0.1, 0.3, 1)) {
  cat(sprintf("  t = %5.2f: objective change %.3g\n", t,
              fit$obj$fn(p + t * v) - f0))
}
