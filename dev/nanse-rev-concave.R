# Reviewer: are the lane's "concave" verdicts real saddles? For each
# fixture: the lost set, the most negative direction of the cached
# unit-diagonal Hessian, the objective along it (a saddle lowers the
# objective on at least one side), and a refit started a step along it
# on each side (a saddle lets the refit reach a higher logLik).
#   Rscript dev/nanse-rev-concave.R lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
probe <- function(lab, mkfit, refit = NULL) {
  cat("\n==", lab, "\n")
  w <- character()
  fit <- withCallingHandlers(mkfit(), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  for (x in w) cat("  warn:", substr(x, 1, 200), "\n")
  sdr <- ns$sdr_of(fit)
  cat("  lost:", paste(sprintf("%s(%s)", names(sdr$se_lost), sdr$se_lost),
                       collapse = " "), "\n")
  p <- fit$opt$par
  nm <- ns$outer_par_names(fit)
  H <- fit$cache$hessian_fixed$H %||% stats::optimHess(p, fit$obj$fn,
                                                       fit$obj$gr)
  H <- (H + t(H)) / 2
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), symmetric = TRUE)
  k <- which.min(e$values)
  cat(sprintf("  min unit-diag eigenvalue %.3g (largest %.3g)\n",
              e$values[k], e$values[1]))
  v <- e$vectors[, k] / D
  v <- v / sqrt(sum(v^2))
  cat("  direction:", paste(sprintf("%s=%.2f", nm, v)[abs(v) > 0.1],
                            collapse = " "), "\n")
  f0 <- fit$obj$fn(p)
  ch <- vapply(c(-1, -0.1, -0.01, 0.01, 0.1, 1), function(t) {
    fit$obj$fn(p + t * v) - f0
  }, 0)
  cat("  objective change at t = -1, -.1, -.01, .01, .1, 1:",
      paste(signif(ch, 3), collapse = " "), "\n")
  for (t in c(-0.3, 0.3)) {
    o <- tryCatch(stats::nlminb(p + t * v, fit$obj$fn, fit$obj$gr),
                  error = function(e) NULL)
    cat(sprintf("  nlminb on the same objective from t = %4.1f: ", t),
        if (is.null(o)) "error" else
          sprintf("logLik change %.3g, code %d", f0 - o$objective,
                  o$convergence), "
")
  }
}
# frmtmb.sample test-gr-by-draws.R's fixture
set.seed(11)
dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
f1 <- bf(y ~ x + (1 + x | gr(g, by = f)))
probe("gr(g, by = f)",
      function() frm(f1, family = gaussian(), data = dd),
      function(st) {
        fit0 <- frm(f1, family = gaussian(), data = dd,
                    control = frmtmb_control(check_se = "ignore"))
        tpl <- fit0$estimates
        # outer parameters in template order: replace and refit
        frm(f1, family = gaussian(), data = dd,
            start = ns$relist_outer(fit0, st),
            control = frmtmb_control(check_se = "ignore"))
      })
# frmtmb.spline test-difference.R's by-factor smooth beside exact gp()
set.seed(8)
n <- 120
d2 <- data.frame(x = sort(stats::runif(n, 0, 6)),
                 fac = factor(rep(c("A", "B"), length.out = n)))
d2$y <- sin(d2$x) + ifelse(d2$fac == "B", 0.3 * d2$x, 0) +
  stats::rnorm(n, 0, 0.3)
f2 <- bf(y ~ fac + s(x, by = fac, k = 6) + gp(x))
probe("s(x, by = fac) + gp(x)",
      function() frm(f2, family = gaussian(), data = d2),
      function(st) NULL)
