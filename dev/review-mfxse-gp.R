# Review of lane mfxse, attack 7: the average slope of an exact gp() field.
#   Rscript dev/review-mfxse-gp.R
# Data seed 8150 (dev/mfxse-gpslope.R's design, y = sin(x) + 0.3 z + noise).
# 1. The same model by four routes, all with z: is 0.0308 the exact field's
#    own value, or one route's?
# 2. The step: marginaleffects' eps and numderiv.
# 3. The split of the delta variance: inner (theta held) against carried.
# 4. A profile: the average slope of the conditional mode f(theta) as one
#    theta moves over +-3 standard errors, inner problem re-solved, against
#    its linearization. This is what the delta method assumes is linear.
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(lane, base, user))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressWarnings(suppressMessages(e))
set.seed(8150)
m <- 90
d <- data.frame(x = runif(m, 0, 6), z = runif(m, 0, 6))
d$y <- sin(d$x) + 0.3 * d$z + rnorm(m, sd = 0.3)
h <- 1e-4
avg_sl <- function(f) {
  mean((marginaleffects::get_predict(f, newdata = transform(d, x = x + h / 2))$estimate -
          marginaleffects::get_predict(f, newdata = transform(d, x = x - h / 2))$estimate) / h)
}
cat("\n1. routes (all y ~ z + ...):\n")
specs <- c("gp(x)", "gp(x, gr = FALSE)", "gp(x, k = 15)", "gp(x, k = 40)", "s(x, k = 20)")
fits <- list()
for (s in specs) {
  f <- quiet(frm(as.formula(paste("y ~ z +", s)), data = d))
  fits[[s]] <- f
  sl <- quiet(avg_slopes(f, variables = "x"))
  th <- f$estimates$theta
  cat(sprintf("  %-18s estimate %.5f se %.4f  theta %s  conv %d\n", s, sl$estimate,
              sl$std.error, paste(sprintf("%.3f", th), collapse = " "), f$opt$convergence))
}
fit <- fits[["gp(x)"]]
cat("\n2. step:\n")
for (e in c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)) {
  a <- quiet(avg_slopes(fit, variables = "x", eps = e))
  b <- quiet(avg_slopes(fit, variables = "x", eps = e, numderiv = "fdcenter"))
  cat(sprintf("  eps %.0e: se %.5f, fdcenter %.5f, estimate %.5f\n", e, a$std.error,
              b$std.error, a$estimate))
}
cat("\n3. split of the delta variance:\n")
jc <- ns$get_joint_cov(fit)
rn <- jc$names; V <- as.matrix(jc$V)
cat("  joint rows:", paste(names(table(rn)), table(rn), collapse = ", "), "\n")
mp <- list(); for (cp in unique(rn)) mp[[cp]] <- seq_len(sum(rn == cp))
put <- function(f, v) { for (cp in names(mp)) f$estimates[[cp]][mp[[cp]]] <- v[rn == cp]; f }
v0 <- numeric(length(rn)); for (cp in names(mp)) v0[rn == cp] <- fit$estimates[[cp]]
g <- vapply(seq_along(v0), function(j) {
  hj <- 1e-5 * max(1, abs(v0[j])); vp <- v0; vm <- v0
  vp[j] <- vp[j] + hj; vm[j] <- vm[j] - hj
  (avg_sl(put(fit, vp)) - avg_sl(put(fit, vm))) / (2 * hj)
}, 0)
Q <- ns$joint_precision(fit); r <- fit$obj$env$random
Wrr <- as.matrix(Matrix::solve(Q[r, r]))
cat(sprintf("  total %.5f; inner only (theta, beta held, Qrr^-1) %.5f\n",
            sqrt(sum(g * (V %*% g))), sqrt(sum(g[r] * (Wrr %*% g[r])))))
thr <- which(rn == "theta")
cat("  theta:", sprintf("%.4f", v0[thr]), " se:", sprintf("%.4f", sqrt(diag(V)[thr])), "\n")
cat("  d avg slope / d theta at fixed z:", sprintf("%.4f", g[thr]), "\n")
cat("\n4. profile over each theta, inner re-solved:\n")
env <- fit$obj$env
par <- env$last.par.best
cat("  last.par.best equals the estimates:", isTRUE(all.equal(unname(par), v0)), "\n")
outer <- setdiff(seq_along(par), r)
for (k in thr) {
  kk <- match(k, outer)
  sek <- sqrt(V[k, k])
  base_sl <- avg_sl(fit)
  for (s in c(-3, -2, -1, -0.5, 0.5, 1, 2, 3)) {
    p <- par[outer]; p[kk] <- p[kk] + s * sek
    invisible(fit$obj$fn(p))
    full <- env$last.par
    f2 <- put(fit, unname(full))
    # the linearization: d slope(f(theta)) / d theta = g_theta + g_b' J
    cat(sprintf("  theta[%d] %+.1f se: avg slope %.5f (change %+.5f)\n", k, s,
                avg_sl(f2), avg_sl(f2) - base_sl))
  }
  invisible(fit$obj$fn(par[outer]))
}
# the slope of the conditional mode along theta, by central difference,
# times theta's se: the carried part of the delta, one theta at a time
cat("\n5. the slope at single positions (slopes() at x = 1, 3, 5, z = 3):\n")
nd <- data.frame(x = c(1, 3, 5), z = 3)
for (s in names(fits)) {
  sl <- quiet(slopes(fits[[s]], newdata = nd, variables = "x"))
  cat(sprintf("  %-18s %s\n", s, paste(sprintf("%.4f (%.4f)", sl$estimate, sl$std.error), collapse = "  ")))
}
