# Diagnose test-difference.R:490 on the base build: where do the two
# routes to the per-row standard error differ, and by how much?
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r6"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(8)
n <- 120
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
fit <- suppressWarnings(
  frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)),
      family = stats::gaussian(), data = d))
print(fit$estimates$theta)
cat("se_lost:\n"); print(frmtmb:::sdr_of(fit)$se_lost)
cat("pdHess:", frmtmb:::sdr_of(fit)$pdHess, "\n")
jc <- frmtmb:::get_joint_cov(fit)
cat("jc$null is NULL:", is.null(jc$null), "\n")
ev <- eigen(jc$V, symmetric = TRUE, only.values = TRUE)$values
cat("V eigen range:", format(range(ev), digits = 17), "\n")
gx <- d$x[-1] - diff(d$x) / 2
A <- data.frame(x = gx, fac = factor("A", levels = levels(d$fac)))
B <- data.frame(x = gx, fac = factor("B", levels = levels(d$fac)))
for (g in list(A = A, B = B)) {
  lb <- frm_lp_basis(fit, newdata = g, extra_cov = TRUE)
  C <- as.matrix(lb$A)
  q1 <- rowSums((C %*% lb$V) * C)
  q2 <- diag(C %*% lb$V %*% t(C))
  cat("coef var (rowSums) range:", format(range(q1), digits = 6), "\n")
  cat("coef var (triple) range:", format(range(q2), digits = 6), "\n")
  cat("extra_var range:", format(range(lb$extra_var), digits = 6), "\n")
  cat("min(q + extra):", format(min(q2 + lb$extra_var), digits = 6), "\n")
  cat("n rows q<0:", sum(q2 < 0), " q+e<0:", sum(q2 + lb$extra_var < 0),
      "\n")
  cat("se_nonest:", sum(lb$se_nonest), "\n")
  ref <- frm_linpred(fit, newdata = g, se.fit = TRUE)$se.fit
  mine <- sqrt(pmax(q2 + lb$extra_var, 0))
  cat("max rel:", format(max(abs(mine / ref - 1)), digits = 6), "\n")
}
r <- try(frm_curve(fit, newdata = A, contrast = B, simultaneous = FALSE))
print(attr(r, "check"))
