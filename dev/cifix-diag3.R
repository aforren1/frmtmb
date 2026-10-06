# On the scan's failing fits, compare the coefficient variance of the
# grid rows under the joint covariance as 0.68.0 forms it (Q^-1, which
# carries the flat outer parameters) and conditional on the lost outer
# parameters (their rows and columns of Q removed before the solve).
# Usage: Rscript dev/cifix-diag3.R <lib or "base"> <n> <seed> [<n> <seed>]
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), "\n")
cases <- matrix(as.integer(args[-1]), ncol = 2, byrow = TRUE)
for (k in seq_len(nrow(cases))) {
  n <- cases[k, 1]; s <- cases[k, 2]
  set.seed(s)
  d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                  fac = factor(rep(c("A", "B"), length.out = n)))
  d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
    stats::rnorm(n, 0, 0.3)
  fit <- suppressWarnings(
    frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)),
        family = stats::gaussian(), data = d))
  cat("\n== n =", n, "seed =", s, " theta:",
      format(fit$estimates$theta, digits = 4), "\n")
  lost <- frmtmb:::sdr_of(fit)$se_lost
  print(lost)
  Q <- frmtmb:::joint_precision(fit)
  rn <- rownames(Q)
  random <- c("b", "miss", if (fit$REML) "beta")
  outer_rows <- which(!rn %in% random)
  om <- frmtmb:::outer_par_names(fit)
  stopifnot(length(outer_rows) == length(om))
  drop <- outer_rows[match(names(lost), om)]
  cat("dropped Q rows:", drop, "(", rn[drop], ")\n")
  Qc <- Q[-drop, -drop]
  ev <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  evc <- eigen(as.matrix(Qc), symmetric = TRUE, only.values = TRUE)$values
  cat("Q eig range:", format(range(ev), digits = 4),
      " Q_cond eig range:", format(range(evc), digits = 4), "\n")
  V <- as.matrix(Matrix::solve(Q))
  Vc <- matrix(0, nrow(Q), ncol(Q))
  Vc[-drop, -drop] <- as.matrix(Matrix::solve(Qc))
  gx <- d$x[-1] - diff(d$x) / 2
  A <- data.frame(x = gx, fac = factor("A", levels = levels(d$fac)))
  lb <- suppressWarnings(frm_lp_basis(fit, newdata = A))
  C <- as.matrix(lb$A)
  p <- lb$coef_pos
  q <- rowSums((C %*% V[p, p]) * C)
  qc <- rowSums((C %*% Vc[p, p]) * C)
  cat("q (0.68.0) range:", format(range(q), digits = 4), "\n")
  cat("q (conditional) range:", format(range(qc), digits = 4), "\n")
  cat("max |q/qc - 1| over rows with q > 0:",
      format(max(abs(q[q > 0] / qc[q > 0] - 1)), digits = 4), "\n")
  cat("extra_var range:", format(range(lb$extra_var), digits = 4), "\n")
  # the same coefficient variance from mgcv's convention: smoothing
  # parameters (and here every variance parameter) held at the estimate
  vp <- which(rn %in% c("theta", "betad"))
  Vfix <- solve(as.matrix(Q)[-vp, -vp])
  pr <- match(p, seq_len(nrow(Q))[-vp])
  qf <- rowSums((C %*% Vfix[pr, pr]) * C)
  cat("q (all outer fixed) range:", format(range(qf), digits = 4), "\n")
}
