# Review of lane nlstop: why nl_limit_message() is silent on two limits of
# a body, at the fits of dev/reviews/nlstop-review-limit.R (start at the
# truth): a * (exp(b * x) - 1) on linear data (n 1000, seed 2) and the Hill
# body on step data (n 1000, seed 1). Prints each gate of the check and
# the objective along the ridge.
# Usage: Rscript dev/reviews/nlstop-review-fndebug.R
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
one <- function(lab, f, d, st, ridge) {
  r <- tt(frm(f, data = d, start = list(beta = st)))
  cat("==", lab, "ll", llv(r), "code", r$opt$convergence, "warn", wsum(r),
      "\n  est", format(r$estimates$beta, digits = 6), "\n")
  cat("  nl_limit_message():", is.null(ns$nl_limit_message(r)), "(TRUE = NULL)\n")
  obj <- r$obj; p0 <- r$opt$par; j <- which(names(p0) == "beta")
  g0 <- drop(obj$gr(p0))[j]
  H <- matrix(NA, length(j), length(j))
  for (i in seq_along(j)) {
    h <- 1e-4 * abs(p0[j[i]]); up <- p0; dn <- p0
    up[j[i]] <- up[j[i]] + h; dn[j[i]] <- dn[j[i]] - h
    H[, i] <- (drop(obj$gr(up))[j] - drop(obj$gr(dn))[j]) / (2 * h)
  }
  D <- p0[j]
  Hu <- (H + t(H)) / 2 * outer(D, D) + diag(D * g0, length(j))
  e <- eigen(Hu, symmetric = TRUE)
  v <- e$vectors[, length(j)]
  big <- abs(v) >= 0.1 * max(abs(v))
  cat("  g0", format(g0, digits = 3), "\n  eigen", format(e$values, digits = 3),
      "\n  v", format(v, digits = 3), "big", big, "\n")
  cat("  cond z |D| sqrt(H_ii):", format(abs(D) * sqrt(abs(diag(H))), digits = 3), "\n")
  lf <- log(1000); t <- lf / max(abs(v))
  cat("  quad gate 0.5 lam t^2 =", 0.5 * e$values[length(j)] * t^2,
      "<= 1e3 tol =", 1e3 * 1e-3, "\n")
  for (s in c(1, -1)) {
    q <- p0; q[j] <- p0[j] * exp(s * t * v)
    cat("  side", s, "objective gain", r$opt$objective - obj$fn(q), "\n")
  }
  for (s in ridge) {
    q <- p0; q[j] <- p0[j] * s
    cat(sprintf("  ridge %s: gain %.4g\n", paste(format(s), collapse = ","),
                r$opt$objective - obj$fn(q)))
  }
  obj$fn(p0)
}
set.seed(2); n <- 1000
d <- data.frame(x = runif(n, 0, 3)); d$y <- 2 * d$x + rnorm(n, 0, 0.5)
one("explin", bf(y ~ a * (exp(b * x) - 1), a + b ~ 1, nl = TRUE), d,
    c(2, 0.5), list(c(10, 0.1), c(1000, 0.001), c(0.1, 10)))
set.seed(1)
d <- data.frame(x = c(runif(n / 2, 0.5, 4.5), runif(n / 2, 5.5, 10)))
d$y <- 10 * (d$x > 5) + rnorm(n, 0, 0.5)
one("hillstep", bf(y ~ Emax * x^h / (EC50^h + x^h), Emax + EC50 + h ~ 1,
                   nl = TRUE), d, c(10, 5, 2),
    list(c(1, 1, 10), c(1, 1, 1000), c(1, 1, 0.5)))
cat("DONE\n")
