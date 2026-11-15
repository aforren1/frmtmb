# Short check, item 7: core NEWS says the margin is "the fit over the
# observed rows with the term's matrix cut to them". True for fcor()
# (M_oo); for sar() the margin of the observed rows is not sar(W_oo).
# Lane only.
source("dev/reviews/arspat2-short-common.R")
W <- rook(7)
W <- W / rowSums(W)
n <- nrow(W)
set.seed(21)
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.6 * W, 1 + 0.5 * d$x + rnorm(n)))
miss <- sort(sample(n, 12))
ok <- !(seq_len(n) %in% miss)
dm <- d
dm$y[miss] <- NA
for (ty in c("lag", "error")) {
  fo <- stats::as.formula(paste0("y | mi() ~ x + sar(W, type = '", ty, "')"))
  fc <- stats::as.formula(paste0("y ~ x + sar(Wo, type = '", ty, "')"))
  a <- quiet(frm(fo, data = dm, data2 = list(W = W)))
  b <- quiet(frm(fc, data = d[ok, ], data2 = list(Wo = W[ok, ok])))
  # dense reference of the margin at a's optimum, gaussian
  p <- a$opt$par
  ac <- a$frame$autocor[[1L]]
  rho <- unname(ns$autocor_natural(p[names(p) == "thetaac"], ac))
  A <- diag(n) - rho * W
  B <- if (ty == "lag") solve(A) else diag(n)
  sig <- exp(p[names(p) == "betad"])
  mu <- as.vector(B %*% (cbind(1, d$x) %*% p[names(p) == "beta"]))
  S <- sig^2 * solve(crossprod(A))[ok, ok]
  L <- chol(S)
  r <- backsolve(L, d$y[ok] - mu[ok], transpose = TRUE)
  ll_ref <- -0.5 * (sum(ok) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(r * r))
  cat(sprintf("sar %s: mi fit logLik %.6f (dense margin at its optimum %.6f); sar(W_oo) on observed rows %.6f; rho %.4f vs %.4f\n",
              ty, as.numeric(logLik(a)), ll_ref, as.numeric(logLik(b)),
              rho, unname(ns$autocor_natural(b$opt$par[names(b$opt$par) == "thetaac"],
                                             b$frame$autocor[[1L]]))))
}
