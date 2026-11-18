# Lane squant: the population (n -> Inf) fit of the extended log-F
# density to a known residual law. For each smoothing h (lambda = h
# sigma) and quantile tau it reports the bias of mu* against the true
# tau-quantile, in units of the residual sd, and the ratio of the
# model-based (Hessian) variance of mu-hat to the sandwich variance.
# Expectations are midpoint sums over the law's quantile function on
# 4e5 points. No package code; no randomness.
#
# Run: Rscript dev/squant-pop.R > dev/squant-out/pop.txt

elf_lpdf <- function(y, mu, sigma, h, tau) {
  lam <- h * sigma
  r <- y - mu
  t <- r / lam
  l1pe <- ifelse(t > 0, t + log1p(exp(-t)), log1p(exp(t)))
  (1 - tau) * r / sigma - h * l1pe - log(lam) -
    lbeta(h * (1 - tau), h * tau)
}

ls_ <- sqrt((exp(1) - 1) * exp(1))
laws <- list(
  normal = list(d = dnorm, q = qnorm),
  t3 = list(d = function(x) dt(x * sqrt(3), 3) * sqrt(3),
            q = function(p) qt(p, 3) / sqrt(3)),
  lognormal = list(d = function(x) dlnorm(x * ls_ + exp(0.5)) * ls_,
                   q = function(p) (qlnorm(p) - exp(0.5)) / ls_)
)

N <- 400000L
u <- (seq_len(N) - 0.5) / N

pop_fit <- function(law, h, tau) {
  q <- law$q(tau)
  yy <- law$q(u)
  obj <- function(p) -mean(elf_lpdf(yy, p[1], exp(p[2]), h, tau))
  o <- optim(c(q, log(0.3)), obj, method = "BFGS",
             control = list(reltol = 1e-15, maxit = 1000))
  p <- o$par
  e <- 1e-5
  sc <- function(p) {
    cbind((elf_lpdf(yy, p[1] + e, exp(p[2]), h, tau) -
             elf_lpdf(yy, p[1] - e, exp(p[2]), h, tau)) / (2 * e),
          (elf_lpdf(yy, p[1], exp(p[2] + e), h, tau) -
             elf_lpdf(yy, p[1], exp(p[2] - e), h, tau)) / (2 * e))
  }
  s0 <- sc(p)
  B <- crossprod(s0) / N
  A <- matrix(0, 2, 2)
  for (j in 1:2) {
    pp <- p
    pm <- p
    pp[j] <- pp[j] + 1e-4
    pm[j] <- pm[j] - 1e-4
    A[, j] <- colMeans(sc(pp) - sc(pm)) / 2e-4
  }
  A <- (A + t(A)) / 2
  Ai <- solve(-A)
  sand <- Ai %*% B %*% Ai
  qr_var <- tau * (1 - tau) / law$d(q)^2
  c(h = h, tau = tau, q = q, mu = p[1], bias = p[1] - q,
    sigma = exp(p[2]), lambda = h * exp(p[2]), var_hess = Ai[1, 1],
    var_sand = sand[1, 1], var_qr = qr_var,
    ratio_hess_sand = Ai[1, 1] / sand[1, 1], grad = max(abs(colMeans(s0))))
}

hs <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4)
taus <- c(0.5, 0.75, 0.9, 0.95)
cat("law\th\ttau\tbias_sd\tsigma\tlambda\tsd_hess\tsd_sand\tsd_qr\t",
    "hess_over_sand\tbias_over_se_n100\tbias_over_se_n1000\tgrad\n", sep = "")
for (nm in names(laws)) {
  for (tau in taus) {
    for (h in hs) {
      r <- pop_fit(laws[[nm]], h, tau)
      cat(sprintf(paste0("%s\t%.2f\t%.2f\t%+.5f\t%.4f\t%.4f\t%.4f\t%.4f\t",
                         "%.4f\t%.4f\t%.3f\t%.3f\t%.1e\n"),
        nm, h, tau, r[["bias"]], r[["sigma"]], r[["lambda"]],
        sqrt(r[["var_hess"]]), sqrt(r[["var_sand"]]), sqrt(r[["var_qr"]]),
        r[["ratio_hess_sand"]],
        r[["bias"]] / sqrt(r[["var_sand"]] / 100),
        r[["bias"]] / sqrt(r[["var_sand"]] / 1000), r[["grad"]]))
    }
  }
}
