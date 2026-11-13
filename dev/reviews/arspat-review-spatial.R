# Review of lane arspat: fcor()/sar() with mi() on student() (Laplace
# against the exact multivariate-t margin), and sar() with a complex
# spectrum (log-determinant against determinant(), rho's interval against
# spatialreg's rule). Output to stdout.
source("dev/reviews/arspat-review-common.R")
suppressPackageStartupMessages(library(frmtmb.spatial))
cat("frmtmb.spatial from", dirname(dirname(system.file(package = "frmtmb.spatial"))), "\n")

sar_rho <- function(fit) {
  ac <- fit$frame$autocor[[1L]]
  unname(ns$autocor_natural(fit$estimates$thetaac[ac$theta_idx], ac))
}

## 1. fcor() + mi() on student(): Laplace vs exact margin
fcor_case <- function(nmiss, seed = 3, n = 60, tdf = 4) {
  set.seed(seed)
  M <- 0.5^abs(outer(1:n, 1:n, "-"))
  d <- data.frame(x = rnorm(n))
  # several independent blocks so nu is identified: one chi-square per 10
  e <- as.vector(t(chol(M)) %*% rnorm(n)) * 0.8 *
    rep(sqrt(tdf / rchisq(n / 10, tdf)), each = 10)
  d$y <- 1 + 0.5 * d$x + e
  miss <- sort(sample(n, nmiss))
  d$y[miss] <- NA
  f <- frm(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M),
           family = student())
  dp <- ns$eval_dpars(f)[[1L]]
  ok <- !is.na(d$y)
  mu <- rep_len(dp$mu, n)
  ex <- mvtnorm::dmvt(d$y[ok], mu[ok], dp$sigma[1]^2 * M[ok, ok],
                      df = dp$nu[1], log = TRUE)
  ll <- as.numeric(logLik(f))
  # the exact maximum: the same model on the observed rows, M cut
  Mo <- M[ok, ok]
  f0 <- frm(y ~ x + fcor(Mo), data = d[ok, ], data2 = list(Mo = Mo),
            family = student())
  dp0 <- ns$eval_dpars(f0)[[1L]]
  cat(sprintf(paste0("fcor student nmiss %2d: at Laplace est: laplace %.6f ",
                     "exact %.6f diff %.3e | max laplace %.6f exact %.6f | ",
                     "nu %.3g vs %.3g, sigma %.5f vs %.5f, b_x %.5f vs %.5f\n"),
              nmiss, ll, ex, ll - ex, ll, as.numeric(logLik(f0)),
              dp$nu[1], dp0$nu[1], dp$sigma[1], dp0$sigma[1],
              fixef(f)[2, 1], fixef(f0)[2, 1]))
  invisible(f)
}
for (k in c(3, 10, 20)) chk(paste("fcor student mi", k), fcor_case(k))
chk("fcor student mi: what the fit says", {
  f <- fcor_case(10)
  out <- capture.output(print(summary(f)))
  cat("summary lines mentioning Laplace or approx:",
      sum(grepl("Laplace|approx", out)), "\n")
  NULL
})

## 2. sar() + mi() on student()
sar_W <- function(n, seed = 5, k = 4, sym = TRUE) {
  set.seed(seed)
  xy <- cbind(runif(n), runif(n))
  D <- as.matrix(dist(xy)); diag(D) <- Inf
  K <- t(apply(D, 1, function(r) as.numeric(rank(r) <= k)))
  if (sym) K <- pmax(K, t(K))
  K / rowSums(K)
}
sar_case <- function(nmiss, type, seed = 4, n = 80, tdf = 4) {
  W <- sar_W(n)
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  A <- diag(n) - 0.4 * W
  e <- rnorm(n, 0, 0.7) * rep(sqrt(tdf / rchisq(n / 10, tdf)), each = 10)
  d$y <- if (type == "lag") as.vector(solve(A, 1 + 0.5 * d$x + e)) else
    1 + 0.5 * d$x + as.vector(solve(A, e))
  miss <- sort(sample(n, nmiss))
  d$y[miss] <- NA
  f <- frm(stats::as.formula(paste0("y | mi() ~ x + sar(W, type = '", type,
                                    "')")),
           data = d, data2 = list(W = W), family = student())
  dp <- ns$eval_dpars(f)[[1L]]
  rho <- sar_rho(f)
  A <- diag(n) - rho * W
  eta <- rep_len(dp$mu, n)
  m <- if (type == "lag") as.vector(solve(A, eta)) else eta
  S <- dp$sigma[1]^2 * solve(crossprod(A))
  ok <- !is.na(d$y)
  ex <- mvtnorm::dmvt(d$y[ok], m[ok], S[ok, ok], df = dp$nu[1], log = TRUE)
  ll <- as.numeric(logLik(f))
  cat(sprintf("sar %-5s student nmiss %2d: laplace %.6f exact %.6f diff %.3e nu %.3g rho %.4f\n",
              type, nmiss, ll, ex, ll - ex, dp$nu[1], rho))
  invisible(f)
}
for (ty in c("lag", "error")) for (k in c(3, 16)) {
  chk(paste("sar student mi", ty, k), sar_case(k, ty))
}

## 3. complex spectrum: log-determinant and rho's interval
chk("complex W: logdet against determinant()", {
  n <- 40
  W <- sar_W(n, seed = 21, k = 3, sym = FALSE)
  ev <- eigen(W, only.values = TRUE)$values
  isreal <- abs(Im(ev)) <= 1e-8 * max(Mod(ev))
  cat("complex eigenvalues:", sum(!isreal), " real:", sum(isreal), "\n")
  cat("min Re over all:", min(Re(ev)), " min real eigenvalue:",
      min(Re(ev[isreal])), "\n")
  cat("lane interval (1/min Re, 1/max Re):", 1 / min(Re(ev)), 1 / max(Re(ev)), "\n")
  cat("spatialreg interval (real eigenvalues only):",
      1 / min(Re(ev[isreal])), 1 / max(Re(ev[isreal])), "\n")
  aux <- frmtmb.spatial:::sar_build(list(M = quote(W), type = "lag"),
                                    seq_len(n), n, function(e, w) W)
  cat("aux lo hi:", aux$lo, aux$hi, "\n")
  for (rho in c(aux$lo * 0.999999, aux$lo * 0.9, -0.5, 0, 0.3, 0.9,
                aux$hi * 0.999999, 1 / min(Re(ev[isreal])) * 0.99)) {
    dt <- determinant(diag(n) - rho * W, logarithm = TRUE)
    ld <- as.numeric(frmtmb.spatial:::sar_logdet(rho, aux))
    cat(sprintf("rho %+.6f logdet %.12f determinant %.12f (sign %d) diff %.2e\n",
                rho, ld, as.numeric(dt$modulus), dt$sign,
                ld - as.numeric(dt$modulus)))
  }
  NULL
})
chk("complex W: a negative rho between the two lower bounds", {
  n <- 40
  W <- sar_W(n, seed = 21, k = 3, sym = FALSE)
  ev <- eigen(W, only.values = TRUE)$values
  isreal <- abs(Im(ev)) <= 1e-8 * max(Mod(ev))
  lo_lane <- 1 / min(Re(ev)); lo_sr <- 1 / min(Re(ev[isreal]))
  rho0 <- if (lo_sr < lo_lane) (lo_lane + max(lo_sr, -3)) / 2 else NA
  cat("true rho", rho0, "\n")
  if (!is.na(rho0)) {
    set.seed(9)
    out <- NULL
    for (r in 1:5) {
      d <- data.frame(x = rnorm(n))
      d$y <- 1 + 0.5 * d$x + as.vector(solve(diag(n) - rho0 * W,
                                             rnorm(n, 0, 0.5)))
      f <- frm(y ~ x + sar(W, type = "error"), data = d, data2 = list(W = W))
      out <- c(out, sar_rho(f))
    }
    cat("fitted rho over 5 draws:", format(out, digits = 5), " lower bound",
        lo_lane, "\n")
  }
  NULL
})
