# Short check of fix round 1, item 3: is a binomial at 1e5 trials reading -444.00
# at 1e-6 above -444.85 at the reference a property of the model, or of one seed?
# Independent of the tape: the field integrated by Laplace written out here, then
# importance-sampled on the Laplace Gaussian, at each nugget, at each model's own
# maximum. Seed 5601 is the row the lane read; 12 more seeds say how often it
# happens. Under the lane's reading, white noise raises the likelihood "where the
# field cannot follow the data at a known noise".
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-binom.R binom.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
marg <- function(p, x, suc, tr, nu, is_draws = 0L) {
  u <- sort(unique(x)); J <- match(x, u); n <- length(x)
  s2 <- exp(2 * p[2])
  K <- s2 * exp(-outer(u, u, "-")^2 / (2 * exp(2 * p[3]))) + diag(nu * s2, length(u))
  R <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(R)) return(1e10)
  # Newton on the field, the mode of the joint, in K's own well-conditioned form
  f <- rep(0, length(u))
  ll <- function(f) {
    eta <- p[1] + f[J]
    sum(dbinom(suc, tr, plogis(eta), log = TRUE))
  }
  jl <- function(f) {
    a <- backsolve(R, f, transpose = TRUE)
    ll(f) - 0.5 * sum(a^2) - sum(log(diag(R))) - 0.5 * length(u) * log(2 * pi)
  }
  for (it in 1:80) {
    pr <- plogis(p[1] + f[J])
    g <- as.numeric(tapply(suc - tr * pr, factor(J, seq_along(u)), sum))
    w <- as.numeric(tapply(tr * pr * (1 - pr), factor(J, seq_along(u)), sum))
    B <- K %*% diag(w, length(u)) + diag(length(u))
    fn <- as.numeric(solve(B, K %*% (w * f + g)))
    st <- 1
    old <- jl(f)
    while (jl(f + st * (fn - f)) < old - 1e-12 && st > 1e-6) st <- st / 2
    d <- max(abs(st * (fn - f)))
    f <- f + st * (fn - f)
    if (d < 1e-10) break
  }
  pr <- plogis(p[1] + f[J])
  w <- as.numeric(tapply(tr * pr * (1 - pr), factor(J, seq_along(u)), sum))
  # log|K^-1 + W| + log|K| = log|I + K W| through the symmetric form
  sw <- sqrt(w)
  S <- diag(length(u)) + outer(sw, sw) * K
  LS <- chol(S)
  lap <- jl(f) + sum(log(diag(R))) - sum(log(diag(LS)))
  if (!is_draws) return(-lap)
  # the posterior covariance K - K sw S^-1 sw K
  V <- K - K %*% (sw * chol2inv(LS)) %*% (sw * K)
  V <- (V + t(V)) / 2
  Lv <- t(chol(V + diag(1e-14 * mean(diag(V)), nrow(V))))
  set.seed(77)
  lw <- vapply(seq_len(is_draws), function(i) {
    z <- rnorm(length(u))
    fs <- f + as.numeric(Lv %*% z)
    jl(fs) - (-0.5 * sum(z^2) - sum(log(diag(Lv))) - 0.5 * length(u) * log(2 * pi))
  }, 0)
  mx <- max(lw)
  w0 <- exp(lw - mx)
  list(is = mx + log(mean(w0)), lap = lap, ess = sum(w0)^2 / sum(w0^2))
}
one <- function(seed, tr, label) {
  set.seed(seed)
  n <- 60
  x <- round(runif(n, 0, 6), 3)
  invisible(rnorm(n))
  z <- round(runif(n, 0, 6), 3)
  suc <- rbinom(n, tr, plogis(sin(2 * x)))
  starts <- list(c(0, 0.7, 0.3), c(0, 1, 0.2), c(0, 0.3, 0))
  fit <- function(nu) {
    best <- NULL
    for (st in starts) {
      o <- nlminb(st, marg, x = x, suc = suc, tr = tr, nu = nu,
                  control = list(iter.max = 300, eval.max = 500))
      if (is.null(best) || o$objective < best$objective) best <- o
    }
    best
  }
  a <- fit(1e-6)
  b <- fit(1e-12)
  ia <- marg(a$par, x, suc, tr, 1e-6, 4000L)
  ib <- marg(b$par, x, suc, tr, 1e-12, 4000L)
  ic <- marg(a$par, x, suc, tr, 1e-12, 4000L)
  cat(sprintf(paste0("BIN %s seed %d tr %g | 1e-6: lap %.4f is %.4f (ess %.0f) sdgp %.3f rho %.3f | ",
                     "1e-12: lap %.4f is %.4f (ess %.0f) sdgp %.3f rho %.3f | ",
                     "1e-6 less 1e-12, lap %.3f is %.3f | at 1e-6's point read at 1e-12: lap %.4f\n"),
              label, seed, tr, -a$objective, ia$is, ia$ess, exp(a$par[2]), exp(a$par[3]),
              -b$objective, ib$is, ib$ess, exp(b$par[2]), exp(b$par[3]),
              b$objective - a$objective, ia$is - ib$is, ic$lap))
}
# the lane's row: the seed, draw order and trials of dev/review-gpgap-silent.R
set.seed(5601)
n <- 60
x0 <- round(runif(n, 0, 6), 3); e0 <- rnorm(n); z0 <- round(runif(n, 0, 6), 3)
for (seed in 6001:6012) one(seed, 1e5, "sim")
for (seed in 6001:6006) one(seed, 1e4, "sim")
for (seed in 6001:6006) one(seed, 1e6, "sim")
