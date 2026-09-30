# Punch 1: the new tests' designs, replayed on the round-1 xbeta density
# (RTMB::dbeta() on a taped argument, and the 50-step Lentz fraction
# with a hard switch at m, both copied here verbatim from the round-1
# source), beside xbeta() on the lane build. This is how the B1 and B2
# pins are seen to fail: the round-1 build itself is no longer installed.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
r1_cf <- function(x, a, b, N) {
  qab <- a + b; qap <- a + 1; qam <- a - 1
  cc <- 1
  d <- 1 / (1 - qab * x / qap)
  lh <- log(abs(d))
  for (k in seq_len(N)) {
    k2 <- 2 * k
    aa <- k * (b - k) * x / ((qam + k2) * (a + k2))
    d <- 1 / (1 + aa * d); cc <- 1 + aa / cc; lh <- lh + log(abs(d * cc))
    aa <- -(a + k) * (qab + k) * x / ((a + k2) * (qap + k2))
    d <- 1 / (1 + aa * d); cc <- 1 + aa / cc; lh <- lh + log(abs(d * cc))
  }
  a * log(x) + b * log1p(-x) - RTMB::lbeta(a, b) - log(a) + lh
}
r1_ibeta <- function(x, a, b, N = 50L) {
  m <- (a + 1) / (a + b + 2)
  t <- m - x
  xd <- x - 0.5 * ((x - m) + abs(x - m))
  xc <- x + 0.5 * ((m - x) + abs(m - x))
  w <- (t + abs(t)) / (2 * abs(t) + 1e-300)
  ld <- r1_cf(xd, a, b, N)
  lc <- r1_cf(1 - xc, b, a, N)
  cap <- log1p(-2^-53)
  lc <- 0.5 * (lc + cap - abs(cap - lc))
  w * ld + (1 - w) * log1p(-exp(lc))
}
r1_xbeta <- frmtmb_family(
  "xbeta_round1", dpars = c("mu", "phi", "kappa"),
  links = list(mu = "logit", phi = "log", kappa = "log"),
  lpdf = function(y, dpars, aterms) {
    i0 <- as.numeric(y <= 0); i1 <- as.numeric(y >= 1); ib <- i0 + i1
    phi <- dpars[["phi"]]; kap <- dpars[["kappa"]]
    mp <- dpar_complement(dpars, "mu", "logit")
    a <- mp$p * phi; b <- mp$q * phi; d <- 1 + 2 * kap
    z <- (y * (1 - ib) + 0.5 * ib + kap) / d
    out <- (1 - ib) * (RTMB::dbeta(z, a, b, log = TRUE) - log(d))
    w0 <- which(i0 == 1); w1 <- which(i1 == 1)
    if (length(w0) || length(w1)) {
      q <- kap / d + 0 * y; a <- a + 0 * y; b <- b + 0 * y
      if (length(w0)) out[w0] <- r1_ibeta(q[w0], a[w0], b[w0])
      if (length(w1)) out[w1] <- r1_ibeta(q[w1], b[w1], a[w1])
    }
    out
  },
  init_dpars = list(mu = function(y, a) min(max(mean(y), 0.05), 0.95),
                    phi = function(y, a) 5, kappa = function(y, a) 0.1))
ref_xbeta <- function(y, mu, phi, kappa) {
  a <- mu * phi; b <- (1 - mu) * phi; d <- 1 + 2 * kappa
  sum(ifelse(y <= 0, stats::pbeta(kappa / d, a, b, log.p = TRUE),
             ifelse(y >= 1, stats::pbeta((1 + kappa) / d, a, b,
                                         lower.tail = FALSE, log.p = TRUE),
                    stats::dbeta((y + kappa) / d, a, b, log = TRUE) -
                      log(d))))
}
fitmsg <- function(expr) {
  tryCatch(suppressWarnings(expr), error = function(e) {
    m <- conditionMessage(e)
    if (grepl("NA/NaN gradient", m)) "ERROR: NA/NaN gradient evaluation" else
      paste("ERROR:", substr(m, 1, 60))
  })
}
cat("== B1: free kappa at phi 1e4 and 2e4 (seed 4200, the new test) ==\n")
for (phi in c(1e4, 2e4)) {
  set.seed(4200)
  z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
  d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
  for (fam in list(r1_xbeta, xbeta())) {
    f <- fitmsg(frm(y ~ 1, family = fam, data = d))
    if (is.character(f)) {
      cat(sprintf("  phi %g %-13s %s\n", phi, fam$family, f))
    } else {
      dp <- frmtmb:::eval_dpars(f)[[1]]
      cat(sprintf("  phi %g %-13s logLik %.8f  reference at the estimates %.8f\n",
                  phi, fam$family, as.numeric(logLik(f)),
                  ref_xbeta(d$y, dp$mu[1], dp$phi[1], dp$kappa[1])))
    }
  }
}
cat("== B2: kappa held at 0.05, phi 2e5 (seed 4200, the new test) ==\n")
set.seed(4200)
z <- stats::rbeta(1000, 0.05 / 1.1 * 2e5, (1 - 0.05 / 1.1) * 2e5)
d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
for (fam in list(r1_xbeta, xbeta())) {
  f <- fitmsg(frm(bf(y ~ 1, kappa = 0.05), family = fam, data = d))
  if (is.character(f)) { cat("  ", fam$family, f, "\n"); next }
  dp <- frmtmb:::eval_dpars(f)[[1]]
  cat(sprintf("  %-13s logLik %.8f  reference at the estimates %.8f  diff %.2e\n",
              fam$family, as.numeric(logLik(f)),
              ref_xbeta(d$y, dp$mu[1], dp$phi[1], 0.05),
              as.numeric(logLik(f)) - ref_xbeta(d$y, dp$mu[1], dp$phi[1], 0.05)))
}
cat("== B2: the incomplete beta at the new test's points, worst over the grid ==\n")
worst <- c(r1 = 0, lane = 0)
for (a in c(0.3, 4, 60, 3e3, 2e5)) for (b in c(0.5, 7, 300, 7e3, 8e5)) {
  s <- a + b; m <- (a + 1) / (s + 2); sd <- sqrt(a * b / (s * s * (s + 1)))
  xs <- c(m, m * (1 - 1e-3), m * (1 + 1e-3), m - 2 * sd, m + 2 * sd,
          m - 5.5 * sd, m + 5.5 * sd, m / 3, 0.45)
  for (x in xs[xs > 0 & xs < 0.5]) {
    r <- stats::pbeta(x, a, b, log.p = TRUE)
    if (!is.finite(r) || r < -600) next
    scale <- abs(a * log(x)) + abs(b * log1p(-x)) + abs(lbeta(a, b))
    tol <- 64 * .Machine$double.eps * max(scale, abs(r))
    worst["r1"] <- max(worst["r1"], abs(r1_ibeta(x, a, b) - r) / tol)
    worst["lane"] <- max(worst["lane"],
                         abs(frmtmb:::log_ibeta_half(x, a, b) - r) / tol)
  }
}
cat(sprintf("  worst error in units of the test's tolerance: round 1 %.3g, lane %.3g\n",
            worst["r1"], worst["lane"]))
