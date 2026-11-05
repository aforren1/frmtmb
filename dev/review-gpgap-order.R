# Short check of fix round 1, item 1: does the ORDER of the factor change how
# much rounding moves the field? At the target (sdgp 2.17, lscale 0.2539 of dmax,
# nugget 1e-12 sdgp^2), one kernel in four orders, double against 300 bits, the same
# z. The unit is the noise: true sigma 1e-3. A field moved by a good part of sigma
# makes the sampler's gradient noise. Also the pivots, and how much of z the
# likelihood can see. Seed 6601 for x; z's seed 45.
#
#   Rscript dev/review-gpgap-order.R
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(Rmpfr))
prec <- 300
mchol <- function(K) {
  n <- nrow(K)
  L <- mpfrArray(0, prec, dim = c(n, n))
  for (j in seq_len(n)) {
    s <- K[j, j]
    if (j > 1) for (q in seq_len(j - 1)) s <- s - L[j, q]^2
    L[j, j] <- sqrt(s)
    if (j < n) for (i in (j + 1):n) {
      s <- K[i, j]
      if (j > 1) for (q in seq_len(j - 1)) s <- s - L[i, q] * L[j, q]
      L[i, j] <- s / L[j, j]
    }
  }
  L
}
maximin <- function(x) {
  p <- which.min(abs(x - mean(x)))
  while (length(p) < length(x)) {
    d <- vapply(seq_along(x), function(i) min(abs(x[i] - x[p])), 0)
    p <- c(p, which.max(d))
  }
  p
}
for (n in c(15L, 40L)) {
  set.seed(6601)
  x <- round(runif(n, 0, 6), 3)
  sdg <- 2.17
  rho <- 0.2539 * max(dist(x))
  xs <- sort(x)
  orders <- list(sorted = order(x), first_appearance = seq_len(n),
                 reversed = rev(order(x)), maximin = maximin(x))
  set.seed(45)
  z0 <- rnorm(n)
  for (nm in names(orders)) {
    p <- orders[[nm]]
    xp <- x[p]
    Kd <- sdg^2 * exp(-outer(xp, xp, "-")^2 / (2 * rho^2))
    diag(Kd) <- diag(Kd) + 1e-12 * sdg^2
    Ld <- t(chol(Kd))
    mx <- mpfr(xp, prec)
    Km <- mpfrArray(0, prec, dim = c(n, n))
    for (i in seq_len(n)) for (j in seq_len(n)) {
      Km[i, j] <- mpfr(sdg, prec)^2 * exp(-(mx[i] - mx[j])^2 / (2 * mpfr(rho, prec)^2))
      if (i == j) Km[i, j] <- Km[i, j] + mpfr(1e-12, prec) * mpfr(sdg, prec)^2
    }
    Lm <- mchol(Km)
    fd <- as.numeric(Ld %*% z0)
    fm <- mpfr(rep(0, n), prec)
    for (i in seq_len(n)) {
      t0 <- mpfr(0, prec)
      for (q in seq_len(i)) t0 <- t0 + Lm[i, q] * mpfr(z0[q], prec)
      fm[i] <- t0
    }
    err <- max(abs(as.numeric(mpfr(fd, prec) - fm)))
    # how the entries of the factor themselves err, largest per column
    dL <- abs(Ld - matrix(as.numeric(Lm), n, n))
    pv <- diag(Ld)
    cat(sprintf(paste0("ORDER n %d %-17s field error %.2e (%.2e of sigma 1e-3) | max factor entry error %.2e | ",
                       "pivots: first 3 %s, smallest %.1e, count below 1e-3 %d | the z that hold the signal: ",
                       "pivot share of z1..5 in the field's variance %.3f\n"),
                n, nm, err, err / 1e-3, max(dL),
                paste(signif(pv[1:3], 2), collapse = ","), min(pv), sum(pv < 1e-3),
                sum(Ld[, 1:5]^2) / sum(Ld^2)))
  }
}
