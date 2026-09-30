# Reviewer: where does an xbeta fit's gradient go NaN at high phi?
# Section B of dev/fams2-rev-dens.R (seed 4200, mu 0.047, kappa 0.05,
# phi 2e4) and A1 (interior only, seed 4101).
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(frmtmb); library(RTMB) })
lib <- frmtmb:::log_ibeta_half
cat("== 1. log_ibeta_half gradient in (x, log a, log b) over large shapes ==\n")
G <- MakeTape(function(p) lib(p[1], exp(p[2]), exp(p[3])), c(0.1, 0, 0))
J <- G$jacfun()
bad <- 0; tot <- 0; ex <- NULL
for (la in seq(log(10), log(1e6), length.out = 40)) {
  for (r in c(0.01, 0.03, 0.05, 0.1, 0.3)) {
    a <- exp(la) * r / (1 - r); b <- exp(la)
    m <- (a + 1) / (a + b + 2)
    for (dx in c(-0.3, -0.1, -0.03, -0.01, 0, 0.01, 0.03, 0.1, 0.3)) {
      x <- m * (1 + dx)
      if (x >= 0.5) next
      tot <- tot + 1
      v <- G(c(x, log(a), log(b))); g <- J(c(x, log(a), log(b)))
      if (!is.finite(v) || !all(is.finite(g))) {
        bad <- bad + 1
        if (is.null(ex) || nrow(ex) < 12) ex <- rbind(ex, c(x = x, a = a, b = b, v = v, g))
      }
    }
  }
}
cat("points", tot, " non-finite value or gradient:", bad, "\n")
if (!is.null(ex)) print(ex)
cat("\n== 2. the pieces at one bad point ==\n")
if (!is.null(ex)) {
  x <- ex[1, "x"]; a <- ex[1, "a"]; b <- ex[1, "b"]
  m <- (a + 1) / (a + b + 2)
  cat("x", x, "m", m, "\n")
  cat("direct cf  ", frmtmb:::log_ibeta_cf(min(x, m), a, b, 50L), "\n")
  cat("complement ", frmtmb:::log_ibeta_cf(1 - max(x, m), b, a, 50L), "\n")
  cat("pbeta      ", pbeta(x, a, b, log.p = TRUE), "\n")
  # the Lentz terms along the way
  cf_trace <- function(x, a, b, N = 50) {
    qab <- a + b; qap <- a + 1; qam <- a - 1; cc <- 1
    d <- 1 / (1 - qab * x / qap); out <- d
    for (k in seq_len(N)) {
      k2 <- 2 * k
      aa <- k * (b - k) * x / ((qam + k2) * (a + k2))
      d <- 1 / (1 + aa * d); cc <- 1 + aa / cc; out <- c(out, d * cc)
      aa <- -(a + k) * (qab + k) * x / ((a + k2) * (qap + k2))
      d <- 1 / (1 + aa * d); cc <- 1 + aa / cc; out <- c(out, d * cc)
    }
    out
  }
  tr <- cf_trace(min(x, m), a, b)
  cat("direct Lentz factors <= 0:", sum(tr <= 0), " min |factor|", min(abs(tr)), "\n")
  tr2 <- cf_trace(1 - max(x, m), b, a)
  cat("complement Lentz factors <= 0:", sum(tr2 <= 0), " min |factor|", min(abs(tr2)), "\n")
}
cat("\n== 3. RTMB::dbeta gradient at large shapes ==\n")
D <- MakeTape(function(p) RTMB::dbeta(p[1], exp(p[2]), exp(p[3]), log = TRUE), c(0.3, 0, 0))
DJ <- D$jacfun()
for (s in 10^(2:8)) {
  p <- c(0.3, log(0.3 * s), log(0.7 * s))
  cat(sprintf("shape sum %.0e: value %.6g grad %s\n", s, D(p),
              paste(format(DJ(p), digits = 5), collapse = " ")))
}
cat("\n== 4. the fit's own objective along a phi path (section B, phi 2e4) ==\n")
set.seed(4200)
n <- 1000; mu <- 0.047; kap <- 0.05; phi <- 2e4
z <- rbeta(n, mu * phi, (1 - mu) * phi)
y <- pmin(pmax((1 + 2 * kap) * z - kap, 0), 1)
# evaluate the likelihood tape directly through a custom RTMB objective that
# calls the family lpdf
fam <- xbeta()
lpdf <- frmtmb:::as_frmtmb_family(fam)$lpdf
obj <- MakeADFun(function(p) {
  -sum(lpdf(y, list(mu = plogis(p$m), phi = exp(p$lp), kappa = exp(p$lk)), list()))
}, list(m = qlogis(mu), lp = log(phi), lk = log(kap)), silent = TRUE)
for (lp in log(c(1e3, 1e4, 2e4, 5e4, 1e5, 1e6))) {
  for (lk in log(c(0.03, 0.05, 0.1))) {
    p <- c(qlogis(mu), lp, lk)
    cat(sprintf("phi %.0e kappa %.2f: nll %.6f grad %s\n", exp(lp), exp(lk),
                obj$fn(p), paste(format(obj$gr(p), digits = 5), collapse = " ")))
  }
}
o <- try(nlminb(c(qlogis(0.2), log(5), log(0.1)), obj$fn, obj$gr,
            control = list(trace = 1, iter.max = 300)))
if (!inherits(o, "try-error")) print(o[c("par", "objective", "convergence", "message")])
cat("\n== 5. which piece: boundary rows alone, interior rows alone ==\n")
for (part in c("boundary", "interior")) {
  yy <- if (part == "boundary") y[y == 0] else y[y > 0 & y < 1]
  ob <- MakeADFun(function(p) {
    -sum(lpdf(yy, list(mu = plogis(p$m), phi = exp(p$lp), kappa = exp(p$lk)), list()))
  }, list(m = qlogis(mu), lp = log(phi), lk = log(kap)), silent = TRUE)
  for (ph in c(1e3, 2e3, 5e3, 2e4)) {
    p <- c(qlogis(mu), log(ph), log(kap))
    cat(sprintf("  %s rows (%d), phi %.0e: nll %.6f grad %s\n", part, length(yy), ph,
                ob$fn(p), paste(format(ob$gr(p), digits = 5), collapse = " ")))
  }
}
cat("\n== 6. Beta() interior density through MakeADFun at the same shapes ==\n")
yi <- y[y > 0 & y < 1]
ob <- MakeADFun(function(p) -sum(RTMB::dbeta(yi, plogis(p$m) * exp(p$lp),
                                           (1 - plogis(p$m)) * exp(p$lp), log = TRUE)),
                list(m = qlogis(mu), lp = log(phi)), silent = TRUE)
for (ph in c(1e3, 2e3, 5e3, 2e4)) cat(sprintf("  phi %.0e: grad %s\n", ph,
  paste(format(ob$gr(c(qlogis(mu), log(ph))), digits = 5), collapse = " ")))
