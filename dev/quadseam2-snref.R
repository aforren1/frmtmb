# Item 4836: is the stall stop of the skew normal under quadrature a
# false maximum of the rule, or a real one of the integral? The exact
# marginal log-likelihood by stats::integrate() per group (rel.tol 1e-12)
# at each fit's estimate, and the rule taped at each estimate read at the
# other one.
#   Rscript dev/quadseam2-snref.R <libpath|base> [seed]
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
seed <- if (length(a) > 1) as.integer(a[2]) else 6L
set.seed(seed); n <- 200
xs <- -abs(rnorm(n)) * 3
y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
m3 <- mean((y - mean(y))^3) / sd(y)^3
bad <- list(betad = c(log(sd(y)), 2 * sign(m3) + 0.5 * m3))
fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
# skew normal by mean and sd (RTMBdist::dskewnorm2), written with
# pnorm(log.p = TRUE): the double version underflows to -Inf in the tail
lsn <- function(y, mu, sigma, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  om <- sigma / sqrt(1 - 2 / pi * dl^2)
  xi <- mu - om * dl * sqrt(2 / pi)
  z <- (y - xi) / om
  log(2) - log(om) + dnorm(z, log = TRUE) + pnorm(alpha * z, log.p = TRUE)
}
exact_ll <- function(beta, betad, theta) {
  mu <- beta[1] + beta[2] * dd$xs
  sg <- exp(betad[1]); al <- betad[2]; sb <- exp(theta[1])
  tot <- 0
  for (lv in levels(dd$g)) {
    i <- which(dd$g == lv)
    lf <- function(b) {
      vapply(b, function(bb) {
        sum(lsn(dd$y[i], mu[i] + bb, sg, al)) +
          dnorm(bb, 0, sb, log = TRUE)
      }, 0)
    }
    m <- optimize(function(b) -lf(b), c(-10 * sb - 5, 10 * sb + 5))$minimum
    if (sb < 1e-6) { tot <- tot + lf(0) - dnorm(0, 0, sb, log = TRUE); next }
    c0 <- lf(m)
    v <- integrate(function(b) exp(lf(b) - c0), -Inf, Inf,
                   rel.tol = 1e-12, subdivisions = 1000L)$value
    tot <- tot + c0 + log(v)
  }
  tot
}
fits <- list()
for (st in c("default", "stall")) {
  fits[[st]] <- suppressWarnings(suppressMessages(frm(fo,
    family = skew_normal(), data = dd, quadrature = TRUE,
    start = if (st == "stall") bad)))
}
lap <- suppressWarnings(suppressMessages(frm(fo, family = skew_normal(),
                                             data = dd)))
fits$laplace <- lap
for (nm in names(fits)) {
  f <- fits[[nm]]
  e <- f$estimates
  ex <- exact_ll(e$beta, e$betad, e$theta)
  cat(sprintf("FIT %-8s alpha %.5f ll %.9f exact %.9f err %.3e\n", nm,
              e$betad[2], as.numeric(logLik(f)), ex, as.numeric(logLik(f)) - ex))
}
# the rules taped at each quadrature estimate, read at both points
rf <- lapply(fits[c("default", "stall")], frmtmb:::quad_reanchored_fn)
for (a1 in names(rf)) for (a2 in c("default", "stall")) {
  cat(sprintf("RULE taped at %-8s read at %-8s nll %.9f\n", a1, a2,
              rf[[a1]](fits[[a2]]$opt$par)))
}
# exact profile along the segment from the stall estimate to the default
p0 <- fits$stall$estimates; p1 <- fits$default$estimates
for (t in c(0, 0.25, 0.5, 0.75, 1)) {
  b <- (1 - t) * p0$beta + t * p1$beta
  bd <- (1 - t) * p0$betad + t * p1$betad
  th <- (1 - t) * p0$theta + t * p1$theta
  cat(sprintf("SEG t %.2f exact %.9f\n", t, exact_ll(b, bd, th)))
}
cat("END\n")
for (nm in names(fits)) {
  cat("EST", nm, "theta", format(fits[[nm]]$estimates$theta, digits = 6),
      "betad", format(fits[[nm]]$estimates$betad, digits = 6), "\n")
}
