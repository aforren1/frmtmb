# Reviewer re-check (punch round 1), B2: frmtmb:::logitnormal_sd() on
# extreme inputs, against an independent reference: integrate() over
# the logit in pieces at fixed breakpoints, and a 1e6-draw Monte Carlo
# (seed 1). Also the switch at logit SD 4, NaN/Inf/0 inputs, log = TRUE,
# and a three-component mixture.
#   Rscript dev/surface-rev2-lnsd.R > dev/surface-rev-out/rev2-lnsd.txt
source("dev/surface-rev-env.R")
rev_env("lane")
suppressPackageStartupMessages(library(frmtmb))
lnsd <- frmtmb:::logitnormal_sd
# reference: E and Var of plogis(Z) (or log plogis), Z ~ N(mu, s^2),
# by integrate() on the standard-normal scale in fixed pieces, so the
# integrator cannot step over the logistic's transition
ref <- function(mu, s, log = FALSE) {
  f <- function(u) {
    z <- mu + s * u
    if (log) stats::plogis(z, log.p = TRUE) else stats::plogis(z)
  }
  br <- c(-40, -10, -6, -3, -1, 0, 1, 3, 6, 10, 40)
  # add the transition point u0 = -mu / s and its neighbourhood
  if (s > 0) {
    u0 <- -mu / s
    br <- sort(unique(c(br, u0 + c(-50, -20, -5, -1, 0, 1, 5, 20, 50) / s)))
    br <- br[br >= -40 & br <= 40]
  }
  I <- function(g) sum(vapply(seq_len(length(br) - 1L), function(j) {
    stats::integrate(function(u) g(u) * stats::dnorm(u), br[j], br[j + 1],
                     rel.tol = 1e-12, abs.tol = 0,
                     subdivisions = 2000L, stop.on.error = FALSE)$value
  }, 0))
  m <- I(f)
  sqrt(I(function(u) (f(u) - m)^2))
}
mc <- function(mu, s, log = FALSE, R = 1e6) {
  set.seed(1)
  z <- rnorm(R, mu, s)
  sd(if (log) stats::plogis(z, log.p = TRUE) else stats::plogis(z))
}
cat("p (smaller side), logit SD: frmtmb | reference integrate | MC 1e6\n")
grid <- expand.grid(p = c(0.5, 0.1, 1e-3, 1e-8, 1e-30),
                    s = c(0, 1e-15, 1e-8, 0.01, 1, 3.999, 4.001, 10, 100,
                          1e4, 1e8))
for (i in seq_len(nrow(grid))) {
  p <- grid$p[i]; s <- grid$s[i]
  mu <- qlogis(p)
  v <- lnsd(p, 1 - p, s)
  r <- if (s < 1e-6) p * (1 - p) * s else ref(mu, s)
  m <- if (s < 1e-6) NA else mc(mu, s)
  cat(sprintf("p %-7g s %-7g: %-12.6g | %-12.6g | %-12.6g rel %s\n", p, s,
              v, r, m, if (r > 0) format(v / r - 1, digits = 3) else "-"))
}
cat("\np near 1 read from the complement q (p = 1 - 1e-12, q = 1e-12):\n")
cat("  s = 1:", lnsd(1 - 1e-12, 1e-12, 1), "| reference",
    ref(qlogis(1e-12), 1), "\n")
cat("  s = 10:", lnsd(1 - 1e-12, 1e-12, 10), "| reference",
    ref(qlogis(1e-12), 10), "\n")
cat("\ndelta limit, s -> 0: sd / (p q s) at s = 1e-8:",
    lnsd(0.3, 0.7, 1e-8) / (0.21 * 1e-8), "; at s = 1e-15:",
    lnsd(0.3, 0.7, 1e-15) / (0.21 * 1e-15), "\n")
cat("\nNaN/Inf/0 inputs (p, q, lse):\n")
cases <- list(c(NaN, 0.5, 1), c(0.5, 0.5, NaN), c(0.5, 0.5, Inf),
              c(0, 1, 1), c(1, 0, 1), c(0.5, 0.5, -1), c(NA, 0.5, 1),
              c(0.2, 0.8, 0))
for (cc in cases) {
  r <- tryCatch(format(lnsd(cc[1], cc[2], cc[3])),
                error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("  p %s q %s lse %s -> %s\n", cc[1], cc[2], cc[3], r))
}
cat("\nlog = TRUE (sd of log p):\n")
for (pp in list(c(1e-3, 1), c(0.999, 1), c(1e-8, 10), c(1 - 1e-9, 6))) {
  p <- pp[1]; s <- pp[2]
  cat(sprintf("  p %-9g s %-3g: %-12.6g | ref %-12.6g | MC %-12.6g\n", p, s,
              lnsd(p, 1 - p, s, log = TRUE),
              ref(qlogis(p), s, log = TRUE), mc(qlogis(p), s, log = TRUE)))
}
t0 <- proc.time()
invisible(lnsd(rep(1e-6, 2000), rep(1 - 1e-6, 2000), rep(7, 2000)))
cat("\n2000 wide rows (integrate path):", (proc.time() - t0)[3], "s\n")

## three components
set.seed(8)
d3 <- data.frame(y = c(rnorm(80, -4), rnorm(80, 0), rnorm(80, 4)))
f3 <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian(), gaussian()),
          data = d3)
pm <- pp_mixture(f3)
cat("\nthree components: dim", dim(pm), "|", dimnames(pm)[[3]], "\n")
cat("row sums of Estimate: range", format(range(rowSums(pm[, "Estimate", ]))),
    "\n")
cat("any NA Est.Error:", anyNA(pm[, "Est.Error", ]), "| range",
    format(range(pm[, "Est.Error", ])), "\n")
cat("Q2.5 <= Estimate <= Q97.5 everywhere:",
    all(pm[, "Q2.5", ] <= pm[, "Estimate", ] &
          pm[, "Estimate", ] <= pm[, "Q97.5", ]), "\n")
pl <- pp_mixture(f3, log = TRUE)
cat("log = TRUE: Estimate equals log(Estimate):",
    isTRUE(all.equal(pl[, "Estimate", ], log(pm[, "Estimate", ]))),
    "| any non-finite Est.Error:", any(!is.finite(pl[, "Est.Error", ])),
    "\n")
# two components: the two columns must report the same spread
set.seed(4)
dm <- data.frame(y = c(rnorm(150, -1.5), rnorm(150, 1.5)))
f2 <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()), data = dm)
p2 <- pp_mixture(f2)
cat("two components: max relative |Est.Error[1] - Est.Error[2]|",
    format(max(abs(p2[, "Est.Error", 1] / p2[, "Est.Error", 2] - 1))), "\n")
