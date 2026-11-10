# Lane support, part B: frmtmb's xbeta log density against brms's
# log_lik_xbeta() (betareg::dxbeta behind it) and against the definition
# taken in 256-bit arithmetic (Rmpfr).
#
# Section 1 rebuilds the prep of brms's block "log_lik for the xbeta
# model runs without errors" (seed 568, ns 50, nobs 8, observation 3),
# as dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R does.
# Section 2 adds the end rows y = 0 and y = 1 and a spread of shapes,
# which the block does not reach.
#
#   Rscript dev/support-xbeta-ref.R   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(Rmpfr)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "betareg",
    format(packageVersion("betareg")), "brms",
    format(packageVersion("brms")), "\n")
prec <- 256
fam <- frmtmb::brmsfamily("xbeta")
fr <- function(y, mu, phi, kappa) {
  n <- length(mu)
  as.numeric(frmtmb::row_lpdf(fam, rep_len(y, n), rep_len(y, n),
                              list(mu = mu, phi = rep_len(phi, n),
                                   kappa = rep_len(kappa, n)),
                              list(), NULL))
}
# log I_x(a, b) in mpfr, from I_x(a, b) = x^a (1 - x)^b / (a B(a, b))
# 2F1(a + b, 1; a + 1; x): every term of the series is positive, so
# nothing cancels, and at x <= 0.4 (kappa / (1 + 2 kappa) with kappa
# below 2) the terms fall below 2^-1000 of the sum well before n = 4000
# for b up to 500. Rmpfr's pbetaI() takes whole shapes only.
log_ibeta_mp <- function(x, a, b, nmax = 4000) {
  n <- mpfr(0:nmax, prec)
  lt <- lgamma(a + b + n) - lgamma(a + b) - lgamma(a + 1 + n) +
    lgamma(a + 1) + n * log(x)
  m <- max(lt)
  stopifnot(asNumeric(lt[length(lt)] - m) < -600)
  a * log(x) + b * log1p(-x) - log(a) -
    (lgamma(a) + lgamma(b) - lgamma(a + b)) + m + log(sum(exp(lt - m)))
}
# the definition: Z ~ Beta(mu phi, (1 - mu) phi), Y = (1 + 2 kappa) Z -
# kappa censored to [0, 1]
ref <- function(y, mu, phi, kappa) {
  n <- max(length(y), length(mu), length(phi), length(kappa))
  y <- rep_len(y, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  kappa <- rep_len(kappa, n)
  vapply(seq_len(n), function(i) {
    a <- mpfr(mu[i], prec) * phi[i]
    b <- (1 - mpfr(mu[i], prec)) * phi[i]
    d <- 1 + 2 * mpfr(kappa[i], prec)
    q <- mpfr(kappa[i], prec) / d
    if (y[i] <= 0) return(asNumeric(log_ibeta_mp(q, a, b)))
    if (y[i] >= 1) return(asNumeric(log_ibeta_mp(q, b, a)))
    z <- (mpfr(y[i], prec) + kappa[i]) / d
    asNumeric((a - 1) * log(z) + (b - 1) * log1p(-z) -
                (lgamma(a) + lgamma(b) - lgamma(a + b)) - log(d))
  }, numeric(1))
}
rel <- function(x, r) max(abs(x - r) / pmax(abs(r), 1e-300))

cat("== section 1: brms's block, seed 568 ==\n")
set.seed(568L)
ns <- 50
nobs <- 8
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(
  mu = matrix(rbeta(ns * nobs, 1.2, 2.3), ncol = nobs),
  phi = rexp(ns, 0.01),
  kappa = rexp(ns, 2)
)
prep$data <- list(Y = rbeta(nobs, 2, 3))
i <- 3
ll_f <- fr(prep$data$Y[i], prep$dpars$mu[, i], prep$dpars$phi,
           prep$dpars$kappa)
ll_b <- brms:::log_lik_xbeta(i, prep = prep)
ll_r <- ref(prep$data$Y[i], prep$dpars$mu[, i], prep$dpars$phi,
            prep$dpars$kappa)
cat(sprintf(paste("n %d; frmtmb vs 256-bit: max rel %.3g; brms vs 256-bit:",
                  "max rel %.3g; frmtmb vs brms: max rel %.3g\n"),
            length(ll_f), rel(ll_f, ll_r), rel(ll_b, ll_r), rel(ll_f, ll_b)))
cat(sprintf("frmtmb vs brms: max abs %.3g; range of values [%.4g, %.4g]\n",
            max(abs(ll_f - ll_b)), min(ll_r), max(ll_r)))

cat("== section 2: end rows and a spread of shapes, seed 581 ==\n")
set.seed(581L)
m <- 100
mu <- rbeta(m, 1, 1)
phi <- exp(runif(m, log(0.5), log(500)))
kappa <- exp(runif(m, log(1e-3), log(2)))
for (yv in c(0, 1, 0.3)) {
  f <- fr(yv, mu, phi, kappa)
  b <- betareg::dxbeta(rep(yv, m), mu = mu, phi = phi, nu = kappa,
                       log = TRUE)
  r <- ref(yv, mu, phi, kappa)
  fin <- is.finite(r)
  cat(sprintf(paste("y = %.1f: n %d (%d with a finite reference, range",
                    "[%.4g, %.4g]); frmtmb vs 256-bit max rel %.3g;",
                    "betareg vs 256-bit max rel %.3g; max abs: frmtmb",
                    "%.3g, betareg %.3g\n"),
              yv, m, sum(fin), min(r[fin]), max(r[fin]),
              rel(f[fin], r[fin]), rel(b[fin], r[fin]),
              max(abs(f[fin] - r[fin])), max(abs(b[fin] - r[fin]))))
}
