# Audit harness for the classic TMB examples.
#
# For each example that frmtmb claims to cover, run the RTMB reference
# script from dev/tmb-examples/ and the one-formula frm() spelling in the
# SAME process, then compare the maximized log-likelihoods. Both sides
# tape the same objective through RTMB, so a disagreement is a
# modeling-surface difference, not numerical noise.
#
#   Rscript dev/tmb-examples-check.R          # every case
#   Rscript dev/tmb-examples-check.R spatial  # one case
#
# The reference scripts source() their data files into the global
# environment, so each case gets its own R subprocess.

self <- normalizePath(sub("^--file=", "",
                          grep("^--file=", commandArgs(), value = TRUE)[1]),
                      winslash = "/")
root <- normalizePath(file.path(dirname(self), ".."), winslash = "/")
EX <- file.path(root, "dev", "tmb-examples")

# Run the reference script the way it expects: cwd on the example folder
# (it source()s its data file by relative path) and assignments landing
# in the global environment.
reference <- function(file) {
  owd <- setwd(EX)
  on.exit(setwd(owd))
  invisible(capture.output(suppressWarnings(
    source(file.path(EX, file), local = FALSE, echo = FALSE))))
}

cases <- list()

cases$linreg <- function() {
  reference("linreg.R")
  d <- data.frame(Y = data$Y, x = data$x)
  list(ref = -opt$value,
       fit = frm(bf(Y ~ x), family = gaussian(), data = d))
}

cases$dataeval <- function() {
  reference("dataeval.R")
  # The example is about tape reuse: a, b and sd are shared across the
  # ten data chunks, so the fitted model is one pooled regression.
  list(ref = -opt$objective,
       fit = frm(bf(y ~ x), family = gaussian(), data = df))
}

cases$tweedie <- function() {
  reference("tweedie.R")
  list(ref = -opt$objective,
       fit = frm(bf(y ~ 1), family = tweedie(), data = data.frame(y = data$y)))
}

cases$compois <- function() {
  reference("compois.R")
  # frmtmb's compois() is RTMB::dcompois2, the MEAN parameterization,
  # which is the example's second fit.
  list(ref = -fit.mean$objective,
       fit = frm(bf(x ~ 1), family = compois(), data = data.frame(x = data$x)))
}

cases$spatial <- function() {
  reference("spatial.R")
  # X[, 1] is the constant 0.1, so column one is an intercept on a
  # rescaled coefficient. exp(pos + 0 | grp) is Sigma = sd^2 exp(-D/rho),
  # the reference's exp(log_sigma) * u with u ~ MVN(0, exp(-a D)) and
  # rho = 1/a.
  d <- data.frame(y = y, x2 = X[, 2],
                  pos = num_factor(Z[, 1], Z[, 2]),
                  grp = factor(rep(1L, length(y))))
  # The reference drops the Poisson normalizing constant (it sums the
  # kernel y*eta - exp(eta)); put it back to compare log-likelihoods.
  list(ref = -opt$objective - sum(lfactorial(y)),
       fit = frm(bf(y ~ 1 + x2 + exp(pos + 0 | grp)),
                 family = poisson(), data = d))
}

cases$adaptive_integration <- function() {
  reference("adaptive_integration.R")
  # One scalar random effect per observation, marginalized by adaptive
  # quadrature on both sides.
  d <- data.frame(x = x, n = n, c1 = c1, c2 = c2, c3 = c3,
                  obs = factor(seq_along(x)))
  list(ref = -fit$objective,
       fit = frm(bf(x | trials(n) ~ 0 + c1 + c2 + c3 + (1 | obs)),
                 family = binomial(), data = d, quadrature = TRUE))
}

cases$spde <- function() {
  reference("spde.R")
  # Weibull survival with right censoring over a Matern SPDE field.
  # frmtmb's spde() assembles the same
  # Q = tau^2 (kappa^4 c0 + 2 kappa^2 g1 + g2).
  d <- data.frame(time = LeukSurv$time,
                  cens = 1 - LeukSurv$cens,
                  sex = LeukSurv$sex, age = LeukSurv$age,
                  wbc = LeukSurv$wbc, tpi = LeukSurv$tpi,
                  node = factor(mesh$idx$loc, levels = seq_len(mesh$n)))
  fem <- list(c0 = spde$c0, g1 = spde$g1, g2 = spde$g2)
  list(ref = -opt$objective,
       fit = frm(bf(time | cens(cens) ~ sex + age + wbc + tpi +
                      spde(fem, gr = node)),
                 family = weibull(), data = d))
}

# NEAR, not a replication. mat() carries the same Matern correlation,
# but the reference field has unit marginal variance and no nugget, and
# frmtmb's spelling has neither switch: the plain formula fits a
# 4-parameter superset (sd, range, shape, residual sigma). Pinning the
# two extra parameters through lower/upper on the internal theta names
# recovers the reference's phi and kappa, which is the evidence that the
# kernel itself agrees.
cases$matern <- function() {
  reference("matern.R")
  set.seed(123)                     # the coordinates the reference used
  xy <- matrix(runif(2 * 100, 0, 10), 100)
  d <- data.frame(z = as.vector(x), pos = num_factor(xy[, 1], xy[, 2]),
                  grp = factor(rep(1L, 100)))
  free <- frm(bf(z ~ 0 + mat(pos + 0 | grp)), family = gaussian(), data = d)
  pin <- c(theta_1 = 0, `sigma_(Intercept)` = -20)
  fixed <- suppressWarnings(
    frm(bf(z ~ 0 + mat(pos + 0 | grp)), family = gaussian(), data = d,
        lower = pin, upper = pin))
  th <- fixed$estimates$theta
  med <- stats::median(D[D > 0])
  cat(sprintf("  matern free spelling  logLik=%.6f  sd=%.4f  sigma=%.2e\n",
              as.numeric(logLik(free)), sqrt(VarCorr(free)[[1]][1, 1]),
              sigma(free)))
  cat(sprintf("  reference phi=%.6f kappa=%.6f\n", fit$par[1], fit$par[2]))
  cat(sprintf("  pinned    phi=%.6f kappa=%.6f\n",
              0.05 * med + exp(th[2]), 0.1 + 4.9 / (1 + exp(-th[3]))))
  # cross-evaluate: the reference objective at frmtmb's estimates
  cat(sprintf("  ref nll at ref par=%.8f  at pinned par=%.8f\n",
              obj$fn(fit$par),
              obj$fn(c(0.05 * med + exp(th[2]), 0.1 + 4.9 / (1 + exp(-th[3]))))))
  list(ref = -fit$objective, fit = fixed)
}

args <- commandArgs(trailingOnly = TRUE)
want <- if (length(args)) args else names(cases)

if (length(args) == 1L && args[1] %in% names(cases)) {
  suppressMessages(pkgload::load_all(root, quiet = TRUE))
  out <- cases[[args[1]]]()
  ll <- if (is.null(out$fit)) NA_real_ else as.numeric(logLik(out$fit))
  cat(sprintf("%-22s ref=%.10f  frm=%.10f  |diff|=%.3e\n",
              args[1], out$ref, ll, abs(out$ref - ll)))
} else {
  rscript <- file.path(R.home("bin"), "Rscript")
  for (nm in want) {
    st <- system2(rscript, c(shQuote(self), nm), stdout = TRUE, stderr = TRUE)
    hit <- grep("\\|diff\\||^  ", st, value = TRUE)
    cat(if (length(hit)) hit else
        sprintf("%-22s FAILED: %s\n", nm, paste(utils::tail(st, 3),
                                                collapse = " ")),
        sep = "\n")
  }
}
