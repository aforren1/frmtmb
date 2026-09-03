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
reference <- function(file, n_at = NULL) {
  owd <- setwd(EX)
  on.exit(setwd(owd))
  src <- file.path(EX, file)
  if (!is.null(n_at)) {
    # Run the reference at a shorter series than it publishes. Through
    # v0.42 the transform rows were verified this way, because ar1()
    # was dense; the substitution is kept so the audit's own numbers
    # stay reproducible next to the full-size ones.
    txt <- sub("^n <- 1000$", paste("n <-", n_at), readLines(src))
    src <- textConnection(txt)
    on.exit(close(src), add = TRUE, after = FALSE)
  }
  invisible(capture.output(suppressWarnings(
    source(src, local = FALSE, echo = FALSE))))
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

# The multivariate stochastic volatility model at its PUBLISHED size,
# n = 945 and p = 3. Reachable since v0.43, when ar1() stopped building
# a dense 945 x 945 covariance.
#
# The two optimizers do not land on the same point (the reference's own
# nlminb returns convergence = 1 from the example's single cold start),
# so the identity is established the way the matern row establishes its
# own: the REFERENCE objective is cross-evaluated at frmtmb's estimates,
# mapped through the parameterization below. Agreement there is a
# statement about the two likelihoods; agreement of the two optima would
# only be a statement about two optimizers.
#
# The map. frmtmb's sigma has a log link, so
# log sigma_y = mu_x / 2 + h / 2: the intercept is mu_x / 2 and the
# ar1() block IS h / 2. Hence
#
#   mu_x_j      = 2 * betad_j
#   phi_j       = t_j / sqrt(1 + t_j^2)          (theta position 2)
#   sigma_init  = 2 * exp(theta_1)               (marginal sd of h)
#   log_sigma_j = log(sigma_init sqrt(1 - phi^2))  (INNOVATION sd)
#   off_diag_x  = thetar                         (same row-normalized L)
sdv_map_par <- function(fit, p) {
  th <- fit$estimates$theta
  bd <- fit$estimates$betad
  t2 <- th[seq_len(p) * 2L]
  phi <- t2 / sqrt(1 + t2^2)
  sigma_init <- 2 * exp(th[seq_len(p) * 2L - 1L])
  c(phi, log(sigma_init * sqrt(1 - phi^2)), 2 * bd,
    fit$estimates$thetar)
}

sdv_case <- function(script) {
  function() {
    reference(script)
    dd <- data.frame(x1 = y[, 1], x2 = y[, 2], x3 = y[, 3],
                     tim = factor(seq_len(n), levels = seq_len(n)),
                     g = factor(rep(1L, n)))
    t0 <- Sys.time()
    fit <- frm(mvbf(bf(x1 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
                    bf(x2 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
                    bf(x3 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
                    rescor = TRUE), data = dd)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ll <- as.numeric(logLik(fit))
    cross <- obj$fn(sdv_map_par(fit, p))
    cat(sprintf("  %s n=%d p=%d  frm fit %.1f s, reference nlminb conv=%d\n",
                script, n, p, el, opt$convergence))
    cat(sprintf("  ref nll at its own optimum = %.8f\n", opt$objective))
    cat(sprintf("  ref nll at frm's estimates = %.8f  (frm -logLik = %.8f)\n",
                cross, -ll))
    cat(sprintf("  identity |ref(frm par) - (-logLik)| = %.3e\n",
                abs(cross + ll)))
    # and the reference, restarted from frmtmb's point, stays there:
    # that is the reference agreeing at its OWN converged optimum
    lo <- c(rep(-.99, p), rep(-3, p), rep(-3, p), rep(-5, p))
    hi <- c(rep(.99, p), rep(3, p), rep(3, p), rep(5, p))
    warm <- suppressWarnings(
      stats::nlminb(sdv_map_par(fit, p), obj$fn, obj$gr,
                    lower = lo, upper = hi))
    cat(sprintf("  ref restarted at frm's point: nll=%.8f conv=%d\n",
                warm$objective, warm$convergence))
    list(ref = -warm$objective, fit = fit)
  }
}

cases$sdv_multi <- sdv_case("sdv_multi.R")
cases$sdv_multi_compact <- sdv_case("sdv_multi_compact.R")

# The two transform rows, at the example's own n = 1000. Through v0.42
# they were verified at n = 200 because ar1() was dense; the O(d)
# density removes the cap. Both need the pin the audit describes: the
# reference field is marginally standard normal (dautoreg's own scale),
# and frmtmb's ar1() carries a free marginal sd, so the plain spelling
# fits a superset with one extra parameter.
# `body` is a quoted bf() call: the cases are built before the package
# is loaded, so it cannot be evaluated here.
transform_case <- function(script, body, n_at = NULL, map_par = identity) {
  function() {
    reference(script, n_at)
    bfo <- eval(body)
    d <- data.frame(y = y, tim = factor(seq_len(n), levels = seq_len(n)),
                    g = factor(rep(1L, n)))
    # the reference's own start: both shape parameters at 1. A
    # nonlinear body with qgamma()/qbeta() in it has no gradient at
    # zero, so this is not optional on either side.
    st <- list(beta = c(0, 1, 1))
    free <- frm(bfo, data = d, start = st)
    pin <- c(theta_1 = 0)
    fixed <- suppressWarnings(frm(bfo, data = d, start = st,
                                  lower = pin, upper = pin))
    th <- fixed$estimates$theta
    cat(sprintf("  %s n=%d  free logLik=%.10f (df=%s)  pinned=%.10f\n",
                script, n, as.numeric(logLik(free)),
                attr(logLik(free), "df"), as.numeric(logLik(fixed))))
    cat(sprintf("  reference phi=%.6f   pinned rho=%.6f\n",
                opt$par[["phi"]], th[2] / sqrt(1 + th[2]^2)))
    # The reference objective at frmtmb's pinned estimates, in the
    # reference's own parameterization (phi, the two shape parameters,
    # the residual sd). Both scripts run nlminb from a cold start and
    # transform2's Laplace inner problem is delicate at sd = 0.005, so
    # the identity is stated here rather than left to two optimizers.
    par_at <- c(th[2] / sqrt(1 + th[2]^2),
                map_par(fixed$estimates$beta[2:3]), sigma(fixed))
    cross <- obj$fn(par_at)
    cat(sprintf("  ref nll at its optimum=%.8f  at frm's pinned par=%.8f",
                opt$objective, cross))
    cat(sprintf("  |identity|=%.3e\n", abs(cross + as.numeric(logLik(fixed)))))
    list(ref = -opt$objective, fit = fixed)
  }
}

# `scale` is the audit's own name for the third argument, and RTMB's
# qgamma() takes a RATE there positionally, so the fitted coefficient is
# 1 / scale. The likelihood is the same one either way (it is the same
# gamma), which is why the row's logLik matched all along; only the
# cross-evaluation below has to undo the reciprocal.
cases$transform <- transform_case(
  "transform.R",
  quote(bf(y ~ RTMB::qgamma(RTMB::pnorm(z), shape, scale),
           z ~ 0 + ar1(tim + 0 | g), shape ~ 1, scale ~ 1, nl = TRUE)),
  map_par = function(b) c(b[1], 1 / b[2]))

cases$transform2 <- transform_case(
  "transform2.R",
  quote(bf(y ~ RTMB::qbeta(RTMB::pnorm(z), shape1, shape2),
           z ~ 0 + ar1(tim + 0 | g), shape1 ~ 1, shape2 ~ 1, nl = TRUE)))

# The v0.42 spellings, at the n = 200 the dense ar1() forced. Kept so
# that the audit's published numbers stay reproducible against the O(d)
# density: they must not move.
cases$transform_n200 <- transform_case(
  "transform.R",
  quote(bf(y ~ RTMB::qgamma(RTMB::pnorm(z), shape, scale),
           z ~ 0 + ar1(tim + 0 | g), shape ~ 1, scale ~ 1, nl = TRUE)),
  n_at = 200L, map_par = function(b) c(b[1], 1 / b[2]))

cases$transform2_n200 <- transform_case(
  "transform2.R",
  quote(bf(y ~ RTMB::qbeta(RTMB::pnorm(z), shape1, shape2),
           z ~ 0 + ar1(tim + 0 | g), shape1 ~ 1, shape2 ~ 1, nl = TRUE)),
  n_at = 200L)

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
