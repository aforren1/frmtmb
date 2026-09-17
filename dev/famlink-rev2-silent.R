## Recheck, priority 1: failures that used to be LOUD. Models whose
## objective is NaN over a region the optimum sits against, or near it.
## For each, every warning frm() raises, the convergence code and
## message, pdHess, logLik, estimates, SEs and (lane) nonfinite_trials.
## Seed 20260916 per design. Usage: Rscript dev/famlink-rev2-silent.R <lane|base>
ARM <- commandArgs(trailingOnly = TRUE)[1]
source("dev/famlink-rev-common.R")
mk <- function(n = 200) { set.seed(20260916); data.frame(x = runif(n, -2, 2),
                                                        g = factor(rep(1:20, length.out = n))) }
designs <- list(
  pois_identity_wall = function() { d <- mk(); d$y <- rpois(200, pmax(0.05, 0.6 + 1.2 * d$x)); list(y ~ x, d, poisson("identity")) },
  bern_log_wall = function() { d <- mk(); d$y <- rbinom(200, 1, pmin(0.99, exp(-0.2 + 0.5 * d$x))); list(y ~ x, d, bernoulli("log")) },
  bern_identity_wall = function() { d <- mk(); d$y <- rbinom(200, 1, pmin(0.98, pmax(0.02, 0.5 + 0.45 * d$x))); list(y ~ x, d, bernoulli("identity")) },
  binom_identity_trials = function() { d <- mk(); d$n <- 10; d$y <- rbinom(200, 10, pmin(0.99, pmax(0.01, 0.5 + 0.4 * d$x))); list(y | trials(n) ~ x, d, binomial("identity")) },
  gamma_inverse_cross = function() { d <- mk(); mu <- 1 / pmax(0.05, 0.6 + 0.35 * d$x); d$y <- rgamma(200, 3, 3 / mu); list(y ~ x, d, Gamma("inverse")) },
  invgauss_default = function() { d <- mk(); d$y <- statmod_like <- rgamma(200, 20, 20 / exp(0.5 + 0.6 * d$x)); list(y ~ x, d, "inverse.gaussian") },
  invgauss_default_re = function() { d <- mk(); d$y <- rgamma(200, 20, 20 / exp(0.5 + 0.6 * d$x + rnorm(20, 0, 0.3)[d$g])); list(y ~ x + (1 | g), d, "inverse.gaussian") },
  pois_identity_re = function() { d <- mk(); d$y <- rpois(200, pmax(0.05, 1 + 1 * d$x + rnorm(20, 0, 0.3)[d$g])); list(y ~ x + (1 | g), d, poisson("identity")) },
  pois_sqrt = function() { d <- mk(); d$y <- rpois(200, pmax(0, 0.5 + 1 * d$x)^2); list(y ~ x, d, poisson("sqrt")) },
  cox_identity = function() { d <- mk(); d$y <- rexp(200, pmax(0.05, 1 + 0.4 * d$x)); list(y ~ x, d, cox("identity")) },
  zib_identity = function() { d <- mk(); d$n <- 8; d$y <- rbinom(200, 8, pmin(0.95, pmax(0.05, 0.5 + 0.4 * d$x))); list(y | trials(n) ~ x, d, zero_inflated_binomial("identity")) },
  gauss_reml_re = function() { d <- mk(); d$y <- 1 + d$x + rnorm(20)[d$g] + rnorm(200); list(y ~ x + (1 | g), d, gaussian(), REML = TRUE) },
  gauss_sigma_bound = function() { d <- mk(); d$y <- 1 + d$x + rnorm(200, 0, 0.5); list(bf(y ~ x, sigma ~ x), d, gaussian(),
      prior = set_prior("", class = "Intercept", dpar = "sigma", lb = -0.2)) },
  exgauss_bound = function() { d <- mk(); d$y <- rnorm(200, d$x) + rexp(200, 2); list(y ~ x, d, exgaussian(), prior = set_prior("", class = "Intercept", lb = -0.5, ub = 0.5)) },
  pois_identity_optim = function() { d <- mk(); d$y <- rpois(200, pmax(0.05, 0.6 + 1.2 * d$x)); list(y ~ x, d, poisson("identity"), control = frmtmb_control(optimizer = "optim")) },
  pois_identity_restarts = function() { d <- mk(); d$y <- rpois(200, pmax(0.05, 0.6 + 1.2 * d$x)); list(y ~ x, d, poisson("identity"), control = frmtmb_control(restarts = 3)) },
  negbin_identity_re = function() { d <- mk(); d$y <- rnbinom(200, mu = pmax(0.1, 2 + 1.5 * d$x + rnorm(20, 0, 0.4)[d$g]), size = 3); list(y ~ x + (1 | g), d, negbinomial("identity")) }
)
out <- list()
for (nm in names(designs)) {
  spec <- tryCatch(designs[[nm]](), error = function(e) e)
  if (inherits(spec, "error")) { cat(nm, "setup error", conditionMessage(spec), "\n"); next }
  args <- list(formula = spec[[1]], data = spec[[2]], family = spec[[3]])
  extra <- spec[-(1:3)]
  args <- c(args, extra)
  warns <- character()
  fit <- withCallingHandlers(tryCatch(do.call(frm, args), error = function(e) e),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  r <- list(warnings = warns)
  if (inherits(fit, "error")) { r$error <- conditionMessage(fit) } else {
    r$conv <- fit$opt$convergence; r$msg <- fit$opt$message
    r$nonfinite <- fit$opt$nonfinite_trials
    r$ll <- as.numeric(logLik(fit)); r$est <- unlist(fit$estimates)
    vw <- character()
    r$se <- withCallingHandlers(tryCatch(sqrt(diag(vcov(fit))), error = function(e) conditionMessage(e)),
                                warning = function(w) { vw <<- c(vw, conditionMessage(w)); invokeRestart("muffleWarning") })
    r$vcov_warnings <- vw
    dg <- character()
    r$diag <- withCallingHandlers(tryCatch(capture.output(diagnose(fit)), error = function(e) conditionMessage(e)),
                                  warning = function(w) { dg <<- c(dg, conditionMessage(w)); invokeRestart("muffleWarning") })
    r$diag_warn <- dg
    r$print <- capture.output(print(fit)); r$summary <- tryCatch(capture.output(summary(fit)), error = function(e) "")
  }
  out[[nm]] <- r
  nan_w <- sum(grepl("NA/NaN", warns))
  cat(sprintf("%-24s %s conv %s nonfinite %s NaN-warnings %d other-warnings %d\n", nm,
              if (is.null(r$error)) "fit" else paste("ERROR", substr(r$error, 1, 60)),
              format(r$conv), format(r$nonfinite), nan_w, length(warns) - nan_w))
}
saveRDS(out, sprintf("dev/famlink-rev2-silent-%s.rds", ARM))
