# Fix round 1, B1: every family with an mi() response no density reads.
# frm() against the complete-case fit (the exact maximum), L1 against brms's
# program on the observed rows, and a read mi() against integrate().
# Usage: Rscript dev/refuse-f1-unread.R <family>. Seed 169 + family index
# (the lane's data) and 9300 + index for the read arm (the review's design).
fam <- commandArgs(trailingOnly = TRUE)[1L]
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
fams <- c("gaussian", "student", "skew_normal", "beta", "Gamma", "weibull",
          "exponential", "frechet", "inverse.gaussian", "lognormal",
          "exgaussian", "gen_extreme_value", "von_mises", "asym_laplace")
fi <- match(fam, fams)
ff <- switch(fam,
  beta = Beta(), Gamma = stats::Gamma(link = "log"),
  inverse.gaussian = stats::inverse.gaussian(link = "log"),
  gaussian = stats::gaussian(),
  get(fam, asNamespace("frmtmb"))())
bf_fam <- function() {
  switch(fam, beta = brms::brmsfamily("beta"),
         Gamma = brms::brmsfamily("gamma", "log"),
         inverse.gaussian = brms::brmsfamily("inverse.gaussian", "log"),
         brms::brmsfamily(fam))
}
gen <- function(n) {
  x <- rnorm(n)
  eta <- 0.3 + 0.7 * x
  y <- switch(fam,
    gaussian = rnorm(n, eta, 0.8),
    student = eta + 0.8 * rt(n, 3),
    skew_normal = eta + 0.5 * rnorm(n) + 1.2 * abs(rnorm(n)),
    beta = rbeta(n, plogis(eta) * 3, (1 - plogis(eta)) * 3),
    Gamma = rgamma(n, 1.3, 1.3 / exp(eta)),
    weibull = rweibull(n, 1.4, exp(eta)),
    exponential = rexp(n, 1 / exp(eta)),
    frechet = 1 / rweibull(n, 2.2, exp(-eta)),
    inverse.gaussian = {
      m <- exp(eta); v <- rnorm(n)^2; lam <- 2
      yy <- m + m^2 * v / (2 * lam) -
        m / (2 * lam) * sqrt(4 * m * lam * v + m^2 * v^2)
      ifelse(runif(n) < m / (m + yy), yy, m^2 / yy)
    },
    lognormal = rlnorm(n, eta, 0.9),
    exgaussian = eta + rnorm(n, 0, 0.3) + rexp(n, 1 / 1.2),
    gen_extreme_value = eta + 0.8 * (rexp(n)^(-0.2) - 1) / 0.2,
    von_mises = ((2 * atan(0.4 * eta) + rnorm(n, 0, 0.7) + pi) %% (2 * pi)) - pi,
    asym_laplace = eta + rexp(n) - rexp(n))
  data.frame(x = x, y = y)
}
cat("family", fam, "\n")
set.seed(169 + fi)
d <- gen(200)
d$y[sample(200, 60)] <- NA
a <- suppressWarnings(frm(bf(y | mi() ~ x), family = ff, data = d))
b <- suppressWarnings(frm(bf(y ~ x), family = ff, data = d[!is.na(d$y), ]))
ea <- unlist(c(a$estimates$beta, a$estimates$betad))
eb <- unlist(c(b$estimates$beta, b$estimates$betad))
cat(sprintf("UNREAD logLik %.6f complete-case %.6f gap %.3g max coef gap %.3g nmiss slots %d BIC gap %.3g\n",
            logLik(a), logLik(b), logLik(a) - logLik(b), max(abs(ea - eb)),
            length(a$estimates$miss), BIC(a) - BIC(b)))
l1 <- tryCatch({
  dd <- E$brms_mi_observed(a, d)
  bform <- brms::bf(y | mi() ~ x)
  pr <- E$brms_flat_prior(bform, data = dd, family = bf_fam())
  code <- brms::make_stancode(bform, data = dd, family = bf_fam(), prior = pr)
  sdat <- E$brms_standata(bform, data = dd, family = bf_fam(), prior = pr)
  mod <- E$brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  pars <- E$stan_pars_from_fit(a, sdat, code, NULL)
  up <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
  g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
  sprintf("const %.3g grad %.3g", lp - as.numeric(logLik(a)), max(abs(g)))
}, error = function(e) paste("L1 error:", substr(conditionMessage(e), 1, 100)))
cat("UNREAD L1", l1, "\n")
set.seed(9300 + fi)
n <- 80
d <- gen(n)
d$z <- 0.5 + 0.6 * d$y + 0.3 * d$x + rnorm(n, 0, 0.7)
mrows <- sample(n, 8)
d$y[mrows] <- NA
fit <- tryCatch(suppressWarnings(frm(
  bf(y | mi() ~ x) + ff + bf(z ~ mi(y) + x) + gaussian(), data = d)),
  error = function(e) e)
if (inherits(fit, "error")) {
  cat("READ refused:", substr(conditionMessage(fit), 1, 120), "\n")
  quit(save = "no")
}
env <- fit$obj$env
pv <- env$last.par.best
mi_i <- which(names(pv) == "miss")
base <- env$f(pv)
err <- numeric(length(mi_i))
for (k in seq_along(mi_i)) {
  u0 <- pv[mi_i[k]]
  h1 <- function(u) {
    vapply(u, function(ui) {
      p <- pv
      p[mi_i[k]] <- ui
      exp(base - env$f(p))
    }, 0)
  }
  eps <- 1e-3
  d2 <- (env$f(replace(pv, mi_i[k], u0 + eps)) +
           env$f(replace(pv, mi_i[k], u0 - eps)) - 2 * base) / eps^2
  w <- 1 / sqrt(max(d2, 1e-6))
  I <- tryCatch(integrate(h1, u0 - 40 * w, u0 + 40 * w, subdivisions = 2000L,
                          rel.tol = 1e-9, stop.on.error = FALSE)$value,
                error = function(e) NA_real_)
  err[k] <- log(sqrt(2 * pi) * w) - log(I)
}
cat(sprintf("READ chart=%s Laplace - exact per value: mean %+.4f range [%+.4f, %+.4f] conv %d\n",
            fit$frame$mi_map$y$chart, mean(err), min(err), max(err),
            fit$opt$convergence))
l2 <- tryCatch({
  fb <- brms::bf(y | mi() ~ x, family = bf_fam()) +
    brms::bf(z ~ mi(y) + x) + brms::set_rescor(FALSE)
  pr <- E$brms_flat_prior(fb, data = d, family = NULL)
  code <- brms::make_stancode(fb, data = d, prior = pr)
  sdat <- E$brms_standata(fb, data = d, family = NULL, prior = pr)
  mod <- E$brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  pars <- E$stan_pars_from_fit(fit, sdat, code, NULL)
  up <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
  ours <- -env$f(pv) + attr(pars, "logJ")
  g <- E$brms_lp_grad(sf, up, fit, pars)
  idx <- E$brms_inner_index(sf, pars)
  sprintf("const %.3g inner grad %.3g", lp - ours, max(abs(g[idx])))
}, error = function(e) paste("L1 error:", substr(conditionMessage(e), 1, 100)))
cat("READ L1", l2, "\n")
