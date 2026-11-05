# D8: L1 identity against brms 2.23.0's log_prob, per family that takes an
# mi() response, on both shapes: the response itself and a response read
# by an mi() predictor of a second gaussian response. Also the Laplace
# error of the chart: frm()'s logLik and estimates against the exact
# marginal, the complete-case fit (a missing value that no predictor
# reads integrates to 1). Seed 169 + family index.
# Usage: Rscript dev/refuse-mi-l1.R <family>
fam_name <- commandArgs(trailingOnly = TRUE)[1L]
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
fam_i <- match(fam_name, fams)
stopifnot(!is.na(fam_i))
set.seed(169 + fam_i)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
eta <- 0.5 + 0.4 * d$x
d$y <- switch(fam_name,
  gaussian = rnorm(n, eta, 1),
  student = eta + rt(n, 4),
  skew_normal = eta + rnorm(n) + abs(rnorm(n)),
  beta = rbeta(n, plogis(eta) * 8, (1 - plogis(eta)) * 8),
  Gamma = rgamma(n, 3, 3 / exp(eta)),
  weibull = rweibull(n, 2, exp(eta)),
  exponential = rexp(n, 1 / exp(eta)),
  frechet = 1 / rweibull(n, 3, exp(-eta)) * 1,
  inverse.gaussian = {
    mu <- exp(eta) ; v <- rnorm(n)^2
    yy <- mu + mu^2 * v / (2 * 4) -
      mu / (2 * 4) * sqrt(4 * mu * 4 * v + mu^2 * v^2)
    ifelse(runif(n) < mu / (mu + yy), yy, mu^2 / yy)
  },
  lognormal = rlnorm(n, eta, 0.6),
  exgaussian = eta + rnorm(n, 0, 0.5) + rexp(n, 1),
  gen_extreme_value = eta + 0.6 * (rexp(n)^(-0.1) - 1) / 0.1,
  von_mises = {
    a <- 2 * atan(0.5 * d$x)
    ((a + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi
  },
  asym_laplace = eta + rexp(n) - rexp(n))
d_full <- d
rows <- sample(n, 15)
d$y[rows] <- NA
d$w <- d$z + rnorm(n, 0, 0.3)
frm_fam <- function() {
  switch(fam_name,
    beta = frmtmb::Beta(),
    Gamma = Gamma(link = "log"),
    inverse.gaussian = inverse.gaussian(link = "log"),
    get(fam_name, asNamespace("frmtmb"))())
}
brm_fam <- function() {
  if (fam_name == "beta") brms::brmsfamily("beta") else
  if (fam_name == "Gamma") brms::brmsfamily("gamma", "log") else
  if (fam_name == "inverse.gaussian") {
    brms::brmsfamily("inverse.gaussian", "log")
  } else brms::brmsfamily(fam_name)
}
out <- list(family = fam_name)
report <- function(tag, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat(tag, ": ERROR ", substr(conditionMessage(r), 1, 300), "\n", sep = "")
  } else {
    cat(tag, ": ", paste(format(r, digits = 8), collapse = " "), "\n", sep = "")
  }
  invisible(r)
}
ctl <- list()
fit <- report("fit uni", {
  f <- frm(bf(y | mi() ~ x), family = frm_fam(), data = d)
  fit_keep <<- f
  as.numeric(logLik(f))
})
if (exists("fit_keep")) {
  f <- fit_keep
  cc <- report("complete-case logLik", {
    g <- frm(bf(y ~ x), family = frm_fam(), data = d[!is.na(d$y), ])
    cc_keep <<- g
    as.numeric(logLik(g))
  })
  if (exists("cc_keep")) {
    b1 <- frmtmb:::fixef_by_dpar(f)$mu
    b2 <- frmtmb:::fixef_by_dpar(cc_keep)$mu
    se <- sqrt(diag(vcov(cc_keep)))[seq_along(b2)]
    cat("laplace error: logLik frm - exact", as.numeric(logLik(f)) -
          as.numeric(logLik(cc_keep)), "per missing",
        (as.numeric(logLik(f)) - as.numeric(logLik(cc_keep))) / length(rows),
        "\nmu coef diff in SE units:", format((b1 - b2) / se, digits = 3), "\n")
  }
  report("L1 uni", {
    bform <- brms::bf(y | mi() ~ x)
    pr <- E$brms_flat_prior(bform, data = d, family = brm_fam())
    code <- brms::make_stancode(bform, data = d, family = brm_fam(), prior = pr)
    sdat <- E$brms_standata(bform, data = d, family = brm_fam(), prior = pr)
    mod <- E$brms_stan_model(code)
    sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
    pars <- E$stan_pars_from_fit(f, sdat, code, NULL)
    up <- rstan::unconstrain_pars(sf, pars)
    lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
    ours <- -f$obj$env$f(f$obj$env$last.par.best) + attr(pars, "logJ")
    g <- E$brms_lp_grad(sf, up, f, pars)
    idx <- E$brms_inner_index(sf, pars)
    c(const = lp - ours, inner_grad = max(abs(g[idx])), ours = ours)
  })
}
# a family's mi() variable read by a gaussian response's predictor
d2 <- d
d2$resp <- d2$z + 0.5 * d_full$y + rnorm(n, 0, 0.5)
d2$y <- d$y
rm(list = intersect("fit_keep", ls()))
report("fit pred", {
  f2 <- frm(mvbf(bf(resp ~ mi(y) + z), bf(y | mi() ~ x, family = frm_fam()),
                 rescor = FALSE), data = d2)
  fit2_keep <<- f2
  as.numeric(logLik(f2))
})
if (exists("fit2_keep")) {
  report("L1 pred", {
    fb <- brms::bf(resp ~ mi(y) + z) + brms::bf(y | mi() ~ x, family = brm_fam()) +
      brms::set_rescor(FALSE)
    pr <- E$brms_flat_prior(fb, data = d2, family = NULL)
    code <- brms::make_stancode(fb, data = d2, prior = pr)
    sdat <- E$brms_standata(fb, data = d2, prior = pr)
    mod <- E$brms_stan_model(code)
    sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
    f2 <- fit2_keep
    pars <- E$stan_pars_from_fit(f2, sdat, code, NULL)
    up <- rstan::unconstrain_pars(sf, pars)
    lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
    ours <- -f2$obj$env$f(f2$obj$env$last.par.best) + attr(pars, "logJ")
    g <- E$brms_lp_grad(sf, up, f2, pars)
    idx <- E$brms_inner_index(sf, pars)
    c(const = lp - ours, inner_grad = max(abs(g[idx])), ours = ours)
  })
}
