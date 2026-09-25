# Lane arcov: frmtmb's cov = FALSE ARMA against brms's COMPILED Stan
# program, through rstan::log_prob() at frmtmb's estimates.
#
#   FRMTMB_STAN_CACHE=/tmp/lanes/arcov/stan-cache Rscript dev/arcov-brms-stan.R
#
# The translator and the checks are test-brms-likelihood.R's own
# (tests/testthat/helper-brms.R): brms_flat_prior(), stan_pars_from_fit()
# and the same two quantities brms_lp_check() asserts, printed rather
# than asserted, for the models dev/arcov-validate.R checks against the
# R transliteration. Ragged groups with interior gaps, rows shuffled.
.libPaths(c("/opt/rlib/stan", "/opt/rlib/lane-arcov", "/opt/rlib/base",
            "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages({
  library(testthat)
  library(brms)
  library(frmtmb)
})
# the helper runs inside the package namespace under testthat, and so
# does it here
h <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-brms.R", envir = h)

set.seed(29)
ng <- 25
nt <- 8
d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                time = rep(seq_len(nt), ng))
d$x <- rnorm(nrow(d))
d$y <- 1 + 0.7 * d$x + unlist(lapply(seq_len(ng), function(i) {
  as.numeric(arima.sim(list(ar = 0.6, ma = 0.3), nt, sd = 0.6))
}))
set.seed(31)
d <- d[-c(3, 17, 18, 40, 77, 150), ]
d$w <- runif(nrow(d), 0.5, 2)
d$cc <- sample(c(0, 0, 0, 1, -1), nrow(d), TRUE)
d$yp <- d$y + 4
d <- d[sample(nrow(d)), ]

one <- function(label, fo_frm, fo_brms, family, joint = FALSE) {
  fit <- frm(fo_frm, data = d, family = family)
  fam_b <- brms::brmsfamily(family[["family"]])
  prior <- brms_flat_prior(fo_brms, data = d, family = fam_b)
  code <- brms::make_stancode(fo_brms, data = d, family = fam_b,
                              prior = prior)
  sdat <- brms_standata(fo_brms, data = d, family = fam_b, prior = prior)
  rtab <- brms_ranef_table(fo_brms, d, fam_b, prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code), data = sdat,
                                         chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code, rtab)
  up <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
  ours <- if (joint) {
    -fit$obj$env$f(fit$obj$env$last.par.best) + attr(pars, "logJ")
  } else {
    as.numeric(logLik(fit))
  }
  gr <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
  if (joint) gr <- gr[brms_inner_index(sf, pars)]
  data.frame(model = label, stan_log_prob = lp, frmtmb = ours,
             rel_diff = abs(lp - ours) / abs(ours),
             max_abs_grad = max(abs(gr)))
}
environment(one) <- h

out <- rbind(
  one("ar(1)", y ~ x + ar(time, g), brms::bf(y ~ x + ar(time, g)),
      gaussian()),
  one("ma(1)", y ~ x + ma(time, g), brms::bf(y ~ x + ma(time, g)),
      gaussian()),
  one("arma(2,1)", y ~ x + arma(time, g, p = 2, q = 1),
      brms::bf(y ~ x + arma(time, g, p = 2, q = 1)), gaussian()),
  one("student arma(1,1)", y ~ x + arma(time, g),
      brms::bf(y ~ x + arma(time, g)), student()),
  one("weights ma(1)", y | weights(w) ~ x + ma(time, g),
      brms::bf(y | weights(w) ~ x + ma(time, g)), gaussian()),
  one("cens arma(1,1)", y | cens(cc) ~ x + arma(time, g),
      brms::bf(y | cens(cc) ~ x + arma(time, g)), gaussian()),
  one("trunc ma(1)", yp | trunc(lb = 0) ~ x + ma(time, g),
      brms::bf(yp | trunc(lb = 0) ~ x + ma(time, g)), gaussian()),
  one("(1 | g) + arma(1,1), joint", y ~ x + (1 | g) + arma(time, g),
      brms::bf(y ~ x + (1 | g) + arma(time, g)), gaussian(), joint = TRUE))
old <- options(width = 120, digits = 12)
print(out, row.names = FALSE)
options(old)
