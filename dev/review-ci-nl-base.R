# nl = TRUE on every core family: the body reaches the likelihood.
#
# Through 0.78.0 an ordinal family's rule that thresholds replace the
# intercept also took the intercept of `b1 ~ 1`, so frm() fitted the
# thresholds alone without a word: the maximum equalled the intercept-only
# fit's (dev/equiv-repro-thresnl.R, seed 49). A body honored has the linear
# fit's maximum, up to a small share of what the predictor gains over the
# intercept-only model. The file holds that for each family the table says
# works, and each row of frm_compat("nl") it covers.

nl_fam_data <- function(seed = 20261029, n = 300) {
  set.seed(seed)
  x <- round(rnorm(n), 3)
  eta <- 0.3 + 0.5 * x
  list(x = x, n = n, eta = eta,
       pos = exp(0.4 + 0.3 * x + rnorm(n, 0, 0.4)),
       cnt = rpois(n, exp(eta)), prob = plogis(eta),
       zero = runif(n) < 0.2,
       ord = cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = 1:4),
       e1 = rnorm(n), e2 = rnorm(n, 0, 0.5), e3 = rexp(n, 2))
}

nl_fam_response <- function(f, z) {
  n <- z$n
  switch(f,
    gaussian = , student = , skew_normal = , gen_extreme_value = z$eta + z$e1,
    exgaussian = z$eta + z$e2 + z$e3,
    lognormal = , Gamma = , weibull = , exponential = , inverse.gaussian = ,
    frechet = z$pos,
    shifted_lognormal = 0.1 + z$pos,
    beta = , xbeta = plogis(z$eta + z$e2),
    zero_inflated_beta = ifelse(z$zero, 0, plogis(z$eta + z$e2)),
    zero_one_inflated_beta = ifelse(z$zero, rbinom(n, 1, 0.5),
                                    plogis(z$eta + z$e2)),
    tweedie = , hurdle_gamma = , hurdle_lognormal = ifelse(z$zero, 0, z$pos),
    poisson = , negbinomial = , negbinomial2 = , nbinom1 = , geometric = ,
    compois = , discrete_weibull = z$cnt,
    zero_inflated_poisson = , zero_inflated_negbinomial =
      ifelse(z$zero, 0L, z$cnt),
    hurdle_poisson = , hurdle_negbinomial = ifelse(z$zero, 0L, z$cnt + 1L),
    binomial = , beta_binomial = rbinom(n, 10, z$prob),
    zero_inflated_binomial = , zero_inflated_beta_binomial =
      ifelse(z$zero, 0L, rbinom(n, 10, z$prob)),
    bernoulli = rbinom(n, 1, z$prob),
    cumulative = , sratio = , cratio = , acat = ordered(z$ord),
    hurdle_cumulative = ordered(ifelse(z$zero, 0L, as.integer(z$ord)),
                                levels = 0:4),
    von_mises = (z$eta + z$e2 + pi) %% (2 * pi) - pi,
    stop("no generator for ", f))
}

nl_fam_of <- function(f) {
  if (identical(f, "beta")) return(Beta())
  g <- get0(f, envir = asNamespace("frmtmb"), mode = "function")
  if (is.null(g)) g <- get0(f, envir = asNamespace("stats"), mode = "function")
  if (is.null(g)) brmsfamily(f) else g()
}

nl_fam_lhs <- function(f) {
  if (f %in% c("binomial", "beta_binomial", "zero_inflated_binomial",
               "zero_inflated_beta_binomial")) "y | trials(ntr)" else "y"
}

nl_fam_ll <- function(form, f, d, ...) {
  fit <- allow_warnings(
    suppressMessages(frm(form, family = nl_fam_of(f), data = d, ...)),
    allowed = c("Standard errors", "did not report convergence",
                "standard error"))
  as.numeric(logLik(fit))
}

# One family per line, so that a family's failure names it. The families
# whose body brms sees differently (several locations) are below.
nl_fams_single <- c(
  "gaussian", "student", "lognormal", "shifted_lognormal", "skew_normal",
  "exgaussian", "Gamma", "weibull", "exponential", "inverse.gaussian", "beta",
  "tweedie", "poisson", "negbinomial", "negbinomial2", "nbinom1", "geometric",
  "compois", "binomial", "bernoulli", "beta_binomial",
  "zero_inflated_beta_binomial", "zero_inflated_poisson",
  "zero_inflated_negbinomial", "zero_inflated_binomial", "zero_inflated_beta",
  "zero_one_inflated_beta", "xbeta", "hurdle_poisson", "hurdle_negbinomial",
  "hurdle_gamma", "hurdle_lognormal", "cumulative", "hurdle_cumulative",
  "sratio", "cratio", "acat", "von_mises", "frechet", "gen_extreme_value",
  "discrete_weibull")

for (f_i in nl_fams_single) local({
  f <- f_i
  test_that(paste0("the nl body reaches the likelihood of ", f), {
    z <- nl_fam_data()
    d <- data.frame(x = z$x, ntr = 10L)
    d$y <- nl_fam_response(f, z)
    L <- nl_fam_lhs(f)
    lin <- nl_fam_ll(bf(stats::as.formula(paste(L, "~ x"))), f, d)
    int <- nl_fam_ll(bf(stats::as.formula(paste(L, "~ 1"))), f, d)
    nl <- nl_fam_ll(bf(stats::as.formula(paste(L, "~ b0 + b1 * x")),
                       b0 + b1 ~ 1, nl = TRUE), f, d)
    expect_gt(lin - int, 0)
    expect_lt(abs(nl - lin) / (lin - int), 1e-3)
    st <- frm_compat("nl", f)
    expect_identical(st$status[st$feature_b == f], "works")
  })
})

test_that("an ordinal nlpar keeps its whole design", {
  # the absent case of the fix: thresholds still replace the location's
  # own intercept, and an nlpar's intercept is an ordinary coefficient
  z <- nl_fam_data()
  d <- data.frame(x = z$x, z = rnorm(z$n))
  d$y <- ordered(z$ord)
  lin <- frm(bf(y ~ x), family = cumulative(), data = d)
  expect_false("(Intercept)" %in% names(lin$estimates[["beta"]]))
  fit <- frm(bf(y ~ b1 * x, b1 ~ 1 + z, nl = TRUE), family = cumulative(),
             data = d)
  expect_named(fit$estimates[["beta"]], c("b1_(Intercept)", "b1_z"))
  # with thres() too
  fit2 <- frm(bf(y | thres(3) ~ b1 * x, b1 ~ 1, nl = TRUE),
              family = cumulative(), data = d)
  expect_lt(abs(as.numeric(logLik(fit2)) - as.numeric(logLik(lin))) /
              abs(as.numeric(logLik(lin))), 1e-8)
})

test_that("the post-fit paths of an nl ordinal fit read the body", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x, zz = rnorm(z$n))
  d$y <- ordered(z$ord)
  nd <- d[c(3, 7, 11, 19), ]
  for (fam in list(cumulative(), sratio(), cratio(), acat())) {
    fl <- frm(bf(y ~ x + zz), family = fam, data = d)
    fn <- frm(bf(y ~ bx * x + bz * zz, bx + bz ~ 1, nl = TRUE), family = fam,
              data = d)
    scale_of <- max(abs(fitted(fl)))
    expect_lt(max(abs(fitted(fn) - fitted(fl))) / scale_of, 1e-5)
    expect_lt(max(abs(fitted(fn, newdata = nd) - fitted(fl, newdata = nd))) /
                scale_of, 1e-5)
    expect_lt(max(abs(frm_linpred(fn) - frm_linpred(fl))) /
                max(abs(frm_linpred(fl))), 1e-5)
    set.seed(5)
    sn <- simulate(fn, nsim = 1)
    set.seed(5)
    sl <- simulate(fl, nsim = 1)
    # the two fits' coefficients agree to the optimizer's stop, so a draw
    # that lands on a category boundary may differ, and only that
    expect_gt(mean(as.character(sn[[1]]) == as.character(sl[[1]])), 0.98)
  }
})

test_that("an nl ordinal fit keeps the na.exclude padding", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x)
  d$y <- ordered(z$ord)
  d$x[c(2, 5)] <- NA
  fit <- suppressMessages(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE),
                              family = cumulative(), data = d,
                              na.action = stats::na.exclude))
  p <- fitted(fit)
  expect_identical(nrow(p), nrow(d))
  expect_true(all(is.na(p[c(2, 5), , ])))
  expect_false(anyNA(p[-c(2, 5), , ]))
  # a body on the plain path pads once too
  d$g <- z$eta + z$e1
  g <- suppressMessages(frm(bf(g ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                            family = gaussian(), data = d,
                            na.action = stats::na.exclude))
  expect_length(frm_linpred(g), nrow(d))
})

test_that("the prior table of an nl ordinal model lists its thresholds", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x)
  d$y <- ordered(z$ord)
  pt <- default_prior(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d,
                      family = cumulative())
  thres <- pt[pt$class == "Intercept", ]
  expect_identical(thres$coef, c("", "1", "2", "3"))
  # and a class = "Intercept" prior reaches them, as on the linear model
  pr <- set_prior("normal(0, 1)", class = "Intercept")
  flat <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(),
              data = d)
  fit <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(),
             data = d, prior = pr)
  expect_false(isTRUE(all.equal(fit$estimates[["tau_raw"]],
                                flat$estimates[["tau_raw"]])))
})

test_that("cs() under a body is refused in brms's words", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x, zz = rnorm(z$n))
  d$y <- ordered(z$ord)
  expect_error(frm(bf(y ~ b1 * x, b1 ~ 1 + cs(zz), nl = TRUE),
                   family = sratio(), data = d),
               "only supported for the main parameter 'mu'", fixed = TRUE)
  expect_error(frm(bf(y ~ b1 * x + cs(zz), b1 ~ 1, nl = TRUE),
                   family = sratio(), data = d),
               "cs() is not supported in a nonlinear formula body",
               fixed = TRUE)
  # the absent case: cs() on a linear location still fits
  expect_s3_class(frm(bf(y ~ x + cs(zz)), family = sratio(), data = d),
                  "frmtmb_fit")
  expect_identical(frm_compat("cs_pred()", "nl")$status, "refused")
})

# ---------------------------------------------------------------------
# Several locations: brms writes the one body into each, with one set of
# nonlinear parameters. The reference is the density written by hand at
# the estimates; test-nl-brms.R holds them against brms's program.
# ---------------------------------------------------------------------

nl_multi_data <- function(seed = 5001, n = 200) {
  set.seed(seed)
  d <- data.frame(x = round(rnorm(n), 3), ntr = 10L)
  d$cat <- factor(vapply(seq_len(n), function(i) {
    sample(c("a", "b", "c"), 1, prob = c(1, exp(0.4 * d$x[i]), 1))
  }, ""))
  m <- t(vapply(seq_len(n), function(i) {
    stats::rmultinom(1, 10, c(1, exp(0.4 * d$x[i]), 1))[, 1]
  }, numeric(3)))
  colnames(m) <- c("a", "b", "c")
  d$cnt <- m
  s <- (m + 0.5) / 11.5
  colnames(s) <- c("a", "b", "c")
  d$sim <- s
  d$yg <- ifelse(runif(n) < 0.4, -1, 1) + 0.5 * d$x + rnorm(n, 0, 0.5)
  d
}

test_that("categorical and multinomial read one body in every location", {
  d <- nl_multi_data()
  fit <- frm(bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
             family = categorical(), data = d)
  b <- fit$estimates[["beta"]]
  eta <- unname(b["b0_(Intercept)"] + b["b1_(Intercept)"] * d$x)
  lp <- cbind(0, eta, eta)
  p <- exp(lp) / rowSums(exp(lp))
  by_hand <- sum(log(p[cbind(seq_len(nrow(d)), as.integer(d$cat))]))
  expect_lt(abs(as.numeric(logLik(fit)) - by_hand) / abs(by_hand), 1e-10)
  # fitted() is the same softmax, against what the probabilities themselves
  # span in this fit
  expect_lt(max(abs(unname(fitted(fit)[, "Estimate", ]) - p)) /
              diff(range(p)), 1e-7)

  fit2 <- frm(bf(cnt | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
              family = multinomial(), data = d)
  b <- fit2$estimates[["beta"]]
  eta <- unname(b["b0_(Intercept)"] + b["b1_(Intercept)"] * d$x)
  lp <- cbind(0, eta, eta)
  p <- exp(lp) / rowSums(exp(lp))
  by_hand <- sum(vapply(seq_len(nrow(d)), function(i) {
    stats::dmultinom(d$cnt[i, ], prob = p[i, ], log = TRUE)
  }, 0))
  expect_lt(abs(as.numeric(logLik(fit2)) - by_hand) / abs(by_hand), 1e-10)
})

test_that("a mixture reads one body in every component", {
  d <- nl_multi_data()
  fit <- allow_warnings(frm(bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                            family = mixture(gaussian(), gaussian()), data = d,
                            start = list(beta = c(0.1, 0.4))),
                        allowed = c("Standard errors", "standard error",
                                    "did not report convergence",
                                    "the likelihood is flat"))
  expect_named(fit$estimates[["beta"]], c("b0_(Intercept)", "b1_(Intercept)"))
  # the two components share a mean, so the body moves both
  mu <- frm_linpred(fit, dpar = "mu1")
  expect_identical(unname(mu), unname(frm_linpred(fit, dpar = "mu2")))
})

test_that("nl = TRUE on several locations without a brms form is refused", {
  # the case the guard on the copy leaves: any family with several
  # locations that brms has no such spelling for. A family made here, not
  # frmtmb.latent's hmm(): core's checks install no extension, and the
  # skip gate fails a test that skips on one
  local_frmtmb_registries()
  two_fam <- frmtmb_family("tfamtwomu", dpars = c("mu1", "mu2", "sigma"),
    links = list(mu1 = "identity", mu2 = "identity", sigma = "log"),
    primary_dpars = c("mu1", "mu2"),
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      RTMB::dnorm(y, dpars$mu1 + dpars$mu2, dpars$sigma, log = TRUE)
    })
  d <- data.frame(y = rnorm(40), x = rnorm(40))
  expect_error(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = two_fam,
                   data = d),
               "single 'mu' location", fixed = TRUE)
})

# ---------------------------------------------------------------------
# The guard: a body the likelihood does not read is refused.
# ---------------------------------------------------------------------

test_that("a family that does not read the body is refused", {
  local_frmtmb_registries()
  skip_fam <- frmtmb_family("tfamskip", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      RTMB::dnorm(y, 0, dpars$sigma, log = TRUE)
    })
  read_fam <- frmtmb_family("tfamread", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      RTMB::dnorm(y, dpars$mu, dpars$sigma, log = TRUE)
    })
  set.seed(7101)
  d <- data.frame(x = rnorm(120))
  d$y <- 0.5 * d$x + rnorm(120)
  expect_error(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                   family = skip_fam, data = d),
               "does not reach the likelihood of family 'tfamskip'",
               fixed = TRUE)
  # the absent case: the same family with mu read fits, and at the linear
  # maximum
  fit <- frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = read_fam,
             data = d)
  lin <- frm(bf(y ~ x), family = read_fam, data = d)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(lin))) /
              abs(as.numeric(logLik(lin))), 1e-8)
  # and a linear formula on the family that skips mu is not its business
  expect_s3_class(allow_warnings(frm(bf(y ~ 1), family = skip_fam, data = d),
                                 allowed = "Standard errors"), "frmtmb_fit")
})

test_that("a parameter whose design came out empty is refused", {
  # Hollow the design as the ordinal rule did through 0.78.0: nothing in
  # this build can reach it, so the frame is edited by hand
  set.seed(7102)
  d <- data.frame(x = rnorm(60))
  d$y <- 0.5 * d$x + rnorm(60)
  fit <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = gaussian(), data = d)
  fr <- fit$frame
  nll <- build_objective(fr)
  tpl <- fr[["par_template"]]
  # the guard passes on the real frame
  expect_null(nl_check_body_read(fr, nll, tpl))
  key <- linpred_key("y", "b1")
  fr[["linpreds"]][[key]][["X"]] <- fr[["linpreds"]][[key]][["X"]][, 0L,
                                                                   drop = FALSE]
  fr[["linpreds"]][[key]][["idx"]] <- integer(0)
  expect_error(nl_check_body_read(fr, nll, tpl),
               "left it with no coefficient", fixed = TRUE)
})

test_that("the guard leaves a parameter with no coefficients on purpose", {
  set.seed(7103)
  d <- data.frame(x = rnorm(80), o = 0.5)
  d$y <- 1 + 0.5 * d$x + rnorm(80)
  # b ~ 0: brms builds it, and the body reads 0
  f0 <- frm(bf(y ~ a + b * x, a ~ 1, b ~ 0, nl = TRUE), family = gaussian(),
            data = d)
  f1 <- frm(bf(y ~ 1), family = gaussian(), data = d)
  expect_lt(abs(as.numeric(logLik(f0)) - as.numeric(logLik(f1))) /
              abs(as.numeric(logLik(f1))), 1e-8)
  # an offset-only parameter
  f2 <- frm(bf(y ~ a + b * x, a ~ 1, b ~ 0 + offset(o), nl = TRUE),
            family = gaussian(), data = d)
  f3 <- frm(bf(y ~ 1 + offset(0.5 * x)), family = gaussian(), data = d)
  expect_lt(abs(as.numeric(logLik(f2)) - as.numeric(logLik(f3))) /
              abs(as.numeric(logLik(f3))), 1e-8)
  # a body flat in its coefficients at the start but not at the moved
  # points passes: a * exp(b * x) at a = 0
  f4 <- suppressWarnings(frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE),
                             family = gaussian(), data = d,
                             start = list(beta = c(1, 0.1))))
  expect_s3_class(f4, "frmtmb_fit")
})

# ---------------------------------------------------------------------
# Default starts
# ---------------------------------------------------------------------

test_that("a body the default link cannot evaluate starts at the mean", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x, y = z$pos)
  body <- bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)
  pt <- par_template(body, data = d, family = Gamma())
  expect_equal(unname(pt[["beta"]]), c(1 / mean(d$y), 0))
  expect_message(frm(body, family = Gamma(), data = d),
                 "Nonlinear starting values moved")
  fit <- suppressMessages(frm(body, family = inverse.gaussian(), data = d))
  lin <- frm(bf(y ~ x), family = inverse.gaussian(), data = d)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(lin))) /
              abs(as.numeric(logLik(lin))), 1e-8)
  pt <- par_template(body, data = d, family = inverse.gaussian())
  expect_equal(unname(pt[["beta"]]), c(1 / mean(d$y)^2, 0))
})

test_that("a start the link can evaluate does not move", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x, y = z$pos)
  body <- bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)
  # a log link evaluates eta = 0
  pt <- par_template(body, data = d, family = Gamma(link = "log"))
  expect_identical(unname(pt[["beta"]]), c(0, 0))
  # the user's start and a located prior stand
  expect_no_message(f1 <- frm(body, family = Gamma(), data = d,
                              start = list(beta = c("b0_(Intercept)" = 0.4))))
  pr <- set_prior("normal(0.4, 1)", nlpar = "b0")
  pt2 <- par_template(body, data = d, family = Gamma(), prior = pr)
  expect_identical(unname(pt2[["beta"]])[1], 0.4)
})

test_that("a body undefined at its start keeps the error that asks for one", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x, y = z$pos)
  # log(a) at a = 0 is no value to solve from
  msg <- tryCatch(frm(bf(y ~ log(a) + b1 * x, a + b1 ~ 1, nl = TRUE),
                      family = Gamma(), data = d), error = conditionMessage)
  expect_match(msg, "failed from its default starting values", fixed = TRUE)
  # b1 * x changes sign, so no single b1 keeps an inverse link in its domain
  msg <- tryCatch(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = Gamma(),
                      data = d), error = conditionMessage)
  expect_match(msg, "failed from its default starting values", fixed = TRUE)
})

test_that("a multivariate model moves only the body that needs it", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x, y = z$pos, g2 = z$eta + z$e1)
  f <- mvbf(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, family = Gamma()),
            bf(g2 ~ c0 + c1 * x, c0 + c1 ~ 1, nl = TRUE, family = gaussian()))
  fit <- suppressMessages(frm(f, data = d))
  a1 <- as.numeric(logLik(frm(bf(y ~ x), family = Gamma(), data = d)))
  a2 <- as.numeric(logLik(frm(bf(g2 ~ x), family = gaussian(), data = d)))
  expect_lt(abs(as.numeric(logLik(fit)) - (a1 + a2)) / abs(a1 + a2), 1e-8)
})

# ---------------------------------------------------------------------
# Residual correlation terms with a body
# ---------------------------------------------------------------------

nl_ac_data <- function(seed = 158, n = 120, ngrp = 6) {
  set.seed(seed)
  d <- data.frame(x = round(rnorm(n), 3),
                  g = factor(sort(rep(seq_len(ngrp), length.out = n))))
  d$t <- stats::ave(seq_len(n), d$g, FUN = seq_along)
  e <- stats::ave(rnorm(n), d$g, FUN = function(v) {
    as.numeric(stats::filter(v, 0.4, method = "recursive"))
  })
  d$y <- 2 + 0.5 * d$x + e
  d
}

test_that("bf(nl = TRUE, autocor =) equals the linear model", {
  d <- nl_ac_data()
  st <- list(beta = c("b0_(Intercept)" = 2, "b1_(Intercept)" = 0.5))
  terms <- list(~ ar(time = t, gr = g, cov = TRUE),
                ~ ar(time = t, gr = g, cov = FALSE),
                ~ ma(time = t, gr = g, cov = TRUE),
                ~ ma(time = t, gr = g, cov = FALSE),
                ~ arma(time = t, gr = g, cov = TRUE),
                ~ arma(time = t, gr = g, cov = FALSE),
                ~ cosy(time = t, gr = g))
  for (tm in terms) {
    lin <- allow_warnings(frm(bf(y ~ x, autocor = tm), data = d),
                          allowed = c("Standard errors", "standard error"))
    fit <- allow_warnings(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE,
                                 autocor = tm), data = d, start = st),
                          allowed = c("Standard errors", "standard error"))
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(lin))) /
                abs(as.numeric(logLik(lin))), 1e-6)
    expect_identical(fit$spec$responses[[1]]$autocor$label,
                     lin$spec$responses[[1]]$autocor$label)
    # the mean a post-fit path reads is the linear model's, the cov = FALSE
    # shift included
    expect_lt(max(abs(frm_linpred(fit) - frm_linpred(lin))) /
                max(abs(frm_linpred(lin))), 1e-3)
  }
})

test_that("fcor and acformula() take a body too", {
  d <- nl_ac_data()
  n <- nrow(d)
  M <- 0.3^abs(outer(seq_len(n), seq_len(n), "-"))
  st <- list(beta = c("b0_(Intercept)" = 2, "b1_(Intercept)" = 0.5))
  lin <- frm(bf(y ~ x, autocor = ~ fcor(M)), data = d, data2 = list(M = M))
  fit <- frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, autocor = ~ fcor(M)),
             data = d, data2 = list(M = M), start = st)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(lin))) /
              abs(as.numeric(logLik(lin))), 1e-6)
  f2 <- bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE) +
    acformula(~ ar(time = t, gr = g, cov = TRUE))
  fit2 <- frm(f2, data = d, start = st)
  f3 <- bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE,
           autocor = ~ ar(time = t, gr = g, cov = TRUE))
  fit3 <- frm(f3, data = d, start = st)
  expect_identical(as.numeric(logLik(fit2)), as.numeric(logLik(fit3)))
  # update() keeps the term
  up <- update(fit3, newdata = d[-120, ])
  expect_identical(up$spec$responses[[1]]$autocor$label,
                   fit3$spec$responses[[1]]$autocor$label)
})

test_that("a term in the body itself stays refused", {
  d <- nl_ac_data()
  expect_error(frm(bf(y ~ b1 * x + ar(time = t, gr = g), b1 ~ 1, nl = TRUE),
                   data = d), "autocor = ~ ar", fixed = TRUE)
  # and brms's own refusals of the term hold on the body route
  expect_error(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE,
                      autocor = ~ b0 + ar(time = t, gr = g)), data = d),
               "belongs in the model formula", fixed = TRUE)
  expect_identical(frm_compat("ar()", "nl")$status, "conditional")
})

test_that("a mixture's location takes no term, with or without a body", {
  d <- nl_ac_data()
  expect_error(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE,
                      autocor = ~ ar(time = t, gr = g)),
                   family = mixture(gaussian(), gaussian()), data = d),
               "only be (written|specified) on 'mu'")
})

# ---------------------------------------------------------------------
# A body constant across the rows at a stationary point: the review of
# 2026-11-02, blocker 1. b1 * exp(b2 * x) at b1 = b2 = 0 under cumulative():
# every gradient of a body coefficient is exactly 0, because the thresholds
# start at their exact maximum, and frm() returned -413.5626 with code 0,
# the logLik of y ~ 1, where brms and a `start` reach -371.7584. Seen on
# an earlier build, which carried the design fix alone.
# ---------------------------------------------------------------------

nl_saddle_data <- function(seed = 9202, n = 300) {
  set.seed(seed)
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                     c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  d
}

test_that("a saddle at the default start is left, on every ordinal family", {
  d <- nl_saddle_data()
  body <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  st <- list(beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.1))
  for (fam in list(cumulative(), cumulative("probit"), sratio(), cratio(),
                   acat())) {
    only <- as.numeric(logLik(frm(y ~ 1, family = fam, data = d)))
    started <- as.numeric(logLik(allow_warnings(
      frm(body, family = fam, data = d, start = st),
      allowed = c("Standard errors", "did not report convergence",
                  "not identified"))))
    fit <- allow_warnings(suppressMessages(frm(body, family = fam, data = d)),
                          allowed = c("Standard errors", "not identified",
                                      "did not report convergence"))
    ll <- as.numeric(logLik(fit))
    # at least what a start in the region reaches: acat() from the default
    # start ran up the ridge to the linear limit, 5.9 nats short
    expect_gt(ll, started - 1e-6 * abs(started))
    expect_gt(ll - only, 10)
  }
  # the first one by name: it says what happened, and a warning that the
  # data do not determine the coefficients does not reach the user
  expect_message(fit <- frm(body, family = cumulative(), data = d),
                 "reached a higher point")
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(
    frm(body, family = cumulative(), data = d, start = st)))) /
      abs(as.numeric(logLik(fit))), 1e-6)
  expect_gt(abs(fit$estimates[["beta"]][["b2_(Intercept)"]]), 0.1)
})

test_that("hurdle_cumulative leaves the same saddle", {
  d <- nl_saddle_data()
  set.seed(1)
  d$y <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)),
                 levels = 0:4)
  fit <- suppressMessages(frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
                              family = hurdle_cumulative(), data = d))
  only <- as.numeric(logLik(frm(y ~ 1, family = hurdle_cumulative(), data = d)))
  expect_gt(as.numeric(logLik(fit)) - only, 10)
})

test_that("a saddle is left under other families and does not need the RNG", {
  set.seed(9203)
  n <- 250
  d <- data.frame(x = round(rnorm(n), 3))
  d$y <- 1 + 0.8 * d$x + rnorm(n, 0, 0.5)
  d$yc <- rpois(n, exp(0.5 + 0.5 * d$x))
  body <- bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE)
  # b1 and b2 are identified only through their product, and the fit says so
  run <- function() {
    allow_warnings(suppressMessages(frm(body, family = gaussian(), data = d)),
                   allowed = c("Standard errors", "not identified"))
  }
  set.seed(1)
  f1 <- run()
  set.seed(2)
  f2 <- run()
  expect_identical(f1$estimates[["beta"]], f2$estimates[["beta"]])
  lin <- as.numeric(logLik(frm(bf(y ~ x), family = gaussian(), data = d)))
  expect_lt(abs(as.numeric(logLik(f1)) - lin) / abs(lin), 1e-6)
  fp <- allow_warnings(suppressMessages(
    frm(bf(yc ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE),
        family = poisson(), data = d)),
    allowed = c("Standard errors", "not identified"))
  lp <- as.numeric(logLik(frm(bf(yc ~ x), family = poisson(), data = d)))
  expect_lt(abs(as.numeric(logLik(fp)) - lp) / abs(lp), 1e-6)
})

test_that("a fit that is not on a saddle comes back as it was", {
  d <- nl_saddle_data()
  # the restart never runs: a body that varies at the stop
  fit <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(),
             data = d)
  expect_length(nl_saddle_find(fit)$saddle, 0L)
  expect_length(nl_saddle_find(fit)$local, 0L)
  # a true maximum where the body is constant: x^2 cannot take the negative
  # slope of these data, and the curvature says the point is a maximum
  set.seed(9204)
  dn <- data.frame(x = round(rnorm(200), 3))
  dn$y <- 1 - 0.8 * dn$x + rnorm(200, 0, 0.5)
  fn <- allow_warnings(suppressMessages(
    frm(bf(y ~ a + b^2 * x, a + b ~ 1, nl = TRUE), family = gaussian(),
        data = dn, start = list(beta = c("a_(Intercept)" = 1,
                                         "b_(Intercept)" = 0)))),
    allowed = c("Standard errors", "not identified"))
  expect_length(nl_saddle_find(fn)$saddle, 0L)
  expect_length(nl_saddle_find(fn)$local, 1L)
  # a body that cannot vary at all: `a` is the model
  fa <- frm(bf(y ~ a, a ~ 1, nl = TRUE), family = gaussian(), data = dn)
  expect_length(nl_saddle_find(fa)$saddle, 0L)
  expect_length(nl_saddle_find(fa)$local, 0L)
  # a constant column gives a body with no covariate to move on
  dn$k <- 2
  fk <- allow_warnings(frm(bf(y ~ a * k, a ~ 1, nl = TRUE),
                           family = gaussian(), data = dn),
                       allowed = "Standard errors")
  expect_length(nl_saddle_find(fk)$saddle, 0L)
  expect_length(nl_saddle_find(fk)$local, 0L)
})

test_that("a point no restart leaves is a warning that names the body", {
  # every restart is mocked to fail, so the fit still sits on the point
  d <- nl_saddle_data()
  body <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  local_mocked_bindings(nl_saddle_starts = function(...) list())
  got <- list()
  fit <- withCallingHandlers(
    frm(body, family = cumulative(), data = d),
    warning = function(w) {
      got[[length(got) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(fit, "frmtmb_fit")
  named <- Filter(function(w) inherits(w, "frmtmb_nl_saddle"), got)
  expect_length(named, 1L)
  msg <- conditionMessage(named[[1L]])
  expect_match(msg, "b1 * exp(b2 * x)", fixed = TRUE)
  expect_match(msg, "thresholds-only model", fixed = TRUE)
  expect_match(msg, "no effect", fixed = TRUE)
  expect_match(msg, "`start`", fixed = TRUE)
  # with a start of the user's: the same warning
  w <- character(0)
  fit <- withCallingHandlers(
    frm(body, family = cumulative(), data = d,
        start = list(beta = c("b1_(Intercept)" = 0, "b2_(Intercept)" = 0))),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("thresholds-only model", w, fixed = TRUE)))
  expect_s3_class(fit, "frmtmb_fit")
})

# ---------------------------------------------------------------------
# Fix round 2 of the review of 2026-11-02: every body is moved, nothing stops,
# and what the fit says is true.
# ---------------------------------------------------------------------

nl_two_data <- function() {
  set.seed(11001)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$yn <- 1 + 0.7 * d$x + 0.6 * d$z + rnorm(n, 0, exp(-0.5 + 0.4 * d$z))
  invisible(runif(1))
  d$ymix <- ifelse(runif(n) < 0.4, -2 + 0.8 * d$x, 2 - 0.6 * d$x) +
    rnorm(n, 0, 0.5)
  d
}

test_that("a second body does not stay on its saddle", {
  d <- nl_two_data()
  fm <- bf(yn ~ 1) + nlf(mu ~ a + b1 * c1 * x) +
    nlf(sigma ~ exp(s0 + s1 * s2 * z)) + lf(a ~ 1) + lf(b1 ~ 1) + lf(c1 ~ 1) +
    lf(s0 ~ 1) + lf(s1 ~ 1) + lf(s2 ~ 1)
  allowed <- c("Standard errors", "not identified",
               "did not report convergence", "gradient")
  fit <- allow_warnings(suppressMessages(frm(fm, data = d)), allowed = allowed)
  st <- allow_warnings(frm(fm, data = d, start = list(
    beta = c("a_(Intercept)" = 1, "b1_(Intercept)" = 1, "c1_(Intercept)" = 0.7),
    betad = c("s0_(Intercept)" = -0.5, "s1_(Intercept)" = 0.5,
              "s2_(Intercept)" = 0.5))), allowed = allowed)
  ref <- as.numeric(logLik(st))
  # what a start in the region reaches: the sigma body read 16.8 nats short with
  # its slopes at 0
  expect_gt(as.numeric(logLik(fit)), ref - 1e-6 * abs(ref))
  expect_gt(diff(range(frm_linpred(fit, dpar = "sigma"))), 0.1)

  fm2 <- bf(ymix ~ 1) + nlf(mu1 ~ a1 + b1 * c1 * x) +
    nlf(mu2 ~ a2 + b2 * c2 * x) + lf(a1 ~ 1) + lf(b1 ~ 1) + lf(c1 ~ 1) +
    lf(a2 ~ 1) + lf(b2 ~ 1) + lf(c2 ~ 1)
  fam <- mixture(gaussian(), gaussian())
  fit2 <- allow_warnings(suppressMessages(frm(fm2, family = fam, data = d)),
                         allowed = allowed)
  lin2 <- as.numeric(logLik(allow_warnings(
    frm(bf(ymix ~ 1, mu1 ~ x, mu2 ~ x), family = fam, data = d),
    allowed = allowed)))
  # the mu2 body read 65.7 nats short
  expect_gt(as.numeric(logLik(fit2)), lin2 - 1e-6 * abs(lin2))
  expect_gt(diff(range(frm_linpred(fit2, dpar = "mu2"))), 0.1)
})

test_that("two ordinal bodies leave the point without luck from a joint start", {
  set.seed(11001)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  cutp <- c(-Inf, 0.5, 1.5, 3, Inf)
  d$y1 <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n), cutp, labels = 1:4))
  d$y2 <- ordered(cut(1.2 * exp(-0.5 * d$z) + rlogis(n), cutp, labels = 1:4))
  f1 <- bf(y1 ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  f2 <- bf(y2 ~ c1 * exp(c2 * z), c1 + c2 ~ 1, nl = TRUE)
  fm <- mvbf(f1, f2, rescor = FALSE)
  allowed <- c("Standard errors", "not identified",
               "did not report convergence", "gradient")
  sum1 <- as.numeric(logLik(allow_warnings(suppressMessages(
    frm(f1, family = cumulative(), data = d)), allowed = allowed))) +
    as.numeric(logLik(allow_warnings(suppressMessages(
      frm(f2, family = cumulative(), data = d)), allowed = allowed)))
  ctl <- frmtmb_control(optCtrl = list(iter.max = 0))
  # a first fit with both bodies on the point, whatever the joint start reaches
  at0 <- allow_warnings(frm(fm, family = cumulative(), data = d, control = ctl),
                        allowed = allowed)
  look <- nl_saddle_find(at0, use_curvature = FALSE)
  expect_length(look$saddle, 2L)
  # one start moves the coefficients of both
  st <- nl_saddle_starts(at0$frame, look$saddle, NULL)[[1L]]
  expect_setequal(names(st$beta), c("y1_b1_(Intercept)", "y1_b2_(Intercept)",
                                    "y2_c1_(Intercept)", "y2_c2_(Intercept)"))
  # flips of the second body reach each of its signs, the first body held
  flips <- nl_saddle_starts(at0$frame, look$saddle, NULL, 2L, at0)
  expect_true(any(vapply(flips, function(s) {
    s$beta[["y2_c2_(Intercept)"]] < 0
  }, NA)))
  expect_identical(flips[[1L]]$beta[["y1_b1_(Intercept)"]],
                   at0$estimates$beta[["y1_b1_(Intercept)"]])
  run <- function(s, announce) {
    if (is.null(s)) return(at0)
    frm(fm, family = cumulative(), data = d, start = s)
  }
  fit <- suppressMessages(suppressWarnings(
    nl_fit_off_saddle(run, at0$frame, NULL, NULL)))
  expect_lt(abs(as.numeric(logLik(fit)) - sum1) / abs(sum1), 1e-6)
})

test_that("a body with no effect in the data gets a warning, not a stop", {
  set.seed(11002)
  n <- 200
  d <- data.frame(x = round(rnorm(n), 3))
  d$x <- d$x - mean(d$x)
  e <- rnorm(n)
  # a residual exactly orthogonal to x: the slope is 0 to rounding
  d$y <- 1 + stats::resid(stats::lm(e ~ x, data = d))
  body <- bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE)
  got <- list()
  fit <- withCallingHandlers(
    suppressMessages(frm(body, data = d)),
    warning = function(w) {
      got[[length(got) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(fit, "frmtmb_fit")
  lin <- as.numeric(logLik(frm(y ~ x, data = d)))
  expect_lt(abs(as.numeric(logLik(fit)) - lin) / abs(lin), 1e-6)
  msgs <- vapply(Filter(function(w) inherits(w, "frmtmb_nl_saddle"), got),
                 conditionMessage, "")
  expect_length(msgs, 1L)
  expect_match(msgs, "a + b1 * b2 * x", fixed = TRUE)
  expect_match(msgs, "no effect", fixed = TRUE)
  expect_match(msgs, "`start`", fixed = TRUE)
})

test_that("a sharp true maximum is kept, however little a step leaves it", {
  # sigma 1e-6: a step of 1e-7 moves the objective by more than 1e-3 of its
  # value there, and a smoothness test on the restart refused the true point
  set.seed(11004)
  n <- 200
  d <- data.frame(x = round(rnorm(n), 3))
  d$y <- 1 + 2 * d$x + rnorm(n, 0, 1e-6)
  body <- bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE)
  fit <- allow_warnings(suppressMessages(frm(body, data = d)),
                        allowed = c("Standard errors", "not identified",
                                    "did not report convergence", "gradient",
                                    "constant across the rows"))
  expect_s3_class(fit, "frmtmb_fit")
  lin <- as.numeric(logLik(frm(y ~ x, data = d)))
  # the same fit from a start at the truth
  st <- allow_warnings(frm(body, data = d, start = list(beta = c(
    "a_(Intercept)" = 1, "b1_(Intercept)" = 1, "b2_(Intercept)" = 2))),
    allowed = c("Standard errors", "not identified",
                "did not report convergence", "gradient"))
  expect_lt(abs(as.numeric(logLik(st)) - lin) / abs(lin), 1e-6)
  # a default start need not reach it, but it warns if it does not
  if (abs(as.numeric(logLik(fit)) - lin) / abs(lin) > 1e-6) {
    expect_warning(suppressMessages(frm(body, data = d)),
                   class = "frmtmb_nl_saddle")
  }
})

test_that("a constant body at a local maximum is said in a message", {
  set.seed(11005)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3))
  d$yo <- ordered(cut(-0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf),
                      labels = 1:4))
  body <- bf(yo ~ b^2 * x * (1 - b^2), b ~ 1, nl = TRUE)
  msgs <- list()
  fit <- withCallingHandlers(
    frm(body, family = cumulative(), data = d),
    message = function(m) {
      msgs[[length(msgs) + 1L]] <<- m
      invokeRestart("muffleMessage")
    })
  only <- as.numeric(logLik(frm(yo ~ 1, family = cumulative(), data = d)))
  expect_lt(abs(as.numeric(logLik(fit)) - only) / abs(only), 1e-8)
  named <- Filter(function(m) inherits(m, "frmtmb_nl_local"), msgs)
  expect_length(named, 1L)
  expect_match(conditionMessage(named[[1L]]), "local maximum", fixed = TRUE)
  expect_match(conditionMessage(named[[1L]]), "thresholds-only", fixed = TRUE)
})

test_that("iter.max = 0 reads the objective at the start", {
  d <- nl_saddle_data()
  body <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  ctl <- frmtmb_control(optCtrl = list(iter.max = 0))
  msgs <- character(0)
  fit <- withCallingHandlers(
    allow_warnings(frm(body, family = cumulative(), data = d, control = ctl),
                   allowed = c("Standard errors", "not identified",
                               "did not report convergence", "gradient")),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_false(any(grepl("restarted", msgs, fixed = TRUE)))
  expect_identical(unname(fit$estimates$beta[c("b1_(Intercept)",
                                               "b2_(Intercept)")]), c(0, 0))
})

test_that("the return code alone does not start a restart", {
  d <- nl_saddle_data()
  # a capped fit away from a saddle: the code is the cap's
  body <- bf(y ~ b1 * x, b1 ~ 1, nl = TRUE)
  n_run <- 0L
  full <- frm(body, family = cumulative(), data = d)
  capped <- full
  capped$opt$convergence <- 1L
  run <- function(st, announce) {
    n_run <<- n_run + 1L
    capped
  }
  fit <- nl_fit_off_saddle(run, full$frame, NULL, NULL)
  expect_identical(n_run, 1L)
  expect_identical(fit$opt$convergence, 1L)
})

test_that("a restart past an agreement finds the flip of the sign", {
  # acat, seed 4, slope -0.6: a stop at the first agreement read -376.2397 at b1
  # -15200, the ridge, at code 0, and the patterns past it read -370.2536
  # (dev/reviews/2026-11-02-nlfix.md, release fix 2; dev/review-nlfix-d7.R)
  set.seed(401)
  n <- 300
  d <- data.frame(x = round(stats::rnorm(n), 3), z = round(stats::rnorm(n), 3))
  lat <- 1.5 * exp(-0.6 * d$x) + 0.3 * d$z
  d$y <- ordered(cut(lat + stats::rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  body <- bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1, nl = TRUE)
  fit <- suppressMessages(suppressWarnings(frm(body, family = acat(), data = d)))
  best <- -Inf
  for (b2 in c(-0.5, 0.5)) for (b1 in c(-1, 1)) {
    r <- suppressMessages(suppressWarnings(frm(
      body, family = acat(), data = d,
      start = list(beta = c(b1_Intercept = b1, b2_Intercept = b2)))))
    best <- max(best, as.numeric(logLik(r)))
  }
  # 6 nats separate the ridge from the flip, on a scale the fit reads to 1e-6
  expect_lt(best - as.numeric(logLik(fit)), 1e-5 * abs(best))
})

test_that("two restarts on one ridge do not end the round", {
  # equal values at different estimates: the sign of b2 comes from the third start
  set.seed(9202)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                     c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  body <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
  allowed <- c("Standard errors", "not identified",
               "did not report convergence", "gradient")
  for (fam in list(cratio(), sratio())) {
    fit <- allow_warnings(suppressMessages(frm(body, family = fam, data = d)),
                          allowed = allowed)
    st <- allow_warnings(suppressMessages(
      frm(body, family = fam, data = d, start = list(
        beta = c("b1_(Intercept)" = 0.16, "b2_(Intercept)" = -0.18)))),
      allowed = allowed)
    ref <- as.numeric(logLik(st))
    expect_gt(as.numeric(logLik(fit)), ref - 1e-6 * abs(ref))
  }
})

test_that("the restart's message says what is known", {
  # two values that print the same at 7 digits are told apart
  ll <- nl_ll_pair(-541.23671, -541.23674)
  expect_false(identical(ll[1L], ll[2L]))
  expect_identical(nl_ll_pair(-413.5, -371.2), c("-413.5", "-371.2"))
  d <- nl_saddle_data()
  body <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  msgs <- character(0)
  fit <- withCallingHandlers(
    frm(body, family = cumulative(), data = d),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  m <- grep("reached a higher point", msgs, value = TRUE, fixed = TRUE)
  expect_length(m, 1L)
  expect_match(m, "return code 0", fixed = TRUE)
  expect_no_match(m, "better maximum", fixed = TRUE)
  expect_no_match(m, "thresholds-only", fixed = TRUE)
  # a disc body beside a mu body that varies is not the thresholds-only model
  set.seed(11001)
  n <- 300
  dd <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  dd$yo <- ordered(cut(1.5 * exp(0.6 * dd$x) + 0.7 * dd$z + rlogis(n),
                       c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  fm <- bf(yo ~ 1) + nlf(mu ~ b1 * exp(b2 * x)) +
    nlf(disc ~ exp(d1 * d2 * z)) + lf(b1 ~ 1) + lf(b2 ~ 1) + lf(d1 ~ 1) +
    lf(d2 ~ 1)
  ctl <- frmtmb_control(optCtrl = list(iter.max = 0))
  at0 <- allow_warnings(frm(fm, family = cumulative(), data = dd, control = ctl),
                        allowed = c("Standard errors", "not identified",
                                    "did not report convergence", "gradient"))
  frame <- at0$frame
  mu <- list(lp = frame$linpreds[[linpred_key("yo", "mu")]],
             seen = c("b1", "b2"))
  disc <- list(lp = frame$linpreds[[linpred_key("yo", "disc")]],
               seen = c("d1", "d2"))
  plain <- "the model without those covariates"
  expect_identical(nl_model_words(frame, list(disc)), plain)
  # nor is a mu body with a disc body that still varies
  expect_identical(nl_model_words(frame, list(mu)), plain)
  expect_identical(nl_model_words(frame, list(mu, disc)),
                   "the thresholds-only model")
})

test_that("a disc body does not break the ordinal fit-end check", {
  set.seed(11001)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$yo <- ordered(cut(0.8 * d$x + rlogis(n, 0, exp(0.3 * d$z)),
                      c(-Inf, -1, 0, 1, Inf), labels = 1:4))
  got <- character(0)
  # a body with no constant term of its own: identified, and no warning about it
  fit <- withCallingHandlers(
    frm(bf(yo ~ 0 + x) + nlf(disc ~ d1 * z) + lf(d1 ~ 1),
        family = cumulative(), data = d),
    warning = function(w) {
      got <<- c(got, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_s3_class(fit, "frmtmb_fit")
  expect_no_match(got, "fit_check", fixed = TRUE)
  expect_identical(ord_fit_check(fit, "yo"), NULL)
  # the same scale problem as `disc ~ 1 + z`, seen by the standard errors, which
  # is what the help says of it
  got <- character(0)
  fit <- withCallingHandlers(
    frm(bf(yo ~ 0 + x) + nlf(disc ~ d0 + d1 * z) + lf(d0 ~ 1) + lf(d1 ~ 1),
        family = cumulative(), data = d),
    warning = function(w) {
      got <<- c(got, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_no_match(got, "fit_check", fixed = TRUE)
  expect_true(any(grepl("Standard errors|not identified", got)))
})

# ---------------------------------------------------------------------
# The guard's false refusals: the review of 2026-11-02, major 1. Each of
# these bodies stopped with "does not reach the likelihood" on the first
# an earlier build.
# ---------------------------------------------------------------------

test_that("the guard leaves a body whose value did not move", {
  set.seed(6201)
  n <- 200
  d <- data.frame(x = rnorm(n), x2 = runif(n, 10, 20))
  d$y <- 0.5 * d$x + rnorm(n)
  # the move is not a null direction of a linear combination
  for (f in list(bf(y ~ (2 * a - b) * x, a + b ~ 1, nl = TRUE),
                 bf(y ~ 2 * a - b, a + b ~ 1, nl = TRUE),
                 bf(y ~ a - 2 * b + c, a + b + c ~ 1, nl = TRUE))) {
    expect_s3_class(allow_warnings(frm(f, data = d),
                                   allowed = c("Standard errors",
                                               "not identified")),
                    "frmtmb_fit")
  }
  expect_s3_class(
    frm(bf(y ~ 2 * a - b, a + b ~ 1, nl = TRUE), data = d,
        prior = c(set_prior("normal(0, 1)", nlpar = "a"),
                  set_prior("normal(0, 1)", nlpar = "b"))),
    "frmtmb_fit")
  # below rounding of the objective
  d$xt <- d$x * 1e-17
  expect_s3_class(allow_warnings(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE),
                                     data = d),
                                 allowed = c("Standard errors", "gradient",
                                             "did not report")), "frmtmb_fit")
  d$yb <- 1e6 + 0.5e6 * d$x
  d$xs <- d$x * 1e-10
  expect_s3_class(allow_warnings(frm(bf(yb ~ b * xs, b ~ 1, nl = TRUE),
                                     data = d),
                                 allowed = c("Standard errors", "gradient",
                                             "did not report")), "frmtmb_fit")
  # flat or saturated around the start
  d$ys <- ifelse(d$x2 > 15, 1, 0) + rnorm(n, 0, 0.3)
  expect_s3_class(allow_warnings(frm(bf(ys ~ step(x2 - c), c ~ 1, nl = TRUE),
                                     data = d),
                                 allowed = c("Standard errors",
                                             "not identified")), "frmtmb_fit")
  d$xh <- d$x2 * 1000
  d$yt <- exp(-0.002 * d$x2 * 100) + rnorm(n, 0, 0.1)
  expect_s3_class(allow_warnings(
    frm(bf(yt ~ exp(-k * xh), k ~ 1, nl = TRUE), data = d,
        start = list(beta = c("k_(Intercept)" = 1))),
    allowed = c("Standard errors", "gradient", "did not report",
                "not identified")), "frmtmb_fit")
})

test_that("the guard's true positives still fire", {
  local_frmtmb_registries()
  set.seed(6202)
  d <- data.frame(x = rnorm(100))
  d$y <- 0.5 * d$x + rnorm(100)
  body <- bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)
  times0 <- frmtmb_family("tfamzero", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      RTMB::dnorm(y, dpars$mu * 0, dpars$sigma, log = TRUE)
    })
  expect_error(frm(body, family = times0, data = d),
               "moves the body and leaves the objective", fixed = TRUE)
  # a density that reads two rows of mu is not refused
  two <- frmtmb_family("tfamtwo", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      m <- dpars$mu * as.numeric(y > stats::quantile(as.numeric(y), 0.98))
      RTMB::dnorm(y, m, dpars$sigma, log = TRUE)
    })
  expect_s3_class(allow_warnings(frm(body, family = two, data = d),
                                 allowed = c("Standard errors",
                                             "not identified")), "frmtmb_fit")
  # softmax of identical locations is constant under refcat = NA
  dc <- d
  dc$cat <- factor(sample(c("a", "b", "c"), 100, TRUE))
  expect_error(frm(bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                   family = categorical(refcat = NA), data = dc),
               "does not reach the likelihood", fixed = TRUE)
})

# ---------------------------------------------------------------------
# Blocker 2: simulate(newdata =) on a body under cov = FALSE, which read the
# one-step mean and then ran its own recursion on top. Slope of the row means
# on the observed shift, measured 1.0032 against 0.0032 on the linear spelling.
# ---------------------------------------------------------------------

test_that("simulate(newdata =) of a cov = FALSE body has an empty past", {
  d <- nl_ac_data(seed = 4401)
  term <- ~ ar(time = t, gr = g, cov = FALSE)
  st <- list(beta = c("b0_(Intercept)" = 2, "b1_(Intercept)" = 0.5))
  fl <- frm(bf(y ~ x, autocor = term), data = d)
  fn <- frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, autocor = term),
            data = d, start = st)
  nd <- d[d$g %in% levels(d$g)[1:2], ]
  N <- 2000
  sl <- as.matrix(simulate(fl, nsim = N, seed = 3, newdata = nd))
  sn <- as.matrix(simulate(fn, nsim = N, seed = 3, newdata = nd))
  mc <- max(apply(sl, 1, stats::sd)) / sqrt(N)
  # the same seed and a body equal to the linear predictor, up to the
  # optimizer's stop; a double-counted AR part is 25 Monte Carlo errors away
  expect_lt(max(abs(rowMeans(sn) - rowMeans(sl))), 5 * mc)
})

# ---------------------------------------------------------------------
# Major 2: conditional_effects() on an ordinal body.
# ---------------------------------------------------------------------

test_that("conditional_effects() of an ordinal body has the linear's band", {
  z <- nl_fam_data()
  d <- data.frame(x = z$x)
  d$y <- ordered(z$ord)
  fl <- frm(bf(y ~ x), family = cumulative(), data = d)
  fn <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d)
  cl <- conditional_effects(fl, effects = "x")[[1L]]
  cn <- conditional_effects(fn, effects = "x")[[1L]]
  expect_lt(max(abs(cl$estimate__ - cn$estimate__)) /
              diff(range(cl$estimate__)), 1e-4)
  expect_lt(max(abs(cl$se__ - cn$se__)) / max(cl$se__), 1e-3)
  kl <- conditional_effects(fl, effects = "x", categorical = TRUE)[[1L]]
  kn <- conditional_effects(fn, effects = "x", categorical = TRUE)[[1L]]
  expect_lt(max(abs(kl$se__ - kn$se__)) / max(kl$se__), 1e-3)
  # and the bootstrap, the remedy on the other displays
  set.seed(7)
  bn <- suppressWarnings(conditional_effects(fn, effects = "x", band = "boot",
                                             boot = 5, resolution = 5))[[1L]]
  expect_true(all(is.finite(bn$estimate__)))
  bk <- suppressWarnings(conditional_effects(fn, effects = "x", band = "boot",
                                             categorical = TRUE, boot = 5,
                                             resolution = 5))[[1L]]
  expect_true(all(is.finite(bk$estimate__)))
  expect_identical(nlevels(bk$cats__), nlevels(d$y))
})

test_that("the refusal on other displays names what the family takes", {
  set.seed(7302)
  d <- data.frame(x = rnorm(150))
  d$yz <- ifelse(runif(150) < 0.3, 0L, rpois(150, exp(0.5 + 0.5 * d$x)))
  fz <- frm(bf(yz ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
            family = zero_inflated_poisson(), data = d)
  msg <- tryCatch(conditional_effects(fz, effects = "x"),
                  error = conditionMessage)
  expect_match(msg, "band = \"boot\"", fixed = TRUE)
  expect_match(msg, "method = \"predict\"", fixed = TRUE)
  expect_match(msg, "dpar = \"mu\"", fixed = TRUE)
  # each of the three works
  set.seed(7)
  for (a in list(list(band = "boot", boot = 5, resolution = 5),
                 list(method = "predict", resolution = 5),
                 list(dpar = "mu", resolution = 5))) {
    r <- suppressWarnings(do.call(conditional_effects,
                                  c(list(fz, effects = "x"), a)))
    expect_true(all(is.finite(r[[1L]]$estimate__)))
  }
  # a category display of a body with several locations reads fitted()
  set.seed(7301)
  d$cat <- factor(sample(c("a", "b", "c"), 150, TRUE, prob = c(1, 2, 3)))
  fc <- frm(bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
            family = categorical(), data = d, start = list(beta = c(0.1, 0.3)))
  ce <- conditional_effects(fc, effects = "x", categorical = TRUE)[[1L]]
  expect_true(all(is.finite(ce$se__)))
})

# ---------------------------------------------------------------------
# Minors 6, 7, 8
# ---------------------------------------------------------------------

test_that("nlf(mu ~ ...) takes a residual term by either route", {
  d <- nl_ac_data()
  st <- list(beta = c("b0_(Intercept)" = 2, "b1_(Intercept)" = 0.5))
  term <- ~ ar(time = t, gr = g, cov = TRUE)
  ref <- frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, autocor = term),
             data = d, start = st)
  f1 <- frm(bf(y ~ 1, autocor = term) + nlf(mu ~ b0 + b1 * x) + lf(b0 ~ 1) +
              lf(b1 ~ 1), data = d, start = st)
  f2 <- frm(bf(y ~ 1) + nlf(mu ~ b0 + b1 * x) + lf(b0 ~ 1) + lf(b1 ~ 1) +
              acformula(term), data = d, start = st)
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(ref)))
  expect_identical(as.numeric(logLik(f2)), as.numeric(logLik(ref)))
  # the right-hand side of that bf() still has nothing else to say
  expect_error(frm(bf(y ~ x, autocor = term) + nlf(mu ~ b0 + b1 * x) +
                     lf(b0 ~ 1) + lf(b1 ~ 1), data = d, start = st),
               "with nothing to do", fixed = TRUE)
})

test_that("a helper called cs() is the user's", {
  set.seed(4401)
  d <- data.frame(x = rnorm(100))
  d$y <- 0.5 * d$x + rnorm(100)
  cs <- function(v) v * 2
  f <- stats::as.formula("y ~ a * cs(x)")
  environment(f) <- environment()
  fit <- frm(bf(f, a ~ 1, nl = TRUE), data = d)
  expect_equal(2 * fit$estimates[["beta"]][[1L]],
               unname(stats::coef(stats::lm(y ~ 0 + x, d))), tolerance = 1e-5)
})

test_that("update(formula. =) carries the term through the stored call", {
  # a dropped term here is silent: the two logLiks differ by 6.6 nats
  d <- nl_ac_data()
  d$z <- rev(d$x)
  st <- list(beta = c("b0_(Intercept)" = 2, "b1_(Intercept)" = 0.5))
  term <- ~ ar(time = t, gr = g, cov = TRUE)
  fit <- frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, autocor = term),
             data = d, start = st)
  up <- suppressMessages(update(fit, formula. = bf(y ~ b0 + b1 * z, nl = TRUE,
                                                   autocor = term), start = st))
  ref <- frm(bf(y ~ b0 + b1 * z, b0 + b1 ~ 1, nl = TRUE, autocor = term),
             data = d, start = st)
  expect_identical(length(up$frame[["autocor"]]), 1L)
  expect_lt(abs(as.numeric(logLik(up)) - as.numeric(logLik(ref))) /
              abs(as.numeric(logLik(ref))), 1e-7)
  up2 <- suppressMessages(update(up, data = d[d$g != levels(d$g)[1L], ],
                                 start = st))
  expect_identical(length(up2$frame[["autocor"]]), 1L)
  # a stored formula with the term comes back from the call
  expect_true("autocor" %in% names(as.list(up$call[[2L]])) ||
                !is.null(up$bform[["autocor"]]))
})

test_that("a zero row of x^b has Stan's value and gradient", {
  set.seed(9302)
  n <- 150
  d <- data.frame(x = c(0, 0, runif(n - 2, 0, 3)))
  d$y <- 2 * d$x^0.5 + rnorm(n, 0, 0.3)
  for (f in list(bf(y ~ a * x^b, a + b ~ 1, nl = TRUE),
                 bf(y ~ a * pow(x, b), a + b ~ 1, nl = TRUE))) {
    fit <- frm(f, data = d, start = list(beta = c(1, 0.5)))
    expect_identical(fit$opt$convergence, 0L)
    # the zero rows do not depend on b, so the fit on the other rows with
    # those two rows' mean fixed at 0 has the same b
    b <- fit$estimates[["beta"]]
    ref <- stats::nls(y ~ a * x^b, data = d[d$x > 0, ],
                      start = list(a = 1, b = 0.5))
    expect_lt(abs(b[[2L]] - stats::coef(ref)[["b"]]) / abs(b[[2L]]), 0.05)
  }
  # value at the zero row: 1 at y = 0, 0 above, Inf below, on the tape too
  tp <- RTMB::MakeTape(function(p) nl_pow(c(0, 2), p), 0.5)
  expect_equal(tp(0.5), c(0, sqrt(2)))
  expect_equal(tp(0), c(1, 1))
  expect_identical(tp(-1)[1L], Inf)
  expect_equal(as.numeric(tp$jacobian(0.5)), c(0, sqrt(2) * log(2)))
})
