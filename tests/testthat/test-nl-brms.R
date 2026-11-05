# Log-density identity against brms for a nonlinear body: the checks of
# test-brms-likelihood.R (brms_lp_check() in helper-brms.R) on the shapes
# that frm() dropped, refused or could not start through 0.78.0.
#
# Stan compiles here. The file is opt-in, as test-brms-likelihood.R is:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")

nlb_ord <- function(seed, hurdle = FALSE, n = 200) {
  set.seed(seed)
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  g = factor(sprintf("g%02d", sort(rep(1:10, length.out = n)))))
  eta <- 0.3 + 0.5 * d$x - 0.3 * d$z
  d$y <- ordered(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                     labels = 1:4))
  if (hurdle) {
    d$y <- ordered(ifelse(runif(n) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
  }
  d
}

nlb_ac <- function(seed, n = 200, ngrp = 10) {
  set.seed(seed)
  d <- data.frame(x = round(rnorm(n), 3),
                  g = factor(sprintf("g%02d", sort(rep(seq_len(ngrp),
                                                       length.out = n)))))
  d$t <- stats::ave(seq_len(n), d$g, FUN = seq_along)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(n, 0, 0.5)
  d
}

test_that("an nl body under each ordinal family is brms's model", {
  skip_unless_brms_fit()
  fams <- list(cumulative = cumulative, sratio = sratio, cratio = cratio,
               acat = acat)
  bfams <- list(cumulative = brms::cumulative, sratio = brms::sratio,
                cratio = brms::cratio, acat = brms::acat)
  for (nm in names(fams)) {
    d <- nlb_ord(4900 + match(nm, names(fams)))
    fit <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = fams[[nm]](),
               data = d)
    brms_lp_check(brms::bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), bfams[[nm]](), d,
                  fit)
    fit <- frm(bf(y | thres(3) ~ b1 * x, b1 ~ 1, nl = TRUE),
               family = fams[[nm]](), data = d)
    brms_lp_check(brms::bf(y | thres(3) ~ b1 * x, b1 ~ 1, nl = TRUE),
                  bfams[[nm]](), d, fit)
    fit <- frm(bf(y ~ b1 * x, b1 ~ 1, disc ~ 0 + z, nl = TRUE),
               family = fams[[nm]](), data = d)
    brms_lp_check(brms::bf(y ~ b1 * x, b1 ~ 1, disc ~ 0 + z, nl = TRUE),
                  bfams[[nm]](), d, fit)
  }
  d <- nlb_ord(4951, hurdle = TRUE)
  fit <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = hurdle_cumulative(),
             data = d)
  brms_lp_check(brms::bf(y ~ b1 * x, b1 ~ 1, nl = TRUE),
                brms::hurdle_cumulative(), d, fit)
})

test_that("a body from the default start reaches brms's model off the saddle", {
  skip_unless_brms_fit()
  set.seed(9202)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3))
  d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                     c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  for (nm in c("cumulative", "sratio", "cratio", "acat")) {
    fit <- suppressWarnings(suppressMessages(
      frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
          family = get(nm, envir = asNamespace("frmtmb"))(), data = d)))
    only <- as.numeric(logLik(frm(y ~ 1,
      family = get(nm, envir = asNamespace("frmtmb"))(), data = d)))
    expect_gt(as.numeric(logLik(fit)) - only, 10)
    brms_lp_check(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
                  get(nm, envir = asNamespace("brms"))(), d, fit)
  }
})

test_that("a body with a zero row of x^b is brms's model", {
  skip_unless_brms_fit()
  set.seed(9302)
  n <- 150
  d <- data.frame(x = c(0, 0, runif(n - 2, 0, 3)))
  d$y <- 2 * d$x^0.5 + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ a * x^b, a + b ~ 1, nl = TRUE), data = d,
             start = list(beta = c(1, 0.5)))
  brms_lp_check(brms::bf(y ~ a * x^b, a + b ~ 1, nl = TRUE), gaussian(), d, fit)
  fit <- frm(bf(y ~ a * pow(x, b), a + b ~ 1, nl = TRUE), data = d,
             start = list(beta = c(1, 0.5)))
  brms_lp_check(brms::bf(y ~ a * pow(x, b), a + b ~ 1, nl = TRUE), gaussian(),
                d, fit)
})

test_that("an nl body under the other threshold structures", {
  skip_unless_brms_fit()
  d <- nlb_ord(4954)
  fit <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE),
             family = cumulative(threshold = "equidistant"), data = d)
  brms_lp_check(brms::bf(y ~ b1 * x, b1 ~ 1, nl = TRUE),
                brms::cumulative(threshold = "equidistant"), d, fit)
  d <- nlb_ord(4956)
  fit <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE),
             family = cumulative(threshold = "sum_to_zero"), data = d)
  brms_lp_check(brms::bf(y ~ b1 * x, b1 ~ 1, nl = TRUE),
                brms::cumulative(threshold = "sum_to_zero"), d, fit)
  d <- nlb_ord(4955)
  d$gg <- factor(ifelse(as.integer(d$g) <= 5, "p", "q"))
  fit <- frm(bf(y | thres(gr = gg) ~ b1 * x, b1 ~ 1, nl = TRUE),
             family = cumulative(), data = d)
  brms_lp_check(brms::bf(y | thres(gr = gg) ~ b1 * x, b1 ~ 1, nl = TRUE),
                brms::cumulative(), d, fit)
})

test_that("an nl body in one response of a multivariate model", {
  skip_unless_brms_fit()
  d <- nlb_ac(5007)
  d$z2 <- 0.3 + 0.5 * d$x + rnorm(nrow(d))
  fit <- frm(mvbf(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
                  bf(z2 ~ x), rescor = FALSE), data = d,
             start = list(beta = c("y_b1_(Intercept)" = 1,
                                   "y_b2_(Intercept)" = 0.1)))
  brms_lp_check(brms::mvbf(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1,
                                    nl = TRUE),
                           brms::bf(z2 ~ x), rescor = FALSE),
                gaussian(), d, fit)
})

test_that("an nl body with a random effect under an ordinal family", {
  skip_unless_brms_fit()
  d <- nlb_ord(4941)
  fit <- frm(bf(y ~ b1 * x, b1 ~ 1 + (1 | g), nl = TRUE),
             family = cumulative(), data = d)
  brms_lp_check(brms::bf(y ~ b1 * x, b1 ~ 1 + (1 | g), nl = TRUE),
                brms::cumulative(), d, fit, joint = TRUE)
})

test_that("an nl body with every core autocor term is brms's model", {
  skip_unless_brms_fit()
  st <- list(beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.1))
  specs <- list(
    ar_cov = quote(ar(time = t, gr = g, cov = TRUE)),
    ar = quote(ar(time = t, gr = g, cov = FALSE)),
    ma_cov = quote(ma(time = t, gr = g, cov = TRUE)),
    ma = quote(ma(time = t, gr = g, cov = FALSE)),
    arma_cov = quote(arma(time = t, gr = g, cov = TRUE)),
    arma = quote(arma(time = t, gr = g, cov = FALSE)))
  for (i in seq_along(specs)) {
    d <- nlb_ac(1580 + i)
    tm <- eval(call("~", specs[[i]]))
    fit <- frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, autocor = tm),
               data = d, start = st)
    brms_lp_check(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
                           autocor = tm), gaussian(), d, fit)
  }
  d <- nlb_ac(1592)
  d$y <- d$y + rnorm(nlevels(d$g), 0, 0.5)[as.integer(d$g)]
  fit <- frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
                autocor = ~ cosy(time = t, gr = g)), data = d, start = st)
  brms_lp_check(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
                         autocor = ~ cosy(time = t, gr = g)),
                gaussian(), d, fit)
  d <- nlb_ac(1590, n = 40)
  fit <- frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
                autocor = ~ unstr(time = t, gr = g)), data = d, start = st)
  brms_lp_check(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
                         autocor = ~ unstr(time = t, gr = g)),
                gaussian(), d, fit)
  d <- nlb_ac(158)
  M <- diag(nrow(d))
  fit <- frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
                autocor = ~ fcor(M)), data = d, data2 = list(M = M))
  brms_lp_check(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
                         autocor = ~ fcor(M)), gaussian(), d, fit,
                data2 = list(M = M))
})

test_that("a body from the default starts under Gamma and inverse.gaussian", {
  skip_unless_brms_fit()
  set.seed(3001)
  d <- data.frame(x = round(rnorm(200), 3))
  d$y <- exp(0.4 + 0.3 * d$x + rnorm(200, 0, 0.4))
  fit <- suppressMessages(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                              family = Gamma(), data = d))
  brms_lp_check(brms::bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), Gamma(), d,
                fit)
  fit <- suppressMessages(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                              family = inverse.gaussian(), data = d))
  brms_lp_check(brms::bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                inverse.gaussian(), d, fit)
})

test_that("the vignette's inv_logit() body is brms's model", {
  skip_unless_brms_fit()
  set.seed(4001)
  d <- data.frame(ability = rnorm(200))
  d$y <- rbinom(200, 1, 0.33 + 0.67 * plogis(1.2 * d$ability))
  fit <- frm(bf(y ~ 0.33 + 0.67 * inv_logit(eta), eta ~ ability, nl = TRUE),
             family = bernoulli("identity"), data = d)
  brms_lp_check(brms::bf(y ~ 0.33 + 0.67 * inv_logit(eta), eta ~ ability,
                         nl = TRUE), brms::bernoulli("identity"), d, fit)
})

test_that("the families with several locations read brms's one body", {
  skip_unless_brms_fit()
  set.seed(5002)
  n <- 200
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
  run <- function(f_txt, fam, bfam) {
    f <- stats::as.formula(f_txt)
    fit <- suppressWarnings(frm(bf(f, b0 + b1 ~ 1, nl = TRUE), family = fam,
                                data = d, start = list(beta = c(0.1, 0.3))))
    brms_lp_check(brms::bf(f, b0 + b1 ~ 1, nl = TRUE), bfam, d, fit)
  }
  run("cat ~ b0 + b1 * x", categorical(), brms::categorical())
  run("cnt | trials(ntr) ~ b0 + b1 * x", multinomial(), brms::multinomial())
  run("sim ~ b0 + b1 * x", dirichlet(), brms::dirichlet())
  run("cnt | trials(ntr) ~ b0 + b1 * x", dirichlet_multinomial(),
      brms::dirichlet_multinomial())
  run("sim ~ b0 + b1 * x", logistic_normal(), brms::logistic_normal())
  run("sim ~ b0 + b1 * x", brmsfamily("dirichlet2"),
      brms::brmsfamily("dirichlet2"))
  # theta1 ~ 1 so that brms carries an intercept and not a simplex, which
  # the translator has no rule for
  fit <- suppressWarnings(frm(bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, theta1 ~ 1,
                                 nl = TRUE),
                              family = mixture(gaussian(), gaussian()),
                              data = d, start = list(beta = c(0.1, 0.3))))
  brms_lp_check(brms::bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, theta1 ~ 1, nl = TRUE),
                brms::mixture(gaussian(), gaussian()), d, fit)
})
