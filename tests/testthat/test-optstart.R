# Starting values from a least-squares fit of mu (start_ls_plan()) and
# the restart rule (restart_result()); dev/optstart-findings.md has
# the measurements behind every case.

# y = 1 + k z1 + ... + b_g + e: a covariate that explains nearly all of
# the response, so that sd(y) is near k and the residual spread is 1
optstart_wide <- function(seed, n = 500, p = 10, k = 1e3) {
  set.seed(seed)
  Z <- matrix(rnorm(n * p), n, p)
  colnames(Z) <- paste0("z", seq_len(p))
  G <- 10
  d <- data.frame(Z, g = factor(rep(seq_len(G), length.out = n)))
  d$y <- 1 + k * d$z1 + rnorm(G, 0, 0.5)[d$g] + rnorm(n)
  fo <- stats::as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                                "+ (1 | g)"))
  list(d = d, fo = fo)
}

test_that("a covariate that explains most of the response reaches lme4", {
  skip_if_not_installed("lme4")
  # 0.70.0 started sigma at sd(y), near 1e3, and every slope at 0; these
  # three seeds stopped 61.6, 32.7 and 33.3 log-likelihood units below
  # lme4 with a warning (dev/optstart-search.R)
  for (s in c(1, 14, 18)) {
    m <- optstart_wide(s)
    f <- frm(m$fo, data = m$d)
    l4 <- lme4::lmer(m$fo, data = m$d, REML = FALSE)
    ll4 <- as.numeric(stats::logLik(l4))
    expect_identical(f$opt$convergence, 0L, info = paste("seed", s))
    expect_lt(abs(as.numeric(stats::logLik(f)) - ll4),
              1e-8 * abs(ll4), label = paste("seed", s, "gap to lme4"))
  }
})

test_that("mu starts at weighted least squares, sigma at its spread", {
  set.seed(3)
  d <- data.frame(x = rnorm(40), w = runif(40, 0.5, 2))
  d$y <- 100 + 50 * d$x + rnorm(40)
  st <- par_template(bf(y | weights(w) ~ x), data = d, family = gaussian())
  ref <- stats::lm(y ~ x, data = d, weights = w)
  expect_equal(unname(st$beta), unname(stats::coef(ref)))
  r <- stats::residuals(ref)
  expect_equal(unname(st$betad), log(sqrt(sum(d$w * r^2) / sum(d$w))))
  # the start is the spread of the residuals, not of the response
  expect_lt(exp(unname(st$betad)), stats::sd(d$y) / 10)
})

test_that("each family's scale starts near its value", {
  # data from the family's own simulator with a covariate that explains
  # most of the response; 0.70.0 started these between 0.19 and 13.4
  # times the true value (dev/optstart-relations.R)
  reg <- frmtmb:::family_registry
  cases <- list(
    list("gaussian", gaussian(), list(sigma = 1.5), identity, c(10, 20)),
    list("student", student(), list(sigma = 1.5, nu = 4), identity,
         c(10, 20)),
    # exgaussian starts each of sigma and beta at half the residual
    # spread, sqrt(2) times smaller than equal true values
    list("exgaussian", exgaussian(), list(sigma = 1.5, beta = 1.5), identity,
         c(10, 20), 1 / sqrt(2)),
    list("asym_laplace", asym_laplace(), list(sigma = 1, quantile = 0.5),
         identity, c(10, 20)),
    list("lognormal", lognormal(), list(sigma = 0.4), identity, c(1, 0.8)),
    list("Gamma", Gamma(link = "log"), list(shape = 3), exp, c(1, 0.8)),
    list("weibull", weibull(), list(shape = 2), exp, c(1, 0.8)),
    list("beta", Beta(), list(phi = 20), stats::plogis, c(0.3, 1)),
    list("negbinomial", negbinomial(), list(shape = 3), exp, c(1, 0.8)))
  n <- 2000
  for (cs in cases) {
    set.seed(1)
    x <- rnorm(n)
    mu <- cs[[4]](cs[[5]][1] + cs[[5]][2] * x)
    y <- reg[[cs[[1]]]]()$sim(c(list(mu = mu), lapply(cs[[3]], rep, n)),
                              list(), n)
    st <- par_template(bf(y ~ x), data = data.frame(x = x, y = y),
                       family = cs[[2]])
    for (dp in intersect(names(cs[[3]]), c("sigma", "shape", "phi", "beta"))) {
      ratio <- exp(st$betad[[paste0(dp, "_(Intercept)")]]) / cs[[3]][[dp]]
      ratio <- ratio / (if (length(cs) > 5) cs[[6]] else 1)
      expect_true(ratio > 0.8 && ratio < 1.25,
                  info = paste(cs[[1]], dp, "start / true", ratio))
    }
    expect_true(abs(st$beta[["x"]] / cs[[5]][2] - 1) < 0.2,
                info = paste(cs[[1]], "slope start / true",
                             st$beta[["x"]] / cs[[5]][2]))
  }
})

test_that("a start or a prior leaves the rest of the start to the rule", {
  set.seed(5)
  d <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  d$y <- 50 + 30 * d$x + rnorm(60)
  fo <- bf(y ~ x + (1 | g))
  cold <- par_template(fo, data = d, family = gaussian())
  # the entry start names is the user's; the others are the defaults
  st <- par_template(fo, data = d, family = gaussian(),
                     start = list(beta = c(x = 5)))
  expect_identical(st$beta[["x"]], 5)
  expect_identical(st$beta[["(Intercept)"]], cold$beta[["(Intercept)"]])
  expect_identical(st$betad, cold$betad)
  # a prior on an ordinary coefficient moves no starting value
  pr <- par_template(fo, data = d, family = gaussian(),
                     prior = set_prior("normal(0, 1)", class = "b"))
  expect_identical(unclass(pr), unclass(cold))
})

test_that("where the rule does not apply the start is the initializer's", {
  set.seed(5)
  d <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)),
                  cn = rbinom(60, 1, 0.2))
  d$y <- 50 + 30 * d$x + rnorm(60)
  d$yc <- stats::rpois(60, exp(1 + 0.5 * d$x))
  old <- c(mean(d$y), 0)
  old_sigma <- log(stats::sd(d$y))
  same <- function(st, beta, betad) {
    expect_identical(unname(st$beta), beta)
    expect_identical(unname(st$betad), betad)
  }
  # a censored response, and a family with no scale
  same(par_template(bf(y | cens(cn) ~ x), data = d, family = gaussian()),
       old, old_sigma)
  same(par_template(bf(yc ~ x), data = d, family = poisson()),
       c(log(mean(d$yc) + 0.1), 0), NULL)
  # an intercept-only mu has nothing to fit
  same(par_template(bf(y ~ 1 + (1 | g)), data = d, family = gaussian()),
       mean(d$y), old_sigma)
  # a nonlinear mu has no design of its own
  st <- par_template(bf(y ~ a + b * x, a ~ 1, b ~ 1, nl = TRUE), data = d,
                     family = gaussian())
  expect_identical(unname(st$beta), c(0, 0))
  expect_identical(unname(st$betad), old_sigma)
  # a smooth in mu
  st <- par_template(bf(y ~ s(x)), data = d, family = gaussian())
  same(st, c(mean(d$y), rep(0, length(st$beta) - 1L)), old_sigma)
  # a design whose condition number is above 1e5 (degree-8 raw powers)
  dp <- data.frame(u = stats::runif(60, 1, 2))
  dp$y <- sin(3 * dp$u) + stats::rnorm(60, 0, 0.2)
  for (k in 1:8) dp[[paste0("p", k)]] <- dp$u^k
  same(par_template(bf(y ~ p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8),
                    data = dp, family = gaussian()),
       c(mean(dp$y), rep(0, 8)), log(stats::sd(dp$y)))
  # and a link the rule is not written for
  same(par_template(bf(yc + 1 ~ x), data = d, family = Gamma()),
       c(1 / mean(d$yc + 1), 0),
       log(max(mean(d$yc + 1)^2 / stats::var(d$yc + 1), 0.1)))
})

test_that("a restart that did not move keeps the converged run's code", {
  run <- function(obj, code, msg, par, it = 10L) {
    list(par = par, objective = obj, convergence = code, message = msg,
         iterations = it, evaluations = c("function" = it + 2L,
                                          gradient = it + 1L))
  }
  conv <- run(10, 0L, "relative convergence (4)", c(1, 2))
  still <- run(10, 1L, "false convergence (8)", c(1, 2), 1L)
  rr <- frmtmb:::restart_result
  # the run's code, message, iterations and evaluations, with the
  # restart's own kept aside
  r <- rr(conv, still)
  expect_identical(r[c("convergence", "message", "iterations",
                       "evaluations")],
                   conv[c("convergence", "message", "iterations",
                          "evaluations")])
  expect_identical(r$restart,
                   still[c("convergence", "message", "iterations",
                           "evaluations")])
  # a gain within the tolerance but a move keeps the restart's code
  moved <- run(10 - 1e-12, 1L, "false convergence (8)", c(1, 2.1))
  expect_identical(rr(conv, moved), moved)
  # a real gain keeps it too, and a restart that ends higher is dropped
  far <- run(9, 1L, "false convergence (8)", c(1, 2))
  expect_identical(rr(conv, far), far)
  expect_null(rr(conv, run(11, 0L, "relative convergence (4)", c(1, 2))))
  # a converged restart, or a first run that did not converge, is kept
  expect_identical(rr(still, conv), conv)
  # the review's case: a caller's loose rel.tol = 1e-4 at an objective
  # near 1e4 allows a gain of 1, but a restart that gains 0.9 and moves
  # 5 units is a different stop and keeps its code
  first <- run(1e4, 0L, "X-convergence (3)", c(0, 0))
  long <- run(1e4 - 0.9, 10L, "iteration limit reached", c(5, 0), 150L)
  expect_identical(rr(first, long, list(optCtrl = list(rel.tol = 1e-4))),
                   long)
  # the tolerance follows the optimizer that ran: L-BFGS-B's factr
  near <- run(10 - 3e-9, 1L, "ABNORMAL_TERMINATION_IN_LNSRCH", c(1, 2))
  expect_identical(rr(conv, near, list(optimizer = "optim"))$convergence,
                   0L)
  expect_identical(rr(conv, near, list(optimizer = "optim",
                                       optCtrl = list(factr = 1)))$convergence,
                   1L)
  expect_identical(rr(conv, near)$convergence, 1L)
  expect_identical(rr(conv, near, list(optCtrl = list(rel.tol = 1e-9)))$
                     convergence, 0L)
})

test_that("a group sd far above the residual sd converges at lme4", {
  skip_if_not_installed("lme4")
  # a group sd 1000 times the residual sd (dev/optrap-rev2-b1.R re1000):
  # from 0.70.0's start (sd(y) for sigma, given here as `start`) seeds 6
  # and 14, and from the least-squares start seed 8, stopped at lme4's
  # maximum at relative convergence, and the restart from there ended in
  # one iteration at "false convergence (8)" with the same objective
  for (s in c(6, 14, 8)) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 1000)[d$g] + rnorm(100)
    st <- if (s != 8) {
      list(beta = c(mean(d$y), 0), betad = log(stats::sd(d$y)))
    }
    f <- frm(y ~ x + (1 | g), data = d, start = st)
    # lme4's own convergence check is unsure at a group sd this large;
    # its maximum is the one frmtmb reaches either way
    l4 <- allow_warnings(lme4::lmer(y ~ x + (1 | g), data = d, REML = FALSE),
                         c("unable to evaluate scaled gradient",
                           "Model failed to converge"))
    ll4 <- as.numeric(stats::logLik(l4))
    expect_identical(f$opt$convergence, 0L, info = paste("seed", s))
    expect_lt(abs(as.numeric(stats::logLik(f)) - ll4), 1e-8 * abs(ll4),
              label = paste("seed", s, "gap to lme4"))
  }
})

test_that("a mi() coefficient is stepped in its variable's units", {
  skip_if_not_installed("mice")
  skip_if_not_installed("numDeriv")
  # brms_missings' fit_imp2 (lane remeasure, R1): the mi(chl) column is a
  # zero placeholder, so its coefficient got bmi's unit (8) where chl
  # spreads 45; at 0.70.0 the fit ended at "false convergence (8)" with
  # every standard error NaN, at the optimum 0.67.0 reached
  data("nhanes", package = "mice", envir = environment())
  f <- frm(bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
             set_rescor(FALSE), data = nhanes)
  expect_identical(f$opt$convergence, 0L)
  nm <- names(par_template(f)$beta)
  u <- f$par_units[seq_along(nm)]
  ub <- u[nm == "bmi_(Intercept)"]
  expect_identical(unname(ub / u[nm == "bmi_michl"]),
                   2^round(log2(stats::sd(nhanes$chl, na.rm = TRUE))))
  # the standard errors against a Richardson Jacobian of the gradient
  H <- numDeriv::jacobian(function(x) as.numeric(f$obj$gr(x)), f$opt$par)
  ref <- sqrt(diag(solve((H + t(H)) / 2)))[seq_along(nm)]
  se <- sqrt(diag(vcov(f)))[gsub("[()]", "", nm)]
  expect_true(all(is.finite(se)))
  expect_lt(max(abs(se / ref - 1)), 1e-4)
})

test_that("an imputed value is not data for the start or the unit", {
  # mi() stores an imputed value as 0 in the frame; 0.70.0 read those
  # zeros as observations for the response unit (lane remeasure, R1)
  set.seed(9)
  d <- data.frame(x = rnorm(60))
  d$y <- 50 + 30 * d$x + rnorm(60)
  d$y[sample(60, 20)] <- NA
  st <- par_template(bf(y | mi() ~ x), data = d, family = gaussian())
  ref <- stats::lm(y ~ x, data = d)
  expect_equal(unname(st$beta), unname(stats::coef(ref)))
  fr <- frm(bf(y | mi() ~ x), data = d, family = gaussian(),
            dry_run = "frame")
  u <- frmtmb:::predictor_units(fr)
  expect_identical(unname(u), numeric(0))
})
