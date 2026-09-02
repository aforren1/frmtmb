# nlf(): a nonlinear formula for ONE parameter.
#
# `nl = TRUE` says the RESPONSE formula is a nonlinear body, so only mu
# can be nonlinear. nlf() names the parameter, so any of them can be,
# and bodies can be chained. The two spellings meet in one place - a
# per-dpar `nl_body` in the spec - and the first test here is the gate
# that says so: the composed spelling and the flag spelling are the same
# fit to the last bit.

nlf_dat <- local({
  set.seed(111)
  n <- 200
  x <- stats::runif(n, 0, 5)
  data.frame(y = 2.5 * exp(-0.7 * x) + stats::rnorm(n, 0, 0.15), x = x)
})

test_that("nlf() composition is the nl = TRUE model exactly", {
  st <- list(beta = c(1, 0.3))
  direct <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) +
                  gaussian(), data = nlf_dat, start = st)
  # brms writes it this way round; the flag is on the bf() and nlf()
  # supplies mu's body through a named parameter
  composed <- frm(bf(y ~ nlmu, nl = TRUE) +
                    nlf(nlmu ~ a * exp(-b * x)) + lf(a ~ 1, b ~ 1) +
                    gaussian(), data = nlf_dat, start = st)
  # and brms's own idiom `bf(y ~ a) + nlf(a ~ ...)`, where nlf() has
  # already declared the name the response formula uses, so the flag is
  # redundant (brms insists on it; we do not)
  flagless <- frm(bf(y ~ nlmu) + nlf(nlmu ~ a * exp(-b * x)) +
                    lf(a ~ 1, b ~ 1) + gaussian(), data = nlf_dat,
                  start = st)

  expect_equal(as.numeric(logLik(composed)), as.numeric(logLik(direct)),
               tolerance = 1e-10)
  expect_equal(as.numeric(logLik(flagless)), as.numeric(logLik(direct)),
               tolerance = 1e-10)
  expect_equal(unlist(composed$estimates), unlist(direct$estimates),
               tolerance = 1e-10)
  expect_equal(vcov(composed), vcov(direct), tolerance = 1e-10)
  expect_identical(dimnames(vcov(composed)), dimnames(vcov(direct)))
  expect_equal(predict(composed), predict(direct), tolerance = 1e-10)
})

test_that("a nonlinear sigma with a linear mu matches a hand-rolled tape", {
  set.seed(20)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n, -1, 1))
  d$y <- 1 + 0.8 * d$x + stats::rnorm(n, 0, exp(-0.5 + 0.6 * d$z))

  # nl = TRUE cannot spell this: mu stays a design, sigma is a body
  fit <- frm(bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1) +
               gaussian(), data = d)

  yv <- d$y; xv <- d$x; zv <- d$z
  nll_ref <- function(p) {
    -sum(RTMB::dnorm(yv, p$b0 + p$b1 * xv,
                     exp(p$a0 + p$bz * zv), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b0 = 0, b1 = 0, a0 = 0, bz = 0),
                         silent = TRUE)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(iter.max = 1000, eval.max = 1000))

  # the structural claim: the two objectives are the same function, so
  # they agree away from the optimum too
  pv <- c(0.7, 0.4, -0.3, 0.9)
  expect_equal(fit$obj$fn(pv), obj$fn(pv), tolerance = 1e-10)

  expect_equal(as.numeric(logLik(fit)), -opt$objective, tolerance = 1e-8)
  expect_equal(unname(c(fixef(fit)$mu, fixef(fit)$a, fixef(fit)$b)),
               unname(opt$par), tolerance = 1e-5)
  # the log link is applied to the BODY's value, as it is in brms
  expect_equal(unname(predict(fit, dpar = "sigma", type = "response")),
               exp(opt$par[3] + opt$par[4] * zv), tolerance = 1e-5)
})

test_that("nonlinear bodies chain to any depth in dependency order", {
  st <- list(beta = 0.5)
  ref <- frm(bf(y ~ exp(bb) * x, bb ~ 1, nl = TRUE) + gaussian(),
             data = nlf_dat, start = st)
  chain <- frm(bf(y ~ a, nl = TRUE) + nlf(a ~ cc * x) +
                 nlf(cc ~ exp(bb)) + lf(bb ~ 1) + gaussian(),
               data = nlf_dat, start = st)
  deep <- frm(bf(y ~ p1, nl = TRUE) + nlf(p1 ~ p2 * x) +
                nlf(p2 ~ p3 + 0) + nlf(p3 ~ p4 * 1) +
                nlf(p4 ~ exp(p5)) + lf(p5 ~ 1) + gaussian(),
              data = nlf_dat, start = st)

  expect_equal(as.numeric(logLik(chain)), as.numeric(logLik(ref)),
               tolerance = 1e-10)
  expect_equal(as.numeric(logLik(deep)), as.numeric(logLik(ref)),
               tolerance = 1e-10)
  expect_equal(predict(deep), predict(ref), tolerance = 1e-10)
  # the parameters a body reads are computed before it
  expect_identical(names(chain$spec$responses$y$dpars),
                   c("bb", "cc", "a", "mu", "sigma"))
})

test_that("a cycle among the bodies is refused by name", {
  expect_error(
    frm(bf(y ~ a, nl = TRUE) + nlf(a ~ b1 * x) + nlf(b1 ~ a + 1) +
          gaussian(), data = nlf_dat, dry_run = "spec"),
    "depend on each other in a cycle"
  )
  expect_error(
    frm(bf(y ~ a, nl = TRUE) + nlf(a ~ a + 1) + gaussian(),
        data = nlf_dat, dry_run = "spec"),
    "depend on each other in a cycle"
  )
})

test_that("random effects live in the nonlinear parameters either way", {
  set.seed(112)
  n_g <- 25; n_per <- 20
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- stats::runif(n_g * n_per, 0, 5)
  a_g <- 2.5 + stats::rnorm(n_g, 0, 0.5)
  d <- data.frame(y = a_g[g] * exp(-0.7 * x) +
                    stats::rnorm(length(x), 0, 0.15), x = x, g = g)
  st <- list(beta = c(2, 0.5))

  direct <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1,
                   nl = TRUE) + gaussian(), data = d, start = st)
  composed <- frm(bf(y ~ nlmu, nl = TRUE) + nlf(nlmu ~ a * exp(-b * x)) +
                    lf(a ~ 1 + (1 | g), b ~ 1) + gaussian(),
                  data = d, start = st)
  expect_equal(as.numeric(logLik(composed)), as.numeric(logLik(direct)),
               tolerance = 1e-10)
  expect_equal(vcov(composed), vcov(direct), tolerance = 1e-10)
  expect_equal(ranef(composed)[[1]], ranef(direct)[[1]],
               tolerance = 1e-10)

  # and the same random effect behind a nonlinear sigma, which has no
  # nl = TRUE spelling at all
  set.seed(21)
  n <- 400
  d2 <- data.frame(x = stats::rnorm(n), z = stats::runif(n, -1, 1),
                   g = factor(rep(1:20, 20)))
  u <- stats::rnorm(20, 0, 0.4)
  d2$y <- 1 + 0.8 * d2$x +
    stats::rnorm(n, 0, exp(-0.5 + 0.6 * d2$z + u[d2$g]))
  fit <- frm(bf(y ~ x) + nlf(sigma ~ a + b * z) +
               lf(a ~ 1 + (1 | g), b ~ 1) + gaussian(), data = d2)
  expect_named(VarCorr(fit), "a: 1 | g")
  expect_identical(dim(ranef(fit)[[1]]), c(20L, 1L))
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("frm_ode() composes inside an nlf() body", {
  skip_if_not_installed("RTMBode")
  pk_dyn <- function(t, y, p) {
    list(c(-p[1] * y[1], p[1] * y[1] - p[2] / p[3] * y[2]))
  }
  set.seed(2026)
  tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
  id <- factor(rep(seq_len(4), each = length(tt)))
  time <- rep(tt, 4)
  ka <- exp(stats::rnorm(4, 0, 0.3))[as.integer(id)]
  ke <- exp(log(0.2) + stats::rnorm(4, 0, 0.25))[as.integer(id)]
  mu <- 100 * ka / (10 * (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
  d <- data.frame(id = id, time = time, dose = 100,
                  conc = mu + stats::rnorm(length(mu), 0, 0.3))
  st <- list(beta = c(0, log(0.25), log(8)))

  direct <- frm(bf(conc ~ frm_ode(pk_dyn, init = list(dose, 0),
                                  times = time,
                                  parms = list(exp(lka), exp(lke),
                                               exp(lV)),
                                  group = id, output = 2L),
                   lka ~ 1, lke ~ 1, lV ~ 1, nl = TRUE) + gaussian(),
                data = d, start = st)
  vianlf <- frm(bf(conc ~ pk, nl = TRUE) +
                  nlf(pk ~ frm_ode(pk_dyn, init = list(dose, 0),
                                   times = time,
                                   parms = list(exp(lka), exp(lke),
                                                exp(lV)),
                                   group = id, output = 2L)) +
                  lf(lka ~ 1, lke ~ 1, lV ~ 1) + gaussian(),
                data = d, start = st)
  expect_equal(as.numeric(logLik(vianlf)), as.numeric(logLik(direct)),
               tolerance = 1e-10)
  expect_equal(predict(vianlf), predict(direct), tolerance = 1e-10)

  # the within-group constancy check reads every body, not just mu's
  d$phase <- factor(ifelse(d$time <= 2, "early", "late"))
  expect_error(
    frm(bf(conc ~ pk, nl = TRUE) +
          nlf(pk ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                           parms = list(exp(lka), exp(lke), exp(lV)),
                           group = id, output = 2L)) +
          lf(lka ~ 1, lke ~ 1 + phase, lV ~ 1) + gaussian(),
        data = d, dry_run = "frame",
        start = list(beta = c(0, log(0.25), 0, log(8)))),
    "not constant within 'id'"
  )
})

test_that("post-processing follows the nonlinear parameter a body names", {
  set.seed(20)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n, -1, 1))
  d$y <- 1 + 0.8 * d$x + stats::rnorm(n, 0, exp(-0.5 + 0.6 * d$z))
  fit <- frm(bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1) +
               gaussian(), data = d)

  expect_equal(predict(fit, newdata = d), predict(fit), tolerance = 1e-8)
  expect_equal(fitted(fit), predict(fit, type = "response"),
               tolerance = 1e-8)
  expect_length(residuals(fit), n)
  # mu is linear, so its delta-method standard error is unaffected
  expect_false(is.null(predict(fit, se.fit = TRUE)$se.fit))
  expect_error(predict(fit, dpar = "sigma", se.fit = TRUE),
               "se.fit is not supported")
  # the effect display finds the body's covariate, not mu's
  ce <- conditional_effects(fit, dpar = "sigma", band = "boot", boot = 5)
  expect_named(ce, "z")
  expect_error(conditional_effects(fit, dpar = "sigma", band = "wald"),
               "cannot put a wald band on a nonlinear predictor")
  # ... while mu keeps its own analytic band
  expect_named(conditional_effects(fit), "x")
  expect_true(all(c("a_Intercept", "b_Intercept") %in% variables(fit)))
  expect_s3_class(hypothesis(fit, "b_Intercept = 0"), "frmtmb_hypothesis")

  pr <- get_prior(bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1) +
                    gaussian(), data = d)
  expect_true(all(c("a", "b") %in% pr$dpar))
})

test_that("an update() delta still reaches a linear mu behind a nonlinear sigma", {
  set.seed(3)
  d <- data.frame(y = stats::rnorm(60), x = stats::rnorm(60),
                  z = stats::rnorm(60), w = stats::rnorm(60))
  fit <- frm(bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1) +
               gaussian(), data = d)
  up <- update(fit, . ~ . + w)
  expect_identical(colnames(up$frame$linpreds$y.mu$X),
                   c("(Intercept)", "x", "w"))
  # but not a mu that is itself a body
  nl <- frm(bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1) + gaussian(),
            data = d, start = list(beta = 0))
  expect_error(update(nl, . ~ . + w), "cannot be applied to a nonlinear")
})

test_that("nlf() refuses what it cannot mean", {
  d <- data.frame(y = stats::rnorm(50), x = stats::rnorm(50),
                  z = stats::rnorm(50), g = factor(rep(1:5, 10)))
  expect_error(nlf(~x), "two-sided formula naming the parameter")
  expect_error(nlf(a + b ~ exp(x)), "one parameter at a time")
  expect_error(nlf(a.b ~ x), "Invalid parameter name")
  expect_error(nlf(a ~ exp(x), a ~ 1), "both a nonlinear body")
  expect_error(bf(y ~ x, sigma ~ z) + nlf(sigma ~ a * z),
               "nlf\\(\\) sets 'sigma'")
  expect_error(bf(y ~ x) + nlf(sigma ~ a * z) + nlf(sigma ~ a + z),
               "nlf\\(\\) sets 'sigma'")
  expect_error(bf(y ~ x) + nlf(sigma ~ a * z) + lf(sigma ~ z),
               "lf\\(\\) sets 'sigma'")
  expect_error(bf(y ~ x) + bf(z ~ x) + nlf(sigma ~ a),
               "does not say which response")
  # nlf(mu ~ ...) would silently discard the bf()'s own right-hand side
  expect_error(frm(bf(y ~ x) + nlf(mu ~ a * exp(b * x)) +
                     lf(a ~ 1, b ~ 1) + gaussian(), data = d,
                   dry_run = "spec"),
               "with nothing to do")
  expect_s3_class(frm(bf(y ~ 1) + nlf(mu ~ a * exp(b * x)) +
                        lf(a ~ 1, b ~ 1) + gaussian(), data = d,
                      dry_run = "spec"), "frmtmb_spec")
  # a body whose free names are all data leaves nothing to estimate
  expect_error(frm(bf(y ~ x) + nlf(sigma ~ 1 + z) + gaussian(),
                   data = d, dry_run = "spec"),
               "nothing is left to estimate")
  # a declared parameter no body uses
  expect_error(frm(bf(y ~ x) + nlf(shape ~ a) + lf(a ~ 1) + gaussian(),
                   data = d, dry_run = "spec"),
               "not used in the model formula")
  # a nonlinear parameter that shadows a data column
  expect_error(frm(bf(y ~ x) + nlf(sigma ~ z * g) + lf(z ~ 1) +
                     gaussian(), data = d, dry_run = "frame"),
               "also name columns of the data")
  # residual correlation terms are terms, not R code, in any body
  expect_error(frm(bf(y ~ x) + nlf(sigma ~ a + ar(x, g)) + lf(a ~ 1) +
                     gaussian(), data = d, dry_run = "spec"),
               "not supported in a ")
})

test_that("nlf() reaches other families and the rest of the grammar", {
  set.seed(5)
  d <- data.frame(x = stats::rnorm(80), z = stats::rnorm(80),
                  g = factor(rep(1:8, 10)))
  d$y <- stats::rnorm(80)
  d$yg <- stats::rgamma(80, 2, 1)
  d$y2 <- stats::rnorm(80)

  # any dpar of any family, not only sigma
  gam <- frm(bf(yg ~ x) + nlf(shape ~ a + b * z) + lf(a ~ 1, b ~ 1) +
               Gamma(), data = d)
  expect_true(is.finite(as.numeric(logLik(gam))))

  # the nonlinear parameters keep the full predictor grammar
  sm <- frm(bf(y ~ x) + nlf(sigma ~ a) + lf(a ~ s(z, k = 5)) +
              gaussian(), data = d)
  expect_true(is.finite(as.numeric(logLik(sm))))

  # one response of a multivariate model can be the nonlinear one
  mv <- frm(mvbf(bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1),
                 bf(y2 ~ x)) + gaussian(), data = d)
  expect_true(is.finite(as.numeric(logLik(mv))))

  # REML integrates mu's fixed effects; a nonlinear sigma is not one
  rml <- frm(bf(y ~ x + (1 | g)) + nlf(sigma ~ a + b * z) +
               lf(a ~ 1, b ~ 1) + gaussian(), data = d, REML = TRUE)
  expect_true(is.finite(as.numeric(logLik(rml))))
  expect_identical(rml$spec$responses$y$primary_dpars, "mu")
})

test_that("a body can read another dpar's value: varPower(~ fitted(.))", {
  skip_if_not_installed("nlme")
  set.seed(42)
  n <- 400
  d <- data.frame(x = stats::runif(n, 1, 5))
  m <- 2 + 1.5 * d$x
  d$y <- m + stats::rnorm(n, 0, 0.3 * abs(m)^0.8)

  # sd_i = sigma0 * |mu_i|^theta, with sigma's log link turning the body
  # into exactly that. nlme spells it varPower(form = ~ fitted(.)); brms
  # has no spelling at all, because a body name there is a column.
  fit <- frm(bf(y ~ x) + nlf(sigma ~ ls + th * log(abs(mu))) +
               lf(ls ~ 1, th ~ 1) + gaussian(), data = d,
             start = list(betad = c(log(0.3), 0.8)))

  # the objective IS the joint Gaussian likelihood of that model
  yv <- d$y; xv <- d$x
  nll_ref <- function(p) {
    mu <- p$b0 + p$b1 * xv
    -sum(RTMB::dnorm(yv, mu, exp(p$ls + p$th * log(abs(mu))),
                     log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b0 = 1, b1 = 1, ls = log(0.3),
                                       th = 0.8), silent = TRUE)
  pv <- c(1.8, 1.4, log(0.35), 0.7)
  expect_equal(fit$obj$fn(pv), obj$fn(pv), tolerance = 1e-12)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(iter.max = 2000, eval.max = 2000))
  expect_equal(as.numeric(logLik(fit)), -opt$objective, tolerance = 1e-8)

  g <- nlme::gls(y ~ x, data = d,
                 weights = nlme::varPower(form = ~ fitted(.)),
                 method = "ML")
  # the two likelihoods are the same function: ours evaluated at gls's
  # own estimates reproduces gls's logLik. The ESTIMATORS differ - gls
  # alternates between the mean fit and the variance function against
  # iteratively updated fitted values, and stops with a gradient it
  # does not check, while this is one joint maximization - so ours ends
  # a little higher.
  gp <- c(stats::coef(g), log(g$sigma),
          stats::coef(g$modelStruct$varStruct, unconstrained = FALSE))
  expect_equal(fit$obj$fn(gp), -as.numeric(stats::logLik(g)),
               tolerance = 1e-8)
  expect_gte(as.numeric(logLik(fit)), as.numeric(stats::logLik(g)))
  expect_equal(as.numeric(logLik(fit)), as.numeric(stats::logLik(g)),
               tolerance = 1e-3)
  expect_equal(unname(fixef(fit)$mu), unname(stats::coef(g)),
               tolerance = 1e-2)
  expect_equal(exp(fixef(fit)$ls[[1]]), g$sigma, tolerance = 1e-2)
  expect_equal(fixef(fit)$th[[1]],
               unname(stats::coef(g$modelStruct$varStruct,
                                  unconstrained = FALSE)),
               tolerance = 1e-2)

  # newdata and eval_dpars take the same route through the reference
  expect_equal(predict(fit, newdata = d[1:20, ], dpar = "sigma",
                       type = "response"),
               predict(fit, dpar = "sigma",
                       type = "response")[1:20],
               tolerance = 1e-10, ignore_attr = TRUE)

  # a column of the data still wins, so a brms body keeps its meaning
  d$mu <- stats::runif(n, 1, 2)
  fit2 <- frm(bf(y ~ x) + nlf(sigma ~ ls + th * log(abs(mu))) +
                lf(ls ~ 1, th ~ 1) + gaussian(), data = d,
              start = list(betad = c(log(0.3), 0.1)))
  lp <- fit2$frame$linpreds[["y.sigma"]]
  expect_identical(names(lp$data_list), "mu")
  expect_length(lp$nl_dpar_refs, 0L)

  # and two dpars cannot each be a function of the other
  expect_error(
    frm(bf(y ~ 1) + nlf(mu ~ a * sigma) + nlf(sigma ~ b * mu) +
          lf(a ~ 1, b ~ 1) + gaussian(), data = d, dry_run = "spec"),
    "depend on each other in a cycle"
  )
})

test_that("nlf() objects print and carry their formulas", {
  x <- nlf(sigma ~ a * z, a ~ 1)
  expect_s3_class(x, "frmtmb_nlf")
  expect_named(x$nlforms, "sigma")
  expect_named(x$pforms, "a")
  expect_output(print(x), "nonlinear")
  bform <- bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1)
  expect_named(bform$nlforms, "a")
  expect_named(bform$pforms, "b")
  expect_output(print(bform), "a ~ exp\\(b \\* x\\) \\(nonlinear\\)")
  # loop = is accepted for brms source compatibility
  expect_s3_class(nlf(sigma ~ a * z, loop = FALSE), "frmtmb_nlf")
})
