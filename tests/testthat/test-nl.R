test_that("nonlinear fixed-effects model matches nls", {
  set.seed(111)
  n <- 200
  x <- runif(n, 0, 5)
  y <- 2.5 * exp(-0.7 * x) + rnorm(n, 0, 0.15)
  dd <- data.frame(y = y, x = x)

  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) + gaussian(),
             data = dd, start = list(beta = c(1, 0.3)))
  ref <- nls(y ~ a * exp(-b * x), data = dd, start = list(a = 1, b = 0.3))

  expect_lt(abs(fixef_by_dpar(fit)$a[[1]] - coef(ref)[["a"]]), 1e-4)
  expect_lt(abs(fixef_by_dpar(fit)$b[[1]] - coef(ref)[["b"]]), 1e-4)
  # ML sigma^2 = RSS/n at the same coefficients
  sig_ml <- sqrt(sum(residuals(ref)^2) / n)
  expect_lt(abs(exp(fixef_by_dpar(fit)$sigma[[1]]) - sig_ml), 1e-4)
})

test_that("nonlinear mixed model matches a hand-rolled reference", {
  set.seed(112)
  n_g <- 25; n_per <- 20
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- runif(n_g * n_per, 0, 5)
  a_g <- 2.5 + rnorm(n_g, 0, 0.5)
  y <- a_g[g] * exp(-0.7 * x) + rnorm(length(x), 0, 0.15)
  dd <- data.frame(y = y, x = x, g = g)

  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE) +
               gaussian(),
             data = dd, start = list(beta = c(2, 0.5)))

  yv <- dd$y; xv <- dd$x; gi <- as.integer(dd$g)
  nll_ref <- function(p) {
    nll <- -sum(RTMB::dnorm(p$u, 0, exp(p$lsd), log = TRUE))
    a <- p$a0 + p$u[gi]
    mu <- a * exp(-p$b0 * xv)
    nll - sum(RTMB::dnorm(yv, mu, exp(p$ls), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(a0 = 2, b0 = 0.5, ls = 0, lsd = 0,
                              u = numeric(n_g)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)

  # nlpar random effects show up in ranef and VarCorr
  expect_length(varcorr_matrices(fit), 1)
  expect_identical(dim(ranef(fit)[[1]]), c(as.integer(n_g), 1L))
})

test_that("nl prediction and post-processing", {
  set.seed(113)
  n <- 150
  x <- runif(n, 0, 5)
  dd <- data.frame(y = 2 * exp(-0.5 * x) + rnorm(n, 0, 0.1), x = x)
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) + gaussian(),
             data = dd, start = list(beta = c(1, 0.3)))

  expect_equal(frm_linpred(fit, newdata = dd), frm_linpred(fit),
               tolerance = 1e-8)
  expect_equal(unname(fitted(fit)[, "Estimate"]),
               unname(frm_linpred(fit, type = "response")), tolerance = 1e-8)
  a_hat <- frm_linpred(fit, dpar = "a")
  expect_lt(stats::sd(a_hat), 1e-10)   # intercept-only nlpar is constant
  nd <- data.frame(x = c(0, 1, 2))
  p <- frm_linpred(fit, newdata = nd)
  expect_equal(p[1], fixef_by_dpar(fit)$a[[1]], tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_error(frm_linpred(fit, se.fit = TRUE), "se.fit is not supported")
  expect_length(residuals(fit)[, "Estimate"], n)
})

test_that("nl validation errors are clear", {
  # refused when the model is assembled, as brms refuses it at brm():
  # the formulas may still arrive with + lf() (lane formrobust)
  expect_error(frm(bf(y ~ a * exp(-b * x), nl = TRUE), data = NULL,
                   dry_run = "spec"),
               "nonlinear-parameter formula")
  expect_error(frm(bf(y ~ a * exp(-b * x), a ~ 1, cc ~ 1, nl = TRUE) +
                     gaussian(),
                   data = NULL, dry_run = "spec"),
               "not used in the model formula")
})

# --- reserved nonlinear parameter names -------------------------------
# Four usability defects of hierarchical nl models, wt-api. A nonlinear
# parameter named after one of the family's own dpars used to die inside
# model.frame() with "object 'mu' not found", naming nothing; one named
# after a par-template component fits, but `start` then means the
# component and reported ITS length.

nl_reserved_data <- function(n_id = 8, seed = 7) {
  set.seed(seed)
  w <- seq(1, 6, by = 0.5)
  d <- do.call(rbind, lapply(seq_len(n_id), function(i) {
    data.frame(id = i, group = i %% 2, w = w, logw = log(w),
               I = rexp(length(w), rate = 1 / exp(1.2 - 1.5 * log(w))))
  }))
  d$id <- factor(d$id)
  d
}

test_that("a nonlinear parameter named after a family dpar is refused by name",
          {
  d <- nl_reserved_data()
  # the reported spelling: `mu ~ 1 + (1 | id)` alongside nl = TRUE
  expect_error(
    frm(bf(I ~ mu - chi * logw, mu ~ 1 + (1 | id), chi ~ 1 + group,
           nl = TRUE),
        family = exponential(link = "log"), data = d),
    "distributional parameter of family 'exponential'")
  expect_error(
    frm(bf(I ~ mu - chi * logw, mu ~ 1 + (1 | id), chi ~ 1 + group,
           nl = TRUE),
        family = exponential(link = "log"), data = d),
    "This family reserves: mu")
  # and at SPEC time, so par_template() refuses it before any fit
  expect_error(
    par_template(bf(I ~ mu - chi * logw, mu ~ 1 + (1 | id),
                    chi ~ 1 + group, nl = TRUE),
                 data = d, family = exponential(link = "log")),
    "cannot also be a nonlinear parameter")
  # any family, not just this one
  expect_error(
    frm(bf(I ~ mu * chi, mu ~ 1, chi ~ 1, nl = TRUE), gaussian(), data = d),
    "distributional parameter of family 'gaussian'")
})

test_that("a body that names its own parameter is refused, with the data checked first", {
  d <- nl_reserved_data()
  # no formula for mu and no column called mu: the old message was R's
  # own "object 'mu' not found" from eval(predvars, data, env)
  expect_error(
    frm(bf(I ~ mu - chi * logw, chi ~ 1 + group, nl = TRUE),
        family = exponential(link = "log"), data = d),
    "refers to 'mu' itself")
  # a REAL column of that name still wins, as it does for a dpar
  # reference, so this model keeps fitting
  d2 <- d
  d2$mu <- 1
  f <- frm(bf(I ~ mu * apo - chi * logw, apo ~ 1, chi ~ 1 + group,
              nl = TRUE),
           family = exponential(link = "log"), data = d2)
  # mu is computed by the body, so it contributes no coefficient block;
  # what matters is that the fit happened at all
  expect_setequal(names(fixef_by_dpar(f)), c("apo", "chi"))
  expect_true(all(is.finite(fixef(f, flatten = TRUE))))
  # an nlf() body that names ITSELF is the same fault under another
  # spelling
  expect_error(
    frm(bf(I ~ apo - chi * logw, apo ~ 1, chi ~ 1 + group, nl = TRUE) +
          nlf(sigma ~ sigma + 1),
        family = gaussian(), data = d),
    "refers to 'sigma' itself")
})

test_that("a body reading ANOTHER dpar's value is untouched", {
  # the deliberate variance-function extension: `sigma` in mu's body is
  # that parameter's per-row value, not a nonlinear parameter, and the
  # refusals above must not reach it
  set.seed(2)
  dd <- data.frame(x = rnorm(120))
  dd$y <- 2 + dd$x + rnorm(120)
  f <- frm(bf(y ~ sigma * x + a, a ~ 1, nl = TRUE), gaussian(), data = dd)
  expect_true(is.finite(fixef_by_dpar(f)$a[[1]]))
})

test_that("a nonlinear parameter named after a template component fits, and start names the collision", {
  d <- nl_reserved_data()
  # one test per par-template component that a plain hierarchical model
  # carries. Each name FITS - the refusal is only about `start`.
  for (nm in c("beta", "b", "theta")) {
    body <- stats::as.formula(paste0("I ~ ", nm, " - chi * logw"))
    par <- stats::as.formula(paste0(nm, " ~ 1 + (1 | id)"))
    f <- frm(bf(body, par, chi ~ 1 + group, nl = TRUE),
             family = exponential(link = "log"), data = d)
    expect_true(paste0(nm, "_(Intercept)") %in% names(fixef(f, flatten = TRUE)),
                info = nm)
    st <- list(1)
    names(st) <- nm
    # `start$<nm>` sets the COMPONENT. Whether that is a length error or
    # a silent success depends on the component's length, so both paths
    # have to carry the explanation.
    seen <- character(0)
    msg <- tryCatch({
      withCallingHandlers(
        frm(bf(body, par, chi ~ 1 + group, nl = TRUE),
            family = exponential(link = "log"), data = d, start = st),
        warning = function(w) {
          seen <<- c(seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      paste(seen, collapse = " ")
    }, error = function(e) conditionMessage(e))
    expect_match(msg, "not the nonlinear parameter", info = nm)
    expect_match(msg, "par_template\\(\\) lists both", info = nm)
  }
})

test_that("every par-template component name is spelled out by the collision message", {
  # thetaac, thetar and miss need an autocorrelation term, a residual
  # correlation and an imputed column to appear in a template at all;
  # the message that names them is unit-tested instead of fitting three
  # more models for one string each
  tpl <- list(beta = c(`z_(Intercept)` = 0), betad = numeric(0),
              b = numeric(0), theta = numeric(0), thetaac = numeric(0),
              thetar = numeric(0), miss = numeric(0))
  expect_match(frmtmb:::nl_start_collision_msg("beta", tpl),
               "fixed-effect coefficients")
  expect_match(frmtmb:::nl_start_collision_msg("betad", tpl),
               "distributional parameters")
  expect_match(frmtmb:::nl_start_collision_msg("b", tpl),
               "random-effect vector")
  expect_match(frmtmb:::nl_start_collision_msg("theta", tpl),
               "covariance parameters")
  expect_match(frmtmb:::nl_start_collision_msg("thetaac", tpl),
               "autocorrelation parameters")
  expect_match(frmtmb:::nl_start_collision_msg("thetar", tpl),
               "residual-correlation parameters")
  expect_match(frmtmb:::nl_start_collision_msg("miss", tpl),
               "imputed missing values")
})

test_that("newparams carries the same collision message, in its own spelling", {
  d <- nl_reserved_data()
  msg <- tryCatch(
    frm_simulate(bf(I ~ b - chi * logw, b ~ 1 + (1 | id), chi ~ 1 + group,
                    nl = TRUE),
                 data = d, family = exponential(link = "log"),
                 newparams = list(beta = c(1, 1.5, 0.2), b = 1,
                                  theta = 0.3),
                 nsim = 1, seed = 1),
    error = function(e) conditionMessage(e))
  expect_match(msg, "not the nonlinear parameter")
  # the message answers the argument the caller actually used
  expect_match(msg, "`newparams[$]b`")
  expect_false(grepl("start$", msg, fixed = TRUE))
})

test_that("the collision message keeps saying start$ for start", {
  d <- nl_reserved_data()
  msg <- tryCatch(
    frm(bf(I ~ b - chi * logw, b ~ 1 + (1 | id), chi ~ 1 + group,
           nl = TRUE),
        family = exponential(link = "log"), data = d, start = list(b = 1)),
    error = function(e) conditionMessage(e))
  expect_match(msg, "`start[$]b`")
  expect_false(grepl("newparams", msg, fixed = TRUE))
})

# Coefficients of the nonlinear parameters with a direction that
# changes no fitted value have a flat likelihood ridge: the update of
# ledger row brmsfit-methods:955 (bf(count ~ a + b) with fit2's
# a ~ Age + ..., b ~ Age + ...) is exactly flat along a_Intercept + t,
# b_Intercept - t (dev/fixes-u955.R). frm() fitted such a model to an
# arbitrary point on the ridge without a word. It now warns at the
# fitted point, from the Hessian block of those coefficients there
# (nl_flat_message()). Two checks before the fit refused identified
# models (a symbolic one on pnorm() bodies, a Jacobian at fixed points
# on natural-unit bodies; punch rounds 1 and 2), which is why the check
# moved to the fitted point and warns rather than refuses.
nl_sum_data <- local({
  set.seed(955)
  n <- 60
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 3 + 0.5 * d$x + 0.3 * d$z + stats::rnorm(n, 0, 0.4)
  d
})

# the warnings a fit raises and the error it stops with, if any
nl_run <- function(expr) {
  w <- character()
  err <- NA_character_
  tryCatch(suppressMessages(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })), error = function(e) err <<- conditionMessage(e))
  list(w = w, err = err)
}

nl_flat <- "are not identified: at the optimum the likelihood is flat"

nl_flat_w <- function(r) {
  w <- r$w[grepl(nl_flat, r$w, fixed = TRUE)]
  if (length(w)) w[1L] else NA_character_
}

test_that("a flat combination of nonlinear coefficients warns", {
  r <- nl_run(frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE),
                  data = nl_sum_data))
  expect_identical(r$err, NA_character_)
  m <- nl_flat_w(r)
  for (cf in c("a_Intercept", "a_x", "b_Intercept", "b_x")) {
    expect_match(m, cf, fixed = TRUE)
  }
  # a body that is a function of the sum, a difference, a weighted sum
  # and a product of two intercepts
  for (fo in list(bf(y ~ exp(a + b) * z, a ~ 1, b ~ 1 + x, nl = TRUE),
                  bf(y ~ a - b, a ~ 1 + x, b ~ 1, nl = TRUE),
                  bf(y ~ a + 2 * b, a ~ 1 + x, b ~ 1, nl = TRUE),
                  bf(y ~ a * b, a ~ 1, b ~ 1, nl = TRUE))) {
    st <- switch(deparse1(fo$formula), "y ~ a * b" = c(1, 1),
                 "y ~ exp(a + b) * z" = c(0, 0, 0), c(3, 0, 0.5))
    m <- nl_flat_w(nl_run(frm(fo, data = nl_sum_data,
                              start = list(beta = st))))
    expect_match(m, "a_Intercept", fixed = TRUE,
                 label = deparse1(fo$formula))
  }
  # only the coefficients that move are named
  m <- nl_flat_w(nl_run(frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + x + z,
                               nl = TRUE), data = nl_sum_data)))
  expect_match(m, "The coefficients a_x, b_x of", fixed = TRUE)
  # a prior on one coefficient leaves the other flat direction open
  m <- nl_flat_w(nl_run(frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE),
                            data = nl_sum_data,
                            prior = set_prior("normal(0, 1)", nlpar = "b",
                                              coef = "Intercept"))))
  expect_match(m, "a_x", fixed = TRUE)
  expect_false(grepl("a_Intercept", m, fixed = TRUE))
  # the update of brmsfit-methods:955's shape, on Gamma("identity"),
  # here from the parent fit's estimates: it reaches the ridge, and
  # warns
  dg <- nl_sum_data
  dg$y <- exp(dg$y / 3)
  f0 <- frm(bf(y ~ 1 / (1 + exp(-a)) * exp(b * z), a ~ 1 + x, b ~ 1 + x,
               nl = TRUE), data = dg, family = Gamma("identity"),
            start = list(beta = c(2, 0, 0.5, 0)))
  r <- nl_run(update(f0, formula. = bf(y ~ a + b, nl = TRUE)))
  expect_identical(r$err, NA_character_)
  expect_match(nl_flat_w(r), "a_Intercept, a_x, b_Intercept, b_x",
               fixed = TRUE)
})

test_that("the flat check names coefficients or stays silent", {
  # an exact ridge a + b beside c0 ~ 1 (lane fixes' final review, item
  # 7's control): with the noise of the differenced Hessian in the
  # naming threshold (lane nanse's RB3) no coefficient cleared it, and
  # the warning read "The coefficients  of the nonlinear parameter ''"
  # (the 0.68.0 merge). The check now stays silent there, and the
  # standard-error check names the parameters
  set.seed(45)
  d5 <- data.frame(x = rnorm(200))
  d5$y <- 1 + 0.5 * d5$x + rnorm(200, 0, 0.3)
  r <- nl_run(frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x, c0 ~ 1,
                     nl = TRUE), data = d5,
                  start = list(beta = c(0.5, 0.25, 0.5, 0.25, 1))))
  expect_identical(r$err, NA_character_)
  expect_false(any(grepl("The coefficients  of", r$w, fixed = TRUE)))
  named <- grepl("are not identified: at the optimum", r$w, fixed = TRUE) |
    grepl("Standard errors are not available", r$w, fixed = TRUE)
  expect_true(any(named))
  expect_true(any(grepl("a_x", r$w[named], fixed = TRUE)))
})

test_that("refit() does not repeat the flat-direction warning", {
  # documented in ?refit at 0.68.0 (lane fixes' consolidation item 2):
  # the design and its flat direction are the original fit's, which
  # warned, so a parametric bootstrap does not warn once per replicate
  r <- nl_run(f0 <- frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE),
                        data = nl_sum_data))
  expect_false(is.na(nl_flat_w(r)))
  set.seed(3)
  ys <- simulate(f0, nsim = 1)[[1]]
  r2 <- nl_run(refit(f0, ys))
  expect_identical(r2$err, NA_character_)
  expect_true(is.na(nl_flat_w(r2)))
})

test_that("an identified nonlinear model fits without that warning", {
  no_flat <- function(r, lab) {
    expect_identical(r$err, NA_character_, label = lab)
    expect_false(any(grepl(nl_flat, r$w, fixed = TRUE)), label = lab)
  }
  # disjoint terms: the joint design has full rank, and the fit is lm's
  f1 <- frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + z, nl = TRUE),
            data = nl_sum_data)
  ref <- stats::lm(y ~ x + z, data = nl_sum_data)
  # within a small fraction of lm's own standard errors
  se_ref <- unname(sqrt(diag(stats::vcov(ref))))
  expect_lt(max(abs(unname(fixef(f1)[c("a_Intercept", "a_x", "b_z"),
                                     "Estimate"]) - unname(coef(ref))) /
                  se_ref), 1e-3)
  no_flat(nl_run(frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                     data = nl_sum_data, start = list(beta = c(3, 0)))),
          "a * exp(b * x)")
  no_flat(nl_run(frm(bf(y ~ a + exp(b), a ~ 1, b ~ 1 + x, nl = TRUE),
                     data = nl_sum_data,
                     start = list(beta = c(1, 0.5, 0.1)))),
          "a + exp(b), b ~ 1 + x")
  # a prior on every coefficient of one parameter identifies the sum
  no_flat(nl_run(frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE),
                     data = nl_sum_data,
                     prior = set_prior("normal(0, 1)", nlpar = "b"))),
          "a + b with a prior on b")
  # a parameter that a second body also reads breaks the tie
  no_flat(nl_run(frm(bf(y ~ a + b, a ~ 1, b ~ 1, nl = TRUE) +
                       nlf(sigma ~ c0 + 0.1 * a, c0 ~ 1), data = nl_sum_data,
                     start = list(beta = c(3, 0)))),
          "a + b, a in sigma's body")
})

# The two identified models the fixed-point Jacobian check refused
# (punch round 2, B2): in their natural units the body saturates at
# coefficients near 0, a column underflowed to 0 and read as flat. Both
# fit on 0.67.0 with finite standard errors
# (dev/fixes-rev2-log/nl2-rellib-r5.txt).
test_that("a logistic growth curve over calendar years fits", {
  set.seed(22)
  yr <- seq(1900, 2000, length.out = 120)
  d <- data.frame(yr = yr, y = 50 / (1 + exp((1950 - yr) / 12)) +
                    stats::rnorm(120, 0, 1.5))
  r <- nl_run(fit <- frm(bf(y ~ Asym / (1 + exp((xmid - yr) / exp(lscal))),
                            Asym ~ 1, xmid ~ 1, lscal ~ 1, nl = TRUE),
                         data = d, start = list(beta = c(50, 1950, log(12)))))
  expect_identical(r$err, NA_character_)
  expect_false(any(grepl(nl_flat, r$w, fixed = TRUE)))
  # nls on SSlogis is the same model; the two optimizers agree within a
  # small fraction of frmtmb's own standard errors
  ref <- stats::nls(y ~ SSlogis(yr, Asym, xmid, scal), data = d)
  fe <- fixef(fit)
  est <- c(fe["Asym_Intercept", "Estimate"], fe["xmid_Intercept", "Estimate"],
           exp(fe["lscal_Intercept", "Estimate"]))
  se <- c(fe["Asym_Intercept", "Est.Error"], fe["xmid_Intercept", "Est.Error"],
          exp(fe["lscal_Intercept", "Estimate"]) *
            fe["lscal_Intercept", "Est.Error"])
  expect_true(all(is.finite(se)))
  expect_lt(max(abs(est - unname(stats::coef(ref))) / se), 1e-2)
})

test_that("a psychometric function with a lapse rate fits, x 200 to 400", {
  set.seed(21)
  x <- stats::runif(800, 200, 400)
  d <- data.frame(x = x, y = stats::rbinom(800, 1, 0.04 * 0.5 + 0.96 *
                                             stats::pnorm(x, 300, 25)))
  r <- nl_run(fit <- frm(bf(y ~ lapse * 0.5 + (1 - lapse) *
                              pnorm(x, m0, exp(ls)),
                            lapse ~ 1, m0 ~ 1, ls ~ 1, nl = TRUE),
                         family = bernoulli(link = "identity"), data = d,
                         start = list(beta = c(0.05, 300, log(25)))))
  expect_identical(r$err, NA_character_)
  expect_false(any(grepl(nl_flat, r$w, fixed = TRUE)))
  expect_true(all(is.finite(fixef(fit)[, "Est.Error"])))
})

# The three identified models the symbolic check refused (punch round
# 1, B1): stats::D() differentiates pnorm() and dnorm() in their first
# argument only, so a mean and an sd had "identical" derivatives, 0.
# Each fits on 0.67.0 with finite standard errors
# (dev/fixes-rev-log/nl-false-base.txt).
test_that("a response-preparation model in pnorm()'s mean and sd fits", {
  set.seed(7)
  n <- 3000
  pt <- stats::runif(n, 0, 0.6)
  f1 <- stats::pnorm(pt, 0.20, 0.045)
  f2 <- stats::pnorm(pt, 0.32, 0.045)
  p <- f2 * 0.95 + f1 * (1 - f2) * 0.2 + (1 - f1) * (1 - f2) * 0.5
  d1 <- data.frame(pt = pt, y = stats::rbinom(n, 1, p))
  prep <- frm(
    bf(y ~ pnorm(pt, m2, exp(ls)) * 0.95 +
         pnorm(pt, m1, exp(ls)) * (1 - pnorm(pt, m2, exp(ls))) * 0.2 +
         (1 - pnorm(pt, m1, exp(ls))) * (1 - pnorm(pt, m2, exp(ls))) * 0.5,
       m1 ~ 1, m2 ~ 1, ls ~ 1, nl = TRUE),
    family = bernoulli(link = "identity"), data = d1,
    start = list(beta = c(0.15, 0.4, -3)))
  expect_true(all(is.finite(fixef(prep)[, "Est.Error"])))
})

test_that("a psychometric pnorm(x, m0, exp(ls)) fits", {
  set.seed(11)
  x <- stats::runif(800, -2, 2)
  d2 <- data.frame(x = x, y = stats::rbinom(800, 1, stats::pnorm(x, 0.3, 0.7)))
  psy <- frm(bf(y ~ pnorm(x, m0, exp(ls)), m0 ~ 1, ls ~ 1, nl = TRUE),
             family = bernoulli(link = "identity"), data = d2,
             start = list(beta = c(0, 0)))
  # the same model written in the first argument is the same fit
  psy2 <- frm(bf(y ~ pnorm((x - m0) / exp(ls)), m0 ~ 1, ls ~ 1, nl = TRUE),
              family = bernoulli(link = "identity"), data = d2,
              start = list(beta = c(0, 0)))
  expect_lt(abs(as.numeric(logLik(psy) - logLik(psy2))),
            1e-8 * abs(as.numeric(logLik(psy2))))
})

test_that("a Gaussian bump h * dnorm(x, c0, exp(ls)) fits", {
  set.seed(12)
  x <- stats::runif(300, -3, 3)
  d3 <- data.frame(x = x, y = 2 * stats::dnorm(x, 0.5, 1.2) +
                     stats::rnorm(300, 0, 0.05))
  bump <- frm(bf(y ~ h * dnorm(x, c0, exp(ls)), h ~ 1, c0 ~ 1, ls ~ 1,
                 nl = TRUE), data = d3, start = list(beta = c(1, 0, 0)))
  expect_true(all(is.finite(fixef(bump)[, "Est.Error"])))
})

test_that("a parameter carried by its random effects alone is not flat", {
  # zz ~ 0 + (1 | g) is zero without its random effects, so with them
  # held at zero exp(lsd) * zz looked flat in lsd; the check refused
  # drmTMB's agreement model on the first Jacobian build of punch round
  # 1. It evaluates the random effects at nonzero values now
  set.seed(505)
  ng <- 40
  g <- factor(rep(seq_len(ng), each = 10))
  wg <- stats::rnorm(ng)
  d <- data.frame(g = g, x = stats::rnorm(400), w = wg[g])
  d$y <- 1 + 0.5 * d$x + stats::rnorm(ng, 0, exp(-0.5 + 0.5 * wg))[g] +
    stats::rnorm(400, 0, 0.7)
  fit <- frm(bf(y ~ b0 + exp(lsd) * zz, b0 ~ x, lsd ~ 0 + w,
                zz ~ 0 + (1 | g), nl = TRUE), family = gaussian(),
             data = d)
  expect_identical(fit$opt$convergence, 0L)
})
