# Fits whose likelihood has no maximum at finite values, and the edge of
# a family parameter's range: the standard-error check's limit rule
# (R/se-check.R, se_edge_limit()). And a covariance interval of zero
# width (R/confint.R, varcorr_trans_rows()). Lane boundfam,
# dev/boundfam-findings.md.

limit_conds <- function(expr) {
  w <- character()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}

limit_se <- function(fit) {
  suppressWarnings(sqrt(diag(vcov(fit, full = TRUE))))
}

test_that("an all-zero poisson response has no maximum, and says so", {
  # 0.81.0 stopped at an intercept of -636.8 with standard errors of
  # 2.7e137 and 2.4e137, under only the optimizer's code 9, whose own
  # text said the stop was at the maximum. glm() reports 44231 in
  # silence. Every rounding of the covariate must end the same way
  set.seed(2)
  x <- stats::rnorm(50)
  for (k in 0:19) {
    d <- data.frame(x = x * (1 + 2e-15 * k), y = 0L)
    r <- limit_conds(frm(y ~ x, family = poisson(), data = d))
    expect_length(r$warnings, 1L)
    expect_match(r$warnings, "maximum likelihood estimate does not exist",
                 fixed = TRUE)
    expect_match(r$warnings, "(Intercept) (at ", fixed = TRUE)
    expect_match(r$warnings, "toward -Inf", fixed = TRUE)
    expect_false(any(grepl("Optimizer did not report", r$warnings,
                           fixed = TRUE)))
    lost <- sdr_of(r$value)$se_lost
    expect_identical(unname(lost), c("limit", "limit"))
    expect_true(all(is.nan(limit_se(r$value))))
  }
})

test_that("a response censored on every row has no maximum", {
  # seed 3: code 0 and no warning at all on 0.81.0, with standard errors
  # of 6.4e8, 1.3e8 and 2.1e8. survreg() runs out of iterations there
  for (s in c(1L, 3L)) {
    set.seed(s)
    d <- data.frame(x = stats::rnorm(60))
    d$y <- 1 + d$x + stats::rnorm(60)
    d$cen <- "right"
    r <- limit_conds(frm(y | cens(cen) ~ x, data = d))
    expect_length(r$warnings, 1L)
    expect_match(r$warnings, "maximum likelihood estimate does not exist",
                 fixed = TRUE)
    expect_true(all(is.nan(limit_se(r$value))))
  }
})

test_that("one poisson count where a maximum exists is untouched", {
  # the guard-absent case: one count at the median covariate value keeps
  # an interior maximum, glm()'s
  set.seed(2)
  d <- data.frame(x = stats::rnorm(50), y = 0L)
  d$y[order(d$x)[25]] <- 1L
  r <- limit_conds(frm(y ~ x, family = poisson(), data = d))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 0L)
  g <- stats::glm(y ~ x, family = stats::poisson(), data = d)
  se_glm <- sqrt(diag(stats::vcov(g)))
  expect_equal(unname(limit_se(r$value) / se_glm), c(1, 1),
               tolerance = 1e-3)
})

test_that("a negbinomial shape at its poisson limit is a boundary fit", {
  # Seed 2 of poisson data runs the log shape to 15.9, where 0.81.0 kept
  # a standard error of 1218 in silence and MASS::glm.nb() warns that
  # its iteration limit was reached. The coefficients' standard errors
  # are the poisson model's
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.4 * d$x))
  r <- limit_conds(frm(y ~ x, family = negbinomial(), data = d))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 1L)
  expect_s3_class(r$messages[[1L]], "frmtmb_boundary_fit")
  expect_match(conditionMessage(r$messages[[1L]]),
               "shape_(Intercept) is at the end of its range", fixed = TRUE)
  se <- limit_se(r$value)
  expect_true(is.nan(se[["shape_(Intercept)"]]))
  g <- stats::glm(y ~ x, family = stats::poisson(), data = d)
  expect_equal(unname(se[1:2] / sqrt(diag(stats::vcov(g)))), c(1, 1),
               tolerance = 1e-3)
  # check_se = "ignore" is silent, "stop" stops, as for a variance at zero
  expect_silent(frm(y ~ x, family = negbinomial(), data = d,
                    control = frmtmb_control(check_se = "ignore")))
  expect_error(frm(y ~ x, family = negbinomial(), data = d,
                   control = frmtmb_control(check_se = "stop")),
               "end of its range", fixed = TRUE)
})

test_that("a negbinomial shape with a maximum keeps its standard error", {
  # the guard-absent case: seed 1 of the same design has its maximum at a
  # log shape of 4.04 (glm.nb(): theta 56.7, standard error 100)
  set.seed(1)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.4 * d$x))
  r <- limit_conds(frm(y ~ x, family = negbinomial(), data = d))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 0L)
  expect_true(all(is.finite(limit_se(r$value))))
})

test_that("a mixture precision at the binomial limit is a boundary fit", {
  # test-bcm-latent-mixtures.R's Malingering_2: log phi2 at 23.4, with a
  # standard error of 1.3e4 and no word on 0.81.0
  d <- data.frame(k = c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30, 20, 6,
                        44, 44, 27, 25, 17, 14, 27, 35, 30), n = 45L)
  r <- limit_conds(frm(bf(k | trials(n) ~ 1),
                       family = mixture(beta_binomial, beta_binomial),
                       data = d,
                       start = list(beta = stats::qlogis(c(0.55, 0.95)),
                                    betad = c(log(c(15, 15)), 0))))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 1L)
  expect_s3_class(r$messages[[1L]], "frmtmb_boundary_fit")
  expect_match(conditionMessage(r$messages[[1L]]),
               "phi2_(Intercept) is at the end of its range", fixed = TRUE)
  se <- limit_se(r$value)
  expect_true(is.nan(se[["phi2_(Intercept)"]]))
  expect_true(all(is.finite(se[names(se) != "phi2_(Intercept)"])))
})

test_that("a family parameter held at its limit keeps fitted() SEs", {
  # backlog item 3269: zi's intercept runs to its limit and the fit is
  # the binomial model. fitted() moves by at most 2.7e-9 per logit unit
  # along zi, and its standard errors are the binomial model's (0.81.0
  # to 2.8e-10); with zi's direction kept as lost, every row was NaN
  set.seed(7)
  d <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
  d$y <- stats::rbinom(200, d$N, stats::plogis(0.3 + 0.5 * d$x))
  r <- limit_conds(frm(bf(y | trials(N) ~ x),
                       family = zero_inflated_binomial(), data = d))
  expect_identical(unname(sdr_of(r$value)$se_lost), "limit")
  fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = d)
  fe <- limit_conds(stats::fitted(r$value))
  expect_length(fe$warnings, 0L)
  eb <- stats::fitted(fb)[, "Est.Error"]
  expect_equal(fe$value[, "Est.Error"] / eb, rep(1, nrow(d)),
               tolerance = 1e-6)
  # the absent case: an undetermined mean coefficient keeps its lost
  # direction, so a fitted value that reads it has no standard error
  set.seed(2)
  de <- data.frame(x = stats::rnorm(50), y = 0L)
  f0 <- limit_conds(frm(y ~ x, family = poisson(), data = de))$value
  fp <- limit_conds(stats::fitted(f0, newdata = data.frame(x = c(-1, 1))))
  expect_true(all(is.nan(fp$value[, "Est.Error"])))
  expect_match(fp$warnings, "move along a direction the fit does not",
               fixed = TRUE)
})

test_that("a fitted value that barely moves along a limit keeps its SE", {
  # zi ~ x (or hu ~ x) at the probability's end: the intercept is at its
  # end and the slope undetermined there, and the rows move by about
  # 1e-9 per unit along them. The cosine test called 195 of 200 (297 of
  # 300) fitted values NaN; 0.81.0 gave finite values, which the lane
  # keeps to 1.7e-8 relative
  set.seed(1)
  d <- data.frame(x = stats::rnorm(200), N = 10L)
  p <- stats::plogis(-0.5 + 0.8 * d$x)
  d$y <- ifelse(stats::runif(200) < 0, 0L, stats::rbinom(200, d$N, p))
  f1 <- limit_conds(frm(bf(y | trials(N) ~ x, zi ~ x),
                        family = zero_inflated_binomial(), data = d))$value
  expect_identical(unname(sdr_of(f1)$se_lost), c("limit", "limit"))
  r1 <- limit_conds(stats::fitted(f1))
  expect_length(r1$warnings, 0L)
  fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = d)
  expect_equal(r1$value[, "Est.Error"] / stats::fitted(fb)[, "Est.Error"],
               rep(1, nrow(d)), tolerance = 1e-6)
  set.seed(12)
  dh <- data.frame(x = stats::rnorm(300))
  dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
  f4 <- limit_conds(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(),
                        data = dh))$value
  expect_identical(unname(sdr_of(f4)$se_lost), c("limit", "limit"))
  r4 <- limit_conds(stats::fitted(f4))
  expect_length(r4$warnings, 0L)
  expect_true(all(is.finite(r4$value[, "Est.Error"])))
  # the absent case: a separated bernoulli's fitted values read their
  # runaway direction, and stay NaN with the warning
  set.seed(514)
  ds <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
  ds$yb <- as.integer(ds$z > 0)
  fs <- limit_conds(frm(yb ~ z + x, family = bernoulli(), data = ds))$value
  rs <- limit_conds(stats::fitted(fs))
  expect_gt(sum(is.nan(rs$value[, "Est.Error"])), 0.9 * nrow(ds))
  expect_match(rs$warnings, "move along a direction the fit does not",
               fixed = TRUE)
})

test_that("a quantity reading a NaN variance is non-estimable on any route", {
  # Where the check keeps sdreport()'s inverse (a stop short, or a
  # separation that explained it) there is no null basis, and a lost
  # parameter's variance is NaN in V. Under OpenBLAS 0.3.26 the
  # separated bernoulli above takes that route (code 1), and its 240 NaN
  # fitted values came with no warning. jc_nonest() reads V itself
  V <- diag(c(1, NaN, 1))
  G <- rbind(c(1, 0, 1), c(1, 2, 0), c(0, 0, 0))
  jc <- list(V = V, null = NULL, lost_pos = 2L)
  expect_identical(jc_nonest(jc, G, 1:3), c(FALSE, TRUE, FALSE))
  # the absent case: a lost parameter held at zero variance (a family
  # limit, or the propagation covariance's zero row) is not one
  jc$V[2, 2] <- 0
  expect_identical(jc_nonest(jc, G, 1:3), c(FALSE, FALSE, FALSE))
  jc$lost_pos <- integer(0)
  jc$V[2, 2] <- NaN
  expect_identical(jc_nonest(jc, G, 1:3), c(FALSE, FALSE, FALSE))
})

test_that("a family limit on a nonlinear fit still triggers the restart", {
  skip_on_cran()
  # lane gpstop's sweep row pow_a_student_n120_s2: from the default start
  # a student() nu runs to exp(22.7) at a log-likelihood of -4.95. 0.81.0
  # called nu flat, and that warning sent frm() to restart, which
  # reaches a maximum at 55.67. The quiet limit message must do the same
  n <- 120L
  set.seed(81100 + 1000 * 4 + 100 * 2 + n + 10 * 6)
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  w <- round(stats::runif(n, 0.2, 5), 3)
  f <- 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  y <- (2 + fs) * w^0.5 + 0.1 * stats::rt(n, 3)
  d <- data.frame(x, z, w, y)
  r <- limit_conds(frm(bf(y ~ a * w^b, a ~ gp(x), b ~ 1, nl = TRUE),
                       data = d, family = student()))
  cls <- vapply(r$messages, function(m) class(m)[1L], "")
  expect_true("frmtmb_nl_restart" %in% cls)
  # the replaced fit's limit message is not released
  expect_false("frmtmb_limit_fit" %in% cls)
  expect_gt(as.numeric(stats::logLik(r$value)), 0)
})

test_that("a raw-unit column beside a smooth keeps its SE at the maximum", {
  skip_on_cran()
  # A smooth gives the fit par_units under autoscale = FALSE, and the
  # first limit rule then probed xs in raw units: at scale 1e-5 and 1e-6
  # a step of 2 read level both ways, and xs got NaN and the "does not
  # exist" warning at the maximum the autoscaled fit finds. The rule
  # now steps each coefficient in its natural unit
  for (k in list(c("poisson", 1e-6), c("bernoulli", 1e-5),
                 c("bernoulli", 1e-6))) {
    fam <- k[1]
    sc <- as.numeric(k[2])
    set.seed(9)
    n <- 300
    d <- data.frame(x = stats::rnorm(n), z = stats::runif(n))
    eta <- 0.5 * (fam == "poisson") + 0.3 * d$x + sin(2 * pi * d$z)
    d$y <- if (fam == "poisson") stats::rpois(n, exp(eta)) else
      stats::rbinom(n, 1, stats::plogis(eta))
    d$xs <- d$x * sc
    famf <- if (fam == "poisson") poisson() else bernoulli()
    f1 <- suppressWarnings(suppressMessages(
      frm(y ~ xs + s(z), family = famf, data = d)))
    st <- list(beta = f1$estimates$beta, theta = f1$estimates$theta)
    r <- limit_conds(frm(y ~ xs + s(z), family = famf, data = d,
                         start = st,
                         control = frmtmb_control(autoscale = FALSE)))
    expect_length(r$warnings, 0L)
    expect_length(r$messages, 0L)
    expect_true(all(is.finite(limit_se(r$value))))
    expect_false("limit" %in% sdr_of(r$value)$se_lost)
  }
})

test_that("a zero-inflation or hurdle slope that separates warns", {
  # zi (or hu) is 1 below x = 0 and 0 above: a separation of that
  # probability's model, not the family's limit. The first limit rule
  # called it quiet, removed the code-7 warning, and left the intercept
  # a standard error of 1063
  for (s in 1:2) {
    set.seed(s)
    d <- data.frame(x = stats::rnorm(200))
    d$y <- ifelse(d$x < 0, 0L, stats::rpois(200, exp(2)))
    r <- limit_conds(frm(bf(y ~ 1, zi ~ x),
                         family = zero_inflated_poisson(), data = d))
    expect_length(r$messages, 0L)
    expect_true(any(grepl("separate the rows of its probability",
                          r$warnings, fixed = TRUE)))
    se <- limit_se(r$value)
    expect_true(all(is.nan(se[c("zi_(Intercept)", "zi_x")])))
    expect_true(is.finite(se[["(Intercept)"]]))
    set.seed(s)
    d <- data.frame(x = stats::rnorm(200))
    d$y <- ifelse(d$x < 0, 0L, 1L + stats::rpois(200, 2))
    r <- limit_conds(frm(bf(y ~ 1, hu ~ x), family = hurdle_poisson(),
                         data = d))
    expect_length(r$messages, 0L)
    expect_true(any(grepl("separate the rows of its probability",
                          r$warnings, fixed = TRUE)))
    se <- limit_se(r$value)
    expect_true(all(is.nan(se[c("hu_(Intercept)", "hu_x")])))
  }
})

test_that("a probability intercept or one group at its limit stays quiet", {
  skip_if_not_installed("MASS")
  # the controls of the separation rule: zi's intercept run to zero on
  # poisson data (seed 2: logit -19.6), and a negbinomial shape ~ g whose
  # group b is poisson
  set.seed(2)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rpois(300, exp(1 + 0.3 * d$x))
  r <- limit_conds(frm(bf(y ~ x, zi ~ 1), family = zero_inflated_poisson(),
                       data = d))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 1L)
  expect_s3_class(r$messages[[1L]], "frmtmb_boundary_fit")
  expect_true(is.nan(limit_se(r$value)[["zi_(Intercept)"]]))
  set.seed(2)
  d <- data.frame(g = factor(rep(c("a", "b"), each = 200)))
  d$y <- ifelse(d$g == "a", MASS::rnegbin(400, 4, 2), stats::rpois(400, 4))
  r <- limit_conds(frm(bf(y ~ 1, shape ~ g), family = negbinomial(),
                       data = d))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 1L)
  expect_match(conditionMessage(r$messages[[1L]]),
               "shape_gb is at the end of its range", fixed = TRUE)
  se <- limit_se(r$value)
  expect_true(is.nan(se[["shape_gb"]]))
  expect_true(is.finite(se[["shape_(Intercept)"]]))
})

test_that("a zero standard error gives no interval", {
  # The delta method gives a zero standard error where the jacobian row
  # is flat or the quadratic form rounds to zero or below. 0.72.0 to
  # 0.81.0 reported such a correlation as an interval of zero width,
  # [-0.122416, -0.122416] on the 0.72.0 review's seed 7. Here the
  # covariance of the theta block is set to zero, so every row's
  # standard error is zero
  set.seed(62)
  n_g <- 30
  g <- factor(rep(seq_len(n_g), each = 8))
  U <- matrix(stats::rnorm(n_g * 2), n_g) %*%
    chol(matrix(c(0.64, 0.2, 0.2, 0.25), 2))
  dd <- data.frame(g = g, x = stats::rnorm(n_g * 8))
  dd$y <- 1 + 0.5 * dd$x + U[g, 1] + U[g, 2] * dd$x +
    stats::rnorm(n_g * 8, 0, 0.6)
  fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
  real <- sdr_of
  zero_theta <- function(fit, ...) {
    s <- real(fit, ...)
    th <- rownames(s$cov.fixed) == "theta"
    s$cov.fixed[th, ] <- 0
    s$cov.fixed[, th] <- 0
    if (!is.null(s$cov_fixed_prop)) {
      s$cov_fixed_prop[th, ] <- 0
      s$cov_fixed_prop[, th] <- 0
    }
    s
  }
  # the absent case first: the real covariance gives a positive se
  tr <- varcorr_trans_rows(fit)
  expect_true(all(tr$se_t > 0))
  local_mocked_bindings(sdr_of = zero_theta)
  tr <- varcorr_trans_rows(fit)
  expect_true(all(is.na(tr$se_t)))
  ci <- allow_warnings(confint_varcorr(fit), "No interval for",
                       require = "No interval for")
  expect_true(all(is.na(ci$lwr) & is.na(ci$upr)))
})

test_that("a us() row combining the rows before it has no interval", {
  # the multiple correlation of a term with the terms before it, from
  # the unit-diagonal Cholesky row, against the clamp; in dimension 2
  # it is the one correlation's own test
  tc <- varcorr_cor_clamp / sqrt(1 - varcorr_cor_clamp^2)
  th <- c(0, 0, 0, 0.3, 0.6 * tc, 0.81 * tc)
  expect_identical(varcorr_us_row_edge(th, 3L), c(FALSE, TRUE, TRUE))
  th[5:6] <- 0.6 * tc
  expect_identical(varcorr_us_row_edge(th, 3L), c(FALSE, FALSE, FALSE))
  r2 <- c(0, 0, tc * (1 + 1e-9))
  expect_identical(varcorr_us_row_edge(r2, 2L),
                   abs(r2[3] / sqrt(1 + r2[3]^2)) >= varcorr_cor_clamp)
})

test_that("confint_varcorr() gives a singular three-term block no interval", {
  skip_on_cran()
  # The z slope is a combination of the intercept and the x slope, so
  # the block has rank 2 and the third Cholesky row runs off with no
  # correlation near +/-1: seed 6 stops at (1771, -2271), where 0.81.0
  # gave cor(z,(Intercept)) a Fisher-z standard error of 9.2e-4. Over
  # ulp-moved responses the row's correlations get NA bounds wherever
  # the row is past the clamp, and the first row keeps its interval
  set.seed(6)
  ng <- 30
  d <- data.frame(x = stats::rnorm(ng * 10), z = stats::rnorm(ng * 10),
                  g = factor(rep(seq_len(ng), each = 10)))
  u0 <- stats::rnorm(ng, 0, 0.8)
  u1 <- stats::rnorm(ng, 0, 0.5)
  u2 <- 0.4 * u0 - 0.6 * u1
  y0 <- 1 + 0.5 * d$x - 0.3 * d$z + u0[d$g] + u1[d$g] * d$x +
    u2[d$g] * d$z + stats::rnorm(ng * 10, 0, 0.5)
  past <- 0L
  tc <- varcorr_cor_clamp / sqrt(1 - varcorr_cor_clamp^2)
  for (k in 0:5) {
    d$y <- y0 * (1 + 2e-15 * k)
    fit <- suppressWarnings(frm(y ~ x + z + (1 + x + z | g), data = d))
    th <- fit$estimates$theta
    # rows 2 and 3 of the unit-diagonal Cholesky factor: theta 4, and
    # theta 5 and 6; correlations in lower.tri() order (2,1), (3,1), (3,2)
    r2 <- th[4]^2 >= tc^2
    r3 <- th[5]^2 + th[6]^2 >= tc^2
    edge <- c(r2, r3, r3)
    tr <- suppressWarnings(varcorr_trans_rows(fit))
    se_cor <- tr$se_t[tr$type == "cor"]
    expect_true(all(is.na(se_cor[edge])))
    expect_true(all(is.na(se_cor) | se_cor > 0))
    past <- past + any(edge)
  }
  expect_gt(past, 0L)
})
