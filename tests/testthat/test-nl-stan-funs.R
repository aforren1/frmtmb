# Stan's functions in a nonlinear body.
#
# brms writes a body into its Stan program, and its nonlinear vignette calls
# inv_logit() there. Through 0.78.0 frm() stopped with "could not find
# function". The set is brms's names_stan_functions() less the names Stan
# 2.32 does not parse with real arguments (dev/nlfix-stanc.R). The gated
# test at the end holds the values against Stan's own.

stan_fun_grid <- function() {
  list(
    logit = list(c(0.001, 0.2, 0.5, 0.9, 0.999)),
    inv_logit = list(c(-700, -40, -1, 0, 1, 40, 700)),
    inv_cloglog = list(c(-40, -2, -0.1, 0, 0.5, 1.5, 3)),
    Phi = list(c(-20, -2, -0.5, 0, 0.5, 2, 7)),
    square = list(c(-3, -0.2, 0, 0.5, 7)),
    cbrt = list(c(-27, -0.5, 0.5, 8)),
    exp2 = list(c(-3, -0.5, 0, 1, 10)),
    pow = list(c(0.5, 2, 3, 10), c(-1, 0.5, 2, -0.3)),
    inv = list(c(-4, -0.5, 0.5, 2)),
    inv_sqrt = list(c(0.25, 1, 4, 9)),
    inv_square = list(c(-2, -0.5, 0.5, 4)),
    hypot = list(c(3, -1, 0.5, 0), c(4, 2, -0.5, 1)),
    log1m = list(c(-3, 0, 0.5, 0.999)),
    step = list(c(-2, -1e-10, 0, 1e-10, 5)),
    log1p_exp = list(c(-800, -40, -1, 0, 1, 40, 800)),
    log1m_exp = list(c(-800, -40, -1, -0.7, -0.1, -1e-10)),
    log_diff_exp = list(c(0, 1, -1, 5), c(-1, -2, -3, 4.999)),
    log_sum_exp = list(c(0, 1, -1, 800), c(-1, 800, -3, 700)),
    log_inv_logit = list(c(-800, -40, -1, 0, 1, 40, 800)),
    log1m_inv_logit = list(c(-800, -40, -1, 0, 1, 40, 800)))
}

test_that("every Stan name in a body is on the tape and on doubles alike", {
  grids <- stan_fun_grid()
  expect_setequal(names(grids), names(nl_stan_funs))
  for (nm in names(nl_stan_funs)) {
    f <- nl_stan_funs[[nm]]
    g <- grids[[nm]]
    dbl <- do.call(f, g)
    k <- length(g)
    n <- length(g[[1L]])
    tp <- RTMB::MakeTape(function(p) {
      if (k == 1L) f(p) else f(p[seq_len(n)], p[n + seq_len(n)])
    }, unlist(g))
    v <- as.numeric(tp(unlist(g)))
    dv <- as.numeric(dbl)
    expect_identical(is.finite(v), is.finite(dv),
                     label = paste(nm, "finite on the tape and on doubles"))
    # the two paths may round differently in the last bits, and no further
    expect_lt(max(c(0, abs(v - dv)[is.finite(dv)] /
                      pmax(abs(dv)[is.finite(dv)], .Machine$double.xmin))),
              1e-13, label = paste(nm, "tape vs doubles"))
    J <- tp$jacobian(unlist(g))
    expect_false(anyNA(J), label = paste(nm, "jacobian"))
  }
})

# the analytic derivative of each name in each argument, written from the
# closed form and independently of nl_stan_funs. log-space where the plain
# formula rounds at the grid's ends.
stan_fun_grads <- function() {
  list(
    logit = list(function(x) 1 / (x * (1 - x))),
    inv_logit = list(stats::dlogis),
    inv_cloglog = list(function(x) exp(x - exp(x))),
    Phi = list(stats::dnorm),
    square = list(function(x) 2 * x),
    cbrt = list(function(x) (1 / 3) * abs(x)^(-2 / 3)),
    exp2 = list(function(x) log(2) * 2^x),
    pow = list(function(x, y) y * x^(y - 1), function(x, y) x^y * log(x)),
    inv = list(function(x) -1 / x^2),
    inv_sqrt = list(function(x) -0.5 * x^(-1.5)),
    inv_square = list(function(x) -2 * x^(-3)),
    hypot = list(function(x, y) x / sqrt(x^2 + y^2),
                 function(x, y) y / sqrt(x^2 + y^2)),
    log1m = list(function(x) -1 / (1 - x)),
    step = list(function(x) 0 * x),
    log1p_exp = list(stats::plogis),
    log1m_exp = list(function(x) -1 / expm1(-x)),
    log_diff_exp = list(function(x, y) -1 / expm1(y - x),
                        function(x, y) -1 / expm1(x - y)),
    log_sum_exp = list(function(x, y) stats::plogis(x - y),
                       function(x, y) stats::plogis(y - x)),
    log_inv_logit = list(function(x) stats::plogis(-x)),
    log1m_inv_logit = list(function(x) -stats::plogis(x)))
}

test_that("every Stan name has the derivative of its closed form", {
  grids <- stan_fun_grid()
  grads <- stan_fun_grads()
  expect_setequal(names(grads), names(nl_stan_funs))
  for (nm in names(nl_stan_funs)) {
    f <- nl_stan_funs[[nm]]
    g <- grids[[nm]]
    k <- length(g)
    n <- length(g[[1L]])
    tp <- RTMB::MakeTape(function(p) {
      if (k == 1L) f(p) else f(p[seq_len(n)], p[n + seq_len(n)])
    }, unlist(g))
    J <- tp$jacobian(unlist(g))
    for (a in seq_len(k)) {
      ours <- diag(J[seq_len(n), (a - 1L) * n + seq_len(n), drop = FALSE])
      ref <- do.call(grads[[nm]][[a]], g)
      keep <- is.finite(ref)
      # a true derivative below 1e-280 may round to 0 in the product that
      # forms it, on the tape and in Stan's own inv_logit
      expect_lt(max(c(0, abs(ours[keep] - ref[keep]) /
                        pmax(abs(ref[keep]), 1e-280))), 1e-9,
                label = paste(nm, "argument", a))
    }
  }
})

test_that("the ends the plain spellings lost keep their value and gradient", {
  grad <- function(f, x) diag(RTMB::MakeTape(function(p) f(p), x)$jacobian(x))
  # 5e-63 at x = 5, where expm1's own derivative rounds to 0
  expect_lt(abs(grad(nl_stan_funs$inv_cloglog, 5) - exp(5 - exp(5))) /
              exp(5 - exp(5)), 1e-9)
  # no overflow and no underflow in the squares
  expect_lt(abs(nl_stan_funs$hypot(1e200, 1e200) / (sqrt(2) * 1e200) - 1),
            1e-14)
  expect_lt(abs(nl_stan_funs$hypot(1e-200, 1e-200) / (sqrt(2) * 1e-200) - 1),
            1e-14)
  expect_identical(nl_stan_funs$hypot(0, 0), 0)
  # Stan's -Inf at 0, not 0 / 0
  expect_identical(grad(nl_stan_funs$inv_square, 0), -Inf)
})

test_that("a column and a call of one name in one body", {
  set.seed(4006)
  d <- data.frame(inv_logit = rnorm(120))
  d$y <- stats::plogis(0.8 * d$inv_logit) + rnorm(120, 0, 0.1)
  fit <- frm(bf(y ~ inv_logit(a * inv_logit), a ~ 1, nl = TRUE), data = d)
  expect_gt(fit$estimates[["beta"]][[1L]], 0.3)
})

test_that("the gradient of a Stan name is the one the formula has", {
  # the ends are where a plain spelling loses it
  grad <- function(f, x) {
    diag(RTMB::MakeTape(function(p) f(p), x)$jacobian(x))
  }
  x <- c(-800, -40, -1, 0, 1, 40, 800)
  expect_equal(grad(nl_stan_funs$inv_logit, x), stats::dlogis(x),
               tolerance = 1e-10)
  expect_equal(grad(nl_stan_funs$log1p_exp, x), stats::plogis(x),
               tolerance = 1e-10)
  expect_equal(grad(nl_stan_funs$log_inv_logit, x), stats::plogis(-x),
               tolerance = 1e-10)
  expect_equal(grad(nl_stan_funs$log1m_inv_logit, x), -stats::plogis(x),
               tolerance = 1e-10)
  expect_identical(grad(nl_stan_funs$step, c(-2, 0, 5)), c(0, 0, 0))
})

test_that("a body can call inv_logit(), the vignette's item response model", {
  set.seed(4001)
  d <- data.frame(ability = rnorm(300))
  d$y <- rbinom(300, 1, 0.33 + 0.67 * plogis(1.2 * d$ability))
  fit <- frm(bf(y ~ 0.33 + 0.67 * inv_logit(eta), eta ~ ability, nl = TRUE),
             family = bernoulli("identity"), data = d)
  # the same model with the logistic spelled out
  ref <- frm(bf(y ~ 0.33 + 0.67 / (1 + exp(-eta)), eta ~ ability, nl = TRUE),
             family = bernoulli("identity"), data = d)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))) /
              abs(as.numeric(logLik(ref))), 1e-9)
  # the post-fit paths evaluate the same body
  mu <- as.numeric(frm_linpred(fit))
  b <- fit$estimates[["beta"]]
  by_hand <- 0.33 + 0.67 * stats::plogis(b[[1L]] + b[[2L]] * d$ability)
  expect_lt(max(abs(mu - by_hand)) / max(abs(by_hand)), 1e-12)
  nd <- data.frame(ability = c(-1, 0, 2))
  expect_length(as.numeric(fitted(fit, newdata = nd)[, "Estimate"]), 3L)
})

test_that("the other names resolve in a body too", {
  set.seed(4003)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + 0.5 * d$x + rnorm(200, 0, 0.5)
  d$p <- stats::plogis(d$x)
  pairs <- list(
    list(y ~ b0 + b1 * inv_cloglog(x), y ~ b0 + b1 * (1 - exp(-exp(x)))),
    list(y ~ b0 + b1 * Phi(x), y ~ b0 + b1 * pnorm(x)),
    list(y ~ b0 + b1 * square(x), y ~ b0 + b1 * x^2),
    list(y ~ b0 + b1 * cbrt(x), y ~ b0 + b1 * sign(x) * abs(x)^(1 / 3)),
    list(y ~ b0 + b1 * exp2(x), y ~ b0 + b1 * 2^x),
    list(y ~ b0 + b1 * pow(p, 2), y ~ b0 + b1 * p^2),
    list(y ~ b0 + b1 * logit(p), y ~ b0 + b1 * x),
    list(y ~ b0 + b1 * log1m(p), y ~ b0 + b1 * log(1 - p)),
    list(y ~ b0 + b1 * inv(p), y ~ b0 + b1 / p),
    list(y ~ b0 + b1 * inv_sqrt(p), y ~ b0 + b1 / sqrt(p)),
    list(y ~ b0 + b1 * inv_square(p), y ~ b0 + b1 / p^2),
    list(y ~ b0 + b1 * hypot(x, p), y ~ b0 + b1 * sqrt(x^2 + p^2)),
    list(y ~ b0 + b1 * log1p_exp(x), y ~ b0 + b1 * log(1 + exp(x))),
    list(y ~ b0 + b1 * log1m_exp(-p), y ~ b0 + b1 * log(1 - exp(-p))),
    list(y ~ b0 + b1 * log_inv_logit(x), y ~ b0 + b1 * log(stats::plogis(x))),
    list(y ~ b0 + b1 * log1m_inv_logit(x),
         y ~ b0 + b1 * log(1 - stats::plogis(x))),
    list(y ~ b0 + b1 * log_sum_exp(x, p), y ~ b0 + b1 * log(exp(x) + exp(p))),
    list(y ~ b0 + b1 * log_diff_exp(x + 1, p - 2),
         y ~ b0 + b1 * log(exp(x + 1) - exp(p - 2))))
  for (pr in pairs) {
    st <- frm(bf(pr[[1L]], b0 + b1 ~ 1, nl = TRUE), family = gaussian(),
              data = d)
    rf <- frm(bf(pr[[2L]], b0 + b1 ~ 1, nl = TRUE), family = gaussian(),
              data = d)
    expect_lt(abs(as.numeric(logLik(st)) - as.numeric(logLik(rf))) /
                abs(as.numeric(logLik(rf))), 1e-8,
              label = deparse1(pr[[1L]]))
  }
})

test_that("a helper the user wrote under a Stan name is found first", {
  set.seed(4004)
  d <- data.frame(x = rnorm(120))
  d$y <- 1 + 0.5 * d$x + rnorm(120, 0, 0.5)
  # a user's inv_logit() is not Stan's, and is the one the body means
  inv_logit <- function(v) 2 * v
  f <- as.formula("y ~ b0 + inv_logit(b1) * x")
  environment(f) <- environment()
  fit <- frm(bf(f, b0 + b1 ~ 1, nl = TRUE), family = gaussian(), data = d)
  expect_equal(unname(fit$estimates[["beta"]])[2L] * 2,
               unname(coef(lm(y ~ x, d))[2L]), tolerance = 1e-5)
  # stats::step() on the search path is not a helper, and a data column
  # called `step` stays a column
  fit2 <- frm(bf(y ~ b0 + b1 * step(x), b0 + b1 ~ 1, nl = TRUE),
              family = gaussian(), data = d)
  ref <- lm(y ~ I(x >= 0), d)
  expect_lt(abs(as.numeric(logLik(fit2)) - as.numeric(logLik(ref))) /
              abs(as.numeric(logLik(ref))), 1e-8)
  d$step <- d$x
  fit3 <- frm(bf(y ~ b0 + b1 * step, b0 + b1 ~ 1, nl = TRUE),
              family = gaussian(), data = d)
  expect_lt(abs(as.numeric(logLik(fit3)) -
                  as.numeric(logLik(lm(y ~ x, d)))) /
              abs(as.numeric(logLik(fit3))), 1e-8)
})

test_that("a name that is not Stan's or the user's keeps its error", {
  set.seed(4005)
  d <- data.frame(x = rnorm(40))
  d$y <- rnorm(40)
  expect_error(frm(bf(y ~ b0 + nosuchfun(x) * b1, b0 + b1 ~ 1, nl = TRUE),
                   family = gaussian(), data = d),
               "nosuchfun", fixed = TRUE)
})

test_that("the values are Stan's", {
  skip_unless_brms_fit()
  grids <- stan_fun_grid()
  # the whole grid in one program, computed by Stan's math library in
  # generated quantities, row by row as a brms body is
  nms <- names(grids)
  one <- function(nm, i) {
    g <- grids[[nm]]
    args <- if (length(g) == 1L) {
      sprintf("x_%s[%d]", nm, i)
    } else {
      sprintf("x_%s[%d], y_%s[%d]", nm, i, nm, i)
    }
    sprintf("  r_%s[%d] = %s(%s);", nm, i, nm, args)
  }
  data_block <- unlist(lapply(nms, function(nm) {
    g <- grids[[nm]]
    n <- length(g[[1L]])
    c(sprintf("  vector[%d] x_%s;", n, nm),
      if (length(g) == 2L) sprintf("  vector[%d] y_%s;", n, nm))
  }))
  gq_decl <- vapply(nms, function(nm) {
    sprintf("  vector[%d] r_%s;", length(grids[[nm]][[1L]]), nm)
  }, "")
  gq <- unlist(lapply(nms, function(nm) {
    vapply(seq_along(grids[[nm]][[1L]]), function(i) one(nm, i), "")
  }))
  code <- paste(c("data {", data_block, "}", "generated quantities {", gq_decl,
                  gq, "}"), collapse = "\n")
  sdat <- list()
  for (nm in nms) {
    g <- grids[[nm]]
    sdat[[paste0("x_", nm)]] <- as.array(g[[1L]])
    if (length(g) == 2L) sdat[[paste0("y_", nm)]] <- as.array(g[[2L]])
  }
  sf <- suppressMessages(suppressWarnings(rstan::sampling(
    brms_stan_model(code), data = sdat, algorithm = "Fixed_param",
    chains = 1, iter = 1, warmup = 0, refresh = 0)))
  ex <- rstan::extract(sf)
  for (nm in nms) {
    ours <- as.numeric(do.call(nl_stan_funs[[nm]], grids[[nm]]))
    stan <- as.numeric(ex[[paste0("r_", nm)]])
    expect_identical(is.finite(ours), is.finite(stan), label = nm)
    fin <- is.finite(stan)
    # relative to the largest value in the grid's own scale, with Stan's
    # rounding of a value that underflows (a denormal against 0) excluded
    big <- abs(stan[fin]) > .Machine$double.xmin
    expect_lt(max(c(0, abs(ours[fin][big] - stan[fin][big]) /
                      abs(stan[fin][big]))), 1e-12, label = nm)
  }
})
