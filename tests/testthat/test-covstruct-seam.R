# frmtmb_register_covstruct(), the public covariance-structure seam.
#
# Each structure here is registered through the exported function alone,
# the way another package's .onLoad() would, fitted, and checked against
# a computation that does not go through the registry: the dense
# marginal likelihood of a gaussian model, and core's own gr(cov = A) for
# the same field. The Laplace approximation is exact for a gaussian
# response, so agreement is to optimizer noise, and every tolerance is a
# ratio to a quantity the run measures. Names no extension uses (xcs,
# toyfld), and local_frmtmb_registries() puts every registry back.

# An exchangeable block with a FIXED correlation of 0.5 and one
# estimated standard deviation: xcs(0 + t | g).
xcs_C <- function(d) matrix(0.5, d, d) + diag(0.5, d)
register_xcs <- function() {
  frmtmb_register_covstruct(
    "xcs",
    lpdf = function(b, theta, blk) {
      d <- blk[["dim"]]
      Sigma <- exp(2 * theta[1]) * xcs_C(d)
      dim(b) <- c(d, length(b) %/% d)
      sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(2 * theta[1]) * xcs_C(blk[["dim"]])
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = 1L, sd_idx = 1L)
}

# A field over the levels of g with a FIXED correlation matrix A over
# the levels and one standard deviation: toyfld(A, gr = g), which is
# core's (1 | gr(g, cov = A)) written through the field seam.
register_toyfld <- function() {
  frmtmb_register_covstruct(
    "toyfld",
    lpdf = function(b, theta, blk) {
      K <- exp(2 * theta[1]) * blk[["aux"]][["A"]]
      sum(RTMB::dmvnorm(b, 0, K, log = TRUE))
    },
    vcov = function(theta, blk) {
      matrix(exp(2 * theta[1]), 1, 1,
             dimnames = list("sd(toyfld)", "sd(toyfld)"))
    },
    npar = 1L, sd_idx = 1L,
    natural = function(theta, blk) {
      list(term = "sd(toyfld)", type = "sd", est = theta[1],
           grad = matrix(1, 1, 1))
    },
    draw = function(theta, blk) {
      K <- exp(2 * theta[1]) * blk[["aux"]][["A"]]
      drop(crossprod(chol(K), stats::rnorm(nrow(K))))
    },
    levels_coupled = TRUE,
    term = list(
      args = c("A", "gr"), usage = "toyfld(A, gr = g)",
      parse = function(args, env) {
        if (is.null(args[["A"]])) stop("toyfld() needs A")
        list(A = args[["A"]])
      },
      build = function(spec, gr, lookup) {
        A <- lookup(spec[["A"]], "toyfld()")
        lv <- levels(droplevels(as.factor(gr)))
        list(levels = lv, index = match(as.character(gr), lv),
             aux = list(A = A[lv, lv]))
      }))
}

# -log of the gaussian marginal y ~ N(X beta, Z S Z' + sigma^2 I), in
# closed form: the independent reference
dense_ll <- function(y, X, beta, Z, S, sigma) {
  V <- Z %*% S %*% t(Z) + sigma^2 * diag(length(y))
  R <- chol(V)
  r <- backsolve(R, y - X %*% beta, transpose = TRUE)
  -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(r^2))
}

xcs_data <- function(seed = 11L, G = 30L, per = 3L) {
  set.seed(seed)
  d <- expand.grid(t = factor(c("a", "b", "c")), rep = seq_len(per),
                   g = factor(sprintf("g%02d", seq_len(G))))
  u <- t(crossprod(chol(0.8^2 * xcs_C(3L)), matrix(rnorm(3 * G), 3)))
  d$x <- rnorm(nrow(d))
  d$y <- 0.4 + 0.7 * d$x + u[cbind(as.integer(d$g), as.integer(d$t))] +
    rnorm(nrow(d), 0, 0.6)
  d
}

test_that("a registered bar structure fits as its density says", {
  local_frmtmb_registries()
  register_xcs()
  d <- xcs_data()
  fit <- frm(bf(y ~ x + xcs(0 + t | g)), family = gaussian(), data = d)
  est <- fit$estimates
  beta <- est[["beta"]]
  s <- exp(est[["theta"]][1])
  sigma <- exp(est[["betad"]][["sigma_(Intercept)"]])
  X <- stats::model.matrix(~x, d)
  # the block, level-major: one level's three coefficients together
  Z <- matrix(0, nrow(d), 3L * nlevels(d$g))
  Z[cbind(seq_len(nrow(d)),
          (as.integer(d$g) - 1L) * 3L + as.integer(d$t))] <- 1
  S <- kronecker(diag(nlevels(d$g)), s^2 * xcs_C(3L))
  ref <- dense_ll(d$y, X, beta, Z, S, sigma)
  ll <- as.numeric(logLik(fit))
  # the same density at the same point: an identity, so rounding only
  expect_lt(abs(ll - ref), 1e-9 * abs(ref))
  # and the fit is at the maximum of the dense likelihood
  obj <- function(p) -dense_ll(d$y, X, p[1:2], Z,
                               kronecker(diag(nlevels(d$g)),
                                         exp(2 * p[3]) * xcs_C(3L)),
                               exp(p[4]))
  o <- stats::optim(c(beta, log(s), log(sigma)), obj, method = "BFGS",
                    control = list(reltol = 1e-14))
  expect_lt(abs(-o$value - ll), 1e-7 * abs(ll))
  expect_lte(-o$value - ll, 1e-7 * abs(ll))
  # VarCorr() reads the registered vcov at the estimate
  vc <- VarCorr(fit)
  sdg <- vc$g$sd[, "Estimate"]
  expect_lt(max(abs(sdg - s)), 1e-12 * s)
  # the SE checks and confint() read sd_idx: one log sd, and the fixed
  # correlation reported at its value with no width, as it is not
  # estimated
  expect_identical(unname(log_sd_theta_index(fit)), 1L)
  cv <- confint_varcorr(fit)
  expect_identical(cv$type, rep(c("sd", "cor"), each = 3L))
  expect_lt(max(cv$upr[4:6] - cv$lwr[4:6]), 1e-8 * cv$estimate[4])
  expect_equal(cv$estimate[4:6], rep(0.5, 3L), tolerance = 1e-12)
  # an unseen level draws from the registered vcov
  nd <- d[1:3, ]
  nd$g <- factor("new")
  p_new <- predict(fit, newdata = nd, allow_new_levels = TRUE,
                   summary = FALSE)
  expect_true(all(is.finite(p_new)))
  expect_identical(dim(simulate(fit, nsim = 2, seed = 3)), c(nrow(d), 2L))
  # and frm_compat() describes the structure it now accepts
  expect_true("xcs" %in% frm_compat_features()$name)
})

test_that("a registered field term fits as core's gr(cov = A) does", {
  local_frmtmb_registries()
  register_toyfld()
  set.seed(5)
  G <- 15L
  lv <- sprintf("s%02d", seq_len(G))
  A <- 0.6^abs(outer(seq_len(G), seq_len(G), "-"))
  dimnames(A) <- list(lv, lv)
  d <- data.frame(g = factor(rep(lv, each = 5L), levels = lv),
                  x = rnorm(5L * G))
  phi <- drop(crossprod(chol(0.9^2 * A), rnorm(G)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$g)] + rnorm(nrow(d), 0, 0.5)
  f_toy <- frm(bf(y ~ x + toyfld(A, gr = g)), family = gaussian(),
               data = d, data2 = list(A = A))
  f_core <- frm(bf(y ~ x + (1 | gr(g, cov = A))), family = gaussian(),
                data = d, data2 = list(A = A))
  ll <- as.numeric(logLik(f_core))
  expect_lt(abs(as.numeric(logLik(f_toy)) - ll), 1e-7 * abs(ll))
  b_toy <- f_toy$estimates[["beta"]]
  se <- sqrt(diag(vcov(f_core)))
  expect_lt(max(abs(b_toy - f_core$estimates[["beta"]]) / se), 1e-3)
  # the field's own row, under the name its natural slot gives
  cv <- confint_varcorr(f_toy)
  expect_identical(cv$term, "sd(toyfld)")
  ref_sd <- exp(f_core$estimates[["theta"]][1])
  expect_lt(abs(cv$estimate - ref_sd), 1e-3 * ref_sd)
  # a field is not a brms grouping factor
  expect_null(ngrps(f_toy))
  # the block is one draw of the whole field, through its draw slot
  expect_identical(dim(simulate(f_toy, nsim = 2, seed = 3)),
                   c(nrow(d), 2L))
  # levels_coupled: an unseen location has no law of its own, so it
  # predicts at the population level, as gr(cov = A) does
  nd <- d[1:2, ]
  nd$g <- factor(c("s01", "elsewhere"))
  p <- fitted(f_toy, newdata = nd, allow_new_levels = TRUE)
  p0 <- fitted(f_toy, newdata = nd, re_formula = NA)
  expect_identical(p[2, "Estimate"], p0[2, "Estimate"])
  expect_false(identical(p[1, "Estimate"], p0[1, "Estimate"]))
  pc <- fitted(f_core, newdata = nd, allow_new_levels = TRUE)
  expect_lt(max(abs(p[, "Estimate"] - pc[, "Estimate"])),
            1e-3 * max(abs(pc[, "Estimate"])))
  # the grammar is the parser's: a bar spelling and an unknown argument
  # are refused by name
  expect_error(frm(bf(y ~ x + (1 | toyfld(A, gr = g))), family = gaussian(),
                   data = d, data2 = list(A = A)),
               "toyfld\\(\\) is not a bar term; write toyfld\\(A, gr = g\\)")
  expect_error(frm(bf(y ~ x + toyfld(A, gr = g, rho = 1)),
                   family = gaussian(), data = d, data2 = list(A = A)),
               "toyfld\\(\\): unknown argument\\(s\\) rho")
})

test_that("the registries are as the tests found them", {
  # each test above registered inside local_frmtmb_registries()
  expect_false(any(c("xcs", "toyfld") %in% cs_names()))
  expect_false(any(c("xcs", "toyfld") %in% frm_compat_features()$name))
  # an unregistered structure around a bar term is refused by name; it
  # used to fit the default us block (frmtmb 0.73.0)
  d <- xcs_data(G = 4L)
  expect_error(frm(bf(y ~ x + xcs(0 + t | g)), family = gaussian(), data = d),
               "xcs() is not a covariance structure", fixed = TRUE)
})

test_that("registration refuses what it could not honor, and commits nothing", {
  local_frmtmb_registries()
  before <- cs_names()
  feats <- frm_compat_features()$name
  ok_lpdf <- function(b, theta, blk) sum(RTMB::dnorm(b, 0, exp(theta), TRUE))
  ok_vcov <- function(theta, blk) diag(exp(2 * theta[1]), blk[["dim"]])
  # a core structure, and a name the formula grammar already holds
  expect_error(frmtmb_register_covstruct("ar1", ok_lpdf, ok_vcov, 1L),
               "one of frmtmb's own covariance structures")
  expect_error(frmtmb_register_covstruct("gr", ok_lpdf, ok_vcov, 1L),
               "already means something else")
  expect_error(frmtmb_register_covstruct("not ok", ok_lpdf, ok_vcov, 1L),
               "syntactic name")
  # a name the compatibility vocabulary holds under another kind
  expect_error(frmtmb_register_covstruct("rescor", ok_lpdf, ok_vcov, 1L),
               "One display name carries one kind")
  expect_error(frmtmb_register_covstruct("tcs", ok_lpdf, ok_vcov, -1L),
               "whole numbers")
  expect_error(frmtmb_register_covstruct("tcs", ok_lpdf, function(x) x, 1L),
               "vcov")
  expect_error(frmtmb_register_covstruct(
    "tcs", ok_lpdf, ok_vcov, 1L,
    term = list(args = "A", parse = function(args, env) list(),
                build = function(spec, gr, lookup) list())),
    "\"gr\"")
  expect_identical(cs_names(), before)
  expect_identical(frm_compat_features()$name, feats)
  # the case the refusals must not reach: a correct registration takes,
  # and a second one from the same caller replaces it
  frmtmb_register_covstruct("tcs", ok_lpdf, ok_vcov, 1L)
  frmtmb_register_covstruct("tcs", ok_lpdf, ok_vcov, 1L, sd_idx = 1L)
  expect_identical(cs_names(), c(before, "tcs"))
  expect_identical(cs_entry("tcs")$sd_idx(list(dim = 2L)), 1L)
})

test_that("every built-in structure has the entry shape the seam builds", {
  # the seam is the path core uses: each consumer reads the same slots
  # off a built-in entry as off a registered one
  bad <- unlist(lapply(names(covstruct_core), function(nm) {
    cs_check_entry(nm, covstruct_core[[nm]])
  }))
  expect_null(bad)
  # and the checker is not vacuous
  expect_match(cs_check_entry("x", list(npar = 1L)), "missing")
})

test_that("a wrapper no structure answers to is refused, not fitted as us", {
  # reformulas finds a bar inside any call, so frmtmb 0.73.0 fitted
  # foo(1 | g) and us_t(1 | g) as (1 | g) with no word
  # (dev/spatial-wrapper-probe.R)
  set.seed(1)
  d <- data.frame(g = gl(10, 5), x = rnorm(50))
  d$y <- rnorm(50) + rnorm(10)[d$g]
  expect_error(frm(bf(y ~ x + foo(1 | g)), family = gaussian(), data = d),
               "foo() is not a covariance structure", fixed = TRUE)
  expect_error(frm(bf(y ~ x + us_t(1 | g)), family = gaussian(), data = d),
               "gr(g, dist = \"student\")", fixed = TRUE)
  # the absent case: a known wrapper, a plain bar, and a call that holds
  # no bar are untouched
  f1 <- frm(bf(y ~ x + diag(1 | g)), family = gaussian(), data = d)
  f2 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  expect_identical(f1$frame$re_blocks[[1]]$covstruct, "diag")
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-8)
  expect_no_error(frm(bf(y ~ log(exp(x)) + (1 | g)), family = gaussian(),
                      data = d, dry_run = "spec"))
})

# A scalar iid block, the plainest bar structure: myiid(1 | g) is
# (1 | g). The slips below are each one slot of it gone wrong.
iid_slots <- function() {
  list(lpdf = function(b, theta, blk) {
         sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
       },
       vcov = function(theta, blk) {
         matrix(exp(2 * theta[1]), 1, 1,
                dimnames = list(blk[["cnms"]], blk[["cnms"]]))
       },
       npar = 1L, sd_idx = 1L)
}
iid_data <- function() {
  set.seed(1)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + d$x + rnorm(20)[d$g] + rnorm(100, 0, 0.5)
  d
}

test_that("a block is checked against its structure when the frame is built", {
  local_frmtmb_registries()
  d <- iid_data()
  reg <- function(nm, ...) {
    s <- iid_slots()
    over <- list(...)
    s[names(over)] <- over
    do.call(frmtmb_register_covstruct, c(list(name = nm), s))
  }
  fit_of <- function(nm) {
    frm(bf(stats::as.formula(paste0("y ~ x + ", nm, "(1 | g)"))),
        family = gaussian(), data = d)
  }
  # the absent case: the correct structure fits as (1 | g)
  reg("myiid")
  ll0 <- as.numeric(logLik(frm(bf(y ~ x + (1 | g)), family = gaussian(),
                              data = d)))
  expect_equal(as.numeric(logLik(fit_of("myiid"))), ll0, tolerance = 1e-8)
  # each slip the review constructed fitted silently before
  reg("slipstart", start = c(0, 0, 0))
  expect_error(fit_of("slipstart"), "'slipstart' structure's start slot")
  reg("slipsd", sd_idx = 4L)
  expect_error(fit_of("slipsd"), "'slipsd' structure's sd_idx slot")
  reg("slipvc", vcov = function(theta, blk) diag(3))
  expect_error(fit_of("slipvc"), "'slipvc' structure's vcov slot")
  reg("slipnames", vcov = function(theta, blk) matrix(1, 1, 1))
  expect_error(fit_of("slipnames"), "gave dimnames")
  reg("slipinf", lpdf = function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE)) + log(0)
  })
  expect_error(fit_of("slipinf"), "did not give one finite number at b = 0")
  reg("slipnp", npar = function(blk) -1)
  expect_error(fit_of("slipnp"), "'slipnp' structure's npar slot")
  # and the registration itself refuses an entry it cannot complete
  expect_error(frmtmb_register_covstruct("slipnolpdf", NULL,
                                         iid_slots()$vcov, 1L), "lpdf")
})

test_that("a structure cannot take a name the rest of frmtmb owns", {
  local_frmtmb_registries()
  s <- iid_slots()
  try_reg <- function(nm, ...) {
    do.call(frmtmb_register_covstruct, c(list(name = nm), s, list(...)))
  }
  before <- cs_names()
  # a function formulas call as a term: poly(x, 2) would stop working
  expect_error(try_reg("poly"), "formulas call poly[(][)] as a term")
  expect_error(try_reg("ns"), "formulas call ns[(][)] as a term")
  expect_error(try_reg("log"), "formulas call log[(][)] as a term")
  expect_error(try_reg("s"), "cannot register 's'")
  # the brms vocabulary, core's families and addition terms
  for (nm in c("gaussian", "cox", "trials", "bhaz", "sar", "gr", "mm")) {
    expect_error(try_reg(nm), "cannot register", info = nm)
  }
  # an extension's names, which only the home package registers
  for (nm in c("dec", "wiener")) {
    expect_error(try_reg(nm), "and only frmtmb.eam registers it", info = nm)
  }
  expect_identical(cs_names(), before)
  # the absent case: poly() in a formula is unharmed by the refusal
  d <- iid_data()
  expect_no_error(frm(bf(y ~ poly(x, 2)), family = gaussian(), data = d,
                      dry_run = "spec"))
  # functions formulas do not call are free names for a structure, and
  # deciding loads no package
  for (nm in c("toeplitz", "kernel", "matern", "rw1")) {
    expect_no_error(try_reg(nm))
    expect_true(nm %in% cs_names(), info = nm)
  }
})

test_that("a correct density of any shape passes the block check", {
  local_frmtmb_registries()
  # An intrinsic RW1 over 30 times, the local-level model with no
  # intercept: its density is flat along the constant vector, which the
  # data identify through the field.
  n <- 30L
  Q1 <- crossprod(diff(diag(n)))
  dimnames(Q1) <- list(as.character(seq_len(n)), as.character(seq_len(n)))
  frmtmb_register_covstruct(
    "rwone",
    lpdf = function(b, theta, blk) {
      Q <- blk[["aux"]][["Q"]]
      k <- nrow(Q) - 1
      -k * theta[1] - 0.5 * exp(-2 * theta[1]) * sum(b * (Q %*% b)) -
        0.5 * k * log(2 * pi)
    },
    vcov = function(theta, blk) {
      matrix(exp(2 * theta[1]), 1, 1, dimnames = list("sd", "sd"))
    },
    npar = 1L, sd_idx = 1L,
    natural = function(theta, blk) {
      list(term = "sd(rwone)", type = "sd", est = theta[1],
           grad = matrix(1, 1, 1))
    },
    levels_coupled = TRUE,
    term = list(args = c("Q", "gr"), usage = "rwone(Q, gr = t)",
                parse = function(args, env) list(Q = args[["Q"]]),
                build = function(spec, gr, lookup) {
                  Q <- lookup(spec[["Q"]], "rwone()")
                  lv <- rownames(Q)
                  list(levels = lv, index = match(as.character(gr), lv),
                       aux = list(Q = Q))
                }))
  set.seed(7)
  d <- data.frame(t = seq_len(n), x = rnorm(n))
  d$y <- 0.5 * d$x + cumsum(rnorm(n, 0, 0.7)) + rnorm(n, 0, 0.5)
  f <- frm(bf(y ~ 0 + x + rwone(Q1, gr = t)), family = gaussian(),
           data = d, data2 = list(Q1 = Q1))
  expect_true(is.finite(as.numeric(logLik(f))))
  # A mean far from zero: the shift goes into the intercept, so the fit
  # is (1 | g)'s.
  dg <- iid_data()
  s <- iid_slots()
  s$lpdf <- function(b, theta, blk) {
    sum(RTMB::dnorm(b, 600, exp(theta[1]), log = TRUE))
  }
  do.call(frmtmb_register_covstruct, c(list(name = "shiftmean"), s))
  ll0 <- as.numeric(logLik(frm(bf(y ~ x + (1 | g)), family = gaussian(),
                               data = dg)))
  f <- frm(bf(y ~ x + shiftmean(1 | g)), family = gaussian(), data = dg)
  expect_equal(as.numeric(logLik(f)), ll0, tolerance = 1e-6)
  # Bimodal, modes at -1 and 1 and a local minimum at 0: no test of the
  # level or the curvature at zero tells it from a negated density.
  s$lpdf <- function(b, theta, blk) {
    sd <- exp(theta[1])
    sum(log(0.5 * exp(RTMB::dnorm(b, -1, sd, log = TRUE)) +
              0.5 * exp(RTMB::dnorm(b, 1, sd, log = TRUE))))
  }
  s$vcov <- function(theta, blk) {
    matrix(1 + exp(2 * theta[1]), 1, 1,
           dimnames = list(blk[["cnms"]], blk[["cnms"]]))
  }
  do.call(frmtmb_register_covstruct, c(list(name = "twomodes"), s))
  th <- log(0.3)
  expect_lt(s$lpdf(numeric(3), th, NULL), s$lpdf(rep(1, 3), th, NULL))
  # fitted to data drawn from it, it recovers its sd (0.3)
  set.seed(3)
  u <- sample(c(-1, 1), 40, TRUE) + rnorm(40, 0, 0.3)
  d2 <- data.frame(g = factor(rep(1:40, each = 8)), x = rnorm(320))
  d2$y <- 1 + d2$x + u[d2$g] + rnorm(320, 0, 0.3)
  f <- frm(bf(y ~ x + twomodes(1 | g)), family = gaussian(), data = d2)
  expect_true(is.finite(as.numeric(logLik(f))))
  expect_lt(abs(log(exp(f$estimates[["theta"]][1]) / 0.3)), log(2))
})

test_that("another package's structure is replaced only on request", {
  local_frmtmb_registries()
  s <- iid_slots()
  do.call(frmtmb_register_covstruct,
          c(list(name = "ownedby"), s, list(package = "frmtmbtoyowner")))
  # the session is not the owner
  expect_error(do.call(frmtmb_register_covstruct, c(list(name = "ownedby"), s)),
               "already registered by the package frmtmbtoyowner")
  # the owner replaces it, and so does an explicit replace = TRUE
  expect_no_error(do.call(frmtmb_register_covstruct,
                          c(list(name = "ownedby"), s,
                            list(package = "frmtmbtoyowner"))))
  expect_no_error(do.call(frmtmb_register_covstruct,
                          c(list(name = "ownedby"), s,
                            list(replace = TRUE))))
  expect_null(cs_entry("ownedby")[["package"]])
})

test_that("frm_compat() refuses quadrature and importance as frm() does", {
  local_frmtmb_registries()
  do.call(frmtmb_register_covstruct, c(list(name = "myiid"), iid_slots()))
  for (m in c("quadrature", "importance")) {
    expect_identical(frm_compat("myiid", m)$status, "refused", info = m)
  }
  d <- iid_data()
  d$yb <- rbinom(100, 1, 0.5)
  expect_error(frm(bf(yb ~ x + myiid(1 | g)), family = bernoulli(),
                   data = d, quadrature = TRUE), "quadrature")
  expect_error(frm(bf(yb ~ x + myiid(1 | g)), family = bernoulli(),
                   data = d, importance = 50L), "importance")
  # the absent case: frmtmb's own structures keep their rows
  expect_identical(frm_compat("us", "quadrature")$status, "conditional")
  expect_identical(frm_compat("ar1", "importance")$status, "conditional")
})

test_that("a nested wrapper and a re_formula wrapper are refused by name", {
  set.seed(1)
  d <- data.frame(g = gl(10, 6), x = rnorm(60))
  d$y <- rnorm(60) + rnorm(10)[d$g]
  # frmtmb 0.73.0 and the first fix fitted this as (1 | g)
  expect_error(frm(bf(y ~ x + foo(bar(1 | g))), family = gaussian(),
                   data = d), "foo() is not a covariance structure",
               fixed = TRUE)
  expect_error(frm(bf(y ~ x + diag(foo(1 | g))), family = gaussian(),
                   data = d), "foo() is not a covariance structure",
               fixed = TRUE)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  expect_error(frm_linpred(fit, re_formula = ~ foo(1 | g)),
               "re_formula: foo() is not a covariance structure",
               fixed = TRUE)
  # the absent case: the spellings a model writes are untouched
  expect_identical(frm_linpred(fit, re_formula = ~ (1 | g)),
                   frm_linpred(fit))
  expect_identical(frm_linpred(fit, re_formula = ~ us(1 | g)),
                   frm_linpred(fit))
  f2 <- frm(bf(y ~ x + diag((1 | g))), family = gaussian(), data = d)
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(fit)),
               tolerance = 1e-8)
})

test_that("a fit whose structure's package is absent is refused by name", {
  d <- iid_data()
  make_fits <- function() {
    local_frmtmb_registries()
    do.call(frmtmb_register_covstruct,
            c(list(name = "awayiid"), iid_slots(),
              list(package = "frmtmbnotinstalled")))
    do.call(frmtmb_register_covstruct,
            c(list(name = "sessiid"), iid_slots()))
    list(away = frm(bf(y ~ x + awayiid(1 | g)), family = gaussian(),
                    data = d),
         sess = frm(bf(y ~ x + sessiid(1 | g)), family = gaussian(),
                    data = d))
  }
  fits <- make_fits()
  # the registrations are gone, as in a session that reads a saved fit
  expect_false(any(c("awayiid", "sessiid") %in% cs_names()))
  msg <- paste0("The awayiid() structure of this model is from the ",
                "extension package frmtmbnotinstalled, which is not ",
                "installed")
  calls <- list(
    print = function(f) utils::capture.output(print(f)),
    summary = function(f) summary(f), vcov = function(f) vcov(f),
    predict = function(f) predict(f), fitted = function(f) fitted(f),
    confint = function(f) confint(f), VarCorr = function(f) VarCorr(f),
    ranef = function(f) ranef(f), update = function(f) update(f),
    simulate = function(f) simulate(f, nsim = 1),
    confint_varcorr = function(f) confint_varcorr(f))
  for (nm in names(calls)) {
    e <- tryCatch(calls[[nm]](fits$away), error = function(e) e)
    expect_s3_class(e, "frmtmb_extension_missing")
    expect_match(conditionMessage(e), msg, fixed = TRUE, info = nm)
    e2 <- tryCatch(calls[[nm]](fits$sess), error = function(e) e)
    expect_match(conditionMessage(e2), "registered in an R session",
                 fixed = TRUE, info = nm)
  }
})

test_that("a field with one row per level is not an observation-level effect", {
  # levels_coupled is the property the OLRE check reads: a fixed matrix
  # correlates the levels, so the field's variance and the residual sd
  # are separately identified, as for gr(cov = A). 0.73.0 named three
  # structures there and warned on every car() fit with one row per area
  local_frmtmb_registries()
  register_toyfld()
  set.seed(8)
  G <- 30L
  lv <- sprintf("s%02d", seq_len(G))
  A <- 0.7^abs(outer(seq_len(G), seq_len(G), "-"))
  dimnames(A) <- list(lv, lv)
  d <- data.frame(g = factor(lv, levels = lv), x = rnorm(G))
  d$y <- 1 + d$x + drop(crossprod(chol(A), rnorm(G))) + rnorm(G, 0, 0.3)
  olre <- "gives every observation its own random effect"
  seen <- character()
  withCallingHandlers(
    frm(bf(y ~ x + toyfld(A, gr = g)), family = gaussian(), data = d,
        data2 = list(A = A)),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl(olre, seen, fixed = TRUE)))
  # the case the check exists for still warns: an iid effect per row
  allow_warnings(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d),
                 olre, require = olre)
})
