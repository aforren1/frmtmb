# A covariate far from zero against its spread: `x + m` is the model on
# `x` reparametrized exactly, so the slope, its standard error and the
# log-likelihood must not move. Up to 0.79.1 the optimizer stopped at its
# start with code 0 at m = 1e6 (poisson with (1 | g): slope 6 to 9
# standard errors off, 20 to 47 log-likelihood units short, mostly
# silent), and the frame dropped the slope as aliased from m = 1e7
# (dev/bigx-findings.md).

offset_data <- function(family, m, seed, re = TRUE, n = 150L, G = 15L) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), length.out = n))
  z <- stats::rnorm(n)
  u <- if (re) stats::rnorm(G, 0, 0.7)[g] else 0
  eta <- 0.5 + 0.4 * z + u
  y <- switch(family,
    poisson = stats::rpois(n, exp(eta)),
    bernoulli = stats::rbinom(n, 1, stats::plogis(eta)),
    gaussian = eta + stats::rnorm(n, 0, 0.5))
  data.frame(y = y, w = z + m, z = z, g = g)
}

# the slope's shift in standard errors of the fit on z, the ratio of the
# two standard errors, and the log-likelihoods' difference relative to
# the one on z
offset_compare <- function(fw, fz) {
  # from about 1e7 the whole matrix says that a product formed from it
  # loses digits (warn_cov_products()); its diagonal is accurate
  sw <- sqrt(diag(allow_warnings(vcov(fw), "is ill-conditioned")))
  sz <- sqrt(diag(vcov(fz)))
  list(shift = (fixef_vec(fw)[[2]] - fixef_vec(fz)[[2]]) / sz[[2]],
       se_ratio = sw[[2]] / sz[[2]],
       ll_rel = (as.numeric(logLik(fw)) - as.numeric(logLik(fz))) /
         abs(as.numeric(logLik(fz))))
}

fixef_vec <- function(f) unname(f$estimates[["beta"]])

expect_same_model <- function(cmp) {
  expect_lt(abs(cmp$shift), 1e-2)
  expect_lt(abs(cmp$se_ratio - 1), 1e-3)
  expect_lt(abs(cmp$ll_rel), 1e-8)
}

test_that("a poisson GLMM on x + 1e6 is the GLMM on x", {
  d <- offset_data("poisson", 1e6, 2)
  fz <- frm(bf(y ~ z + (1 | g)), family = poisson(), data = d)
  expect_no_warning(
    fw <- frm(bf(y ~ w + (1 | g)), family = poisson(), data = d))
  expect_equal(fw$opt$convergence, 0L)
  expect_same_model(offset_compare(fw, fz))
})

test_that("a bernoulli GLMM on x + 1e5 is the GLMM on x", {
  d <- offset_data("bernoulli", 1e5, 2)
  fz <- frm(bf(y ~ z + (1 | g)), family = bernoulli(), data = d)
  expect_no_warning(
    fw <- frm(bf(y ~ w + (1 | g)), family = bernoulli(), data = d))
  expect_same_model(offset_compare(fw, fz))
})

test_that("a column 1e8 from zero against sd 1 is kept, not dropped", {
  d <- offset_data("gaussian", 1e8, 3, re = FALSE)
  fz <- frm(bf(y ~ z), family = gaussian(), data = d)
  fw <- NULL
  # a stop at "false convergence (8)" at the maximum is a message of its
  # own class, and only that message is allowed
  expect_no_message(
    fw <- allow_at_max(allow_warnings(frm(bf(y ~ w), family = gaussian(),
                                          data = d),
                                      "false convergence (8)")))
  expect_length(fixef_vec(fw), 2L)
  cmp <- offset_compare(fw, fz)
  expect_lt(abs(cmp$shift), 1e-2)
  expect_lt(abs(cmp$se_ratio - 1), 1e-3)
  # an exact alias is still dropped: a constant at the same distance
  d$k <- 1e8
  expect_message(frm(bf(y ~ w + k), family = gaussian(), data = d,
                     dry_run = "frame"),
                 "dropping column(s): k", fixed = TRUE)
})

test_that("only an ill-conditioned design gets the design's coordinates", {
  d <- offset_data("poisson", 1e6, 2)
  fz <- frm(bf(y ~ z + (1 | g)), family = poisson(), data = d)
  fw <- frm(bf(y ~ w + (1 | g)), family = poisson(), data = d)
  expect_null(fz$obj$env$.frmtmb_opt_basis)
  B <- fw$obj$env$.frmtmb_opt_basis
  expect_false(is.null(B))
  expect_equal(B$pos, 1:2)
  # the inverse is the inverse, to the rounding of the product
  rnd <- 64 * .Machine$double.eps * max(abs(B$M)) * max(abs(B$Mi))
  expect_lt(max(abs(B$M %*% B$Mi - diag(3))), rnd)
  opt_coords <- getFromNamespace("opt_coords", "frmtmb")
  open <- list(lower = -Inf, upper = Inf)
  expect_false(is.null(opt_coords(fw$obj, "nlminb", open, NULL)$basis))
  # a caller's optimizer keeps the coordinates it was handed; frmtmb's
  # own wrappers do not
  expect_null(opt_coords(fw$obj, function(...) NULL, open, NULL)$basis)
  own <- structure(function(...) NULL, frmtmb_builtin = TRUE)
  expect_false(is.null(opt_coords(fw$obj, own, open, NULL)$basis))
  # a bound on the slope: the centering map, in which the slope is its
  # own coordinate (so the box stays a box) and the intercept centered
  held <- list(lower = c(-Inf, 0, -Inf), upper = Inf)
  cs <- opt_coords(fw$obj, "nlminb", held, NULL)
  expect_false(is.null(cs$basis))
  p <- fw$opt$par
  q <- cs$to_q(p)
  expect_identical(q[2:3], unname(p[2:3]))
  expect_equal(q[1], unname(p[1] + mean(d$w) * p[2]))
  expect_identical(cs$lower[2], 0)
  expect_null(fw$obj$env$.frmtmb_basis_note)
  # a bound on a parameter outside the block leaves the block's
  # coordinates in place
  th <- list(lower = c(-Inf, -Inf, -3), upper = Inf)
  expect_identical(opt_coords(fw$obj, "nlminb", th, NULL)$basis$M,
                   B$M)
  # a bound on the intercept leaves nothing to center on: natural
  # coordinates, and the fit says so by class
  ic <- list(lower = c(-1e9, -Inf, -Inf), upper = Inf)
  expect_null(opt_coords(fw$obj, "nlminb", ic, NULL)$basis)
  expect_warning(getFromNamespace("opt_coords_warn", "frmtmb")(fw$obj),
                 class = "frmtmb_ill_design")
  fw$obj$env$.frmtmb_basis_note <- NULL
})

test_that("a bound on a shifted slope fits the centered model", {
  # any lb or ub on class b used to send the fit back to natural
  # coordinates, where it stopped 6 to 9 standard errors off in silence
  d <- offset_data("poisson", 1e6, 1)
  se_ref <- unname(sqrt(diag(vcov(
    frm(bf(y ~ z + (1 | g)), family = poisson(), data = d))))[2])
  for (lb in c(-5, 0.45)) {
    pr <- set_prior("normal(0, 10)", class = "b", lb = lb)
    fz <- allow_warnings(frm(bf(y ~ z + (1 | g)), family = poisson(),
                             data = d, prior = pr),
                         "Standard errors are not available")
    w <- character(0)
    fw <- withCallingHandlers(
      frm(bf(y ~ w + (1 | g)), family = poisson(), data = d, prior = pr),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    # lb = 0.45 binds near the estimate: the slope sits on its bound, and
    # the only warning allowed is that its standard error is not there
    expect_true(all(grepl("Standard errors are not available", w,
                          fixed = TRUE)))
    if (lb < 0) expect_length(w, 0L)
    expect_equal(fw$opt$convergence, 0L)
    expect_lt(abs(fixef_vec(fw)[[2]] - fixef_vec(fz)[[2]]) / se_ref, 1e-3)
    expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
              1e-8)
  }
})

test_that("a nonlinear body linear in a shifted column fits the centered one", {
  d <- offset_data("gaussian", 1e6, 1)
  fz <- frm(bf(y ~ b0 + b1 * z, b0 ~ 1 + (1 | g), b1 ~ 1, nl = TRUE),
            family = gaussian(), data = d)
  expect_no_warning(
    fw <- frm(bf(y ~ b0 + b1 * w, b0 ~ 1 + (1 | g), b1 ~ 1, nl = TRUE),
              family = gaussian(), data = d))
  b1 <- function(f) unname(f$estimates[["beta"]][2])
  se1 <- function(f) unname(sqrt(diag(vcov(f)))[2])
  expect_lt(abs(b1(fw) - b1(fz)) / se1(fz), 1e-2)
  expect_lt(abs(se1(fw) / se1(fz) - 1), 1e-3)
  expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1), 1e-8)
})

test_that("a raw polynomial far from zero keeps only columns it can estimate", {
  # x on [1000, 1010], degree 8: the columns aliased at the rank
  # tolerance are dropped with lm()'s message, those that take the
  # centered design past the condition number the optimizer and the
  # standard errors can use with a warning of their own (the model is
  # not the one written), and the rest have finite standard errors
  set.seed(501)
  n <- 200
  d <- data.frame(xx = stats::runif(n, 1000, 1010))
  d$yy <- sin(d$xx / 2) + stats::rnorm(n, 0, 0.3)
  fit <- NULL
  expect_message(
    fit <- allow_warnings(frm(bf(yy ~ poly(xx, 8, raw = TRUE)),
                              family = gaussian(), data = d),
                          "but they are dropped",
                          require = "The model fitted is not the one written"),
    "rank deficient; dropping column(s)", fixed = TRUE)
  X <- fit$frame$linpreds[[1]]$X
  Xc <- X
  Xc[, -1] <- X[, -1] - rep(colMeans(X[, -1]), each = n)
  kap <- getFromNamespace("lp_design_kappa", "frmtmb")(list(X = Xc))
  expect_lte(kap, getFromNamespace("design_kappa_max", "frmtmb"))
  se <- sqrt(diag(allow_warnings(vcov(fit), "is ill-conditioned")))
  expect_true(all(is.finite(se)))
  g <- stats::glm.fit(X, d$yy)
  llg <- -n / 2 * (log(2 * pi * sum(g$residuals^2) / n) + 1)
  expect_lt(abs(as.numeric(logLik(fit)) / llg - 1), 1e-6)
})

test_that("products and powers of shifted columns are kept and fitted", {
  # at an offset of 1e5 the raw design's condition number is above 1e10
  # and the centered one's is not: the frame keeps w1:w2 and I(w1^2),
  # which the base build dropped, 40 and 5.7 log-likelihood units short.
  # A stop at the maximum may say "false convergence (8)", as a warning or
  # as the message of class frmtmb_converged_at_max
  set.seed(1)
  n <- 200
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  d <- data.frame(z1, z2, w1 = z1 + 1e5, w2 = z2 + 1e5)
  d$y <- 1 + 0.5 * z1 - 0.3 * z2 + 0.4 * z1 * z2 + stats::rnorm(n, 0, 0.5)
  d$yp <- stats::rpois(n, exp(0.3 + 0.4 * z1 - 0.2 * z1^2))
  cases <- list(
    list(bf(y ~ w1 * w2), bf(y ~ z1 * z2), gaussian(), 4L),
    list(bf(yp ~ w1 + I(w1^2)), bf(yp ~ z1 + I(z1^2)), poisson(), 3L))
  for (cc in cases) {
    fz <- frm(cc[[2]], family = cc[[3]], data = d)
    fw <- NULL
    expect_no_message(
      fw <- allow_at_max(allow_warnings(frm(cc[[1]], family = cc[[3]],
                                            data = d),
                                        "false convergence (8)")))
    expect_length(fw$estimates[["beta"]], cc[[4]])
    expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
              1e-6)
  }
  # at 1e6 the data still identify both columns: the frame keeps them
  # and fits them, or drops them with the warning that says the model
  # is not the one written, never with "rank deficient"
  d$w1 <- d$z1 + 1e6
  d$w2 <- d$z2 + 1e6
  for (cc in cases) {
    fz <- frm(cc[[2]], family = cc[[3]], data = d)
    cls <- character(0)
    fw <- NULL
    expect_no_message(fw <- allow_at_max(withCallingHandlers(
      frm(cc[[1]], family = cc[[3]], data = d),
      warning = function(x) {
        cls <<- c(cls, if (inherits(x, "frmtmb_design_capped")) "capped"
                  else conditionMessage(x))
        invokeRestart("muffleWarning")
      })))
    if (length(fw$estimates[["beta"]]) == cc[[4]]) {
      expect_true(all(grepl("false convergence (8)", cls, fixed = TRUE)))
      expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
                1e-6)
    } else {
      expect_true("capped" %in% cls)
    }
  }
})

test_that("the profile of a shifted slope is the centered one", {
  d <- offset_data("poisson", 1e6, 2)
  fz <- frm(bf(y ~ z + (1 | g)), family = poisson(), data = d)
  fw <- frm(bf(y ~ w + (1 | g)), family = poisson(), data = d)
  expect_no_warning(cw <- confint(fw, parm = "w", method = "profile"))
  cz <- confint(fz, parm = "z", method = "profile")
  se <- unname(sqrt(diag(vcov(fz)))[2])
  expect_true(all(is.finite(cw[1, 1:2])))
  expect_lt(max(abs(cw[1, 1:2] - cz[1, 1:2])) / se, 1e-3)
  rw <- confint(fw, parm = "w", method = "uniroot")
  rz <- confint(fz, parm = "z", method = "uniroot")
  expect_lt(max(abs(rw[1, 1:2] - rz[1, 1:2])) / se, 1e-3)
})

test_that("an iteration limit keeps its code and its warning", {
  d <- offset_data("poisson", 1e6, 2)
  w <- character(0)
  f <- withCallingHandlers(
    frm(bf(y ~ w + (1 | g)), family = poisson(), data = d,
        control = frmtmb_control(optCtrl = list(iter.max = 3,
                                                eval.max = 400))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_false(identical(f$opt$convergence, 0L))
  expect_null(f$opt$natural_finish)
  expect_true(any(grepl("iteration limit", w, fixed = TRUE)))
})

test_that("natural_finish() keeps a stop the natural run cannot improve", {
  natural_finish <- getFromNamespace("natural_finish", "frmtmb")
  obj <- RTMB::MakeADFun(function(p) (p$x - 3)^2 + 100 * (p$y - p$x^2)^2,
                         list(x = -1, y = 1), silent = TRUE)
  opt <- list(par = c(x = -1, y = 1), objective = obj$fn(c(-1, 1)),
              convergence = 1L, message = "false convergence (8)")
  ctl <- list(optCtrl = list(iter.max = 1, eval.max = 2))
  out <- natural_finish(obj, opt, "nlminb", obj$fn, obj$gr,
                        list(lower = -Inf, upper = Inf), ctl, NULL, NULL)
  expect_identical(out, opt)
  # from far enough the natural run converges, and is kept
  ctl$optCtrl <- list(iter.max = 1000, eval.max = 1000)
  out2 <- natural_finish(obj, opt, "nlminb", obj$fn, obj$gr,
                         list(lower = -Inf, upper = Inf), ctl, NULL, NULL)
  expect_equal(out2$convergence, 0L)
  expect_lt(out2$objective, opt$objective)
})

test_that("autoscale = TRUE on x + 1e5 keeps its standard errors", {
  d <- offset_data("poisson", 1e5, 1)
  fz <- frm(bf(y ~ z + (1 | g)), family = poisson(), data = d)
  expect_no_warning(
    fw <- frm(bf(y ~ w + (1 | g)), family = poisson(), data = d,
              control = frmtmb_control(autoscale = TRUE)))
  expect_same_model(offset_compare(fw, fz))
})

test_that("a restart in the design's coordinates is judged by its gain", {
  restart_result <- getFromNamespace("restart_result", "frmtmb")
  opt <- list(par = c(-4e5, 0.4), objective = 250, convergence = 0L,
              message = "relative convergence (4)", iterations = 8L,
              evaluations = c(12L, 9L))
  # moved within the objective's rounding, gained nothing, and stopped
  # at a code of its own
  opt2 <- list(par = c(-4e5 + 0.1, 0.4 + 1e-7), objective = 250 - 1e-12,
               convergence = 1L, message = "false convergence (8)",
               iterations = 4L, evaluations = c(6L, 4L))
  expect_equal(restart_result(opt, opt2)$convergence, 1L)
  kept <- restart_result(opt, opt2, any_move = TRUE)
  expect_equal(kept$convergence, 0L)
  expect_equal(kept$restart$convergence, 1L)
  # a restart that gained is its own stop either way
  opt3 <- opt2
  opt3$objective <- 249
  expect_equal(restart_result(opt, opt3, any_move = TRUE)$convergence, 1L)
})

test_that("a product of shifted columns beside `0 + f` is judged centered", {
  # `0 + f` spans the constant with no intercept column, and the rank
  # test read such a design uncentered: it dropped w1:w2 at 1e5 (capped)
  # and at 1e6 with the message "rank deficient", 84 to 101
  # log-likelihood units short (release review of 0.80.0)
  set.seed(7)
  n <- 300
  d <- data.frame(z1 = stats::rnorm(n), z2 = stats::rnorm(n),
                  f = factor(sample(letters[1:3], n, TRUE)))
  d$y <- c(a = 0, b = 1, c = -1)[as.character(d$f)] + 0.5 * d$z1 -
    0.4 * d$z2 + 0.6 * d$z1 * d$z2 + stats::rnorm(n, 0, 0.5)
  fz <- frm(y ~ 0 + f + z1 * z2, data = d)
  for (m in c(1e5, 1e6)) {
    d$w1 <- d$z1 + m
    d$w2 <- d$z2 + m
    cls <- character(0)
    fw <- NULL
    expect_no_message(fw <- allow_at_max(withCallingHandlers(
      frm(y ~ 0 + f + w1 * w2, data = d),
      warning = function(x) {
        cls <<- c(cls, if (inherits(x, "frmtmb_design_capped")) "capped"
                  else class(x)[1L])
        invokeRestart("muffleWarning")
      })))
    if (m == 1e5) {
      # the centered design is well conditioned here: kept and fitted
      expect_length(fw$estimates[["beta"]], 6L)
      expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
                1e-6)
    } else if (length(fw$estimates[["beta"]]) == 6L) {
      expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
                1e-6)
    } else {
      expect_true("capped" %in% cls)
    }
  }
})
