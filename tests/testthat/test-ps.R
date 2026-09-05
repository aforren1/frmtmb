# ps(): a penalized coefficient block whose VALUE a nonlinear body
# consumes. dev/spline-seam-proposal.md, Part 2.

test_that("the branch-free basis agrees with splines::splineDesign", {
  knots <- frmtmb:::ps_knots(0, 1.5, 12L, 4L)
  W <- frmtmb:::ps_weights(knots, 4L)
  x <- seq(0.01, 1.49, length.out = 37)
  B <- frmtmb:::ps_design(x, knots, 4L, W)
  R <- splines::splineDesign(knots, x, ord = 4L, outer.ok = TRUE)
  expect_equal(dim(B), dim(R))
  expect_lt(max(abs(B - R)), 1e-11)
  # B-splines are a partition of unity on the interior
  expect_lt(max(abs(rowSums(B) - 1)), 1e-11)

  # the same values with x as an AD input, and the taped derivative
  tp <- RTMB::MakeTape(function(p) as.vector(
    frmtmb:::ps_design(p, knots, 4L, W)), x)
  expect_lt(max(abs(matrix(tp(x), length(x), 12L) - R)), 1e-11)
  J <- tp$jacobian(x)
  Dref <- splines::splineDesign(knots, x, ord = 4L, derivs = 1,
                                outer.ok = TRUE)
  Dtape <- vapply(seq_len(12L), function(i)
    diag(J[(i - 1L) * length(x) + seq_along(x), ]), numeric(length(x)))
  expect_lt(max(abs(Dtape - Dref)), 1e-10)

  # exactly zero outside the knot span, which is why the fit reports
  # its coverage
  out <- frmtmb:::ps_design(c(-1, 3), knots, 4L, W)
  expect_lt(max(abs(out)), 1e-9)
})

test_that("ps() refuses to be called as a function and validates its args", {
  expect_error(ps(1), "is a term, not a function")
  d <- data.frame(y = rnorm(40), t = runif(40))
  expect_error(frm(bf(y ~ a * ps(t, k = 2), a ~ 1, nl = TRUE), d, gaussian()),
               "at least")
  expect_error(frm(bf(y ~ a * ps(t, k = 80), a ~ 1, nl = TRUE), d, gaussian()),
               "above 50")
  expect_error(frm(bf(y ~ a * ps(t, degree = 9), a ~ 1, nl = TRUE), d,
                   gaussian()), "between 1 and 5")
  # a literal, not an expression that happens to evaluate: `length(t)`
  # resolves `t` to base::t in any environment the parser could use, so
  # accepting it would accept a basis size of 1
  expect_error(frm(bf(y ~ a * ps(t, k = length(t)), a ~ 1, nl = TRUE), d,
                   gaussian()), "must be a literal constant")
  expect_error(frm(bf(y ~ a * ps(t, k = kk), a ~ 1, nl = TRUE), d,
                   gaussian()), "must be a literal constant")
  expect_error(frm(bf(y ~ a * ps(ps(t)), a ~ 1, nl = TRUE), d, gaussian()),
               "ps\\(\\) inside ps\\(\\)")
  # the three the proposal listed and this does not implement: refused
  # by name rather than through match.call()'s "unused argument"
  d$g <- factor(rep(1:4, 10))
  for (bad in c("ps(t, by = g)", "ps(t, id = 1)", "ps(t, penalty = 1)")) {
    f <- eval(parse(text = paste0("bf(y ~ a * ", bad,
                                  ", a ~ 1, nl = TRUE)")))
    expect_error(frm(f, d, gaussian()), "unknown argument", label = bad)
  }
})

test_that("the penalty eigensplit is Wood (2004)'s, and centred", {
  skip_on_cran()
  set.seed(7071)
  d <- data.frame(t = runif(150, 0, 1))
  d$y <- 2 + sin(2 * pi * d$t) + rnorm(150, 0, 0.15)
  fit <- frm(bf(y ~ lev + ps(t, k = 12), lev ~ 1, nl = TRUE), d, gaussian())
  pt <- fit$frame$linpreds[["y.mu"]]$ps_terms[[1]]

  m <- pt$k
  D2 <- matrix(0, m - 2L, m)
  for (i in seq_len(m - 2L)) D2[i, i + 0:2] <- c(1, -2, 1)
  S <- crossprod(D2)
  # the null space goes to the fixed coefficients, so S U0 is zero
  expect_lt(max(abs(S %*% pt$U0)), 1e-10)
  # the range space is rescaled so that the penalty IS the sum of
  # squares of the random-effect block: Us' S Us = I
  expect_lt(max(abs(crossprod(pt$Us, S %*% pt$Us) - diag(ncol(pt$Us)))), 1e-9)
  # the sum-to-zero constraint, which is what makes the curve's level
  # identifiable against the intercept in the body
  expect_lt(max(abs(c(colSums(pt$U0), colSums(pt$Us)))), 1e-12)
  expect_identical(pt$n_fixed + pt$n_pen, m - 1L)

  # the block is one variance in theta, reported under the term's label
  vc <- VarCorr(fit)
  expect_true(pt$label %in% names(vc))
  expect_identical(length(fit$frame$re_blocks[[pt$block_id]]$theta_idx), 1L)
})

test_that("a ps() curve recovers a known shape and predicts on newdata", {
  skip_on_cran()
  set.seed(7072)
  n_id <- 40
  d <- expand.grid(t = seq(0, 1, length.out = 10), id = factor(1:n_id))
  sh <- rnorm(n_id, 0, 0.06)[d$id]
  d$y <- 3 + sin(2 * pi * (d$t + sh)) + rnorm(nrow(d), 0, 0.15)
  fit <- suppressWarnings(
    frm(bf(y ~ lev + ps(t + shift, k = 10, pad = 0.4),
           lev ~ 1, shift ~ 0 + (1 | id), nl = TRUE), d, gaussian()))

  expect_equal(fixef(fit)$lev[[1]], 3, tolerance = 0.15)
  expect_equal(exp(fixef(fit)$sigma[[1]]), 0.15, tolerance = 0.05)
  expect_equal(sqrt(as.numeric(VarCorr(fit)[["shift: 1 | id"]])[1]), 0.06,
               tolerance = 0.04)

  grid <- seq(0.05, 0.95, length.out = 15)
  nd <- data.frame(t = grid, id = factor(1, levels = levels(d$id)))
  pop <- predict(fit, newdata = nd, re.form = NA)
  expect_lt(sqrt(mean((pop - (3 + sin(2 * pi * grid)))^2)), 0.15)

  # the frozen basis: prediction on new data uses the FITTED knots, so
  # a grid over a sub-range gives the same curve as the full grid
  nd2 <- nd[5:10, ]
  expect_equal(predict(fit, newdata = nd2, re.form = NA), pop[5:10])

  # in sample, and through simulate()
  expect_equal(length(fitted(fit)), nrow(d))
  s <- simulate(fit, nsim = 2, seed = 5)
  expect_identical(dim(s), c(nrow(d), 2L))

  # eval_dpars(b = NULL) drops the random-effect contribution, and a
  # penalized block IS a random effect: what is left is the null space
  ev0 <- frmtmb::eval_dpars(fit, b = NULL)
  expect_true(all(is.finite(ev0[["y"]]$mu)))
  expect_false(isTRUE(all.equal(ev0[["y"]]$mu,
                                frmtmb::eval_dpars(fit)[["y"]]$mu)))
})

test_that("frm_lp_basis() reaches a ps() block's coefficients", {
  skip_on_cran()
  set.seed(7073)
  n_id <- 25
  d <- expand.grid(t = seq(0, 1, length.out = 8), id = factor(1:n_id))
  sh <- rnorm(n_id, 0, 0.05)[d$id]
  d$y <- 2 + sin(2 * pi * (d$t + sh)) + rnorm(nrow(d), 0, 0.2)
  fit <- suppressWarnings(
    frm(bf(y ~ lev + ps(t + shift, k = 8, pad = 0.4),
           lev ~ 1, shift ~ 0 + (1 | id), nl = TRUE), d, gaussian()))
  nd <- data.frame(t = seq(0.1, 0.9, length.out = 9),
                   id = factor(1, levels = levels(d$id)))
  lb <- frm_lp_basis(fit, newdata = nd, re.form = NA)
  expect_equal(lb$eta, unname(predict(fit, newdata = nd, re.form = NA)))
  # the block's own coefficients are in the design
  pt <- fit$frame$linpreds[["y.mu"]]$ps_terms[[1]]
  expect_identical(ncol(lb$A), 1L + pt$n_fixed + pt$n_pen)
  expect_true(any(grepl("ps\\(", lb$coef_names)))

  # against a central difference in each coefficient
  jc <- frm_joint_cov(fit)
  comp <- jc$names
  idx <- unlist(lapply(unique(comp), function(cp) seq_len(sum(comp == cp))))
  worst <- 0
  for (k in seq_along(lb$coef_pos)) {
    pos <- lb$coef_pos[k]
    h <- 1e-6
    fp <- fit; fp$estimates[[comp[pos]]][idx[pos]] <-
      fp$estimates[[comp[pos]]][idx[pos]] + h
    fm <- fit; fm$estimates[[comp[pos]]][idx[pos]] <-
      fm$estimates[[comp[pos]]][idx[pos]] - h
    fd <- (predict(fp, newdata = nd, re.form = NA) -
             predict(fm, newdata = nd, re.form = NA)) / (2 * h)
    worst <- max(worst, max(abs(fd - lb$A[, k])))
  }
  expect_lt(worst, 1e-5)
})

test_that("ps() refuses the fitting options it cannot answer", {
  skip_on_cran()
  set.seed(7074)
  d <- data.frame(t = runif(120), id = factor(rep(1:12, 10)))
  d$y <- 2 + sin(2 * pi * d$t) + rnorm(120, 0, 0.2)
  f <- bf(y ~ lev + ps(t, k = 8), lev ~ 1 + (1 | id), nl = TRUE)
  expect_error(frm(f, d, gaussian(), REML = TRUE), "REML = TRUE cannot")
  expect_error(frm(f, d, gaussian(), quadrature = TRUE),
               "quadrature = TRUE cannot")
  expect_error(frm(f, d, gaussian(), control = frmtmb_control(profile = TRUE)),
               "profile = TRUE\\) cannot")
  # the importance correction refuses every nonlinear predictor, so a
  # ps() block inside one is refused by that guard rather than a new one
  expect_error(frm(f, d, gaussian(), importance = 20),
               "cannot correct a nonlinear predictor")
  # mvbf
  d$y2 <- d$y + rnorm(120, 0, 0.1)
  expect_error(
    frm(mvbf(bf(y ~ lev + ps(t, k = 8), lev ~ 1, nl = TRUE), bf(y2 ~ t)),
        d, gaussian()),
    "multivariate")
})

test_that("predict() at newdata past the knot span says so", {
  skip_on_cran()
  set.seed(7077)
  n_id <- 30
  d <- expand.grid(t = seq(0, 1, length.out = 8), id = factor(1:n_id))
  sh <- rnorm(n_id, 0, 0.05)[d$id]
  d$y <- 2 + sin(2 * pi * (d$t + sh)) + rnorm(nrow(d), 0, 0.15)
  fit <- suppressWarnings(
    frm(bf(y ~ lev + ps(t + shift, k = 8, pad = 0.3),
           lev ~ 1, shift ~ 0 + (1 | id), nl = TRUE), d, gaussian()))
  pt <- fit$frame$linpreds[["y.mu"]]$ps_terms[[1]]
  span <- pt$knot_range

  # inside the span: silent, and that is the ordinary case
  nd_in <- data.frame(t = seq(0.1, 0.9, length.out = 9),
                      id = factor(1, levels = levels(d$id)))
  expect_no_warning(p_in <- predict(fit, newdata = nd_in, re.form = NA))

  # outside it: the basis is a partial sum and then exactly zero, so the
  # curve decays and the prediction bends to the rest of the body. The
  # arithmetic is right and silence about it is the defect.
  nd_out <- data.frame(t = c(0.5, span[2] + 0.5, span[2] + 2),
                       id = factor(1, levels = levels(d$id)))
  expect_warning(p_out <- predict(fit, newdata = nd_out, re.form = NA),
                 "outside the frozen knot span")
  # the message carries the span itself, which is what a reader needs
  w <- tryCatch(predict(fit, newdata = nd_out, re.form = NA),
                warning = function(e) conditionMessage(e))
  expect_true(grepl(format(span[1], digits = 4), w, fixed = TRUE))
  expect_true(grepl(format(span[2], digits = 4), w, fixed = TRUE))
  expect_true(grepl("2 of 3", w, fixed = TRUE))

  # past the last knot the curve reads as exactly zero, so the
  # prediction is the rest of the body and nothing else
  far <- data.frame(t = 50, id = factor(1, levels = levels(d$id)))
  expect_equal(suppressWarnings(predict(fit, newdata = far,
                                        re.form = NA))[[1]],
               fixef(fit)$lev[[1]], tolerance = 1e-8)

  # in sample it stays quiet: the fit-end report has already said it
  expect_no_warning(predict(fit))
  expect_no_warning(fitted(fit))
})

test_that("a ps() fit reports the knot span it left", {
  skip_on_cran()
  set.seed(7075)
  n_id <- 30
  d <- expand.grid(t = seq(0, 1, length.out = 8), id = factor(1:n_id))
  sh <- rnorm(n_id, 0, 0.25)[d$id]
  d$y <- 2 + sin(2 * pi * (d$t + sh)) + rnorm(nrow(d), 0, 0.15)
  # pad = 0 puts the knots on the data range alone, so the fitted shifts
  # push rows outside it
  expect_warning(
    frm(bf(y ~ lev + ps(t + shift, k = 8, pad = 0), lev ~ 1,
           shift ~ 0 + (1 | id), nl = TRUE), d, gaussian()),
    "outside the knot span")
})

test_that("the SMOCC warped-growth likelihood is an identity", {
  skip_on_cran()
  skip_if_not_installed("brokenstick")
  # The model of D'Alessandro, Thoresen and Sorensen (2026),
  # arXiv:2603.11728, section 4:
  #   hgt = b0 + b1 sex + u1 + exp(b2 sex) f(age + b3 GA + u2) + e
  # A subsample keeps the test quick; the whole 200 subjects are fitted
  # in vignette("case-studies").
  raw <- brokenstick::smocc_200
  d <- data.frame(id = factor(raw$id), age = 52.1775 * raw$age,
                  sex = as.numeric(raw$sex == "male"), ga = raw$ga - 40,
                  hgt = raw$hgt)
  d <- d[!is.na(d$hgt), ]
  d <- d[d$id %in% levels(d$id)[1:60], ]
  d$id <- droplevels(d$id)

  fit <- suppressWarnings(frm(
    bf(hgt ~ int + exp(amp) * ps(age + shift, k = 15, pad = 0.25),
       int ~ sex + (1 | id), amp ~ 0 + sex, shift ~ 0 + ga + (1 | id),
       nl = TRUE),
    d, gaussian(), start = list(beta = c(68, 2, 0, 1, 0))))

  pt <- fit$frame$linpreds[["hgt.mu"]]$ps_terms[[1]]
  bk <- fit$frame$re_blocks[[pt$block_id]]
  bn <- names(fit$frame$par_template$beta)

  # the paper's joint density, written out, with the basis taken from
  # splines::splineDesign rather than from frmtmb's own construction
  ref_nll <- function(p) {
    gamma <- as.vector(pt$U0 %*% p$beta[pt$beta_idx] +
                         pt$Us %*% p$b[bk$b_idx])
    arg <- d$age + p$beta[match("shift_ga", bn)] * d$ga +
      p$b[fit$frame$re_blocks[[2]]$b_idx][as.integer(d$id)]
    f <- as.vector(splines::splineDesign(pt$knots, arg, ord = pt$ord,
                                         outer.ok = TRUE) %*% gamma)
    mu <- p$beta[match("int_(Intercept)", bn)] +
      p$beta[match("int_sex", bn)] * d$sex +
      p$b[fit$frame$re_blocks[[1]]$b_idx][as.integer(d$id)] +
      exp(p$beta[match("amp_sex", bn)] * d$sex) * f
    -(sum(stats::dnorm(d$hgt, mu, exp(p$betad[1]), log = TRUE)) +
        sum(stats::dnorm(p$b[fit$frame$re_blocks[[1]]$b_idx], 0,
                         exp(p$theta[fit$frame$re_blocks[[1]]$theta_idx]),
                         log = TRUE)) +
        sum(stats::dnorm(p$b[fit$frame$re_blocks[[2]]$b_idx], 0,
                         exp(p$theta[fit$frame$re_blocks[[2]]$theta_idx]),
                         log = TRUE)) +
        sum(stats::dnorm(p$b[bk$b_idx], 0, exp(p$theta[bk$theta_idx]),
                         log = TRUE)))
  }

  nll <- frmtmb:::build_objective(fit$frame)
  expect_equal(nll(fit$estimates), ref_nll(fit$estimates), tolerance = 1e-11)

  set.seed(7076)
  p2 <- fit$estimates
  p2$beta <- p2$beta + rnorm(length(p2$beta), 0, 0.15)
  p2$betad <- p2$betad + rnorm(length(p2$betad), 0, 0.15)
  p2$b <- p2$b + rnorm(length(p2$b), 0, 0.3)
  p2$theta <- p2$theta + rnorm(length(p2$theta), 0, 0.15)
  expect_equal(nll(p2), ref_nll(p2), tolerance = 1e-11)
})
