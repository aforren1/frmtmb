# Multi-membership random effects: (x | mm(g1, g2)).
#
# The design claim is narrow and the tests pin it down: mm() changes the
# Z MATRIX and nothing else. The block is an ordinary us/diag block over
# the pooled level set, so the likelihood, the Laplace approximation and
# every post-fit method are the single-membership ones. What has to be
# checked is therefore the Z construction (against brms, exactly), that
# the fit it produces is the one a hand-rolled RTMB objective over the
# same Z gives, and that prediction rebuilds the same rows on newdata.

mm_data <- function(seed = 4711, n = 150) {
  set.seed(seed)
  d <- data.frame(
    x = rnorm(n),
    g1 = factor(sample(letters[1:6], n, TRUE)),
    g2 = factor(sample(letters[3:9], n, TRUE)),
    w1 = runif(n, 0.2, 1),
    w2 = runif(n, 0.2, 2),
    c1 = rnorm(n),
    c2 = rnorm(n),
    fc = factor(sample(c("p", "q"), n, TRUE))
  )
  d$one <- 1
  d$zero <- 0
  # the DGP is the multi-membership model itself: a pooled intercept and
  # a pooled slope, each observation seeing the average of its two
  # members, so the mmc() variance is away from the boundary
  lv <- letters[1:9]
  u <- stats::setNames(rnorm(9, 0, 0.8), lv)
  v <- stats::setNames(rnorm(9, 0, 0.5), lv)
  i1 <- as.character(d$g1)
  i2 <- as.character(d$g2)
  d$y <- 1 + 0.5 * d$x + 0.5 * (u[i1] + u[i2]) +
    0.5 * (v[i1] * d$c1 + v[i2] * d$c2) + rnorm(n, 0, 0.6)
  d
}

# ------------------------------------------------------------ parsing

test_that("mm() is read as a grouping call, not a design term", {
  sp <- parse_spec(bf(y ~ x + (1 | mm(g1, g2))) + gaussian())
  rt <- sp$responses$y$dpars$mu$re[[1L]]
  expect_identical(rt$covstruct, "us")
  expect_identical(rt$mm$gvars, c("g1", "g2"))
  expect_null(rt$mm$weights_expr)
  expect_true(rt$mm$scale)
  expect_identical(deparse1(rt$mm$lhs), "1")
  expect_length(rt$mm$mmc, 0L)
  # the bar keeps the user's spelling, which is the term label
  expect_identical(deparse1(rt$bar), "1 | mm(g1, g2)")
})

test_that("mmc() terms are split off the left of the bar", {
  sp <- parse_spec(bf(y ~ (1 + c1 + mmc(c1, c2) | mm(g1, g2))) +
                     gaussian())
  mm <- sp$responses$y$dpars$mu$re[[1L]]$mm
  expect_identical(deparse1(mm$lhs), "1 + c1")
  expect_length(mm$mmc, 1L)
  expect_identical(mm$mmc[[1L]]$label, "mmc(c1, c2)")
  # an mmc()-only left side keeps the implicit intercept, as (x | g) does
  sp0 <- parse_spec(bf(y ~ (mmc(c1, c2) | mm(g1, g2))) + gaussian())
  expect_identical(deparse1(sp0$responses$y$dpars$mu$re[[1L]]$mm$lhs), "1")
  sp1 <- parse_spec(bf(y ~ (0 + mmc(c1, c2) | mm(g1, g2))) + gaussian())
  expect_identical(deparse1(sp1$responses$y$dpars$mu$re[[1L]]$mm$lhs), "0")
})

test_that("a column called mm is not mistaken for the special", {
  dd <- mm_data(n = 60)
  dd$mm <- rnorm(nrow(dd))
  fr <- frm(bf(y ~ mm + (1 | g1)) + gaussian(), data = dd,
            dry_run = "frame")
  expect_true("mm" %in% colnames(fr$linpreds[["y.mu"]]$X))
})

# ------------------------------------------------------- Z construction

test_that("the level set is the members' levels pooled in order", {
  dd <- mm_data(n = 80)
  fr <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd,
            dry_run = "frame")
  bk <- fr$re_blocks[[1L]]
  # brms frame_re(): unique(ulapply(groups, extract_levels)), so the
  # members' own level orders are concatenated and deduplicated
  expect_identical(bk$levels, unique(c(levels(dd$g1), levels(dd$g2))))
  expect_identical(bk$n_levels, length(bk$levels))
  expect_identical(bk$covstruct, "us")
  expect_identical(bk$group_name, "mm(g1, g2)")
  # one variance parameter: the block is an ordinary scalar block
  expect_length(bk$theta_idx, 1L)
})

test_that("equal default weights put 1/J on each member's column", {
  dd <- mm_data(n = 60)
  fr <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd,
            dry_run = "frame")
  bk <- fr$re_blocks[[1L]]
  Z <- as.matrix(fr$linpreds[["y.mu"]]$Z)
  expect_equal(unname(rowSums(Z)), rep(1, nrow(dd)), tolerance = 1e-12)
  j1 <- match(as.character(dd$g1), bk$levels)
  j2 <- match(as.character(dd$g2), bk$levels)
  Zb <- matrix(0, nrow(dd), ncol(Z))
  Zb[cbind(seq_len(nrow(dd)), j1)] <- Zb[cbind(seq_len(nrow(dd)), j1)] + 0.5
  Zb[cbind(seq_len(nrow(dd)), j2)] <- Zb[cbind(seq_len(nrow(dd)), j2)] + 0.5
  expect_lt(max(abs(Z - Zb)), 1e-14)
})

test_that("scale = TRUE divides supplied weights by their row sums", {
  dd <- mm_data(n = 60)
  frs <- frm(bf(y ~ (1 | mm(g1, g2, weights = cbind(w1, w2)))) +
               gaussian(), data = dd, dry_run = "frame")
  Zs <- as.matrix(frs$linpreds[["y.mu"]]$Z)
  expect_equal(unname(rowSums(Zs)), rep(1, nrow(dd)), tolerance = 1e-12)
  frr <- frm(bf(y ~ (1 | mm(g1, g2, weights = cbind(w1, w2),
                            scale = FALSE))) + gaussian(),
             data = dd, dry_run = "frame")
  Zr <- as.matrix(frr$linpreds[["y.mu"]]$Z)
  expect_equal(unname(rowSums(Zr)), dd$w1 + dd$w2, tolerance = 1e-12)
})

test_that("mmc() gives one coefficient with member-specific values", {
  dd <- mm_data(n = 60)
  fr <- frm(bf(y ~ (1 + mmc(c1, c2) | mm(g1, g2))) + gaussian(),
            data = dd, dry_run = "frame")
  bk <- fr$re_blocks[[1L]]
  expect_identical(bk$dim, 2L)
  expect_identical(bk$cnms, c("(Intercept)", "mmc(c1, c2)"))
  Z <- as.matrix(fr$linpreds[["y.mu"]]$Z)
  j1 <- match(as.character(dd$g1), bk$levels)
  j2 <- match(as.character(dd$g2), bk$levels)
  n <- nrow(dd)
  Zb <- matrix(0, n, ncol(Z))
  add <- function(M, j, col, v) {
    idx <- cbind(seq_len(n), (j - 1L) * 2L + col)
    M[idx] <- M[idx] + v
    M
  }
  Zb <- add(Zb, j1, 1L, 0.5)
  Zb <- add(Zb, j2, 1L, 0.5)
  Zb <- add(Zb, j1, 2L, 0.5 * dd$c1)   # member 1 uses mmc()'s 1st arg
  Zb <- add(Zb, j2, 2L, 0.5 * dd$c2)   # member 2 uses its 2nd
  expect_lt(max(abs(Z - Zb)), 1e-14)
})

test_that("a degenerate mm(g, g) with weights c(1, 0) IS (1 | g)", {
  dd <- mm_data(n = 120)
  fa <- frm(bf(y ~ x + (1 | mm(g1, g1, weights = cbind(one, zero)))) +
              gaussian(), data = dd)
  fb <- frm(bf(y ~ x + (1 | g1)) + gaussian(), data = dd)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fb)),
               tolerance = 1e-10)
  expect_equal(unlist(fixef(fa)), unlist(fixef(fb)), tolerance = 1e-8)
  expect_equal(unname(unlist(ranef(fa))), unname(unlist(ranef(fb))),
               tolerance = 1e-8)
  # and so is mm(g, g) with the default 1/2 + 1/2 weights, which lands
  # the whole weight on the same column
  fc <- frm(bf(y ~ x + (1 | mm(g1, g1))) + gaussian(), data = dd)
  expect_equal(as.numeric(logLik(fc)), as.numeric(logLik(fb)),
               tolerance = 1e-10)
})

# ------------------------------------------ brms structural agreement

test_that("mm() designs match brms make_standata J/W/Z arrays", {
  skip_unless_brms()
  dd <- mm_data(n = 150)
  n <- nrow(dd)
  # brms parses mm()/mmc() by evaluating them in the formula environment
  mm <- brms::mm
  mmc <- brms::mmc

  check <- function(brmform, frmform, dim) {
    sd <- brms_standata(brmform, data = dd, family = gaussian())
    fr <- frm(frmform, data = dd, dry_run = "frame")
    bk <- fr$re_blocks[[1L]]
    expect_identical(bk$n_levels, as.integer(sd$N_1))
    expect_identical(bk$dim, as.integer(sd$M_1))
    Z <- as.matrix(fr$linpreds[["y.mu"]]$Z)
    # brms keeps the membership implicit: one level index J_1_k, one
    # weight W_1_k and one covariate column Z_1_c_k per member, summed
    # into the linear predictor by the Stan program. Our Z is that sum.
    Zb <- matrix(0, n, ncol(Z))
    for (cc in seq_len(dim)) {
      for (k in 1:2) {
        idx <- cbind(seq_len(n),
                     (sd[[paste0("J_1_", k)]] - 1L) * bk$dim + cc)
        Zb[idx] <- Zb[idx] +
          sd[[paste0("W_1_", k)]] * sd[[paste0("Z_1_", cc, "_", k)]]
      }
    }
    expect_lt(max(abs(Z - Zb)), 1e-12)
  }

  check(brms::bf(y ~ x + (1 | mm(g1, g2))),
        bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), 1L)
  check(brms::bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2)))),
        bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2)))) +
          gaussian(), 1L)
  check(brms::bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2),
                                 scale = FALSE))),
        bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2),
                           scale = FALSE))) + gaussian(), 1L)
  check(brms::bf(y ~ x + (1 + mmc(c1, c2) | mm(g1, g2))),
        bf(y ~ x + (1 + mmc(c1, c2) | mm(g1, g2))) + gaussian(), 2L)
  check(brms::bf(y ~ x + (0 + mmc(c1, c2) | mm(g1, g2))),
        bf(y ~ x + (0 + mmc(c1, c2) | mm(g1, g2))) + gaussian(), 1L)
})

test_that("the pooled level order matches brms's own", {
  skip_unless_brms()
  mm <- brms::mm
  # deliberately non-alphabetical factor levels, and a member level that
  # only one of the two factors carries
  dd <- data.frame(
    y = rnorm(24),
    g1 = factor(rep(c("z", "y", "x"), 8), levels = c("z", "y", "x")),
    g2 = factor(rep(c("b", "c", "a"), 8)))
  sd <- brms_standata(brms::bf(y ~ (1 | mm(g1, g2))), data = dd,
                      family = gaussian())
  fr <- frm(bf(y ~ (1 | mm(g1, g2))) + gaussian(), data = dd,
            dry_run = "frame")
  bk <- fr$re_blocks[[1L]]
  expect_identical(bk$levels, c("z", "y", "x", "a", "b", "c"))
  expect_identical(match(as.character(dd$g1), bk$levels),
                   as.integer(sd$J_1_1))
  expect_identical(match(as.character(dd$g2), bk$levels),
                   as.integer(sd$J_1_2))
})

# ---------------------------------------------- hand-rolled likelihood

# A plain RTMB objective over the SAME weighted Z, with its own
# parameterization of the covariance. It shares no code with the package
# beyond the design matrix, so agreement at the optimum says the
# registry nll and the Laplace machinery read an mm block exactly as
# they read any other block.
mm_reference_loglik <- function(Zmat, y, X, nlev, dim) {
  dat <- list(y = y, X = X, Z = Zmat)
  par <- list(beta = rep(0, ncol(X)), logsd = 0,
              logtheta = numeric(dim),
              cortheta = numeric(dim * (dim - 1L) / 2L),
              b = numeric(nlev * dim))
  f <- function(p) {
    RTMB::getAll(p, dat)
    nll <- 0
    if (dim == 1L) {
      nll <- nll - sum(RTMB::dnorm(b, 0, exp(logtheta[1]), log = TRUE))
    } else {
      # written out rather than built with matrix helpers, which strip
      # the AD class off their arguments
      s1 <- exp(logtheta[1]); s2 <- exp(logtheta[2])
      rho <- cortheta[1] / sqrt(1 + cortheta[1]^2)
      i1 <- seq(1, by = 2, length.out = nlev)
      z1 <- b[i1] / s1
      z2 <- b[i1 + 1L] / s2
      om <- 1 - rho^2
      nll <- nll + nlev * (log(2 * pi) + log(s1) + log(s2) +
                             0.5 * log(om)) +
        sum(z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * om)
    }
    eta <- as.vector(X %*% beta) + as.vector(Z %*% b)
    nll - sum(RTMB::dnorm(y, eta, exp(logsd), log = TRUE))
  }
  obj <- RTMB::MakeADFun(f, par, random = "b", silent = TRUE)
  op <- suppressWarnings(
    stats::nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 5000, iter.max = 5000,
                                 rel.tol = 1e-14, x.tol = 1e-12)))
  -op$objective
}

test_that("an equal-weight mm fit equals a hand-rolled RTMB objective", {
  dd <- mm_data(n = 150)
  form <- bf(y ~ x + (1 | mm(g1, g2))) + gaussian()
  fr <- frm(form, data = dd, dry_run = "frame")
  fit <- frm(form, data = dd)
  ref <- mm_reference_loglik(as.matrix(fr$linpreds[["y.mu"]]$Z), dd$y,
                             stats::model.matrix(~ x, dd),
                             fr$re_blocks[[1L]]$n_levels, 1L)
  expect_equal(as.numeric(logLik(fit)), ref, tolerance = 1e-8)
})

test_that("an mmc() random slope equals a hand-rolled RTMB objective", {
  dd <- mm_data(n = 150)
  form <- bf(y ~ x + (1 + mmc(c1, c2) | mm(g1, g2))) + gaussian()
  fr <- frm(form, data = dd, dry_run = "frame")
  fit <- frm(form, data = dd)
  ref <- mm_reference_loglik(as.matrix(fr$linpreds[["y.mu"]]$Z), dd$y,
                             stats::model.matrix(~ x, dd),
                             fr$re_blocks[[1L]]$n_levels, 2L)
  expect_equal(as.numeric(logLik(fit)), ref, tolerance = 1e-6)
})

test_that("weighted memberships equal a hand-rolled RTMB objective", {
  dd <- mm_data(n = 150)
  form <- bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2)))) +
    gaussian()
  fr <- frm(form, data = dd, dry_run = "frame")
  fit <- frm(form, data = dd)
  ref <- mm_reference_loglik(as.matrix(fr$linpreds[["y.mu"]]$Z), dd$y,
                             stats::model.matrix(~ x, dd),
                             fr$re_blocks[[1L]]$n_levels, 1L)
  expect_equal(as.numeric(logLik(fit)), ref, tolerance = 1e-8)
})

# ------------------------------------------------------------ post-fit

test_that("ranef, VarCorr, ngrps and simulate see an ordinary block", {
  dd <- mm_data(n = 120)
  fit <- frm(bf(y ~ x + (1 + mmc(c1, c2) | mm(g1, g2))) + gaussian(),
             data = dd)
  re <- ranef(fit)
  expect_length(re, 1L)
  expect_identical(dim(re[[1L]]), c(9L, 2L))
  expect_identical(rownames(re[[1L]]),
                   unique(c(levels(dd$g1), levels(dd$g2))))
  expect_identical(colnames(re[[1L]]),
                   c("(Intercept)", "mmc(c1, c2)"))
  vc <- VarCorr(fit)
  expect_identical(dim(vc[[1L]]), c(2L, 2L))
  expect_identical(dimnames(vc[[1L]])[[1L]],
                   c("(Intercept)", "mmc(c1, c2)"))
  expect_identical(unname(ngrps(fit)), 9L)
  sm <- simulate(fit, nsim = 3)
  expect_identical(dim(as.matrix(sm)), c(nrow(dd), 3L))
  expect_false(anyNA(sm))
  expect_false(anyNA(residuals(fit)))
})

test_that("REML runs over a multi-membership block", {
  dd <- mm_data(n = 120)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd,
             REML = TRUE)
  expect_true(is.finite(as.numeric(logLik(fit))))
  # REML integrates the fixed effects only, so the block is untouched
  expect_identical(fit$frame$re_blocks[[1L]]$covstruct, "us")
  expect_gt(as.numeric(VarCorr(fit)[[1L]]), 0)
})

test_that("diag() and || give uncorrelated multi-membership effects", {
  dd <- mm_data(n = 120)
  fd <- frm(bf(y ~ x + diag(1 + c1 | mm(g1, g2))) + gaussian(),
            data = dd)
  expect_identical(fd$frame$re_blocks[[1L]]$covstruct, "diag")
  expect_length(fd$frame$re_blocks[[1L]]$theta_idx, 2L)
  # (x || g) expands the way lme4 does, so it becomes two diag blocks
  fb <- frm(bf(y ~ x + (1 + c1 || mm(g1, g2))) + gaussian(), data = dd)
  expect_length(fb$frame$re_blocks, 2L)
  expect_equal(as.numeric(logLik(fd)), as.numeric(logLik(fb)),
               tolerance = 1e-6)
})

# --------------------------------------------------------- prediction

test_that("predict on the training rows reproduces fitted()", {
  dd <- mm_data(n = 120)
  for (form in list(
    bf(y ~ x + (1 | mm(g1, g2))) + gaussian(),
    bf(y ~ x + (1 + mmc(c1, c2) | mm(g1, g2))) + gaussian(),
    bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2)))) + gaussian()
  )) {
    fit <- frm(form, data = dd)
    expect_lt(max(abs(predict(fit, newdata = dd) - fitted(fit))), 1e-10)
  }
})

test_that("a partially-new member set predicts its known members", {
  dd <- mm_data(n = 120)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd)
  bk <- fit$frame$re_blocks[[1L]]
  nd <- dd[1:6, ]
  nd$g2 <- factor(c("zz", "zz", as.character(nd$g2[3:6])),
                  levels = c("zz", levels(dd$g2)))
  expect_error(predict(fit, newdata = nd), "New levels")
  p <- predict(fit, newdata = nd, allow_new_levels = TRUE)
  bhat <- ranef(fit)[[1L]][, 1L]
  j1 <- match(as.character(nd$g1), bk$levels)
  j2 <- match(as.character(nd$g2), bk$levels)
  manual <- as.numeric(stats::model.matrix(~ x, nd) %*%
                         unlist(fixef(fit)$mu)) +
    0.5 * bhat[j1] + 0.5 * ifelse(is.na(j2), 0, bhat[j2])
  # the unknown member contributes the population value (zero); the
  # known one still contributes its own weighted effect
  expect_lt(max(abs(p - manual)), 1e-10)
  se <- predict(fit, newdata = nd, allow_new_levels = TRUE,
                se.fit = TRUE)$se.fit
  expect_true(all(is.finite(se)))
})

test_that("one new level reached through both members is ONE draw", {
  # The members of a row are independent draws of the block only when
  # they name DIFFERENT unseen levels. The SAME unseen label through
  # both members is one draw, so the weights add before the variance is
  # taken: (w1 + w2)^2 S, not w1^2 S + w2^2 S. Rows whose members are
  # ALL unseen carry no b column at all, so the difference from the
  # population-level prediction is exactly the new-level term.
  dd <- mm_data(n = 120)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd)
  S <- as.numeric(VarCorr(fit)[[1L]])

  nd <- dd[rep(1L, 4L), ]          # one row four times: eta is common
  nd$g1 <- factor(rep("zz", 4L), levels = c("zz", levels(dd$g1)))
  # rows 1-2 name the same unseen level twice, rows 3-4 two distinct ones
  nd$g2 <- factor(c("zz", "zz", "qq", "qq"),
                  levels = c("zz", "qq", levels(dd$g2)))

  p <- predict(fit, newdata = nd, allow_new_levels = TRUE, se.fit = TRUE)
  base <- predict(fit, newdata = nd, se.fit = TRUE, re.form = NA)
  extra <- unname(p$se.fit^2 - base$se.fit^2)
  expect_equal(extra[1:2], rep((0.5 + 0.5)^2 * S, 2), tolerance = 1e-8)
  expect_equal(extra[3:4], rep(2 * 0.5^2 * S, 2), tolerance = 1e-8)
  # the two conventions differ by a factor of two, which is the point
  expect_equal(extra[1] / extra[3], 2, tolerance = 1e-8)
})

test_that("the expected-response path groups new levels the same way", {
  # A lognormal mean is not the mu dpar, so this is the OTHER standard
  # error path (predict_mean_se). All four rows share one eta, hence one
  # dpar gradient, so the ratio of the new-level terms is the ratio of
  # the two conventions and nothing else.
  dd <- mm_data(n = 120)
  dd$yp <- exp(dd$y / 4)
  fit <- frm(bf(yp ~ x + (1 | mm(g1, g2))) + lognormal(), data = dd)

  nd <- dd[rep(1L, 4L), ]
  nd$g1 <- factor(rep("zz", 4L), levels = c("zz", levels(dd$g1)))
  nd$g2 <- factor(c("zz", "zz", "qq", "qq"),
                  levels = c("zz", "qq", levels(dd$g2)))

  p <- predict(fit, newdata = nd, type = "response",
               allow_new_levels = TRUE, se.fit = TRUE)
  base <- predict(fit, newdata = nd, type = "response", se.fit = TRUE,
                  re.form = NA)
  extra <- unname(p$se.fit^2 - base$se.fit^2)
  expect_equal(extra[1], extra[2], tolerance = 1e-10)
  expect_equal(extra[3], extra[4], tolerance = 1e-10)
  expect_equal(extra[1] / extra[3], 2, tolerance = 1e-6)
})

test_that("newdata that drops a factor level keeps the block's columns", {
  dd <- mm_data(n = 120)
  fit <- frm(bf(y ~ x + diag(1 + fc | mm(g1, g2))) + gaussian(),
             data = dd)
  rows <- which(dd$fc == "p")[1:5]
  nd <- dd[rows, ]
  nd$fc <- droplevels(nd$fc)   # xlev has to restore the dropped column
  expect_lt(max(abs(predict(fit, newdata = nd) - fitted(fit)[rows])),
            1e-10)
})

# --------------------------------------------------------- refusals

test_that("mm() refuses the arguments it does not implement", {
  dd <- mm_data(n = 60)
  f <- function(form) frm(bf(form) + gaussian(), data = dd,
                          dry_run = "frame")
  expect_error(f(y ~ (1 | mm(g1))), "at least two membership")
  expect_error(f(y ~ (1 | mm(factor(g1), g2))), "bare column name")
  expect_error(f(y ~ (1 | mm(g1, g2, cor = FALSE))), "is not supported")
  expect_error(f(y ~ (1 | mm(g1, g2, wieghts = w1))),
               "unknown argument")
  expect_error(f(y ~ (1 | mm(g1, g2, weights = cbind(w1, w2),
                             scale = 3))), "must be TRUE or FALSE")
})

test_that("mm() refuses structures the weighted design cannot carry", {
  dd <- mm_data(n = 60)
  f <- function(form) frm(bf(form) + gaussian(), data = dd,
                          dry_run = "frame")
  expect_error(f(y ~ cs(1 + c1 | mm(g1, g2))),
               "default \\(us\\) and diag structures only")
  expect_error(f(y ~ ar1(0 + fc | mm(g1, g2))),
               "default \\(us\\) and diag structures only")
  expect_error(f(y ~ (1 | gr(mm(g1, g2), cov = A))),
               "relationship matrix indexes one level")
  expect_error(f(y ~ (1 | q | mm(g1, g2)) + (0 + x | q | mm(g1, g2))),
               "cannot share an \\|ID\\| key")
})

test_that("mmc() is refused wherever it has no member to index", {
  dd <- mm_data(n = 60)
  f <- function(form) frm(bf(form) + gaussian(), data = dd,
                          dry_run = "frame")
  expect_error(f(y ~ (1 + mmc(c1) | mm(g1, g2))),
               "one variable per membership variable")
  expect_error(f(y ~ (1 + mmc(a = c1, c2) | mm(g1, g2))),
               "unnamed variables")
  expect_error(f(y ~ (1 + mmc(c1, c2):x | mm(g1, g2))),
               "term of its own")
  expect_error(f(y ~ (1 + mmc(c1, c2) | g1)),
               "one covariate value per MEMBER")
  expect_error(f(y ~ mmc(c1, c2)),
               "not a population-level predictor")
  expect_error(f(y ~ mm(g1, g2)),
               "not a population-level predictor")
  expect_error(f(y ~ (0 + mmc(fc, c2) | mm(g1, g2))),
               "requires numeric variables")
})

test_that("mm() weights are checked for shape, sign and scale", {
  dd <- mm_data(n = 60)
  dd$nw <- -dd$w1
  f <- function(form) frm(bf(form) + gaussian(), data = dd,
                          dry_run = "frame")
  expect_error(f(y ~ (1 | mm(g1, g2, weights = w1))),
               "one column per membership variable")
  expect_error(f(y ~ (1 | mm(g1, g2, weights = cbind(nw, w2)))),
               "cannot scale negative weights")
  expect_error(f(y ~ (1 | mm(g1, g2, weights = cbind(zero, zero)))),
               "sum to zero")
  # scale = FALSE takes them as they are, negatives included
  expect_s3_class(f(y ~ (1 | mm(g1, g2, weights = cbind(nw, w2),
                                scale = FALSE))), "frmtmb_frame")
  expect_error(f(y ~ (0 | mm(g1, g2))), "at least one coefficient")
})

test_that("mm blocks coexist with the rest of the grammar", {
  dd <- mm_data(n = 120)
  # a plain factor block alongside the membership block: two blocks,
  # each with its own levels and its own variance
  fr <- frm(bf(y ~ x + (1 | mm(g1, g2)) + (1 | fc)) + gaussian(),
            data = dd, dry_run = "frame")
  expect_length(fr$re_blocks, 2L)
  expect_identical(fr$re_blocks[[1L]]$group_name, "mm(g1, g2)")
  expect_identical(fr$re_blocks[[2L]]$levels, levels(dd$fc))
  # a membership block on a dpar of its own
  frs <- frm(bf(y ~ x, sigma ~ (1 | mm(g1, g2))) + gaussian(),
             data = dd, dry_run = "frame")
  expect_identical(frs$re_blocks[[1L]]$dpar, "sigma")
  # sparse_x stores the fixed design differently and must not disturb it
  frd <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd,
             dry_run = "frame")
  frx <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd,
             dry_run = "frame",
             control = frmtmb_control(sparse_x = TRUE))
  expect_lt(max(abs(as.matrix(frx$linpreds[["y.mu"]]$Z) -
                      as.matrix(frd$linpreds[["y.mu"]]$Z))), 1e-14)
  # writing the members the other way round fits the same model with a
  # permuted level order
  fa <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd)
  fb <- frm(bf(y ~ x + (1 | mm(g2, g1))) + gaussian(), data = dd)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fb)),
               tolerance = 1e-8)
  ra <- ranef(fa)[[1L]][, 1L]
  rb <- ranef(fb)[[1L]][, 1L]
  expect_equal(ra[order(names(ra))], rb[order(names(rb))],
               tolerance = 1e-6)
})

test_that("mm blocks are declared in the compatibility registry", {
  ft <- frm_compat_features()
  expect_true("mm()" %in% ft$name)
  expect_true("mmc()" %in% ft$name)
  pairs <- frm_compat()
  hit <- pairs[pairs$feature_a == "mm()" | pairs$feature_b == "mm()", ]
  expect_gt(nrow(hit), 0L)
  # the refusals the parser really raises must be declared as refusals
  refused <- c("gr_cov", "gr_prec", "|ID|", "ar1", "cs")
  for (f2 in refused) {
    row <- hit[hit$feature_a == f2 | hit$feature_b == f2, ]
    expect_identical(row$status, "refused")
  }
})
