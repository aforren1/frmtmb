# frmtmb_control(autoscale = TRUE): internal predictor standardization
# via the two-stage warm start (scaled pre-fit, exact back-transform,
# ordinary unscaled fit). ML is invariant under affine column
# reparameterization, so hand-standardized references are exact.

# Collect a fit together with every warning it raised.
fit_catching <- function(expr) {
  w <- character(0)
  fit <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warnings = w)
}

sim_badscale_glmm <- function(seed = 202, n_g = 30, n_per = 20) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_g), each = n_per))
  x1 <- stats::rnorm(n_g * n_per) * 1e6
  x2 <- stats::rnorm(n_g * n_per) * 1e-6
  b0 <- stats::rnorm(n_g, 0, 0.5)
  eta <- 0.3 + 4e-7 * x1 + 3e5 * x2 + b0[g]
  data.frame(y = stats::rpois(length(eta), exp(eta)),
             x1 = x1, x2 = x2, g = g)
}

test_that("autoscale reproduces the plain fit on a well-scaled model", {
  dd <- sim_pois_glmm()
  form <- bf(y ~ x + (1 | g)) + poisson()
  f0 <- frm(form, data = dd)
  fa <- frm(form, data = dd,
            control = frmtmb_control(autoscale = TRUE))
  expect_loglik_equal(fa, f0, tol = 1e-8)
  expect_vector_equal(fixef_by_dpar(fa)$mu, fixef_by_dpar(f0)$mu, tol = 1e-6)
  expect_vector_equal(fa$estimates$theta, f0$estimates$theta, tol = 1e-5)
})

test_that("autoscale is a silent no-op when nothing qualifies", {
  set.seed(303)
  dd <- data.frame(f = factor(rep(letters[1:3], 100)),
                   g = factor(rep(1:30, each = 10)))
  b0 <- stats::rnorm(30, 0, 0.5)
  dd$y <- stats::rpois(300, exp(0.5 + c(0, 0.3, -0.2)[dd$f] + b0[dd$g]))
  form <- bf(y ~ f + (1 | g)) + poisson()
  f0 <- frm(form, data = dd)
  expect_null(autoscale_plan(f0$frame))
  fa <- NULL
  expect_silent(fa <- frm(form, data = dd,
                          control = frmtmb_control(autoscale = TRUE)))
  # no qualifying column: exactly the single ordinary fit
  expect_identical(fa$opt$par, f0$opt$par)
  expect_null(fa$par_units)
})

test_that("autoscale rescues a badly scaled poisson GLMM", {
  db <- sim_badscale_glmm()
  form <- bf(y ~ x1 + x2 + (1 | g)) + poisson()

  auto <- fit_catching(frm(form, data = db,
                           control = frmtmb_control(autoscale = TRUE)))
  expect_s3_class(auto$fit, "frmtmb_fit")
  expect_identical(auto$warnings, character(0))
  expect_equal(auto$fit$opt$convergence, 0L)
  fa <- auto$fit

  # hand-standardized reference with the standardization undone by
  # hand: exact under ML reparameterization invariance
  m1 <- mean(db$x1); s1 <- stats::sd(db$x1)
  m2 <- mean(db$x2); s2 <- stats::sd(db$x2)
  db2 <- db
  db2$z1 <- (db$x1 - m1) / s1
  db2$z2 <- (db$x2 - m2) / s2
  fr <- frm(bf(y ~ z1 + z2 + (1 | g)) + poisson(), data = db2)
  expect_loglik_equal(fa, fr, tol = 1e-6)
  cz <- fixef_by_dpar(fr)$mu
  ref <- c(cz[["(Intercept)"]] - cz[["z1"]] * m1 / s1 -
             cz[["z2"]] * m2 / s2,
           cz[["z1"]] / s1, cz[["z2"]] / s2)
  expect_equal(unname(fixef_by_dpar(fa)$mu), ref, tolerance = 1e-6)
  expect_vector_equal(fa$estimates$theta, fr$estimates$theta, tol = 1e-4)

  # SEs transform by the same column scales (finite-difference Hessians
  # on both sides, so agreement is at step-size accuracy)
  se_a <- summary(fa)$coefficients$mu[, "Std. Error"]
  se_r <- summary(fr)$coefficients$mu[, "Std. Error"]
  expect_equal(unname(se_a[c("x1", "x2")]),
               unname(se_r[c("z1", "z2")]) / c(s1, s2),
               tolerance = 1e-3)
  # prediction SEs ride on the joint precision and match too
  pa <- frm_linpred(fa, se.fit = TRUE)
  pr <- frm_linpred(fr, se.fit = TRUE)
  expect_true(all(is.finite(pa$se.fit)))
  expect_vector_equal(pa$se.fit, pr$se.fit, tol = 1e-3)

  # the plain fit either matches or (as here, typically) fails where
  # autoscale did not
  plain <- fit_catching(frm(form, data = db))
  if (!inherits(plain$fit, "error") && !length(plain$warnings)) {
    expect_loglik_equal(plain$fit, fa, tol = 1e-6)
  } else {
    expect_gt(length(plain$warnings), 0)
  }

  # warm-started refit skips the pre-fit (template short-circuit)
  rf <- refit(fa, db$y)
  expect_loglik_equal(rf, fa, tol = 1e-6)
})

test_that("autoscale leaves smooth and mo() columns untouched", {
  skip_if_not_installed("mgcv")
  set.seed(404)
  n <- 300
  g <- factor(rep(1:30, each = 10))
  xbig <- stats::rnorm(n) * 1e6
  xs <- stats::runif(n)
  m <- factor(sample(1:4, n, replace = TRUE), ordered = TRUE)
  b0 <- stats::rnorm(30, 0, 0.4)
  mo_eff <- c(0, 0.3, 0.5, 0.9)[as.integer(m)]
  dd <- data.frame(xbig = xbig, xs = xs, m = m, g = g)
  dd$y <- 1 + 4e-7 * xbig + sin(2 * pi * xs) + mo_eff + b0[g] +
    stats::rnorm(n, 0, 0.5)

  form <- bf(y ~ xbig + mo(m) + s(xs) + (1 | g)) + gaussian()
  auto <- fit_catching(frm(form, data = dd,
                           control = frmtmb_control(autoscale = TRUE)))
  expect_identical(auto$warnings, character(0))
  fa <- auto$fit

  # only the parametric xbig column is planned; the smooth basis and
  # the mo() placeholder stay untouched in the scaled frame
  plan <- autoscale_plan(fa$frame)
  expect_named(plan, "y.mu")
  lp <- fa$frame$linpreds$y.mu
  expect_identical(plan$y.mu$cols, match("xbig", colnames(lp$X)))
  sf <- autoscale_frame(fa$frame, plan)
  other <- setdiff(seq_len(ncol(lp$X)), plan$y.mu$cols)
  expect_identical(sf$linpreds$y.mu$X[, other], lp$X[, other])

  # hand-standardized reference
  mb <- mean(dd$xbig); sb <- stats::sd(dd$xbig)
  dd2 <- dd
  dd2$xz <- (dd$xbig - mb) / sb
  fr <- frm(bf(y ~ xz + mo(m) + s(xs) + (1 | g)) + gaussian(),
            data = dd2)
  expect_loglik_equal(fa, fr, tol = 1e-6)
  ca <- fixef_by_dpar(fa)$mu
  cr <- fixef_by_dpar(fr)$mu
  expect_equal(ca[["xbig"]], cr[["xz"]] / sb, tolerance = 1e-6)
  expect_equal(ca[["(Intercept)"]],
               cr[["(Intercept)"]] - cr[["xz"]] * mb / sb,
               tolerance = 1e-6)
  expect_equal(unname(ca[["mom"]]), unname(cr[["mom"]]),
               tolerance = 1e-4)
  sm <- grep("s(xs)", names(ca), fixed = TRUE)
  expect_vector_equal(ca[sm], cr[sm], tol = 1e-4)
  expect_vector_equal(fa$estimates$theta, fr$estimates$theta, tol = 1e-4)
})

test_that("autoscale scales without centering when there is no intercept", {
  db <- sim_badscale_glmm(seed = 505)
  form <- bf(y ~ 0 + x1 + x2 + (1 | g)) + poisson()
  auto <- fit_catching(frm(form, data = db,
                           control = frmtmb_control(autoscale = TRUE)))
  expect_identical(auto$warnings, character(0))
  fa <- auto$fit
  plan <- autoscale_plan(fa$frame)
  expect_identical(plan$y.mu$center, c(0, 0))

  s1 <- stats::sd(db$x1); s2 <- stats::sd(db$x2)
  db2 <- db
  db2$z1 <- db$x1 / s1
  db2$z2 <- db$x2 / s2
  fr <- frm(bf(y ~ 0 + z1 + z2 + (1 | g)) + poisson(), data = db2)
  expect_loglik_equal(fa, fr, tol = 1e-6)
  cz <- fixef_by_dpar(fr)$mu
  expect_equal(unname(fixef_by_dpar(fa)$mu),
               c(cz[["z1"]] / s1, cz[["z2"]] / s2),
               tolerance = 1e-6)
})

test_that("autoscale covers dpar formulas and combines with profile", {
  set.seed(606)
  n <- 300
  g <- factor(rep(1:30, each = 10))
  xbig <- stats::rnorm(n) * 1e5
  b0 <- stats::rnorm(30, 0, 0.4)
  dd <- data.frame(xbig = xbig, g = g)
  dd$y <- stats::rnorm(n, 1 + 5e-6 * xbig + b0[g],
                       exp(-0.5 + 2e-6 * xbig))
  form <- bf(y ~ xbig + (1 | g), sigma ~ xbig) + gaussian()
  auto <- fit_catching(frm(form, data = dd,
                           control = frmtmb_control(autoscale = TRUE)))
  expect_identical(auto$warnings, character(0))
  fa <- auto$fit
  # both the beta (mu) and betad (sigma) columns are planned
  plan <- autoscale_plan(fa$frame)
  expect_setequal(vapply(plan, `[[`, "", "par"), c("beta", "betad"))

  mb <- mean(dd$xbig); sb <- stats::sd(dd$xbig)
  dd2 <- dd
  dd2$xz <- (dd$xbig - mb) / sb
  fr <- frm(bf(y ~ xz + (1 | g), sigma ~ xz) + gaussian(), data = dd2)
  expect_loglik_equal(fa, fr, tol = 1e-6)
  # within 1e-3 of each slope's own standard error: both fits step mu in
  # the response's unit (0.5 here, predictor_units()), and nlminb's
  # stopping tests leave the two 1.5e-6 and 2.1e-6 apart relative,
  # far inside an SE (dev/optrap-p1-autoscale.R)
  co <- summary(fr)$coefficients
  expect_lt(abs(fixef_by_dpar(fa)$mu[["xbig"]] * sb -
                  fixef_by_dpar(fr)$mu[["xz"]]) /
              co$mu["xz", "Std. Error"], 1e-3)
  expect_lt(abs(fixef_by_dpar(fa)$sigma[["xbig"]] * sb -
                  fixef_by_dpar(fr)$sigma[["xz"]]) /
              co$sigma["xz", "Std. Error"], 1e-3)

  # profile = TRUE moves beta into the inner problem; autoscale rides
  # on top (both stages profiled)
  db <- sim_badscale_glmm(seed = 707)
  fp <- fit_catching(frm(bf(y ~ x1 + x2 + (1 | g)) + poisson(),
                         data = db,
                         control = frmtmb_control(autoscale = TRUE,
                                                  profile = TRUE)))
  expect_identical(fp$warnings, character(0))
  m1 <- mean(db$x1); s1 <- stats::sd(db$x1)
  m2 <- mean(db$x2); s2 <- stats::sd(db$x2)
  db2 <- db
  db2$z1 <- (db$x1 - m1) / s1
  db2$z2 <- (db$x2 - m2) / s2
  frp <- frm(bf(y ~ z1 + z2 + (1 | g)) + poisson(), data = db2,
             control = frmtmb_control(profile = TRUE))
  expect_loglik_equal(fp$fit, frp, tol = 1e-6)
  expect_equal(fixef_by_dpar(fp$fit)$mu[["x1"]],
               fixef_by_dpar(frp)$mu[["z1"]] / s1,
               tolerance = 1e-5)
})

# a diagonal structure registered from outside, with a log sd per
# column; `bump` adds a constant to each variance, so that the sds no
# longer scale with the columns
register_test_diag <- function(name, bump = 0) {
  frmtmb_register_covstruct(
    name,
    lpdf = function(b, theta, blk) {
      d <- blk[["dim"]]
      B <- RTMB::matrix(b, d)
      s <- sqrt(exp(2 * theta[seq_len(d)]) + bump)
      sum(RTMB::dnorm(as.vector(B / s), 0, 1, log = TRUE)) -
        ncol(B) * sum(log(s))
    },
    vcov = function(theta, blk) {
      V <- diag(exp(2 * theta[seq_len(blk[["dim"]])]) + bump,
                blk[["dim"]])
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = function(blk) blk[["dim"]],
    sd_idx = function(blk) seq_len(blk[["dim"]]))
}

# the messages a call gives, muffled
messages_of <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(expr, message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  msgs
}

test_that("autoscale names every structure whose slope it leaves", {
  # through 0.74.0 a random slope in a structure other than us() or
  # diag() was skipped without a word, frmtmb's own homdiag() or cs() as
  # much as one another package registers. Now a slope is rescaled
  # wherever its structure's density scales with the column, and named
  # where it does not
  local_frmtmb_registries()
  register_test_diag("mydiag")
  register_test_diag("mybump", bump = 0.5)
  set.seed(31)
  d <- data.frame(g = factor(rep(1:15, each = 8)))
  d$x <- stats::rnorm(nrow(d), 0, 1e-4)
  d$y <- 1 + 2000 * d$x + stats::rnorm(15)[d$g] + stats::rnorm(nrow(d))
  ctl <- frmtmb_control(autoscale = TRUE)
  msgs <- messages_of(frm(y ~ x + mybump(1 + x | g), data = d,
                          control = ctl))
  expect_true(any(grepl("random slopes of mybump() on x", msgs,
                        fixed = TRUE)))
  # one sd over an intercept and a slope cannot scale with the slope
  msgs <- messages_of(frm(y ~ x + homdiag(1 + x | g), data = d,
                          control = ctl))
  expect_true(any(grepl("random slopes of homdiag() on x", msgs,
                        fixed = TRUE)))
  # but not under the default for the column's size alone, as through
  # 0.82.0: one sd shared with the intercept does not stall there, and
  # autoscale = FALSE reached the maximum on each such design the review
  # of lane slopescale built (F2)
  msgs <- messages_of(frm(y ~ 1 + homdiag(1 + x | g), data = d))
  expect_false(any(grepl("random slopes of", msgs)))
  # the absent cases: a structure whose density scales with the column
  # is rescaled, frmtmb's own or a registered one, and without
  # autoscale nothing is planned or named
  for (fo in list(y ~ x + (1 + x || g), y ~ x + mydiag(1 + x | g),
                  y ~ 1 + homdiag(0 + x | g))) {
    msgs <- messages_of(frm(fo, data = d, control = ctl))
    expect_false(any(grepl("random slopes of", msgs)), info = deparse(fo))
  }
  msgs <- messages_of(frm(y ~ x + mybump(1 + x | g), data = d,
                          control = frmtmb_control(autoscale = FALSE)))
  expect_false(any(grepl("random slopes of", msgs)))
})

# Item 1018 of the backlog. Seen to fail on 0.82.0 (rellib-r21,
# dev/slopescale-findings.md): the slope's column was rescaled only when
# a fixed effect read it too, so these fits stayed at their start (sd
# 1.000) at code 0, 9.2 to 12.7 log-likelihood units short on gaussian
# and up to 111 on poisson, most with no warning.
slope_data <- function(seed, s, family = "gaussian") {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- stats::rnorm(300)
  eta <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
    stats::rnorm(20, 0, 0.5)[g] * x0
  y <- if (family == "gaussian") eta + stats::rnorm(300) else
    stats::rpois(300, exp(eta))
  data.frame(g = g, x = x0 * s, y = y)
}

# the same model on the column times `s` and on the column itself agree
# where they should: the log-likelihood, every sd with the slope's in
# units of 1/s, ranef() and coef() likewise, and the fitted values
expect_same_model <- function(fs, f1, s, info = NULL) {
  l1 <- as.numeric(logLik(f1))
  expect_lt(abs(as.numeric(logLik(fs)) - l1), 1e-8 * abs(l1))
  slope <- function(nm) grepl("x$", nm)
  sds <- function(f, gn) {
    m <- VarCorr(f)[[gn]]$sd
    stats::setNames(m[, "Estimate"], rownames(m))
  }
  for (gn in names(VarCorr(f1))) {
    s1 <- sds(f1, gn)
    ss <- sds(fs, gn)
    ss[slope(names(ss))] <- ss[slope(names(ss))] * s
    expect_lt(max(abs(ss / s1 - 1)), 1e-4, label = paste(info, gn, "sd"))
  }
  for (acc in list(ranef, coef)) {
    a1 <- acc(f1)
    as_ <- acc(fs)
    for (gn in names(a1)) {
      m1 <- as.matrix(a1[[gn]])
      ms <- as.matrix(as_[[gn]])
      ms <- sweep(ms, 2, ifelse(slope(colnames(ms)), s, 1), "*")
      expect_lt(max(abs(ms - m1)) / max(abs(m1)), 1e-4,
                label = paste(info, gn, "effects"))
    }
  }
  y <- f1$frame$y[[1L]]
  expect_lt(max(abs(fitted(fs)[, "Estimate"] - fitted(f1)[, "Estimate"])) /
              stats::sd(y), 1e-4, label = paste(info, "fitted"))
}

test_that("a random slope with no fixed effect on its column is rescaled", {
  for (fam in c("gaussian", "poisson")) {
    for (fo in list(y ~ 1 + (0 + x | g), y ~ 1 + (1 + x | g))) {
      s <- 1e-4
      ds <- slope_data(2, s, fam)
      d1 <- ds
      d1$x <- ds$x / s
      fm <- if (fam == "gaussian") gaussian() else poisson()
      fs <- frm(fo, family = fm, data = ds)
      f1 <- frm(fo, family = fm, data = d1)
      expect_same_model(fs, f1, s, info = paste(fam, deparse(fo)))
      expect_true(isTRUE(fs$autoscaled))
    }
  }
  # predict() draws the same values, seed for seed
  set.seed(9)
  p1 <- predict(f1)[, "Estimate"]
  set.seed(9)
  ps <- predict(fs)[, "Estimate"]
  expect_lt(max(abs(ps - p1)) / stats::sd(d1$y), 1e-3)
})

test_that("a random slope is rescaled in every covariance structure", {
  skip_on_cran()
  # dev/slopescale-structs.R's design, seed 1
  set.seed(1)
  ng <- 20
  np <- 12
  g <- factor(rep(seq_len(ng), each = np))
  h <- factor(sample(seq_len(ng), ng * np, TRUE))
  a <- factor(rep(1:5, each = ng * np / 5))
  tf <- factor(rep(1:4, length.out = ng * np))
  x <- stats::rnorm(ng * np)
  eta <- 0.5 + stats::rnorm(ng, 0, 0.5)[g] + stats::rnorm(ng, 0, 0.6)[g] * x +
    0.5 * stats::rnorm(5, 0, 0.4)[a] * x +
    stats::rnorm(ng * 4, 0, 0.5)[(as.integer(g) - 1) * 4 + as.integer(tf)] *
      x
  sig <- exp(-0.3 + stats::rnorm(ng, 0, 0.4)[g] * x)
  dd <- data.frame(y = stats::rnorm(length(eta), eta, sig), g = g, h = h,
                   a = a, tf = tf, x = x)
  forms <- list(
    cs = bf(y ~ 1 + cs(1 + x | g)),
    ar1 = bf(y ~ 1 + (1 | g) + ar1(0 + tf:x | g)),
    homdiag = bf(y ~ 1 + (1 | g) + homdiag(0 + x | g)),
    gr = bf(y ~ 1 + (1 + x | gr(g, cor = FALSE))),
    mm = bf(y ~ 1 + (1 + x | mm(g, h))),
    dpar = bf(y ~ 1 + (1 | g), sigma ~ 1 + (0 + x | g)),
    id = bf(y ~ 1 + (1 + x | p | g), sigma ~ 1 + (1 | p | g)))
  s <- 1e-4
  ds <- dd
  ds$x <- dd$x * s
  for (nm in names(forms)) {
    # ar1's correlation is at its limit in this design, on both scales
    fs <- allow_boundary(frm(forms[[nm]], family = gaussian(), data = ds))
    f1 <- allow_boundary(frm(forms[[nm]], family = gaussian(), data = dd))
    expect_same_model(fs, f1, s, info = nm)
  }
  # a nested term: this design puts the outer level's correlation at its
  # limit, so only the likelihood is compared
  fo <- bf(y ~ 1 + (1 + x | a/g))
  fs <- allow_boundary(frm(fo, family = gaussian(), data = ds))
  f1 <- allow_boundary(frm(fo, family = gaussian(), data = dd))
  l1 <- as.numeric(logLik(f1))
  expect_lt(abs(as.numeric(logLik(fs)) - l1), 1e-8 * abs(l1))
})

test_that("each structure autoscale trusts has a density that scales", {
  # autoscale_exact_core lists the structures whose slopes are rescaled
  # without a check at fit time; this is that check, on every one, so a
  # change of parameterization fails here
  set.seed(3)
  n_g <- 12
  dd <- data.frame(g = factor(rep(seq_len(n_g), each = 5)),
                   tim = factor(rep(1:5, n_g)),
                   tn = num_factor(rep(c(0, 0.3, 0.4, 1.1, 1.5), n_g)),
                   pos = num_factor(rep(round(stats::runif(5), 2) * 10,
                                        n_g),
                                    rep(round(stats::runif(5), 2) * 10,
                                        n_g)),
                   x = stats::runif(5 * n_g), z = stats::rnorm(5 * n_g))
  dd$y <- stats::rnorm(5 * n_g)
  A <- diag(n_g) * 0.5 + 0.5
  dimnames(A) <- list(levels(dd$g), levels(dd$g))
  Q <- solve(A)
  forms <- list(
    y ~ 1 + (1 + x + z | g), y ~ 1 + (1 + x || g),
    y ~ 1 + (1 + x | gr(g, dist = "student")),
    y ~ 1 + (1 + x | gr(g, dist = "student", cor = FALSE)),
    y ~ 1 + cs(tim + 0 | g), y ~ 1 + toep(tim + 0 | g),
    y ~ 1 + hetar1(tim + 0 | g), y ~ 1 + (1 + x | gr(g, cov = A)),
    y ~ 1 + (1 + x | gr(g, prec = Q)), y ~ 1 + homdiag(1 + x | g),
    y ~ 1 + ar1(tim + 0 | g), y ~ 1 + homcs(tim + 0 | g),
    y ~ 1 + homtoep(tim + 0 | g), y ~ 1 + ou(tn + 0 | g),
    y ~ 1 + exp(pos + 0 | g), y ~ 1 + gau(pos + 0 | g),
    y ~ 1 + mat(pos + 0 | g))
  seen <- character(0)
  for (fo in forms) {
    fr <- frm(fo, family = gaussian(), data = dd, data2 = list(A = A, Q = Q),
              dry_run = "frame")
    for (bk in fr$re_blocks) {
      cs <- bk$covstruct
      D <- bk$dim
      col <- cs %in% frmtmb:::autoscale_exact_core$col
      expect_true(col || cs %in% frmtmb:::autoscale_exact_core$block,
                  info = cs)
      cpos <- if (col) 2^(seq_len(D) - 2) else rep(3, D)
      expect_true(frmtmb:::autoscale_block_exact(bk, cpos), info = cs)
      # the inverse: one scale cannot follow columns scaled apart
      if (!col && D > 1L) {
        expect_false(frmtmb:::autoscale_block_exact(bk, 2^seq_len(D)),
                     info = cs)
      }
      seen <- c(seen, cs)
    }
  }
  expect_setequal(unique(seen), unlist(frmtmb:::autoscale_exact_core))
  # a registered structure is checked at fit time, without the random
  # number generator
  local_frmtmb_registries()
  register_test_diag("mydiag")
  register_test_diag("mybump", bump = 0.5)
  for (nm in c("mydiag", "mybump")) {
    fr <- frm(stats::as.formula(paste0("y ~ 1 + ", nm, "(1 + x | g)")),
              family = gaussian(), data = dd, dry_run = "frame")
    seed0 <- .Random.seed
    ok <- frmtmb:::autoscale_block_exact(fr$re_blocks[[1L]], c(1, 0.01))
    expect_identical(.Random.seed, seed0)
    expect_identical(ok, nm == "mydiag")
  }
})

test_that("a slope column of spread 1 leaves the fit as autoscale = FALSE", {
  # the guard's condition absent: the plan holds the slope, the default
  # does not engage, and the fit is the plain one bit for bit
  d <- slope_data(2, 1)
  fd <- frm(y ~ 1 + (1 + x | g), data = d)
  ff <- frm(y ~ 1 + (1 + x | g), data = d,
            control = frmtmb_control(autoscale = FALSE))
  plan <- frmtmb:::autoscale_plan(fd$frame)
  expect_length(plan$y.mu$z, 1L)
  expect_null(frmtmb:::autoscale_decide(fd$frame, frmtmb_control(), FALSE))
  expect_identical(fd$opt$par, ff$opt$par)
  expect_identical(as.numeric(logLik(fd)), as.numeric(logLik(ff)))
  expect_false(isTRUE(fd$autoscaled))
})

test_that("an uncentered slope column of small sd engages the default", {
  # x = c0 + s z: the root mean square is near c0, the sd is s. The
  # default engages on the sd, as 0.82.0 did through the fixed twin;
  # engaging on the root mean square left these three fits 9.84, 6.05
  # and 6.09 units short (review of lane slopescale, F1,
  # dev/reviews/slopescale-review-uncentered.R)
  mk <- function(seed, s, c0) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 15))
    z <- stats::rnorm(300)
    eta <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
      stats::rnorm(20, 0, 0.5)[g] * z + 0.3 * z
    data.frame(g = g, x = c0 + s * z, y = eta + stats::rnorm(300))
  }
  for (cs in list(c(5, 2), c(5, 4), c(7, 4))) {
    d <- mk(cs[2], 0.03, cs[1])
    # the same model, the column centered and standardized
    dz <- d
    dz$x <- (d$x - mean(d$x)) / stats::sd(d$x)
    fd <- frm(y ~ x + (1 + x | g), data = d)
    fz <- frm(y ~ x + (1 + x | g), data = dz)
    lz <- as.numeric(logLik(fz))
    expect_lt(lz - as.numeric(logLik(fd)), 1e-8 * abs(lz),
              label = paste("c0", cs[1], "seed", cs[2]))
    expect_true(isTRUE(fd$autoscaled))
  }
})

test_that("a slope stopped at its start is not called at zero", {
  # autoscale = FALSE keeps the stall of item 1018: the slope's sd stays
  # at 1, where its maximum is near 5e4. 0.82.0 called it "at zero" in a
  # boundary message and gave no warning (rellib-r21); the upward probe
  # now counts the gain it finds out there, so the stop is called short
  d <- slope_data(2, 1e-5)
  conds <- character(0)
  f <- withCallingHandlers(
    frm(y ~ 1 + (0 + x | g), data = d,
        control = frmtmb_control(autoscale = FALSE)),
    message = function(m) {
      conds <<- c(conds, paste(class(m)[1], conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      conds <<- c(conds, paste(class(w)[1], conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("frmtmb_boundary_fit", conds, fixed = TRUE)))
  expect_false(any(grepl("is at zero", conds, fixed = TRUE)))
  expect_true(any(grepl("stopped short of the maximum", conds,
                        fixed = TRUE)))
})

test_that("a slope column far from 0 against its spread engages the default", {
  # x = c0 + s z: a calendar year (2000 + 3 z) and 5 + 0.1 z. Their sd
  # is above autoscale_slope_sd, so 0.83.0 did not engage, and the plain
  # fit stopped 41.1, 9.81, 15.3 and 11.1 units short with no warning
  # (dev/slopetrig-repro.R). The ratio of root mean square to sd engages
  # it. The fit is the model's maximum, and every accessor reads the
  # model as written: x = m + sx * xz maps the centered fit's slope by
  # 1 / sx and moves its intercept by the slope times m / sx
  mk <- function(seed, s, c0) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 15))
    z <- stats::rnorm(300)
    eta <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
      stats::rnorm(20, 0, 0.5)[g] * z + 0.3 * z
    data.frame(g = g, x = c0 + s * z, y = eta + stats::rnorm(300))
  }
  rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  # In the year's own units the correlation is near -1, and its
  # parameter's curvature is below what the Hessian can measure, so on
  # the first case it has no standard error, or (6 of 21 roundings,
  # dev/slopetrig-ulp.R) nlminb reports false convergence at the
  # maximum; filed in dev/test-backlog.md (round of 2026-11-16, lane
  # slopetrig)
  allow_warnings({
    for (cs in list(list(2000, 3, 2, y ~ 1 + (1 + x | g)),
                    list(2000, 3, 2, y ~ x + (1 + x | g)),
                    list(5, 0.1, 1, y ~ 1 + (1 + x | g)),
                    list(5, 0.1, 3, y ~ x + (1 + x | g)))) {
      info <- paste("c0", cs[[1]], "seed", cs[[3]], deparse(cs[[4]]))
      d <- mk(cs[[3]], cs[[2]], cs[[1]])
      m <- mean(d$x)
      sx <- stats::sd(d$x)
      dz <- d
      dz$x <- (d$x - m) / sx
      fd <- frm(cs[[4]], data = d)
      fz <- frm(cs[[4]], data = dz)
      expect_true(isTRUE(fd$autoscaled), label = info)
      lz <- as.numeric(logLik(fz))
      expect_lt(lz - as.numeric(logLik(fd)), 1e-8 * abs(lz), label = info)
      sd_x <- function(f) VarCorr(f)$g$sd["x", "Estimate"]
      expect_lt(rel(sd_x(fd) * sx, sd_x(fz)), 1e-4,
                label = paste(info, "sd"))
      rd <- as.matrix(ranef(fd)$g)
      rz <- as.matrix(ranef(fz)$g)
      expect_lt(rel(rd[, "x"] * sx, rz[, "x"]), 1e-4,
                label = paste(info, "re"))
      expect_lt(rel(rd[, "Intercept"],
                    rz[, "Intercept"] - rz[, "x"] * m / sx), 1e-4,
                label = paste(info, "re intercept"))
      cd <- as.matrix(coef(fd)$g)
      cz <- as.matrix(coef(fz)$g)
      expect_lt(rel(cd[, "x"] * sx, cz[, "x"]), 1e-4,
                label = paste(info, "coef"))
      if ("x" %in% rownames(fixef(fd))) {
        expect_lt(rel(fixef(fd)["x", c("Estimate", "Est.Error")] * sx,
                      fixef(fz)["x", c("Estimate", "Est.Error")]), 1e-4,
                  label = paste(info, "fixef"))
      }
      expect_lt(rel(fitted(fd)[, "Estimate"], fitted(fz)[, "Estimate"]),
                1e-4, label = paste(info, "fitted"))
      # Standard errors are invariant too, in exact arithmetic. The year
      # with no fixed x is left out: there the correlation's direction
      # is lost and the others are skewed (fixef 0.63 percent, fitted 2.9
      # percent; filed in dev/test-backlog.md, lane slopetrig). With no
      # fixed x the fixed intercept is the same parameter in both fits
      fx <- "x" %in% rownames(fixef(fd))
      if (fx || cs[[1]] != 2000) {
        if (!fx) {
          expect_lt(rel(fixef(fd)["Intercept", "Est.Error"],
                        fixef(fz)["Intercept", "Est.Error"]), 1e-3,
                    label = paste(info, "fixef SE"))
        }
        expect_lt(rel(fitted(fd)[, "Est.Error"], fitted(fz)[, "Est.Error"]),
                  1e-3, label = paste(info, "fitted SE"))
      }
      nd <- data.frame(g = d$g[c(1, 16)],
                       x = c(1990, 2010) * (cs[[1]] / 2000))
      nz <- nd
      nz$x <- (nd$x - m) / sx
      expect_lt(rel(fitted(fd, newdata = nd)[, "Estimate"],
                    fitted(fz, newdata = nz)[, "Estimate"]), 1e-4,
                label = paste(info, "new data"))
    }
  }, c("Standard errors are not available for 1 of 5 parameters. theta_3",
       "VarCorr() entries move along a direction",
       "Optimizer did not report convergence: false convergence (8)"))
})

test_that("a slope column near its spread's scale leaves the fit alone", {
  # the guard's condition absent: Days 0 to 9 (a ratio of root mean
  # square to sd of 1.85 over every row) and a dose given on half the
  # rows (0, or 1 + 0.2 z: 5.1 over its nonzero rows, 1.39 over every
  # row) are not near the intercept's column, so the default does not
  # engage and the fit is autoscale = FALSE's bit for bit
  set.seed(11)
  g <- factor(rep(1:18, each = 10))
  days <- rep(0:9, 18)
  d1 <- data.frame(g = g, x = days,
                   y = 250 + stats::rnorm(18, 0, 25)[g] +
                     (10 + stats::rnorm(18, 0, 6)[g]) * days +
                     stats::rnorm(180, 0, 25))
  set.seed(12)
  g <- factor(rep(1:20, each = 15))
  trt <- rep(c(TRUE, FALSE), length.out = 300)
  x <- ifelse(trt, 1 + 0.2 * stats::rnorm(300), 0)
  d2 <- data.frame(g = g, x = x,
                   y = 0.5 + stats::rnorm(20, 0, 0.5)[g] +
                     (0.3 + stats::rnorm(20, 0, 0.5)[g]) * x +
                     stats::rnorm(300))
  for (d in list(d1, d2)) {
    fd <- frm(y ~ x + (1 + x | g), data = d)
    ff <- frm(y ~ x + (1 + x | g), data = d,
              control = frmtmb_control(autoscale = FALSE))
    cl <- autoscale_plan(fd$frame)$y.mu$z[[1]]
    expect_lt(cl$ratio, autoscale_slope_ratio)
    expect_gt(cl$spread, autoscale_slope_sd)
    expect_null(autoscale_decide(fd$frame, frmtmb_control(), FALSE))
    expect_false(isTRUE(fd$autoscaled))
    expect_identical(fd$opt$par, ff$opt$par)
    expect_identical(as.numeric(logLik(fd)), as.numeric(logLik(ff)))
  }
  # and the ratio alone decides: the same Days moved to 2 to 11 engage
  d1$x <- d1$x + 2
  fr <- frm(y ~ x + (1 + x | g), data = d1)
  expect_gt(autoscale_plan(fr$frame)$y.mu$z[[1]]$ratio, autoscale_slope_ratio)
  expect_true(isTRUE(fr$autoscaled))
})

test_that("a gr(by =) block reads its ratio over the rows it covers", {
  # (1 + x | gr(g, by = h)) builds one block per level of h, each 0 on
  # the other level's rows. Over every row a calendar year there read
  # 1.41, the default did not engage, and the fit was 5.8 to 41.4 units
  # short on 10 of 10 seeds, 9 silent (review of lane slopetrig, finding
  # 3). Over the rows each block covers it reads about 690
  set.seed(1)
  g <- factor(rep(1:20, each = 15))
  z <- stats::rnorm(300)
  d <- data.frame(g = g, h = factor(ifelse(as.integer(g) <= 10, "a", "b")),
                  x = 2000 + 3 * z)
  d$y <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
    (0.3 + stats::rnorm(20, 0, 0.5)[g]) * z + stats::rnorm(300)
  fr <- frm(y ~ 1 + (1 + x | gr(g, by = h)), data = d, dry_run = "frame")
  zz <- autoscale_plan(fr)$y.mu$z
  expect_length(zz, 2L)
  want <- vapply(levels(d$h), function(l) {
    on <- d$h == l
    sqrt(mean(d$x[on]^2)) / stats::sd(d$x[on])
  }, 0)
  got <- vapply(zz, function(e) e$ratio, 0)
  expect_lt(max(abs(sort(got) / sort(want) - 1)), 1e-12)
  expect_gt(min(got), autoscale_slope_ratio)
  expect_false(is.null(autoscale_decide(fr, frmtmb_control(), FALSE)))
})


test_that("a correlation run out at a short stop is not called a boundary", {
  # autoscale = FALSE keeps the stall of a slope on a column spread 1e-4
  # or 1e-5: the slope's sd moves, but its correlation parameter runs to
  # 266 to 657 while the fit stops 3.9 to 35.1 units short. 0.83.0
  # called it "at the limit of its range" in a boundary message and gave
  # no warning (dev/slopetrig-cor.R). Moving it back toward 0 gains, so
  # it is not at its edge, and the fit warns. Which warning depends on
  # the rounding: "stopped short" where the row is flat (21 of 21 ulp
  # draws per case with the reference BLAS, dev/slopetrig-ulp.R), and
  # on the poisson case under OpenBLAS 0.3.26 a Hessian row with no
  # usable curvature; either is right, the boundary verdict is not
  for (cs in list(list(4, 1e-4, "gaussian"), list(5, 1e-4, "gaussian"),
                  list(5, 1e-5, "poisson"))) {
    info <- paste(cs, collapse = " ")
    d <- slope_data(cs[[1]], cs[[2]], cs[[3]])
    msgs <- character(0)
    warns <- character(0)
    f <- withCallingHandlers(
      frm(y ~ 1 + (1 + x | g), data = d,
          family = if (cs[[3]] == "gaussian") gaussian() else poisson(),
          control = frmtmb_control(autoscale = FALSE)),
      message = function(m) {
        msgs <<- c(msgs, paste(class(m)[1], conditionMessage(m)))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_false(any(grepl("frmtmb_boundary_fit", msgs, fixed = TRUE)),
                 label = info)
    expect_gt(length(warns), 0L, label = info)
    lost <- sdr_of(f)$se_lost
    expect_false(any(lost %in% c("boundary", "limit")), label = info)
    j <- which(names(f$opt$par) == "theta")[3]
    expect_gt(abs(f$opt$par[[j]]), se_cor_clamp_theta(), label = info)
    expect_false(se_edge_cor(f, length(f$opt$par))[j], label = info)
    expect_gt(se_cor_gain_back(f, j), se_edge_tol, label = info)
  }
})

test_that("a correlation at its limit keeps the boundary message", {
  # the guard's condition absent: no slope variance in the data, so the
  # correlation runs to +/-1 and its parameter past the clamp; every
  # move back toward 0 loses, and the fit is the boundary fit it was
  # (test-review-v29.R, "a lost correlation leaves the sd intervals")
  for (s in c(7L, 8L)) {
    set.seed(s)
    g <- factor(rep(seq_len(20), each = 10))
    d <- data.frame(g = g, x = stats::rnorm(200))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.7)[g] + stats::rnorm(200)
    conds <- character(0)
    f <- withCallingHandlers(
      frm(y ~ x + (1 + x | g), data = d),
      message = function(m) {
        conds <<- c(conds, paste(class(m)[1], conditionMessage(m)))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        conds <<- c(conds, paste(class(w)[1], conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
    j <- which(names(f$opt$par) == "theta")[3]
    expect_gt(abs(f$opt$par[[j]]), se_cor_clamp_theta())
    expect_identical(se_cor_gain_back(f, j), 0)
    expect_identical(sdr_of(f)$se_lost, c(theta_3 = "boundary"))
    expect_true(any(grepl("frmtmb_boundary_fit", conds, fixed = TRUE)))
    expect_false(any(grepl("stopped short", conds, fixed = TRUE)))
  }
})
