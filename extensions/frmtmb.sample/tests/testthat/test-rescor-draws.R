# The predictive draws and the pointwise density of a set_rescor(TRUE)
# model on draws, against brms 2.23.0's posterior_predict_gaussian_mv(),
# posterior_predict_student_mv() and log_lik(resp =).
#
# SEEN TO FAIL on frmtmb.sample 0.24.0 (rellib-r15): posterior_predict()
# without `resp` returned the FIRST response alone, so no joint draw
# existed, and resp = c(...) was an error; log_lik(resp = "y1") of a
# rescor model returned the joint density of every response, and
# log_lik(combine = FALSE) was refused (dev/rescorpp-findings.md has
# the failing output). The z bounds are on a Monte Carlo error the run
# measures from the spread over rows.
#
# Fix round 1 of the lane's review added blocks for defects of the
# lane's own first build, each seen to fail on that build by a mutation
# of the loaded namespace (dev/rescorpp-mutate.R): the Student-t tails
# (a normal draw), the three-response permutations (a correlation left
# in model order), predict()'s category refusal and the repeated name.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

rc_case <- local({
  cache <- list()
  function(which) {
    skip_sampler()
    if (is.null(cache[[which]])) {
      set.seed(switch(which, g3 = 11, gs = 12, st = 13, g2f = 14, ga = 15,
                      ord = 16))
      n <- 120
      x <- stats::rnorm(n)
      Z <- matrix(stats::rnorm(3 * n), n) %*%
        chol(matrix(c(1, .8, .5, .8, 1, .3, .5, .3, 1), 3))
      dd <- data.frame(x = x)
      if (which == "st") {
        Z <- Z / sqrt(stats::rchisq(n, 6) / 6)
      }
      s1 <- if (which == "gs") exp(0.5 * x) else 1
      dd$y1 <- 1 + x + s1 * Z[, 1]
      dd$y2 <- -1 + 0.5 * x + 1.5 * Z[, 2]
      dd$y3 <- 0.3 * x + 0.5 * Z[, 3]
      # names with an underscore, which brms spells without it
      dd$y_a <- dd$y1
      dd$y_b <- dd$y2
      dd$o <- cut(x + stats::rlogis(n), c(-Inf, -0.5, 0.5, Inf),
                  labels = FALSE)
      form <- switch(which,
        g3 = mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x)) + set_rescor(TRUE),
        gs = mvbf(bf(y1 ~ x, sigma ~ x), bf(y2 ~ x)) + set_rescor(TRUE),
        st = mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE),
        g2f = mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(FALSE),
        ga = mvbf(bf(y_a ~ x), bf(y_b ~ x)) + set_rescor(TRUE),
        ord = mvbf(bf(o ~ x, family = cumulative()),
                   bf(y2 ~ x, family = gaussian()), rescor = FALSE))
      fit <- if (which == "ord") {
        frm(form, data = dd)
      } else {
        frm(form, family = if (which == "st") student() else gaussian(),
            data = dd)
      }
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 1200, refresh = 0, seed = 4)))
      cache[[which]] <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache[[which]]
  }
})

# per row, the correlation over draws of two responses' predictive
# residuals, and the mean of `a - b` over rows in units of its own
# standard error
rc_rowcor <- function(e1, e2, method = "pearson") {
  vapply(seq_len(ncol(e1)), function(i) {
    stats::cor(e1[, i], e2[, i], method = method)
  }, 1)
}
rc_z <- function(a, b) {
  dlt <- a - b
  mean(dlt) / (stats::sd(dlt) / sqrt(length(dlt)))
}
# the expected correlation over draws given the draws themselves:
# Cov(e_j, e_k | draw) = r sigma_j sigma_k
rc_target <- function(ds, rj, rk) {
  dr <- as.data.frame(posterior::as_draws_df(ds))
  r <- dr[[paste0("rescor__", rj, "__", rk)]]
  sj <- posterior_epred(ds, dpar = "sigma", resp = rj)
  sk <- posterior_epred(ds, dpar = "sigma", resp = rk)
  colMeans(r * sj * sk) / sqrt(colMeans(sj^2) * colMeans(sk^2))
}

test_that("posterior_predict() of a rescor model is brms's joint draw", {
  cs <- rc_case("g3")
  ds <- cs$ds
  set.seed(1)
  pp <- posterior_predict(ds)
  # brms's layout: draws x observations x responses, named as brms
  # names the responses
  expect_identical(dim(pp), c(ndraws(ds), nrow(cs$dd), 3L))
  expect_identical(dimnames(pp)[[3L]], c("y1", "y2", "y3"))
  ep <- posterior_epred(ds)
  expect_identical(dim(ep), dim(pp))
  e <- pp - ep
  for (pr in list(c(1L, 2L), c(1L, 3L), c(2L, 3L))) {
    rj <- c("y1", "y2", "y3")[pr[1L]]
    rk <- c("y1", "y2", "y3")[pr[2L]]
    cf <- rc_rowcor(e[, , pr[1L]], e[, , pr[2L]])
    tg <- rc_target(ds, rj, rk)
    expect_lt(abs(rc_z(cf, tg)), 4.5, label = paste(rj, rk))
  }
  # the independent draw this replaced is far outside that bound: its
  # correlation is 0 where the target is near 0.8
  cf0 <- rc_rowcor(posterior_predict(ds, resp = "y1") - ep[, , 1L],
                   posterior_predict(ds, resp = "y2") - ep[, , 2L])
  expect_gt(abs(rc_z(cf0, rc_target(ds, "y1", "y2"))), 20)
  # the spread of a sum carries the correlation: Var(e1 + e2) against
  # its expectation given the draws
  dr <- as.data.frame(posterior::as_draws_df(ds))
  ev <- mean(dr$sigma_y1^2 + dr$sigma_y2^2 +
               2 * dr$rescor__y1__y2 * dr$sigma_y1 * dr$sigma_y2)
  vs <- apply(e[, , 1L] + e[, , 2L], 2L, stats::var)
  expect_lt(abs(rc_z(vs, rep(ev, length(vs)))), 4.5)
})

test_that("resp = narrows the joint draw, and one response is its margin", {
  cs <- rc_case("g3")
  ds <- cs$ds
  set.seed(2)
  p13 <- posterior_predict(ds, resp = c("y1", "y3"))
  expect_identical(dim(p13), c(ndraws(ds), nrow(cs$dd), 2L))
  expect_identical(dimnames(p13)[[3L]], c("y1", "y3"))
  ep <- posterior_epred(ds, resp = c("y1", "y3"))
  expect_identical(dimnames(ep)[[3L]], c("y1", "y3"))
  cf <- rc_rowcor(p13[, , 1L] - ep[, , 1L], p13[, , 2L] - ep[, , 2L])
  expect_lt(abs(rc_z(cf, rc_target(ds, "y1", "y3"))), 4.5)
  # the order named is the order of the layers, and of the correlation.
  # Three responses, because a 2 x 2 correlation matrix is the same in
  # either order: a matrix left in model order would give the layers
  # (y3, y1) the correlation of (y1, y2), 0.8 where 0.5 is wanted
  perm <- c("y3", "y1", "y2")
  pp3 <- posterior_predict(ds, resp = perm)
  expect_identical(dimnames(pp3)[[3L]], perm)
  ep3 <- posterior_epred(ds, resp = perm)
  for (pr in list(c(1L, 2L), c(1L, 3L), c(2L, 3L))) {
    cf <- rc_rowcor(pp3[, , pr[1L]] - ep3[, , pr[1L]],
                    pp3[, , pr[2L]] - ep3[, , pr[2L]])
    tg <- rc_target(ds, sort(perm[pr])[1L], sort(perm[pr])[2L])
    expect_lt(abs(rc_z(cf, tg)), 4.5, label = paste(perm[pr], collapse = ""))
  }
  # one response: a matrix, its marginal normal with sd sigma
  p1 <- posterior_predict(ds, resp = "y1")
  expect_true(is.matrix(p1))
  dr <- as.data.frame(posterior::as_draws_df(ds))
  sdr <- apply(p1 - ep[, , 1L], 2L, stats::sd)
  expect_lt(abs(rc_z(sdr, rep(sqrt(mean(dr$sigma_y1^2)), length(sdr)))),
            4.5)
})

test_that("a draw reads the draws draw_ids names, through core's draw", {
  cs <- rc_case("g3")
  ds <- cs$ds
  set.seed(9)
  pp <- posterior_predict(ds, draw_ids = c(5L, 9L))
  set.seed(9)
  ok <- rep(TRUE, nrow(cs$dd))
  for (k in 1:2) {
    sh <- frmtmb.sample:::draws_fit_at(ds, c(5L, 9L)[k])
    ys <- frmtmb::predict_rescor_draw(sh, sh$spec$responses,
                                      frmtmb::eval_dpars(sh), ok)
    for (j in 1:3) expect_identical(pp[k, , j], ys[[j]])
  }
})

test_that("a joint draw at newdata, with sigma varying by row", {
  cs <- rc_case("gs")
  ds <- cs$ds
  set.seed(3)
  pp <- posterior_predict(ds)
  ep <- posterior_epred(ds)
  e <- pp - ep
  cf <- rc_rowcor(e[, , 1L], e[, , 2L])
  expect_lt(abs(rc_z(cf, rc_target(ds, "y1", "y2"))), 4.5)
  # each row's spread follows its own sigma
  s1 <- posterior_epred(ds, dpar = "sigma", resp = "y1")
  expect_gt(stats::sd(colMeans(s1)) / mean(colMeans(s1)), 0.2)
  sdr <- apply(e[, , 1L], 2L, stats::sd)
  expect_lt(abs(rc_z(sdr / sqrt(colMeans(s1^2)), rep(1, length(sdr)))),
            4.5)
  nd <- data.frame(x = seq(-2, 2, length.out = 40))
  pn <- posterior_predict(ds, newdata = nd)
  expect_identical(dim(pn), c(ndraws(ds), 40L, 2L))
  en <- posterior_epred(ds, newdata = nd)
  sn <- posterior_epred(ds, newdata = nd, dpar = "sigma", resp = "y1")
  s2 <- posterior_epred(ds, newdata = nd, dpar = "sigma", resp = "y2")
  r <- as.data.frame(posterior::as_draws_df(ds))$rescor__y1__y2
  tg <- colMeans(r * sn * s2) / sqrt(colMeans(sn^2) * colMeans(s2^2))
  cn <- rc_rowcor(pn[, , 1L] - en[, , 1L], pn[, , 2L] - en[, , 2L])
  expect_lt(abs(rc_z(cn, tg)), 4.5)
})

test_that("a Student-t rescor draw is the multivariate t, one nu per row", {
  skip_if_not_installed("mvtnorm")
  cs <- rc_case("st")
  ds <- cs$ds
  set.seed(4)
  pp <- posterior_predict(ds)
  ep <- posterior_epred(ds)
  e <- pp - ep
  # the reference: mvtnorm's multivariate t at each draw's parameters,
  # whose one chi-square per row the responses share
  dr <- as.data.frame(posterior::as_draws_df(ds))
  n <- nrow(cs$dd)
  ref <- array(NA_real_, dim(e))
  for (s in seq_len(ndraws(ds))) {
    r <- dr$rescor__y1__y2[s]
    z <- mvtnorm::rmvt(n, sigma = matrix(c(1, r, r, 1), 2), df = dr$nu[s])
    ref[s, , 1L] <- dr$sigma_y1[s] * z[, 1L]
    ref[s, , 2L] <- dr$sigma_y2[s] * z[, 2L]
  }
  for (m in c("spearman", "kendall")) {
    expect_lt(abs(rc_z(rc_rowcor(e[, , 1L], e[, , 2L], m),
                       rc_rowcor(ref[, , 1L], ref[, , 2L], m))), 4.5,
              label = m)
  }
  # a shared chi-square couples the SIZES of the residuals beyond what
  # their correlation alone gives
  expect_lt(abs(rc_z(rc_rowcor(abs(e[, , 1L]), abs(e[, , 2L]), "spearman"),
                     rc_rowcor(abs(ref[, , 1L]), abs(ref[, , 2L]),
                               "spearman"))), 4.5)
  # The rank checks above barely tell a multivariate t from a normal, so
  # the tails decide: each margin standardized by its draw's sigma is
  # t(nu) at that draw's nu, and lands beyond the draw's 99.5 percent
  # t quantile in 1 percent of the cells. A normal draw puts about 0.03
  # percent there at nu near 6; each cell is an independent draw, so
  # the count is binomial
  for (j in 1:2) {
    u <- e[, , j] / dr[[paste0("sigma_y", j)]]
    hit <- abs(u) > stats::qt(0.995, dr$nu)
    z <- (mean(hit) - 0.01) / sqrt(0.01 * 0.99 / length(hit))
    expect_lt(abs(z), 4.5, label = paste("tail share of y", j))
  }
})

test_that("log_lik(resp =) of a rescor model is brms's: a margin or a joint", {
  skip_if_not_installed("mvtnorm")
  cs <- rc_case("g3")
  ds <- cs$ds
  dd <- cs$dd
  dr <- as.data.frame(posterior::as_draws_df(ds))
  S <- ndraws(ds)
  n <- nrow(dd)
  mu <- function(r) {
    outer(dr[[paste0("b_", r, "_Intercept")]], rep(1, n)) +
      outer(dr[[paste0("b_", r, "_x")]], dd$x)
  }
  # one response: its own normal density, as brms's log_lik(resp = "y1")
  l1 <- log_lik(ds, resp = "y1")
  ref1 <- stats::dnorm(matrix(dd$y1, S, n, byrow = TRUE), mu("y1"),
                       dr$sigma_y1, log = TRUE)
  expect_lt(max(abs(l1 - ref1)), 64 * .Machine$double.eps * max(abs(ref1)))
  # two of three: their bivariate density, in either order
  rows <- c(1L, 7L, S)
  m1 <- mu("y1")
  m3 <- mu("y3")
  l13 <- log_lik(ds, resp = c("y1", "y3"))
  expect_identical(dim(l13), c(S, n))
  ref13 <- t(vapply(rows, function(s) {
    r <- dr$rescor__y1__y3[s]
    D <- diag(c(dr$sigma_y1[s], dr$sigma_y3[s]))
    Sg <- D %*% matrix(c(1, r, r, 1), 2) %*% D
    vapply(seq_len(n), function(i) {
      mvtnorm::dmvnorm(c(dd$y1[i], dd$y3[i]), c(m1[s, i], m3[s, i]), Sg,
                       log = TRUE)
    }, 1)
  }, numeric(n)))
  expect_lt(max(abs(l13[rows, ] - ref13)),
            1e3 * .Machine$double.eps * max(abs(ref13)))
  # all three in another order, against the density written in that
  # order: a 2 x 2 correlation matrix is the same either way round, so
  # only a permutation of three can catch a matrix left in model order
  perm <- c("y3", "y1", "y2")
  mm <- list(y1 = m1, y2 = mu("y2"), y3 = m3)
  lp <- log_lik(ds, resp = perm)
  refp <- t(vapply(rows, function(s) {
    C <- diag(3)
    dimnames(C) <- list(c("y1", "y2", "y3"), c("y1", "y2", "y3"))
    C["y1", "y2"] <- C["y2", "y1"] <- dr$rescor__y1__y2[s]
    C["y1", "y3"] <- C["y3", "y1"] <- dr$rescor__y1__y3[s]
    C["y2", "y3"] <- C["y3", "y2"] <- dr$rescor__y2__y3[s]
    D <- diag(vapply(perm, function(r) dr[[paste0("sigma_", r)]][s], 1))
    Sg <- D %*% C[perm, perm] %*% D
    vapply(seq_len(n), function(i) {
      mvtnorm::dmvnorm(vapply(perm, function(r) dd[[r]][i], 1),
                       vapply(perm, function(r) mm[[r]][s, i], 1), Sg,
                       log = TRUE)
    }, 1)
  }, numeric(n)))
  expect_lt(max(abs(lp[rows, ] - refp)),
            1e3 * .Machine$double.eps * max(abs(refp)))
  # brms returns the joint matrix whatever combine says
  expect_identical(log_lik(ds, combine = FALSE), log_lik(ds))
})

test_that("a Student-t response alone is the t margin with the shared nu", {
  # brms 2.23.0 reads nu = 1 here (its one-response prep finds no nu of
  # its own), which is not the margin of its multivariate t; frmtmb
  # uses the shared nu. dev/rescorpp-findings.md has brms's numbers
  cs <- rc_case("st")
  ds <- cs$ds
  dd <- cs$dd
  dr <- as.data.frame(posterior::as_draws_df(ds))
  n <- nrow(dd)
  m1 <- outer(dr$b_y1_Intercept, rep(1, n)) + outer(dr$b_y1_x, dd$x)
  z <- (matrix(dd$y1, ndraws(ds), n, byrow = TRUE) - m1) / dr$sigma_y1
  ref <- stats::dt(z, dr$nu, log = TRUE) - log(dr$sigma_y1)
  l1 <- log_lik(ds, resp = "y1")
  expect_lt(max(abs(l1 - ref)), 1e3 * .Machine$double.eps * max(abs(ref)))
})

test_that("without rescor each layer is the response's own draw", {
  # the guard-absent case: an independent multivariate model keeps its
  # per-response simulator, so the first layer is bit for bit the draw
  # resp = "y1" gives from the same stream
  cs <- rc_case("g2f")
  ds <- cs$ds
  set.seed(6)
  pp <- posterior_predict(ds)
  set.seed(6)
  p1 <- posterior_predict(ds, resp = "y1")
  expect_identical(dim(pp), c(ndraws(ds), nrow(cs$dd), 2L))
  expect_identical(unname(pp[, , 1L]), p1)
  e <- pp - posterior_epred(ds)
  expect_lt(abs(rc_z(rc_rowcor(e[, , 1L], e[, , 2L]),
                     rep(0, nrow(cs$dd)))), 4.5)
  expect_identical(unname(posterior_epred(ds)[, , 2L]),
                   posterior_epred(ds, resp = "y2"))
})

test_that("the summaries and the error follow brms's multivariate layout", {
  cs <- rc_case("g3")
  ds <- cs$ds
  n <- nrow(cs$dd)
  set.seed(7)
  pe <- predictive_error(ds, ndraws = 50)
  set.seed(7)
  pp <- posterior_predict(ds, ndraws = 50)
  for (j in 1:3) {
    y <- cs$dd[[c("y1", "y2", "y3")[j]]]
    expect_identical(pe[, , j], rep(y, each = 50L) - pp[, , j])
  }
  expect_identical(dim(fitted(ds, ndraws = 20)), c(n, 4L, 3L))
  expect_identical(dim(predict(ds, ndraws = 20)), c(n, 4L, 3L))
  expect_identical(dim(residuals(ds, ndraws = 20)), c(n, 4L, 3L))
  expect_identical(dim(residuals(ds, ndraws = 20, type = "pearson")),
                   c(n, 4L, 3L))
  expect_identical(dimnames(fitted(ds, ndraws = 20))[[3L]],
                   c("y1", "y2", "y3"))
  expect_identical(dim(posterior_linpred(ds, ndraws = 20)), c(20L, n, 3L))
  # brms stops here too: its posterior_predict() is not a matrix
  expect_error(predictive_interval(ds), "needs one response")
  expect_identical(dim(predictive_interval(ds, resp = "y2")), c(n, 2L))
  expect_error(posterior_predict(ds, resp = "nope"),
               "Valid response variables")
})

test_that("resp is read the same way by log_lik() and the draw methods", {
  # brms's spelling of y_a is ya; a name given twice is refused, as brms
  # stops on it; an empty resp is every response, as brms's
  # validate_resp() reads it. On frmtmb.sample 0.24.0 log_lik() refused
  # "ya", returned Inf in every cell for a repeated name (a singular
  # correlation matrix) and a matrix with no responses for character(0)
  cs <- rc_case("ga")
  ds <- cs$ds
  l_a <- log_lik(ds, resp = "y_a")
  expect_identical(log_lik(ds, resp = "ya"), l_a)
  expect_identical(log_lik(ds, resp = c("yb", "y_a")),
                   log_lik(ds, resp = c("y_b", "y_a")))
  expect_identical(log_lik(ds, resp = character(0)), log_lik(ds))
  expect_error(log_lik(ds, resp = c("y_a", "y_a")), "twice")
  expect_error(log_lik(ds, resp = c("ya", "y_a")), "twice")
  expect_error(posterior_predict(ds, resp = c("y_b", "yb")), "twice")
  expect_identical(dim(posterior_predict(ds, resp = character(0),
                                         ndraws = 3)),
                   c(3L, nrow(cs$dd), 2L))
  # loo() reads its resp through log_lik()
  la <- loo(ds, resp = "ya")
  lb <- loo(ds, resp = "y_a")
  expect_identical(la$estimates, lb$estimates)
  expect_identical(attr(la, "yhash"), attr(lb, "yhash"))
})

test_that("predict() refuses to summarize a category beside another response", {
  # Without resp, predict() summarized an ordinal response's category
  # CODES as numbers (a mean of 1.9 for categories 1 to 3), so
  # predict(ds)[, , "o"] differed from predict(ds, resp = "o"), the
  # category shares; brms stops on both calls
  cs <- rc_case("ord")
  ds <- cs$ds
  expect_error(predict(ds, ndraws = 20), "Response 'o' is a category")
  expect_error(predict(ds, resp = c("y2", "o"), ndraws = 20),
               "Response 'o' is a category")
  # one response, and the draws themselves, still answer
  expect_identical(dim(predict(ds, resp = "o", ndraws = 20)),
                   c(nrow(cs$dd), 3L))
  expect_identical(dim(predict(ds, resp = "y2", ndraws = 20)),
                   c(nrow(cs$dd), 4L))
  expect_identical(dim(predict(ds, ndraws = 20, summary = FALSE)),
                   c(20L, nrow(cs$dd), 2L))
})

test_that("the loo_*() functions take several responses, as brms does", {
  skip_if_not_installed("loo")
  cs <- rc_case("g3")
  ds <- cs$ds
  quiet <- function(e) suppressMessages(e)
  le <- quiet(loo_epred(ds, resp = c("y3", "y1")))
  expect_identical(dim(le), c(nrow(cs$dd), 1L, 2L))
  expect_identical(dimnames(le)[[3L]], c("y3", "y1"))
  # the weights are those of the set's joint density
  ps <- frmtmb.sample:::draws_psis_of(log_lik(ds, resp = c("y3", "y1")), ds)
  expect_identical(le, quiet(loo_epred(ds, resp = c("y3", "y1"),
                                       psis_object = ps)))
  set.seed(1)
  lpi <- quiet(loo_predictive_interval(ds, resp = c("y1", "y2")))
  expect_identical(dim(lpi), c(nrow(cs$dd), 2L, 2L))
})
