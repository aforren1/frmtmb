# me() on new data against brms 2.23.0 at the same parameter values (a
# Fixed_param brmsfit, 4000 equal draws, so the spread over draws is the
# latent value's and the observation's alone). Seen to fail on 0.78.0
# (rellib-r16): predict() had the spread of sigma alone, sd ratio 0.55 at
# sdx = 0.8 (dev/predread-me-compare.R).

skip_unless_brms()
skip_on_cran()
if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
  skip("set FRMTMB_BRMS_FIT_TESTS=true to run the fits against brms")
}

mnb_data <- function() {
  set.seed(41)
  n <- 120
  g <- rep(1:12, each = 10)
  d <- data.frame(g = factor(g), z = stats::rnorm(n), x = stats::rnorm(n),
                  sx = 0.4, sw = 0.25, sg = 0.3)
  lat <- stats::rnorm(12)[g]
  d$xg <- lat + stats::rnorm(12, 0, 0.3)[g]
  d$w <- 0.6 * d$x + stats::rnorm(n)
  d$y <- 1 + 0.8 * d$x + 0.5 * d$w + 0.6 * lat + 0.3 * d$z +
    stats::rnorm(n, 0, 0.5)
  d
}
mnb_grid <- function() {
  data.frame(x = c(-1, 0, 1), sx = c(0.2, 0.4, 0.8), w = 0.5, sw = 0.25,
             z = c(0, 1, -1), xg = c(0.5, 0.5, -1), sg = 0.3,
             g = factor(c("1", "1", "2"), levels = as.character(1:12)))
}
# per-row z-scores of mean and variance against their Monte Carlo errors,
# and the two-sided Bonferroni point at 0.001, the survey's rule
mnb_agree <- function(a, b) {
  va <- apply(a, 2, stats::var)
  vb <- apply(b, 2, stats::var)
  zm <- (colMeans(a) - colMeans(b)) / sqrt(va / nrow(a) + vb / nrow(b))
  m4 <- function(s) colMeans(sweep(s, 2, colMeans(s))^4)
  zv <- (va - vb) / sqrt((m4(a) - va^2) / nrow(a) + (m4(b) - vb^2) / nrow(b))
  crit <- stats::qnorm(1 - 0.001 / (4 * ncol(a)))
  expect_lt(max(abs(zm)), crit)
  expect_lt(max(abs(zv)), crit)
}

test_that("predict(newdata =) has brms's spread on every shape", {
  forms <- list(
    list(f = bf(y ~ me(x, sx) + z), b = brms::bf(y ~ me(x, sx) + z)),
    list(f = bf(y ~ me(xg, sg, gr = g) + z),
         b = brms::bf(y ~ me(xg, sg, gr = g) + z)),
    list(f = bf(y ~ me(x, sx) + me(w, sw)),
         b = brms::bf(y ~ me(x, sx) + me(w, sw))),
    list(f = bf(y ~ me(x, sx) + me(x, sx):z),
         b = brms::bf(y ~ me(x, sx) + me(x, sx):z)))
  d <- mnb_data()
  nd <- mnb_grid()
  with_brms_me(function() {
    for (m in forms) {
      fit <- frm(m$f, data = d, family = gaussian())
      bfit <- brms_fixed_fit(m$b, stats::gaussian(), d, fit, ndraws = 4000)
      set.seed(7)
      fp <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000L,
                    propagate_error = FALSE)
      bp <- brms::posterior_predict(bfit, newdata = nd)
      mnb_agree(fp, bp)
      s <- simulate(fit, nsim = 4000L, seed = 7, newdata = nd)
      mnb_agree(t(vapply(s, function(v) as.numeric(v), numeric(3L))), bp)
      # a summary at the estimate: the mean of brms's draw, whose spread is
      # the part of Est.Error that is not coefficient uncertainty
      be <- brms::posterior_epred(bfit, newdata = nd)
      lb <- frm_lp_basis(fit, newdata = nd)
      mcse <- apply(be, 2, stats::sd) / sqrt(nrow(be))
      expect_lt(max(abs(as.numeric(lb$eta) - colMeans(be)) / mcse), 4.5)
      ratio <- sqrt(lb$extra_var) / apply(be, 2, stats::sd)
      expect_lt(max(abs(ratio - 1)), 4 / sqrt(2 * nrow(be)))
    }
  })
})

test_that("rows of one me(gr = ) level share a latent value, as in brms", {
  d <- mnb_data()
  nd <- mnb_grid()
  fit <- frm(bf(y ~ me(xg, sg, gr = g) + z), data = d, family = gaussian())
  with_brms_me(function() {
    bfit <- brms_fixed_fit(brms::bf(y ~ me(xg, sg, gr = g) + z),
                           stats::gaussian(), d, fit, ndraws = 2000)
    be <- brms::posterior_epred(bfit, newdata = nd)
    expect_equal(stats::cor(be[, 1], be[, 2]), 1)
    expect_lt(abs(stats::cor(be[, 1], be[, 3])), 4 / sqrt(2000))
    E <- as.matrix(frm_lp_basis(fit, newdata = nd,
                                extra_cov = TRUE)$extra_cov)
    expect_equal(E[1, 2], E[1, 1])
    expect_equal(E[1, 3], 0)
    # both refuse an x that varies within a level, and a level the fit did
    # not see unless it is allowed
    bad <- nd
    bad$xg[2L] <- 9
    expect_error(brms::posterior_epred(bfit, newdata = bad),
                 "unique for each group")
    expect_error(fitted(fit, newdata = bad), "unique for each group")
    nw <- nd
    nw$g <- as.character(nw$g)
    nw$g[3L] <- "new"
    expect_error(brms::posterior_epred(bfit, newdata = nw), "allow_new_levels")
    expect_error(fitted(fit, newdata = nw), "allow_new_levels")
  })
})

test_that("two responses read one latent; brms draws one per predictor", {
  d <- mnb_data()
  d$y2 <- -0.5 + 0.7 * d$x + stats::rnorm(nrow(d), 0, 0.6)
  nd <- mnb_grid()
  fit <- frm(mvbf(bf(y ~ me(x, sx)), bf(y2 ~ me(x, sx)), rescor = FALSE),
             data = d, family = gaussian())
  with_brms_me(function() {
    bfit <- brms_fixed_fit(brms::mvbf(brms::bf(y ~ me(x, sx)),
                                      brms::bf(y2 ~ me(x, sx)),
                                      rescor = FALSE),
                           stats::gaussian(), d, fit, ndraws = 2000)
    be <- brms::posterior_epred(bfit, newdata = nd)
    bound <- 4 / sqrt(2000)
    for (i in 1:3) expect_lt(abs(stats::cor(be[, i, 1], be[, i, 2])), bound)
    set.seed(9)
    fp <- t(replicate(2000, {
      sh <- fit
      sh[["me_draw"]] <- new.env(parent = emptyenv())
      c(frm_linpred(sh, newdata = nd, resp = "y"),
        frm_linpred(sh, newdata = nd, resp = "y2"))
    }))
    # one latent, a stated divergence: brms's Stan program has one Xme_1, and its
    # replicates contradict it (dev/review-predread-mejoint.R)
    for (i in 1:3) expect_gt(stats::cor(fp[, i], fp[, 3L + i]), 0.999)
  })
})

test_that("fitted() on a new row is brms's mean, and its interval brms's", {
  # dev/review-predread-mefam.R, seed 7971: brms over frm read 1.41 to 1.76 at
  # sdx 1.2 on log links, and fitted()'s interval left the support
  set.seed(7971)
  n <- 300
  d <- data.frame(lat = stats::rnorm(n), sx = 0.4)
  d$xo <- d$lat + stats::rnorm(n, 0, 0.4)
  eta <- 0.4 + 0.7 * d$lat
  d$cnt <- stats::rpois(n, exp(eta))
  d$pos <- exp(eta + stats::rnorm(n, 0, 0.5))
  d$hp <- ifelse(d$cnt == 0, 1L, d$cnt)
  d$hp[stats::runif(n) < 0.3] <- 0L
  nd <- data.frame(xo = 1, sx = c(0.3, 0.8, 1.2))
  specs <- list(
    list(f = bf(cnt ~ me(xo, sx)), b = brms::bf(cnt ~ me(xo, sx)),
         ff = poisson(), bb = stats::poisson()),
    list(f = bf(pos ~ me(xo, sx)), b = brms::bf(pos ~ me(xo, sx)),
         ff = lognormal(), bb = brms::lognormal()),
    list(f = bf(hp ~ me(xo, sx)), b = brms::bf(hp ~ me(xo, sx)),
         ff = hurdle_poisson(), bb = brms::hurdle_poisson()))
  with_brms_me(function() {
    for (m in specs) {
      fit <- frm(m$f, family = m$ff, data = d)
      bfit <- brms_fixed_fit(m$b, m$bb, d, fit, ndraws = 4000)
      be <- brms::posterior_epred(bfit, newdata = nd)
      fi <- fitted(fit, newdata = nd)
      # heavy tails: the mean within 5 Monte Carlo errors
      mcse <- apply(be, 2, stats::sd) / sqrt(nrow(be))
      expect_lt(max(abs(fi[, "Estimate"] - colMeans(be)) / mcse), 5)
      # brms's draws hold no coefficient uncertainty, so this is a floor
      expect_true(all(fi[, "Est.Error"] > 0.8 * apply(be, 2, stats::sd)))
      expect_true(all(fi[, "Q2.5"] > 0))
      lo <- apply(be, 2, stats::quantile, 0.025)
      hi <- apply(be, 2, stats::quantile, 0.975)
      expect_true(all(abs(log(fi[, "Q2.5"] / lo)) < 0.25))
      expect_true(all(abs(log(fi[, "Q97.5"] / hi)) < 0.25))
    }
  })
})
