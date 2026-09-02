# Closing the ordinal prediction surface (v0.32): fitted() returns the
# category-probability matrix predict(type = "response") returns, and
# every consumer that needed one number per row now says which number it
# uses. Before this, fitted() stayed on the latent predictor, which
# broke the predict(type = "response") == fitted() invariant, made
# plot(fit) and residuals() report latent-scale quantities without
# saying so, and left dharma_residuals() and conditional_effects()
# refusing or plotting the wrong scale.

skip_on_cran()

ordfit_data <- function(seed, n = 250, tau = c(-0.8, 0.6), beta = 0.9,
                        levels = c("lo", "mid", "hi")) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  eta <- beta * dd$x
  p <- cbind(stats::plogis(tau[1] - eta),
             stats::plogis(tau[2] - eta) - stats::plogis(tau[1] - eta),
             1 - stats::plogis(tau[2] - eta))
  dd$y <- factor(levels[apply(p, 1L, function(pr) sample(3L, 1L,
                                                         prob = pr))],
                 levels = levels, ordered = TRUE)
  dd
}

## fitted() == predict(type = "response") ------------------------------

test_that("fitted() is the category matrix on all four ordinal families", {
  dd <- ordfit_data(101)
  for (fam in list(cumulative(), sratio(), cratio(), acat())) {
    fit <- frm(bf(y ~ x) + fam, data = dd)
    ft <- fitted(fit)
    expect_true(is.matrix(ft), info = fam$family)
    expect_equal(dim(ft), c(nrow(dd), 3L), info = fam$family)
    expect_equal(colnames(ft), levels(dd$y), info = fam$family)
    expect_equal(rownames(ft), rownames(dd), info = fam$family)
    expect_equal(unname(rowSums(ft)), rep(1, nrow(dd)),
                 tolerance = 1e-12, info = fam$family)
    expect_true(all(ft > 0 & ft < 1), info = fam$family)
    # the invariant the v0.31 note carved an exception out of
    expect_identical(ft, predict(fit, type = "response"),
                     info = fam$family)
    # the latent predictor is still reachable, by name
    expect_true(is.numeric(predict(fit, type = "link")),
                info = fam$family)
  }
})

test_that("fitted() honors cs() and pads under na.exclude", {
  dd <- ordfit_data(102)
  fit <- frm(bf(y ~ x + cs(z)) + sratio(), data = dd)
  ft <- fitted(fit)
  expect_identical(ft, predict(fit, type = "response"))
  expect_equal(unname(rowSums(ft)), rep(1, nrow(dd)), tolerance = 1e-12)
  # cs() really moves the answer: the column is not constant in z
  expect_gt(stats::sd(ft[, 1]), 0.02)

  dna <- dd
  dna$x[c(4L, 30L)] <- NA
  fna <- frm(bf(y ~ x) + cumulative(), data = dna,
             na.action = stats::na.exclude)
  fn <- fitted(fna)
  expect_equal(dim(fn), c(nrow(dna), 3L))
  expect_true(all(is.na(fn[c(4L, 30L), ])))
  expect_false(anyNA(fn[-c(4L, 30L), ]))
  expect_identical(fn, predict(fna, type = "response"))
})

## residuals -----------------------------------------------------------

test_that("ordinal residuals score the categories by their own codes", {
  dd <- ordfit_data(103)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  P <- fitted(fit)
  k <- seq_len(ncol(P))
  m <- as.numeric(P %*% k)
  v <- as.numeric(P %*% (k^2)) - m^2
  yc <- as.integer(dd$y)

  expect_equal(as.numeric(residuals(fit)), yc - m, tolerance = 1e-12)
  expect_equal(as.numeric(residuals(fit, type = "pearson")),
               (yc - m) / sqrt(v), tolerance = 1e-12)
  # a residual on a score, not on the latent scale: the old path
  # returned y - eta, which is a different number entirely
  expect_false(isTRUE(all.equal(as.numeric(residuals(fit)),
                                yc - as.numeric(predict(fit)))))
  # pearson residuals are standardized, so their spread is about 1
  expect_lt(abs(stats::sd(residuals(fit, type = "pearson")) - 1), 0.2)
  # the order-only residuals stay available and stay refused where they
  # were refused
  expect_true(is.numeric(residuals(fit, type = "osa")))
  expect_error(residuals(fit, type = "deviance"), "cumulative")
})

## plot ----------------------------------------------------------------

test_that("plot() runs on an ordinal fit and uses the expected category", {
  dd <- ordfit_data(104, n = 150)
  fit <- frm(bf(y ~ x) + cratio(), data = dd)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(fit, which = 1L, ask = FALSE))
  expect_silent(plot(fit, ask = FALSE))
})

## DHARMa --------------------------------------------------------------

test_that("dharma_residuals() runs and calibrates on an ordinal fit", {
  skip_if_not_installed("DHARMa")
  dd <- ordfit_data(105, n = 400)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  res <- dharma_residuals(fit, nsim = 250, seed = 7)
  expect_s3_class(res, "DHARMa")
  expect_equal(length(res$scaledResiduals), nrow(dd))
  # the fitted axis is the expected category index, on the response's
  # own 1..K scale rather than the latent one
  fpr <- res$fittedPredictedResponse
  expect_true(all(fpr > 1 & fpr < 3))
  expect_equal(as.numeric(fpr),
               as.numeric(fitted(fit) %*% seq_len(3)),
               tolerance = 1e-12)
  # a correctly specified model: the uniformity test must not reject
  u <- DHARMa::testUniformity(res, plot = FALSE)
  expect_gt(u$p.value, 0.01)

  # and it DOES detect a misspecified one: dropping the only predictor
  # leaves the residuals sorted by x
  fit0 <- frm(bf(y ~ 1) + cumulative(), data = dd)
  res0 <- dharma_residuals(fit0, nsim = 250, seed = 7)
  sr <- res0$scaledResiduals
  expect_lt(stats::cor.test(sr, dd$x, method = "spearman")$p.value, 1e-6)
})

## conditional_effects -------------------------------------------------

test_that("conditional_effects() draws one curve per ordinal category", {
  dd <- ordfit_data(106)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  ce <- conditional_effects(fit, resolution = 20L)
  df <- ce[["x"]]
  expect_true("cats__" %in% names(df))
  expect_equal(levels(df$cats__), levels(dd$y))
  expect_equal(nrow(df), 20L * 3L)

  # the curves are the category probabilities predict() gives on the
  # same grid, and they sum to one at every grid point
  grid <- unique(df$x)
  P <- predict(fit, newdata = data.frame(x = grid), type = "response")
  expect_equal(df$estimate__, as.numeric(P), tolerance = 1e-10)
  s <- as.numeric(tapply(df$estimate__, df$x, sum))
  expect_equal(s, rep(1, 20L), tolerance = 1e-10)

  # bands are on the probability scale and cannot leave it
  expect_true(all(df$lower__ > 0 & df$upper__ < 1))
  expect_true(all(df$lower__ < df$estimate__ &
                    df$estimate__ < df$upper__))
  expect_true(all(df$se__ > 0))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(ce, ask = FALSE))
})

test_that("the ordinal effect standard errors include the thresholds", {
  dd <- ordfit_data(107)
  fit <- frm(bf(y ~ x + cs(z)) + sratio(), data = dd)
  rspec <- fit$spec$responses[[1L]]
  lp <- fit$frame$linpreds[[1L]]
  nd <- data.frame(x = c(-1, 0, 1.5), z = c(0.5, 0, -1))
  ed <- frmtmb:::lp_eta_design(fit, lp, nd, FALSE, FALSE)
  ps <- frmtmb:::ord_prob_se(fit, rspec, lp, ed, nd, FALSE)
  expect_equal(ps$P, predict(fit, newdata = nd, type = "response"),
               tolerance = 1e-12, ignore_attr = TRUE)

  # reference: a full numeric Jacobian of the probabilities with
  # respect to EVERY estimated parameter, against the joint covariance
  jc <- frmtmb:::get_joint_cov(fit)
  nms <- jc$names
  th0 <- numeric(length(nms))
  for (nm in unique(nms)) th0[nms == nm] <- fit$estimates[[nm]]
  pf <- function(th) {
    f2 <- fit
    for (nm in unique(nms)) f2$estimates[[nm]] <- th[nms == nm]
    frmtmb:::ord_probs(f2, rspec, nd, FALSE, FALSE)
  }
  ref <- matrix(0, nrow(nd), 3L)
  for (i in seq_len(nrow(nd))) {
    for (kk in seq_len(3L)) {
      g <- vapply(seq_along(th0), function(j) {
        h <- 1e-5 * max(1, abs(th0[j]))
        hi <- th0; hi[j] <- hi[j] + h
        lo <- th0; lo[j] <- lo[j] - h
        (pf(hi)[i, kk] - pf(lo)[i, kk]) / (2 * h)
      }, 0)
      ref[i, kk] <- sqrt(drop(t(g) %*% jc$V %*% g))
    }
  }
  expect_vector_equal(as.vector(ps$se), as.vector(ref), tol = 1e-8)
  # dropping the threshold rows would understate the bands, so they are
  # really contributing
  expect_true(any(grepl("^tau_raw$", nms)))
})

test_that("conditional_effects() keeps the latent route and refuses draws", {
  dd <- ordfit_data(108, n = 150)
  fit <- frm(bf(y ~ x) + acat(), data = dd)
  # naming the dpar opts back into the ordinary latent-scale display
  ce <- conditional_effects(fit, dpar = "mu", resolution = 10L)
  expect_false("cats__" %in% names(ce[["x"]]))
  expect_equal(ce[["x"]]$estimate__,
               predict(fit, newdata = data.frame(x = ce[["x"]]$x),
                       type = "link"),
               tolerance = 1e-10, ignore_attr = TRUE)
  # a prediction interval on a probability curve is not a thing
  expect_error(conditional_effects(fit, method = "predict"),
               "predictive distribution")
})

test_that("an ordinal effect with a second predictor gets a panel each", {
  dd <- ordfit_data(109)
  dd$f <- factor(rep(c("a", "b"), length.out = nrow(dd)))
  fit <- frm(bf(y ~ x + f) + cumulative(), data = dd)
  ce <- conditional_effects(fit, effects = "x:f", resolution = 12L)
  df <- ce[["x:f"]]
  expect_equal(nrow(df), 12L * 2L * 3L)
  expect_setequal(as.character(df$cats__), levels(dd$y))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(ce, ask = FALSE))

  # a discrete varied predictor draws points and error bars instead
  ce2 <- conditional_effects(fit, effects = "f")
  expect_equal(nrow(ce2[["f"]]), 2L * 3L)
  expect_silent(plot(ce2, ask = FALSE))
})

## pp_check ------------------------------------------------------------

test_that("pp_check() compares ordinal draws with the response codes", {
  skip_if_not_installed("bayesplot")
  dd <- ordfit_data(113, n = 150)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  # simulate() returns ordered factors; handing bayesplot the factor
  # labels made it stop with "is.numeric(predictions) is not TRUE"
  expect_s3_class(pp_check(fit, type = "bars", ndraws = 20), "ggplot")
  expect_s3_class(pp_check(fit, ndraws = 5), "ggplot")
})

## posterior_epred array convention ------------------------------------

# v0.31 returned each draw's n x K prediction flattened column by column
# into a draws x (n * K) matrix with "<obs>.<cat>" column names. That is
# not what brms does: ?brms::posterior_epred.brmsfit documents an
# S x N x C array for categorical and ordinal models and an S x N matrix
# otherwise, and frmtmb follows brms spelling for brms-origin functions.

test_that("posterior_epred returns a draws x obs x category array", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  dd <- ordfit_data(114, n = 120)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 11))

  nd <- data.frame(x = c(-1.2, 0, 0.8, 1.5))
  ep <- posterior_epred(ds, newdata = nd, ndraws = 8)
  expect_true(is.array(ep) && length(dim(ep)) == 3L)
  expect_equal(dim(ep), c(8L, 4L, 3L))
  expect_null(dimnames(ep)[[1]])
  expect_equal(dimnames(ep)[[3]], levels(dd$y))
  # ep[, , "hi"] must be addressable by the response's own level names
  expect_equal(dim(ep[, , "hi"]), c(8L, 4L))

  # each draw x observation slice is a distribution
  for (k in seq_len(dim(ep)[1])) {
    expect_vector_equal(rowSums(ep[k, , ]), rep(1, 4), tol = 1e-12)
  }
  expect_true(all(ep > 0 & ep < 1))

  # continuity with the v0.31 flattened matrix: the array is exactly
  # that matrix reshaped, so category slice k is the k-th block of n
  # columns and the old column names still describe the block order
  flat <- matrix(as.vector(ep), nrow = dim(ep)[1])
  expect_equal(dim(flat), c(8L, 12L))
  for (k in seq_len(3L)) {
    expect_identical(as.vector(ep[, , k]),
                     as.vector(flat[, (k - 1L) * 4L + seq_len(4L)]))
  }

  # and each draw's slice is the matrix predict(type = "response")
  # returns, so the posterior mean tracks the MLE probabilities
  P <- predict(fit, newdata = nd, type = "response")
  epf <- posterior_epred(ds, newdata = nd, ndraws = 60)
  expect_lt(max(abs(apply(epf, c(2L, 3L), mean) - P)), 0.1)

  # without newdata the observation margin carries the data rownames
  ep0 <- posterior_epred(ds, ndraws = 4)
  expect_equal(dim(ep0), c(4L, nrow(dd), 3L))
  expect_equal(dimnames(ep0)[[2]], rownames(dd))
  expect_equal(dimnames(ep0)[[3]], levels(dd$y))

  # the other two draws methods are statements about one number per
  # observation and keep the plain draws x observations matrix
  pp <- posterior_predict(ds, newdata = nd, ndraws = 8)
  expect_true(is.matrix(pp))
  expect_equal(dim(pp), c(8L, 4L))
  expect_true(all(pp %in% 1:3))
  pl <- posterior_linpred(ds, newdata = nd, ndraws = 8)
  expect_true(is.matrix(pl))
  expect_equal(dim(pl), c(8L, 4L))

  # pp_check() reads the predictive side, which the change does not move
  skip_if_not_installed("bayesplot")
  expect_s3_class(pp_check(ds, type = "bars", ndraws = 10), "ggplot")
})

test_that("a scalar-response family keeps the draws x obs matrix", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(115)
  gd <- data.frame(x = stats::rnorm(80))
  gd$y <- stats::rpois(80, exp(0.3 + 0.4 * gd$x))
  fp <- frm(bf(y ~ x) + poisson(), data = gd)
  dp <- suppressWarnings(frm_sample(fp, chains = 1, iter = 300,
                                    refresh = 0, seed = 12))
  nd <- data.frame(x = c(-1, 0, 1))
  ep <- posterior_epred(dp, newdata = nd, ndraws = 6)
  expect_true(is.matrix(ep))
  expect_equal(length(dim(ep)), 2L)
  expect_equal(dim(ep), c(6L, 3L))
  expect_null(dimnames(ep))
})

## emmeans / insight ---------------------------------------------------

test_that("emmeans works on the latent scale of an ordinal fit", {
  skip_if_not_installed("emmeans")
  dd <- ordfit_data(110)
  dd$f <- factor(rep(c("a", "b"), length.out = nrow(dd)))
  fit <- frm(bf(y ~ x + f) + cumulative(), data = dd)
  # the ordinal design has no intercept (the thresholds take its place),
  # so the emmeans basis has to be rebuilt by NAME or the reference grid
  # is non-conformable
  em <- emmeans::emmeans(fit, ~ f)
  s <- data.frame(summary(em))
  expect_equal(nrow(s), 2L)
  expect_true(all(is.finite(s$emmean)))
  # the contrast is the fixed-effect difference, on the latent scale
  ct <- data.frame(summary(pairs(em)))
  expect_equal(ct$estimate[1], -unname(fixef(fit)$mu["fb"]),
               tolerance = 1e-6)
})

test_that("insight accessors survive the ordinal fitted() change", {
  skip_if_not_installed("insight")
  dd <- ordfit_data(111, n = 150)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  expect_equal(insight::find_response(fit), "y")
  gp <- insight::get_predicted(fit)
  # insight stays on the latent predictor, the clm-like convention
  expect_equal(length(as.numeric(gp)), nrow(dd))
  expect_equal(as.numeric(gp), unname(predict(fit, type = "link")),
               tolerance = 1e-8)
  expect_s3_class(insight::get_parameters(fit), "data.frame")
})

test_that("marginaleffects still keys ordinal categories by group", {
  skip_if_not_installed("marginaleffects")
  dd <- ordfit_data(112, n = 150)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  p <- marginaleffects::avg_predictions(fit)
  expect_equal(nrow(p), 3L)
  expect_setequal(as.character(p$group), levels(dd$y))
  expect_equal(sum(p$estimate), 1, tolerance = 1e-8)
})
