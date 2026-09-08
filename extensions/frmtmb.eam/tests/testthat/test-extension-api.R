## The one thing this package promises another extension package.
##
## `wiener_lpdf()` exists so that frmtmb.learn::rlddm() can reach the
## Wiener density without a colon. What has to hold for that promise to
## be worth making is here: it is the same function wiener() evaluates,
## it tapes, and it is checked against something outside this package.

test_that("wiener_lpdf is the density wiener() itself evaluates", {
  # the family's own lpdf, reached the way frm() reaches it, against the
  # export. Same numbers or the export is documenting a different
  # function from the one the package fits with.
  fam <- wiener()
  t <- c(0.35, 0.7, 1.4, 3.0)
  ndt <- 0.15
  dp <- list(mu = 1.2, bs = 1.5, ndt = ndt, bias = 0.3)
  for (up in c(0, 1)) {
    theirs <- fam$lpdf(t, dp, list(dec = rep(up, length(t))))
    ours <- wiener_lpdf(t - ndt, dp$mu, dp$bs, dp$bias, up)
    expect_equal(as.numeric(theirs), as.numeric(ours), tolerance = 1e-14)
  }
})

test_that("wiener_lpdf integrates to the diffusion's exit probability", {
  # an outside check: the probability of ever reaching a boundary is a
  # closed form, and the density over that boundary must integrate to it
  v <- 1.2
  a <- 1.5
  w <- 0.3
  t <- seq(1e-5, 12, length.out = 200000)
  dt <- diff(t)[[1L]]
  p_lo_analytic <- (exp(-2 * v * a) - exp(-2 * v * a * w)) /
    (exp(-2 * v * a) - 1)
  expect_equal(sum(exp(wiener_lpdf(t, v, a, w, 0))) * dt, p_lo_analytic,
               tolerance = 1e-6)
  expect_equal(sum(exp(wiener_lpdf(t, v, a, w, 1))) * dt,
               1 - p_lo_analytic, tolerance = 1e-6)
})

test_that("wiener_lpdf agrees with RWiener on both boundaries", {
  skip_if_not_installed("RWiener")
  x <- c(0.35, 0.7, 1.4, 3.0)
  a <- 1.5
  v <- 1.2
  w <- 0.3
  t0 <- 0.1
  for (resp in c("upper", "lower")) {
    # resp repeated to x's length: RWiener warns rather than errors when
    # it has to recycle, and a warning in a suite is noise that hides
    # the next one
    ref <- RWiener::dwiener(x, alpha = a, tau = t0, beta = w, delta = v,
                            resp = rep(resp, length(x)), give_log = TRUE)
    got <- wiener_lpdf(x - t0, v, a, w, as.numeric(resp == "upper"))
    expect_equal(as.numeric(got), as.numeric(ref), tolerance = 1e-13)
  }
})

test_that("wiener_lpdf differentiates on a tape", {
  # the whole reason the density is written the way it is. A numeric
  # answer alone would not say the export is usable by a family.
  tp <- RTMB::MakeTape(function(p) {
    sum(wiener_lpdf(c(0.4, 0.9, 1.6), p[1], p[2], p[3], c(1, 0, 1)))
  }, c(1.2, 1.5, 0.3))
  j <- tp$jacobian(c(1.2, 1.5, 0.3))
  expect_length(as.numeric(j), 3L)
  expect_true(all(is.finite(j)))
  # against a finite difference, so the tape is checked rather than
  # merely exercised
  f <- function(p) {
    sum(wiener_lpdf(c(0.4, 0.9, 1.6), p[[1L]], p[[2L]], p[[3L]],
                    c(1, 0, 1)))
  }
  p <- c(1.2, 1.5, 0.3)
  fd <- vapply(seq_along(p), function(k) {
    h <- 1e-6
    pu <- p
    pd <- p
    pu[[k]] <- pu[[k]] + h
    pd[[k]] <- pd[[k]] - h
    (f(pu) - f(pd)) / (2 * h)
  }, 0)
  expect_equal(as.numeric(j), fd, tolerance = 1e-6)
})

test_that("wiener_lpdf refuses a boundary indicator that is not 0 or 1", {
  expect_error(wiener_lpdf(0.5, 1, 1.5, 0.5, 2),
               "says which boundary")
  expect_error(wiener_lpdf(0.5, 1, 1.5, 0.5, NA),
               "says which boundary")
})

test_that("a decision time at or below zero is -Inf and not NaN", {
  got <- wiener_lpdf(c(-1, 0, 0.5), 1.2, 1.5, 0.5, 1)
  expect_true(all(is.infinite(got[1:2])))
  expect_true(all(got[1:2] < 0))
  expect_true(is.finite(got[[3L]]))
})
