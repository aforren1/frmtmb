# The variance of a one-scale dense block where exp() underflows (lane
# surface, 2026-10-06; the gpby review's frm_sample(fit) on an exact
# y ~ gp(x) fit, which did not move).
#
# Seen to fail on 0.68.1 (rellib-r6): at log sd -1137.64, the point the
# sampler's first long leapfrog step reached (dev/surface-gp-diag4.R),
# the exact gp() field's log density was +Inf: exp(2 * log sd) is 0,
# the covariance is the zero matrix, and RTMB's dmvnorm() of a nonzero
# field against it returned +Inf. The chain read that as a point of
# infinite density and never left it: accept 0, stepsize NaN, 300 of
# 300 transitions divergent (dev/surface-sample-repros.R item 9).

test_that("the floor is the variance itself wherever a fit can be", {
  th <- seq(-326, 50, length.out = 200001)
  expect_identical(frmtmb:::sd2_floored(th), exp(2 * th))
  # and the floor where exp() underflows
  expect_identical(frmtmb:::sd2_floored(-1137.64), 1e-300)
  expect_identical(frmtmb:::sd2_floored(-Inf), 1e-300)
})

gp_fit <- function() {
  set.seed(5)
  n <- 60
  d <- data.frame(x = round(stats::runif(n, 0, 6), 1))
  d$y <- 0.5 + sin(d$x) + stats::rnorm(n, 0, 0.3)
  frm(bf(y ~ gp(x)), family = gaussian(), data = d)
}

test_that("an exact gp() field has no +Inf log density at a tiny sd", {
  fit <- gp_fit()
  blk <- fit$frame$re_blocks[[1L]]
  expect_identical(blk$covstruct, "gp")
  nll <- frmtmb:::covstruct_registry[["gp"]]$nll
  b <- fit$estimates$b
  th <- fit$estimates$theta
  # the control: at the estimates the density is finite
  expect_true(is.finite(nll(b, th, blk)))
  for (lsd in c(-400, -1137.64, -5000)) {
    v <- nll(b, c(lsd, th[2L]), blk)
    expect_false(isTRUE(v == Inf), info = paste("log sd", lsd))
    expect_lt(v, nll(b, th, blk))
  }
})

test_that("the floor leaves an exact gp() fit's objective unchanged", {
  fit <- gp_fit()
  blk <- fit$frame$re_blocks[[1L]]
  b <- fit$estimates$b
  th <- fit$estimates$theta
  # the density as 0.68.1 wrote it, beside the one the fit uses
  old <- sum(RTMB::dmvnorm(t(matrix(b, nrow = 1L)), 0,
                           exp(2 * th[1L]) * frmtmb:::gp_corr(th, blk),
                           log = TRUE))
  expect_identical(frmtmb:::covstruct_registry[["gp"]]$nll(b, th, blk), old)
})

test_that("the other one-scale dense blocks take the same floor", {
  # test-v11.R's construction, which converges
  set.seed(2)
  n_g <- 50
  times <- c(0, 0.3, 0.4, 1.1, 1.5, 2.7)
  Sig <- 0.9^2 * exp(-1.2 * abs(outer(times, times, "-")))
  U <- matrix(stats::rnorm(n_g * 6), n_g) %*% chol(Sig)
  d <- data.frame(y = 1 + as.vector(t(U)) + stats::rnorm(n_g * 6, 0, 0.4),
                  g = factor(rep(seq_len(n_g), each = 6)),
                  tim = num_factor(rep(times, n_g)))
  fit <- frm(bf(y ~ 1 + ou(tim + 0 | g)), family = gaussian(), data = d)
  blk <- fit$frame$re_blocks[[1L]]
  nll <- frmtmb:::covstruct_registry[[blk$covstruct]]$nll
  b <- fit$estimates$b + 0.1
  th <- fit$estimates$theta
  expect_true(is.finite(nll(b, th, blk)))
  expect_false(isTRUE(nll(b, c(-1137.64, th[-1L]), blk) == Inf))
})

test_that("every structure the floor touches has no +Inf at a tiny sd", {
  # the six nll()s that take sd2_floored(): gp, ou, homcs, homtoep, the
  # spatial entries (exp, gau, mat) and gr(cov =) with one coefficient.
  # Each block is read from its frame at the start values with the
  # field moved off zero, and its log sd put where exp() underflows.
  # Seen on 0.68.1 (rellib-r6, dev/surface-p1-tests/floor-base.log):
  # all eight gave +Inf there.
  set.seed(3)
  n_g <- 12
  dd <- data.frame(g = factor(rep(seq_len(n_g), each = 5)),
                   tim = factor(rep(1:5, n_g)),
                   tn = num_factor(rep(c(0, 0.3, 0.4, 1.1, 1.5), n_g)),
                   pos = num_factor(rep(round(stats::runif(5), 2) * 10,
                                        n_g),
                                    rep(round(stats::runif(5), 2) * 10,
                                        n_g)),
                   x = stats::runif(5 * n_g))
  dd$y <- stats::rnorm(5 * n_g)
  A <- diag(n_g) * 0.5 + 0.5
  dimnames(A) <- list(levels(dd$g), levels(dd$g))
  forms <- list(
    gp = y ~ gp(x),
    ou = y ~ 1 + ou(tn + 0 | g),
    homcs = y ~ 1 + homcs(tim + 0 | g),
    homtoep = y ~ 1 + homtoep(tim + 0 | g),
    exp = y ~ 1 + exp(pos + 0 | g),
    gau = y ~ 1 + gau(pos + 0 | g),
    mat = y ~ 1 + mat(pos + 0 | g),
    gr_cov = y ~ 1 + (1 | gr(g, cov = A)))
  seen <- character()
  for (nm in names(forms)) {
    fr <- frm(forms[[nm]], family = gaussian(), data = dd,
              data2 = list(A = A), dry_run = "frame")
    blk <- fr$re_blocks[[1L]]
    seen <- c(seen, blk$covstruct)
    nll <- frmtmb:::covstruct_registry[[blk$covstruct]]$nll
    tpl <- fr$par_template
    b <- tpl$b + seq_along(tpl$b) / length(tpl$b)
    th <- tpl$theta
    v0 <- nll(b, th, blk)
    expect_true(is.finite(v0), info = nm)
    v <- nll(b, c(-1137.64, th[-1L]), blk)
    expect_false(isTRUE(v == Inf), info = nm)
    expect_lt(v, v0)
  }
  expect_setequal(unique(seen), c("gp", "ou", "homcs", "homtoep", "exp",
                                  "gau", "mat", "gr_cov"))
})
