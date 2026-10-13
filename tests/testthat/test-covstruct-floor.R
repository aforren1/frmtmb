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

# The blocks with a scale per coefficient (dev/moparity-findings.md,
# item 5). Seen on 0.69.0 (rellib-r7, dev/moparity-floor-sweep.R): at
# log sd -1137.64 on the tape, us with two and three coefficients, cs,
# toep and gr(cov =) with two gave +Inf, ar1, hetar1, us_t with two and
# gr(prec =) with two gave NaN.
multi_scale_frames <- function() {
  set.seed(3)
  n_g <- 12
  dd <- data.frame(g = factor(rep(seq_len(n_g), each = 5)),
                   tim = factor(rep(1:5, n_g)),
                   x = stats::runif(5 * n_g), z = stats::rnorm(5 * n_g))
  dd$y <- stats::rnorm(5 * n_g)
  A <- diag(n_g) * 0.5 + 0.5
  dimnames(A) <- list(levels(dd$g), levels(dd$g))
  Q <- solve(A)
  forms <- list(
    us2 = y ~ 1 + (1 + x | g),
    us3 = y ~ 1 + (1 + x + z | g),
    ust2 = y ~ 1 + (1 + x | gr(g, dist = "student")),
    ar1 = y ~ 1 + ar1(tim + 0 | g),
    hetar1 = y ~ 1 + hetar1(tim + 0 | g),
    cs = y ~ 1 + cs(tim + 0 | g),
    toep = y ~ 1 + toep(tim + 0 | g),
    grcov2 = y ~ 1 + (1 + x | gr(g, cov = A)),
    grprec2 = y ~ 1 + (1 + x | gr(g, prec = Q)))
  lapply(forms, function(f) {
    frm(f, family = gaussian(), data = dd, data2 = list(A = A, Q = Q),
        dry_run = "frame")
  })
}

# the block's log density on RTMB's tape, the one a sampler evaluates
block_tape <- function(fr, nll = NULL) {
  blk <- fr$re_blocks[[1L]]
  nll <- nll %||% frmtmb:::covstruct_registry[[blk$covstruct]]$nll
  tpl <- fr$par_template
  b <- tpl$b + seq_along(tpl$b) / length(tpl$b)
  th <- tpl$theta[blk$theta_idx]
  obj <- RTMB::MakeADFun(function(p) -nll(p$b, p$theta, blk),
                         list(b = b, theta = th), silent = TRUE)
  list(obj = obj, b = b, th = th, blk = blk)
}

test_that("a block with several scales has no +Inf or NaN at a tiny sd", {
  frs <- multi_scale_frames()
  seen <- character()
  for (nm in names(frs)) {
    tp <- block_tape(frs[[nm]])
    reg <- frmtmb:::covstruct_registry[[tp$blk$covstruct]]
    seen <- c(seen, tp$blk$covstruct)
    v0 <- -tp$obj$fn(c(tp$b, tp$th))
    expect_true(is.finite(v0), info = nm)
    sdi <- reg$sd_idx(tp$blk$dim)
    expect_gt(length(sdi), 0L)
    for (k in c(as.list(sdi), list(sdi))) {
      t2 <- tp$th
      t2[k] <- -1137.64
      v <- -tp$obj$fn(c(tp$b, t2))
      lab <- paste(nm, "log sd", paste(k, collapse = ","))
      expect_false(is.nan(v), info = lab)
      expect_false(isTRUE(v == Inf), info = lab)
      expect_lt(v, v0)
    }
  }
  expect_setequal(unique(seen), c("us", "us_t", "ar1", "hetar1", "cs",
                                  "toep", "gr_cov", "gr_prec"))
})

test_that("the per-scale floor leaves the tape unchanged above it", {
  # the densities as 0.69.0 wrote them: the same functions with exp()
  # where the floor is, so the comparison reads no copy of the code
  ns <- asNamespace("frmtmb")
  old_env <- new.env(parent = ns)
  old_env$sd_floored <- exp
  old_us_sigma <- frmtmb:::us_sigma
  environment(old_us_sigma) <- old_env
  old_env$us_sigma <- old_us_sigma
  # the inverse, built in: the old density does differ at the floor
  frs <- multi_scale_frames()
  for (nm in names(frs)) {
    blk <- frs[[nm]]$re_blocks[[1L]]
    new_nll <- frmtmb:::covstruct_registry[[blk$covstruct]]$nll
    old_nll <- new_nll
    environment(old_nll) <- old_env
    tn <- block_tape(frs[[nm]], new_nll)
    to <- block_tape(frs[[nm]], old_nll)
    sdi <- frmtmb:::covstruct_registry[[blk$covstruct]]$sd_idx(blk$dim)
    for (lsd in c(seq(-308.4, 4, length.out = 41), tn$th[sdi[1L]])) {
      t2 <- tn$th
      t2[sdi] <- lsd
      p <- c(tn$b, t2)
      expect_identical(tn$obj$fn(p), to$obj$fn(p), info = paste(nm, lsd))
      expect_identical(tn$obj$gr(p), to$obj$gr(p), info = paste(nm, lsd))
    }
    t2 <- tn$th
    t2[sdi] <- -1137.64
    p <- c(tn$b, t2)
    old <- to$obj$fn(p)
    expect_true(is.nan(old) || isTRUE(old == -Inf), info = nm)
    expect_false(identical(tn$obj$fn(p), old), info = nm)
  }
})
