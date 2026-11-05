# An exact gp()'s nugget is white noise of sd 1e-3 * sdgp on the field, one
# draw per position. The likelihood cannot tell it from the noise at that
# position, so once its variance reaches that noise frm() and brms read different
# models, and on the base build frm() mostly said nothing
# (dev/gpgap-findings.md, dev/reviews/2026-11-02-gpgap.md M2, M5).

gp_lownoise <- function(s, n = 60, rep = 1L, a = 1) {
  set.seed(20261102)
  x <- round(stats::runif(n %/% rep, 0, 6), 3)
  x <- rep(x, each = rep)
  x <- x[sample(length(x))]
  data.frame(x = x, y = a * sin(2 * x) + a * s * stats::rnorm(length(x)),
             z = round(stats::runif(length(x), 0, 6), 3),
             g = factor(rep(seq_len(6), length.out = length(x))))
}

# The marginal of y integrates the field in closed form: a covariance
# conditioned by sigma^2, whatever the kernel's conditioning, so it reads
# either nugget without a tape.
gp_marginal_nll <- function(p, d, nugget) {
  u <- sort(unique(d$x))
  j <- match(d$x, u)
  s2 <- exp(2 * p[2])
  k <- s2 * exp(-outer(u, u, "-")^2 / (2 * exp(2 * p[3])))
  diag(k) <- diag(k) + nugget(s2)
  v <- k[j, j]
  diag(v) <- diag(v) + exp(2 * p[4])
  r <- chol(v)
  z <- backsolve(r, d$y - p[1], transpose = TRUE)
  0.5 * sum(z^2) + sum(log(diag(r))) + 0.5 * nrow(d) * log(2 * pi)
}
at_frm <- function(s2) 1e-6 * s2
at_brms <- function(s2) 1e-12

# frm() with the warnings of the regime let through: a fit that the guard reads
# stops there, so an SE check's own note on a sigma at its floor is allowed
quiet_fit <- function(expr) {
  allow_warnings(expr, c("Standard errors are not available",
                         "Optimizer did not report convergence",
                         "not small against the noise"))
}
fires <- function(fit) !is.null(gp_nugget_message(fit))
gp_blocks <- function(fit) {
  Filter(function(b) identical(b$covstruct, "gp"), fit$frame$re_blocks)
}

test_that("a nugget that reaches the noise warns by class, and it is real", {
  d <- gp_lownoise(1e-3)
  seen <- NULL
  fit <- withCallingHandlers(
    frm(bf(y ~ gp(x)), family = gaussian(), data = d),
    warning = function(w) {
      if (inherits(w, "frmtmb_gp_nugget")) seen <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(seen, c("frmtmb_gp_nugget", "frmtmb_warning"))
  expect_match(conditionMessage(seen), "gp(x)", fixed = TRUE)
  # says what moves, and does not promise the remedy
  expect_match(conditionMessage(seen), "an estimated sigma reads low")
  expect_match(conditionMessage(seen), "sdgp and the length scale read low")
  expect_match(conditionMessage(seen), "k must be large enough")
  # tryCatch reaches it by name
  hit <- tryCatch(frm(bf(y ~ gp(x)), family = gaussian(), data = d),
                  frmtmb_gp_nugget = function(w) "caught")
  expect_identical(hit, "caught")

  # frm() is the maximum of the closed form at its own nugget; brms's nugget has
  # a maximum above it, by more than a quarter of the one-parameter 95 percent
  # boundary, with a sigma that is another order of magnitude
  bk <- gp_blocks(fit)[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  p0 <- c(mean(d$y), th[[1]], th[[2]], log(1e-3))
  a <- stats::nlminb(p0, gp_marginal_nll, d = d, nugget = at_frm)
  b <- stats::nlminb(p0, gp_marginal_nll, d = d, nugget = at_brms)
  q <- stats::qchisq(0.95, 1) / 8
  expect_gt(a$objective - b$objective, q)
  expect_lt(abs(fit$opt$objective - a$objective), q)
  expect_lt(as.numeric(stats::sigma(fit)), 0.1 * exp(b$par[4]))
})

test_that("the noise at a position is the likelihood's curvature", {
  # gaussian with one row per position, where it is sigma^-2 exactly, and with
  # m rows per position, where it is m sigma^-2. The bound is the rounding of the
  # difference of two large numbers, the prior's diagonal at cond(K) * eps.
  for (rep in c(1L, 5L)) {
    d <- gp_lownoise(0.1, n = 60, rep = rep)
    fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
    bk <- gp_blocks(fit)[[1L]]
    info <- gp_position_info(fit, bk)
    want <- rep / as.numeric(stats::sigma(fit))^2
    expect_length(info, 60L %/% rep)
    v <- cs_entry("gp")$vcov(fit$estimates$theta[bk$theta_idx], bk)
    bound <- kappa(v, exact = TRUE) * .Machine$double.eps *
      max(diag(chol2inv(chol(v)))) / want
    expect_lt(max(abs(info / want - 1)), bound)
  }
})

test_that("the threshold is one standard error of sigma, either side", {
  # m = 60: the rule's sqrt(2 / m) is 0.183 of sigma^2 at nugget 1e-6 sdgp^2,
  # about sigma / sdgp = 2.3e-3. sdgp is near 2 on this curve.
  low <- frm(bf(y ~ gp(x)), family = gaussian(), data = gp_lownoise(8e-3))
  up <- quiet_fit(frm(bf(y ~ gp(x)), family = gaussian(),
                      data = gp_lownoise(3e-3)))
  ratio <- function(fit) {
    bk <- gp_blocks(fit)[[1L]]
    1e-6 * exp(2 * fit$estimates$theta[bk$theta_idx][[1L]]) /
      as.numeric(stats::sigma(fit))^2
  }
  thr <- sqrt(2 / 60)
  expect_lt(ratio(low), thr)
  expect_gt(ratio(up), thr)
  expect_false(fires(low))
  expect_true(fires(up))
  # a factor of 2 in the rule moves it past one of these rows
  near_low <- frm(bf(y ~ gp(x)), family = gaussian(),
                  data = gp_lownoise(5.5e-3))
  near_up <- quiet_fit(frm(bf(y ~ gp(x)), family = gaussian(),
                           data = gp_lownoise(4.2e-3)))
  expect_gt(ratio(near_low), 0.5 * thr)
  expect_lt(ratio(near_low), thr)
  expect_gt(ratio(near_up), thr)
  expect_lt(ratio(near_up), 2 * thr)
  expect_false(fires(near_low))
  expect_true(fires(near_up))
})

test_that("an ordinary field says nothing", {
  # the guarded thing is absent, at sigma 0.3 and 0.1 of the curve, and on the
  # trends where sigma is right to a few percent and the maxima coincide
  for (s in c(0.3, 0.1)) {
    fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = gp_lownoise(s))
    expect_false(fires(fit))
  }
  set.seed(5701)
  x <- round(stats::runif(100, 0, 10), 3)
  e <- stats::rnorm(100)
  for (f in list(x, 100 * stats::plogis(x - 5))) {
    fit <- frm(bf(y ~ gp(x)), family = gaussian(),
               data = data.frame(x = x, y = f + 0.3 * e))
    expect_false(fires(fit))
  }
  # one low row in 60, a group of one, is not a low-noise fit
  d <- gp_lownoise(0.1)
  d$w <- 0.1
  d$w[1L] <- 1e-5
  fit <- frm(bf(y ~ gp(x), sigma ~ 0 + offset(log(w))), data = d)
  expect_false(fires(fit))
})

test_that("a share of the positions at low noise is read, one row is not", {
  # sigma 1e-3 on 15 rows of 60: the median of the ratio does not see them and
  # the field moves; the row of the test above does not move it
  d <- gp_lownoise(1)
  set.seed(5601)
  e <- stats::rnorm(nrow(d))
  d$low <- factor(ifelse(seq_len(nrow(d)) <= 15L, "a", "b"))
  d$y <- sin(2 * d$x) + ifelse(d$low == "a", 1e-3, 0.3) * e
  fit <- quiet_fit(frm(bf(y ~ gp(x), sigma ~ low), data = d))
  expect_true(fires(fit))
})

test_that("a trend on ordinary noise that the floor moves warns", {
  # the regime is the field's range against the noise, not noise-free data
  set.seed(5701)
  x <- round(stats::runif(100, 0, 10), 3)
  d <- data.frame(x = x, y = x^2 + 0.3 * stats::rnorm(100))
  fit <- quiet_fit(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
  expect_true(fires(fit))
})

test_that("sigma with no coefficient of its own: the fit stands, offset read", {
  # the first form of this diagnostic threw on the first spelling, and lost a
  # fit that was complete; the last has a coefficient and an offset, and it read
  # the coefficient alone
  d <- gp_lownoise(1e-3)
  d$wsm <- 1e-3
  fit <- quiet_fit(frm(bf(y ~ gp(x), sigma ~ 0 + offset(log(wsm))), data = d))
  expect_s3_class(fit, "frmtmb_fit")
  expect_true(fires(fit))
  d$wsm <- 0.3
  d$y <- sin(2 * d$x) + 0.3 * stats::rnorm(nrow(d))
  expect_false(fires(frm(bf(y ~ gp(x), sigma ~ 0 + offset(log(wsm))), data = d)))
  d <- gp_lownoise(1e-3)
  d$wlow <- 1e-3
  fit <- quiet_fit(frm(bf(y ~ gp(x), sigma ~ 1 + offset(log(wlow))), data = d))
  expect_true(fires(fit))
})

test_that("a failure to read the tape cannot lose the fit", {
  d <- gp_lownoise(0.1)
  local_mocked_bindings(gp_position_info = function(fit, bk) stop("no tape"))
  fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
  expect_s3_class(fit, "frmtmb_fit")
  expect_null(gp_nugget_message(fit))
})

test_that("every route with noise at a position is read", {
  noise <- 1e-3
  d <- gp_lownoise(noise)
  rows <- list(
    # known noise: the noise is not biased, the field is
    se = quote(frm(bf(y | se(se) ~ gp(x)), data = cbind(d, se = noise))),
    re_sigma = quote(frm(bf(y ~ gp(x), sigma ~ (1 | g)), data = d)),
    smooth_sigma = quote(frm(bf(y ~ gp(x), sigma ~ s(z, k = 5)), data = d)),
    fixed_sigma = quote(frm(bf(y ~ gp(x), sigma ~ z), data = d)),
    mixture = quote(frm(bf(y ~ 1, mu1 ~ gp(x), mu2 ~ 1),
                        family = mixture(gaussian(), gaussian()), data = d)),
    nl = quote(frm(bf(y ~ b1 + b2 * z, b1 ~ gp(x), b2 ~ 1, nl = TRUE),
                   data = d)),
    student = quote(frm(bf(y ~ gp(x)), family = student(), data = d)))
  for (nm in names(rows)) {
    fit <- suppressWarnings(eval(rows[[nm]]))
    expect_true(fires(fit), info = nm)
  }
})

test_that("a family's own information is read, on both sides", {
  set.seed(5601)
  n <- 60
  x <- round(stats::runif(n, 0, 6), 3)
  for (lm in c(4, 16)) {
    d <- data.frame(x = x, cnt = stats::rpois(n, exp(lm + sin(2 * x))))
    fit <- suppressWarnings(frm(bf(cnt ~ gp(x)), family = poisson(), data = d))
    expect_identical(fires(fit), lm == 16, info = lm)
  }
  for (tr in c(1e3, 1e7)) {
    d <- data.frame(x = x, tr = tr)
    d$suc <- stats::rbinom(n, d$tr, stats::plogis(sin(2 * x)))
    fit <- suppressWarnings(frm(bf(suc | trials(tr) ~ gp(x)),
                                family = binomial(), data = d))
    expect_identical(fires(fit), tr == 1e7, info = tr)
  }
})

test_that("replicated positions use the replication", {
  # 20 positions of 30 rows at sigma 5e-3 has r = nugget / (sigma^2 / 30): sigma
  # is right, and sdgp moves. Unreplicated at nearly that sigma stays quiet.
  rp <- gp_lownoise(5e-3, n = 600, rep = 30L)
  fit <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = rp))
  expect_true(fires(fit))
  un <- gp_lownoise(8e-3, n = 60)
  expect_false(fires(frm(bf(y ~ gp(x)), family = gaussian(), data = un)))
})

test_that("each sub-GP is read on its own, and only the low one named", {
  d <- gp_lownoise(0.01 / 3)
  set.seed(20261103)
  # sigma 0.01: the field in x of sd about 3 has noise 0.3 of sigma, the field in
  # z of sd about 0.1 has noise 0.01 of it
  d$z <- round(stats::runif(nrow(d), 0, 6), 3)
  d$y <- 3 * d$y + 0.1 * sin(2 * d$z)
  # either way round: every block is read
  for (fo in list(y ~ gp(x) + gp(z), y ~ gp(z) + gp(x))) {
    fit <- quiet_fit(frm(do.call(bf, list(fo)), family = gaussian(), data = d))
    msg <- gp_nugget_message(fit)
    expect_match(msg, "gp(x)", fixed = TRUE)
    expect_false(grepl("gp(z)", msg, fixed = TRUE))
  }
})

test_that("a field with no nugget, or one the noise dwarfs, is left alone", {
  d <- gp_lownoise(1e-3)
  # the Hilbert-space form has no nugget, and at this k it finds brms's model
  # within one standard error of log sigma
  fk <- frm(bf(y ~ gp(x, k = 30)), family = gaussian(), data = d)
  expect_false(fires(fk))
  s_brms <- stats::nlminb(c(mean(d$y), 0.5, 0.2, log(1e-3)), gp_marginal_nll,
                          d = d, nugget = at_brms)$par[4]
  expect_lt(abs(log(as.numeric(stats::sigma(fk))) - s_brms),
            1 / sqrt(2 * nrow(d)))
  # a gp() of sigma itself, whose curvature is about 2 per row
  d$w <- stats::rnorm(nrow(d), 0, exp(-1 + sin(2 * d$x)))
  fs <- frm(bf(w ~ 1, sigma ~ gp(x)), family = gaussian(), data = d)
  expect_false(fires(fs))
})

test_that("a block carries its nugget to the density and the kriging", {
  d <- gp_lownoise(0.3)
  fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
  bk <- gp_blocks(fit)[[1L]]
  expect_identical(gp_block_nugget(bk), gp_nugget)
  th <- fit$estimates$theta[bk$theta_idx]
  bk2 <- bk
  bk2$gp_nugget <- 1e-12
  expect_identical(gp_block_nugget(bk2), 1e-12)
  sd2 <- exp(2 * th[[1L]])
  expect_equal(unname(diag(cs_entry("gp")$vcov(th, bk2))),
               rep(sd2 * (1 + 1e-12), bk$dim))
  pos <- matrix(sort(unique(d$x)))
  cross <- gp_cross_cov(th, bk2, pos[1L, , drop = FALSE], pos)
  expect_equal(cross[1L, 1L], sd2 * (1 + 1e-12))

  # an unseen position: the kriging variance carries the block's nugget
  fit2 <- fit
  ib <- which(vapply(fit2$frame$re_blocks, function(b) {
    identical(b$covstruct, "gp")
  }, NA))
  fit2$frame$re_blocks[[ib]]$gp_nugget <- 1e-12
  nd <- data.frame(x = 2000)
  v1 <- frm_lp_basis(fit, newdata = nd)$extra_var
  v2 <- frm_lp_basis(fit2, newdata = nd)$extra_var
  # far from the data the kriging takes nothing off the prior variance, so the
  # difference is the two nuggets, against rounding of sd^2
  expect_equal(as.numeric(v1 - v2), sd2 * (1e-6 - 1e-12),
               tolerance = 1e4 * .Machine$double.eps / 1e-6)

  # a block's factor order moves no number the kriging reports, only its rounding
  ord <- rev(seq_len(bk$dim))
  fit3 <- fit
  fit3$frame$re_blocks[[ib]]$gp_order <- ord
  nd2 <- data.frame(x = c(1.2345, 3.21, 2000))
  a1 <- frm_lp_basis(fit, newdata = nd2)
  a3 <- frm_lp_basis(fit3, newdata = nd2)
  scale_xw <- max(abs(a1$A))
  expect_lt(max(abs(a1$A - a3$A)), 1e-4 * scale_xw)
  expect_lt(max(abs(a1$extra_var - a3$extra_var)), 1e-6 * sd2)
  # the factor: a row-permuted triangle whose product is the kernel
  b3 <- fit3$frame$re_blocks[[ib]]
  l3 <- cs_entry("gp")$chol_L(th, b3)
  v3 <- unname(as.matrix(cs_entry("gp")$vcov(th, b3)))
  expect_lt(max(abs(tcrossprod(l3) - v3)), 1e-10 * sd2)
  lp <- l3[ord, ]
  expect_lt(max(abs(lp[upper.tri(lp)])), 1e-12 * sqrt(sd2))
})

test_that("a numeric kernel past LAPACK's reach factors, barely jittered", {
  # two rows at one position leave K singular in double: the factor a draw needs
  # comes back, with the smallest jitter that works, and a factor LAPACK can take
  # is not touched
  x <- c(0, 0, 1, 2)
  k <- exp(-outer(x, x, "-")^2 / 2)
  expect_error(chol(k))
  r <- gp_chol_num(k)
  expect_lt(max(abs(crossprod(r) - k)), 10 * gp_nugget)
  ok <- exp(-outer(x[-1L], x[-1L], "-")^2 / 2) + diag(1e-6, 3L)
  expect_identical(gp_chol_num(ok), chol(ok))
  # past what the jitter may reach it stops, as LAPACK does: no jitter of at most
  # the nugget makes a negative eigenvalue
  bad <- ok
  bad[1L, 2L] <- bad[2L, 1L] <- 2
  expect_error(gp_chol_num(bad))
})

