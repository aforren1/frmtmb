# asym_laplace's log density is linear in mu on each side of its mode,
# so with mu linear in a random effect the data give the effect no
# curvature: the Laplace approximation is max_b (l(b) - b' Sigma^-1 b /
# 2), which never falls as the block's scale grows. Through 0.82.0 such
# fits ended at whatever scale TMB's inner search left (all 40 of the
# review's fits at "false convergence", the approximation at exact
# inner modes 2.9 to 49 above the reported logLik in median,
# dev/reviews/wcheck2-review-al.R); frm() now refuses them.

ald_data <- function() {
  set.seed(20261114)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n),
                  g = factor(rep(1:20, each = 4)))
  d$y <- 1 + 0.5 * d$x + sin(3 * d$z) +
    stats::rnorm(20, 0, 0.7)[as.integer(d$g)] + stats::rnorm(n)
  d$yp <- exp(d$y / 3)
  d$yz <- ifelse(seq_len(n) %% 5 == 0, 0, d$y)
  q8 <- stats::quantile(d$y, 0.85)
  d$cc <- ifelse(d$y > q8, "right", "none")
  d$yc <- pmin(d$y, q8)
  d
}

test_that("asym_laplace with a random effect in mu is refused", {
  d <- ald_data()
  al <- asym_laplace()
  for (fo in list(bf(y ~ x + (1 | g)), bf(y ~ x + (1 + x | g)),
                  bf(y ~ x + (1 + x || g)), bf(y ~ x + s(z)))) {
    err <- tryCatch(frm(fo, family = al, data = d), error = function(e) e)
    expect_s3_class(err, "frmtmb_ald_no_laplace")
    expect_match(conditionMessage(err), "frm_sample()", fixed = TRUE)
  }
  err <- tryCatch(frm(y ~ x + s(z), family = al, data = d),
                  error = conditionMessage)
  expect_match(err, "frm() cannot fit `s(z)` under the Laplace",
               fixed = TRUE)
  expect_error(frm(bf(yz ~ x + (1 | g)),
                   family = zero_inflated_asym_laplace(), data = d),
               class = "frmtmb_ald_no_laplace")
  # a prior on b or sigma leaves the scale free: the approximation still
  # never falls in it
  expect_error(frm(bf(y ~ x + (1 | g)), family = al, data = d,
                   prior = set_prior("normal(0, 1)", class = "b")),
               class = "frmtmb_ald_no_laplace")
})

test_that("the refusal stays open where the proof does not hold", {
  d <- ald_data()
  al <- asym_laplace()
  frame_of <- function(fo, fam = al) {
    frm(fo, family = fam, data = d, dry_run = "frame")
  }
  open <- function(fr, prior = NULL, exact = FALSE) {
    expect_null(check_ald_laplace(fr$spec, fr, prior, exact = exact))
  }
  fg <- frame_of(bf(y ~ x + (1 | g)))
  # the guard's positive case on the same frame
  expect_error(check_ald_laplace(fg$spec, fg, NULL),
               class = "frmtmb_ald_no_laplace")
  # quadrature or importance: the marginal is an integral
  open(fg, exact = TRUE)
  # a prior on the scale can bound it
  open(fg, prior = set_prior("normal(0, 1)", class = "sd"))
  # a log link curves mu in the effect
  open(frame_of(bf(yp ~ x + (1 | g)), asym_laplace(link = "log")))
  # a random effect in sigma: the density is smooth there, and it
  # couples with mu's effects
  open(frame_of(bf(y ~ x, sigma ~ (1 | g))))
  open(frame_of(bf(y ~ x + (1 | g), sigma ~ (1 | g))))
  # a censored row's tail probability is curved
  open(frame_of(bf(yc | cens(cc) ~ x + (1 | g))))
  # an exact gp(): its objective was not seen to rise all the way
  # (TMB's inner search failed on its covariance at a large sd)
  open(frame_of(bf(y ~ x + gp(z))))
  # a Student-t block beside a gaussian one: its curvature depends on
  # its mode, which moves with the gaussian block's sd, so the
  # approximation can fall in that sd (by 10.4 at nu 2.5,
  # dev/reviews/wcheck2-short-student.R)
  d$h <- d$g
  open(frame_of(bf(y ~ x + (1 | g) + (1 | gr(h, dist = "student")))))
  # a smooth density: gaussian gives the effect its curvature
  open(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d,
           dry_run = "frame"))
  # no random effect at all
  open(frame_of(bf(y ~ x)))
})

test_that("the scaled structures profiled later are refused too", {
  # backlog item 5449: each structure's covariance is its sds times a
  # correlation, so scaling the sds by c scales it by c^2, and its
  # Laplace objective rose along that scale in every profile
  # (dev/alapse-ald-profile.R). Data with two rows per group and time,
  # so no effect is observation-level
  set.seed(20261116)
  d <- data.frame(g = factor(rep(1:20, each = 8)),
                  tt = factor(rep(rep(1:4, each = 2), 20)))
  d$x <- stats::rnorm(160)
  d$tn <- num_factor(c(0, 0.3, 0.9, 1.5)[as.integer(d$tt)])
  d$y <- 1 + 0.5 * d$x + stats::rnorm(160)
  d$cc <- ifelse(d$y > 2, "right", "none")
  d$yc <- pmin(d$y, 2)
  Q <- diag(20) * 2
  for (i in 1:19) Q[i, i + 1] <- Q[i + 1, i] <- -0.9
  dimnames(Q) <- list(levels(d$g), levels(d$g))
  terms <- list(hetar1 = "hetar1(tt + 0 | g)", ou = "ou(tn + 0 | g)",
                toep = "toep(tt + 0 | g)", homtoep = "homtoep(tt + 0 | g)",
                homcs = "homcs(tt + 0 | g)", exp = "exp(tn + 0 | g)",
                gau = "gau(tn + 0 | g)", mat = "mat(tn + 0 | g)",
                gr_prec = "(1 | gr(g, prec = Q))")
  frame_of <- function(lhs, tm, fam) {
    fo <- stats::as.formula(paste(lhs, "~ x +", tm))
    frm(bf(fo), family = fam, data = d, data2 = list(Q = Q),
        dry_run = "frame")
  }
  for (nm in names(terms)) {
    fr <- frame_of("y", terms[[nm]], asym_laplace())
    err <- tryCatch(check_ald_laplace(fr$spec, fr, NULL), error = identity)
    expect_s3_class(err, "frmtmb_ald_no_laplace")
    expect_match(conditionMessage(err), "frm_sample()", fixed = TRUE)
    # the guard's absent cases: a gaussian block on a smooth density,
    # and censored rows, whose tail probabilities are curved
    fg <- frame_of("y", terms[[nm]], gaussian())
    expect_null(check_ald_laplace(fg$spec, fg, NULL), label = nm)
    fc <- frame_of("yc | cens(cc)", terms[[nm]], asym_laplace())
    expect_null(check_ald_laplace(fc$spec, fc, NULL), label = nm)
  }
})
