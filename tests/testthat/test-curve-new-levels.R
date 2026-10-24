## `allow_new_levels`, and the population curve of a factor-smooth model.
##
## From frmtmb 0.65.0 `re_formula = NA` KEEPS an `s(t, g, bs = "fs")`
## term, as brms does, so the grid has to name a level of `g` and the
## curve is that level's. The documented way to read the POPULATION
## curve is an UNSEEN level with `allow_new_levels = TRUE`: mgcv's fs
## basis gives it a zero row. Before this argument existed the three
## curve functions could not pass the flag to core, and the vignette's
## population curve failed `R CMD build`.
##
## What is settled here is that the band at the unseen level describes
## the population smooth alone, measured against the band at a SEEN
## level, whose extra width is the between-subject part.

sp_fs_fixture <- local({
  memo <- NULL
  function() {
    if (!is.null(memo)) return(memo)
    set.seed(4)
    n_sub <- 12
    n_rep <- 4
    n_t <- 25
    peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
    sub <- rep(seq_len(n_sub), each = n_rep * n_t)
    d <- data.frame(subject = factor(sub),
                    t = rep(seq(0, 1, length.out = n_t),
                            times = n_sub * n_rep))
    h <- stats::rnorm(n_sub, 1, 0.12)
    sv <- stats::rnorm(n_sub, 0.5, 0.04)
    d$v <- peak(d$t, h[sub], sv[sub]) + stats::rnorm(nrow(d), 0, 0.06)
    # suppressWarnings: a gradient warning on some platforms is not what
    # this file measures, and the covariance check below would refuse a
    # fit whose covariance could not be read
    fit <- suppressWarnings(frmtmb::frm(
      frmtmb::bf(v ~ s(t, k = 10) + s(t, subject, bs = "fs", k = 5)),
      family = stats::gaussian(), data = d))
    tt <- seq(0.05, 0.95, length.out = 30)
    lev <- c(levels(d$subject), "population")
    memo <<- list(
      d = d, fit = fit,
      unseen = data.frame(t = tt,
                          subject = factor("population", levels = lev)),
      seen = data.frame(t = tt,
                        subject = factor("1", levels = levels(d$subject))))
    memo
  }
})

test_that("an unseen fs level gives the population prediction, checked", {
  o <- sp_fs_fixture()
  cv <- frm_curve(o$fit, newdata = o$unseen, re_formula = NA,
                  allow_new_levels = TRUE, simultaneous = FALSE)
  ref <- frmtmb::frm_linpred(o$fit, newdata = o$unseen, type = "link",
                             re_formula = NA, allow_new_levels = TRUE,
                             se.fit = TRUE)
  expect_equal(cv$.estimate, as.numeric(ref$fit))
  expect_equal(cv$.se, as.numeric(ref$se.fit))
  # the self-check ran, on the same rows, and passed: a disagreement
  # above `tol` would have refused rather than returned
  ck <- attr(cv, "check")
  expect_identical(ck$n_predict, 1L)
  expect_true(is.finite(ck$cov_rel_error))
  expect_true(attr(cv, "spec")$allow_new_levels)
})

test_that("the population band carries no between-subject variance", {
  o <- sp_fs_fixture()
  lb <- frmtmb::frm_lp_basis(o$fit, newdata = o$unseen, re_formula = NA,
                             allow_new_levels = TRUE)
  # an fs level core does not know adds no variance of its own; a bar
  # term's unseen level would, and the band would then be wider than the
  # population curve's
  expect_true(all(lb$extra_var == 0))
  # the band a reader wants for the population curve: the seam restricted
  # by hand to the intercept and the s(t) coefficients
  keep <- !grepl("subject", lb$coef_names, fixed = TRUE)
  expect_gt(sum(keep), 1L)
  expect_lt(sum(keep), length(keep))
  A <- as.matrix(lb$A)[, keep, drop = FALSE]
  se_pop <- sqrt(diag(A %*% lb$V[keep, keep, drop = FALSE] %*% t(A)))
  un <- frm_curve(o$fit, newdata = o$unseen, re_formula = NA,
                  allow_new_levels = TRUE, simultaneous = FALSE)
  sn <- frm_curve(o$fit, newdata = o$seen, re_formula = NA,
                  simultaneous = FALSE)
  # the scale is what a subject's own curve adds to the band, which this
  # run measures; the population band must be closer to the hand
  # restriction than that by many orders of magnitude
  gap_seen <- max(abs(sn$.se - se_pop))
  gap_unseen <- max(abs(un$.se - se_pop))
  expect_gt(gap_seen, 0)
  expect_lt(gap_unseen / gap_seen, 1e-8)
})

test_that("the default still refuses an unseen fs level, by name", {
  o <- sp_fs_fixture()
  expect_error(frm_curve(o$fit, newdata = o$unseen, re_formula = NA,
                         simultaneous = FALSE),
               "New levels in the factor-smooth term.*population",
               class = "frmtmb_new_levels")
  expect_error(frm_curve_deriv(o$fit, var = "t", newdata = o$unseen,
                               re_formula = NA, simultaneous = FALSE),
               "population", class = "frmtmb_new_levels")
  expect_error(frm_curve_feature(o$fit, var = "t", newdata = o$unseen,
                                 re_formula = NA),
               "population", class = "frmtmb_new_levels")
  # and the grouping column is still required, whatever the flag says
  expect_error(frm_curve(o$fit, newdata = o$unseen[, "t", drop = FALSE],
                         re_formula = NA, allow_new_levels = TRUE,
                         simultaneous = FALSE),
               "needs the grouping column `subject`")
})

test_that("the flag is validated like the package's other flags", {
  o <- sp_fs_fixture()
  for (bad in list(NA, "yes", c(TRUE, TRUE), 1)) {
    expect_error(frm_curve(o$fit, newdata = o$unseen,
                           allow_new_levels = bad),
                 "`allow_new_levels` must be TRUE or FALSE")
    expect_error(frm_curve_deriv(o$fit, var = "t", newdata = o$unseen,
                                 allow_new_levels = bad),
                 "`allow_new_levels` must be TRUE or FALSE")
    expect_error(frm_curve_feature(o$fit, var = "t", newdata = o$unseen,
                                   allow_new_levels = bad),
                 "`allow_new_levels` must be TRUE or FALSE")
  }
})

test_that("the derivative and the feature read the unseen level too", {
  o <- sp_fs_fixture()
  cv <- frm_curve(o$fit, newdata = o$unseen, re_formula = NA,
                  allow_new_levels = TRUE, simultaneous = FALSE)
  d1 <- frm_curve_deriv(o$fit, var = "t", newdata = o$unseen,
                        re_formula = NA, allow_new_levels = TRUE,
                        simultaneous = FALSE)
  expect_true(is.finite(attr(d1, "check")$cov_rel_error))
  # a stored curve carries the flag, so its derivative needs no restating
  d1c <- frm_curve_deriv(cv, var = "t", simultaneous = FALSE)
  expect_identical(d1c$.estimate, d1$.estimate)
  expect_identical(d1c$.se, d1$.se)
  # a curve stored at a SEEN level with the default, then a new grid at
  # the unseen level with TRUE passed here: the caller's TRUE is honored
  # rather than overridden by the stored FALSE
  cs <- frm_curve(o$fit, newdata = o$seen, re_formula = NA,
                  simultaneous = FALSE)
  d1s <- frm_curve_deriv(cs, var = "t", newdata = o$unseen,
                         allow_new_levels = TRUE, simultaneous = FALSE)
  expect_identical(d1s$.estimate, d1$.estimate)
  pk <- frm_curve_feature(o$fit, var = "t", type = "maximum",
                          newdata = o$unseen, re_formula = NA,
                          allow_new_levels = TRUE)
  expect_identical(nrow(pk), 1L)
  expect_true(is.finite(attr(pk, "check")$cov_rel_error))
  # the peak of the population curve, located on the curve frm_curve()
  # drew: the estimate there is the grid's maximum, to within the grid
  expect_lt(abs(pk$.estimate - cv$t[which.max(cv$.estimate)]),
            diff(cv$t[1:2]))
})

test_that("a difference across two unseen levels adds two draws", {
  # An unseen level of a BAR term carries its block's covariance and
  # loads no design column. Which rows share a draw is decided by the
  # level's label (frmtmb's new_key), and frm_lp_basis(extra_cov = TRUE)
  # on the stacked grid returns the cross term that decides it: two
  # DIFFERENT unseen levels are independent draws and add, the SAME
  # unseen level in both grids is one draw and cancels. Before that
  # seam this pair was refused, because sp_same_latent() could not tell
  # the two apart and would have cancelled two independent draws.
  o <- sp_fs_fixture()
  fit <- frmtmb::frm(frmtmb::bf(v ~ s(t, k = 10) + (1 | subject)),
                     family = stats::gaussian(), data = o$d)
  lev <- c(levels(o$d$subject), "A", "B")
  gA <- data.frame(t = o$unseen$t, subject = factor("A", levels = lev))
  gB <- transform(gA, subject = factor("B", levels = lev))
  la <- frmtmb::frm_lp_basis(fit, newdata = gA, re_formula = NULL,
                             allow_new_levels = TRUE)
  expect_true(all(la$extra_var > 0))
  dif <- frm_curve(fit, newdata = gA, contrast = gB, re_formula = NULL,
                   allow_new_levels = TRUE, simultaneous = FALSE)
  lb <- frmtmb::frm_lp_basis(fit, newdata = gB, re_formula = NULL,
                             allow_new_levels = TRUE)
  C <- as.matrix(la$A) - as.matrix(lb$A)
  v_coef <- rowSums((C %*% la$V) * C)
  expect_lt(max(abs(dif$.se^2 / (v_coef + la$extra_var + lb$extra_var) -
                      1)), 1e-10)
  # one unseen level on both sides, the curve moved along t: its draw
  # cancels and the difference is the smooth's alone
  gA2 <- transform(gA, t = rev(t))
  same <- frm_curve(fit, newdata = gA, contrast = gA2, re_formula = NULL,
                    allow_new_levels = TRUE, simultaneous = FALSE)
  l2 <- frmtmb::frm_lp_basis(fit, newdata = gA2, re_formula = NULL,
                             allow_new_levels = TRUE)
  C2 <- as.matrix(la$A) - as.matrix(l2$A)
  expect_lt(max(abs(same$.se^2 - rowSums((C2 %*% la$V) * C2))),
            1e-10 * max(la$extra_var))
  # the same pair with the grouping term dropped is the ordinary
  # population difference
  ok <- frm_curve(fit, newdata = gA, contrast = gB, re_formula = NA,
                  allow_new_levels = TRUE, simultaneous = FALSE)
  expect_true(all(ok$.estimate == 0))
})

test_that("an unseen bar-term level carries its draw into every route", {
  # One unseen level is ONE draw shared by every row, so the grid's
  # covariance is C V C' + Z S Z'. frm_lp_basis(extra_cov = TRUE)
  # returns the second part whole, and the simultaneous band, the
  # derivative and the feature read it. Before, they drew from C V C'
  # alone and were refused: the simultaneous critical value came out
  # BELOW qnorm(0.975) and covered 35 and 15 percent
  # (dev/reviews/2026-09-29-splinecurve.md, Finding A).
  o <- sp_fs_fixture()
  lev <- c(levels(o$d$subject), "new")
  gu <- data.frame(t = o$unseen$t, subject = factor("new", levels = lev))
  for (slope in c(FALSE, TRUE)) {
    form <- if (slope) v ~ s(t, k = 10) + (1 + t | subject) else
      v ~ s(t, k = 10) + (1 | subject)
    fit <- suppressWarnings(frmtmb::frm(frmtmb::bf(form),
                                        family = stats::gaussian(),
                                        data = o$d))
    lb <- frmtmb::frm_lp_basis(fit, newdata = gu, re_formula = NULL,
                               allow_new_levels = TRUE, extra_cov = TRUE)
    expect_true(all(lb$extra_var > 0))
    A <- as.matrix(lb$A)
    cv <- frm_curve(fit, newdata = gu, re_formula = NULL,
                    allow_new_levels = TRUE, nsim = 2000, seed = 1)
    expect_equal(cv$.se^2,
                 unname(diag(A %*% lb$V %*% t(A))) + lb$extra_var)
    Sig <- unname(A %*% lb$V %*% t(A)) + unname(lb$extra_cov)
    expect_lt(max(abs(attr(cv, "Sigma") - Sig)), 1e-12 * max(abs(Sig)))
    # the band now standardizes by the variance it draws from, so its
    # critical value is at least the pointwise one
    expect_gt(cv$.crit_sim[1L], cv$.crit[1L])
    # the derivative: a random intercept moves the curve and not its
    # slope, and a random slope adds its own variance to the slope's
    dv <- frm_curve_deriv(fit, var = "t", newdata = gu, re_formula = NULL,
                          allow_new_levels = TRUE, simultaneous = FALSE)
    dp <- frm_curve_deriv(fit, var = "t", newdata = gu, re_formula = NA,
                          allow_new_levels = TRUE, simultaneous = FALSE)
    if (!slope) {
      expect_lt(max(abs(dv$.se / dp$.se - 1)), 1e-6)
    } else {
      # the slope's variance, read off the new level's own variance
      # s11 + 2 t s12 + t^2 s22 at three values of t
      g3 <- data.frame(t = c(0, 1, 2), subject = gu$subject[1:3])
      e3 <- frmtmb::frm_lp_basis(fit, newdata = g3, re_formula = NULL,
                                 allow_new_levels = TRUE)$extra_var
      s22 <- (e3[3] - 2 * e3[2] + e3[1]) / 2
      expect_lt(max(abs(dv$.se^2 / (dp$.se^2 + s22) - 1)), 1e-4)
    }
    pk <- frm_curve_feature(fit, var = "t", type = "maximum",
                            newdata = gu, re_formula = NULL,
                            allow_new_levels = TRUE)
    expect_identical(nrow(pk), 1L)
    expect_true(is.finite(pk$.se))
    # core fills an ABSENT grouping column with NA under TRUE, which is
    # the same route
    cv2 <- frm_curve(fit, newdata = gu[, "t", drop = FALSE],
                     re_formula = NULL, allow_new_levels = TRUE,
                     nsim = 2000, seed = 1)
    expect_equal(cv2$.se, cv$.se)
  }
})

test_that("an explicit allow_new_levels = FALSE beats a stored TRUE", {
  # The stored value applies only when the caller says nothing. An OR of
  # the two could not express "refuse", and it once answered an explicit
  # FALSE on an unseen grid.
  o <- sp_fs_fixture()
  st <- frm_curve(o$fit, newdata = o$seen, re_formula = NA,
                  allow_new_levels = TRUE, simultaneous = FALSE)
  expect_error(frm_curve_deriv(st, var = "t", newdata = o$unseen,
                               allow_new_levels = FALSE,
                               simultaneous = FALSE),
               "population", class = "frmtmb_new_levels")
  expect_error(frm_curve_feature(st, var = "t", newdata = o$unseen,
                                 allow_new_levels = FALSE),
               "population", class = "frmtmb_new_levels")
  # said nothing: the stored TRUE is reused, and the answer is the one a
  # direct call with TRUE gives
  d_st <- frm_curve_deriv(st, var = "t", newdata = o$unseen,
                          simultaneous = FALSE)
  d_fit <- frm_curve_deriv(o$fit, var = "t", newdata = o$unseen,
                           re_formula = NA, allow_new_levels = TRUE,
                           simultaneous = FALSE)
  expect_identical(d_st$.estimate, d_fit$.estimate)
  expect_identical(d_st$.se, d_fit$.se)
})
