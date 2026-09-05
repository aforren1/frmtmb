# frm_joint_cov() and frm_lp_basis(): the two exported seams of
# dev/spline-seam-proposal.md, Parts 1a and 1b. The property that
# matters is that predict(se.fit = TRUE) IS a consumer of the seam, so
# the two can never report different numbers.

sp_fit <- function() {
  set.seed(5051)
  dd <- data.frame(x = rnorm(150), g = factor(rep(1:15, each = 10)))
  dd$y <- rnorm(150, 1 + 2 * dd$x + rnorm(15, 0, 0.5)[dd$g], 0.4)
  list(d = dd, fit = frm(bf(y ~ s(x, k = 8) + (1 | g)), dd, gaussian()))
}

test_that("frm_joint_cov() spans beta, b and theta and labels every row", {
  skip_on_cran()
  o <- sp_fit()
  jc <- frm_joint_cov(o$fit)
  expect_true(all(c("beta", "b", "theta") %in% jc$names))
  expect_identical(dim(jc$V), c(length(jc$names), length(jc$names)))
  expect_identical(length(jc$labels), length(jc$names))
  expect_false(anyDuplicated(jc$labels) > 0)
  # the b rows the fixed-effect covariance cannot reach
  expect_gt(sum(jc$names == "b"), 0)
  expect_false(any(grepl("^b\\.", rownames(vcov(o$fit, full = TRUE)))))
  expect_error(frm_joint_cov(list()), "needs a model fitted by frm")
})

test_that("frm_lp_basis() reproduces predict(se.fit = TRUE) exactly", {
  skip_on_cran()
  o <- sp_fit()
  nd <- data.frame(x = seq(-2, 2, length.out = 12),
                   g = factor(1, levels = levels(o$d$g)))
  for (rf in list(NA, NULL)) {
    lb <- frm_lp_basis(o$fit, newdata = nd, re.form = rf)
    pr <- predict(o$fit, newdata = nd, re.form = rf, se.fit = TRUE)
    expect_equal(lb$eta, pr$fit)
    se <- sqrt(rowSums((lb$A %*% lb$V) * lb$A) + lb$extra_var)
    expect_equal(se, pr$se.fit, tolerance = 1e-12)
    expect_identical(length(lb$coef_pos), ncol(lb$A))
    expect_identical(length(lb$coef_names), ncol(lb$A))
    expect_identical(dim(lb$V), c(ncol(lb$A), ncol(lb$A)))
  }
  # the whole grid covariance, which predict() reduces to its diagonal
  lb <- frm_lp_basis(o$fit, newdata = nd, re.form = NA)
  Sigma <- lb$A %*% lb$V %*% t(lb$A)
  expect_equal(sqrt(diag(Sigma)),
               predict(o$fit, newdata = nd, re.form = NA,
                       se.fit = TRUE)$se.fit,
               tolerance = 1e-12)
  expect_true(isSymmetric(unname(Sigma), tol = 1e-10))
})

test_that("frm_lp_basis() works in sample and on a distributional dpar", {
  skip_on_cran()
  set.seed(5052)
  dd <- data.frame(x = rnorm(120))
  dd$y <- rnorm(120, 1 + 2 * dd$x, exp(-0.5 + 0.3 * dd$x))
  fit <- frm(bf(y ~ x, sigma ~ x), dd, gaussian())
  lb <- frm_lp_basis(fit, dpar = "sigma")
  pr <- predict(fit, dpar = "sigma", se.fit = TRUE)
  expect_equal(unname(lb$eta), unname(pr$fit))
  expect_equal(unname(sqrt(rowSums((lb$A %*% lb$V) * lb$A))),
               unname(pr$se.fit), tolerance = 1e-12)
  expect_true(all(grepl("^betad\\.", lb$coef_names)))
})

test_that("frm_lp_basis() gives a nonlinear body an exact Jacobian", {
  skip_on_cran()
  set.seed(5053)
  d3 <- data.frame(t = rep(seq(0, 3, length.out = 12), 20),
                   id = factor(rep(1:20, each = 12)))
  ri <- rnorm(20, 0, 0.3)[d3$id]
  d3$y <- (5 + ri) * (1 - exp(-1.2 * d3$t)) + rnorm(nrow(d3), 0, 0.2)
  fit <- frm(bf(y ~ ult * (1 - exp(-exp(lrc) * t)),
                ult ~ 1 + (1 | id), lrc ~ 1, nl = TRUE), d3, gaussian())
  nd <- data.frame(t = seq(0, 3, length.out = 7),
                   id = factor(1, levels = levels(d3$id)))
  lb <- frm_lp_basis(fit, newdata = nd, re.form = NA)
  expect_equal(lb$eta, unname(predict(fit, newdata = nd, re.form = NA)))

  # A is d eta / d coef; check every column against a central difference
  jc <- frm_joint_cov(fit)
  comp <- jc$names
  idx <- unlist(lapply(unique(comp), function(cp) seq_len(sum(comp == cp))))
  worst <- 0
  for (k in seq_along(lb$coef_pos)) {
    pos <- lb$coef_pos[k]
    h <- 1e-6
    fp <- fit; fp$estimates[[comp[pos]]][idx[pos]] <-
      fp$estimates[[comp[pos]]][idx[pos]] + h
    fm <- fit; fm$estimates[[comp[pos]]][idx[pos]] <-
      fm$estimates[[comp[pos]]][idx[pos]] - h
    fd <- (predict(fp, newdata = nd, re.form = NA) -
             predict(fm, newdata = nd, re.form = NA)) / (2 * h)
    worst <- max(worst, max(abs(fd - lb$A[, k])))
  }
  expect_lt(worst, 1e-6)

  # se.fit stays refused for a nonlinear predictor; this is the route
  expect_error(predict(fit, se.fit = TRUE), "nonlinear predictor")
  expect_error(frm_lp_basis(fit, newdata = nd, allow_new_levels = TRUE),
               "allow_new_levels")
})

test_that("frm_lp_basis() refuses arguments it cannot interpret", {
  skip_on_cran()
  o <- sp_fit()
  expect_error(frm_lp_basis(o$fit, newdata = 1:3), "must be a data frame")
  expect_error(frm_lp_basis(o$fit, re.form = "g"), "must be NULL")
  expect_error(frm_lp_basis(o$fit, dpar = "nope"), "unknown dpar")
  expect_error(frm_lp_basis(o$fit, resp = "nope"), "Unknown response")
})
