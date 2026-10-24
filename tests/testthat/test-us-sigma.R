# us_sigma() and us_theta(): the unstructured covariance a us()
# parameter vector encodes, and its inverse, exported so that a family
# outside frmtmb can carry a free covariance matrix as an extra
# parameter. frmtmb.latent's mixture_mvn() is the first such family;
# these tests pin the contract it relies on.

test_that("us_sigma() follows the documented parameterization", {
  # sds 2 and 0.5; the one correlation entry 1 gives L = [1 0; 1 1],
  # whose normalized second row is (1, 1) / sqrt(2)
  S <- as.matrix(us_sigma(c(log(2), log(0.5), 1), 2))
  expect_equal(sqrt(diag(S)), c(2, 0.5))
  expect_equal(stats::cov2cor(S)[1, 2], 1 / sqrt(2))
  expect_equal(S, t(S))
  # zeros give the identity, at every dimension
  for (d in 1:4) {
    expect_equal(as.matrix(us_sigma(numeric(d + d * (d - 1) / 2), d)),
                 diag(d))
  }
})

test_that("the correlation entries fill L's lower triangle by column", {
  # d = 4 is the first dimension where the order is visible: by column
  # (frmtmb) and by row (glmmTMB) put the six entries in different
  # cells. The reference assigns each cell by name rather than through
  # lower.tri(), so it does not share the code under test.
  sds <- c(1.5, 0.7, 2, 1.1)
  r <- c(0.1, -0.2, 0.3, 0.4, -0.5, 0.6)
  build <- function(cells) {
    L <- diag(4)
    for (i in seq_along(cells)) {
      L[cells[[i]][1], cells[[i]][2]] <- r[i]
    }
    L <- L / sqrt(rowSums(L^2))
    (L %*% t(L)) * (sds %o% sds)
  }
  by_col <- list(c(2, 1), c(3, 1), c(4, 1), c(3, 2), c(4, 2), c(4, 3))
  by_row <- list(c(2, 1), c(3, 1), c(3, 2), c(4, 1), c(4, 2), c(4, 3))
  S <- as.matrix(us_sigma(c(log(sds), r), 4))
  expect_equal(S, build(by_col))
  # the guard on the guard: the two orders really differ at d = 4, so
  # the assertion above could fail
  expect_false(isTRUE(all.equal(build(by_col), build(by_row))))
})

test_that("us_sigma() is what a us() random-effect block reports", {
  # the same map VarCorr() reads for a correlated-slopes block, so a
  # family that calls it parameterizes its covariance as us() does
  set.seed(8)
  g <- factor(rep(1:40, each = 6))
  x <- rnorm(240)
  b <- matrix(rnorm(80), 40) %*% chol(matrix(c(1, 0.4, 0.4, 0.5), 2))
  y <- 1 + b[g, 1] + (0.5 + b[g, 2]) * x + rnorm(240, 0, 0.5)
  fit <- frm(bf(y ~ x + (x | g)), family = gaussian(),
             data = data.frame(y, x, g))
  th <- fit$estimates[["theta"]]
  # brms layout: the estimate is the first statistic of the cov array
  vc <- VarCorr(fit)[["g"]][["cov"]][, 1L, ]
  S <- as.matrix(us_sigma(th, 2, floored = TRUE))
  expect_equal(unname(S), unname(vc), tolerance = 1e-10)
})

test_that("us_sigma() runs on the tape, with its derivatives", {
  th <- c(0.3, -0.2, 0.5, 0.7, -0.4, 0.2)
  f <- function(p) as.vector(us_sigma(p, 3))
  tp <- RTMB::MakeTape(f, numeric(6))
  expect_equal(tp(th), f(th))
  # the AD Jacobian against central differences of the plain function,
  # held as a ratio to the Jacobian's own scale
  J <- tp$jacobian(th)
  h <- 1e-6
  Jfd <- vapply(seq_along(th), function(j) {
    e <- replace(numeric(length(th)), j, h)
    (f(th + e) - f(th - e)) / (2 * h)
  }, numeric(9))
  expect_equal(dim(J), dim(Jfd))
  expect_lt(max(abs(J - Jfd)) / max(abs(J)), 1e-6)
})

test_that("us_sigma() refuses a vector of the wrong length by name", {
  # five entries cannot be a 3 x 3 block, which takes six
  expect_error(us_sigma(numeric(5), 3), "takes 6 parameters, not 5")
  expect_error(us_sigma(numeric(4), 2), "takes 3 parameters, not 4")
  expect_error(us_sigma(numeric(6), 3), NA)
  expect_error(us_sigma(0, 0), "positive whole number")
  expect_error(us_sigma(numeric(3), 2.5), "positive whole number")
  expect_error(us_sigma(numeric(3), NA), "positive whole number")
  # a d that is not whole by a hair says so in the digits it prints
  expect_error(us_sigma(numeric(3), 2 + 1e-12), "2.000000000001",
               fixed = TRUE)
})

test_that("us_sigma(floored = TRUE) floors as the random-effect blocks do", {
  # exp(-400) underflows to 0. A 1 x 1 random-effect block floors the
  # variance (sd2_floored()), a larger one each standard deviation
  # (sd_floored()); us_sigma() does the same in each case.
  expect_identical(as.matrix(us_sigma(-400, 1))[1, 1], 0)
  expect_identical(as.matrix(us_sigma(-400, 1, floored = TRUE))[1, 1],
                   sd2_floored(-400))
  expect_gt(sd2_floored(-400), 0)
  S2 <- as.matrix(us_sigma(c(-400, -400, 0), 2, floored = TRUE))
  expect_identical(diag(S2), rep(sd_floored(-400)^2, 2))
  expect_gt(min(diag(S2)), 0)
})

test_that("us_theta() inverts us_sigma()", {
  set.seed(21)
  for (d in 1:5) {
    A <- matrix(stats::rnorm(d * d), d)
    S <- crossprod(A) + diag(d)
    th <- us_theta(S)
    expect_length(th, d + d * (d - 1) / 2)
    expect_null(names(th))
    expect_equal(as.matrix(us_sigma(th, d)), S, tolerance = 1e-12)
    th2 <- stats::rnorm(d + d * (d - 1) / 2)
    expect_equal(us_theta(as.matrix(us_sigma(th2, d))), th2,
                 tolerance = 1e-10)
  }
  # dimnames on the input do not reach the parameter vector
  S <- matrix(c(2, 0.5, 0.5, 1), 2, dimnames = list(c("a", "b"),
                                                     c("a", "b")))
  expect_null(names(us_theta(S)))
})

test_that("us_theta() refuses what is not a covariance matrix", {
  expect_error(us_theta(1:3), "square numeric matrix")
  expect_error(us_theta(matrix(1, 2, 3)), "square numeric matrix")
  expect_error(us_theta(matrix(c(1, NA, NA, 1), 2)),
               "square numeric matrix")
  expect_error(us_theta(matrix(c(1, 0.5, 0.2, 1), 2)), "symmetric")
  expect_error(us_theta(matrix(c(0, 0, 0, 1), 2)), "positive definite")
  # positive diagonal, not positive definite: refused by us_theta() in
  # its own words, not by the internal correlation step
  m <- tryCatch(us_theta(matrix(c(1, 2, 2, 1), 2)),
                error = conditionMessage)
  expect_match(m, "us_theta(): `S` must be a positive definite covariance",
               fixed = TRUE)
  expect_no_match(m, "requested correlation matrix", fixed = TRUE)
})
