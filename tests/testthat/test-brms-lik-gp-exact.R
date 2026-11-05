# The exact gp() against brms's Stan program, on data whose positions are
# NOT sorted. frmtmb sorts a term's distinct positions and brms keeps them
# in order of first appearance, so the translator must reorder the field
# before it reads z off a factor. Before that fix every row
# of this file read tens to thousands of nats. The suite never saw it:
# test-brms-likelihood.R's exact row has sorted x
# (dev/gpgap-findings.md).
#
# Stan compiles here, once per spelling, cached under FRMTMB_STAN_CACHE:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")

gpx_data <- function(freq = 0.7) {
  set.seed(20261102)
  n <- 80
  d <- data.frame(x = sample(0:12, n, TRUE), z = sample(0:12, n, TRUE),
                  f = factor(sample(c("a", "b"), n, TRUE)))
  e <- stats::rnorm(n, 0, 0.3)
  d$y <- 2 * sin(freq * d$x) + e
  d$ybz <- 2 * sin(freq * d$x) + 1.5 * cos(freq * d$z) + e
  d$ybf <- ifelse(d$f == "a", 2 * sin(freq * d$x), 1.5 * cos(freq * d$x)) + e
  d$y2 <- 2 * sin(freq * d$z) + stats::rnorm(n, 0, 0.3)
  # distinct on 0.001 steps, for a term that keeps one value per row
  d$xu <- round(d$x + stats::runif(n, -0.04, 0.04), 3)
  d
}


test_that("brms_gp_order() finds brms's positions in frmtmb's sort", {
  xb <- matrix(c(0.5, 0.1, 0.9, 0.3), ncol = 1L)
  sd <- list(Xgp_1 = xb, Jgp_1 = 1:4)
  io <- brms_gp_order(sd, "zgp_1")
  # frmtmb sorts: 0.1, 0.3, 0.5, 0.9. brms's first value, 0.5, is frmtmb's
  # third
  expect_identical(io, c(3L, 1L, 4L, 2L))
  expect_identical(sort(xb[, 1])[io], xb[, 1])
  # two coordinates sort on the first, then the second, as the frame does
  sd2 <- list(Xgp_1 = cbind(c(1, 0, 1), c(2, 9, 1)), Jgp_1 = 1:3)
  expect_identical(brms_gp_order(sd2, "zgp_1"), c(3L, 1L, 2L))
  # a by-level's and a response's data carry the z's own suffix
  sd3 <- list(Xgp_y_1_2 = xb, Jgp_y_1_2 = 1:4)
  expect_identical(brms_gp_order(sd3, "zgp_y_1_2"), io)
  # gr = FALSE has no Jgp and holds the rows in row order on both sides
  expect_identical(brms_gp_order(list(Xgp_1 = xb), "zgp_1"), 1:4)
})

test_that("exact gp() on unsorted positions", {
  skip_unless_brms_fit()
  dd <- gpx_data()
  specs <- list(
    list(bf(y ~ gp(x)), brms::bf(y ~ gp(x)), gaussian()),
    list(bf(ybz ~ gp(x, z, iso = FALSE)), brms::bf(ybz ~ gp(x, z, iso = FALSE)),
         gaussian()),
    list(bf(ybf ~ gp(x, by = f)), brms::bf(ybf ~ gp(x, by = f)), gaussian()),
    list(mvbf(bf(y ~ gp(x)), bf(y2 ~ gp(z)), rescor = FALSE),
         brms::mvbf(brms::bf(y ~ gp(x)), brms::bf(y2 ~ gp(z)), rescor = FALSE),
         NULL)
  )
  for (s in specs) {
    fit <- if (is.null(s[[3]])) {
      frm(s[[1]], data = dd)
    } else {
      frm(s[[1]], family = s[[3]], data = dd)
    }
    # the kernels are not diagonal: the order is what is being read
    bk <- fit$frame$re_blocks[[1L]]
    ev <- eigen(brms_block_cov(fit, bk), symmetric = TRUE,
                only.values = TRUE)$values
    expect_gt(max(ev) / min(ev), 1e3)
    expect_true(fit$opt$convergence == 0L)
    brms_lp_check(s[[2]], s[[3]], dd, fit, joint = TRUE)
  }
})

# The row below reads the order, on a kernel where the two nuggets' sizes differ by
# 0.026 nats and 2.3 in the gradient. brms's one line is read at frmtmb's relative
# size for that row alone; brms's own line is read in the other rows here and in
# test-brms-likelihood.R.
gpx_rel_nugget <- function(code) {
  new <- gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);",
              code, fixed = TRUE)
  testthat::expect_false(identical(new, code))
  new
}

test_that("exact gp() under gr = FALSE keeps row order on both sides", {
  # the case where brms_gp_order() must not sort
  skip_unless_brms_fit()
  dd <- gpx_data()
  fit <- frm(bf(y ~ gp(xu, gr = FALSE)), family = gaussian(), data = dd)
  bform <- brms::bf(y ~ gp(xu, gr = FALSE))
  expect_true(is.unsorted(dd$xu))
  brms_lp_check(bform, gaussian(), dd, fit, joint = TRUE,
                code_fix = gpx_rel_nugget)
})

test_that("exact gp() on a kernel near the identity is brms's own program", {
  # at a length scale short against the spacing no eigenvalue is near
  # either nugget, so the two factors agree and brms's program needs no
  # rewriting
  skip_unless_brms_fit()
  dd <- gpx_data(freq = 2.5)
  fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = dd)
  ev <- eigen(brms_block_cov(fit, fit$frame$re_blocks[[1L]]),
              symmetric = TRUE, only.values = TRUE)$values
  expect_lt(max(ev) / min(ev), 100)
  brms_lp_check(brms::bf(y ~ gp(x)), gaussian(), dd, fit, joint = TRUE)
})

test_that("the approximate gp() on unsorted positions needs no reorder", {
  # its map has no positions: z is a coefficient on a fixed basis, so the plain
  # term has to read brms's own program as it stands (the sorted row of
  # test-brms-likelihood.R could not tell). `by` is read in test-gp-by.R.
  skip_unless_brms_fit()
  dd <- gpx_data()
  stopifnot(is.unsorted(dd$x))
  fit <- frm(bf(y ~ gp(x, k = 10)), family = gaussian(), data = dd)
  expect_identical(fit$frame$re_blocks[[1L]]$covstruct, "hsgp")
  brms_lp_check(brms::bf(y ~ gp(x, k = 10)), gaussian(), dd, fit, joint = TRUE)
})
