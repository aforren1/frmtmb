test_that("frame assembly: dimensions and index bookkeeping", {
  dd <- sim_pois_glmm()
  fr <- frm(bf(y ~ x + (x | g)) + poisson(), data = dd,
               dry_run = "frame")
  expect_s3_class(fr, "frmtmb_frame")
  expect_identical(fr$n_obs, nrow(dd))

  lp <- fr$linpreds[["y.mu"]]
  expect_identical(dim(lp$X), c(nrow(dd), 2L))
  expect_identical(dim(lp$Z), c(nrow(dd), 2L * nlevels(dd$g)))

  expect_length(fr$re_blocks, 1)
  bk <- fr$re_blocks[[1]]
  expect_identical(bk$covstruct, "us")
  expect_identical(bk$dim, 2L)
  expect_identical(bk$n_levels, nlevels(dd$g))
  expect_identical(bk$theta_idx, 1:3)
  expect_length(bk$b_idx, 2L * nlevels(dd$g))

  expect_named(fr$par_template, c("beta", "b", "theta"))
  expect_length(fr$par_template$beta, 2)
  expect_length(fr$par_template$b, 2L * nlevels(dd$g))
  expect_length(fr$par_template$theta, 3)
})

test_that("gaussian frame adds betad for sigma", {
  dd <- data.frame(y = rnorm(50), x = rnorm(50),
                   g = factor(rep(1:5, each = 10)))
  fr <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
               dry_run = "frame")
  expect_named(fr$par_template, c("beta", "betad", "b", "theta"))
  expect_identical(names(fr$par_template$betad), "sigma_(Intercept)")
  expect_length(fr$par_template$theta, 1)
})

test_that("NA rows are dropped consistently across design and aterms", {
  dd <- sim_pois_glmm(n_g = 5, n_per = 6)
  dd$w <- runif(nrow(dd), 0.5, 2)
  dd$x[3] <- NA
  dd$w[7] <- NA
  fr <- frm(bf(y | weights(w) ~ x + (1 | g)) + poisson(),
               data = dd, dry_run = "frame")
  expect_identical(fr$n_obs, nrow(dd) - 2L)
  expect_length(fr$aterm_values$y$weights, fr$n_obs)
  expect_length(fr$y$y, fr$n_obs)
})

test_that("factor binomial responses convert to 0/1", {
  dd <- data.frame(y = factor(sample(c("no", "yes"), 40, replace = TRUE)),
                   x = rnorm(40))
  fr <- frm(bf(y ~ x) + binomial(), data = dd, dry_run = "frame")
  expect_true(all(fr$y$y %in% c(0, 1)))
})

test_that("invalid responses are rejected at assembly time", {
  dd <- data.frame(y = c(-1, 2, 3), x = 1:3)
  expect_error(frm(bf(y ~ x) + poisson(), data = dd,
                      dry_run = "frame"),
               "non-negative integers")
})
