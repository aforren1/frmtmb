# frm_sample() on mm() with an |ID| key and with cov = A, which core builds
# as of this round. brms's names and the prior table are read here; the
# log density is the one core's L1 identities pin.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

mm_s_data <- function(seed = 70, n = 100) {
  set.seed(seed)
  lv <- letters[1:6]
  d <- data.frame(x = stats::rnorm(n),
                  g1 = factor(sample(lv, n, TRUE), levels = lv),
                  g2 = factor(sample(lv, n, TRUE), levels = lv))
  u <- stats::setNames(stats::rnorm(6, 0, 0.8), lv)
  d$y <- 0.5 * d$x + 0.5 * (u[as.character(d$g1)] + u[as.character(d$g2)]) +
    stats::rnorm(n, 0, 0.5)
  d
}

test_that("a keyed mm() block samples under brms's names", {
  skip_sampler()
  dd <- mm_s_data()
  s <- suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2))),
    data = dd, chains = 1, iter = 300, seed = 2, refresh = 0)))
  cn <- colnames(as.matrix(as_draws_matrix(s)))
  expect_true(all(c("sd_mmg1g2__Intercept", "sd_mmg1g2__sigma_Intercept",
                    "cor_mmg1g2__Intercept__sigma_Intercept") %in% cn))
  ep <- posterior_epred(s, newdata = dd[1:3, ])
  expect_identical(dim(ep)[2L], 3L)
})

test_that("mm(cov = ) samples, with its default priors", {
  skip_sampler()
  skip_if_not_installed("brms")
  dd <- mm_s_data()
  lv <- letters[1:6]
  A <- 0.5^abs(outer(1:6, 1:6, "-"))
  dimnames(A) <- list(lv, lv)
  tab <- as.data.frame(default_prior(bf(y ~ x + (1 | mm(g1, g2, cov = A))),
                                     data = dd, data2 = list(A = A),
                                     route = "sample"))
  bp <- as.data.frame(brms::default_prior(
    brms::bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = dd,
    data2 = list(A = A)))
  expect_identical(tab$prior[tab$class == "sd" & !nzchar(tab$group)],
                   bp$prior[bp$class == "sd" & !nzchar(bp$group)])
  s <- suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = dd, data2 = list(A = A),
    chains = 1, iter = 300, seed = 3, refresh = 0)))
  dm <- as.matrix(as_draws_matrix(s))
  expect_true(all(dm[, "sd_mmg1g2__Intercept"] > 0))
})
