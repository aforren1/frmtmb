# posterior_linpred() with no dpar on a family with several location
# predictors stacks them, as brms 2.23.0 does: S x N x (K - 1) for a
# categorical or multinomial fit with a reference category, every
# margin unnamed (dev/catmulti-out/brms-fits.txt, from
# dev/catmulti-brms-fits.R). 0.19.0 returned the first predictor's
# S x N matrix.

lps_cat_data <- function(seed, n = 200) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x, -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$f <- factor(apply(P, 1, function(p) {
    sample(c("a", "b", "c"), 1, prob = p)
  }), levels = c("a", "b", "c"))
  Yc <- t(vapply(seq_len(n), function(i) {
    as.numeric(stats::rmultinom(1, 10, P[i, ]))
  }, numeric(3)))
  colnames(Yc) <- c("y1", "y2", "y3")
  dd$Yc <- Yc
  dd$N <- 10L
  dd
}

lps_sample <- function(fit) {
  allow_warnings(
    ds <- suppressMessages(frm_sample(fit, chains = 1, iter = 200,
                                      warmup = 100, refresh = 0, seed = 4)),
    c("Effective Samples Size", "R-hat"))
  ds
}

test_that("posterior_linpred() stacks a categorical fit's predictors", {
  skip_sampler()
  dd <- lps_cat_data(1)
  fit <- frm(bf(f ~ x), family = categorical(), data = dd)
  ds <- lps_sample(fit)
  pl <- posterior_linpred(ds, ndraws = 5)
  expect_identical(dim(pl), c(5L, 200L, 2L))
  expect_null(dimnames(pl))
  # each layer is that category's predictor, the dpar = call's matrix
  for (k in 1:2) {
    one <- posterior_linpred(ds, ndraws = 5, dpar = c("mub", "muc")[k])
    expect_identical(dim(one), c(5L, 200L))
    expect_identical(unname(pl[, , k]), unname(one))
  }
  # transform = TRUE asks brms for dpar "mu", which is not one of them
  expect_error(posterior_linpred(ds, transform = TRUE),
               paste("Invalid argument 'dpar'. Valid distributional",
                     "parameters are: 'mub', 'muc'"),
               fixed = TRUE)
})

test_that("multinomial draws stack, and one predictor stays a matrix", {
  skip_sampler()
  dd <- lps_cat_data(2)
  fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dd)
  ds <- lps_sample(fm)
  pl <- posterior_linpred(ds, ndraws = 4)
  expect_identical(dim(pl), c(4L, 200L, 2L))
  # brms's posterior_epred() is S x N x K, named by the columns
  ep <- posterior_epred(ds, ndraws = 4)
  expect_identical(dim(ep), c(4L, 200L, 3L))
  expect_identical(dimnames(ep)[[3L]], c("y1", "y2", "y3"))
  # the expected counts of a draw sum to the row's trials
  expect_lt(max(abs(apply(ep, 1:2, sum) - 10)) / 10, 1e-12)
  # a two-category categorical has one predictor: brms's S x N
  d2 <- dd
  d2$f <- factor(ifelse(dd$f == "a", "a", "b"))
  f2 <- suppressMessages(frm(bf(f ~ x), family = categorical(), data = d2))
  ds2 <- lps_sample(f2)
  expect_identical(dim(posterior_linpred(ds2, ndraws = 3)), c(3L, 200L))
})

test_that("transform = TRUE is refused on a two-category categorical too", {
  # one location dpar, mub, and no mu: the first build returned the
  # identity-link predictor as if transformed (review finding 3); brms
  # 2.23.0 refuses it, "Valid distributional parameters are: 'mub'"
  skip_sampler()
  dd <- lps_cat_data(3)
  dd$f <- factor(ifelse(dd$f == "a", "a", "b"))
  f2 <- suppressMessages(frm(bf(f ~ x), family = categorical(), data = dd))
  ds <- lps_sample(f2)
  expect_error(posterior_linpred(ds, transform = TRUE),
               "Valid distributional parameters are: 'mub'", fixed = TRUE)
  expect_identical(dim(posterior_linpred(ds, transform = TRUE, dpar = "mub",
                                         ndraws = 2)), c(2L, 200L))
})
