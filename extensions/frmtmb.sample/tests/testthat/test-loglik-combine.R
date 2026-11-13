# log_lik(combine = FALSE) on a multivariate model without rescor: brms
# 2.23.0 returns a draws x observations x responses array, the layers
# named by response and summing to the combined matrix
# (dev/sampfront-brms-probe.R: dim 200 48 2, dimnames NULL NULL "y" "z").
# Through frmtmb.sample 0.28.0 it was refused (backlog 4584). The draws
# are built by hand around the ML estimates, so no sampler is needed.

cb_draws <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit, data = fit$data),
            class = "frmtmb_draws")
}

# the matrix without its attributes, as a layer of the array holds it
bare <- function(m) {
  attributes(m) <- list(dim = dim(m))
  m
}

cb_data <- function(seed = 4584L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(40), k = stats::rpois(40, 3))
  dd$y_a <- 1 + 0.5 * dd$x + stats::rnorm(40)
  dd$z <- -1 + 0.3 * dd$x + stats::rnorm(40)
  dd
}

test_that("combine = FALSE stacks one layer per response, as brms does", {
  dd <- cb_data()
  fit <- frm(mvbf(bf(y_a ~ x), bf(z ~ x), bf(k ~ x, family = poisson()),
                  rescor = FALSE), data = dd)
  ds <- cb_draws(fit)
  a <- log_lik(ds, combine = FALSE)
  expect_identical(dim(a), c(6L, 40L, 3L))
  # brms's response names, ya for y_a
  expect_identical(dimnames(a), list(NULL, NULL, c("ya", "z", "k")))
  for (j in seq_len(3L)) {
    r <- c("y_a", "z", "k")[j]
    expect_identical(a[, , j], bare(log_lik(ds, resp = r)),
                     info = r)
  }
  # the layers sum to the combined matrix, up to the order of the sum
  tot <- log_lik(ds)
  expect_lt(max(abs(apply(a, c(1, 2), sum) - tot)),
            64 * .Machine$double.eps * max(abs(tot)))
  # two responses named: their two layers
  a2 <- log_lik(ds, combine = FALSE, resp = c("z", "ya"))
  expect_identical(dimnames(a2)[[3L]], c("z", "ya"))
  expect_identical(a2[, , 2L], a[, , 1L])
  # the guard's condition absent: one response has nothing to combine,
  # and the draws argument is read as in the matrix
  expect_identical(log_lik(ds, combine = FALSE, resp = "z"),
                   log_lik(ds, resp = "z"))
  expect_identical(dim(log_lik(ds, combine = FALSE, ndraws = 3)),
                   c(3L, 40L, 3L))
  # and a subset of the fitted rows is cut in every layer
  expect_identical(log_lik(ds, combine = FALSE, newdata = dd[5:9, ]),
                   a[, 5:9, , drop = FALSE])
})

test_that("combine = FALSE on responses fitted on different rows refuses", {
  dd <- cb_data()
  dd$keep <- rep(c(TRUE, FALSE), 20L)
  fit <- frm(mvbf(bf(y_a | subset(keep) ~ x), bf(z ~ x), rescor = FALSE),
             data = dd)
  ds <- cb_draws(fit)
  expect_error(log_lik(ds, combine = FALSE),
               "fitted on different rows", fixed = TRUE)
  expect_identical(dim(log_lik(ds, combine = FALSE, resp = "y_a")),
                   c(6L, 20L))
})
