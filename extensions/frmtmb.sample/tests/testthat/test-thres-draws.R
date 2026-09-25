# Grouped ordinal thresholds, thres(gr = ), on sampled draws: each row
# on newdata has to reach the simulator with its group, or it has no
# thresholds to be drawn from.

test_that("posterior_predict() on newdata reads each row's group", {
  skip_on_cran()
  set.seed(11)
  n <- 240
  d <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE))
  u <- rlogis(n, 0.7 * d$x)
  tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6),
              c = c(-1.5, -0.4, 0.5, 1.6))
  d$y <- vapply(seq_len(n), function(i) 1L + sum(u[i] > tau[[d$g[i]]]),
                1L)
  fit <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative())
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 2)))
  nd <- data.frame(x = c(3, 3, 3), g = c("a", "b", "c"))
  top <- c(4L, 3L, 5L)
  pp <- posterior_predict(ds, newdata = nd, ndraws = 50)
  expect_true(all(t(pp) <= top))
  ep <- posterior_epred(ds, newdata = nd)
  expect_identical(dim(ep)[3L], 5L)
  expect_true(all(ep[, 2L, 4:5] == 0))
  expect_true(all(c("Intercept[b,2]", "Intercept[c,4]") %in%
                    rownames(fixef(ds))))
})
