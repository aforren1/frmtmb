# Where the brms-parity lanes of 2026-09-25 meet. Each lane tested its
# feature alone; these blocks cover the combinations the merge created,
# including two defects that neither lane could have seen: the grouped-
# threshold category count read the FIRST response's family, and a
# cs() term sized its coefficients from a threshold vector that a
# multivariate frame stores under the response's own name.

skip_on_cran()

parity_data <- function(n = 400, seed = 11) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(c("a", "b"), each = n / 2)),
                  x = stats::rnorm(n))
  lat <- 0.8 * d$x + stats::rlogis(n)
  d$o <- ifelse(d$g == "a", findInterval(lat, c(-1, 0, 1.2)) + 1,
                findInterval(lat, c(-0.5, 0.8)) + 1)
  d$y <- 0.5 * d$x + stats::rnorm(n)
  d$s <- findInterval(0.6 * d$x + stats::rlogis(n), c(-1, 0.3, 1.5)) + 1
  d
}

# The two sides of each identity below are the same sum of the same
# densities, so the residual is rounding; the ratio is to the size of
# the log-likelihood itself.
ll_num <- function(f) as.numeric(stats::logLik(f))

test_that("grouped thresholds on the second response of an mv model", {
  d <- parity_data()
  f_mv <- suppressWarnings(frm(bf(y ~ x) + gaussian() +
                                 bf(o | thres(gr = g) ~ x) + cumulative(),
                               data = d))
  f_o <- suppressWarnings(frm(bf(o | thres(gr = g) ~ x) + cumulative(),
                              data = d))
  f_y <- frm(bf(y ~ x) + gaussian(), data = d)
  expect_lt(abs(ll_num(f_mv) - ll_num(f_o) - ll_num(f_y)) /
              abs(ll_num(f_mv)), sqrt(.Machine$double.eps))
  expect_true(all(c("o_Intercept[a,3]", "o_Intercept[b,2]") %in%
                    rownames(fixef(f_mv))))
  # four categories: the largest group's three thresholds plus one. The
  # merge first read the gaussian response's family here and counted
  # every threshold of both groups as one vector, six categories
  p_mv <- fitted(f_mv, resp = "o")
  p_o <- fitted(f_o)
  expect_identical(dim(p_mv)[3], dim(p_o)[3])
  expect_lt(max(abs(p_mv[, "Estimate", ] - p_o[, "Estimate", ])),
            1e3 * sqrt(.Machine$double.eps) * max(p_o[, "Estimate", ]))
})

test_that("cs() on an ordinal response inside an mv model", {
  d <- parity_data()
  f_mv <- frm(bf(s ~ cs(x)) + sratio() + bf(y ~ x) + gaussian(), data = d)
  f_s <- frm(bf(s ~ cs(x)) + sratio(), data = d)
  f_y <- frm(bf(y ~ x) + gaussian(), data = d)
  expect_lt(abs(ll_num(f_mv) - ll_num(f_s) - ll_num(f_y)) /
              abs(ll_num(f_mv)), sqrt(.Machine$double.eps))
  # one coefficient per threshold: three. The merge first sized them
  # from a threshold vector it could not find, zero coefficients
  lp <- f_mv$frame$linpreds[[linpred_key("s", "mu")]]
  expect_length(f_mv$estimates[[lp$cs[[1L]]$par]], 3L)
})

test_that("0 + Intercept and center = FALSE inside an mv model", {
  d <- parity_data()
  f1 <- frm(bf(s ~ x, center = FALSE) + cumulative() +
              bf(y ~ 0 + Intercept + x) + gaussian(), data = d)
  f0 <- frm(bf(s ~ x) + cumulative() + bf(y ~ x) + gaussian(), data = d)
  expect_lt(abs(ll_num(f1) - ll_num(f0)) / abs(ll_num(f0)),
            sqrt(.Machine$double.eps))
  pt <- prior_summary(f1)
  expect_false(any(pt$class == "Intercept" & pt$resp == "y"))
})
