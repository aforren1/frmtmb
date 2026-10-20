# An end category of an ordinal family at tied thresholds, on doubles.
# frmtmb 0.71.0's cumulative() gave NaN there where the category has a
# positive probability: the masked interior pair of an end row was
# log(0), and 0 * -Inf is NaN (found by lane bin2port, log_lik:462 of
# brms's suite). The tape folded the masked product away, so fits were
# not affected; log_lik() of draws evaluates on doubles and was.
# dev/catmulti-tie-identity.R shows the taped value and gradient are
# bit for bit 0.71.0's at distinct thresholds.

ot_exact <- function(y, mu, tau, cdf) {
  vapply(seq_along(mu), function(s) {
    F <- c(0, cdf(tau - mu[s]), 1)
    log(F[y[s] + 1L] - F[y[s]])
  }, 0)
}

ot_check <- function(got, want) {
  fin <- is.finite(want)
  # zero-width categories are -Inf on both sides
  expect_identical(got[!fin], want[!fin])
  expect_true(all(is.finite(got[fin])))
  if (!any(fin)) return(invisible())
  expect_lt(max(abs(got[fin] - want[fin]) / abs(want[fin])), 1e-12)
}

test_that("cumulative() end categories are finite at tied thresholds", {
  mu <- c(-0.7, 0, 1.3)
  ties <- list(c(0, 0, 0), c(-1, -1, 1), c(-1, 1, 1))
  for (lnk in c("logit", "probit", "cloglog")) {
    lp <- cumulative(link = lnk)[["ord_lpdf_make"]](identity)
    cdf <- switch(lnk, logit = stats::plogis, probit = stats::pnorm,
                  cloglog = function(x) 1 - exp(-exp(x)))
    for (tau in ties) {
      for (y in 1:4) {
        got <- as.numeric(lp(rep(y, 3), list(mu = mu, disc = 1), list(),
                             list(tau_raw = tau)))
        ot_check(got, ot_exact(rep(y, 3), mu, tau, cdf))
      }
    }
  }
})

test_that("the tape and doubles agree at tied thresholds", {
  lp <- cumulative()[["ord_lpdf_make"]](identity)
  f <- function(p) {
    "c" <- RTMB::ADoverload("c")
    tau <- c(p[1], p[1], p[2])
    sum(lp(c(1, 1, 4), list(mu = c(-0.7, 0, 1.3), disc = 1), list(),
           list(tau_raw = tau)))
  }
  tp <- RTMB::MakeTape(f, c(0, 1))
  want <- sum(ot_exact(c(1, 1, 4), c(-0.7, 0, 1.3), c(0, 0, 1),
                       stats::plogis))
  expect_lt(abs(f(c(0, 1)) - want) / abs(want), 1e-12)
  expect_lt(abs(tp(c(0, 1)) - want) / abs(want), 1e-12)
})

test_that("hurdle_cumulative() and the sequential families stay finite", {
  mu <- c(-0.7, 0, 1.3)
  tau <- c(0, 0, 0)
  lph <- hurdle_cumulative()[["ord_lpdf_make"]](identity)
  for (y in c(1, 4)) {
    v <- as.numeric(lph(rep(y, 3), list(mu = mu, hu = 0.2, disc = 1),
                        list(), list(tau_raw = tau)))
    ot_check(v, log(0.8) + ot_exact(rep(y, 3), mu, tau, stats::plogis))
  }
  # a zero row reads category 1 in its masked ordinal term, so it was
  # NaN at a tie too (found by lane bin2port's reviewer,
  # dev/catmulti-bin2port-rev-cumfix.R); it is log(hu)
  v0 <- as.numeric(lph(rep(0, 3), list(mu = mu, hu = 0.2, disc = 1),
                       list(), list(tau_raw = tau)))
  ot_check(v0, rep(log(0.2), 3))
  # sratio, cratio and acat form no difference of two distribution
  # values, so a tie leaves every category finite (pinned, not changed)
  for (fam in list(sratio(), cratio(), acat())) {
    lp <- fam[["ord_lpdf_make"]](identity)
    for (y in 1:4) {
      v <- as.numeric(lp(rep(y, 3), list(mu = mu, disc = 1), list(),
                         list(tau_raw = tau)))
      expect_true(all(is.finite(v)), label = fam[["family"]])
    }
  }
})
