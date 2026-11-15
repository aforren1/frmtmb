# The stationary-point escape moves one predictor off its declared point
# per start, with the others held. With two on their points at once one
# pass moved only one, and the fit kept the other there in silence: both
# alphas of mixture(skew_normal(), skew_normal()) started at 0 ended 7 to
# 9 below the fit from the default start (dev/quadseam2-findings.md,
# item 5222). The escape is now asked again from the kept optimum.

sn_draw <- function(n, mu, sigma, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  mu + sigma * (dl * abs(stats::rnorm(n)) +
                  sqrt(1 - dl^2) * stats::rnorm(n))
}

test_that("two skew normal components on alpha = 0 both escape", {
  set.seed(101)
  n <- 400
  z <- stats::rbinom(n, 1, 0.5)
  dd <- data.frame(y = ifelse(z == 1, sn_draw(n, -3, 1, 5),
                              sn_draw(n, 3, 1, -5)))
  mf <- mixture(skew_normal(), skew_normal())
  fo <- bf(y ~ 1)
  u <- frm(fo, family = mf, data = dd, dry_run = "objective")
  free <- frm(fo, family = mf, data = dd)
  bd <- u$estimates$betad
  al <- grep("^alpha", names(bd))
  expect_length(al, 2L)
  bd[al] <- 0
  stuck <- frm(fo, family = mf, data = dd, start = list(betad = bd))
  a <- stuck$estimates$betad[al]
  expect_true(all(abs(a) > 1))
  ll <- as.numeric(logLik(free))
  expect_lt(abs(as.numeric(logLik(stuck)) - ll), abs(ll) * 1e-8)
  esc <- stuck$opt[["stationary_escape"]]
  # two passes, four starts then two from the kept optimum; or one pass
  # of four when a start that moves one alpha carries the other off its
  # point on the way, which ends at the same fit
  expect_true(unname(esc[["starts"]]) %in% c(4, 6))
})
