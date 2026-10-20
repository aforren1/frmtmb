# A state's draws read that state's rows of every addition term. The
# hmm simulator handed the state-dependent family the whole term list,
# so the j-th row of a state read the trials() of the j-th row of the
# data: 26% of the draws of a binomial hmm exceeded N with frmtmb.latent
# 0.6.1 (dev/ndaterm-findings.md in the frmtmb repository; the
# reviewer's construction, dev/ndaterm-rev-probe3.R). The tolerance is
# a multiple of a Monte Carlo standard error measured from the draws.

skip_on_cran()

test_that("an hmm draws each state at its own rows' trials()", {
  set.seed(7)
  d <- data.frame(id = rep(1:20, each = 30), t = rep(1:30, 20),
                  N = sample(c(2, 50), 600, TRUE))
  st <- unlist(lapply(1:20, function(i) {
    s <- numeric(30)
    s[1] <- 1
    for (j in 2:30) {
      s[j] <- if (stats::runif(1) < 0.9) s[j - 1] else 3 - s[j - 1]
    }
    s
  }))
  d$y <- stats::rbinom(600, d$N, c(0.2, 0.8)[st])
  # the starting values take a logit of rows at 0 or N; that warning is
  # the start's, and only it is let through
  fit <- withCallingHandlers(
    frm(bf(y | trials(N) ~ 1),
        family = hmm(K = 2, binomial(), time = t, group = id), data = d),
    warning = function(w) {
      if (grepl("NaNs produced", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  s <- as.matrix(simulate(fit, nsim = 50, seed = 1))
  # 0.6.1: share above N 0.26, and a mean of 12.2 on the N = 2 rows
  # where fitted() says 0.91
  expect_identical(sum(sweep(s, 1, d$N, ">")), 0L)
  small <- s[d$N == 2, , drop = FALSE]
  se <- stats::sd(as.numeric(small)) / sqrt(length(small))
  expect_lt(abs(mean(small) - mean(fitted(fit)[d$N == 2, "Estimate"])),
            6 * se)
})
