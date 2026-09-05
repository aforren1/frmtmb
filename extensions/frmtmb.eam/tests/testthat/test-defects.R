## Regressions for the defects the feasibility study queued and the
## API notes in dev-findings.md recorded. Each test reproduces the
## reported behavior first and then pins what replaced it.

test_that("the density is -Inf below the non-decision time, not NaN", {
  # Queued defect: a decision time at or below zero came out NaN,
  # because log(u) of a negative u is NaN and NaN spreads through every
  # term after it. Zero is the right likelihood there and -Inf is the
  # right log likelihood; NaN is not an answer at all, and a mixture's
  # log-sum-exp cannot recover from one.
  for (t in c(-2, -0.1, -1e-9, 0)) {
    v <- ddm_lpdf_lower(t, 1, 1.4, 0.5)
    expect_false(is.nan(v))
    expect_identical(v, -Inf)
  }
  # both boundaries, through the data-carried reflection
  expect_identical(ddm_lpdf_both(-0.1, 1.2, 1.4, 0.4, 1), -Inf)
  expect_identical(ddm_lpdf_both(-0.1, 1.2, 1.4, 0.4, 0), -Inf)

  # and it is silent about it: the NaN used to arrive with a warning
  # from log(), which is noise in a fit that evaluates the density
  # thousands of times
  expect_silent(ddm_lpdf_lower(-0.1, 1, 1.4, 0.5))

  # vectorized, mixing rows either side of the edge, because that is
  # how it actually arrives
  got <- ddm_lpdf_both(c(-0.5, 0.4, 0, 1.2), 1, 1.4, 0.5, c(0, 1, 1, 0))
  expect_identical(is.finite(got), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("the floor on the normalized time is inert on the support", {
  # The fix must not perturb the density anywhere it was already right.
  # The floor is the smallest positive double, so adding it to any
  # normalized time a fit can reach is a no-op, and this asserts that
  # rather than trusting it.
  skip_if_not_installed("RWiener")
  gr <- expand.grid(t = c(0.01, 0.05, 0.4, 2.5, 8),
                    a = c(0.4, 1.4, 3.5), w = c(0.2, 0.5, 0.8))
  worst <- 0
  for (i in seq_len(nrow(gr))) {
    r <- log(RWiener::dwiener(gr$t[i] + 1e-9, gr$a[i], 1e-9, gr$w[i], 1.1,
                              resp = "lower"))
    # the reference itself underflows at the extremes of this grid; the
    # claim is about the rows where there is something to compare
    if (!is.finite(r)) next
    worst <- max(worst, abs(
      ddm_lpdf_lower(gr$t[i], 1.1, gr$a[i], gr$w[i]) - r) / abs(r))
  }
  expect_lt(worst, 1e-11)
})

test_that("max_ndt above min(rt) is refused alone and allowed in a mixture", {
  skip_if_not_installed("RWiener")
  # Queued defect: the guard is right for a bare Wiener model and wrong
  # inside a mixture, where the contaminant component is precisely what
  # covers the trials the diffusion cannot reach. Before the fix,
  # mixture(wiener(max_ndt = 0.4), ...) could not be fitted at all when
  # any response time fell below 0.4.
  set.seed(3)
  dat <- ddm_simulate(300, mu = 0.9, bs = 1.4, ndt = 0.30)
  k <- sample(300, 18)
  dat$rt[k] <- stats::runif(18, 0.12, 2.5)   # fast and slow guesses
  expect_lt(min(dat$rt), 0.4)

  # the refusal stands for the family on its own, and now says what to
  # do about it
  expect_error(frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
                   family = wiener(max_ndt = 0.4), data = dat),
               "above the smallest response time")
  expect_error(frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
                   family = wiener(max_ndt = 0.4), data = dat),
               "allow_unreachable")

  # and lifts inside a mixture
  fit <- frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
             family = mixture(wiener(max_ndt = 0.4,
                                     allow_unreachable = TRUE),
                              lognormal()),
             data = dat)
  expect_true(is.finite(as.numeric(logLik(fit))))
  # the point of the exercise: the non-decision time is free to sit
  # above the fastest response time, which a bare fit cannot do
  e <- unlist(fixef(fit))
  ndt <- 0.4 / (1 + exp(-e[["ndt1.(Intercept)"]]))
  expect_gt(ndt, min(dat$rt))

  # every row still has a likelihood, because the contaminant carries
  # the ones the Wiener component gives zero
  expect_true(all(is.finite(stats::residuals(fit, type = "response"))))
})

test_that("an unreachable row is a clean zero, not a NaN gradient", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")
  # Why the family's density floors the decision time at a margin below
  # each row's response time rather than returning -Inf: a true -Inf
  # exponentiates to zero in the mixture's sum, but differentiates to
  # NaN, and one NaN in a gradient stops the fit. A finite log density
  # of about -1/delta does both jobs.
  nd <- ddm_nodes(c("st"), c(sz = 1L, st = 21L))
  y <- c(0.5, 0.2)                     # the second row is unreachable
  lp <- function(p) {
    ddm_lpdf_var(y, p[1], 1.4, 0.5, p[2], 0, 0, p[3], c(1, 0), nd,
                 TRUE, 1e-9)
  }
  p <- c(1.0, 0.30, 0.05)
  v <- lp(p)
  expect_true(is.finite(v[1]))
  expect_lt(v[2], -1e6)                # a likelihood of zero, in double
  expect_equal(exp(v[2]), 0)
  tp <- RTMB::MakeTape(function(q) sum(lp(q)), p)
  g <- as.numeric(tp$jacobian(p))
  expect_true(all(is.finite(g)))
})

test_that("dec() carries the indicator the way brms spells it", {
  skip_if_not_installed("RWiener")
  # The ledger's finding 1: frmtmb's addition terms used to be a closed
  # set with no way in, so the indicator had to travel as vint() and the
  # user had to hand-code a factor to 0/1. frmtmb 0.49.0 added the
  # registry and this package registers `dec` when it loads.
  set.seed(21)
  dat <- ddm_simulate(400, mu = 0.9, bs = 1.4, ndt = 0.28)
  dat$resp <- factor(ifelse(dat$upper == 1, "upper", "lower"),
                     levels = c("lower", "upper"))
  # a factor, read on its levels with the second one as the upper
  # boundary, which is brms's rule
  f1 <- frm(bf(rt | dec(resp) ~ 1, bias = 0.5), family = wiener(),
            data = dat)
  f2 <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(),
            data = dat)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-10)

  # logical too, since FALSE sorts below TRUE
  dat$lg <- dat$upper == 1
  f3 <- frm(bf(rt | dec(lg) ~ 1, bias = 0.5), family = wiener(),
            data = dat)
  expect_equal(as.numeric(logLik(f3)), as.numeric(logLik(f2)),
               tolerance = 1e-10)

  # a column with one level is ambiguous and says so rather than
  # silently coding everything to the lower boundary
  dat$one <- factor("upper")
  expect_error(frm(bf(rt | dec(one) ~ 1, bias = 0.5), family = wiener(),
                   data = dat),
               "has two levels")
})

test_that("the missing decision indicator is refused, naming both spellings", {
  skip_if_not_installed("RWiener")
  set.seed(5)
  dat <- ddm_simulate(80, mu = 0.6, bs = 1.3, ndt = 0.2)
  # frmtmb_family(required_aterms =) declares the terms a density needs
  # ALL of. This family needs EITHER of two, so the refusal stays
  # written out; see ?wiener and dev-findings.md.
  expect_error(frm(bf(rt ~ 1, bias = 0.5), family = wiener(), data = dat),
               "decision indicator is missing")
  expect_error(frm(bf(rt ~ 1, bias = 0.5), family = wiener(), data = dat),
               "vint\\(upper\\)")
})

test_that("the family derives its link from the data through family_finalize", {
  skip_if_not_installed("RWiener")
  # The ledger's finding 3: the bound on the non-decision time is a
  # property of the response, and the family object is built before
  # frm() sees any data. This package used to have valid_y() write the
  # bound into an environment the link closures read at call time, which
  # worked only for as long as the undocumented slot order held.
  # frmtmb 0.49.0 added family_finalize() and the environment is gone.
  fam <- wiener()
  expect_true(is.function(fam[["family_finalize"]]))
  # the unfinalized family says its bound is missing rather than
  # returning a silent NA
  expect_error(fam$links[["ndt"]]$linkinv(0), "bound is not set yet")

  # and the finalized one carries a real link, on the family the fit
  # actually holds
  set.seed(6)
  dat <- ddm_simulate(200, mu = 0.9, bs = 1.3, ndt = 0.25)
  fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5), family = wiener(),
             data = dat)
  lk <- family(fit)$links[["ndt"]]
  expect_equal(lk$linkinv(0), min(dat$rt) / 2, tolerance = 1e-12)
  expect_no_error(lk$linkinv(3))
})

test_that("ddm_floor() holds its floor where a smooth maximum does not", {
  # The idiom this file's density depends on, asserted on its own so
  # that the two spellings cannot be confused again. A floor written as
  # the smooth maximum (x + lo + |x - lo|) / 2 looks equivalent and is
  # not: when `lo` is far below the rounding of `x`, both occurrences of
  # `lo` vanish and it collapses to the positive part, which is exactly
  # zero for a negative x. Zero is the value the floor exists to
  # prevent, because log(0) is -Inf and its gradient is NaN.
  lo <- 1e-300
  maxform <- function(x) 0.5 * (x + lo + abs(x - lo))
  for (x in c(0.3, 0, -1e-17, -1e-300)) {
    expect_gte(ddm_floor(x, lo), lo)
    expect_true(is.finite(log(ddm_floor(x, lo))))
  }
  # inert above the floor, to the last bit
  expect_identical(ddm_floor(0.3, lo), 0.3)
  expect_identical(ddm_floor(1, lo), 1)
  # and the counter-example that makes the helper worth having
  expect_identical(maxform(-1e-17), 0)
  expect_identical(ddm_floor(-1e-17, lo), lo)

  # the other floor in the density uses a denormal, where the same
  # argument is sharper still
  expect_identical(ddm_floor(0.5, ddm_u_floor), 0.5)
  expect_gt(ddm_floor(-0.05, ddm_u_floor), 0)
})

test_that("a range past the cap whose span rounds NEGATIVE stays finite", {
  skip_if_not_installed("RTMB")
  # `span` is the difference of two independently rounded evaluations of
  # the same cap, so on a row whose whole non-decision-time range has
  # moved past it, the surviving share is not reliably zero: it is zero
  # or a rounding residue of either sign. These three numbers produce
  # the negative on the platforms measured (Windows and Ubuntu), which
  # is the case a smooth maximum turns into log(0). The sign itself is
  # not asserted, because which way a given platform's arithmetic rounds
  # it is not this package's to promise; the helper's own
  # test above pins the negative input deterministically, and every
  # assertion below must hold whichever way this one rounds.
  y <- 0.152
  t0 <- 0.300
  st <- 0.098
  delta <- 1e-9 * y
  cap <- y - delta
  span <- ddm_smin(t0 + st / 2, cap) - ddm_smin(t0 - st / 2, cap)
  expect_lt(abs(span), 1e-12)
  expect_gt(t0 - st / 2, cap)          # the whole range really is past

  nd <- ddm_nodes("st", c(sz = 1L, st = 21L))
  v <- ddm_lpdf_var(y, 1.0, 1.4, 0.5, t0, 0, 0, st, 0, nd, TRUE, delta)
  expect_true(is.finite(v))
  expect_lt(v, -1e6)                   # a likelihood of zero, in double
  expect_equal(exp(v), 0)

  g <- as.numeric(RTMB::MakeTape(function(p) {
    ddm_lpdf_var(y, p[1], 1.4, 0.5, p[2], 0, 0, p[3], 0, nd, TRUE, delta)
  }, c(1.0, t0, st))$jacobian(c(1.0, t0, st)))
  expect_true(all(is.finite(g)))
})

test_that("a mixture survives a negative-span row end to end", {
  skip_if_not_installed("RWiener")
  # The end-to-end companion to the pinned point above. Coverage here
  # has to be arranged rather than hoped for: whether any one
  # unreachable row rounds its span negative depends on the exact
  # non-decision time the optimizer is standing on, so a single planted
  # row would land on the defect only by luck. Forty fast guesses
  # spread below the true non-decision time make it near-certain at
  # every evaluation, and the last assertion checks that the fit really
  # did end up in that regime rather than trusting it.
  set.seed(3)
  dat <- ddm_simulate(320, mu = 0.9, bs = 1.4, ndt = 0.34)
  k <- sample(320, 40)
  dat$rt[k] <- stats::runif(40, 0.15, 0.26)
  dat$upper[k] <- stats::rbinom(40, 1, 0.5)

  fit <- frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
             family = mixture(wiener(max_ndt = 0.4, variability = "st",
                                     allow_unreachable = TRUE),
                              lognormal()),
             data = dat)
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_true(all(is.finite(unlist(fixef(fit)))))

  e <- unlist(fixef(fit))
  ndt <- 0.4 / (1 + exp(-e[["ndt1.(Intercept)"]]))
  st <- 0.8 / (1 + exp(-e[["st1.(Intercept)"]]))
  # the non-decision time settled above the fastest response times,
  # which is the whole reason a mixture is here
  expect_gt(ndt, min(dat$rt))
  # and the fit ends where the defect lives: rows whose whole range is
  # past the cap. Whether any of their spans rounds NEGATIVE rather than
  # to exactly zero is a property of the platform's arithmetic, not of
  # the model: Windows produced negatives here and Ubuntu did not, so the
  # sign is not asserted: a past-the-cap span is zero in exact arithmetic
  # and a rounding residue of either sign in double. What is asserted is
  # the invariant the floor exists for: the residue is at rounding scale
  # and the floored share stays finite and negligible, whichever way it
  # rounded.
  cap <- dat$rt - 1e-9 * min(dat$rt)
  past <- (ndt - st / 2) > cap
  expect_gt(sum(past), 0)
  span <- ddm_smin(ndt + st / 2, cap) - ddm_smin(ndt - st / 2, cap)
  expect_true(all(abs(span[past]) < 1e-12))
  share <- ddm_floor(span[past] / st, ddm_share_floor)
  expect_true(all(is.finite(share)))
  expect_true(all(share <= 1e-12))
})
