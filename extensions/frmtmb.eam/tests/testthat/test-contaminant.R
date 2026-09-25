## wiener(contaminant = TRUE): a uniform contaminant over the observed
## response range, mixed in at a proportion `lambda` on a logit link
## (item 3.5 of dev/extension-gaps-plan.md).
##
## Every value a new behavior produces is computed inside tryCatch(),
## so on a build without it each assertion FAILS rather than the block
## stopping at its first error with the rest unrun.

ct_try <- function(expr) tryCatch(expr, error = function(e) NULL)

# the finalized family, as frm() would hand it to the objective
# The observed range is given explicitly: decision (a), 2026-09-24,
# makes it the default only under a declared deadline.
ct_family <- function(d, ..., contaminant_range = range(d$rt)) {
  fit <- ct_try(suppressWarnings(frm(
    bf(rt | dec(upper) ~ 1, bias = 0.5, lambda = 0.1),
    family = wiener(contaminant = TRUE, ...,
                    contaminant_range = contaminant_range), data = d)))
  if (is.null(fit)) NULL else stats::family(fit)
}

ct_data <- function(n = 400L, seed = 41L, lam = 0, clo = 0.35) {
  set.seed(seed)
  d <- ddm_simulate(n, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  hit <- stats::runif(n) < lam
  d$rt[hit] <- stats::runif(sum(hit), clo, 4)
  d$upper[hit] <- stats::rbinom(sum(hit), 1, 0.5)
  d
}

test_that("the mixture integrates to one over time and both boundaries", {
  d <- ct_data()
  fam <- ct_family(d)
  tot <- NA_real_
  if (!is.null(fam)) {
    cr <- fam$contaminant_range
    dp <- list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5, lambda = 0.2)
    dens <- function(t, bnd) {
      exp(fam$lpdf(t, dp, list(dec = rep(bnd, length(t)))))
    }
    tot <- 0
    for (up in 0:1) {
      # split at the range's edges, where the uniform steps
      for (lim in list(c(0.3 + 1e-12, cr[1]), cr, c(cr[2], 60))) {
        tot <- tot + stats::integrate(dens, lim[1], lim[2], bnd = up,
                                      rel.tol = 1e-12,
                                      subdivisions = 4000L)$value
      }
    }
  }
  expect_lt(abs(tot - 1), 1e-9)
})

test_that("at lambda = 0 the density is the plain family's", {
  d <- ct_data()
  fam <- ct_family(d)
  plain <- frmtmb.eam:::ddm_lpdf_both
  dp <- list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5, lambda = 0)
  y <- d$rt
  lp <- plain(y - 0.3, 0.8, 1.4, 0.5, d$upper)
  lc <- if (is.null(fam)) rep(NA_real_, length(y)) else
    fam$lpdf(y, dp, list(dec = d$upper))
  # the fallback keeps `up` a real selection on a build without the
  # contaminant, so the identity below cannot pass on NA against NA
  g <- if (is.null(fam)) -log(2 * diff(range(y))) else
    -log(2 * diff(fam$contaminant_range))
  up <- lp >= g
  # BITWISE where the diffusion's density is at least the contaminant's,
  # and to the rounding of one log-sum anchor elsewhere
  expect_gt(sum(up), 100L)
  expect_identical(lc[up], lp[up])
  expect_lt(max(abs(lc[!up] - lp[!up]) / abs(lp[!up])),
            64 * .Machine$double.eps)
})

test_that("a fit on clean data is the plain fit, and lambda collapses", {
  # THE COLLAPSE CASE. Maximum likelihood has no prior to hold the
  # mixing proportion off its edge, so on data with no contaminant it
  # runs down the logit until the likelihood stops changing. What must
  # be true there: the log-likelihood is the plain family's, the other
  # estimates are the plain fit's, and diagnose() says lambda is at the
  # end of its link rather than reporting an estimate.
  # Seed 1 is one that collapses: over seeds 1 to 40 of this design, 4
  # do and 36 stop inside the link (dev/phase3b-log/collapse-scan.txt)
  d <- ct_data(n = 600L, seed = 1L)
  f <- bf(rt | dec(upper) ~ 1, bias = 0.5)
  fp <- frm(f, family = wiener(), data = d)
  fc <- ct_try(suppressWarnings(frm(f, family = wiener(
    contaminant = TRUE, contaminant_range = range(d$rt)),
                                    data = d)))
  lc <- if (is.null(fc)) NA_real_ else as.numeric(logLik(fc))
  lp <- as.numeric(logLik(fp))
  expect_lt(abs(lc - lp) / abs(lp), 1e-6)
  eta <- if (is.null(fc)) NA_real_ else
    frmtmb::fixef_by_dpar(fc)$lambda[[1]]
  expect_lt(eta, -8)
  dg <- if (is.null(fc)) "" else
    paste(utils::capture.output(frmtmb::diagnose(fc)), collapse = "\n")
  expect_match(dg, "end of its link: lambda")
  mu <- if (is.null(fc)) NA_real_ else frmtmb::fixef_by_dpar(fc)$mu[[1]]
  expect_lt(abs(mu - frmtmb::fixef_by_dpar(fp)$mu[[1]]), 1e-3)
})

test_that("lambda recovers on contaminated data", {
  # the window starts at the non-decision time, below which no
  # diffusion row falls; a window from 0.35 left 11 rows outside it
  d <- ct_data(n = 3000L, seed = 43L, lam = 0.08, clo = 0.3)
  fc <- ct_try(frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
                   family = wiener(contaminant = TRUE,
                                   contaminant_range = c(0.3, 4)),
                   data = d))
  ci <- if (is.null(fc)) matrix(NA, 1, 3) else
    stats::confint(fc)["lambda_(Intercept)", , drop = FALSE]
  # a single draw: the Wald interval on the logit covers the truth
  expect_true(isTRUE(ci[1, 1] < stats::qlogis(0.08) &&
                       stats::qlogis(0.08) < ci[1, 2]))
})

test_that("a stated bound above the fastest response is allowed", {
  # a guess can be faster than the non-decision time, and under the
  # contaminant such a row has a likelihood, so max_ndt above min(rt)
  # is a model rather than a mistake
  d <- ct_data(n = 400L, seed = 44L)
  d$rt[1:5] <- c(0.12, 0.15, 0.18, 0.2, 0.22)
  fc <- ct_try(suppressWarnings(frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
                                    family = wiener(
                                      contaminant = TRUE, max_ndt = 0.5,
                                      contaminant_range = range(d$rt)),
                                    data = d)))
  ll <- if (is.null(fc)) NA_real_ else as.numeric(logLik(fc))
  expect_true(is.finite(ll))
  t0 <- if (is.null(fc)) NA_real_ else ndt_time(fc)[1]
  expect_gt(t0, 0.22)
  # without the contaminant the same bound is still refused
  msg <- tryCatch({
    frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
        family = wiener(max_ndt = 0.5), data = d)
    ""
  }, error = conditionMessage)
  expect_match(msg, "above the smallest")
})

test_that("censoring composes with the contaminant", {
  d <- ct_data(n = 400L, seed = 45L, lam = 0.05)
  d$code <- as.integer(d$rt > 1.5)
  d$y <- pmin(d$rt, 1.5)
  fc <- ct_try(frm(bf(y | dec(upper) + cens(code) ~ 1, bias = 0.5),
                   family = wiener(contaminant = TRUE,
                                   contaminant_range = range(d$y)),
                   data = d))
  hand <- NA_real_
  if (!is.null(fc)) {
    fam <- stats::family(fc)
    fe <- frmtmb::fixef_by_dpar(fc)
    lam <- stats::plogis(fe$lambda[[1]])
    a <- exp(fe$bs[[1]])
    t0 <- fam$links$ndt$linkinv(fe$ndt[[1]])
    cr <- fam$contaminant_range
    r <- d$code == 1
    f <- exp(frmtmb.eam:::ddm_lpdf_both(d$y[!r] - t0, fe$mu[[1]], a, 0.5,
                                        d$upper[!r]))
    g <- 0.5 / diff(cr)
    S <- exp(frmtmb.eam:::ddm_rt_lcdf2(d$y[r] - t0, fe$mu[[1]], a, 0.5)$lS)
    G <- (d$y[r] - cr[1]) / diff(cr)
    hand <- sum(log((1 - lam) * f + lam * g)) +
      sum(log((1 - lam) * S + lam * pmax(1 - G, 1e-300)))
  }
  ll <- if (is.null(fc)) NA_real_ else as.numeric(logLik(fc))
  expect_lt(abs(ll - hand) / abs(hand), 1e-12)
})

test_that("a left-censored row under the contaminant is scored at its edge", {
  # decision (b): the diffusion's defective function at the row's own
  # boundary, and half of the contaminant's G, the coin flip landing there
  d <- ct_data(n = 400L, seed = 49L, lam = 0.05)
  d$code <- -as.integer(d$rt < 0.55)
  d$y <- pmax(d$rt, 0.55)
  fc <- ct_try(frm(bf(y | dec(upper) + cens(code) ~ 1, bias = 0.5),
                   family = wiener(contaminant = TRUE,
                                   contaminant_range = c(0.3, 4)),
                   data = d))
  hand <- NA_real_
  if (!is.null(fc)) {
    fam <- stats::family(fc)
    fe <- frmtmb::fixef_by_dpar(fc)
    lam <- stats::plogis(fe$lambda[[1]])
    a <- exp(fe$bs[[1]])
    t0 <- fam$links$ndt$linkinv(fe$ndt[[1]])
    l <- d$code == -1
    f <- exp(frmtmb.eam:::ddm_lpdf_both(d$y[!l] - t0, fe$mu[[1]], a, 0.5,
                                        d$upper[!l]))
    g <- ifelse(d$y[!l] >= 0.3 & d$y[!l] <= 4, 0.5 / 3.7, 0)
    fb <- tryCatch(get("ddm_rt_lcdf_b", asNamespace("frmtmb.eam")),
                   error = function(e) function(...) NA_real_)
    Fb <- exp(fb(d$y[l] - t0, fe$mu[[1]], a, 0.5,
                  d$upper[l]))
    G <- (d$y[l] - 0.3) / 3.7
    hand <- sum(log((1 - lam) * f + lam * g)) +
      sum(log((1 - lam) * Fb + lam * 0.5 * G))
  }
  expect_gt(sum(d$code == -1), 10)
  ll <- if (is.null(fc)) NA_real_ else as.numeric(logLik(fc))
  expect_lt(abs(ll - hand) / abs(hand), 1e-12)
})

test_that("the conditional mean and the draw follow the contaminant", {
  skip_if_not_installed("RWiener")
  d <- ct_data(n = 400L, seed = 46L)
  fam <- ct_family(d)
  dp <- list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5, lambda = 0.3)
  n <- 5000L
  av <- list(dec = rep(1, n))
  m <- if (is.null(fam)) NA_real_ else fam$post$mean_fn(dp, list(dec = 1))
  set.seed(47)
  s <- if (is.null(fam)) NA_real_ else fam$sim(dp, av, n)
  # the simulated mean agrees with the closed form to its own Monte
  # Carlo error, four standard errors wide
  expect_lt(abs(mean(s) - m), 4 * stats::sd(s) / sqrt(n))
  # and neither is the plain family's
  plain <- frmtmb.eam:::ddm_mean_rt(dp, list(dec = 1))
  expect_gt(abs(m - plain), 0.05)
})

test_that("contaminant = TRUE refuses what it cannot mean", {
  m1 <- tryCatch({ wiener(contaminant = TRUE, allow_unreachable = TRUE); "" },
                 error = conditionMessage)
  expect_match(m1, "cannot be combined")
  m2 <- tryCatch({ wiener(contaminant = NA); "" }, error = conditionMessage)
  expect_match(m2, "must be TRUE or FALSE")
})

test_that("the default window needs a deadline", {
  # Decision (a), 2026-09-24. Without a deadline the observed range
  # is refused, classed; under trunc(ub =) it is the window from the
  # fastest response to the deadline
  d <- ct_data()
  f <- bf(rt | dec(upper) ~ 1, bias = 0.5)
  e <- tryCatch({ frm(f, family = wiener(contaminant = TRUE), data = d)
    NULL }, error = function(e) e)
  expect_s3_class(e, "frmtmb_eam_contaminant_range_error")
  expect_match(if (is.null(e)) "" else conditionMessage(e),
               "declares no response deadline")
  dd <- d[d$rt < 2.5, ]
  fit <- ct_try(suppressWarnings(frm(
    bf(rt | dec(upper) + trunc(ub = 2.5) ~ 1, bias = 0.5),
    family = wiener(contaminant = TRUE), data = dd)))
  cr <- if (is.null(fit)) c(NA, NA) else stats::family(fit)$contaminant_range
  expect_identical(cr, c(min(dd$rt), 2.5))
  dd$ub <- ifelse(seq_len(nrow(dd)) %% 2 == 0, 2.5, 3)
  m <- tryCatch({ frm(bf(rt | dec(upper) + trunc(ub = ub) ~ 1, bias = 0.5),
                      family = wiener(contaminant = TRUE), data = dd)
    "" }, error = conditionMessage)
  expect_match(m, "differs between rows")
})

test_that("a fitted range is kept", {
  d <- ct_data()
  fam <- ct_family(d)
  cr <- if (is.null(fam)) c(NA, NA) else fam$contaminant_range
  expect_identical(cr, range(d$rt))
})

test_that("a stated range replaces the observed one", {
  d <- ct_data()
  fam <- ct_family(d, contaminant_range = c(0, 6))
  cr <- if (is.null(fam)) c(NA, NA) else fam$contaminant_range
  expect_identical(cr, c(0, 6))
  # a response outside the stated range, scored directly, has the
  # diffusion's density alone, times the diffusion's share
  dp <- list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5, lambda = 0.2)
  fam2 <- ct_family(d[d$rt >= 0.4 & d$rt <= 2, ],
                    contaminant_range = c(0.4, 2))
  y <- c(0.35, 3)
  got <- if (is.null(fam2)) c(NA, NA) else
    fam2$lpdf(y, dp, list(dec = c(1, 0)))
  want <- log(0.8) + frmtmb.eam:::ddm_lpdf_both(y - 0.3, 0.8, 1.4, 0.5,
                                                c(1, 0))
  expect_equal(got, want, tolerance = 64 * .Machine$double.eps)
})

test_that("a response outside contaminant_range is refused, classed", {
  # Punch round 2: a response past the window's end was fitted without a
  # word, with contaminant density zero on that row
  d <- ct_data()
  e <- tryCatch({
    frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
        family = wiener(contaminant = TRUE,
                        contaminant_range = c(0, max(d$rt) - 0.01)),
        data = d)
    NULL
  }, error = function(e) e)
  expect_s3_class(e, "frmtmb_eam_contaminant_range_error")
  expect_match(if (is.null(e)) "" else conditionMessage(e),
               "outside contaminant_range")
})

test_that("contaminant_range is refused where it means nothing", {
  m1 <- tryCatch({ wiener(contaminant_range = c(0, 5)); "" },
                 error = conditionMessage)
  expect_match(m1, "there is none without contaminant = TRUE")
  m2 <- tryCatch({ wiener(contaminant = TRUE, contaminant_range = c(2, 1))
    "" }, error = conditionMessage)
  expect_match(m2, "two finite response")
})

test_that("the objective is smooth in lambda when rows sit below ndt", {
  # Punch round 1, B1. Under a stated max_ndt a row below the
  # non-decision time has a log density near -1.8e9, and the first
  # spelling formed log(lambda) + d with d near 1.8e9, which rounds
  # log(lambda) to 2.4e-7. The objective was a staircase in lambda. A
  # central difference at h = 1e-8 then disagrees with one at h = 1e-4
  # by orders of magnitude; on a smooth objective the two agree.
  d <- ct_data(n = 400L, seed = 48L)
  d$rt[1:5] <- c(0.12, 0.15, 0.18, 0.2, 0.22)
  fam <- ct_family(d, max_ndt = 0.5)
  obj <- function(eta) {
    if (is.null(fam)) return(NA_real_)
    dp <- list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5,
               lambda = stats::plogis(eta), .eta_lambda = eta)
    sum(fam$lpdf(d$rt, dp, list(dec = d$upper)))
  }
  fd <- function(h) (obj(-3 + h) - obj(-3 - h)) / (2 * h)
  wide <- fd(1e-4)
  narrow <- fd(1e-8)
  # rows below ndt are really there, or this would test nothing
  expect_gt(sum(d$rt < 0.3), 3L)
  # the narrow difference carries the rounding of a smooth sum of 400
  # terms of order 5, about 400 * 5 * eps / 1e-8; a staircase step of
  # 2.4e-7 per row would put it near 4 * 2.4e-7 / 2e-8 = 48 away
  expect_lt(abs(narrow - wide), 400 * 5 * .Machine$double.eps / 1e-8)
})
