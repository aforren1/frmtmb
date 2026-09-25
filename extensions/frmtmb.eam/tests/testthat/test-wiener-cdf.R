## wiener(): the distribution of the response time, and cens() and
## trunc() on top of it (items 3.4, and punch round 1: the survival at
## large |v| a, and the known boundary on left and interval censoring).
##
## Two references, both stored as fixtures so that no Rmpfr is needed
## here: a 400-bit one on the grid test-density.R pins the density on
## (dev/phase3b-cdf-reference-v2.R), and a 700-bit one over v in
## -20..20, a in 0.3..6, w in 0.001..0.999 and u from 1e-3 to 10
## (dev/phase3b-cdf-reference-v3.R), which carries the survival and
## both defective distribution functions as logs. RWiener and WienR are
## the established implementations; each is checked at the accuracy it
## has.
##
## Every value a new behavior produces is computed inside tryCatch(),
## so on a build without it each assertion FAILS rather than the block
## stopping at its first error with the rest unrun.

cdf_try <- function(expr) tryCatch(expr, error = function(e) NULL)
ns_fun <- function(nm) cdf_try(get(nm, asNamespace("frmtmb.eam")))
rt_cdf <- function(t, v, a, w) {
  f <- ns_fun("ddm_rt_lcdf2")
  if (is.null(f)) {
    return(list(lF = rep(NA_real_, length(t)), lS = rep(NA_real_, length(t))))
  }
  f(t, v, a, w)
}
rt_cdf_b <- function(t, v, a, w, up) {
  f <- ns_fun("ddm_rt_lcdf_b")
  if (is.null(f)) return(rep(NA_real_, length(t)))
  f(t, v, a, w, up)
}

cdf_ref <- function() {
  x <- utils::read.csv(test_path("fixtures", "wiener-cdf-ref.csv"),
                       colClasses = c(F_ref = "character",
                                      S_ref = "character"))
  x$F_ref <- as.numeric(x$F_ref)
  x$S_ref <- as.numeric(x$S_ref)
  x
}
cdf_ref3 <- function() {
  x <- utils::read.csv(test_path("fixtures", "wiener-rtcdf-ref.csv"),
                       colClasses = c(lS = "character", lF = "character",
                                      lFl = "character", lFu = "character"))
  for (nm in c("lS", "lF", "lFl", "lFu", "v", "a", "w", "t")) {
    x[[nm]] <- as.numeric(x[[nm]])
  }
  x
}

test_that("F and S match the 400-bit reference on the density's grid", {
  ref <- cdf_ref()
  k <- ref$grid == "density"
  expect_identical(sum(k), 1125L)
  r <- rt_cdf(ref$t[k], ref$v[k], ref$a[k], ref$w[k])
  eF <- abs(exp(r$lF) - ref$F_ref[k]) / ref$F_ref[k]
  eS <- abs(exp(r$lS) - ref$S_ref[k]) / ref$S_ref[k]
  # RELATIVE, down to F = 1e-78 and S = 1e-50 on this grid
  expect_lt(max(eF), 1e-10)
  expect_lt(max(eS), 1e-10)
})

test_that("S, F and both defective functions hold out to |v| a = 120", {
  # Punch round 1, M1. The first version formed the small-time survival
  # as 1 - F, floored at 1e-300, and returned log S = -382.4 for a true
  # -73.97 at large |v| a. Errors are on the LOG, which is the relative
  # error of the quantity; the worst measured on the shipped build are
  # 8.0e-11 (S), 4.3e-14 (F) and 9.4e-13 (either boundary).
  ref <- cdf_ref3()
  expect_identical(nrow(ref), 3150L)
  expect_gt(max(abs(ref$v) * ref$a), 100)
  r <- rt_cdf(ref$t, ref$v, ref$a, ref$w)
  bl <- rt_cdf_b(ref$t, ref$v, ref$a, ref$w, 0)
  bu <- rt_cdf_b(ref$t, ref$v, ref$a, ref$w, 1)
  expect_lt(max(abs(r$lS - ref$lS)), 1e-9)
  expect_lt(max(abs(r$lF - ref$lF)), 1e-12)
  expect_lt(max(abs(bl - ref$lFl)), 1e-11)
  expect_lt(max(abs(bu - ref$lFu)), 1e-11)
  # a survival far below the smallest double is a finite, right log:
  # the reference goes to log S = -72063
  deep <- ref$lS < -800
  expect_gt(sum(deep), 50L)
  expect_lt(max(abs(r$lS[deep] - ref$lS[deep]) / abs(ref$lS[deep])), 1e-12)
})

test_that("the two defective functions sum to the marginal one", {
  ref <- cdf_ref3()
  bl <- rt_cdf_b(ref$t, ref$v, ref$a, ref$w, 0)
  bu <- rt_cdf_b(ref$t, ref$v, ref$a, ref$w, 1)
  r <- rt_cdf(ref$t, ref$v, ref$a, ref$w)
  m <- pmax(bl, bu)
  lsum <- m + log(exp(bl - m) + exp(bu - m))
  # measured 1.4e-14 over the 3150 rows
  expect_lt(max(abs(lsum - r$lF)), 1e-12)
})

test_that("RWiener::pwiener() agrees to 1e-10 on the density's grid", {
  skip_if_not_installed("RWiener")
  ref <- cdf_ref()
  k <- ref$grid == "density"
  rw <- mapply(function(t, v, a, w) {
    RWiener::pwiener(t + 1e-9, a, 1e-9, w, v, resp = "both")
  }, ref$t[k], ref$v[k], ref$a[k], ref$w[k])
  r <- rt_cdf(ref$t[k], ref$v[k], ref$a[k], ref$w[k])
  # ABSOLUTE, because that is the accuracy RWiener has: its own error
  # against the 400-bit truth is the larger part of this difference
  expect_lt(max(abs(exp(r$lF) - rw)), 1e-10)
})

test_that("WienR agrees per boundary, at the accuracy WienR has", {
  skip_if_not_installed("WienR")
  ref <- cdf_ref3()
  set.seed(36)
  k <- sample(which(ref$u >= 5e-3), 200)
  wl <- mapply(function(t, v, a, w) {
    WienR::pWDM(t, "lower", a = a, v = v, w = w, precision = 1e-12)$value
  }, ref$t[k], ref$v[k], ref$a[k], ref$w[k])
  bl <- rt_cdf_b(ref$t[k], ref$v[k], ref$a[k], ref$w[k], 0)
  # absolute, as WienR's precision argument is; its own error against
  # the 700-bit reference is 2e-13 over the whole grid
  expect_lt(max(abs(exp(bl) - wl)), 1e-11)
})

test_that("an interval's mass matches 1200-bit Rmpfr to 1e-10 in the log", {
  # Punch round 1: late in the distribution F(y2) and F(y) are both
  # within 1e-13 of the boundary probability and a difference of the two
  # came out negative. Punch round 2: the review's 511 Rmpfr intervals
  # found the first repair 25 units out in the log where the interval
  # held 1e-40 of the mass before it. The fixture is those 511 and 240
  # more to u1 = 1e-3 and |v| a = 80 (dev/phase3b-interval-fixture.R);
  # the worst error measured on the build is 4.9e-11.
  ref <- utils::read.csv(test_path("fixtures", "wiener-interval-ref.csv"))
  f <- ns_fun("ddm_rt_linterval_b")
  got <- if (is.null(f)) NA_real_ else
    cdf_try(f(ref$t1, ref$t2, ref$v, ref$a, ref$w, 0))
  if (is.null(got)) got <- NA_real_
  expect_lt(max(abs(got - ref$ref)), 1e-10)
})

test_that("the second derivatives are right at zero drift and far out", {
  # Punch round 1: the Laplace approximation differentiates twice. The
  # inner Hessian was -1e15 at a drift of exactly zero, every fit's
  # start, and NaN wherever a blend's tanh saturated past 700
  f <- ns_fun("ddm_rt_lcdf_b")
  g <- ns_fun("ddm_rt_linterval_b")
  ok <- FALSE
  if (!is.null(f) && !is.null(g)) {
    t <- c(0.05, 0.3, 1)
    tp <- RTMB::MakeTape(function(p) {
      sum(f(t, p[1], exp(p[2]), 0.5, c(0, 1, 1))) +
        sum(g(t, t + 0.1, p[1], exp(p[2]), 0.5, c(1, 0, 1)))
    }, c(0, 0))
    hf <- tp$jacfun()
    H <- hf$jacobian(c(0, log(1.5)))
    fd <- (tp$jacobian(c(1e-5, log(1.5))) - tp$jacobian(c(-1e-5, log(1.5)))) /
      2e-5
    ok <- abs(H[1, 1] - fd[1]) / max(1, abs(fd[1])) < 1e-6
    for (v in c(-4, -0.3, 0.3, 4)) for (la in log(c(0.5, 1.5, 4))) {
      ok <- ok && all(is.finite(hf$jacobian(c(v, la))))
    }
  }
  expect_true(ok)
})

test_that("a hierarchical fit with left and interval censoring converges", {
  # the design every fit of the first recovery arm failed on, smaller
  set.seed(28)
  d <- do.call(rbind, lapply(1:5, function(s) {
    x <- ddm_simulate(150, mu = 0.8 + stats::rnorm(1, 0, 0.3), bs = 1.4,
                      ndt = 0.25, bias = 0.5)
    x$s <- factor(s)
    x
  }))
  d$code <- 0L
  d$y2 <- d$rt
  fast <- d$rt < 0.45
  d$code[fast] <- -1L
  d$rt[fast] <- 0.45
  pick <- which(!fast)[seq(5, sum(!fast), by = 5)]
  lo <- floor(d$rt[pick] * 10) / 10
  d$code[pick] <- 2L
  d$y2[pick] <- lo + 0.1
  d$rt[pick] <- pmax(lo, 0.45)
  fit <- cdf_try(suppressWarnings(frm(
    bf(rt | dec(upper) + cens(code, y2) ~ 1 + (1 | s), bias = 0.5),
    family = wiener(), data = d)))
  ll <- if (is.null(fit)) NA_real_ else as.numeric(logLik(fit))
  expect_true(is.finite(ll))
})

test_that("the slots compute the censored rows only, and the same numbers", {
  # Punch round 1: computed on every row, the tape did not fit in memory
  # at 12,000 rows. Each slot now works on the rows frmtmb reads; on
  # those rows the values are the full computation's.
  set.seed(29)
  d <- ddm_simulate(200, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  d$code <- 0L
  d$code[d$rt > 1.1] <- 1L
  d$code[d$rt < 0.45] <- -1L
  d$y <- pmin(pmax(d$rt, 0.45), 1.1)
  fit <- cdf_try(frm(bf(y | dec(upper) + cens(code) ~ 1, bias = 0.5),
                     family = wiener(), data = d))
  ok <- FALSE
  if (!is.null(fit)) {
    fam <- stats::family(fit)
    av <- fit$frame$aterm_values$y
    dp <- list(mu = 0.7, bs = 1.3, ndt = 0.28, bias = 0.5)
    lS <- fam$lccdf(d$y, dp, av)
    Fv <- fam$lcdf(d$y, dp, av)
    r <- d$code == 1
    l <- d$code == -1
    full_S <- frmtmb.eam:::ddm_rt_lcdf2(d$y - 0.28, 0.7, 1.3, 0.5)$lS
    full_F <- exp(frmtmb.eam:::ddm_rt_lcdf_b(d$y - 0.28, 0.7, 1.3, 0.5,
                                             d$upper))
    ok <- sum(r) > 10 && sum(l) > 10 &&
      identical(lS[r], full_S[r]) && all(lS[!r] == 0) &&
      max(abs(Fv[l] - full_F[l]) / full_F[l]) < 1e-15 && all(Fv[!l] == 1)
  }
  expect_true(ok)
})

test_that("F integrates the density, over both boundaries", {
  set.seed(34)
  worst <- 0
  for (i in 1:12) {
    v <- runif(1, -3, 3); a <- runif(1, 0.5, 3); w <- runif(1, 0.2, 0.8)
    t1 <- runif(1, 0.01, 0.5); t2 <- t1 + runif(1, 0.05, 3)
    dens <- function(t) {
      exp(frmtmb.eam:::ddm_lpdf_both(t, v, a, w, 0)) +
        exp(frmtmb.eam:::ddm_lpdf_both(t, v, a, w, 1))
    }
    I <- stats::integrate(dens, t1, t2, rel.tol = 1e-12,
                          subdivisions = 2000L)$value
    Fd <- diff(exp(rt_cdf(c(t1, t2), v, a, w)$lF))
    worst <- max(worst, abs(Fd - I) / I)
  }
  expect_lt(worst, 1e-9)
})

test_that("at or below the non-decision time F is 0 and S is 1", {
  r <- rt_cdf(c(-0.5, 0, 1e-14), 1, 1.4, 0.5)
  expect_true(all(exp(r$lF) < 1e-250))
  expect_true(all(abs(r$lS) < 1e-250))
})

test_that("the gradients are finite where a fit can step", {
  # Punch round 1: a left-censored row at the fastest response puts its
  # decision time at zero when ndt reaches its bound, and the first
  # version returned NaN gradients there
  f <- ns_fun("ddm_rt_lcdf_b")
  g <- ns_fun("ddm_rt_lcdf2")
  q <- c(0.45, 0.6, 1.0)
  up <- c(0, 1, 1)
  ok <- FALSE
  if (!is.null(f) && !is.null(g)) {
    tp <- RTMB::MakeTape(function(p) {
      sum(f(q - p[3], p[1], exp(p[2]), 0.5, up)) +
        sum(g(q - p[3], p[1], exp(p[2]), 0.5)$lS)
    }, c(0, 0, 0))
    ok <- TRUE
    for (v in c(-60, -5, 0, 5, 60)) for (la in c(-4, 0, 6)) {
      for (t0 in c(0, 0.44, 0.45, 0.4500001, 2)) {
        ok <- ok && all(is.finite(tp$jacobian(c(v, la, t0))))
      }
    }
  }
  expect_true(ok)
})

# the likelihood written out by hand, at the fitted parameters, with
# RWiener as the density and the distribution functions: a right-
# censored row over both boundaries, a left- or interval-censored one
# at its own boundary
hand_ll <- function(d, y, code, y2, v, a, w, t0) {
  pb <- function(q, resp) {
    if (q <= t0) return(0)
    RWiener::pwiener(q, a, t0, w, v, resp = resp)
  }
  vapply(seq_len(nrow(d)), function(i) {
    b <- if (d$upper[i] == 1) "upper" else "lower"
    if (code[i] == 0) {
      log(RWiener::dwiener(y[i], a, t0, w, v, resp = b))
    } else if (code[i] == 1) {
      log(1 - pb(y[i], "both"))
    } else if (code[i] == -1) {
      log(pb(y[i], b))
    } else {
      log(pb(y2[i], b) - pb(y[i], b))
    }
  }, numeric(1))
}

fit_pars <- function(fit) {
  fe <- frmtmb::fixef_by_dpar(fit)
  lk <- stats::family(fit)$links
  list(v = fe$mu[[1]], a = exp(fe$bs[[1]]),
       t0 = lk$ndt$linkinv(fe$ndt[[1]]))
}

test_that("a right-censored fit scores what the hand-written likelihood says", {
  skip_if_not_installed("RWiener")
  set.seed(21)
  d <- ddm_simulate(400, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  cut <- 1.1
  d$code <- as.integer(d$rt > cut)
  d$y <- pmin(d$rt, cut)
  fit <- cdf_try(frm(bf(y | dec(upper) + cens(code) ~ 1, bias = 0.5),
                     family = wiener(), data = d))
  ll <- if (is.null(fit)) NA_real_ else as.numeric(logLik(fit))
  p <- if (is.null(fit)) list(v = NA, a = NA, t0 = NA) else fit_pars(fit)
  hand <- if (is.null(fit)) NA_real_ else
    sum(hand_ll(d, d$y, d$code, NULL, p$v, p$a, 0.5, p$t0))
  expect_gt(sum(d$code), 40)
  expect_lt(abs(ll - hand) / abs(hand), 1e-10)
})

test_that("a right-censored row's boundary is not read, exactly", {
  set.seed(21)
  d <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  d$code <- as.integer(d$rt > 1.1)
  d$y <- pmin(d$rt, 1.1)
  d2 <- d
  d2$upper[d2$code == 1] <- 1 - d2$upper[d2$code == 1]
  f <- bf(y | dec(upper) + cens(code) ~ 1, bias = 0.5)
  a <- cdf_try(frm(f, family = wiener(), data = d))
  b <- cdf_try(frm(f, family = wiener(), data = d2))
  la <- if (is.null(a)) NA_real_ else as.numeric(logLik(a))
  lb <- if (is.null(b)) -Inf else as.numeric(logLik(b))
  expect_identical(la, lb)
})

test_that("a left-censored row's boundary IS read", {
  # Decision (b), 2026-09-24: a left- or interval-censored trial reached
  # a boundary, and dec() says which
  set.seed(25)
  d <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  d$code <- as.integer(d$rt < 0.5) * -1L
  d$y <- pmax(d$rt, 0.5)
  d2 <- d
  d2$upper[d2$code == -1] <- 1 - d2$upper[d2$code == -1]
  f <- bf(y | dec(upper) + cens(code) ~ 1, bias = 0.5)
  a <- cdf_try(frm(f, family = wiener(), data = d))
  b <- cdf_try(frm(f, family = wiener(), data = d2))
  la <- if (is.null(a)) NA_real_ else as.numeric(logLik(a))
  lb <- if (is.null(b)) NA_real_ else as.numeric(logLik(b))
  expect_gt(sum(d$code == -1), 20)
  expect_true(is.finite(la) && is.finite(lb) && abs(la - lb) > 1)
})

test_that("all four censoring codes match the hand-written likelihood", {
  skip_if_not_installed("RWiener")
  set.seed(22)
  N <- 300L
  d <- ddm_simulate(N, mu = 0.6, bs = 1.5, ndt = 0.25, bias = 0.5)
  d <- d[d$rt < 2.5, ]
  N <- nrow(d)
  d$code <- 0L
  d$code[1:60] <- 2L
  d$code[61:100] <- 1L
  d$code[101:140] <- -1L
  d$y <- ifelse(d$code == 2L, d$rt - 0.04, d$rt)
  d$y2 <- d$rt + 0.04
  fit <- cdf_try(frm(bf(y | dec(upper) + cens(code, y2) ~ 1, bias = 0.5),
                     family = wiener(), data = d))
  ll <- if (is.null(fit)) NA_real_ else as.numeric(logLik(fit))
  p <- if (is.null(fit)) list(v = NA, a = NA, t0 = NA) else fit_pars(fit)
  hand <- if (is.null(fit)) NA_real_ else
    sum(hand_ll(d, d$y, d$code, d$y2, p$v, p$a, 0.5, p$t0))
  expect_lt(abs(ll - hand) / abs(hand), 1e-10)
})

test_that("trunc() divides by the window's mass over both boundaries", {
  skip_if_not_installed("RWiener")
  set.seed(22)
  d <- ddm_simulate(300, mu = 0.6, bs = 1.5, ndt = 0.25, bias = 0.5)
  d <- d[d$rt < 2.5, ]
  N <- nrow(d)
  ft <- cdf_try(frm(bf(rt | dec(upper) + trunc(ub = 2.5) ~ 1, bias = 0.5),
                    family = wiener(), data = d))
  llt <- if (is.null(ft)) NA_real_ else as.numeric(logLik(ft))
  pt <- if (is.null(ft)) list(v = NA, a = NA, t0 = NA) else fit_pars(ft)
  hand_t <- if (is.null(ft)) NA_real_ else
    sum(hand_ll(d, d$rt, rep(0L, N), NULL, pt$v, pt$a, 0.5, pt$t0)) -
      N * log(RWiener::pwiener(2.5, pt$a, pt$t0, 0.5, pt$v, resp = "both"))
  expect_lt(abs(llt - hand_t) / abs(hand_t), 1e-10)
})

test_that("left or interval censoring with trunc() is refused by name", {
  set.seed(26)
  d <- ddm_simulate(200, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  d <- d[d$rt < 2, ]
  d$code <- as.integer(d$rt < 0.5) * -1L
  d$y <- pmax(d$rt, 0.5)
  msg <- tryCatch({
    frm(bf(y | dec(upper) + cens(code) + trunc(ub = 2) ~ 1, bias = 0.5),
        family = wiener(), data = d)
    ""
  }, error = conditionMessage)
  expect_match(msg, "cannot be combined with trunc")
})

test_that("censoring reaches a per-group bound through ndt_group()", {
  set.seed(23)
  d <- do.call(rbind, lapply(1:4, function(s) {
    x <- ddm_simulate(150, mu = 0.8, bs = 1.4, ndt = 0.2 + 0.05 * s,
                      bias = 0.5)
    x$s <- s
    x
  }))
  d$code <- as.integer(d$rt > 1.2)
  d$y <- pmin(d$rt, 1.2)
  fit <- cdf_try(frm(bf(y | dec(upper) + cens(code) + ndt_group(s) ~ 1,
                        ndt ~ 1, bias = 0.5),
                     family = wiener(), data = d))
  hand <- NA_real_
  nb <- NA_integer_
  if (!is.null(fit)) {
    t0 <- ndt_time(fit)
    fe <- frmtmb::fixef_by_dpar(fit)
    r <- d$code == 1
    lS <- frmtmb.eam:::ddm_rt_lcdf2(d$y[r] - t0[r], fe$mu[[1]],
                                    exp(fe$bs[[1]]), 0.5)$lS
    lf <- frmtmb.eam:::ddm_lpdf_both(d$y[!r] - t0[!r], fe$mu[[1]],
                                     exp(fe$bs[[1]]), 0.5, d$upper[!r])
    hand <- sum(lS) + sum(lf)
    nb <- length(unique(round(t0, 12)))
  }
  # the per-group bounds are really different, or this would be the
  # scalar case under another name
  expect_identical(nb, 4L)
  ll <- if (is.null(fit)) NA_real_ else as.numeric(logLik(fit))
  expect_lt(abs(ll - hand) / abs(hand), 1e-12)
})

test_that("variability with cens() is refused by name", {
  set.seed(24)
  d <- ddm_simulate(100, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  d$code <- as.integer(d$rt > 1.1)
  msg <- tryCatch({
    frm(bf(rt | dec(upper) + cens(code) ~ 1, bias = 0.5),
        family = wiener(variability = "sv"), data = d)
    ""
  }, error = conditionMessage)
  expect_match(msg, "plain model only")
})

test_that("lba() and gddm() refuse cens() and trunc() by name", {
  set.seed(27)
  dl <- lba_simulate(200, v = c(2.4, 1.6), A = 0.5, k = 0.4, ndt = 0.2)
  dl$code <- as.integer(dl$rt > 1.2)
  refusal <- function(expr) {
    tryCatch({
      expr
      ""
    }, error = conditionMessage)
  }
  not_built <- "censoring with cens\\(\\) or trunc\\(\\) is not built"
  m1 <- refusal(frm(bf(rt | vint(choice) + cens(code) ~ 1),
                    family = lba(2), data = dl))
  expect_match(m1, paste0("^lba\\(\\): ", not_built))
  m2 <- refusal(frm(bf(rt | vint(choice) + trunc(ub = 5) ~ 1),
                    family = lba(2), data = dl))
  expect_match(m2, paste0("^lba\\(\\): ", not_built))
  set.seed(28)
  dg <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                      control = gddm_control(t_max = 2))
  dg$cond <- 1L
  dg$code <- as.integer(dg$rt > 1)
  small <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
  m3 <- refusal(frm(bf(rt | vint(upper, cond) + cens(code) ~ 1,
                       bias = 0.5),
                    family = gddm(control = small), data = dg))
  expect_match(m3, paste0("^gddm\\(\\): ", not_built))
  m4 <- refusal(frm(bf(rt | vint(upper, cond) + trunc(ub = 2) ~ 1,
                       bias = 0.5),
                    family = gddm(control = small), data = dg))
  expect_match(m4, paste0("^gddm\\(\\): ", not_built))
  m5 <- refusal(lba(2, contaminant = TRUE))
  expect_match(m5, "^lba\\(\\): contaminant = TRUE is not built")
  m6 <- refusal(rdm(2, contaminant = TRUE))
  expect_match(m6, "^rdm\\(\\): contaminant = TRUE is not built")
  # punch round 2: NA is not TRUE, and was reported as if it were
  m7 <- refusal(lba(2, contaminant = NA))
  expect_match(m7, "^lba\\(\\): `contaminant` must be FALSE")
})

test_that("the compatibility table says cens() and trunc() work", {
  expect_identical(frmtmb::frm_compat("wiener", "cens()")$status, "works")
  expect_identical(frmtmb::frm_compat("wiener", "trunc()")$status, "works")
})
