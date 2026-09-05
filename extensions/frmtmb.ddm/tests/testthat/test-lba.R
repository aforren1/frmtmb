## The linear ballistic accumulator. rtdists is the canonical R
## implementation and is the external reference; where a claim cannot be
## made against it (the survival function in its far tail, the joint
## distribution of choice and time) the reference is a quadrature of the
## density or the generative process itself.

## The single-accumulator law, reached without going through a fit.
acc <- function(v, A = 0.6, k = 0.6, s = 1, posdrift = TRUE) {
  list(v = v, A = A, b = A + k, s = s, posdrift = posdrift)
}
law <- frmtmb.ddm:::lba_law
race <- frmtmb.ddm:::lba_race_lpdf

# ------------------------------------------------------------------ (a)
test_that("the defective density matches rtdists where rtdists is accurate", {
  skip_if_not_installed("rtdists")
  gr <- expand.grid(A = c(0.2, 0.5, 1, 2), k = c(0.1, 0.4, 1.5),
                    v = c(-0.5, 0.2, 1, 2.5, 4), s = c(0.5, 1, 1.7),
                    t = c(0.02, 0.1, 0.35, 0.8, 2, 6, 20))
  # h is the standardized drift that puts the accumulator at the top of
  # the start-point range. rtdists writes the normal difference as a
  # subtraction of two lower tails, which holds full precision only
  # while pnorm(h) is far from one. Adjudicated at 200 bits, rtdists is
  # already 1.35e-4 wrong at h = 7.6 where this family is 4.6e-15 wrong,
  # so parity is asserted where rtdists is accurate and the rest is
  # checked against Rmpfr in the next test.
  h <- (gr$k - gr$v * gr$t) / (gr$s * gr$t)
  for (pd in c(TRUE, FALSE)) {
    ref <- rtdists::dlba_norm(rt = gr$t, A = gr$A, b = gr$A + gr$k, t0 = 0,
                              mean_v = gr$v, sd_v = gr$s, posdrift = pd)
    mine <- exp(law$ldens(gr$t, list(v = gr$v, A = gr$A, b = gr$A + gr$k,
                                     s = gr$s, posdrift = pd)))
    keep <- ref > 1e-290 & h < 4.5
    expect_gt(sum(keep), 900)
    expect_lt(max(abs(mine[keep] - ref[keep]) / ref[keep]), 1e-11)
  }
})

test_that("in the fast tail we beat rtdists, against a 200-bit reference", {
  skip_if_not_installed("Rmpfr")
  gr <- expand.grid(A = c(0.2, 0.5, 1, 2), k = c(0.1, 0.4, 1.5),
                    v = c(-0.5, 0.2, 1, 2.5, 4), s = c(0.5, 1, 1.7),
                    t = c(0.02, 0.1, 0.35, 0.8, 2, 6, 20))
  hi <- (gr$A + gr$k - gr$v * gr$t) / (gr$s * gr$t)
  lo <- (gr$k - gr$v * gr$t) / (gr$s * gr$t)
  # Phi(hi) - Phi(lo) at 200 bits, through erfc so that neither tail
  # saturates in the reference either
  rt2 <- sqrt(Rmpfr::mpfr(2, 200))
  truth <- vapply(seq_along(hi), function(i) {
    a <- Rmpfr::mpfr(hi[i], 200); b <- Rmpfr::mpfr(lo[i], 200)
    as.numeric((Rmpfr::erfc(-a / rt2) - Rmpfr::erfc(-b / rt2)) / 2)
  }, numeric(1))
  mine <- frmtmb.ddm:::lba_phidiff(hi, lo)
  subtractive <- stats::pnorm(hi) - stats::pnorm(lo)
  pos <- truth > 0
  bulk <- pos & lo < 8
  tl <- pos & lo >= 8
  expect_gt(sum(tl), 40)
  # the form rtdists uses returns exactly zero on most of the tail
  expect_gt(sum(subtractive[tl] == 0), 40)
  # ours returns none there, and stays within a part per thousand
  expect_identical(sum(mine[tl] == 0), 0L)
  expect_lt(max(abs(mine[tl] - truth[tl]) / truth[tl]), 1e-3)
  # and in the bulk it is orders better than the subtractive form
  expect_lt(max(abs(mine[bulk] - truth[bulk]) / truth[bulk]), 1e-12)
  expect_gt(max(abs(subtractive[bulk] - truth[bulk]) / truth[bulk]), 1e-4)
  # and where the two implementations disagree inside the nominally
  # safe band, the 200-bit reference says it is rtdists that is wrong
  near <- pos & lo > 7 & lo < 8
  if (any(near)) {
    em <- abs(mine[near] - truth[near]) / truth[near]
    es <- abs(subtractive[near] - truth[near]) / truth[near]
    expect_lt(max(em), 1e-12)
    expect_gt(max(es), 1e-5)
  }
  # the caveat, pinned rather than hidden: a handful of rows are only
  # good to parts per million, so this is an improvement, not a proof
  rel <- abs(mine[pos] - truth[pos]) / truth[pos]
  expect_lt(sum(rel > 1e-6), 10L)
  expect_lt(max(rel), 1e-3)
})

test_that("the survival function matches rtdists wherever rtdists has digits", {
  skip_if_not_installed("rtdists")
  gr <- expand.grid(A = c(0.3, 1), k = c(0.2, 0.9), v = c(0.3, 1.5, 3),
                    s = c(0.7, 1), t = c(0.05, 0.3, 1, 4, 15))
  ref <- 1 - rtdists::plba_norm(rt = gr$t, A = gr$A, b = gr$A + gr$k, t0 = 0,
                                mean_v = gr$v, sd_v = gr$s)
  mine <- exp(law$lsurv(gr$t, list(v = gr$v, A = gr$A, b = gr$A + gr$k,
                                   s = gr$s, posdrift = TRUE)))
  keep <- ref > 1e-12
  expect_gt(sum(keep), 100)
  expect_lt(max(abs(mine[keep] - ref[keep]) / ref[keep]), 1e-9)
})

test_that("the survival function keeps its digits where 1 - plba_norm loses them", {
  skip_if_not_installed("rtdists")
  # a fast competitor that has almost certainly finished: rtdists's
  # distribution function rounds to exactly one, so 1 - F is exactly zero
  # and a loser's log-contribution would be -Inf on an ordinary row
  p <- acc(v = 10, A = 0.6, k = 0.6)
  tt <- c(1, 2, 5)
  expect_equal(1 - rtdists::plba_norm(tt, A = 0.6, b = 1.2, t0 = 0,
                                      mean_v = 10, sd_v = 1),
               rep(0, 3))
  mine <- exp(law$lsurv(tt, p))
  expect_true(all(mine > 0))
  # the reference here is a high-accuracy quadrature of the density
  quad <- vapply(tt, function(q) {
    stats::integrate(function(u) exp(law$ldens(u, p)), q, Inf,
                     rel.tol = 1e-12, subdivisions = 5000L)$value
  }, numeric(1))
  expect_lt(max(abs(mine - quad) / quad), 1e-9)
})

test_that("the race matches rtdists n1PDF for two, three and four accumulators", {
  skip_if_not_installed("rtdists")
  set.seed(7)
  A <- 0.7; k <- 0.6; s <- 1; b <- A + k
  for (nacc in c(2L, 3L, 4L)) {
    m <- 40L
    tt <- runif(m, 0.25, 2.5)
    V <- matrix(runif(m * nacc, -0.2, 3), ncol = nacc)
    win <- rep_len(seq_len(nacc), m)
    pars <- lapply(seq_len(nacc), function(j)
      list(v = V[, j], A = A, b = b, s = s, posdrift = TRUE))
    mine <- exp(race(tt, win, law, pars))
    ref <- vapply(seq_len(m), function(i) {
      ord <- c(win[i], setdiff(seq_len(nacc), win[i]))
      rtdists::n1PDF(tt[i], A = A, b = b, t0 = 0, mean_v = V[i, ord],
                     sd_v = rep(s, nacc), silent = TRUE)
    }, numeric(1))
    expect_lt(max(abs(mine - ref) / ref), 1e-11)
  }
})

test_that("the family's log-density is the race, through a fitted object", {
  skip_if_not_installed("rtdists")
  set.seed(21)
  dat <- lba_simulate(150, v = c(2.2, 1.4, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  fit <- frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = dat)
  ll <- as.numeric(logLik(fit))
  e <- fixef(fit)
  A <- exp(e$A[[1]]); k <- exp(e$k[[1]])
  ndt <- stats::family(fit)$links$ndt$linkinv(e$ndt[[1]])
  vv <- c(e$v1[[1]], e$v2[[1]], e$v3[[1]])
  ref <- vapply(seq_len(nrow(dat)), function(i) {
    ord <- c(dat$choice[i], setdiff(1:3, dat$choice[i]))
    log(rtdists::n1PDF(dat$rt[i] - ndt, A = A, b = A + k, t0 = 0,
                       mean_v = vv[ord], sd_v = rep(1, 3), silent = TRUE))
  }, numeric(1))
  expect_lt(abs(ll - sum(ref)), 1e-6)
})

# ------------------------------------------------------------------ (b)
test_that("the analytic joint distribution reproduces the generative process", {
  set.seed(99)
  A <- 0.6; k <- 0.6; ndt <- 0.15; vv <- c(2.2, 1.4, 0.7)
  N <- 200000L
  d <- lba_simulate(N, v = vv, A = A, k = k, ndt = ndt)
  pars <- lapply(vv, function(v) acc(v, A, k))
  edges <- c(0.15, 0.4, 0.6, 0.9, 1.4, 3, Inf)
  worst <- 0
  for (j in 1:3) {
    for (m in seq_len(length(edges) - 1L)) {
      lo <- edges[m]; hi <- min(edges[m + 1L], 80)
      emp <- mean(d$choice == j & d$rt > lo & d$rt <= hi)
      ana <- stats::integrate(function(q) {
        exp(race(q - ndt, rep(j, length(q)), law, pars))
      }, lo, hi, rel.tol = 1e-10, subdivisions = 3000L)$value
      se <- sqrt(max(emp * (1 - emp), 1e-12) / N)
      worst <- max(worst, abs(emp - ana) / se)
    }
  }
  # 18 cells; the analytic value is exact, so the only error is Monte
  # Carlo and the largest standardized gap should look like a normal max
  expect_lt(worst, 4.5)
})

test_that("choice probabilities sum to one under truncation and less without", {
  A <- 0.6; k <- 0.6; vv <- c(1.5, 0.4, -0.6)
  for (pd in c(TRUE, FALSE)) {
    pars <- lapply(vv, function(v) acc(v, A, k, posdrift = pd))
    tot <- sum(vapply(1:3, function(j) {
      stats::integrate(function(q) exp(race(q, rep(j, length(q)), law, pars)),
                       1e-8, Inf, rel.tol = 1e-10,
                       subdivisions = 5000L)$value
    }, numeric(1)))
    if (pd) {
      # every accumulator arrives eventually, so the race always ends
      expect_equal(tot, 1, tolerance = 1e-6)
    } else {
      # the missing mass is trials on which no accumulator would respond
      expect_lt(tot, 1 - 1e-3)
      expect_equal(tot, 1 - prod(stats::pnorm(-vv)), tolerance = 1e-6)
    }
  }
})

test_that("the untruncated survival adds the never-arriving mass exactly", {
  # S_untrunc(t) - S_trunc-numerator = P(drift <= 0) = Phi(-v/s). Spelled
  # 1 - Phi(v/s) it is a subtraction near one and loses everything in the
  # same place the survival was written out to protect: 7 percent wrong
  # at v/s = 8 and total loss at v/s = 10.
  for (v in c(-1, 0.5, 2, 5, 6, 8, 10, 12)) {
    p <- acc(v, posdrift = FALSE)
    pt <- acc(v, posdrift = TRUE)
    tt <- c(0.3, 1, 1.5, 2, 2.5, 3)
    su <- exp(law$lsurv(tt, p))
    # sf, the finite part, is the truncated survival times Phi(v/s)
    sf <- exp(law$lsurv(tt, pt)) * stats::pnorm(v)
    got <- su - sf
    want <- stats::pnorm(-v)
    expect_equal(got, rep(want, length(tt)), tolerance = 1e-9,
                 info = paste("v =", v))
  }
  # the spelling that fails: at v/s = 10 it is a total loss
  expect_identical(1 - stats::pnorm(10), 0)
  expect_gt(stats::pnorm(-10), 0)
})

test_that("the model is identified only up to scale, which sd_v fixes", {
  # multiplying A, b, every drift and the drift sd by one constant leaves
  # the distribution unchanged, which is why sd_v is not estimated
  tt <- c(0.2, 0.6, 1.5)
  base <- lapply(c(1.8, 1.0), function(v) acc(v, A = 0.5, k = 0.5, s = 1))
  cc <- 2.7
  scaled <- lapply(c(1.8, 1.0) * cc, function(v)
    acc(v, A = 0.5 * cc, k = 0.5 * cc, s = cc))
  expect_equal(race(tt, c(1, 2, 1), law, base),
               race(tt, c(1, 2, 1), law, scaled), tolerance = 1e-9)
})

# ------------------------------------------------------------------ (c)
test_that("the taped gradient matches numDeriv at several parameter points", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")
  set.seed(5)
  d <- lba_simulate(120, v = c(2.4, 1.5, 0.9), A = 0.5, k = 0.4, ndt = 0.2)
  nll <- function(par) {
    b <- exp(par$lA) + exp(par$lk)
    pars <- lapply(1:3, function(j)
      list(v = par$v[j], A = exp(par$lA), b = b, s = 1, posdrift = TRUE))
    -sum(race(d$rt - exp(par$lndt), d$choice, law, pars))
  }
  init <- list(v = c(2.4, 1.5, 0.9), lA = log(0.5), lk = log(0.4),
               lndt = log(0.2))
  obj <- RTMB::MakeADFun(nll, init, silent = TRUE)
  # every point keeps the non-decision time below the fastest response:
  # above it the decision time is negative and there is no density to
  # differentiate, which is the state the bounded ndt link exists to
  # make unreachable in a fit
  expect_gt(min(d$rt), 0.28)
  pts <- list(unlist(init),
              unlist(init) + c(0.6, -0.4, 0.3, 0.2, -0.3, -0.25),
              c(1.7, 1.05, 0.6, log(0.62), log(0.53), log(0.15)),
              c(3.5, 0.2, -0.4, log(0.9), log(0.15), log(0.05)))
  for (p in pts) {
    expect_lt(max(abs(as.numeric(obj$gr(p)) - numDeriv::grad(obj$fn, p))), 1e-6)
  }
})

test_that("value and gradient agree in the far fast tail", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")
  # This is the regime the subtractive form could not survive. Once
  # pnorm(h) rounds to one it returns exactly zero, so the density loses
  # a term while the tape keeps differentiating the term that was
  # written, and value and gradient stop describing the same surface.
  # At this exact point that gap was 239.7 in absolute terms, 6 percent
  # of the gradient. The log-space form in lba_phidiff() removes it.
  set.seed(5)
  d <- lba_simulate(120, v = c(2.4, 1.5, 0.9), A = 0.5, k = 0.4, ndt = 0.2)
  nll <- function(par) {
    b <- exp(par$lA) + exp(par$lk)
    accs <- lapply(1:3, function(j)
      list(v = par$v[j], A = exp(par$lA), b = b, s = 1, posdrift = TRUE))
    -sum(race(d$rt - exp(par$lndt), d$choice, law, accs))
  }
  obj <- RTMB::MakeADFun(nll, list(v = c(2.4, 1.5, 0.9), lA = log(0.5),
                                   lk = log(0.4), lndt = log(0.2)),
                         silent = TRUE)
  p <- c(1.7, 1.05, 0.6, log(0.62), log(0.53), log(0.25))
  # the point really is in the defective regime: some rows are past 8.3
  h <- (0.53 / (d$rt - 0.25)) - 1.7
  expect_gt(max(h), 8.3)
  expect_gt(sum(h > 8.3), 0)
  expect_lt(max(abs(as.numeric(obj$gr(p)) - numDeriv::grad(obj$fn, p))), 1e-5)
})

test_that("the fitted objective's gradient is zero at the optimum", {
  set.seed(6)
  d <- lba_simulate(600, v = c(2.4, 1.5, 0.9), A = 0.5, k = 0.4, ndt = 0.2)
  fit <- frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = d)
  expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-3)
})

test_that("a decision time at or below zero gives a wall, not a NaN", {
  pars <- lapply(c(2, 1.2), function(v) acc(v))
  lp <- race(c(-1, 0, 1e-9, 0.4), c(1, 1, 2, 2), law, pars)
  expect_true(all(is.finite(lp)))
  expect_true(all(lp[1:3] < -100))
  # predict() on new data faster than anything in the training set holds
  # the training bound, so the row can land there and must not take the
  # whole prediction with it
  set.seed(17)
  d <- lba_simulate(120, v = c(2.2, 1.3, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  fit <- frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = d)
  nd <- d[1:5, ]
  nd$rt <- min(d$rt) / 4
  expect_true(all(is.finite(logLik(fit))))
  expect_silent(p <- predict(fit, newdata = nd, type = "link"))
})

# ------------------------------------------------------------------ (d)
test_that("three accumulators and a covariate on one drift recover", {
  skip_on_cran()
  set.seed(2024)
  truth <- c(v1 = 2.4, v2 = 1.5, v3 = 0.9, b_x = 0.9, A = 0.5, k = 0.4,
             ndt = 0.2)
  R <- 12L
  N <- 1200L
  est <- matrix(NA_real_, R, 4,
                dimnames = list(NULL, c("v1", "v2_int", "v2_x", "v3")))
  for (r in seq_len(R)) {
    x <- stats::rnorm(N)
    V <- cbind(truth[["v1"]],
               truth[["v2"]] + truth[["b_x"]] * x,
               truth[["v3"]])
    d <- lba_simulate(N, v = V, A = truth[["A"]], k = truth[["k"]],
                      ndt = truth[["ndt"]])
    d$x <- x
    f <- try(frm(bf(rt | vint(choice) ~ 1, v2 ~ x), family = lba(3),
                 data = d), silent = TRUE)
    if (inherits(f, "try-error")) next
    e <- fixef(f)
    est[r, ] <- c(e$v1[[1]], e$v2[[1]], e$v2[[2]], e$v3[[1]])
  }
  ok <- stats::complete.cases(est)
  expect_gte(sum(ok), 10L)
  est <- est[ok, , drop = FALSE]
  tru <- c(truth[["v1"]], truth[["v2"]], truth[["b_x"]], truth[["v3"]])
  for (j in seq_len(4)) {
    mcse <- stats::sd(est[, j]) / sqrt(nrow(est))
    # the estimator is consistent, so the mean over replicates should sit
    # within a few Monte Carlo standard errors of the truth
    expect_lt(abs(mean(est[, j]) - tru[j]), max(4 * mcse, 0.06))
  }
})

test_that("a covariate on one drift moves that drift and not the others", {
  set.seed(31)
  N <- 1500L
  x <- stats::rnorm(N)
  V <- cbind(2.4, 1.5 + 1.2 * x, 0.9)
  d <- lba_simulate(N, v = V, A = 0.5, k = 0.4, ndt = 0.2)
  d$x <- x
  fit <- frm(bf(rt | vint(choice) ~ x), family = lba(3), data = d)
  e <- fixef(fit)
  # every drift got its own copy of the formula; only the second should
  # have found a slope, which is the capability the family exists for
  expect_gt(e$v2[["x"]], 0.8)
  expect_lt(abs(e$v1[["x"]]), 0.35)
  expect_lt(abs(e$v3[["x"]]), 0.35)
})

# ------------------------------------------------------------------ (e)
test_that("two accumulators behave consistently with the Wiener family", {
  skip_on_cran()
  # The LBA and the diffusion are different models and will not agree on
  # a likelihood. What must agree is the qualitative behaviour: raising
  # one drift raises that response's share and shortens its mean time.
  set.seed(77)
  N <- 30000L
  share <- numeric(0); mt_win <- numeric(0); mt_lose <- numeric(0)
  gaps <- c(0, 0.4, 0.9, 1.6)
  for (g in gaps) {
    d <- lba_simulate(N, v = c(1.4 + g, 1.4), A = 0.5, k = 0.4, ndt = 0.2)
    share <- c(share, mean(d$choice == 1))
    mt_win <- c(mt_win, mean(d$rt[d$choice == 1]))
    mt_lose <- c(mt_lose, mean(d$rt[d$choice == 2]))
  }
  # at a zero gap the two accumulators are exchangeable
  expect_equal(share[1], 0.5, tolerance = 0.02)
  expect_equal(mt_win[1], mt_lose[1], tolerance = 0.02)
  # the favoured response gets commoner and faster as the gap opens
  expect_true(all(diff(share) > 0))
  expect_true(all(diff(mt_win) < 0))
  # and the errors get slower than the correct responses, the ordering a
  # diffusion with a start-point bias also produces
  expect_true(all(mt_lose[-1] > mt_win[-1]))

  # the same ordering out of the two-boundary diffusion, for comparison
  dd <- ddm_simulate(N, mu = 0.8, bs = 1.4, ndt = 0.2, bias = 0.5)
  expect_gt(mean(dd$upper == 1), 0.5)
  expect_lt(mean(dd$rt[dd$upper == 1]), mean(dd$rt[dd$upper == 0]))
})

test_that("a two-accumulator LBA fit is not a Wiener fit", {
  # stated rather than assumed: the same data give different likelihoods,
  # so the cross-check above is about behaviour, not about parity
  set.seed(12)
  d <- lba_simulate(800, v = c(2.2, 1.2), A = 0.5, k = 0.4, ndt = 0.2)
  d$upper <- as.integer(d$choice == 1)
  fl <- frm(bf(rt | vint(choice) ~ 1), family = lba(2), data = d)
  fw <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(), data = d)
  expect_false(isTRUE(all.equal(as.numeric(logLik(fl)),
                                as.numeric(logLik(fw)), tolerance = 1e-3)))
})

# ------------------------------------------------------------------ (f)
test_that("the constructor refuses an impossible accumulator count", {
  expect_error(lba(), "number of accumulators")
  expect_error(lba(1), "number of accumulators")
  expect_error(lba(2.5), "number of accumulators")
  expect_error(lba(c(2, 3)), "number of accumulators")
})

test_that("the constructor refuses a drift scale that cannot identify the model", {
  expect_error(lba(3, sd_v = 0), "fixes the scale")
  expect_error(lba(3, sd_v = -1), "fixes the scale")
  expect_error(lba(3, sd_v = c(1, 1)), "fixes the scale")
  expect_silent(lba(3, sd_v = c(1, 1.2, 0.8)))
})

test_that("the constructor refuses a non-flag posdrift and a bad ndt bound", {
  expect_error(lba(3, posdrift = 1), "truncated at zero")
  expect_error(lba(3, max_ndt = 0), "bounds the non-decision time")
  expect_error(lba(3, max_ndt = c(0.2, 0.3)), "bounds the non-decision time")
})

test_that("a choice outside 1..n is refused by name", {
  set.seed(8)
  d <- lba_simulate(50, v = c(2, 1, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  d$bad <- d$choice
  d$bad[3] <- 4L
  expect_error(frm(bf(rt | vint(bad) ~ 1), family = lba(3), data = d),
               "whole number from 1 to 3")
  d$bad[3] <- 0L
  expect_error(frm(bf(rt | vint(bad) ~ 1), family = lba(3), data = d),
               "whole number from 1 to 3")
  # a fractional choice never reaches the family: vint() is an integer
  # term and frmtmb refuses it first, which is the better place for it
  d$bad <- as.numeric(d$choice) + 0.5
  expect_error(frm(bf(rt | vint(bad) ~ 1), family = lba(3), data = d),
               "must be integers")
})

test_that("a non-positive response time is refused", {
  set.seed(9)
  d <- lba_simulate(50, v = c(2, 1, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  d$rt[5] <- 0
  expect_error(frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = d),
               "strictly positive, finite")
})

test_that("a missing choice indicator is refused rather than silently dropped", {
  set.seed(10)
  d <- lba_simulate(50, v = c(2, 1, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  # declared through required_aterms, so the refusal comes from frmtmb and
  # names both the family and the spelling that supplies the term
  expect_error(frm(bf(rt ~ 1), family = lba(3), data = d), "vint1")
})

test_that("a threshold inside the start-point range cannot be reached", {
  set.seed(11)
  d <- lba_simulate(50, v = c(2, 1, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  # b = A + k with a log link on k, so b > A at every linear predictor;
  # the only way to ask for b <= A is to pin k, and the link's range
  # check refuses that
  expect_error(frm(bf(rt | vint(choice) ~ 1, k = 0), family = lba(3),
                   data = d), "log link")
  # log(-0.2) warns on its way to the refusal, which is the link range
  # check doing its job rather than anything this family should catch
  expect_error(suppressWarnings(
    frm(bf(rt | vint(choice) ~ 1, k = -0.2), family = lba(3), data = d)),
    "log link")
})

test_that("a non-decision-time bound above the fastest response is refused", {
  set.seed(13)
  d <- lba_simulate(80, v = c(2, 1, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  expect_error(frm(bf(rt | vint(choice) ~ 1),
                   family = lba(3, max_ndt = max(d$rt)), data = d),
               "above the fastest response")
  expect_silent(fit <- frm(bf(rt | vint(choice) ~ 1),
                           family = lba(3, max_ndt = min(d$rt)), data = d))
})

test_that("lba_simulate refuses arguments it cannot draw from", {
  expect_error(lba_simulate(0, v = c(1, 2)), "number of trials")
  expect_error(lba_simulate(10, v = 1), "at least")
  expect_error(lba_simulate(10, v = c(1, 2), A = 0), "must both be positive")
  expect_error(lba_simulate(10, v = c(1, 2), k = -1), "must both be positive")
  expect_error(lba_simulate(10, v = c(1, 2), sd_v = 0), "positive number")
  expect_error(lba_simulate(10, v = matrix(1, 3, 2)), "row of drift means")
})

# ------------------------------------------------------- family plumbing
test_that("the fixed drift scale is visible on the family", {
  f <- lba(3, sd_v = 1.5)
  expect_equal(f[["lba_sd_v"]], rep(1.5, 3))
  expect_equal(f[["lba_n"]], 3L)
  expect_true(f[["lba_posdrift"]])
  expect_equal(f$dpars, c("v1", "v2", "v3", "A", "k", "ndt"))
  expect_equal(f$primary_dpars, c("v1", "v2", "v3"))
})

test_that("the non-decision-time link is derived from the response", {
  set.seed(14)
  d <- lba_simulate(60, v = c(2, 1, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  fit <- frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = d)
  lk <- stats::family(fit)$links$ndt
  expect_match(lk$name, "^scaled_logit")
  # the link cannot reach its own bound at any finite linear predictor,
  # so the non-decision time stays strictly below the fastest response
  expect_lt(lk$linkinv(1e4), min(d$rt) + 1e-12)
  expect_lt(lk$linkinv(fixef(fit)$ndt[[1]]), min(d$rt))
})

test_that("simulate() redraws times holding each row's observed choice", {
  set.seed(15)
  d <- lba_simulate(200, v = c(2.4, 1.4, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
  fit <- frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = d)
  s <- simulate(fit, nsim = 2)
  expect_equal(dim(s), c(200L, 2L))
  expect_true(all(is.finite(as.matrix(s))))
  expect_true(all(as.matrix(s) > 0))
  # conditioning on the choice makes the per-choice mean times track the
  # observed ones, which an unconditional draw would not
  for (j in 1:3) {
    expect_equal(mean(s[[1]][d$choice == j]), mean(d$rt[d$choice == j]),
                 tolerance = 0.15)
  }
})

test_that("the untruncated convention is a different model, not a rescaling", {
  skip_if_not_installed("rtdists")
  set.seed(16)
  d <- lba_simulate(400, v = c(2.2, 1.3), A = 0.5, k = 0.4, ndt = 0.2)
  ft <- frm(bf(rt | vint(choice) ~ 1), family = lba(2), data = d)
  fu <- frm(bf(rt | vint(choice) ~ 1), family = lba(2, posdrift = FALSE),
            data = d)
  expect_false(isTRUE(all.equal(as.numeric(logLik(ft)),
                                as.numeric(logLik(fu)), tolerance = 1e-4)))
  # and the untruncated fit agrees with rtdists's own posdrift = FALSE
  e <- fixef(fu)
  A <- exp(e$A[[1]]); k <- exp(e$k[[1]])
  ndt <- stats::family(fu)$links$ndt$linkinv(e$ndt[[1]])
  vv <- c(e$v1[[1]], e$v2[[1]])
  ref <- sum(vapply(seq_len(nrow(d)), function(i) {
    ord <- c(d$choice[i], setdiff(1:2, d$choice[i]))
    log(rtdists::n1PDF(d$rt[i] - ndt, A = A, b = A + k, t0 = 0,
                       mean_v = vv[ord], sd_v = c(1, 1),
                       args.dist = list(posdrift = FALSE), silent = TRUE))
  }, numeric(1)))
  expect_lt(abs(as.numeric(logLik(fu)) - ref), 1e-6)
})
