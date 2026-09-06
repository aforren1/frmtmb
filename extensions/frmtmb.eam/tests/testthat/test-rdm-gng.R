## The racing diffusion model and the go/no-go diffusion model.
##
## Neither family has a canonical R implementation the way lba() has
## rtdists, so the references here are built rather than borrowed, and
## each is chosen to be independent of the thing it is checking:
##
##   RDM     statmod's inverse Gaussian, averaged over the start point
##           with the package's own Gauss-Legendre rule. statmod knows
##           nothing about races and this file's algebra is not in it,
##           so agreement is evidence about the derivation and not
##           about one implementation agreeing with itself.
##   go/no-go  WienR's Wiener distribution function, which is what EMC2
##           itself calls, and RWiener's, which is a second and
##           independent one.
##
## EMC2 is the package these two families exist to match, but every
## function that computes either likelihood in it is internal. The
## comparison against it therefore lives in
## dev/rdm-gng-emc2-reference.R and its results are written up in
## dev/rdm-gng-findings.md; nothing here reaches into another package's
## namespace.

law <- frmtmb.eam:::rdm_law
race <- frmtmb.eam:::lba_race_lpdf
acc <- function(v, A = 0.5, k = 0.5) list(v = v, A = A, k = k)

# The reference average over the distance to travel. The rule is the
# package's own, but the LAW is statmod's, which is the half that
# matters: a start-point average of an inverse Gaussian is what the
# model says, and this computes it without any of the log-space algebra
# the family uses.
rdm_ref <- function(t, v, A, k, fn, n = 60L) {
  gl <- frmtmb.eam:::ddm_gauss_legendre(n)
  vapply(seq_along(t), function(i) {
    d <- k[i] + A[i] * gl$x
    sum(gl$w * fn(t[i], d, v[i]))
  }, numeric(1))
}
ref_dens <- function(t, d, v) statmod::dinvgauss(t, mean = d / v, shape = d^2)
ref_surv <- function(t, d, v) {
  statmod::pinvgauss(t, mean = d / v, shape = d^2, lower.tail = FALSE)
}

# ------------------------------------------------------------ RDM (a)

test_that("the single-accumulator density is the averaged inverse Gaussian", {
  skip_if_not_installed("statmod")
  gr <- expand.grid(t = c(0.05, 0.2, 0.5, 1, 2, 5), v = c(0.3, 1, 2, 4),
                    k = c(0.3, 1, 2.5), A = c(0.1, 0.5, 1.5))
  mine <- exp(law$ldens(gr$t, acc(gr$v, gr$A, gr$k)))
  ref <- rdm_ref(gr$t, gr$v, gr$A, gr$k, ref_dens)
  keep <- ref > 1e-250
  expect_gt(sum(keep), 200)
  expect_lt(max(abs(mine[keep] - ref[keep]) / ref[keep]), 1e-9)
  # and nothing on the grid underflowed to a hole
  expect_identical(sum(mine == 0), 0L)
})

test_that("the survival is the averaged inverse-Gaussian upper tail", {
  skip_if_not_installed("statmod")
  gr <- expand.grid(t = c(0.05, 0.2, 0.5, 1, 2, 5), v = c(0.3, 1, 2, 4),
                    k = c(0.3, 1, 2.5), A = c(0.1, 0.5, 1.5))
  mine <- exp(law$lsurv(gr$t, acc(gr$v, gr$A, gr$k)))
  ref <- rdm_ref(gr$t, gr$v, gr$A, gr$k, ref_surv)
  keep <- ref > 1e-12
  expect_gt(sum(keep), 200)
  expect_lt(max(abs(mine[keep] - ref[keep]) / ref[keep]), 1e-9)
})

test_that("the survival keeps its digits where one minus the CDF loses them", {
  skip_if_not_installed("statmod")
  # A fast accumulator that has almost certainly finished. This is an
  # ordinary row of a race, not an extreme one: a loser that has
  # probably-but-not-certainly finished is what every trial with a
  # clear winner has three of.
  #
  # The failure of the subtractive form is sharper than returning zero.
  # It STICKS at the rounding of one, about 8.9e-16, and stays there
  # however far into the tail it is asked, so it reports the same
  # survival at t = 3 as at t = 20 while the true value falls by 238
  # orders of magnitude.
  p <- acc(v = 8, A = 0.1, k = 0.3)
  tt <- c(3, 5, 10, 20)
  gl <- frmtmb.eam:::ddm_gauss_legendre(60L)
  subtractive <- vapply(tt, function(q) {
    d <- 0.3 + 0.1 * gl$x
    1 - sum(gl$w * statmod::pinvgauss(q, mean = d / 8, shape = d^2))
  }, numeric(1))
  # stuck: every entry is the same, and it is the rounding of one
  expect_lt(stats::sd(subtractive), 1e-20)
  expect_lt(max(subtractive), 1e-14)

  # the reference that still has digits: the upper tail asked for
  # directly, so that nothing is subtracted from one anywhere
  direct <- rdm_ref(tt, rep(8, 4), rep(0.1, 4), rep(0.3, 4), ref_surv)
  expect_true(all(direct > 0))
  expect_lt(min(direct), 1e-280)
  # and the subtractive form is wrong by more than 250 orders of
  # magnitude at the far end, which is the point
  expect_gt(max(subtractive) / min(direct), 1e250)

  mine <- exp(law$lsurv(tt, p))
  expect_true(all(mine > 0))
  # 1.4e-08 as measured, on values down to 1e-282. The closed form loses
  # relative accuracy that far out and the rate is recorded in
  # dev/rdm-gng-findings.md; what matters here is that it is still
  # tracking the tail rather than sitting on a floor.
  expect_lt(max(abs(mine - direct) / direct), 1e-7)
})

test_that("the race is the winner's density times the losers' survivals", {
  skip_if_not_installed("statmod")
  # The composition EMC2 uses for every race family it has, written out
  # against the independent law rather than against this package's.
  set.seed(7)
  A <- 0.4; k <- 0.6
  for (nacc in c(2L, 3L, 4L)) {
    m <- 30L
    tt <- stats::runif(m, 0.2, 2.5)
    V <- matrix(stats::runif(m * nacc, 0.4, 3.5), ncol = nacc)
    win <- rep_len(seq_len(nacc), m)
    pars <- lapply(seq_len(nacc), function(j) acc(V[, j], A, k))
    mine <- race(tt, win, law, pars)
    ref <- vapply(seq_len(m), function(i) {
      lw <- log(rdm_ref(tt[i], V[i, win[i]], A, k, ref_dens))
      ls <- sum(vapply(setdiff(seq_len(nacc), win[i]), function(j) {
        log(rdm_ref(tt[i], V[i, j], A, k, ref_surv))
      }, numeric(1)))
      lw + ls
    }, numeric(1))
    expect_lt(max(abs(mine - ref)), 1e-8)
  }
})

test_that("the family's log-density is the race, through a fitted object", {
  skip_if_not_installed("statmod")
  set.seed(21)
  dat <- rdm_simulate(200, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5,
                      ndt = 0.2)
  fit <- frm(bf(rt | vint(choice) ~ 1), family = rdm(3), data = dat)
  e <- fixef(fit)
  A <- exp(e$A[[1]]); k <- exp(e$k[[1]])
  ndt <- stats::family(fit)$links$ndt$linkinv(e$ndt[[1]])
  vv <- exp(c(e$v1[[1]], e$v2[[1]], e$v3[[1]]))
  ref <- vapply(seq_len(nrow(dat)), function(i) {
    t <- dat$rt[i] - ndt
    j <- dat$choice[i]
    log(rdm_ref(t, vv[j], A, k, ref_dens)) +
      sum(vapply(setdiff(1:3, j), function(m) {
        log(rdm_ref(t, vv[m], A, k, ref_surv))
      }, numeric(1)))
  }, numeric(1))
  expect_lt(abs(as.numeric(logLik(fit)) - sum(ref)), 1e-6)
})

# ------------------------------------------------------------ RDM (b)

test_that("the choice probabilities sum to one", {
  # The race always ends: every drift is positive, so every accumulator
  # arrives eventually and no mass escapes. This is the defective
  # density's mass check, and it is the one that would catch a missing
  # or double-counted factor in the survival.
  A <- 0.5; k <- 0.6; vv <- c(2.5, 1.4, 0.8)
  pars <- lapply(vv, function(v) acc(v, A, k))
  tot <- sum(vapply(seq_along(vv), function(j) {
    stats::integrate(function(q) exp(race(q, rep(j, length(q)), law, pars)),
                     1e-9, Inf, rel.tol = 1e-10,
                     subdivisions = 5000L)$value
  }, numeric(1)))
  expect_equal(tot, 1, tolerance = 1e-6)
})

test_that("the analytic joint distribution reproduces the generative process", {
  skip_on_cran()
  set.seed(99)
  A <- 0.5; k <- 0.6; ndt <- 0.15; vv <- c(3.0, 2.0, 1.2)
  N <- 200000L
  d <- rdm_simulate(N, v = vv, A = A, k = k, ndt = ndt)
  pars <- lapply(vv, function(v) acc(v, A, k))
  edges <- c(0.15, 0.3, 0.45, 0.7, 1.2, Inf)
  worst <- 0
  for (j in 1:3) {
    for (m in seq_len(length(edges) - 1L)) {
      lo <- edges[m]; hi <- min(edges[m + 1L], 60)
      emp <- mean(d$choice == j & d$rt > lo & d$rt <= hi)
      ana <- stats::integrate(function(q) {
        exp(race(q - ndt, rep(j, length(q)), law, pars))
      }, lo, hi, rel.tol = 1e-10, subdivisions = 3000L)$value
      se <- sqrt(max(emp * (1 - emp), 1e-12) / N)
      worst <- max(worst, abs(emp - ana) / se)
    }
  }
  # 15 cells, an exact analytic value, so the only error is Monte Carlo
  expect_lt(worst, 4.5)
})

test_that("fixing the diffusion coefficient at one identifies the scale", {
  skip_if_not_installed("statmod")
  # The invariance the documentation claims, checked in both directions.
  # Scaling the distances and the drifts together leaves the model alone
  # ONLY if the diffusion coefficient scales with them; this family holds
  # it at one, so the scale is identified rather than free.
  tt <- c(0.3, 0.8, 1.6)
  base <- acc(2.2, A = 0.5, k = 0.5)
  for (cc in c(1.5, 2, 3)) {
    scaled <- acc(2.2 * cc, A = 0.5 * cc, k = 0.5 * cc)
    # not a reparameterization: with the diffusion held at one these are
    # different models, and increasingly so as the constant grows
    expect_gt(max(abs(law$lsurv(tt, base) - law$lsurv(tt, scaled))), 1)
  }
  # and the invariance really is the one claimed, checked through a law
  # that carries the diffusion coefficient explicitly. The first passage
  # to a distance d under drift v and diffusion s is inverse Gaussian
  # with mean d/v and shape (d/s)^2, so scaling d, v and s by one
  # constant leaves both alone.
  gl <- frmtmb.eam:::ddm_gauss_legendre(60L)
  avg <- function(t, A, k, v, s) {
    d <- k + A * gl$x
    sum(gl$w * statmod::pinvgauss(t, mean = d / v, shape = (d / s)^2,
                                  lower.tail = FALSE))
  }
  a1 <- vapply(tt, function(q) avg(q, 0.5, 0.5, 2.2, 1), numeric(1))
  a2 <- vapply(tt, function(q) avg(q, 1.0, 1.0, 4.4, 2), numeric(1))
  expect_equal(a1, a2, tolerance = 1e-12)
  # and the s = 1 member of that family is this one
  expect_equal(log(a1), as.numeric(law$lsurv(tt, base)), tolerance = 1e-9)
})

# ------------------------------------------------------------ RDM (c)

test_that("the taped gradient matches numDeriv at several points", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")
  set.seed(5)
  d <- rdm_simulate(150, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
  nll <- function(par) {
    pars <- lapply(1:3, function(j)
      acc(exp(par$lv[j]), exp(par$lA), exp(par$lk)))
    -sum(race(d$rt - exp(par$lndt), d$choice, law, pars))
  }
  init <- list(lv = log(c(3.0, 2.0, 1.2)), lA = log(0.5), lk = log(0.5),
               lndt = log(0.2))
  obj <- RTMB::MakeADFun(nll, init, silent = TRUE)
  expect_gt(min(d$rt), 0.21)
  pts <- list(unlist(init),
              unlist(init) + c(0.3, -0.25, 0.2, 0.15, -0.2, -0.3),
              c(log(4.5), log(1.1), log(0.6), log(0.8), log(0.3),
                log(0.12)),
              c(log(1.2), log(0.9), log(0.5), log(0.2), log(1.2),
                log(0.05)))
  for (p in pts) {
    expect_lt(max(abs(as.numeric(obj$gr(p)) - numDeriv::grad(obj$fn, p))),
              1e-6)
  }
})

test_that("the fitted objective's gradient is zero at the optimum", {
  set.seed(6)
  d <- rdm_simulate(800, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
  fit <- frm(bf(rt | vint(choice) ~ 1), family = rdm(3), data = d)
  expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-3)
})

test_that("a decision time at or below zero gives a wall, not a NaN", {
  pars <- lapply(c(2.5, 1.4), function(v) acc(v))
  lp <- race(c(-1, 0, 1e-9, 0.4), c(1, 1, 2, 2), law, pars)
  expect_true(all(is.finite(lp)))
  expect_true(all(lp[1:3] < -100))
  # predict() on new data faster than anything in training holds the
  # training bound, so the row can land there and must not take the
  # whole prediction with it
  set.seed(17)
  d <- rdm_simulate(120, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
  fit <- frm(bf(rt | vint(choice) ~ 1), family = rdm(3), data = d)
  nd <- d[1:5, ]
  nd$rt <- min(d$rt) / 4
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_silent(p <- predict(fit, newdata = nd, type = "link"))
})

# ------------------------------------------------------------ RDM (d)

test_that("three accumulators and a covariate on one drift recover", {
  skip_on_cran()
  set.seed(2024)
  truth <- c(v1 = 3.0, v2 = 2.0, v3 = 1.2, b_x = 0.35, A = 0.5, k = 0.5,
             ndt = 0.2)
  R <- 12L
  N <- 1500L
  est <- matrix(NA_real_, R, 4,
                dimnames = list(NULL, c("v1", "v2_int", "v2_x", "v3")))
  for (r in seq_len(R)) {
    x <- stats::rnorm(N)
    V <- cbind(truth[["v1"]],
               truth[["v2"]] * exp(truth[["b_x"]] * x),
               truth[["v3"]])
    d <- rdm_simulate(N, v = V, A = truth[["A"]], k = truth[["k"]],
                      ndt = truth[["ndt"]])
    d$x <- x
    f <- try(frm(bf(rt | vint(choice) ~ 1, v2 ~ x), family = rdm(3),
                 data = d), silent = TRUE)
    if (inherits(f, "try-error")) next
    e <- fixef(f)
    est[r, ] <- c(e$v1[[1]], e$v2[[1]], e$v2[[2]], e$v3[[1]])
  }
  ok <- stats::complete.cases(est)
  expect_gte(sum(ok), 10L)
  est <- est[ok, , drop = FALSE]
  # every drift has a log link, so the intercepts are compared on it and
  # the slope is already there
  tru <- c(log(truth[["v1"]]), log(truth[["v2"]]), truth[["b_x"]],
           log(truth[["v3"]]))
  for (j in seq_len(4)) {
    mcse <- stats::sd(est[, j]) / sqrt(nrow(est))
    expect_lt(abs(mean(est[, j]) - tru[j]), max(4 * mcse, 0.08))
  }
})

test_that("a covariate on one drift moves that drift and not the others", {
  set.seed(31)
  N <- 2000L
  x <- stats::rnorm(N)
  V <- cbind(3.0, 2.0 * exp(0.5 * x), 1.2)
  d <- rdm_simulate(N, v = V, A = 0.5, k = 0.5, ndt = 0.2)
  d$x <- x
  fit <- frm(bf(rt | vint(choice) ~ x), family = rdm(3), data = d)
  e <- fixef(fit)
  expect_gt(e$v2[["x"]], 0.3)
  expect_lt(abs(e$v1[["x"]]), 0.2)
  expect_lt(abs(e$v3[["x"]]), 0.2)
})

test_that("rdm refuses what it cannot read", {
  set.seed(3)
  d <- rdm_simulate(60, v = c(3, 2), A = 0.5, k = 0.5, ndt = 0.2)
  expect_error(rdm(1), "2 or more")
  expect_error(rdm(3, max_ndt = -1), "positive finite number")
  # the winner is mandatory, and named
  expect_error(frm(bf(rt ~ 1), family = rdm(2), data = d), "vint")
  # a two-level decision indicator is not a winner index
  d$bad <- d$choice - 1L
  expect_error(frm(bf(rt | vint(bad) ~ 1), family = rdm(2), data = d),
               "whole number from 1 to 2")
  expect_error(rdm_simulate(10, v = c(1, -1)), "positive and finite")

  # dec() is refused rather than silently dropped.
  d$two <- (d$choice > 1) + 0L
  expect_error(frm(bf(rt | dec(two) + vint(choice) ~ 1), family = rdm(3),
                   data = d), "cannot mean anything here")
})

test_that("both race families refuse dec(), and used not to agree", {
  # The regression. Before the punch round lba() FITTED a model written
  # with dec() alongside vint(), dropping the term with no warning and
  # estimates bit-identical to the model without it, while rdm() on the
  # same formula refused. Two families with the same dpars and the same
  # formula interface should not disagree about whether a term means
  # anything, and the direction of the old disagreement was the harmful
  # one: dec() IS the spelling under wiener(), so a ported model fitted
  # while quietly ignoring half of what its author wrote.
  set.seed(11)
  d <- lba_simulate(120, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.4, ndt = 0.2)
  d$two <- (d$choice > 1) + 0L
  expect_error(frm(bf(rt | dec(two) + vint(choice) ~ 1), family = lba(3),
                   data = d), "cannot mean anything here")
  # one shared template, so both refusals name their own family
  expect_error(frm(bf(rt | dec(two) + vint(choice) ~ 1), family = lba(3),
                   data = d), "^lba:")
  r <- rdm_simulate(120, v = c(3, 2, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
  r$two <- (r$choice > 1) + 0L
  expect_error(frm(bf(rt | dec(two) + vint(choice) ~ 1), family = rdm(3),
                   data = r), "^rdm:")
  # and neither refusal costs the ordinary model anything
  expect_no_error(frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = d))
  expect_no_error(frm(bf(rt | vint(choice) ~ 1), family = rdm(3), data = r))
})

# ------------------------------------------------------- go/no-go (a)

nogo <- function(t, v, a, w) exp(frmtmb.eam:::ddm_nogo_lprob(t, v, a, w))

test_that("the no-go probability matches an independent Wiener CDF", {
  skip_if_not_installed("WienR")
  # w = 0.9 is here because of the punch round. The first version of
  # this grid stopped at 0.7 and filtered to ref > 1e-6, and the blend
  # constant's worst corner needs a large drift AND a start point
  # strongly biased toward the go boundary AND a wide boundary, so
  # nothing in the suite could see it.
  gr <- expand.grid(t = c(0.05, 0.2, 0.6, 1.5, 4),
                    a = c(0.3, 0.8, 1.4, 2.5, 4),
                    v = c(-3, -1, 0, 1.5, 4),
                    w = c(0.25, 0.3, 0.5, 0.7, 0.9))
  ref <- 1 - WienR::pWDM(gr$t, rep("upper", nrow(gr)), a = gr$a, v = gr$v,
                         w = gr$w, t0 = 0, precision = 1e-12)$value
  mine <- nogo(gr$t, gr$v, gr$a, gr$w)
  # every value is a probability, on every row of the grid and not just
  # the ones the reference can score
  expect_true(all(mine >= 0 & mine <= 1))

  # Stratified rather than filtered. WienR at precision 1e-12 is itself
  # only good to about 4 percent once the no-go probability falls to
  # 1e-13 -- measured, and the reason the sharp corner is pinned
  # against a 260-bit constant in the next test rather than here.
  tight <- ref > 1e-6
  expect_gt(sum(tight), 400)
  expect_lt(max(abs(mine[tight] - ref[tight]) / ref[tight]), 1e-9)

  # Below that the assertion is on the REFERENCE's terms. Only 7 rows of
  # this grid land there at all, and the worst of them is a no-go
  # probability of 3.0e-13 where WienR is 1.4e-04 away from this
  # family; a part in a thousand is what WienR has left, not what this
  # family has. The corner is really guarded by the 260-bit constants
  # in the next test, which need no reference package.
  loose <- ref > 0 & ref <= 1e-6
  expect_gte(sum(loose), 5)
  expect_lt(max(abs(mine[loose] - ref[loose]) / ref[loose]), 1e-3)
})

test_that("the blend holds up where the small-time route collapses", {
  # The punch-round regression, and the sharpest assertion in this file.
  #
  # These reference values are 260-bit computations, checked by two
  # INDEPENDENT high-precision routes that agree with each other to
  # between 8e-64 and 7e-77 on these very rows, so they are constants
  # rather than another implementation of the thing under test. No
  # Rmpfr at test time.
  #
  # The corner they pin needs a large drift AND a start point strongly
  # biased toward the go boundary AND a wide boundary, all at once. It
  # is where the small-time route stops being merely inaccurate and
  # becomes hopeless: at the first row it is 12.1 RELATIVE, because it
  # forms 1 - F_upper and there is nothing left to subtract from. With
  # the blend centred at u0 = 0.06 that route still carried a share of
  # about 3e-05 there and dragged the answer to 7.75e-05 relative; at
  # 0.02 the weight has saturated and the answer is the large-time
  # route's own.
  ref <- data.frame(
    t = c(2.5, 2.0, 2.5, 1.0, 0.6, 4.0),
    v = c(5, 5, 5, 5, 3, 2),
    a = c(4, 4, 4, 4, 2.5, 1.4),
    w = c(0.90, 0.90, 0.75, 0.90, 0.85, 0.80),
    truth = c(8.2014416716939391e-16, 4.2147569533687343e-13,
              1.1883283058299342e-13, 2.9321985692509448e-07,
              0.0092538804834898583, 0.0076638946359617807))
  mine <- nogo(ref$t, ref$v, ref$a, ref$w)
  rel <- abs(mine - ref$truth) / ref$truth
  expect_lt(max(rel), 1e-11)
  # and the worst of them, the one the review found, individually
  expect_lt(rel[1], 1e-11)
  expect_equal(mine[1], ref$truth[1], tolerance = 1e-11)

  # The mechanism, asserted rather than described: on that row the
  # large-time route alone is right and the small-time route is not,
  # so the blend is only correct because the weight has saturated onto
  # the large one.
  lrg <- frmtmb.eam:::ddm_nogo_large(ref$t[1], ref$v[1], ref$a[1], ref$w[1])
  sml <- frmtmb.eam:::ddm_nogo_small(ref$t[1], ref$v[1], ref$a[1], ref$w[1])
  expect_lt(abs(lrg - ref$truth[1]) / ref$truth[1], 1e-13)
  expect_gt(abs(sml - ref$truth[1]) / ref$truth[1], 1)
})

test_that("the blend centre is on the correct side of the crossing", {
  # u0 is a measured constant, not a taste. Pinned so that a later edit
  # has to re-measure rather than nudge: below the crossing near
  # u = 0.025 the SMALL route is the accurate one, so an over-low
  # centre hands over too early and is as wrong as an over-high one.
  expect_equal(frmtmb.eam:::ddm_cdf_u0, 0.02)
  expect_equal(frmtmb.eam:::ddm_cdf_us, 0.12)
  # the weight really has saturated at the pinned corner
  u <- 2.5 / 4^2
  lam <- 0.5 * (1 + tanh((log(u) - log(frmtmb.eam:::ddm_cdf_u0)) /
                           frmtmb.eam:::ddm_cdf_us))
  expect_gt(lam, 1 - 1e-13)
  # and it has NOT saturated the other way at a genuinely small-time row
  u2 <- 0.05 / 4^2
  lam2 <- 0.5 * (1 + tanh((log(u2) - log(frmtmb.eam:::ddm_cdf_u0)) /
                            frmtmb.eam:::ddm_cdf_us))
  expect_lt(lam2, 1e-13)
})

test_that("RWiener agrees too, at the accuracy RWiener has", {
  skip_if_not_installed("RWiener")
  q <- c(0.1, 0.3, 0.8, 1.5, 3)
  for (w in c(0.35, 0.5, 0.65)) {
    ref <- 1 - vapply(q, function(x) {
      RWiener::pwiener(x, 1.4, 1e-9, w, 1, resp = "upper")
    }, numeric(1))
    expect_lt(max(abs(nogo(q, 1, 1.4, w) - ref) / ref), 1e-7)
  }
})

test_that("the two series agree with each other in the overlap band", {
  # A reference-free check that the DERIVATION is right rather than that
  # one implementation agrees with itself: the small-time image sum and
  # the large-time eigenfunction sum are separate derivations of the
  # same quantity and share no algebra.
  gr <- expand.grid(t = c(0.3, 0.8, 2), a = c(0.8, 1.4, 2.5),
                    v = c(-2, 0, 1.5, 4), w = c(0.3, 0.5, 0.7))
  u <- gr$t / gr$a^2
  keep <- u > 0.05 & u < 5              # where both CONVERGE
  expect_gt(sum(keep), 30)
  # At the CONVERGED truncation. The image sum's terms carry
  # exp(-2 j (w + j) / u), which decays in j only for a small u: at
  # u = 5 the j = 5 term is still 1.7e-05, so a wide band is a statement
  # about the derivation and needs enough terms to be one. The shipped
  # truncation is 4, and the next assertion is what covers it.
  s <- frmtmb.eam:::ddm_nogo_small(gr$t[keep], gr$v[keep], gr$a[keep],
                                   gr$w[keep], K = 12L)
  l <- frmtmb.eam:::ddm_nogo_large(gr$t[keep], gr$v[keep], gr$a[keep],
                                   gr$w[keep])
  # measured at 1.47e-08, on a row whose no-go probability is 8.5e-07;
  # the tolerance is set just above what two independent series of a
  # transcendental function actually reach, not at a round number
  expect_lt(max(abs(s - l) / l), 1e-7)

  # The shipped truncation is 4, and the band it has to be right on is
  # not this one. The blend gives the small route a non-zero weight only
  # below u = 0.197, where tanh has not yet saturated; above that the
  # weight is exactly one and the small route's value is multiplied by
  # exactly zero. Inside its own band four terms and twelve agree to
  # every bit, which is the assertion that matters for the likelihood.
  band <- u > 0.001 & u < 0.197
  expect_gt(sum(band), 5)
  s4 <- frmtmb.eam:::ddm_nogo_small(gr$t[band], gr$v[band], gr$a[band],
                                    gr$w[band], K = 4L)
  s12 <- frmtmb.eam:::ddm_nogo_small(gr$t[band], gr$v[band], gr$a[band],
                                     gr$w[band], K = 12L)
  expect_identical(s4, s12)
})

test_that("the no-go probability is finite at every edge a fit can reach", {
  edge <- list(c(1e-12, 1, 1.4, 0.5), c(1e-9, 0, 1.4, 0.5),
               c(50, 5, 0.3, 0.5), c(1, 0, 1, 0.5),
               c(1, -20, 5, 0.99), c(1, 20, 5, 0.01),
               c(1, 1e-9, 1.4, 0.5), c(1e-6, 3, 0.2, 0.2))
  for (p in edge) {
    v <- nogo(p[1], p[2], p[3], p[4])
    expect_true(is.finite(v))
    expect_gte(v, 0)
    expect_lte(v, 1 + 1e-12)
  }
  # zero drift is the one the gambler's-ruin ratio is 0/0 at, and it has
  # a known answer: with an unbiased start point the two boundaries are
  # exchangeable, so the process is as likely to be waiting or below as
  # it is to have gone up
  expect_equal(nogo(1e6, 0, 1, 0.5), 0.5, tolerance = 1e-9)
  expect_equal(nogo(1e6, 0, 1, 0.25), 0.75, tolerance = 1e-9)
})

# ------------------------------------------------------- go/no-go (b)

test_that("the go density and the no-go probability sum to one", {
  # The mass check for this family. The go branch is defective on
  # (ndt, deadline] and the no-go branch carries the rest, so the two
  # together are a probability distribution over the pair
  # (responded, time). A missing reflection or a CDF at the wrong
  # boundary would break this and nothing else would.
  fam <- wiener_gng(deadline = 1.5)
  dp <- list(mu = 1.1, bs = 1.4, ndt = 0.25, bias = 0.45)
  fin <- fam[["family_finalize"]](fam, c(0.4, 1.2), list(dec = c(1, 1)))
  dens <- function(q) {
    exp(fin[["lpdf"]](q, lapply(dp, function(v) rep(v, length(q))),
                      list(dec = rep(1, length(q)))))
  }
  go_mass <- stats::integrate(dens, dp[["ndt"]] + 1e-10, 1.5,
                              rel.tol = 1e-11, subdivisions = 4000L)$value
  no_mass <- exp(fin[["lpdf"]](1.5, dp, list(dec = 0)))
  expect_equal(go_mass + no_mass, 1, tolerance = 1e-8)
  # and each piece is a probability in its own right
  expect_gt(go_mass, 0.5)
  expect_gt(no_mass, 0.05)
})

test_that("the mass check holds across parameters and deadlines", {
  dpars <- list(c(mu = 1.5, bs = 1.2, ndt = 0.2, bias = 0.5),
                c(mu = -0.5, bs = 1.8, ndt = 0.15, bias = 0.35),
                c(mu = 0, bs = 1, ndt = 0.1, bias = 0.6),
                c(mu = 3, bs = 2.2, ndt = 0.3, bias = 0.5))
  for (td in c(0.8, 1.5, 4)) {
    for (p in dpars) {
      fam <- wiener_gng(deadline = td)
      dp <- as.list(p)
      fin <- fam[["family_finalize"]](fam, c(p[["ndt"]] + 0.05, td),
                                      list(dec = c(1, 1)))
      go <- stats::integrate(function(q) {
        exp(fin[["lpdf"]](q, lapply(dp, function(v) rep(v, length(q))),
                          list(dec = rep(1, length(q)))))
      }, p[["ndt"]] + 1e-10, td, rel.tol = 1e-11,
      subdivisions = 4000L)$value
      no <- exp(fin[["lpdf"]](td, dp, list(dec = 0)))
      expect_equal(go + no, 1, tolerance = 1e-7,
                   label = paste0("mass at deadline ", td, ", mu ",
                                  p[["mu"]]))
    }
  }
})

test_that("a go row is wiener()'s own upper-boundary density", {
  # The go branch is not a new density and must not have become one.
  set.seed(4)
  y <- c(0.4, 0.6, 0.9, 1.3)
  dp <- list(mu = 1.1, bs = 1.4, ndt = 0.25, bias = 0.45)
  fam <- wiener_gng(deadline = 2)
  fin <- fam[["family_finalize"]](fam, y, list(dec = rep(1, 4)))
  mine <- fin[["lpdf"]](y, lapply(dp, function(v) rep(v, 4)),
                        list(dec = rep(1, 4)))
  wf <- wiener()
  wfin <- wf[["family_finalize"]](wf, y, list(dec = rep(1, 4)))
  ref <- wfin[["lpdf"]](y, lapply(dp, function(v) rep(v, 4)),
                        list(dec = rep(1, 4)))
  expect_equal(as.numeric(mine), as.numeric(ref), tolerance = 1e-12)
})

# ------------------------------------------------------- go/no-go (c)

test_that("the taped gradient matches numDeriv, both branches", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("numDeriv")
  set.seed(8)
  d <- wiener_gng_simulate(300, mu = 1.0, bs = 1.4, ndt = 0.25,
                           bias = 0.45, deadline = 1.5)
  # both branches really are exercised
  expect_gt(sum(d$responded == 1), 50)
  expect_gt(sum(d$responded == 0), 20)
  lp <- frmtmb.eam:::gng_lpdf
  nll <- function(par) {
    dp <- list(mu = par$mu, bs = exp(par$lbs), ndt = exp(par$lndt),
               bias = 1 / (1 + exp(-par$lbias)))
    -sum(lp(d$rt, dp, list(dec = d$responded), 1.5))
  }
  init <- list(mu = 1.0, lbs = log(1.4), lndt = log(0.25), lbias = 0)
  obj <- RTMB::MakeADFun(nll, init, silent = TRUE)
  pts <- list(unlist(init),
              c(0.4, log(1.9), log(0.18), 0.4),
              c(2.2, log(0.9), log(0.3), -0.5),
              c(-0.6, log(2.4), log(0.12), 0.8))
  for (p in pts) {
    expect_lt(max(abs(as.numeric(obj$gr(p)) - numDeriv::grad(obj$fn, p))),
              1e-6)
  }
})

test_that("the gradient stays finite where a no-go trial is very surprising", {
  skip_if_not_installed("RTMB")
  # The regime the large-time route exists for: a strong drift to the go
  # boundary and a late deadline make a no-go trial almost impossible,
  # and 1 - F would be a subtraction with nothing left in it. The value
  # must stay finite and the gradient with it.
  lp <- frmtmb.eam:::gng_lpdf
  y <- c(0.5, 2.0)
  dec <- c(1, 0)
  f <- function(p) {
    dp <- list(mu = p[1], bs = exp(p[2]), ndt = 0.2, bias = 0.5)
    sum(lp(y, dp, list(dec = dec), 2.5))
  }
  p <- c(8, log(3))
  expect_true(is.finite(f(p)))
  g <- as.numeric(RTMB::MakeTape(f, p)$jacobian(p))
  expect_true(all(is.finite(g)))
  expect_true(all(g != 0))              # a wall, not a flat floor
  # the no-go probability really is tiny there: 3.8e-11, which is well
  # past where the 1 - F route has anything left
  expect_lt(nogo(2.3, 8, 3, 0.5), 1e-9)
})

test_that("the fitted objective's gradient is zero at the optimum", {
  set.seed(9)
  d <- wiener_gng_simulate(1200, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
  fit <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
             family = wiener_gng(deadline = 1.5), data = d)
  expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-3)
})

# ------------------------------------------------------- go/no-go (d)

test_that("the simulator's go rate is the analytic go probability", {
  # An independent check on the whole no-go branch: the simulator runs
  # the diffusion and never touches the distribution function, so the
  # share of trials it records as responses is a Monte Carlo estimate of
  # 1 - nogo.
  set.seed(12)
  N <- 60000L
  for (p in list(c(mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.5, td = 1.5),
                 c(mu = -0.5, bs = 1.2, ndt = 0.2, bias = 0.6, td = 1.0),
                 c(mu = 2.0, bs = 1.8, ndt = 0.15, bias = 0.4, td = 2.0))) {
    d <- wiener_gng_simulate(N, mu = p[["mu"]], bs = p[["bs"]],
                             ndt = p[["ndt"]], bias = p[["bias"]],
                             deadline = p[["td"]])
    emp <- mean(d$responded)
    ana <- 1 - nogo(p[["td"]] - p[["ndt"]], p[["mu"]], p[["bs"]],
                    p[["bias"]])
    se <- sqrt(emp * (1 - emp) / N)
    expect_lt(abs(emp - ana) / se, 4.5,
              label = paste0("go rate at mu ", p[["mu"]]))
  }
})

test_that("the go response times follow the go branch of the density", {
  skip_on_cran()
  set.seed(13)
  N <- 60000L
  d <- wiener_gng_simulate(N, mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.45,
                           deadline = 1.5)
  rt <- d$rt[d$responded == 1]
  fam <- wiener_gng(deadline = 1.5)
  dp <- list(mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.45)
  fin <- fam[["family_finalize"]](fam, rt, list(dec = rep(1, length(rt))))
  dens <- function(q) {
    exp(fin[["lpdf"]](q, lapply(dp, function(v) rep(v, length(q))),
                      list(dec = rep(1, length(q)))))
  }
  edges <- c(0.25, 0.4, 0.5, 0.65, 0.85, 1.1, 1.5)
  obs <- as.numeric(table(cut(rt, edges)))
  ex <- vapply(seq_len(length(edges) - 1L), function(i) {
    stats::integrate(dens, edges[i], edges[i + 1L],
                     rel.tol = 1e-10)$value
  }, numeric(1)) * N
  expect_lt(sum((obs - ex)^2 / ex), stats::qchisq(1 - 1e-6, length(obs)))
})

test_that("wiener_gng recovers its parameters", {
  skip_on_cran()
  set.seed(2025)
  truth <- c(mu = 1.0, bs = 1.4, ndt = 0.25)
  R <- 10L
  N <- 2000L
  est <- matrix(NA_real_, R, 3, dimnames = list(NULL, names(truth)))
  for (r in seq_len(R)) {
    d <- wiener_gng_simulate(N, mu = truth[["mu"]], bs = truth[["bs"]],
                             ndt = truth[["ndt"]], bias = 0.5,
                             deadline = 1.5)
    f <- try(frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
                 family = wiener_gng(deadline = 1.5), data = d),
             silent = TRUE)
    if (inherits(f, "try-error")) next
    e <- fixef(f)
    est[r, ] <- c(e$mu[[1]], exp(e$bs[[1]]),
                  stats::family(f)$links$ndt$linkinv(e$ndt[[1]]))
  }
  ok <- stats::complete.cases(est)
  expect_gte(sum(ok), 8L)
  est <- est[ok, , drop = FALSE]
  for (j in seq_len(3)) {
    mcse <- stats::sd(est[, j]) / sqrt(nrow(est))
    expect_lt(abs(mean(est[, j]) - truth[[j]]), max(4 * mcse, 0.05),
              label = paste0("recovery of ", names(truth)[j]))
  }
})

test_that("a covariate on the drift recovers", {
  set.seed(41)
  N <- 3000L
  x <- stats::rbinom(N, 1, 0.5)
  d <- wiener_gng_simulate(N, mu = 0.4 + 1.2 * x, bs = 1.4, ndt = 0.25,
                           bias = 0.5, deadline = 1.5)
  d$x <- x
  fit <- frm(bf(rt | dec(responded) ~ x, bias = 0.5),
             family = wiener_gng(deadline = 1.5), data = d)
  e <- fixef(fit)
  expect_equal(e$mu[["(Intercept)"]], 0.4, tolerance = 0.2)
  expect_equal(e$mu[["x"]], 1.2, tolerance = 0.3)
  # the covariate really does change the go rate, which is the thing a
  # go/no-go design measures
  expect_gt(mean(d$responded[x == 1]) - mean(d$responded[x == 0]), 0.05)
})

# ------------------------------------------------------- go/no-go (e)

test_that("the two deadline spellings are one model", {
  set.seed(14)
  d <- wiener_gng_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
  f1 <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = 1.5), data = d)
  f2 <- frm(bf(rt | dec(responded) + vreal(deadline) ~ 1, bias = 0.5),
            family = wiener_gng(), data = d)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-10)
  expect_equal(unlist(fixef(f1)), unlist(fixef(f2)), tolerance = 1e-6)
})

test_that("a deadline that varies by row is read per row", {
  set.seed(15)
  n <- 1200L
  td <- rep(c(0.9, 2.0), length.out = n)
  d <- wiener_gng_simulate(n, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = td)
  # the longer deadline really does produce more responses
  expect_gt(mean(d$responded[td == 2.0]), mean(d$responded[td == 0.9]))
  fit <- frm(bf(rt | dec(responded) + vreal(deadline) ~ 1, bias = 0.5),
             family = wiener_gng(), data = d)
  expect_true(is.finite(as.numeric(logLik(fit))))
  e <- fixef(fit)
  expect_equal(e$mu[["(Intercept)"]], 1.0, tolerance = 0.3)
  # and reading one deadline for all rows is a DIFFERENT model, so the
  # per-row spelling is doing something
  f2 <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = 2.0), data = d)
  expect_gt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(f2))), 1)
})

test_that("the non-decision-time bound comes from the go trials alone", {
  set.seed(16)
  d <- wiener_gng_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
  fastest_go <- min(d$rt[d$responded == 1])
  # a no-go row's entry is a placeholder; put an absurdly small one there
  # and the bound must not move
  d$rt[d$responded == 0] <- 0.001
  fit <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
             family = wiener_gng(deadline = 1.5), data = d)
  lk <- stats::family(fit)$links$ndt
  expect_equal(lk$linkinv(0), fastest_go / 2, tolerance = 1e-12)
  # and the placeholder really did not change the answer
  d2 <- d; d2$rt[d2$responded == 0] <- 1.5
  f2 <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = 1.5), data = d2)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(f2)),
               tolerance = 1e-10)
})

test_that("wiener_gng refuses what it cannot read", {
  set.seed(17)
  d <- wiener_gng_simulate(200, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
  expect_error(wiener_gng(deadline = -1), "positive finite number")
  expect_error(wiener_gng(max_ndt = 0), "positive finite number")

  # the indicator is mandatory, and so is a deadline from somewhere
  expect_error(frm(bf(rt ~ 1, bias = 0.5),
                   family = wiener_gng(deadline = 1.5), data = d), "dec")
  expect_error(frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
                   family = wiener_gng(), data = d), "vreal")

  # vint() is the wrong spelling here and says so
  d$up <- d$responded
  expect_error(frm(bf(rt | dec(responded) + vint(up) ~ 1, bias = 0.5),
                   family = wiener_gng(deadline = 1.5), data = d),
               "travels through dec\\(\\)")

  # a response after its own deadline has no probability
  d2 <- d
  d2$rt[which(d2$responded == 1)[1]] <- 1.9
  expect_error(frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
                   family = wiener_gng(deadline = 1.5), data = d2),
               "after their own deadline")

  # a data set with no responses identifies nothing
  d3 <- d; d3$responded <- 0
  expect_error(frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
                   family = wiener_gng(deadline = 1.5), data = d3),
               "needs go trials")
})

test_that("dec() reads a logical the way it reads a 0/1 column", {
  set.seed(18)
  d <- wiener_gng_simulate(300, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
  d$lg <- d$responded == 1
  f1 <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = 1.5), data = d)
  f2 <- frm(bf(rt | dec(lg) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = 1.5), data = d)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-10)
})

test_that("weights scale both families' log likelihoods", {
  set.seed(19)
  d <- rdm_simulate(300, v = c(3, 2, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
  d$w1 <- 1
  f0 <- frm(bf(rt | vint(choice) ~ 1), family = rdm(3), data = d)
  f1 <- frm(bf(rt | vint(choice) + weights(w1) ~ 1), family = rdm(3),
            data = d)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f0)),
               tolerance = 1e-6)

  g <- wiener_gng_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
  g$w1 <- 1
  h0 <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = 1.5), data = g)
  h1 <- frm(bf(rt | dec(responded) + weights(w1) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = 1.5), data = g)
  expect_equal(as.numeric(logLik(h1)), as.numeric(logLik(h0)),
               tolerance = 1e-6)
})

test_that("both families refuse the post-fit methods they have no mean for", {
  set.seed(20)
  d <- rdm_simulate(200, v = c(3, 2, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
  f <- frm(bf(rt | vint(choice) ~ 1), family = rdm(3), data = d)
  # Without a refusing mean these did not error: frmtmb fell back to
  # "the first primary dpar on the response scale is the mean", and
  # returned a DRIFT RATE of 3.51 for data whose response times average
  # 0.36. That is the regression these three pin.
  expect_error(fitted(f), "no mean response time")
  expect_error(predict(f, type = "response"), "no mean response time")
  expect_error(residuals(f, type = "response"), "no mean response time")
  # pearson needs the mean before it needs the variance, so the refusal
  # a user meets is the mean's, not the variance function's
  expect_error(residuals(f, type = "pearson"), "no mean response time")
  expect_error(residuals(f, type = "deviance"), "unit deviance")

  g <- wiener_gng_simulate(300, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
  h <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
           family = wiener_gng(deadline = 1.5), data = g)
  expect_error(fitted(h), "no mean response to report")
  expect_error(predict(h, type = "response"), "no mean response to report")
  expect_error(residuals(h, type = "response"), "no mean response to report")
  expect_error(residuals(h, type = "pearson"), "no mean response to report")
  expect_error(residuals(h, type = "deviance"), "unit deviance")
  # the link scale is always available and is what a user reads instead
  expect_no_error(predict(f, type = "link"))
  expect_no_error(predict(h, type = "link"))

  # The refusal comes from the FAMILY, not from a core guard, and that
  # is the property that matters. Core's own protection is backwards:
  # mean_is_mu() reads a MISSING mean_fn as agreement that the mean is
  # mu, so fitted() refuses only for a family with no dpar called mu.
  # rdm() would be protected by an accident of naming and
  # wiener_gng(), whose drift IS called mu, not at all. Declaring a
  # mean that stops makes both refusals independent of which way that
  # core guard happens to fall, so they survive the core fix unchanged.
  expect_match(tryCatch(fitted(f), error = conditionMessage), "^rdm:")
  expect_match(tryCatch(fitted(h), error = conditionMessage), "^wiener_gng:")
  # rdm has no dpar called mu and wiener_gng does; both refuse anyway
  expect_false("mu" %in% names(stats::family(f)$links))
  expect_true("mu" %in% names(stats::family(h)$links))
})

# =================================== across-trial variability, go/no-go

test_that("the go branch IS wiener()'s density at the same variability", {
  skip_if_not_installed("RWiener")
  # Not a tolerance, an identity. gng_lpdf_var() calls the same
  # ddm_lpdf_var() with up = 1, the same node sets built from the same
  # `nodes` argument, and the same `delta`. Getting the two `delta`s to
  # agree is why wiener_gng() has a gng_family(cfg, ub, delta) builder:
  # at 0.3.0 it floored the decision time at a flat 1e-12 while wiener()
  # used 1e-9 * min(y), and no identity survives two margins.
  set.seed(21)
  td <- 3.0
  d <- ddm_simulate(300L, mu = 1.0, bs = 1.4, ndt = 0.30, bias = 0.5,
                    sv = 0.5, sz = 0.2, st = 0.1)
  d <- d[d$upper == 1 & d$rt < td, ]
  d$go <- 1L
  at <- list(dec = rep(1, nrow(d)))
  # family_finalize() rather than frm(): it is the seam that resolves the
  # bound and the margin, it is what a fit would call, and it does not
  # build a tape. Two full-variability tapes in one process is what the
  # optimizer cannot afford here.
  fin <- function(fam) fam[["family_finalize"]](fam, d$rt, at)
  for (vv in list("sv", "sz", "st", c("sv", "sz"), c("sv", "sz", "st"))) {
    fw <- fin(wiener(variability = vv))
    fg <- fin(wiener_gng(deadline = td, variability = vv))
    dp <- list(mu = 1.1, bs = 1.35, ndt = 0.22, bias = 0.47,
               sv = 0.6, sz = 0.18, st = 0.09)[c("mu", "bs", "ndt",
                                                 "bias", vv)]
    expect_identical(fg[["lpdf"]](d$rt, dp, at),
                     fw[["lpdf"]](d$rt, dp, at))
  }
})

test_that("the go mass and the no-go probability still add to one", {
  td <- 1.5
  dp <- list(mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.5)
  for (vv in list(character(0), "sv", "sz", "st", c("sv", "sz", "st"))) {
    nd <- ddm_nodes(vv, c(sz = 7L, st = 21L))
    ndc <- gng_nodes(vv, c(sv = 31L, sz = 9L, st = 9L))
    sv <- if ("sv" %in% vv) 0.5 else 0
    sz <- if ("sz" %in% vv) 0.2 else 0
    st <- if ("st" %in% vv) 0.1 else 0
    # the support starts at ndt - st / 2, not at ndt: an integral begun
    # at ndt misses mass the model has and reports a defect of 8e-04
    # that belongs to the integration limit
    lo <- dp$ndt - st / 2
    f <- function(y) {
      exp(ddm_lpdf_var(y, dp$mu, dp$bs, dp$bias, dp$ndt, sv, sz, st, 1,
                       nd, "st" %in% vv, 1e-9 * lo))
    }
    brk <- sort(unique(c(lo, lo + 1e-6, 0.4, 0.7, 1.0, td)))
    mass <- 0
    for (i in seq_len(length(brk) - 1L)) {
      mass <- mass + stats::integrate(f, brk[i], brk[i + 1L],
                                      rel.tol = 1e-13,
                                      subdivisions = 4000L)$value
    }
    ng <- exp(ddm_nogo_lprob_var(td, dp$mu, dp$bs, dp$bias, dp$ndt,
                                 sv, sz, st, ndc))
    expect_lt(abs(1 - (mass + ng)), 1e-10)
  }
})

test_that("the variability average is tape-safe and differentiates", {
  ndc <- gng_nodes(c("sv", "sz", "st"), c(sv = 11L, sz = 7L, st = 7L))
  lp <- function(p) {
    ddm_nogo_lprob_var(1.5, p[1], p[2], p[3], p[4], p[5], p[6], p[7], ndc)
  }
  for (p in list(c(1.0, 1.4, 0.50, 0.25, 0.5, 0.20, 0.10),
                 c(-0.5, 2.0, 0.35, 0.30, 0.3, 0.10, 0.05),
                 c(3.0, 2.5, 0.80, 0.20, 1.5, 0.25, 0.08))) {
    tp <- RTMB::MakeTape(lp, p)
    g <- tp$jacobian(p)
    expect_true(all(is.finite(g)))
    skip_if_not_installed("numDeriv")
    fd <- numDeriv::grad(lp, p)
    expect_lt(max(abs(g - fd)) / max(1, max(abs(fd))), 1e-7)
  }
})

test_that("a start point pushed past a boundary is a limit, not a NaN", {
  # wiener() says a wide sz at a biased start can push the uniform range
  # past a boundary, and that the density there is a barrier. True of the
  # density. NOT true of the no-go probability, which takes log1p(-w) and
  # returns NaN above one and a probability ABOVE ONE below zero. A NaN
  # is not a barrier, it is the end of the tape.
  # the NaN arrives with log1p(-w)'s own warning, which is the point
  expect_warning(bad <- ddm_nogo_lprob(1.25, 1.0, 1.4, 1.05), "NaN")
  expect_true(is.nan(bad))
  expect_gt(ddm_nogo_lprob(1.25, 1.0, 1.4, -0.05), 0)
  expect_true(is.finite(ddm_lpdf_lower_sv(1.25, 1.0, 1.4, 1.05, 0.5)))

  # the clamp turns both into the boundary case they are
  ndc <- gng_nodes("sz", c(sv = 11L, sz = 7L, st = 7L))
  v <- ddm_nogo_lprob_var(1.5, 1.0, 1.4, 0.80, 0.25, 0, 0.50, 0, ndc)
  expect_true(is.finite(v))
  expect_lt(v, 0)
  # It clamps the DEVIATION, so it is EXACTLY inert at a deviation of
  # zero: the expression collapses to lo + (-lo). Clamping the value
  # instead cost an ulp on every row and with it two exact identities
  # this family is entitled to.
  for (w in c(0.5, 0.45, 0.8, 0.999)) {
    expect_identical(w + ddm_wclamp_dev(0, w), w)
  }
  # and it bites only where the range leaves the boundaries
  expect_identical(0.5 + ddm_wclamp_dev(0.2, 0.5), 0.7)
  expect_lt(0.8 + ddm_wclamp_dev(0.25, 0.8), 1)
})

test_that("the plain go/no-go family is untouched by any of this", {
  skip_if_not_installed("RWiener")
  # The no-variability path must still be the 0.3.0 path, bit for bit:
  # the node sets are the one-node rules, the widths are zero, and a
  # single node at the middle of a zero-width interval is an evaluation
  # at the middle, exactly.
  set.seed(1)
  d <- wiener_gng_simulate(200L, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
  dp <- lapply(list(mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.45),
               function(v) rep(v, nrow(d)))
  at <- list(dec = d$responded)
  plain <- gng_lpdf(d$rt, dp, at, 1.5)
  ndc <- gng_nodes(character(0), c(sv = 15L, sz = 7L, st = 7L))
  viaq <- ddm_nogo_lprob_var(1.5, dp$mu, dp$bs, dp$bias, dp$ndt,
                             0, 0, 0, ndc)
  nogo <- d$responded == 0
  expect_identical(plain[nogo], viaq[nogo])
})

test_that("wiener_gng refuses variability it cannot read", {
  expect_error(wiener_gng(variability = "sw"), "names the across-trial")
  expect_error(wiener_gng(variability = "sw"), "^wiener_gng\\(\\)")
  expect_error(wiener_gng(nogo_nodes = c(sv = 0)), "node counts")
  expect_error(wiener_gng(nogo_nodes = c(bogus = 7)), "node counts")
  expect_error(wiener_gng(nodes = c(sv = 7)), "node counts")
})

# ============================================== cens() and trunc() seams

test_that("rdm scores a right-censored race as the product of survivals", {
  set.seed(7)
  N <- 400L
  d <- rdm_simulate(N, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
  cut <- stats::quantile(d$rt, 0.75)
  d$cens <- as.integer(d$rt > cut)
  d$obs <- pmin(d$rt, cut)
  d$choice[d$cens == 1] <- 1L
  fit <- frm(bf(obs | vint(choice) + cens(cens) ~ 1), family = rdm(3),
             data = d)

  # the same likelihood written out by hand at the fitted parameters
  fe <- fixef(fit)
  lk <- stats::family(fit)$links
  dp <- list(v1 = exp(fe$v1[[1]]), v2 = exp(fe$v2[[1]]),
             v3 = exp(fe$v3[[1]]), A = exp(fe$A[[1]]), k = exp(fe$k[[1]]),
             ndt = lk$ndt$linkinv(fe$ndt[[1]]))
  vp <- paste0("v", 1:3)
  accs <- rdm_pars(dp, vp)
  hand <- vapply(seq_len(N), function(i) {
    if (d$cens[i] == 1) rdm_lccdf(d$obs[i], dp, vp)
    else lba_race_lpdf(d$obs[i] - dp$ndt, d$choice[i], rdm_law, accs)
  }, numeric(1))
  expect_lt(abs(as.numeric(logLik(fit)) - sum(hand)) /
              abs(as.numeric(logLik(fit))), 1e-12)

  # the winner on a censored row is not read, and is EXACTLY not read
  d2 <- d
  d2$choice[d2$cens == 1] <- 3L
  f2 <- frm(bf(obs | vint(choice) + cens(cens) ~ 1), family = rdm(3),
            data = d2)
  expect_identical(as.numeric(logLik(f2)), as.numeric(logLik(fit)))
})

test_that("rdm takes all four censoring codes and a truncation bound", {
  set.seed(11)
  N <- 300L
  d <- rdm_simulate(N, v = c(3.0, 2.0), A = 0.5, k = 0.5, ndt = 0.2)
  d$y2 <- d$rt + 0.05
  d$code <- 0L
  d$code[1:100] <- 2L
  d$code[101:150] <- 1L
  d$code[151:200] <- -1L
  d$yy <- ifelse(d$code == 2L, pmax(d$rt - 0.05, 0.21), d$rt)
  fit <- frm(bf(yy | vint(choice) + cens(code, y2) ~ 1), family = rdm(2),
             data = d)
  fe <- fixef(fit)
  lk <- stats::family(fit)$links
  dp <- list(v1 = exp(fe$v1[[1]]), v2 = exp(fe$v2[[1]]),
             A = exp(fe$A[[1]]), k = exp(fe$k[[1]]),
             ndt = lk$ndt$linkinv(fe$ndt[[1]]))
  vp <- c("v1", "v2")
  accs <- rdm_pars(dp, vp)
  Fq <- function(q) 1 - exp(rdm_lccdf(q, dp, vp))
  hand <- vapply(seq_len(N), function(i) {
    if (d$code[i] == 0L) {
      lba_race_lpdf(d$yy[i] - dp$ndt, d$choice[i], rdm_law, accs)
    } else if (d$code[i] == 1L) {
      rdm_lccdf(d$yy[i], dp, vp)
    } else if (d$code[i] == -1L) {
      log(Fq(d$yy[i]))
    } else {
      log(Fq(d$y2[i]) - Fq(d$yy[i]))
    }
  }, numeric(1))
  expect_lt(abs(as.numeric(logLik(fit)) - sum(hand)) /
              abs(as.numeric(logLik(fit))), 1e-12)
})

test_that("a right-censored go/no-go trial IS a no-go trial", {
  skip_if_not_installed("RWiener")
  # Not a tolerance either. A trial whose clock stopped before it
  # responded and a no-go trial are the same statement about the same
  # process, so the two spellings are the same arithmetic.
  set.seed(9)
  td <- 1.5
  g <- wiener_gng_simulate(400L, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = td)
  fa <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = td), data = g)
  g$cens <- as.integer(g$responded == 0)
  g$allgo <- 1L
  fb <- frm(bf(rt | dec(allgo) + cens(cens) ~ 1, bias = 0.5),
            family = wiener_gng(deadline = td), data = g)
  expect_identical(as.numeric(logLik(fb)), as.numeric(logLik(fa)))
})

test_that("go/no-go refuses the censoring it cannot mean", {
  skip_if_not_installed("RWiener")
  # The likelihood is a defective density plus a point mass. A window
  # normalizer on the response scale renormalizes the density and says
  # nothing about the mass, so lcdf is deliberately absent and frmtmb
  # refuses the rest by name.
  set.seed(9)
  td <- 1.5
  g <- wiener_gng_simulate(200L, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = td)
  g$allgo <- 1L
  g$cl <- ifelse(g$responded == 0, -1L, 0L)
  expect_error(frm(bf(rt | dec(allgo) + cens(cl) ~ 1, bias = 0.5),
                   family = wiener_gng(deadline = td), data = g),
               "need a family with a CDF")
  expect_error(frm(bf(rt | dec(allgo) + trunc(lb = 0.2) ~ 1, bias = 0.5),
                   family = wiener_gng(deadline = td), data = g),
               "need a family with a CDF")
  expect_null(wiener_gng(deadline = td)[["lcdf"]])
  expect_true(is.function(wiener_gng(deadline = td)[["lccdf"]]))
  expect_true(is.function(rdm(2)[["lcdf"]]))
  expect_true(is.function(rdm(2)[["lccdf"]]))
})

test_that("the image sum's shortened truncation changed nothing", {
  # ddm_cdf_ks went from 12 to 4. The blend gives the small-time route a
  # non-zero weight only below u = 0.197, and at that u the |j| = 2 term
  # is already exp(-2 * 2 * 2.5 / 0.197), which is 9e-23. Measured over
  # this grid every truncation from 2 to 12 is bit-identical.
  gr <- expand.grid(t = c(0.05, 0.2, 0.6, 1.5, 4, 8, 15),
                    v = c(-2, -0.5, 0.5, 1, 2, 5),
                    a = c(0.8, 1.4, 2.5, 4),
                    w = c(0.25, 0.45, 0.5, 0.75, 0.9))
  blend <- function(K) {
    small <- 1 - ddm_lower_cdf_small(gr$t, -gr$v, gr$a, 1 - gr$w, K = K)
    large <- ddm_nogo_large(gr$t, gr$v, gr$a, gr$w)
    u <- ddm_floor(gr$t / (gr$a * gr$a), ddm_u_floor)
    lam <- 0.5 * (1 + tanh((log(u) - log(ddm_cdf_u0)) / ddm_cdf_us))
    (1 - lam) * log(ddm_floor(small, ddm_share_floor)) +
      lam * log(ddm_floor(large, ddm_share_floor))
  }
  base <- blend(12L)
  for (K in c(8L, 6L, 4L, 3L, 2L)) expect_identical(blend(K), base)
  expect_identical(ddm_cdf_ks, 4L)
})
