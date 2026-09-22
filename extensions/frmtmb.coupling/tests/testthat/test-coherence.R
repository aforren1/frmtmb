## The extractors, and the hierarchical fit they are for. This is where the
## package's claim lives: averaging per-subject coherences is biased and
## concatenating subjects is worse, so a hierarchical fit is what is left.

cp_form <- function(rhs = "1") {
  stats::as.formula(paste0("w11 | vreal(w22, w12r, w12i) + vint(n) ~ ", rhs))
}
cp_bf <- function(rhs = "1") {
  frmtmb::bf(cp_form(rhs),
             stats::as.formula(paste0("pow2 ~ ", rhs)),
             stats::as.formula(paste0("coh ~ ", rhs)),
             stats::as.formula(paste0("phase ~ ", rhs)))
}

## One frequency per unit, drawn from a named spectral matrix per unit.
cp_units <- function(eta, phase, n, la = NULL, lb = NULL) {
  N <- length(eta)
  if (is.null(la)) la <- numeric(N)
  if (is.null(lb)) lb <- numeric(N)
  out <- vapply(seq_len(N), function(i) {
    a <- exp(la[i]); b <- exp(lb[i]); cm <- b * exp(eta[i])
    z1 <- complex(real = stats::rnorm(n, 0, sqrt(0.5)),
                  imaginary = stats::rnorm(n, 0, sqrt(0.5)))
    z2 <- complex(real = stats::rnorm(n, 0, sqrt(0.5)),
                  imaginary = stats::rnorm(n, 0, sqrt(0.5)))
    d1 <- a * z1
    d2 <- cm * complex(modulus = 1, argument = -phase[i]) * z1 + b * z2
    cr <- sum(d1 * Conj(d2))
    c(sum(Mod(d1)^2), sum(Mod(d2)^2), Re(cr), Im(cr))
  }, numeric(4))
  data.frame(id = factor(seq_len(N)), w11 = out[1, ], w22 = out[2, ],
             w12r = out[3, ], w12i = out[4, ], n = n)
}

test_that("the coherence interval stays inside the unit interval", {
  # A Wald interval on the coherence scale walks off the end near a
  # boundary, which is exactly where few segments put it. This one is
  # formed on the logit scale, so it cannot.
  set.seed(201)
  # eta = 3 was the old top of this loop, a coherence of 0.998. 16 is a
  # coherence of 1 to ten decimals. The ceiling here is the DRAW routine
  # rather than the family: above about eta = 17 the simulated matrix
  # loses its own determinant to rounding and valid_y() rightly refuses
  # it. The density itself is checked directly to eta = 40 in
  # test-cross-wishart.R, where no draw is needed.
  # The invariant: at every eta the answer is either an interval of
  # POSITIVE width inside (0, 1), or a refusal by name. What must never
  # come back is a 95 percent interval whose lower, estimate and upper
  # are the same number, which is what subtracting from 1 used to give.
  for (eta in c(-3, 0, 3, 16)) {
    d <- cp_units(rep(eta, 15L), rep(0.4, 15L), n = 4L)
    fit <- suppressWarnings(
      frmtmb::frm(cp_bf(), family = cross_wishart(), data = d))
    expect_true(is.finite(as.numeric(stats::logLik(fit))),
                info = paste("logLik at eta", eta))
    co <- try(frm_coherence(fit, newdata = d[1, ]), silent = TRUE)
    if (inherits(co, "try-error")) {
      # only legitimate where the fit has run out of curvature
      expect_gt(eta, 3)
      expect_match(conditionMessage(attr(co, "condition")),
                   "no interval exists")
      next
    }
    expect_true(co$.lower > 0 && co$.upper < 1, info = paste("eta", eta))
    expect_gt(co$.upper - co$.lower, 0)
    expect_true(co$.estimate > co$.lower && co$.estimate < co$.upper,
                info = paste("eta", eta))
    expect_equal(co$.estimate, stats::plogis(co$.eta))
  }
})

test_that("an interval near a coherence of 1 has positive width", {
  # Two signals differing by 1e-5 of noise put the estimate at a
  # coherence of 1 to eight decimals. The interval must still be an
  # interval: a 95 percent interval of width zero is a wrong answer
  # rather than a narrow one.
  set.seed(211)
  src <- stats::rnorm(4096)
  d <- frm_cross_spectrum(src, src + stats::rnorm(4096, 0, 1e-5),
                          segments = 16L)
  # a deliberately hard fit: it sits against the boundary and the
  # optimizer says so, which is not what is under test here
  fit <- suppressWarnings(
    frmtmb::frm(cp_bf(), family = cross_wishart(), data = d))
  co <- try(frm_coherence(fit, newdata = d[1, ]), silent = TRUE)
  if (inherits(co, "try-error")) {
    # the fit carries no curvature there, and that is refused by name
    # rather than returned as a point
    expect_match(conditionMessage(attr(co, "condition")),
                 "no interval exists")
  } else {
    expect_gt(co$.upper - co$.lower, 0)
    expect_gt(co$.se, 0)
    expect_lt(co$.lower, co$.estimate)
    expect_gt(co$.upper, co$.estimate)
  }
  # either way the fit itself is finite, which it was not before the
  # complement moved to the log scale
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})

test_that("a group at a coherence of 1 does not poison the fit", {
  # Eight groups at a moderate coherence and one approaching 1. The
  # extreme group used to drive eta past 34 and return NaN standard
  # errors while valid_y() passed every row.
  set.seed(212)
  N <- 9L
  eta <- c(rep(qlogis(0.3) / 2, N - 1L), 16)
  d <- cp_units(eta, rep(0.3, N), n = 8L)
  fit <- suppressWarnings(
    frmtmb::frm(cp_bf("1 + (1 | id)"), family = cross_wishart(), data = d))
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  co <- try(frm_coherence(fit, re_formula = NA), silent = TRUE)
  if (!inherits(co, "try-error")) {
    expect_true(all(is.finite(co$.estimate)))
    expect_true(all(co$.upper - co$.lower > 0))
  }
})

test_that("the interval width follows the confidence level", {
  set.seed(202)
  d <- cp_units(rep(0.2, 30L), rep(0.5, 30L), n = 8L)
  fit <- frmtmb::frm(cp_bf(), family = cross_wishart(), data = d)
  w <- vapply(c(0.5, 0.8, 0.95, 0.99), function(l) {
    co <- frm_coherence(fit, newdata = d[1, ], level = l)
    co$.upper - co$.lower
  }, 0)
  expect_true(all(diff(w) > 0))
  expect_error(frm_coherence(fit, level = 1), "strictly between 0 and 1")
  expect_error(frm_coherence(fit, level = 0), "strictly between 0 and 1")
})

test_that("both extractors refuse a fit of another family", {
  set.seed(203)
  d <- data.frame(y = stats::rgamma(60, 3, 1), x = stats::rnorm(60))
  g <- frmtmb::frm(frmtmb::bf(y ~ x), family = stats::Gamma(link = "log"),
                   data = d)
  expect_error(frm_coherence(g), "cross_wishart")
  expect_error(frm_phase(g), "cross_wishart")
  expect_error(frm_cross_simulate(g), "cross_wishart")
})

test_that("phase comes back in radians and is recovered", {
  set.seed(204)
  for (ph in c(-2.0, 0, 1.5)) {
    d <- cp_units(rep(0.8, 40L), rep(ph, 40L), n = 8L)
    fit <- frmtmb::frm(cp_bf(), family = cross_wishart(), data = d)
    p <- frm_phase(fit, newdata = d[1, ])
    expect_lt(abs(p$.estimate - ph), 0.15, label = paste("phase", ph))
    expect_true(p$.lower < ph && p$.upper > ph)
  }
})

test_that("re_formula separates the group answer from the population one", {
  set.seed(205)
  N <- 30L
  eta <- stats::rnorm(N, 0.3, 0.6)
  d <- cp_units(eta, stats::rnorm(N, 0.5, 0.3), n = 12L)
  fit <- frmtmb::frm(cp_bf("1 + (1 | id)"), family = cross_wishart(),
                     data = d)
  per <- frm_coherence(fit, re_formula = NULL)
  pop <- frm_coherence(fit, re_formula = NA)
  expect_equal(nrow(per), N)
  expect_equal(length(unique(round(pop$.estimate, 10))), 1L)
  expect_gt(stats::sd(per$.estimate), 0)
  # shrinkage: the per-group spread is narrower than the raw per-group one
  raw <- (d$w12r^2 + d$w12i^2) / (d$w11 * d$w22)
  expect_lt(stats::sd(per$.estimate), stats::sd(raw))
})

test_that("a coherence spectrum is read off a smooth over frequency", {
  # The case a smooth exists for: coherence that varies with frequency,
  # and cannot leave (0, 1) at any point on the smooth because the link
  # is the constraint.
  set.seed(206)
  nf <- 80L
  f <- seq_len(nf) / (2 * nf)
  eta_true <- 2.5 * exp(-((f - 0.15)^2) / (2 * 0.03^2)) - 1.5
  d <- cp_units(eta_true, rep(0.3, nf), n = 10L)
  d$freq <- f
  fit <- frmtmb::frm(cp_bf("s(freq, k = 12)"), family = cross_wishart(),
                     data = d)
  co <- frm_coherence(fit)
  expect_equal(nrow(co), nf)
  expect_true(all(co$.lower > 0 & co$.upper < 1))
  # the fitted peak is near the true one
  expect_lt(abs(f[which.max(co$.estimate)] - f[which.max(eta_true)]), 0.04)
  # and the smooth beats a flat model
  flat <- frmtmb::frm(cp_bf(), family = cross_wishart(), data = d)
  expect_lt(stats::AIC(fit), stats::AIC(flat))
})

test_that("the hierarchical fit beats averaging per-subject coherences", {
  # The package's reason to exist, in one comparison. Subjects share a
  # coherence and a phase spread; the naive estimate averages each
  # subject's own coherence, the pooled one concatenates their segments.
  skip_on_cran()
  set.seed(207)
  N <- 40L; n <- 8L; R <- 60L
  mu_eta <- 0; sd_eta <- 0.4; sd_phi <- 0.8
  truth <- stats::plogis(2 * mu_eta)
  naive <- numeric(R); pooled <- numeric(R); model <- numeric(R)
  for (r in seq_len(R)) {
    eta <- stats::rnorm(N, mu_eta, sd_eta)
    ph <- stats::rnorm(N, 0.8, sd_phi)
    d <- cp_units(eta, ph, n, la = stats::rnorm(N, 0, 0.3),
                  lb = stats::rnorm(N, 0, 0.3))
    naive[r] <- mean((d$w12r^2 + d$w12i^2) / (d$w11 * d$w22))
    pooled[r] <- (sum(d$w12r)^2 + sum(d$w12i)^2) / (sum(d$w11) * sum(d$w22))
    fit <- try(frmtmb::frm(cp_bf("1 + (1 | id)"), family = cross_wishart(),
                           data = d), silent = TRUE)
    model[r] <- if (inherits(fit, "try-error")) NA_real_ else
      frm_coherence(fit, newdata = d[1, ], re_formula = NA)$.estimate
  }
  expect_gt(sum(!is.na(model)), 0.9 * R)
  b_naive <- mean(naive) - truth
  b_pooled <- mean(pooled) - truth
  b_model <- mean(model, na.rm = TRUE) - truth
  # naive is biased UP, pooling is biased DOWN by phase cancellation, and
  # the model is closer to the truth than either
  expect_gt(b_naive, 0.005)
  expect_lt(b_pooled, -0.15)
  expect_lt(abs(b_model), abs(b_naive))
  expect_lt(abs(b_model), abs(b_pooled))
})

test_that("a random effect on every dpar keeps the variance component alive", {
  # The failure that looks most like success: with the two channel powers
  # left as free per-subject parameters they are incidental, they starve
  # the coherence variance component, and it collapses to zero without a
  # warning. Random effects on all four keep it.
  skip_on_cran()
  set.seed(208)
  N <- 40L; n <- 4L; R <- 20L
  collapsed <- 0L
  for (r in seq_len(R)) {
    eta <- stats::rnorm(N, 0, 0.4)
    d <- cp_units(eta, stats::rnorm(N, 0.8, 0.8), n,
                  la = stats::rnorm(N, 0, 0.3), lb = stats::rnorm(N, 0, 0.3))
    fit <- try(frmtmb::frm(cp_bf("1 + (1 | id)"), family = cross_wishart(),
                           data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    # VarCorr returns one variance MATRIX per block, so the standard
    # deviation is the square root of its diagonal
    sds <- unlist(lapply(frmtmb::varcorr_matrices(fit),
                         function(m) sqrt(diag(m))))
    if (!length(sds) || all(is.na(sds)) ||
          max(sds, na.rm = TRUE) < 1e-3) collapsed <- collapsed + 1L
  }
  # with all four dpars carrying the random effect the variance component
  # survives in the large majority of runs
  expect_lt(collapsed, R / 2)
})

## Item 2.6's design, small enough for the ordinary suite: one
## subject-by-condition cell per (id, cond) with `n_rep` replicate rows
## drawn from the same spectral matrix, no frequency axis, and only the
## coherence carrying structure.
##
## The id-by-condition deviations are drawn at sd 1 and SCALED, because
## `rnorm(n, 0, 0)` returns without touching the random-number stream:
## the two arms below have to differ in that one component and in
## nothing else, and drawing at sd 0 would move every Wishart draw
## after it as well.
##
## cp_units() reads HALF the logit coherence, since its coherence is
## `plogis(2 * eta)`, so a contrast of `b_cond` on the logit scale is
## passed as `b_cond / 2`.
##
## `lb` is not free. cp_units() builds the second channel as
## `cm z1 + b z2` with `cm = b exp(eta)`, so its power is
## `b^2 (1 + exp(2 eta))`: leaving `lb` at 0 would tie the second
## channel's POWER to the coherence and give a 2.6-fold power spread
## across these rows, which is the trap `?cross_wishart` warns about and
## which `pow2 ~ 1` would then be the wrong model for. Choosing `lb` to
## cancel that factor holds both powers at 1 and leaves the coherence
## alone, since coherence is invariant to scaling a channel.
cp_cells <- function(seed, n_sub, n_rep, sd_idcond, sd_id = 0.35,
                     b0 = -0.6, b_cond = 0.5, n = 8L) {
  set.seed(seed)
  d <- expand.grid(rep = seq_len(n_rep), cond = factor(c("a", "b")),
                   id = factor(seq_len(n_sub)))
  u_id <- stats::rnorm(n_sub) * sd_id
  u_ic <- stats::rnorm(n_sub * 2L) * sd_idcond
  ic <- as.integer(d$id) + n_sub * (as.integer(d$cond) - 1L)
  lc <- b0 + b_cond * (d$cond == "b") + u_id[as.integer(d$id)] + u_ic[ic]
  w <- cp_units(lc / 2, rep(0.4, nrow(d)), n = n,
                la = numeric(nrow(d)), lb = -0.5 * log1p(exp(lc)))
  cbind(d[, c("id", "cond")],
        w[, c("w11", "w22", "w12r", "w12i", "n")])
}

## The powers and the phase are constants in that design, so an
## intercept is the correct model for them and only the coherence
## carries the ladder.
cp_bf_coh <- function(rhs) {
  frmtmb::bf(cp_form("1"), pow2 ~ 1,
             stats::as.formula(paste0("coh ~ ", rhs)), phase ~ 1)
}

cp_width <- function(rhs, d) {
  fit <- suppressWarnings(
    frmtmb::frm(cp_bf_coh(rhs), family = cross_wishart(), data = d,
                se = TRUE))
  ci <- suppressWarnings(stats::confint(fit))
  unname(ci["coh_condb", "upr"] - ci["coh_condb", "lwr"])
}

test_that("a within-condition contrast needs the within-subject term", {
  # Item 2.6 of dev/extension-gaps-plan.md. `cond` varies INSIDE a
  # subject, so `(1 | id)` shifts both of that subject's conditions
  # together and cancels out of the contrast: `(1 | id:cond)` is the
  # term that carries the contrast's between-subject spread. Drop it
  # and the estimate hardly moves while its interval collapses, which
  # is the failure one fit cannot see.
  #
  # At the realistic design, 40 subjects x 2 conditions x 60
  # frequencies, dev/coh-findings.md measures the interval of
  # `cond + s(freq, by = cond) + (1 | id)` at 0.37 times the correct
  # model's over 148 replicates, covering a truth of 0.5 on 78 of them
  # against the correct model's 141.
  #
  # The two arms. TREATMENT has the subject-by-condition effect; NULL
  # does not, which makes `(1 | id)` correct there too, so a widening
  # that were really "an added term always widens an interval" would
  # show up in both.
  #
  # The statistic is PAIRED, seed by seed, and that is not a
  # refinement. The unpaired version, the mean treatment ratio against
  # the largest null ratio, PASSED on data whose truth had no
  # subject-by-condition effect at all, 1.0889 against 1.0772
  # (dev/coh-absent.R): two blocks of four from one distribution land
  # either way often enough. Paired at one seed the two arms differ
  # only by the scaling of the deviations, so with the effect absent
  # they are the SAME data and the quotient is exactly 1, which fails
  # this assertion rather than flipping a coin. Over 40 seeds of this
  # generator the quotient runs 1.422 to 2.596 with a median of 1.893
  # (dev/coh-paired.R).
  #
  # How often it would fire on correct data, observation first: 0 of
  # those 40 fired, which bounds the rate below 0.072 a seed and below
  # 0.26 for the minimum of four, at 95 percent confidence and assuming
  # no shape. A lognormal fit to the same 40 reads far rarer, about 1
  # in 2,200, and a normal fit to them reads 1 in 90, so treat the
  # small number as the fit's and not the run's (dev/coh-falsealarm2.R).
  #
  # Where it stops discriminating, so that a failure there is read as
  # the limit it is rather than as a defect: with the effect PRESENT
  # but SMALL, `sd_idcond` at 0.15 on seeds 2700:2703, the assertion
  # FAILS at min(treat / null) = 0.999999956, because the correct
  # model's component collapses to zero on two of those four seeds. The
  # floor is between 0.15 and 0.5.
  skip_on_cran()
  # The pairing is load-bearing and costs no fit to check. The two arms
  # are the same data at `sd_idcond = 0` only because the generator
  # writes `rnorm(n) * sd`: `rnorm(n, 0, 0)` returns without drawing at
  # all, and that spelling would leave the arms with different Wishart
  # draws and quietly turn the quotient below into the unpaired
  # statistic that failed open.
  stream_after <- function(sd_ic) {
    cp_cells(2610L, n_sub = 4L, n_rep = 2L, sd_idcond = sd_ic)
    .Random.seed
  }
  expect_identical(stream_after(0), stream_after(0.5))
  ratio <- function(sd_ic, seeds) {
    vapply(seeds, function(s) {
      d <- cp_cells(s, n_sub = 16L, n_rep = 6L, sd_idcond = sd_ic)
      cp_width("cond + (1 | id) + (1 | id:cond)", d) /
        cp_width("cond + (1 | id)", d)
    }, numeric(1))
  }
  seeds <- 2610:2613
  treat <- ratio(0.5, seeds)
  null <- ratio(0, seeds)
  expect_gt(min(treat / null), 1)
  # A sign check, and only that: adding the term must never NARROW the
  # interval. It does not discriminate, and dev/coh-absent.R shows it
  # passing with the effect absent, so it is here for the direction and
  # the line above is the one that carries the claim.
  expect_gt(min(treat), 1)
})

test_that("frm_cross_simulate draws matrices with the fitted first moment", {
  set.seed(209)
  d <- cp_units(rep(0.5, 25L), rep(0.7, 25L), n = 12L)
  fit <- frmtmb::frm(cp_bf(), family = cross_wishart(), data = d)
  sims <- frm_cross_simulate(fit, nsim = 3L, seed = 2, newdata = d[1, ])
  expect_length(sims, 3L)
  expect_named(sims[[1]], c("w11", "w22", "w12r", "w12i", "n"))
  expect_equal(sims[[1]]$n, d$n[1])
  # every draw is positive definite, which is the property that makes it
  # refittable rather than merely plausible
  for (s in sims) {
    expect_true(s$w11 > 0 && s$w22 > 0)
    expect_true(s$w11 * s$w22 - s$w12r^2 - s$w12i^2 > 0)
  }
  # and the same seed gives the same draw
  expect_equal(frm_cross_simulate(fit, nsim = 1L, seed = 2, newdata = d[1, ]),
               sims[1])
})

test_that("the simulated first moment matches the fitted matrix", {
  skip_on_cran()
  set.seed(210)
  d <- cp_units(rep(0.5, 25L), rep(0.7, 25L), n = 16L)
  fit <- frmtmb::frm(cp_bf(), family = cross_wishart(), data = d)
  sims <- frm_cross_simulate(fit, nsim = 4000L, seed = 5, newdata = d[1, ])
  m <- rowMeans(vapply(sims, function(z)
    c(z$w11, z$w22, z$w12r, z$w12i), numeric(4)))
  e <- vapply(frmtmb::fixef_by_dpar(fit), function(z) z[["(Intercept)"]], 0)
  s11 <- exp(e[["mu"]]); s22 <- exp(e[["pow2"]])
  ch <- stats::plogis(e[["coh"]]); ph <- e[["phase"]]
  tgt <- 16 * c(s11, s22, sqrt(ch * s11 * s22) * cos(ph),
                sqrt(ch * s11 * s22) * sin(ph))
  expect_equal(m, tgt, tolerance = 0.05, ignore_attr = TRUE)
})
