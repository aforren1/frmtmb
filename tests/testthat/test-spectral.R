# The Whittle likelihood is an ordinary exponential fit on the
# periodogram, so nothing here tests an estimator. What is tested is the
# data preparation, because every one of its details is a correctness
# detail a user gets wrong once and never notices, and the two refusals
# that the response alone can carry.

ar1 <- function(n, phi, seed) {
  set.seed(seed)
  as.numeric(stats::filter(stats::rnorm(n), phi, method = "recursive",
                           init = stats::rnorm(1, 0, 1 / sqrt(1 - phi^2))))
}

test_that("frm_periodogram agrees with spec.pgram where they overlap", {
  y <- ar1(1024, 0.7, 11)
  pg <- frm_periodogram(y)
  sp <- stats::spec.pgram(y, taper = 0, detrend = FALSE, fast = FALSE,
                          plot = FALSE)
  keep <- seq_len(nrow(pg))
  # the scaling and the frequency axis are R's own, so a user comparing
  # against spec.pgram is comparing the same quantity
  expect_equal(pg$freq, sp$freq[keep])
  expect_equal(pg$pgram, sp$spec[keep])
  # and spec.pgram keeps one ordinate this deliberately does not
  expect_equal(length(sp$freq), nrow(pg) + 1L)
})

test_that("frequency zero and the Nyquist ordinate are dropped", {
  for (m in c(8L, 9L, 16L, 17L, 64L)) {
    set.seed(m)
    pg <- frm_periodogram(stats::rnorm(m), fs = 4)
    # floor((m - 1) / 2): index 0 always goes, index m / 2 goes when even
    expect_identical(nrow(pg), (m - 1L) %/% 2L)
    expect_true(all(pg$freq > 0))
    # the Nyquist frequency is fs / 2, and it is not in the grid
    expect_true(all(abs(pg$freq - 2) > 1e-12))
    expect_lt(max(pg$freq), 2)
  }
})

test_that("fs puts the frequencies in Hz and the ordinates per Hz", {
  set.seed(2)
  w <- stats::rnorm(4096, 0, 3)
  for (fs in c(1, 256)) {
    pg <- frm_periodogram(w, fs = fs)
    expect_equal(pg$freq, seq_len(nrow(pg)) * fs / 4096)
    # E I = var / fs, which is what makes the model a density per Hz
    expect_equal(mean(pg$pgram), stats::var(w) / fs, tolerance = 0.02)
  }
})

test_that("mean removal moves nothing that is returned, linear does", {
  y <- ar1(512, 0.5, 3)
  a <- frm_periodogram(y, detrend = "mean")
  b <- frm_periodogram(y + 1000, detrend = "none")
  # not a coincidence to document away: the kept frequencies are exactly
  # the ones where a constant contributes nothing, so this also fails if
  # the zero ordinate ever comes back
  expect_equal(a$pgram, b$pgram)
  ramp <- frm_periodogram(y + 0.05 * seq_along(y), detrend = "linear")
  expect_equal(ramp$pgram, a$pgram, tolerance = 0.05)
  expect_false(isTRUE(all.equal(
    frm_periodogram(y + 0.05 * seq_along(y), detrend = "mean")$pgram,
    a$pgram, tolerance = 0.05)))
})

test_that("a taper keeps E I and is the window spec.taper spells", {
  set.seed(4)
  w <- stats::rnorm(4096)
  none <- frm_periodogram(w)
  for (tp in c("hann", "split_cosine")) {
    pg <- frm_periodogram(w, taper = tp)
    # dividing by sum(h^2) rather than m is what keeps the mean model
    expect_equal(mean(pg$pgram), stats::var(w), tolerance = 0.03)
    expect_false(isTRUE(all.equal(pg$pgram, none$pgram)))
  }
  # the split cosine bell is R's own, at R's own proportion
  expect_equal(frmtmb:::taper_weights(100, "split_cosine", 0.1),
               stats::spec.taper(rep(1, 100), p = 0.1))
  expect_identical(frmtmb:::taper_weights(8, "none", 0.1), rep(1, 8))
  expect_equal(frmtmb:::taper_weights(8, "hann", 0.1)[1L], 0)
})

test_that("segment averaging trades resolution for dispersion", {
  set.seed(5)
  w <- stats::rnorm(4096)
  for (k in c(1L, 4L, 8L)) {
    pg <- frm_periodogram(w, segments = k)
    expect_identical(nrow(pg), (4096L %/% k - 1L) %/% 2L)
    # the average of k independent ordinates has CV 1 / sqrt(k)
    expect_equal(stats::sd(pg$pgram) / mean(pg$pgram), 1 / sqrt(k),
                 tolerance = 0.1)
  }
})

test_that("many series come back from one call", {
  y <- ar1(512, 0.5, 6)
  m <- frm_periodogram(cbind(a = y, b = rev(y)))
  expect_identical(names(m), c("freq", "pgram", "series"))
  expect_identical(levels(m$series), c("a", "b"))
  expect_identical(nrow(m), 2L * 255L)
  # a single series stays two columns: the shape follows the INPUT, so
  # it can be relied on without looking at the data
  expect_identical(names(frm_periodogram(y)), c("freq", "pgram"))
  # groups may have different lengths
  g <- rep(c("s1", "s2"), times = c(300, 212))
  gm <- frm_periodogram(y, group = g)
  expect_identical(nrow(gm), 149L + 105L)
  expect_identical(levels(gm$series), c("s1", "s2"))
  expect_error(frm_periodogram(cbind(y, y), group = g), "columns are already")
  g[3] <- NA
  expect_error(frm_periodogram(y, group = g), "has NA")
  expect_error(frm_periodogram(data.frame(a = factor(letters[1:8]))),
               "not numeric")
})

test_that("a ts brings its own sampling rate", {
  y <- stats::ts(ar1(512, 0.3, 7), frequency = 128)
  expect_equal(frm_periodogram(y)$freq, frm_periodogram(as.numeric(y),
                                                        fs = 128)$freq)
})

test_that("frm_periodogram refuses what it cannot transform", {
  expect_error(frm_periodogram(c(1, NA, 3, 4, 5)), "numeric with no NA")
  expect_error(frm_periodogram(stats::rnorm(2)), "usable ordinate")
  expect_error(frm_periodogram(stats::rnorm(64), segments = 40),
               "usable ordinate")
  expect_error(frm_periodogram(stats::rnorm(64), fs = -1), "positive number")
  expect_error(frm_periodogram(stats::rnorm(64), p = 0.9), "in \\(0, 0.5\\)")
  expect_error(frm_periodogram(stats::rnorm(64), segments = 0),
               "positive integer")
})

test_that("a whittle fit recovers an AR(1) that arima() agrees with", {
  y <- ar1(512, 0.7, 8)
  pg <- frm_periodogram(y)
  pg$w <- 2 * pi * pg$freq
  # tanh, not a bare coefficient: the spectrum cannot tell phi from
  # 1 / phi (the two log spectra differ by a constant to 1e-15), so an
  # unconstrained fit is free to return the reciprocal root
  fit <- frm(bf(pgram ~ ls - log(1 - 2 * tanh(z) * cos(w) + tanh(z)^2),
                ls ~ 1, z ~ 1, nl = TRUE),
             family = whittle(), data = pg)
  ref <- stats::arima(y, order = c(1, 0, 0), method = "ML")
  expect_equal(tanh(fixef(fit)[["z"]][[1L]]), unname(ref$coef[["ar1"]]),
               tolerance = 0.02)
  se <- sqrt(diag(vcov(fit)))
  expect_equal(unname(se[[2L]]), sqrt(ref$var.coef[1L, 1L]), tolerance = 0.05)
  # the innovation variance is the exponentiated intercept
  expect_equal(exp(fixef(fit)[["ls"]][[1L]]), ref$sigma2, tolerance = 0.05)
})

test_that("whittle(tapers = k) fixes the shape rather than estimating it", {
  y <- ar1(2048, 0.5, 9)
  pg <- frm_periodogram(y, segments = 4)
  fit <- frm(bf(pgram ~ freq), family = whittle(tapers = 4), data = pg)
  free <- frm(bf(pgram ~ freq), family = Gamma(link = "log"), data = pg)
  expect_equal(fixef(fit)[["shape"]][[1L]], log(4))
  expect_identical(as.character(fit$frame[["map"]][["betad"]]), NA_character_)
  expect_identical(length(fit$obj$par) + 1L, length(free$obj$par))
  expect_equal(unname(predict(fit, type = "disp")[1L]), 4)
  # the mean model is a Gamma GLM with the dispersion held at 1 / k, so
  # glm() with the same link is the reference for the coefficients
  ref <- stats::glm(pgram ~ freq, family = stats::Gamma(link = "log"),
                    data = pg)
  expect_equal(unname(fixef(fit)[["mu"]]), unname(stats::coef(ref)),
               tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit)),
               sum(stats::dgamma(pg$pgram, shape = 4,
                                 scale = stats::fitted(ref) / 4, log = TRUE)),
               tolerance = 1e-4)
  expect_identical(whittle()[["dpars"]], "mu")
  expect_error(whittle(tapers = 2.5), "positive whole number")
  # the family's value is a default: an explicit one in bf() still wins
  user <- frm(bf(pgram ~ freq, shape = 2), family = whittle(tapers = 4),
              data = pg)
  expect_equal(fixef(user)[["shape"]][[1L]], log(2))
})

test_that("whittle refuses a response that cannot be a periodogram", {
  # long enough that the dispersion check applies: it needs 200 ordinates
  # before the null is narrow enough to say anything
  y <- ar1(8192, 0.5, 10)
  pg <- frm_periodogram(y)
  neg <- pg
  neg$pgram <- log(pg$pgram)          # power in dB, the usual mistake
  expect_error(frm(bf(pgram ~ freq), family = whittle(), data = neg),
               "strictly positive")
  # an already-averaged periodogram declared raw. var(diff(log I)) is
  # 2 trigamma(k) whatever the spectrum is, and only a SMALLER k can
  # push it below the threshold
  sm <- frm_periodogram(y, segments = 8)
  expect_error(frm(bf(pgram ~ freq), family = whittle(), data = sm),
               "already averaged")
  expect_silent(frm(bf(pgram ~ freq), family = whittle(tapers = 8),
                    data = sm))
  # and the raw one is not refused, at a spectrum that is far from flat
  expect_silent(frm(bf(pgram ~ freq), family = whittle(), data = pg))
})

test_that("whittle refuses addition terms an ordinate has no meaning for", {
  pg <- frm_periodogram(ar1(256, 0.4, 12))
  pg$cn <- rep("none", nrow(pg))
  expect_error(frm(bf(pgram | cens(cn) ~ freq), family = whittle(), data = pg),
               "cens")
})

test_that("simulate() is refused by name and names the way out", {
  pg <- frm_periodogram(ar1(512, 0.4, 13))
  fit <- frm(bf(pgram ~ freq), family = whittle(), data = pg)
  err <- tryCatch(stats::simulate(fit, nsim = 1), error = conditionMessage)
  expect_match(err, "whittle")
  expect_match(err, "ORDINATES")
  expect_match(err, "frm_series_draw", fixed = TRUE)
  # the entry points built on the same draw refuse for the same reason
  expect_error(frm_bootstrap(fit, nsim = 2), "ORDINATES")
})

test_that("frm_series_draw inverts the transform the fit is built on", {
  y <- ar1(1024, 0.6, 14)
  pg <- frm_periodogram(y, fs = 128)
  fit <- frm(bf(pgram ~ s(freq, k = 8)), family = whittle(), data = pg)
  z <- frm_series_draw(fit, nsim = 3, seed = 1)
  # odd length on purpose: no Nyquist ordinate has to be invented
  expect_identical(dim(z), c(2L * nrow(pg) + 1L, 3L))
  expect_equal(stats::frequency(z), pg$freq[1L] * (2 * nrow(pg) + 1))
  expect_equal(mean(z[, 1L]), 0, tolerance = 1e-8)
  # the round trip: the drawn series lands back on the fit's own grid
  expect_equal(frm_periodogram(z[, 1L], fs = stats::frequency(z))$freq,
               pg$freq)
  # and its expected periodogram IS the fitted spectrum
  set.seed(2)
  many <- frm_series_draw(fit, nsim = 200)
  rat <- vapply(seq_len(200), function(i) {
    mean(frm_periodogram(many[, i],
                         fs = stats::frequency(many))$pgram / fitted(fit))
  }, 0)
  expect_equal(mean(rat), 1, tolerance = 0.02)
})

test_that("frm_series_draw refuses what has no single spectrum", {
  y <- ar1(512, 0.5, 15)
  st <- frm_periodogram(cbind(a = y, b = rev(y)))
  fit <- frm(bf(pgram ~ freq + (1 | series)), family = whittle(), data = st)
  expect_error(frm_series_draw(fit), "one evenly spaced grid")
  one <- st[st$series == "a", ]
  expect_silent(frm_series_draw(fit, newdata = one, nsim = 1))
  pg <- frm_periodogram(y)
  gfit <- frm(bf(pgram ~ freq), family = Gamma(), data = pg)
  expect_error(frm_series_draw(gfit), "needs a whittle")
  # a model written in a transformed frequency keeps no freq column, and
  # the refusal says what to pass instead
  pg$w <- 2 * pi * pg$freq
  wfit <- frm(bf(pgram ~ w), family = whittle(), data = pg)
  expect_error(frm_series_draw(wfit), "pass the frequencies")
  expect_silent(frm_series_draw(wfit, freq = pg$freq))
})

test_that("leakage-dominated ordinates are refused too", {
  # An untapered periodogram of a spectrum steeper than f^-2 is leakage,
  # not signal: the ordinates stop being independent exponentials and
  # the same dispersion check sees it. Measured in dev/freq-findings.md:
  # refused 97% of the time at exponent 3, 0% at exponent 1.
  set.seed(21)
  fs <- 256
  nf <- 4088
  fl <- seq_len(nf) * fs / (2 * nf + 1)
  amp <- sqrt((2 * nf + 1) * fs * 10 * fl^(-3) / 2)
  z <- complex(real = stats::rnorm(nf, 0, amp),
               imaginary = stats::rnorm(nf, 0, amp))
  # cut a segment out, so it is not periodic and does leak
  x <- (Re(stats::fft(c(0 + 0i, z, Conj(rev(z))), inverse = TRUE)) /
          (2 * nf + 1))[1:1023]

  raw <- frm_periodogram(x, fs = fs)
  raw$logf <- log(raw$freq)
  err <- tryCatch(frm(bf(pgram ~ logf), family = whittle(), data = raw),
                  error = conditionMessage)
  expect_type(err, "character")
  expect_match(err, "leakage")
  expect_match(err, "taper", fixed = TRUE)

  # the taper is the remedy, and then the exponent comes back
  hann <- frm_periodogram(x, fs = fs, taper = "hann")
  hann$logf <- log(hann$freq)
  fit <- frm(bf(pgram ~ logf), family = whittle(), data = hann)
  expect_equal(-fixef(fit)[["mu"]][[2L]], 3, tolerance = 0.15)
})

test_that("the smoothness check reaches short responses", {
  # It used to return early below 200 ordinates, which made it inert at
  # the length this family is aimed at: a one-second epoch at 256 Hz is
  # 127 ordinates. The threshold now shrinks with the row count instead.
  set.seed(31)
  # the review's own reproduction line, which used to be accepted
  sm <- frm_periodogram(stats::rnorm(1024), segments = 4)
  expect_identical(nrow(sm), 127L)
  expect_error(frm(bf(pgram ~ freq), family = whittle(), data = sm),
               "already averaged")
  expect_silent(frm(bf(pgram ~ freq), family = whittle(tapers = 4),
                    data = sm))

  # and it still does not fire on a legitimate raw periodogram of the
  # same length
  expect_silent(frm(bf(pgram ~ freq), family = whittle(),
                    data = frm_periodogram(stats::rnorm(256))))

  # below 24 ordinates there is nothing to test with, so it does not run
  tiny <- frm_periodogram(stats::rnorm(1024), segments = 64)
  expect_lt(nrow(tiny), 24L)
  expect_silent(whittle()[["valid_y"]](tiny$pgram, NULL))
})

test_that("the smoothness threshold shrinks with the row count", {
  f <- frmtmb:::whittle_smooth_frac
  expect_lt(f(24), f(64))
  expect_lt(f(64), f(256))
  # flat outside the knots, and never above the half the old rule used
  expect_equal(f(512), 0.5)
  expect_equal(f(50000), 0.5)
  expect_equal(f(1), f(24))
  expect_true(all(f(c(24, 32, 48, 64, 96, 128, 192, 256)) < 0.5))
})

test_that("the smoothness check needs frequency order, and says so", {
  # the one way to make it refuse something legitimate: a MONOTONE
  # reordering removes the step-to-step variation the statistic reads
  set.seed(32)
  raw <- frm_periodogram(stats::rnorm(4096))$pgram
  vy <- whittle()[["valid_y"]]
  expect_silent(vy(raw, NULL))
  expect_silent(vy(sample(raw), NULL))          # shuffling only adds
  expect_error(vy(sort(raw), NULL), "too smooth")
})
