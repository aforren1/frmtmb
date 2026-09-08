## What frm_cross_spectrum() has to get right is not the transform, which
## is stats::fft(), but the bookkeeping around it: which bins survive, what
## the degrees of freedom are, and whether the matrix it hands over is
## positive definite. Every one of those is a correctness detail a user
## gets wrong once and never notices, because the fit converges regardless.

cp_pair <- function(n = 4096L, alpha = 0.8) {
  src <- stats::rnorm(n)
  list(x = src, y = alpha * src + sqrt(1 - alpha^2) * stats::rnorm(n),
       coh = alpha^2)
}

test_that("frequency zero and Nyquist are dropped", {
  set.seed(1)
  p <- cp_pair(1024L)
  # segments = 4 gives segments of 256, an even length, so both ends exist
  xs <- frm_cross_spectrum(p$x, p$y, segments = 4L)
  expect_equal(nrow(xs), 127L)
  expect_gt(min(xs$freq), 0)
  expect_lt(max(xs$freq), 0.5)
  # an odd segment length has no Nyquist bin to drop, only zero
  xs2 <- frm_cross_spectrum(p$x[1:1020], p$y[1:1020], segments = 4L)
  expect_equal(nrow(xs2), 127L)
  expect_gt(min(xs2$freq), 0)
})

test_that("sfreq puts the frequencies in Hz", {
  set.seed(2)
  p <- cp_pair(2048L)
  a <- frm_cross_spectrum(p$x, p$y, sfreq = 1, segments = 8L)
  b <- frm_cross_spectrum(p$x, p$y, sfreq = 256, segments = 8L)
  expect_equal(b$freq, a$freq * 256)
  expect_equal(a$w11, b$w11)
})

test_that("the degrees of freedom multiply across the three routes", {
  set.seed(3)
  p <- cp_pair(4096L)
  expect_equal(frm_cross_spectrum(p$x, p$y, segments = 8L)$n[1], 8)
  expect_equal(frm_cross_spectrum(p$x, p$y, segments = 4L, tapers = 3L)$n[1], 12)
  expect_equal(frm_cross_spectrum(p$x, p$y, segments = 4L, smooth = 3L)$n[1], 12)
  # but not both at once: they widen the same spectral window
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 4L, tapers = 2L,
                                  smooth = 3L),
               "cannot both be above 1")
})

test_that("tapers = 1 means no taper", {
  # The first sine taper is not a no-op: it correlates neighbouring bins,
  # and smoothing over them then buys fewer draws than `n` claims. What
  # made this visible was the effective degrees of freedom test below
  # reading 6.1 where it should read 8.
  set.seed(31)
  p <- cp_pair(1024L)
  xs <- frm_cross_spectrum(p$x, p$y, segments = 4L, tapers = 1L)
  manual <- vapply(1:4, function(s) {
    idx <- ((s - 1L) * 256L + 1L):(s * 256L)
    Mod(stats::fft(p$x[idx] - mean(p$x[idx]))[6L])^2
  }, 0)
  # the plain periodogram with the 1/N convention, summed over segments
  expect_equal(xs$w11[5], sum(manual) / 256)
  # and the power scale does not jump when tapers goes from 1 to 2
  t2 <- frm_cross_spectrum(p$x, p$y, segments = 4L, tapers = 2L)
  expect_equal(mean(xs$w11 / xs$n), mean(t2$w11 / t2$n), tolerance = 0.15)
})

test_that("every returned matrix is positive definite", {
  set.seed(4)
  p <- cp_pair(2048L)
  for (a in list(list(segments = 2L), list(segments = 8L, tapers = 2L),
                 list(segments = 2L, smooth = 4L),
                 list(segments = 1L, tapers = 2L))) {
    xs <- do.call(frm_cross_spectrum, c(list(p$x, p$y), a))
    expect_true(all(xs$w11 > 0))
    expect_true(all(xs$w22 > 0))
    expect_true(all(xs$w11 * xs$w22 - xs$w12r^2 - xs$w12i^2 > 0),
                info = paste(names(a), unlist(a), collapse = " "))
  }
})

test_that("the matrix is the SUM over draws, not the average", {
  # The density is written for the sum, so smoothing must add rather than
  # average: a smoothed row is exactly the total of the rows it replaces.
  set.seed(5)
  p <- cp_pair(4096L)
  raw <- frm_cross_spectrum(p$x, p$y, segments = 4L)
  sm <- frm_cross_spectrum(p$x, p$y, segments = 4L, smooth = 4L)
  expect_equal(sm$n[1], 4 * raw$n[1])
  for (k in c(1L, 5L, 20L)) {
    idx <- ((k - 1L) * 4L + 1L):(k * 4L)
    expect_equal(sm$w11[k], sum(raw$w11[idx]))
    expect_equal(sm$w12r[k], sum(raw$w12r[idx]))
    expect_equal(sm$freq[k], mean(raw$freq[idx]))
  }
  # and w11 / n estimates a spectral density of order 1 for data of order 1
  expect_lt(abs(log(mean(raw$w11 / raw$n))), 1.5)
})

test_that("a matrix pair gives one block per column with an id", {
  set.seed(6)
  X <- matrix(stats::rnorm(3 * 1024), 1024, 3,
              dimnames = list(NULL, c("s1", "s2", "s3")))
  Y <- 0.7 * X + matrix(stats::rnorm(3 * 1024), 1024, 3)
  xs <- frm_cross_spectrum(X, Y, segments = 4L)
  expect_true("id" %in% names(xs))
  expect_s3_class(xs$id, "factor")
  expect_equal(levels(xs$id), c("s1", "s2", "s3"))
  expect_equal(nrow(xs), 3L * 127L)
  # and each block matches the vector call on that column
  one <- frm_cross_spectrum(X[, 2], Y[, 2], segments = 4L)
  expect_equal(xs[xs$id == "s2", names(one)], one, ignore_attr = TRUE)
})

test_that("the effective degrees of freedom are the nominal ones", {
  # At a true coherence of zero the naive estimate has mean exactly 1/n,
  # so 1/mean(coherence) reads the real degrees of freedom off a
  # simulation. This is the property the whole `n` bookkeeping rests on.
  skip_on_cran()
  set.seed(7)
  R <- 1500L
  for (a in list(list(segments = 8L), list(segments = 2L, tapers = 4L),
                 list(segments = 2L, smooth = 4L),
                 list(segments = 1L, smooth = 8L),
                 list(segments = 4L, tapers = 2L))) {
    v <- numeric(R)
    for (r in seq_len(R)) {
      xs <- do.call(frm_cross_spectrum,
                    c(list(stats::rnorm(512), stats::rnorm(512)), a))
      k <- max(1L, nrow(xs) %/% 4L)
      v[r] <- (xs$w12r[k]^2 + xs$w12i[k]^2) / (xs$w11[k] * xs$w22[k])
    }
    n_eff <- 1 / mean(v)
    expect_equal(n_eff, 8, tolerance = 0.12,
                 info = paste(names(a), unlist(a), collapse = " "))
  }
})

test_that("bad arguments are refused by name", {
  set.seed(8)
  p <- cp_pair(512L)
  expect_error(frm_cross_spectrum(p$x, p$y[1:10]), "same length")
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 0), "at least 1")
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 1L),
               "1 degree of freedom")
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 2.5), "whole number")
  expect_error(frm_cross_spectrum(p$x, p$y, sfreq = -1), "positive finite")
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 400L),
               "samples per segment")
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 4L, tapers = 200L),
               "sine tapers past that are zero")
  expect_error(frm_cross_spectrum(p$x, p$y, tapers = 2L, smooth = 2L),
               "cannot both be above 1")
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 4L, smooth = 500L),
               "wider than")
  expect_error(frm_cross_spectrum(p$x, p$y, frange = c(1, 0)),
               "low below high")
  expect_error(frm_cross_spectrum(p$x, p$y, frange = c(10, 20)),
               "keeps no frequency")
  expect_error(frm_cross_spectrum(c(p$x[-1], NA), p$y), "missing or non-finite")
  expect_error(frm_cross_spectrum(matrix(p$x, ncol = 2), p$y),
               "both be matrices")
})

test_that("frange trims after everything else", {
  set.seed(9)
  p <- cp_pair(2048L)
  full <- frm_cross_spectrum(p$x, p$y, sfreq = 100, segments = 4L)
  cut <- frm_cross_spectrum(p$x, p$y, sfreq = 100, segments = 4L,
                            frange = c(5, 20))
  expect_true(all(cut$freq >= 5 & cut$freq <= 20))
  expect_equal(cut, full[full$freq >= 5 & full$freq <= 20, ],
               ignore_attr = TRUE)
})
