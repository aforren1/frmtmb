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
                 list(segments = 1L, tapers = 2L),
                 list(segments = 4L, window = "hann"),
                 list(segments = 2L, window = "hann"))) {
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
  expect_error(frm_cross_spectrum(c(p$x[-1], Inf), p$y),
               "holds an infinite value")
  expect_error(frm_cross_spectrum(rep(NA_real_, 512), p$y),
               "recorded in both")
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 4L,
                                  window = "hann", smooth = 2L),
               "cannot be combined")
  expect_error(frm_cross_spectrum(p$x, p$y, segments = 4L,
                                  window = "hann", tapers = 2L),
               "two tapers of the same segment")
  expect_error(frm_cross_spectrum(p$x, p$y, window = "hamming"),
               "should be one of")
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

## ---- gaps in the record -----------------------------------------------
## An artifact-rejected recording arrives with NA runs in it. The
## property that makes one call correct is that it equals the hand
## assembly a user would otherwise write: the clean spans transformed
## separately and their accumulators added.

cp_gapped <- function() {
  # 1200 + 600 + 600 usable samples, so at segments = 8 the segment
  # length is 300 and every span holds a whole number of them
  n <- 2500L
  src <- stats::rnorm(n)
  x <- src + stats::rnorm(n)
  y <- 0.8 * src + stats::rnorm(n)
  list(x = x, y = y, gap = c(1201:1250, 1851:1900),
       spans = list(1:1200, 1251:1850, 1901:2500))
}

test_that("a record with two rejected spans is the sum of the three", {
  set.seed(11)
  d <- cp_gapped()
  xg <- d$x; yg <- d$y
  xg[d$gap] <- NA
  yg[d$gap] <- NA
  got <- frm_cross_spectrum(xg, yg, segments = 8L)
  # the same three spans, each cut into segments of the same 300 samples
  parts <- Map(function(ix, k) frm_cross_spectrum(d$x[ix], d$y[ix],
                                                  segments = k),
               d$spans, c(4L, 2L, 2L))
  expect_equal(got$freq, parts[[1]]$freq)
  for (col in c("w11", "w22", "w12r", "w12i")) {
    tot <- parts[[1]][[col]] + parts[[2]][[col]] + parts[[3]][[col]]
    expect_equal(got[[col]], tot, info = col)
  }
  expect_equal(got$n[1], sum(vapply(parts, function(z) z$n[1], 0)))
  expect_equal(got$n[1], 8)
})

test_that("no transform crosses a gap", {
  # A step across a gap would leak into every frequency if a segment
  # ever straddled one. The answer must not move at all.
  set.seed(12)
  d <- cp_gapped()
  xg <- d$x; yg <- d$y
  xg[d$gap] <- NA; yg[d$gap] <- NA
  a <- frm_cross_spectrum(xg, yg, segments = 8L)
  shifted <- d$x
  shifted[1251:2500] <- shifted[1251:2500] + 500
  xg2 <- shifted; xg2[d$gap] <- NA
  b <- frm_cross_spectrum(xg2, yg, segments = 8L)
  expect_equal(a$w11, b$w11)
  expect_equal(a$w12r, b$w12r)
})

test_that("a gap in either signal rejects the pair", {
  set.seed(13)
  p <- cp_pair(2048L)
  xg <- p$x; xg[1000:1099] <- NA
  yg <- p$y; yg[1000:1099] <- NA
  # the same 100 samples removed from x alone, from y alone, or from
  # both give one answer: usability is a property of the pair
  a <- frm_cross_spectrum(xg, p$y, segments = 4L)
  b <- frm_cross_spectrum(p$x, yg, segments = 4L)
  ab <- frm_cross_spectrum(xg, yg, segments = 4L)
  expect_equal(a, b)
  expect_equal(a, ab)
})

test_that("n reports the segments the record supplied, not the request", {
  set.seed(14)
  p <- cp_pair(2048L)
  full <- frm_cross_spectrum(p$x, p$y, segments = 8L)
  expect_equal(full$n[1], 8)
  # one span rejected at the head: the usable count falls and the
  # segment length falls with it
  xg <- p$x; xg[1:200] <- NA
  cut <- frm_cross_spectrum(xg, p$y, segments = 8L)
  expect_lte(cut$n[1], 8)
  expect_gte(cut$n[1], 2)
  # and the answer is the hand assembly of the one surviving span
  one <- frm_cross_spectrum(p$x[201:2048], p$y[201:2048], segments = 8L)
  expect_equal(cut$n[1], one$n[1])
  expect_equal(cut$w11, one$w11)
})

test_that("a record with no whole segment left refuses by name", {
  set.seed(15)
  p <- cp_pair(512L)
  xg <- p$x; xg[10:512] <- NA
  expect_error(frm_cross_spectrum(xg, p$y, segments = 4L),
               "recorded in both")
  # eight spans of 30 samples: asked for 2 segments the length is 120,
  # no span holds one, and the refusal is the rank-one one
  keep <- unlist(lapply(0:7, function(k) k * 64L + 1:30))
  xg2 <- rep(NA_real_, 512); xg2[keep] <- p$x[keep]
  expect_error(frm_cross_spectrum(xg2, p$y, segments = 2L),
               "1 degree of freedom")
})

## ---- the Hann window --------------------------------------------------

test_that("a Hann segment is the windowed periodogram, exactly", {
  # stats::spec.pgram() on the windowed segment is the independent
  # reference. The only thing this package adds is the sqrt(N/sum(w^2))
  # scaling that keeps w11 / n a spectral density, and the identity pins
  # that constant along with the rest of the arithmetic.
  set.seed(16)
  p <- cp_pair(512L)
  xs <- frm_cross_spectrum(p$x, p$y, segments = 2L, window = "hann")
  seglen <- 256L
  w <- 0.5 * (1 - cos(2 * pi * (seq_len(seglen) - 1) / seglen))
  cw <- seglen / sum(w^2)
  ref <- numeric(nrow(xs))
  for (s in 1:2) {
    ix <- ((s - 1L) * seglen + 1L):(s * seglen)
    sp <- stats::spec.pgram(w * (p$x[ix] - mean(p$x[ix])), taper = 0,
                            detrend = FALSE, demean = FALSE, fast = FALSE,
                            plot = FALSE)
    ref <- ref + cw * sp$spec[seq_len(nrow(xs))]
  }
  expect_equal(xs$w11, ref)
  # and the power scale does not move when the window is turned on,
  # within the untapered arm's own standard error on the same data
  box <- frm_cross_spectrum(p$x, p$y, segments = 2L)
  rel <- stats::sd(box$w11 / box$n) / sqrt(nrow(box)) /
    mean(box$w11 / box$n)
  expect_lt(abs(log(mean(xs$w11 / xs$n) / mean(box$w11 / box$n))),
            5 * rel)
})

test_that("a Hann window buys no degrees of freedom and costs none", {
  # Segments are disjoint, so a window inside each one leaves the number
  # of independent complex draws alone. Measured the way every other row
  # of the degrees-of-freedom table is measured, against the untapered
  # arm on the same draws rather than against a written-down number.
  skip_on_cran()
  set.seed(17)
  R <- 1200L
  vb <- numeric(R); vh <- numeric(R)
  coh1 <- function(xs) {
    k <- max(1L, nrow(xs) %/% 4L)
    (xs$w12r[k]^2 + xs$w12i[k]^2) / (xs$w11[k] * xs$w22[k])
  }
  for (r in seq_len(R)) {
    a <- stats::rnorm(512); b <- stats::rnorm(512)
    vb[r] <- coh1(frm_cross_spectrum(a, b, segments = 8L))
    vh[r] <- coh1(frm_cross_spectrum(a, b, segments = 8L,
                                     window = "hann"))
  }
  # the delta-method standard error of 1/mean(v), from the run itself
  se <- function(v) stats::sd(v) / sqrt(length(v)) / mean(v)^2
  expect_lt(abs(1 / mean(vh) - 1 / mean(vb)),
            4 * sqrt(se(vb)^2 + se(vh)^2))
  expect_lt(abs(1 / mean(vh) - 8), 4 * se(vh))
})

test_that("a Hann window correlates adjacent frequencies", {
  # This is why it cannot be smoothed over, and the refusal that says so
  # rests on it. Untapered ordinates are independent across frequency
  # and Hann ones are not; three apart they are independent again.
  skip_on_cran()
  set.seed(18)
  R <- 200L
  lag_of <- function(xs, k) {
    stats::acf(log(xs$w11), lag.max = k, plot = FALSE)$acf[k + 1L]
  }
  m <- vapply(seq_len(R), function(r) {
    a <- stats::rnorm(1024); b <- stats::rnorm(1024)
    bx <- frm_cross_spectrum(a, b, segments = 2L)
    hn <- frm_cross_spectrum(a, b, segments = 2L, window = "hann")
    c(none1 = lag_of(bx, 1L), hann1 = lag_of(hn, 1L),
      hann3 = lag_of(hn, 3L))
  }, numeric(3))
  mu <- rowMeans(m)
  se <- apply(m, 1L, function(v) stats::sd(v) / sqrt(length(v)))
  # the untapered arm is indistinguishable from zero, the Hann arm is
  # far from it, and three bins apart the window is gone again. All
  # three in units of the run's own standard error.
  expect_lt(abs(mu[["none1"]]), 4 * se[["none1"]])
  expect_gt(mu[["hann1"]] / se[["hann1"]], 20)
  expect_lt(abs(mu[["hann3"]]), 4 * se[["hann3"]])
})
