#' The cross-spectrum of a pair of signals, as rows a model can read
#'
#' Turns two time series into one row per frequency holding the
#' Hermitian cross-spectral matrix and the number of independent complex
#' draws that went into it. Those rows are the response
#' [cross_wishart()] consumes.
#'
#' @param x,y The two signals. Either two numeric vectors of equal
#'   length, or two matrices with the same dimensions whose COLUMNS are
#'   UNITS (subjects, trials, sessions) and never channels; a matrix pair
#'   adds an `id` column naming the column. A three-column pair is three
#'   units of the same two signals, not a three-channel recording:
#'   [cross_wishart()] models a channel PAIR, and a third channel has no
#'   route in. To relate three channels, fit the three pairs separately
#'   and say so.
#' @param sfreq Sampling rate in Hz. `freq` comes back in Hz; the
#'   default of 1 returns cycles per sample.
#' @param segments Number of disjoint, non-overlapping blocks the record
#'   is cut into. Each contributes one complex draw per frequency.
#' @param tapers Number of orthogonal sine tapers applied within each
#'   segment. Each contributes one further draw. `1`, the default, means
#'   no taper at all rather than the first sine taper.
#' @param smooth Number of adjacent Fourier bins averaged together. The
#'   returned `freq` is the middle of each group. `tapers` and `smooth`
#'   may not both be above 1; see "How many degrees of freedom".
#' @param frange Optional `c(low, high)` in the same units as `freq`,
#'   applied after everything else.
#'
#' @return A data frame with one row per retained frequency and columns
#'   `freq`, `w11`, `w22`, `w12r`, `w12i` and `n`, plus `id` for the
#'   matrix form. `w11` and `w22` are the two auto-spectra summed over
#'   draws, `w12r` and `w12i` the real and imaginary parts of the summed
#'   cross-spectrum, and `n` the degrees of freedom. The matrix is
#'   `[[w11, w12r + 1i w12i], [w12r - 1i w12i, w22]]`, and it is the SUM
#'   over draws rather than the average, which is the scaling the
#'   complex Wishart density is written for. Each draw carries the `1/N`
#'   periodogram convention, so `w11 / n` estimates a spectral density
#'   rather than a squared transform. Any constant rescaling of the pair
#'   is absorbed by the two fitted powers and leaves coherence and phase
#'   untouched, so the convention matters for reading `mu` and `pow2`
#'   and for nothing else.
#'
#' @section What is dropped, and why:
#' Frequency zero and the Nyquist frequency are removed. Their Fourier
#' coefficients are real rather than complex, so the pair's periodogram
#' there is real Wishart with half the degrees of freedom and a
#' different density. Keeping them is a mistake a user makes once and
#' never notices, because the fit still converges.
#'
#' The mean is removed from every segment before its transform, which is
#' what makes frequency zero uninformative in the first place.
#'
#' @section How many degrees of freedom, and where to get them:
#' The density needs `n >= 2`, because a cross-periodogram from a single
#' draw is rank one: its coherence is exactly 1 by arithmetic, whatever
#' the signals did. Below 4 the estimate is nearly worthless even though
#' it exists.
#'
#' Segments, tapers and smoothing all buy real degrees of freedom. At a
#' true coherence of zero the naive estimate has mean exactly `1/n`, so
#' `1 / mean(coherence)` reads the effective count straight off a
#' simulation. Every accepted configuration was measured that way, 4000
#' replicates on white noise, nominal `n = 8` throughout:
#'
#' | configuration | measured n |
#' | --- | --- |
#' | `segments = 8` | 7.92 |
#' | `segments = 4, smooth = 2` | 8.17 |
#' | `segments = 2, smooth = 4` | 7.99 |
#' | `segments = 1, smooth = 8` | 7.98 |
#' | `segments = 4, tapers = 2` | 8.03 |
#' | `segments = 2, tapers = 4` | 7.76 |
#' | `segments = 1, tapers = 8` | 7.98 |
#'
#' Two combinations are refused rather than counted wrong.
#'
#' **Tapers and smoothing together.** Both widen the same spectral
#' window, so their product is not the degrees of freedom: 4 tapers
#' smoothed over 4 bins measures 5.80 against a nominal 16, and 2 tapers
#' over 2 bins measures 2.70 against 4. Raise `segments` instead.
#'
#' **Overlapping segments.** Overlap raises the nominal count without
#' raising the independent one, so `n` would be a lie and every standard
#' error downstream would be too small.
#'
#' Frequency smoothing assumes the spectrum is flat across the bins it
#' averages, so it is the route with a cost that a white-noise
#' measurement cannot see. Measured on an AR(1) spectrum with
#' `phi = 0.9`, series length 4096, every retained bin of 250 replicates
#' pooled, **that cost is not detectable at these widths**: nominal 2,
#' 4, 8 and 16 deliver 1.994, 3.993, 8.008 and 16.069, against 2.002,
#' 3.992, 8.011 and 15.965 on white noise. It appears only at
#' `smooth = 32`, which delivers 31.5.
#'
#' An earlier draft of this package reported a 5 percent loss at
#' `smooth = 16` on that spectrum. That number came from reading one
#' frequency bin per replicate instead of all of them, and it did not
#' survive replication. Smoothing is safe at the widths anyone uses;
#' what is not safe is smoothing across a band where the COHERENCE
#' varies, which is a modelling error rather than a degrees-of-freedom
#' one and no count will catch it.
#'
#' @examples
#' set.seed(1)
#' src <- rnorm(2048)
#' a <- src + rnorm(2048)
#' b <- 0.8 * src + rnorm(2048)
#' xs <- frm_cross_spectrum(a, b, sfreq = 256, segments = 8)
#' head(xs)
#' @export
frm_cross_spectrum <- function(x, y, sfreq = 1, segments = 8L, tapers = 1L,
                               smooth = 1L, frange = NULL) {
  segments <- cp_count(segments, "segments")
  tapers <- cp_count(tapers, "tapers")
  smooth <- cp_count(smooth, "smooth")
  if (is.matrix(x) || is.matrix(y)) {
    if (!is.matrix(x) || !is.matrix(y) || !identical(dim(x), dim(y))) {
      stop("`x` and `y` must both be matrices of the same dimensions ",
           "when either is a matrix; their columns are the units.",
           call. = FALSE)
    }
    ids <- colnames(x)
    if (is.null(ids)) ids <- as.character(seq_len(ncol(x)))
    out <- lapply(seq_len(ncol(x)), function(j) {
      d <- frm_cross_spectrum(x[, j], y[, j], sfreq = sfreq,
                              segments = segments, tapers = tapers,
                              smooth = smooth, frange = frange)
      cbind(id = factor(ids[j], levels = ids), d)
    })
    return(do.call(rbind, out))
  }
  x <- cp_series(x, "x")
  y <- cp_series(y, "y")
  if (length(x) != length(y)) {
    stop("`x` and `y` must have the same length; they are ", length(x),
         " and ", length(y), ".", call. = FALSE)
  }
  if (!is.numeric(sfreq) || length(sfreq) != 1L || !is.finite(sfreq) ||
        sfreq <= 0) {
    stop("`sfreq` must be one positive finite number.", call. = FALSE)
  }
  seglen <- length(x) %/% segments
  if (seglen < 4L) {
    stop("`segments` = ", segments, " leaves ", seglen,
         " samples per segment; at least 4 are needed for one usable ",
         "Fourier frequency.", call. = FALSE)
  }
  if (segments * tapers * smooth < 2L) {
    stop("`segments`, `tapers` and `smooth` multiply to 1 degree of ",
         "freedom. A cross-periodogram from one complex draw is rank one: ",
         "its coherence is exactly 1 by arithmetic and its determinant is ",
         "zero, so no model can read it. Raise one of the three.",
         call. = FALSE)
  }
  if (tapers > 1L && smooth > 1L) {
    stop("`tapers` and `smooth` cannot both be above 1. Both widen the ",
         "same spectral window, so the degrees of freedom would not be ",
         "their product: measured at a true coherence of zero over 4000 ",
         "replicates, 4 tapers smoothed over 4 bins buys 5.8 independent ",
         "draws rather than 16. Pick one, and raise `segments` for the ",
         "rest.", call. = FALSE)
  }
  if (tapers >= seglen) {
    stop("`tapers` = ", tapers, " is not fewer than the ", seglen,
         " samples in a segment; sine tapers past that are zero.",
         call. = FALSE)
  }
  ## Retained bins: 1 .. floor((seglen - 1) / 2), which drops bin 0 and,
  ## for even seglen, the Nyquist bin.
  keep <- seq_len((seglen - 1L) %/% 2L)
  H <- cp_taper_basis(seglen, tapers)
  acc11 <- numeric(length(keep)); acc22 <- numeric(length(keep))
  accr <- numeric(length(keep)); acci <- numeric(length(keep))
  for (s in seq_len(segments)) {
    idx <- ((s - 1L) * seglen + 1L):(s * seglen)
    xs <- x[idx] - mean(x[idx]); ys <- y[idx] - mean(y[idx])
    Fx <- stats::mvfft(H * xs)[keep + 1L, , drop = FALSE]
    Fy <- stats::mvfft(H * ys)[keep + 1L, , drop = FALSE]
    acc11 <- acc11 + rowSums(Mod(Fx)^2)
    acc22 <- acc22 + rowSums(Mod(Fy)^2)
    cr <- rowSums(Fx * Conj(Fy))
    accr <- accr + Re(cr); acci <- acci + Im(cr)
  }
  ## The 1/N periodogram convention, so that w11 / n estimates a spectral
  ## density of order 1 for data of order 1 rather than of order N. It
  ## changes no coherence and no phase, which are ratios, and it shifts the
  ## fitted log powers by a constant.
  acc11 <- acc11 / seglen; acc22 <- acc22 / seglen
  accr <- accr / seglen; acci <- acci / seglen
  freq <- keep * sfreq / seglen
  if (smooth > 1L) {
    g <- length(keep) %/% smooth
    if (g < 1L) {
      stop("`smooth` = ", smooth, " is wider than the ", length(keep),
           " retained frequencies.", call. = FALSE)
    }
    grp <- rep(seq_len(g), each = smooth)
    take <- seq_len(g * smooth)
    bin <- function(v) as.numeric(tapply(v[take], grp, sum))
    freq <- as.numeric(tapply(freq[take], grp, mean))
    acc11 <- bin(acc11); acc22 <- bin(acc22)
    accr <- bin(accr); acci <- bin(acci)
  }
  out <- data.frame(freq = freq, w11 = acc11, w22 = acc22,
                    w12r = accr, w12i = acci,
                    n = segments * tapers * smooth)
  if (!is.null(frange)) {
    if (!is.numeric(frange) || length(frange) != 2L || anyNA(frange) ||
          frange[1L] >= frange[2L]) {
      stop("`frange` must be `c(low, high)` with low below high.",
           call. = FALSE)
    }
    out <- out[out$freq >= frange[1L] & out$freq <= frange[2L], ,
               drop = FALSE]
    if (!nrow(out)) {
      stop("`frange` keeps no frequency; the retained band runs from ",
           signif(freq[1L], 4), " to ", signif(freq[length(freq)], 4), ".",
           call. = FALSE)
    }
  }
  rownames(out) <- NULL
  out
}

#' The taper basis for one segment: a boxcar at k = 1, sine tapers above.
#'
#' `tapers = 1` means no taper. That is what a user means by it, and it
#' is also what keeps the degrees of freedom honest: a tapered transform
#' has correlated neighbouring bins, so applying the first sine taper
#' and then smoothing over 4 of them buys 6.1 independent draws rather
#' than 8. Measured, 3000 replicates at a true coherence of zero.
#'
#' Sine tapers rather than Slepian sequences above that: they are a
#' closed-form expression instead of an eigenproblem, they are exactly
#' orthogonal, and section 4 of `dev/xspec-findings.md` measures their
#' effective degrees of freedom at the nominal count. They are scaled by
#' `sqrt(N)` so that the returned spectra are on the same scale as the
#' untapered ones rather than jumping when `tapers` goes from 1 to 2.
#'
#' @noRd
cp_taper_basis <- function(N, tapers) {
  if (tapers == 1L) return(matrix(1, N, 1L))
  vapply(seq_len(tapers),
         function(k) sqrt(2 * N / (N + 1)) * sin(pi * k * seq_len(N) / (N + 1)),
         numeric(N))
}

#' @noRd
cp_count <- function(v, nm) {
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 ||
        v != round(v)) {
    stop("`", nm, "` must be one whole number of at least 1.", call. = FALSE)
  }
  as.integer(v)
}

#' @noRd
cp_series <- function(v, nm) {
  if (!is.numeric(v) || !length(v)) {
    stop("`", nm, "` must be a non-empty numeric vector.", call. = FALSE)
  }
  if (anyNA(v) || any(!is.finite(v))) {
    stop("`", nm, "` holds a missing or non-finite value; a transform ",
         "cannot skip a sample the way a regression skips a row.",
         call. = FALSE)
  }
  as.numeric(v)
}
