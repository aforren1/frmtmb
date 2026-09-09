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
#'   and say so. `NA` marks a sample the record does not have; see
#'   "Gaps in the record".
#' @param sfreq Sampling rate in Hz. `freq` comes back in Hz; the
#'   default of 1 returns cycles per sample.
#' @param segments How many disjoint, non-overlapping blocks the record
#'   is cut into, at most. Each contributes one complex draw per
#'   frequency. A record with gaps in it usually supplies fewer, and `n`
#'   says how many it supplied.
#' @param tapers Number of orthogonal sine tapers applied within each
#'   segment. Each contributes one further draw. `1`, the default, means
#'   no taper at all rather than the first sine taper.
#' @param smooth Number of adjacent Fourier bins averaged together. The
#'   returned `freq` is the middle of each group. `tapers` and `smooth`
#'   may not both be above 1; see "How many degrees of freedom".
#' @param window The data window applied to each segment before its
#'   transform. `"none"`, the default, is a boxcar. `"hann"` is the
#'   raised cosine, for a spectrum steep enough that the untapered
#'   transform's high frequencies are leakage from its low ones. It buys
#'   no degrees of freedom and costs none, and it cannot be combined with
#'   `tapers` or `smooth`; see "How many degrees of freedom".
#' @param frange Optional `c(low, high)` in the same units as `freq`,
#'   applied after everything else.
#'
#' @return A data frame with one row per retained frequency and columns
#'   `freq`, `w11`, `w22`, `w12r`, `w12i` and `n`, plus `id` for the
#'   matrix form. `w11` and `w22` are the two auto-spectra summed over
#'   draws, `w12r` and `w12i` the real and imaginary parts of the summed
#'   cross-spectrum, and `n` the degrees of freedom the record actually
#'   supplied. The matrix is
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
#' @section Gaps in the record:
#' `NA` in either signal marks a sample the pair does not have, which is
#' what artifact rejection leaves behind. Such a sample is not
#' interpolated and not skipped: the record is cut at it, and the
#' segments are laid inside the clean spans that remain, so no transform
#' ever crosses a gap. An artifact-rejected recording is then one call
#' rather than one call per surviving span and a hand-written sum.
#'
#' A sample is usable only where BOTH signals have it, because the
#' cross-spectrum is a property of the pair. `NaN` counts as `NA`; an
#' infinity is refused, because it is a value the arithmetic cannot use
#' rather than a value the record is missing.
#'
#' The segment LENGTH comes from the usable sample count divided by
#' `segments`, and each clean span then supplies as many whole segments
#' as fit in it, in record order, up to `segments` in total. The samples
#' left over at the end of each span are dropped, and a span shorter
#' than one segment supplies nothing at all.
#'
#' So a record with ONE clean span always yields the full `segments`,
#' whatever was rejected from its ends, because the length divides down.
#' A record broken into SEVERAL usually yields fewer, since each span
#' drops its own remainder. On such a record, asking for FEWER segments
#' can be refused where more would be accepted: the length comes from
#' the global usable count, so a small `segments` can ask for a segment
#' longer than any surviving span. Three spans of 700, 1043 and 699 are
#' refused at `segments = 2`, which wants 1221 samples in a row, and
#' give 3, 7 and 14 segments at 4, 8 and 16.
#'
#' That is not hidden: `n` is the count the
#' record supplied, it is per-row data the density reads, and every
#' standard error downstream is formed from it. Nothing warns, because
#' on a record with gaps the shortfall is the normal case rather than a
#' mistake.
#'
#' In the matrix form each column is its own record, so a column with
#' more rejected samples than its neighbors gets a shorter segment, its
#' own frequency grid and its own `n`. That is correct and it is why
#' `freq` is a column of the frame rather than an attribute of it: a
#' model reads `freq` per row.
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
#' | `segments = 8, window = "hann"` | 8.02 |
#'
#' Three combinations are refused rather than counted wrong.
#'
#' **Tapers and smoothing together.** Both widen the same spectral
#' window, so their product is not the degrees of freedom: 4 tapers
#' smoothed over 4 bins measures 5.80 against a nominal 16, and 2 tapers
#' over 2 bins measures 2.70 against 4. Raise `segments` instead.
#'
#' **A window and smoothing together**, for the same reason and by the
#' same arithmetic. A Hann window makes each ordinate a weighted sum of
#' its two neighbors and itself, so adjacent bins are correlated and
#' averaging them adds less than it claims. Measured at a true coherence
#' of zero, 3000 replicates a cell, nominal `n = 8` throughout: 4
#' segments smoothed over 2 bins delivers 5.66, 2 segments over 4 bins
#' delivers 4.98 and 1 segment over 8 bins delivers 4.68, against 8.19,
#' 8.02 and 8.05 for the same three without the window.
#'
#' Declaring the shortfall instead of refusing it was considered and
#' rejected. There is no one number to declare: the same window loses a
#' different amount at each smoothing width, and frmtmb's own
#' `whittle()` measured the same thing from the other side, where an
#' honestly declared equivalent degrees of freedom does not rescue a
#' smoothed spectrum from its raw-periodogram check at kernel widths of
#' 7 and above (`dev/reviews/2026-09-08-spectral2.md`). The correlation
#' is a property of the estimate, not of the number attached to it.
#'
#' **A window and tapers together.** Both are a taper on the same
#' segment, and applying one over the other is neither.
#'
#' **Overlapping segments** stay refused. Overlap raises the nominal
#' count without raising the independent one, so `n` would be a lie and
#' every standard error downstream would be too small.
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
#' varies, which is a modeling error rather than a degrees-of-freedom
#' one and no count will catch it.
#'
#' @section What a window does to the rows, which `n` does not say:
#' `n` counts the draws behind ONE row and a window leaves that count
#' alone, which the table above measures. What a window and what
#' `tapers` both change is the relation BETWEEN rows, and `n` says
#' nothing about that.
#'
#' Measured on the log of `w11` from this function, 400 replicates of
#' white noise at 1024 samples, at calls this function accepts. The last
#' column is the variance inflation factor of a mean over the rows,
#' `1 + 2 * sum` of the positive correlations, which is the number a
#' smooth in frequency actually pays. It is summed over twelve lags,
#' which matters for one row: at `tapers = 8` the correlation is still
#' positive past lag 4, and stopping there gives 6.45 instead of 8.02.
#' Every other row has died by lag 3 and is the same either way.
#'
#' | configuration | lag 1 | lag 2 | lag 3 | inflation |
#' | --- | --- | --- | --- | --- |
#' | `segments = 8` | -0.010 | -0.015 | -0.021 | 1.00 |
#' | `segments = 8, window = "hann"` | 0.397 | -0.009 | -0.027 | 1.79 |
#' | `segments = 8, tapers = 2` | 0.525 | 0.104 | -0.035 | 2.24 |
#' | `segments = 4, tapers = 4` | 0.745 | 0.502 | 0.271 | 4.07 |
#' | `segments = 1, tapers = 8` | 0.876 | 0.744 | 0.615 | 8.02 |
#'
#' So a Hann-windowed frame carries about one independent frequency in
#' 1.8, and `tapers` costs MORE rather than less: a `tapers = 4` frame
#' carries about one in four, and eight tapers cost a factor of eight,
#' with ordinates still correlated 0.6 three bins apart. An `s(freq)`
#' smooth fitted on either has fewer effective points than it has rows,
#' so its band is optimistic. By how much is not measured here: the
#' inflation factor says how much information the rows carry, not what
#' a penalized smooth does with it.
#'
#' None of that touches `n`, which is right in every one of these
#' configurations: the per-row degrees of freedom and the correlation
#' between rows are different quantities, and `n` claims only the first.
#'
#' Use a window or a taper where leakage is the larger error, which is a
#' steep spectrum, and not by default.
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
                               smooth = 1L, window = c("none", "hann"),
                               frange = NULL) {
  segments <- cp_count(segments, "segments")
  tapers <- cp_count(tapers, "tapers")
  smooth <- cp_count(smooth, "smooth")
  window <- match.arg(window)
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
                              smooth = smooth, window = window,
                              frange = frange)
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
  if (window != "none" && tapers > 1L) {
    stop("`window` = \"", window, "\" and `tapers` = ", tapers,
         " are two tapers of the same segment, and one applied over the ",
         "other is neither. Sine tapers buy degrees of freedom and a ",
         "window does not; pick the one the record needs.",
         call. = FALSE)
  }
  if (window != "none" && smooth > 1L) {
    stop("`window` = \"", window, "\" and `smooth` = ", smooth,
         " cannot be combined. A window makes each ordinate a weighted ",
         "sum of itself and its neighbors, so adjacent bins are ",
         "correlated and averaging them buys less than `smooth` claims: ",
         "measured at a true coherence of zero over 3000 replicates, a ",
         "Hann window smoothed over 4 bins from 2 segments buys 4.98 ",
         "independent draws rather than 8. Raise `segments` instead.",
         call. = FALSE)
  }
  ## Gaps first, because the segment length is a property of the samples
  ## that survive rather than of the record's nominal length.
  sp <- cp_spans(!is.na(x) & !is.na(y))
  usable <- sum(sp[["len"]])
  if (usable < 4L) {
    stop("`x` and `y` have ", usable, " samples that are recorded in ",
         "both, out of ", length(x), ". At least 4 are needed for one ",
         "usable Fourier frequency, so there is nothing here to ",
         "transform.", call. = FALSE)
  }
  seglen <- usable %/% segments
  gap_note <- if (usable < length(x)) {
    paste0(", from the ", usable, " of ", length(x),
           " samples that are recorded in both signals")
  } else ""
  if (seglen < 4L) {
    stop("`segments` = ", segments, " leaves ", seglen,
         " samples per segment", gap_note,
         "; at least 4 are needed for one usable Fourier frequency.",
         call. = FALSE)
  }
  if (tapers >= seglen) {
    stop("`tapers` = ", tapers, " is not fewer than the ", seglen,
         " samples in a segment; sine tapers past that are zero.",
         call. = FALSE)
  }
  starts <- cp_blocks(sp, seglen, segments)
  nseg <- length(starts)
  if (nseg * tapers * smooth < 2L) {
    stop("this record supplies ", nseg, " segment(s) of ", seglen,
         " samples", gap_note, ", and `segments`, `tapers` and `smooth` ",
         "then multiply to 1 degree of freedom. A cross-periodogram ",
         "from one complex draw is rank one: its coherence is exactly 1 ",
         "by arithmetic and its determinant is zero, so no model can ",
         "read it. Raise one of the three.", call. = FALSE)
  }
  if (tapers > 1L && smooth > 1L) {
    stop("`tapers` and `smooth` cannot both be above 1. Both widen the ",
         "same spectral window, so the degrees of freedom would not be ",
         "their product: measured at a true coherence of zero over 4000 ",
         "replicates, 4 tapers smoothed over 4 bins buys 5.8 independent ",
         "draws rather than 16. Pick one, and raise `segments` for the ",
         "rest.", call. = FALSE)
  }
  ## Retained bins: 1 .. floor((seglen - 1) / 2), which drops bin 0 and,
  ## for even seglen, the Nyquist bin.
  keep <- seq_len((seglen - 1L) %/% 2L)
  H <- cp_taper_basis(seglen, tapers, window)
  acc11 <- numeric(length(keep)); acc22 <- numeric(length(keep))
  accr <- numeric(length(keep)); acci <- numeric(length(keep))
  for (s0 in starts) {
    idx <- s0:(s0 + seglen - 1L)
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
                    n = nseg * tapers * smooth)
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

#' The maximal runs of samples both signals have.
#'
#' Returned as a start and a length rather than as index vectors, so
#' that a long record with few gaps costs a handful of integers rather
#' than a copy of itself.
#'
#' @noRd
cp_spans <- function(ok) {
  r <- rle(ok)
  ends <- cumsum(r[["lengths"]])
  k <- which(r[["values"]])
  list(start = ends[k] - r[["lengths"]][k] + 1L, len = r[["lengths"]][k])
}

#' Where each segment starts, laid inside the clean spans.
#'
#' Whole segments only, from the beginning of each span, in record
#' order, and never more than `cap` of them. The cap is what keeps an
#' ungapped record identical to what it produced before gaps were
#' admitted: `usable %/% segments` rounds down, so the one span can
#' hold more than `segments` whole blocks and would otherwise return
#' more degrees of freedom than the user asked for.
#'
#' @noRd
cp_blocks <- function(sp, seglen, cap) {
  out <- integer(0)
  for (i in seq_along(sp[["start"]])) {
    k <- min(sp[["len"]][i] %/% seglen, cap - length(out))
    if (k > 0L) {
      out <- c(out, sp[["start"]][i] + (seq_len(k) - 1L) * seglen)
    }
    if (length(out) >= cap) break
  }
  as.integer(out)
}

#' The taper basis for one segment: a boxcar at k = 1, sine tapers above.
#'
#' `tapers = 1` means no taper. That is what a user means by it, and it
#' is also what keeps the degrees of freedom honest: a tapered transform
#' has correlated neighboring bins, so applying the first sine taper
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
#' A `window` is the same object at one column and is scaled the same
#' way, by `sqrt(N / sum(w^2))`, which for the periodic Hann window is
#' exactly `sqrt(8/3)`. That is what makes `w11 / n` estimate the same
#' spectral density windowed or not: measured on white noise of unit
#' variance, 400 replicates, the mean is 1.0013 (se 0.0022) with no
#' window and 0.9963 (se 0.0032) with Hann.
#'
#' @noRd
cp_taper_basis <- function(N, tapers, window = "none") {
  if (window == "hann") {
    # The PERIODIC raised cosine, which is the one whose transform has
    # three non-zero taps; the symmetric form's does not, and the
    # neighbor correlation the docs measure is a property of the
    # periodic one.
    w <- 0.5 * (1 - cos(2 * pi * (seq_len(N) - 1) / N))
    return(matrix(w * sqrt(N / sum(w^2)), N, 1L))
  }
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
  if (any(is.infinite(v))) {
    stop("`", nm, "` holds an infinite value. NA marks a sample the ",
         "record does not have, and the transform is cut at it; an ",
         "infinity is a value the arithmetic cannot use and says ",
         "nothing about which samples are trustworthy.", call. = FALSE)
  }
  as.numeric(v)
}
