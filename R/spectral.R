# Frequency-domain models. There is no new likelihood here, and that is
# the point: a periodogram ordinate of a stationary Gaussian series is
# exponential about the spectral density, and frmtmb's exponential() is
# mean-parameterized with a log link, so the Whittle likelihood IS an
# ordinary frm() fit whose response is the periodogram and whose linear
# predictor is the log spectrum. What the user cannot get right by hand
# is the data preparation (which ordinates are exponential and which are
# not) and the assumption behind a smoothed periodogram (which changes
# the shape, not the mean). Those two are what this file supplies.

#' Taper weights, as `stats::spec.taper()` and `stats::spec.pgram()`
#' spell them, so a user comparing against R's own spectral estimate is
#' comparing the same window.
#'
#' @noRd
taper_weights <- function(m, taper, p) {
  switch(taper,
    none = rep(1, m),
    # the periodic Hann window: h[1] is exactly zero, which is why a
    # tapered fit throws away a little of each end of the series
    hann = 0.5 - 0.5 * cos(2 * pi * seq.int(0L, m - 1L) / m),
    split_cosine = stats::spec.taper(rep(1, m), p = p),
    stop("Unknown taper: ", taper, call. = FALSE)
  )
}

#' Remove what the model does not describe. Only `"linear"` changes a
#' returned ordinate: the retained frequencies are exactly the ones at
#' which the DFT of a constant vanishes, so mean removal moves the zero
#' ordinate and nothing else.
#'
#' @noRd
detrend_series <- function(y, detrend) {
  switch(detrend,
    none = y,
    mean = y - mean(y),
    linear = {
      # the two-parameter fit in closed form: no lm() object per segment
      m <- length(y)
      tt <- seq_len(m) - (m + 1) / 2
      y - mean(y) - tt * (sum(tt * y) / sum(tt * tt))
    },
    stop("Unknown detrend: ", detrend, call. = FALSE)
  )
}

#' One series to its ordinates. `segments` non-overlapping blocks are
#' averaged (Bartlett), which is what makes the result Gamma with the
#' shape at the block count rather than exponential.
#'
#' @noRd
pgram_one <- function(y, fs, taper, p, detrend, segments, label) {
  n <- length(y)
  if (!is.numeric(y) || anyNA(y)) {
    stop("frm_periodogram(): series ", label,
         " must be numeric with no NA", call. = FALSE)
  }
  m <- n %/% segments
  # floor((m - 1) / 2) ordinates survive: index 0 always goes, and index
  # m / 2 goes as well when m is even
  nf <- (m - 1L) %/% 2L
  if (nf < 1L) {
    stop("frm_periodogram(): series ", label, " has length ", n,
         " and segments = ", segments, ", which leaves ", nf,
         " usable ordinate(s). A segment needs at least 3 points, and ",
         "far more than that to say anything", call. = FALSE)
  }
  h <- taper_weights(m, taper, p)
  # sum(h^2) rather than m: the taper's own power is divided out, so
  # E I(w) stays the (smoothed) spectral density at any window
  scale <- segments * fs * sum(h * h)
  acc <- numeric(nf)
  keep <- seq.int(2L, nf + 1L)
  for (s in seq_len(segments)) {
    seg <- detrend_series(y[((s - 1L) * m + 1L):(s * m)], detrend)
    d <- stats::fft(h * seg)[keep]
    acc <- acc + Re(d)^2 + Im(d)^2
  }
  list(freq = seq_len(nf) * (fs / m), pgram = acc / scale)
}

#' Periodogram for a Whittle-likelihood fit
#'
#' Turns one or more time series into the data frame a spectral model is
#' fitted to: one row per retained Fourier frequency, with the ordinate
#' that [whittle()] treats as the response.
#'
#' The zero frequency and the Nyquist frequency are **dropped**. Their
#' ordinates are chi-square with one degree of freedom rather than two,
#' so they are not exponential and do not belong in the likelihood; a
#' fit that keeps them is wrong in a way nothing downstream reports. A
#' series of length `m` therefore returns `floor((m - 1) / 2)` rows.
#'
#' @section Scaling:
#' The ordinate is `|sum_t h_t (x_t - trend) exp(-2i pi f t / fs)|^2 /
#' (fs sum_t h_t^2)` and `freq` is in cycles per unit time (Hz when `fs`
#' is in Hz), which is the scaling and the frequency axis of
#' [stats::spec.pgram()]. So the ordinates of white noise average to
#' `var(x) / fs`, and a spectrum written for the model must be a density
#' per Hz. With `fs = 1` the frequency axis runs to just under `0.5`,
#' in cycles per sample.
#'
#' @section What a taper does to the likelihood:
#' A taper trades bias for dependence. Multiplying the series by a
#' window before the transform suppresses the leakage that makes a steep
#' spectrum's high frequencies too large, which is the main small-sample
#' bias of the Whittle likelihood; in exchange, neighboring ordinates
#' are no longer close to independent, so the likelihood is a working
#' one and standard errors from it are optimistic. The exponential
#' marginal survives: dividing by `sum(h^2)` keeps `E I(f)` at the
#' spectral density smoothed by the window, so the mean model is
#' unchanged and only the correlation between rows is not.
#'
#' @section Segment averaging:
#' `segments = k` splits the series into `k` non-overlapping blocks of
#' `floor(n / k)` points each (a remainder at the end is dropped),
#' periodograms each block, and averages them. The average of `k` independent
#' exponential ordinates is Gamma with shape `k` and the same mean, so
#' the matching family is `whittle(tapers = k)`: the two arguments are
#' one decision written twice, and this is the seam to get right. The
#' price is resolution, since the grid belongs to the block length
#' `floor(n / k)` rather than to `n`.
#'
#' @param x A numeric vector (one series), a `ts`, or a matrix or data
#'   frame whose **columns** are series of equal length.
#' @param fs Sampling rate in Hz. Defaults to 1 (cycles per sample), or
#'   to `stats::frequency(x)` when `x` is a `ts`.
#' @param taper `"none"`, `"hann"`, or `"split_cosine"` (the split
#'   cosine bell of [stats::spec.taper()], applied to a proportion `p`
#'   of each end).
#' @param p Proportion tapered at each end for `"split_cosine"`.
#' @param detrend `"mean"` (the default), `"none"`, or `"linear"`.
#'   `"mean"` and `"none"` return identical ordinates, because the
#'   retained frequencies are exactly those where a constant contributes
#'   nothing; `"linear"` does change them.
#' @param segments Number of non-overlapping blocks to average. See
#'   Segment averaging.
#' @param group Optional grouping vector as long as `x` (a vector only),
#'   splitting it into series that need not have equal lengths.
#' @return A data frame with columns `freq` (Hz), `pgram` (the
#'   ordinate), and, when the input names more than one series, `series`
#'   (a factor). Rows are ordered by series and then by frequency.
#' @seealso [whittle()] for the family, and `vignette("spectral")`.
#' @examples
#' set.seed(1)
#' y <- as.numeric(arima.sim(list(ar = 0.7), 512))
#' pg <- frm_periodogram(y, fs = 256)
#' str(pg)
#'
#' # the ordinates average to the variance, in per-Hz units
#' c(mean(pg$pgram), var(y) / 256)
#'
#' # one call for many series: columns are series
#' pg2 <- frm_periodogram(cbind(a = y, b = rev(y)), fs = 256)
#' table(pg2$series)
#' @export
frm_periodogram <- function(x, fs = NULL,
                            taper = c("none", "hann", "split_cosine"),
                            p = 0.1, detrend = c("mean", "none", "linear"),
                            segments = 1L, group = NULL) {
  taper <- match.arg(taper)
  detrend <- match.arg(detrend)
  if (is.null(fs)) fs <- if (stats::is.ts(x)) stats::frequency(x) else 1
  if (!is.numeric(fs) || length(fs) != 1L || !is.finite(fs) || fs <= 0) {
    stop("frm_periodogram(fs =) must be one positive number", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p <= 0 || p >= 0.5) {
    stop("frm_periodogram(p =) must be one number in (0, 0.5)", call. = FALSE)
  }
  segments <- as.integer(segments)
  if (length(segments) != 1L || is.na(segments) || segments < 1L) {
    stop("frm_periodogram(segments =) must be one positive integer",
         call. = FALSE)
  }

  # three input shapes, one internal shape: a named list of series
  multi <- TRUE
  if (is.matrix(x) || is.data.frame(x)) {
    if (!is.null(group)) {
      stop("frm_periodogram(): give `group` with a vector `x`, not with a ",
           "matrix, whose columns are already the series", call. = FALSE)
    }
    nms <- colnames(x) %||% paste0("V", seq_len(ncol(x)))
    series <- stats::setNames(
      lapply(seq_len(ncol(x)), function(j) {
        col <- if (is.data.frame(x)) x[[j]] else x[, j]
        # as.numeric() on a factor returns level CODES, which would
        # transform quietly and mean nothing
        if (!is.numeric(col)) {
          stop("frm_periodogram(): column ", nms[j], " is ", class(col)[1L],
               ", not numeric", call. = FALSE)
        }
        as.numeric(col)
      }), nms)
  } else if (!is.null(group)) {
    x <- as.numeric(x)
    if (length(group) != length(x)) {
      stop("frm_periodogram(): `group` has length ", length(group),
           " and `x` has length ", length(x), call. = FALSE)
    }
    if (anyNA(group)) {
      stop("frm_periodogram(): `group` has NA, and the rows it labels ",
           "would be dropped from every series without saying so",
           call. = FALSE)
    }
    g <- factor(group, levels = unique(as.character(group)))
    series <- split(x, g)
  } else {
    multi <- FALSE
    series <- list(as.numeric(x))
  }

  parts <- lapply(seq_along(series), function(i) {
    nm <- names(series)[i] %||% as.character(i)
    pgram_one(series[[i]], fs, taper, p, detrend, segments, nm)
  })
  out <- data.frame(
    freq = unlist(lapply(parts, `[[`, "freq"), use.names = FALSE),
    pgram = unlist(lapply(parts, `[[`, "pgram"), use.names = FALSE)
  )
  if (multi) {
    lens <- vapply(parts, function(z) length(z[["freq"]]), 0L)
    nms <- names(series)
    out[["series"]] <- factor(rep(nms, lens), levels = nms)
  }
  out
}

#' The one refusal the response alone can carry.
#'
#' Two errors reach here. A response that is not a periodogram at all,
#' which is almost always log power or power in dB and so goes negative.
#' And ordinates that are not independent exponentials, which happens
#' two ways: they were already averaged (Welch, Bartlett, multitaper)
#' and the family says they are raw, or they are LEAKAGE rather than
#' signal, which is what an untapered periodogram returns for a spectrum
#' falling faster than f^-2. Both leave the same fingerprint, so one
#' check catches both and the message names both.
#'
#' The second is detectable, and nearly one-sidedly. The log ratio of
#' two ordinates with the same mean has variance `2 trigamma(k)` at
#' shape k and nothing else enters it, so for neighbouring frequencies
#' `var(diff(log y))` is that variance plus the variance the SPECTRUM
#' contributes across one grid step. A spectrum that moves can only ADD
#' to it, and so can a random reordering of the rows (measured: a
#' shuffled raw periodogram gives 3.46 against a model 3.29, and is
#' accepted).
#'
#' The exception, and the reason this says "nearly": a reordering that
#' is MONOTONE in the response subtracts, and heavily. A raw
#' periodogram sorted by its own power gives 0.002 and is refused every
#' time. So the statistic assumes the rows are in frequency order,
#' which is what `frm_periodogram()` returns; sorting a periodogram by
#' power before fitting it breaks the check, and is one way to make
#' it refuse something legitimate.
#'
#' The threshold is a FRACTION of the model value, and the fraction
#' shrinks with the number of ordinates, because the sampling spread of
#' the statistic does the opposite. A flat half was the first rule and
#' it needed a hard floor at 200 ordinates to stay honest, which made
#' the check inert at exactly the length this feature is aimed at: a
#' one-second epoch at 256 Hz is 127 ordinates, and nothing fired
#' there. `whittle_smooth_frac()` replaces the floor with a curve
#' calibrated at the flat spectrum, the case that MINIMIZES the
#' statistic, so what is bought at small counts is bought by asking for
#' more evidence rather than by refusing to look.
#'
#' Measured in dev/freq-findings.md. At 127 ordinates the untapered
#' exponent-3 power law is now caught 92% of the time (it was 0%), a
#' k-averaged periodogram declared raw is caught 90 to 100% at k = 4
#' and 8 down to 31 ordinates, and no legitimate spectrum tested false-
#' fired: AR(1) at phi 0, 0.9 and 0.99, a 1/f background with a strong
#' alpha peak, and ten stacked 64-point series are all 0%. The two
#' known limits are k = 2 below about 200 ordinates, where Gamma(2) and
#' Gamma(1) are too close to separate (48%), and AR(1) at phi = 0.999,
#' which fires about 8% and is arguably a true positive because that
#' spectrum really is leakage-contaminated at that length.
#'
#' @noRd
whittle_smooth_frac <- function(nf) {
  # calibrated so that no more than about one flat-spectrum sample in
  # 10000 reaches it, with a 15% margin, monotone in nf and capped at
  # the half that the 200-ordinate rule used. Interpolated on log(nf)
  # and flat outside the knots.
  knots_nf <- c(24, 32, 48, 64, 96, 128, 192, 256, 384, 512)
  knots_fr <- c(0.123, 0.171, 0.236, 0.279, 0.359,
                0.400, 0.455, 0.494, 0.500, 0.500)
  stats::approx(log(knots_nf), knots_fr, xout = log(nf), rule = 2)$y
}

whittle_valid_y <- function(tapers) {
  force(tapers)
  function(y, aterms) {
    if (!is.numeric(y) || any(!is.finite(y))) {
      stop("whittle: the response must be finite periodogram ordinates",
           call. = FALSE)
    }
    if (any(y <= 0)) {
      stop("whittle: the response must be strictly positive. A ",
           "periodogram ordinate is a squared modulus; a response that ",
           "goes negative is log power or power in dB, which is not what ",
           "this likelihood is about. Model log power by fitting the ",
           "ordinates themselves - the log link already puts the linear ",
           "predictor on the log-spectrum scale", call. = FALSE)
    }
    # below this the null is too wide for any threshold to separate a
    # smoothed periodogram from a rough one, and a series that short
    # has no spectrum worth fitting either
    if (length(y) < 24L) return(invisible(NULL))
    v <- stats::var(diff(log(y)))
    thr <- whittle_smooth_frac(length(y)) * 2 * trigamma(tapers)
    if (v < thr) {
      stop("whittle(tapers = ", tapers, "): the response is far too ",
           "smooth to be that. var(diff(log(y))) is ",
           format(signif(v, 3)), ", the refusal triggers below ",
           format(signif(thr, 3)), ", and ordinates of shape ", tapers,
           " have an expected ", format(signif(2 * trigamma(tapers), 3)),
           " whatever the spectrum is, since only the spectrum's own ",
           "step-to-step variation adds to it. Two things look like ",
           "this. (1) The ordinates were already averaged - Welch, ",
           "Bartlett, multitaper - and this family says they are raw: ",
           "pass tapers = the number of periodograms that were ",
           "averaged. (2) They are spectral leakage rather than signal, ",
           "which is what an untapered periodogram returns for a ",
           "spectrum falling faster than f^-2, and the estimate from it ",
           "would be the leakage floor rather than the spectrum: pass ",
           "taper = \"hann\" to frm_periodogram(). If the ordinates ",
           "ALREADY carry a hann taper, that taper correlates ",
           "neighbouring ordinates and can trip this check on its ",
           "own: use taper = \"split_cosine\", which does not",
           call. = FALSE)
    }
    invisible(NULL)
  }
}

#' Whittle likelihood for a periodogram response
#'
#' Names the assumption behind a spectral fit. A periodogram ordinate of
#' a stationary Gaussian series is exponential about the spectral
#' density, so with a log link the linear predictor IS the log spectrum
#' and nothing else in the fit changes. `whittle()` is that exponential
#' family, with the periodogram-specific refusals attached.
#'
#' `tapers = k` for `k > 1` is the Gamma likelihood with the shape held
#' at `k`, which is the distribution of an average of `k` independent
#' ordinates: a Welch or Bartlett estimate from
#' `frm_periodogram(segments = k)`, or a `k`-taper multitaper estimate.
#' The shape is data preparation, not a parameter, so it is fixed rather
#' than estimated - through the same mechanism as `bf(y ~ x, shape = k)`,
#' and an explicit `shape` in [bf()] still wins.
#'
#' @section Post-processing speaks about the periodogram:
#' `fitted()` returns the fitted spectral density at each row's
#' frequency, `predict()` the same for new frequencies, and
#' `residuals()` compares an ordinate with that density. None of them
#' is about the series, and none of them can be: the phases are not in
#' a periodogram.
#'
#' `simulate()` is **refused by name** rather than documented, and so
#' are the four entry points built on it: `pp_check()`,
#' `frm_simulate()`, `frm_bootstrap()`, and
#' `conditional_effects(method = "predict")`.
#' Draws from this family are ordinates, a vector of them under the
#' response's name reads as a simulated series, and that is a mistake a
#' user makes once and never notices. Each refusal names the two things
#' to do instead: draw ordinates in the open, with
#' `rexp(nobs(fit), 1 / fitted(fit))` (or `rgamma(nobs(fit), shape = k,`
#' `scale = fitted(fit) / k)` at `tapers = k`), which is also the
#' one-line replacement for a DHARMa-style check; or draw a series
#' whose spectrum is the fitted one with [frm_series_draw()].
#'
#' @section What the family refuses, and what it cannot see:
#' Two things are checked, and both look only at the response.
#'
#' A response that is not positive is refused outright: a periodogram
#' ordinate is a squared modulus, so a negative one means log power or
#' power in dB, which this likelihood is not about.
#'
#' A response too SMOOTH to have the declared shape is also refused.
#' `var(diff(log(y)))` is `2 trigamma(tapers)` for ordinates of that
#' shape whatever the spectrum is, and only the spectrum's own
#' step-to-step variation adds to it, so a value far below that is
#' evidence the ordinates were averaged more than `tapers` says, or
#' that they are spectral leakage rather than signal. The threshold is
#' a calibrated fraction of the model value that shrinks as the number
#' of ordinates falls.
#'
#' **What it cannot see.** It needs the rows in frequency order, which
#' is what [frm_periodogram()] returns: a periodogram sorted by its own
#' power is refused every time, wrongly. It does not run at all below
#' 24 ordinates. It separates `tapers = 2` from a raw periodogram only
#' above about 200 ordinates, because Gamma(2) and Gamma(1) are close.
#' A Hann-tapered periodogram is legitimately raw, but the taper
#' correlates neighbouring ordinates: the lag-one correlation of
#' `log I` is about 0.3, which pulls the expected statistic from
#' 3.29 down to 2.28 against a trigger of 1.65, so such a response is
#' refused about one to two percent of the time. A split-cosine taper
#' does not do this (correlation 0.009). If a refusal names a response
#' you already tapered with Hann, that is this, and the remedy is not
#' another taper.
#' And it is a backstop, not a test to rely on: at 127 ordinates (a
#' one-second epoch at 256 Hz) an untapered exponent-3 power law is
#' caught about 92% of the time and an exponent-2 one much less often.
#' Choose the taper from the shape of the spectrum, not from whether an
#' error appeared. `vignette("spectral")` has the measured rates.
#'
#' @section What it is not:
#' The likelihood treats ordinates as independent, which is exact only
#' in the limit. It is biased for a short series or a spectrum near a
#' unit root; `vignette("spectral")` measures how much and says where it
#' is safe. Tapering reduces that bias and correlates the rows.
#'
#' @param tapers Number of independent periodograms averaged into each
#'   ordinate. `1` (the default) is the raw periodogram and gives the
#'   exponential likelihood; `k > 1` gives Gamma with the shape fixed at
#'   `k`.
#' @param link Link for the mean. The default `"log"` is what makes the
#'   linear predictor the log spectrum, and is very nearly always what
#'   is wanted.
#' @return A `frmtmb_family`.
#' @seealso [frm_periodogram()], [frm_series_draw()],
#'   `vignette("spectral")`.
#' @examples
#' set.seed(1)
#' y <- as.numeric(arima.sim(list(ar = 0.6), 512))
#' pg <- frm_periodogram(y)
#'
#' # a nonparametric log spectrum: the smoothing parameter is estimated
#' # by the same Laplace marginal likelihood as any other smooth
#' fit <- frm(bf(pgram ~ s(freq, k = 8)), family = whittle(), data = pg)
#' fixef(fit)
#'
#' # an averaged periodogram needs the shape it was averaged with
#' pg4 <- frm_periodogram(y, segments = 4)
#' frm(bf(pgram ~ log(freq)), family = whittle(tapers = 4), data = pg4)
#' @export
whittle <- function(tapers = 1, link = "log") {
  if (!is.numeric(tapers) || length(tapers) != 1L || is.na(tapers) ||
        tapers < 1 || tapers != round(tapers)) {
    stop("whittle(tapers =) must be one positive whole number: the count ",
         "of independent periodograms averaged into each ordinate",
         call. = FALSE)
  }
  tapers <- as.integer(tapers)
  fam <- if (tapers == 1L) fam_exponential(link) else fam_Gamma(link)
  fam[["family"]] <- "whittle"
  # cens() and trunc() have no meaning for an ordinate, and inheriting
  # them from exponential() would advertise a censored periodogram
  fam[["accepts_aterms"]] <- check_accepts_aterms("weights")
  fam[["valid_y"]] <- whittle_valid_y(tapers)
  fam[["tapers"]] <- tapers
  if (tapers > 1L) fam[["fixed_dpars"]] <- list(shape = tapers)
  # A draw here is a vector of ORDINATES. Handed back under the response
  # name it is read as a simulated time series, which it is not and
  # cannot be turned into: the phases are not in a periodogram. The
  # refusal is the point of the family (see the sim_refusal note), and
  # everything it costs is one line away, which the note supplies.
  fam[["sim"]] <- NULL
  fam[["sim_refusal"]] <- paste0(
    "This one is deliberate. A draw from a Whittle fit is a vector of ",
    "periodogram ORDINATES, and returned under the response's name it ",
    "is read as a simulated series, which it is not: the phases are not ",
    "in a periodogram. Draw ordinates in the open - rexp(nobs(fit), ",
    "1 / fitted(fit)), or rgamma(nobs(fit), shape = k, ",
    "scale = fitted(fit) / k) at tapers = k - and draw a SERIES with ",
    "the fitted spectrum using frm_series_draw(fit)"
  )
  fam
}

#' The frequency grid a fitted spectrum lives on, recovered from the
#' data the fit was given. A grid is `j * delta` for `j = 1..nf`; a
#' response with any other pattern is either several series stacked (the
#' common case, and the reason for the message) or not a periodogram.
#'
#' @noRd
pgram_grid <- function(f) {
  nf <- length(f)
  if (nf < 2L || is.unsorted(f) || any(!is.finite(f)) || f[1L] <= 0) {
    stop("frm_series_draw(): the frequencies must be an increasing ",
         "Fourier grid with no zero ordinate", call. = FALSE)
  }
  delta <- f[1L]
  if (max(abs(f - delta * seq_len(nf))) > 1e-8 * delta * nf) {
    stop("frm_series_draw(): the frequencies are not one evenly spaced ",
         "grid starting at the spacing. Several series stacked in one ",
         "data frame look like this; draw from one of them by passing ",
         "its rows as `newdata`", call. = FALSE)
  }
  delta
}

#' Draw a time series from a fitted spectrum
#'
#' The inverse of the transform a spectral fit is built on, and the verb
#' to reach for when `simulate()` gave ordinates where a series was
#' wanted. It draws a stationary Gaussian series whose spectral density
#' is the one the fit estimated, by giving each Fourier frequency an
#' independent complex Gaussian amplitude with the fitted variance and
#' transforming back.
#'
#' A simulated PERIODOGRAM cannot be inverted: the phases are not in it.
#' What is well defined is a draw from the fitted MODEL, which is what
#' this is, and it is a different object - two calls with the same
#' fitted spectrum give two unrelated series.
#'
#' The length is fixed by the grid at `2 * nf + 1`, an odd number, so
#' that the drawn series has no Nyquist ordinate to invent and exactly
#' zero mean. The sampling rate follows from the grid spacing, and the
#' result carries it as a `ts` frequency, so
#' `frm_periodogram(frm_series_draw(fit))` lands back on the same grid.
#'
#' @param object A `frmtmb_fit` with a [whittle()] family.
#' @param nsim Number of series to draw.
#' @param seed Optional integer seed.
#' @param newdata Rows to take the spectrum from, for a fit that pools
#'   several series. Defaults to the fitted rows.
#' @param freq The frequencies, either as a column name in the data or
#'   as the numbers themselves. A model written in terms of a
#'   transformed frequency keeps no `freq` column, so passing
#'   `pg$freq` is the general route.
#' @param ... Passed to [predict()], for example `re.form`.
#' @return A `ts` matrix with `2 * nf + 1` rows and `nsim` columns, at
#'   the sampling rate the grid implies.
#' @seealso [whittle()], [frm_periodogram()].
#' @examples
#' set.seed(1)
#' y <- as.numeric(arima.sim(list(ar = 0.6), 512))
#' pg <- frm_periodogram(y, fs = 128)
#' fit <- frm(bf(pgram ~ s(freq, k = 8)), family = whittle(), data = pg)
#'
#' # a new series with the same estimated spectrum, not the same series
#' z <- frm_series_draw(fit, nsim = 2, seed = 1)
#' dim(z)
#' stats::frequency(z)
#' @export
frm_series_draw <- function(object, nsim = 1, seed = NULL, newdata = NULL,
                            freq = "freq", ...) {
  if (!inherits(object, "frmtmb_fit")) {
    stop("frm_series_draw(object =) must be a frmtmb_fit", call. = FALSE)
  }
  rspec <- single_response(object, "frm_series_draw()")
  if (!identical(rspec$family[["family"]], "whittle")) {
    stop("frm_series_draw() needs a whittle() fit: it reads the fitted ",
         "values as a spectral density, and for family '",
         rspec$family[["family"]], "' they are not one", call. = FALSE)
  }
  check_count(nsim, "nsim", min = 1L)
  dat <- newdata %||% object$frame[["data_frame"]]
  if (is.numeric(freq) && length(freq) > 1L) {
    f <- as.numeric(freq)
  } else if (is.character(freq) && length(freq) == 1L) {
    if (!freq %in% names(dat)) {
      stop("frm_series_draw(): no column '", freq, "' in the data the ",
           "spectrum comes from. A model written in terms of a ",
           "transformed frequency (angular frequency, log frequency) ",
           "keeps no frequency column, so pass the frequencies ",
           "themselves: freq = pg$freq", call. = FALSE)
    }
    f <- as.numeric(dat[[freq]])
  } else {
    stop("frm_series_draw(freq =) must be one column name or the ",
         "frequencies themselves", call. = FALSE)
  }
  s <- as.numeric(stats::predict(object, newdata = newdata,
                                 type = "response", ...))
  if (length(s) != length(f)) {
    stop("frm_series_draw(): ", length(s), " fitted values for ",
         length(f), " frequencies", call. = FALSE)
  }
  ord <- order(f)
  f <- f[ord]
  s <- s[ord]
  delta <- pgram_grid(f)
  nf <- length(f)
  m <- 2L * nf + 1L
  fs <- delta * m
  if (!is.null(seed)) set.seed(seed)
  out <- matrix(0, m, nsim)
  # E|Z_j|^2 = m fs S_j is what makes |Z_j|^2 / (m fs) - the periodogram
  # of the result at this scaling - average to S_j
  amp <- sqrt(m * fs * s / 2)
  for (i in seq_len(nsim)) {
    z <- complex(real = stats::rnorm(nf, 0, amp),
                 imaginary = stats::rnorm(nf, 0, amp))
    full <- c(0 + 0i, z, Conj(rev(z)))
    out[, i] <- Re(stats::fft(full, inverse = TRUE)) / m
  }
  stats::ts(out, frequency = fs)
}
