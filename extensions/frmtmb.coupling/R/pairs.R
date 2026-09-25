#' The cross-spectra of several channel pairs, stacked with a `pair` factor
#'
#' [cross_wishart()] models ONE channel pair. A recording with more
#' channels has a pair for every two of them, and the questions asked of
#' it are about all of those pairs at once. `frm_cross_pairs()` runs
#' [frm_cross_spectrum()] on each pair and stacks the results, one block
#' per pair, with a `pair` factor that a model can use.
#'
#' @section Many pairs in one model:
#' Fitting each pair on its own and correcting the p-values afterwards
#' treats the pairs as unrelated. One model over the stacked frame lets
#' them share information instead:
#'
#' * `coh ~ 1 + (1 | pair)` gives each pair its own coherence, shrunk
#'   toward the pairs' common level by as much as the data say the pairs
#'   differ. A pair that stands out after shrinkage stands out against
#'   the others, which is the comparison a multiple-comparison correction
#'   tries to make after the fact.
#' * `coh ~ s(freq, by = pair)` with `pair` also a fixed effect gives
#'   each pair its own coherence spectrum, and [frm_coherence()] reads
#'   each one off with a band.
#' * `coh ~ 0 + pair`, with `mu`, `pow2` and `phase` also by `pair`, is
#'   the separate fits in one call: the likelihood is then a sum of one
#'   term per pair, so each pair's estimate is the one its own fit gives.
#'
#' The power of a channel appears in every pair that includes it, so the
#' blocks share channels. Each block is still a correct likelihood for
#' its own pair; what the stacked likelihood does not model is the
#' dependence between two blocks that share a channel. Read the
#' standard errors of a stacked fit with that in mind.
#'
#' @param X The recording: a numeric matrix or a data frame of numeric
#'   columns whose COLUMNS ARE CHANNELS and whose rows are samples, or a
#'   list of such matrices, one per epoch, each with the same channels.
#'   Channel names come from the column names, or are `ch1`, `ch2`, and
#'   so on when there are none.
#' @param pairs Which pairs to form. `NULL`, the default, forms every
#'   pair of two different channels, in column order. Otherwise a
#'   two-column matrix or data frame, one row per pair, naming each
#'   channel by name or by column number.
#' @param ... Passed to [frm_cross_spectrum()]: `sfreq`, `segments`,
#'   `tapers`, `smooth`, `window`, `frange` and, for a list of epochs,
#'   `group`.
#'
#' @return The rows of [frm_cross_spectrum()] for every pair, with the
#'   factors `pair` (`"a-b"` for channels `a` and `b`), `ch1` and `ch2`
#'   in front. `ch1` is the signal `frm_cross_spectrum()` takes as `x`,
#'   so the fitted phase is the lead of `ch1` over `ch2`.
#'
#' @seealso [frm_cross_spectrum()], [cross_wishart()].
#'
#' @examples
#' set.seed(2)
#' src <- rnorm(2048)
#' X <- cbind(Fz = src + rnorm(2048), Cz = 0.8 * src + rnorm(2048),
#'            Pz = rnorm(2048))
#' xp <- frm_cross_pairs(X, sfreq = 256, segments = 8)
#' table(xp$pair)
#' @export
frm_cross_pairs <- function(X, pairs = NULL, ...) {
  dots <- list(...)
  allowed <- c("sfreq", "segments", "tapers", "smooth", "window", "frange",
               "group")
  dn <- names(dots)
  if (length(dots) && (is.null(dn) || any(!nzchar(dn)))) {
    frm_stop("frm_cross_pairs(): every argument after `pairs` must be ",
             "named. The ones it passes on are: ",
             paste(allowed, collapse = ", "), ".", call. = FALSE)
  }
  bad <- setdiff(dn, allowed)
  if (length(bad)) {
    frm_stop("frm_cross_pairs(): unknown argument",
             if (length(bad) > 1L) "s " else " ",
             paste0("`", bad, "`", collapse = ", "),
             ". It passes ", paste(allowed, collapse = ", "),
             " to frm_cross_spectrum(); the two signals of each pair come ",
             "from `X` and `pairs`.", call. = FALSE)
  }
  epochs <- if (is.list(X) && !is.data.frame(X)) X else list(X)
  if (!length(epochs)) {
    frm_stop("frm_cross_pairs(): `X` is an empty list.", call. = FALSE)
  }
  epochs <- lapply(seq_along(epochs), function(i) {
    cp_channels(epochs[[i]], if (is.list(X) && !is.data.frame(X)) i)
  })
  chans <- colnames(epochs[[1L]])
  for (i in seq_along(epochs)[-1L]) {
    if (!identical(colnames(epochs[[i]]), chans)) {
      frm_stop("frm_cross_pairs(): epoch ", i, " has channels ",
               paste(colnames(epochs[[i]]), collapse = ", "),
               " and epoch 1 has ", paste(chans, collapse = ", "),
               ". Every epoch needs the same channels in the same order.",
               call. = FALSE)
    }
  }
  if (length(chans) < 2L) {
    frm_stop("frm_cross_pairs(): `X` has ", length(chans), " channel. A ",
             "pair needs two.", call. = FALSE)
  }
  pr <- cp_pairs(pairs, chans)
  is_list <- is.list(X) && !is.data.frame(X)
  labels <- paste(chans[pr[, 1L]], chans[pr[, 2L]], sep = "-")
  out <- lapply(seq_len(nrow(pr)), function(k) {
    a <- pr[k, 1L]
    b <- pr[k, 2L]
    xa <- lapply(epochs, function(m) m[, a])
    xb <- lapply(epochs, function(m) m[, b])
    if (!is_list) {
      xa <- xa[[1L]]
      xb <- xb[[1L]]
    }
    d <- do.call(frm_cross_spectrum, c(list(xa, xb), dots))
    cbind(pair = factor(labels[k], levels = labels),
          ch1 = factor(chans[a], levels = chans),
          ch2 = factor(chans[b], levels = chans), d)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' One epoch of a multichannel recording, as a numeric matrix with
#' channel names.
#'
#' @noRd
cp_channels <- function(m, i) {
  what <- if (is.null(i)) "`X`" else paste0("epoch ", i, " of `X`")
  if (is.data.frame(m)) {
    num <- vapply(m, is.numeric, TRUE)
    if (!all(num)) {
      frm_stop("frm_cross_pairs(): ", what, " has the non-numeric column",
               if (sum(!num) > 1L) "s " else " ",
               paste(names(m)[!num], collapse = ", "),
               ". Every column is a channel.", call. = FALSE)
    }
    m <- as.matrix(m)
  }
  if (!is.matrix(m) || !is.numeric(m)) {
    frm_stop("frm_cross_pairs(): ", what, " must be a numeric matrix or ",
             "data frame whose columns are channels.", call. = FALSE)
  }
  if (is.null(colnames(m))) colnames(m) <- paste0("ch", seq_len(ncol(m)))
  cn <- colnames(m)
  if (anyNA(cn) || any(!nzchar(cn)) || anyDuplicated(cn)) {
    frm_stop("frm_cross_pairs(): the channel names of ", what, " must be ",
             "unique and non-empty.", call. = FALSE)
  }
  m
}

#' The pairs as a two-column matrix of channel positions.
#'
#' @noRd
cp_pairs <- function(pairs, chans) {
  if (is.null(pairs)) {
    n <- length(chans)
    return(do.call(rbind, lapply(seq_len(n - 1L), function(i) {
      cbind(i, seq.int(i + 1L, n))
    })))
  }
  if (is.data.frame(pairs)) pairs <- as.matrix(pairs)
  if (!is.matrix(pairs) || ncol(pairs) != 2L || !nrow(pairs)) {
    frm_stop("frm_cross_pairs(): `pairs` must be a two-column matrix or ",
             "data frame with one row per pair.", call. = FALSE)
  }
  if (is.numeric(pairs)) {
    if (anyNA(pairs) || any(pairs != round(pairs)) || any(pairs < 1) ||
          any(pairs > length(chans))) {
      frm_stop("frm_cross_pairs(): a channel number in `pairs` is not a ",
               "column of `X`, which has ", length(chans), ".",
               call. = FALSE)
    }
    pos <- matrix(as.integer(pairs), ncol = 2L)
  } else {
    pos <- matrix(match(as.character(pairs), chans), ncol = 2L)
    if (anyNA(pos)) {
      frm_stop("frm_cross_pairs(): `pairs` names ",
               paste(unique(as.character(pairs)[is.na(pos)]),
                     collapse = ", "),
               ", which ", if (sum(is.na(pos)) > 1L) "are" else "is",
               " not a channel of `X`. The channels are ",
               paste(chans, collapse = ", "), ".", call. = FALSE)
    }
  }
  if (any(pos[, 1L] == pos[, 2L])) {
    i <- which(pos[, 1L] == pos[, 2L])[1L]
    frm_stop("frm_cross_pairs(): row ", i, " of `pairs` pairs channel ",
             chans[pos[i, 1L]], " with itself. Its coherence is 1 by ",
             "arithmetic.", call. = FALSE)
  }
  key <- paste(pmin(pos[, 1L], pos[, 2L]), pmax(pos[, 1L], pos[, 2L]))
  if (anyDuplicated(key)) {
    i <- anyDuplicated(key)
    frm_stop("frm_cross_pairs(): row ", i, " of `pairs` repeats the pair ",
             chans[pos[i, 1L]], " and ", chans[pos[i, 2L]], ", in one order ",
             "or the other. A pair counted twice counts its data twice.",
             call. = FALSE)
  }
  pos
}
