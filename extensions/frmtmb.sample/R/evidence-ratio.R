# The Savage-Dickey evidence ratio brms reports as `Evid.Ratio`, and
# the posterior probability beside it.
#
# WHAT THE NUMBER IS. For a POINT hypothesis the ratio is the Bayes
# factor of the point null against the alternative, which for a nested
# point null is the posterior density of the hypothesis quantity at the
# point divided by its PRIOR density there (Dickey and Lientz 1970).
# For a DIRECTIONAL hypothesis there is no density involved at all: it
# is the posterior odds of the claim, and the prior never enters.
#
# WHERE THE TWO DENSITIES COME FROM. brms's `evidence_ratio()`: both
# are kernel density estimates at the point, the numerator from the
# posterior draws of the hypothesis quantity and the denominator from
# the same quantity evaluated on the prior draws `frm_sample(sample_prior
# = "yes")` stored (`prior_draws()`). Without prior draws for every
# parameter the hypothesis names, the ratio is NA, with no message, as
# in brms. Through frmtmb.sample 0.23.0 the denominator was evaluated
# from the prior specification instead, and only for affine hypotheses
# of coefficients with a proper prior; that is gone, because it gave a
# number where brms gives NA and refused where brms gives one (a
# standard deviation, sigma).

#' The posterior density of a drawn quantity at one point.
#'
#' brms's construction, rebuilt from what `brms:::density_ratio()`
#' does: widen the evaluation range so it contains the point, take
#' `stats::density()` on 4096 grid points, read the grid with a spline
#' at the point, and floor the result at zero. Matching it is the
#' point - two Savage-Dickey ratios computed with different smoothers
#' are not comparable numbers, so a validation against brms has to use
#' brms's smoother.
#'
#' THE BANDWIDTH is `stats::density()`'s default, `bw = "nrd0"`:
#' Silverman's rule of thumb, `0.9 min(sd, IQR / 1.34) n^(-1/5)`. It is
#' named here rather than left implicit because the ratio is only as
#' reproducible as its bandwidth rule, and a reader comparing this
#' number to one from another package needs to know which rule made it.
#'
#' @noRd
er_kde_at <- function(v, point = 0, n = 4096L) {
  v <- v[is.finite(v)]
  if (length(v) < 2L || stats::sd(v) <= 0) return(NA_real_)
  from <- min(v)
  to <- max(v)
  # density() evaluates on [from, to] only, so a point outside the
  # drawn range would be read off the end of the grid instead of
  # estimated there
  if (from > point) {
    from <- point - stats::sd(v) / 4
  } else if (to < point) {
    to <- point + stats::sd(v) / 4
  }
  d <- stats::density(v, from = from, to = to, n = n)
  max(stats::spline(d$x, d$y, xout = point)$y, 0)
}

#' The hypothesis quantity evaluated on the prior draws of the
#' parameters it names, or `NULL` where any of them has none: brms's
#' `eval_hypothesis()` reads `prior_draws(x, variable = parsH)` and uses
#' it only when it has a column for every name.
#'
#' @noRd
er_prior_values <- function(x, ex) {
  vars <- hyp_expr_vars(ex)
  if (!length(vars)) return(NULL)
  pd <- prior_draws.frmtmb_draws(x, variable = vars)
  if (is.null(pd) || ncol(pd) != length(vars)) return(NULL)
  as.numeric(hyp_eval_in(ex, as.list(pd)))
}

#' The consecutive blocks of draws the Monte Carlo error is measured
#' across: the sampler's own chains, but never fewer than four.
#'
#' Between-BLOCK spread rather than a bootstrap over rows: the draws
#' are autocorrelated, and resampling rows independently would report a
#' Monte Carlo error smaller than the one the chains actually have.
#'
#' FOUR IS A FLOOR because two chains give `sd()` of two numbers over
#' `sqrt(2)`, which is one degree of freedom and not an error bar at
#' all: on a measured pair of 2-chain runs of one model it understated
#' a real 0.46 gap as 0.18. Splitting each chain in half keeps the
#' blocks consecutive, so within-block autocorrelation is still
#' respected - it is the split-chain construction the sampler
#' diagnostics already use. More than four chains keep their own
#' boundaries.
#'
#' @noRd
er_blocks <- function(x, n) {
  nc <- tryCatch(nchains(x), error = function(e) 1L)
  if (is.na(nc)) nc <- 1L
  nc <- max(nc, 4L)
  if (n < 8L * nc) return(NULL)
  split(seq_len(n), rep(seq_len(nc), each = n %/% nc,
                        length.out = n))
}

#' Evidence ratios and posterior probabilities for one set of parsed
#' hypotheses, with the between-block Monte Carlo error of each ratio
#' and the prior draws of each hypothesis quantity (`NA` where there are
#' none), which brms returns as `prior_samples`.
#'
#' `point` says which rows were written as a POINT null. The Monte
#' Carlo error of a point ratio is the spread of its numerator between
#' blocks of the chains, over the one denominator of all the prior
#' draws, which are independent.
#'
#' @noRd
er_evidence <- function(x, exs, dirs, draws, labels, point) {
  n <- nrow(draws)
  k_n <- length(exs)
  er <- rep(NA_real_, k_n)
  mcse <- rep(NA_real_, k_n)
  # logical NA where no hypothesis has prior draws, as brms's frame is
  prior <- matrix(NA, n, k_n)
  blocks <- er_blocks(x, n)
  for (k in seq_len(k_n)) {
    v <- draws[, k]
    # brms reads the prior draws of every hypothesis, directional or not,
    # and returns them with the result
    pv <- er_prior_values(x, exs[[k]])
    if (!is.null(pv) && length(pv) == n) prior[, k] <- pv
    if (identical(dirs[k], "less") || identical(dirs[k], "greater")) {
      # brms's posterior odds of the claim; the prior never enters a
      # directional ratio
      hit <- if (identical(dirs[k], "less")) sum(v < 0) else sum(v > 0)
      er[k] <- hit / (n - hit)
      if (!is.null(blocks)) {
        per <- vapply(blocks, function(ix) {
          vb <- v[ix]
          h <- if (identical(dirs[k], "less")) sum(vb < 0) else sum(vb > 0)
          h / (length(vb) - h)
        }, 1)
        per <- per[is.finite(per)]
        if (length(per) > 1L) mcse[k] <- stats::sd(per) / sqrt(length(per))
      }
      next
    }
    if (!isTRUE(point[k]) || is.null(pv)) next
    den <- er_kde_at(pv)
    er[k] <- er_kde_at(v) / den
    if (!is.null(blocks)) {
      per <- vapply(blocks, function(ix) er_kde_at(v[ix]), 1) / den
      per <- per[is.finite(per)]
      if (length(per) > 1L) mcse[k] <- stats::sd(per) / sqrt(length(per))
    }
  }
  pp <- ifelse(is.infinite(er), 1, er / (1 + er))
  list(evid_ratio = er, post_prob = pp, mcse = mcse, prior = prior)
}
