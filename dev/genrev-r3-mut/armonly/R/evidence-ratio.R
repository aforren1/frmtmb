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
# WHERE THE TWO DENSITIES COME FROM, and why they are not symmetric.
# The numerator is a kernel density of the drawn quantity, because that
# is the only thing the draws can give. The denominator is EVALUATED
# from the prior specification the model was sampled under, because
# that specification is known exactly and estimating it from prior
# draws would put Monte Carlo error in a place where none is needed.
# brms estimates both by kernel density, from `sample_prior = "yes"`
# draws; the two agree to the accuracy of brms's denominator, which is
# the noisier half (a weakly-informative prior spreads its draws over a
# range where a fixed-bandwidth estimate at one point is thin).
#
# WHAT REFUSES, AND WHY BY NAME. A denominator needs a proper prior on
# the quantity being tested. frm_sample() leaves class "b" FLAT by
# default, exactly as brms does, so the commonest point hypothesis of
# all has no denominator unless the user wrote a prior. That is a fact
# about the model, not a failure of the method, and saying which
# hypothesis lacks which prior is the only useful thing to report; brms
# reaches the same place from the other side, reporting NA when no
# prior draws were saved.

#' The density of one resolved prior entry at one value, or `NULL` for
#' a `dist` whose density is not one number per parameter.
#'
#' The kinds are core's `prior_base_logdens()` list; this reads the
#' resolved specification rather than reaching into core's tape-side
#' evaluator, which is written for `RTMB` values and carries the change
#' of variables an outer parameter needs and this quantity does not.
#'
#' @noRd
er_prior_density <- function(dist, v) {
  switch(dist$kind,
    normal = stats::dnorm(v, dist$location, dist$scale),
    t = stats::dt((v - dist$location) / dist$scale, dist$df) / dist$scale,
    exponential = stats::dexp(v, dist$rate),
    logistic = stats::dlogis(v, dist$location, dist$scale),
    gamma = stats::dgamma(v, shape = dist$shape, rate = dist$rate),
    # written out because R parameterizes the inverse gamma nowhere;
    # Stan's second argument (and so brms's) is the SCALE
    inv_gamma = if (v <= 0) 0 else {
      exp(dist$shape * log(dist$scale) - lgamma(dist$shape) -
            (dist$shape + 1) * log(v) - dist$scale / v)
    },
    beta = stats::dbeta(v, dist$shape1, dist$shape2),
    NULL)
}

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

#' The resolved prior entries the model was sampled under, or the
#' reason there are none to read.
#'
#' @noRd
er_entries <- function(fit) {
  pl <- fit[["prior"]]
  if (is.null(pl) || !length(unclass(pl))) {
    return(list(entries = NULL,
                why = paste("the model carries no proper prior (it was",
                            "sampled flat), so there is no density to",
                            "divide by")))
  }
  if (length(attr(pl, "overrides") %||% list())) {
    return(list(entries = NULL,
                why = paste("the retired prior = list(name = ...)",
                            "spelling took over parameters this",
                            "specification no longer describes, so the",
                            "density that was sampled cannot be read",
                            "back from it")))
  }
  e <- tryCatch(resolve_prior_input(fit, pl)$entries,
                error = function(err) NULL)
  if (!length(e)) {
    return(list(entries = NULL,
                why = paste("its prior specification does not resolve",
                            "against this model")))
  }
  list(entries = e, why = NULL)
}

#' Position of each name in the `(beta, estimated betad)` coefficient
#' vector `hyp_vals_only()` lays out first, `NA` for a name that is not
#' a coefficient.
#'
#' A random-effect standard deviation, a correlation and a
#' natural-scale summary are all excluded on purpose. Their priors sit
#' on the INTERNAL parameter (a log standard deviation, a transformed
#' correlation) with a change of variables in between, so the prior of
#' the reported quantity is not the entry's density; and a point null
#' at zero on a variance component is a boundary hypothesis that a
#' density ratio does not answer anyway.
#'
#' @noRd
er_coef_pos <- function(fit, nms) {
  raw <- estimated_coef_names(fit)
  bare <- gsub("[()]", "", raw)
  pos <- match(nms, bare)
  ifelse(is.na(pos), match(nms, raw), pos)
}

#' The template slot one coefficient position occupies.
#'
#' @noRd
er_coef_slot <- function(fit, i) {
  tpl <- fit$frame[["par_template"]]
  nb <- length(tpl[["beta"]])
  if (i <= nb) return(list(comp = "beta", idx = i))
  fx <- fit$frame[["betad_fixed_idx"]]
  keep <- setdiff(seq_along(tpl[["betad"]]), fx)
  list(comp = "betad", idx = keep[i - nb])
}

#' The entry covering one slot, or `NULL`.
#'
#' @noRd
er_entry_for <- function(entries, slot) {
  for (e in entries) {
    if (identical(e$comp, slot$comp) && slot$idx %in% e$idx) return(e)
  }
  NULL
}

#' The hypothesis quantity as `a + sum(b * coefficients)`, or `NULL`
#' when it is not affine in them.
#'
#' MEASURED, not parsed. The expression engine already evaluates, so
#' `m + 1` evaluations give the intercept and the slopes and two more
#' check the claim at points none of them used. Parsing the expression
#' instead would have to know every function a hypothesis may contain.
#'
#' BOTH SIGNS are probed, and the negative one is the load-bearing
#' half. The fitting points are `0` and the unit vectors, all
#' non-negative, so `abs()` is indistinguishable from the identity on
#' every one of them and on any non-negative check point. It slipped
#' through, and the prior of `|X|` for `X ~ N(0, 1)` is a half-normal
#' whose density at zero is twice the normal's - so `abs(x) = 0`
#' reported a Bayes factor exactly two times too large, silently, on a
#' hypothesis a brms user would plausibly write. (It is also a boundary
#' null, which is a second reason to refuse it.)
#'
#' @noRd
er_affine <- function(fit, ex, vo, pos) {
  at <- function(v) {
    vv <- vo$vals
    vv[pos] <- v
    hyp_eval(fit, ex, vv, vo$comp)
  }
  m <- length(pos)
  a <- at(rep(0, m))
  b <- numeric(m)
  for (j in seq_len(m)) {
    ej <- numeric(m)
    ej[j] <- 1
    b[j] <- at(ej) - a
  }
  chk <- seq_len(m) * 0.7 - 0.3
  for (probe in list(chk, -chk)) {
    if (!isTRUE(all.equal(at(probe), a + sum(b * probe),
                          tolerance = 1e-8))) {
      return(NULL)
    }
  }
  list(a = a, b = b)
}

#' The prior density of one hypothesis quantity at zero, or the reason
#' there is none.
#'
#' @noRd
er_prior_at_zero <- function(fit, ex, vo, entries) {
  nms <- all.vars(ex)
  if (!length(nms)) {
    return(list(why = paste("it names no parameter, so there is no",
                            "prior density to read")))
  }
  pos <- er_coef_pos(fit, nms)
  if (anyNA(pos)) {
    # every name that reaches here resolved in the hypothesis
    # environment, so a name that is not a coefficient is one of the
    # NATURAL-SCALE summaries that environment also carries. Saying
    # "not a population-level coefficient" about `sigma` was true of
    # the name and false of the model - sigma IS a coefficient here,
    # reported on a different scale - and a user who reads that stops
    # believing the package.
    return(list(why = paste0("it names ",
                             paste(nms[is.na(pos)], collapse = ", "),
                             ", a natural-scale summary rather than a ",
                             "population-level coefficient. A standard ",
                             "deviation, a correlation and a dispersion ",
                             "are each reported on a scale their prior ",
                             "is not written on: the prior sits on the ",
                             "internal parameter with a change of ",
                             "variables in between, so its density is ",
                             "not the density of the tested quantity")))
  }
  slots <- lapply(pos, function(i) er_coef_slot(fit, i))
  ent <- lapply(slots, function(s) er_entry_for(entries, s))
  flat <- vapply(ent, is.null, TRUE)
  if (any(flat)) {
    return(list(why = paste0(paste(nms[flat], collapse = ", "),
                             " has no proper prior (frm_sample() leaves ",
                             "class \"b\" flat, as brms does); write one ",
                             "with set_prior(class = \"b\") and resample ",
                             "to get a Bayes factor for it")))
  }
  bad <- vapply(ent, function(e) {
    !identical(e$scale, "internal") || !is.null(e$offset)
  }, TRUE)
  if (any(bad)) {
    return(list(why = paste0("the prior on ",
                             paste(nms[bad], collapse = ", "),
                             " is not written about that coefficient ",
                             "itself (an Intercept prior is about the ",
                             "intercept at the predictor means, and a ",
                             "class \"sd\" prior is about a standard ",
                             "deviation), so its density is not the ",
                             "density of the tested quantity")))
  }
  aff <- er_affine(fit, ex, vo, pos)
  if (is.null(aff)) {
    return(list(why = paste("it is not an affine function of the",
                            "coefficients it names, and the prior of a",
                            "nonlinear function of them has no closed",
                            "form here")))
  }
  dists <- lapply(ent, function(e) e$dist)
  if (length(pos) == 1L) {
    # a + b v = 0 at v = -a / b, and a density transforms by |b|
    if (abs(aff$b[1L]) < .Machine$double.eps) {
      return(list(why = "the tested quantity does not depend on it"))
    }
    d <- er_prior_density(dists[[1L]], -aff$a / aff$b[1L])
    if (is.null(d)) {
      return(list(why = paste0("its prior is a ", dists[[1L]]$kind,
                               " density, which is not one number per ",
                               "coefficient")))
    }
    return(list(density = d / abs(aff$b[1L])))
  }
  # several coefficients: only the normal convolves in closed form, and
  # only because the entries are independent by construction
  if (!all(vapply(dists, function(d) identical(d$kind, "normal"), TRUE))) {
    return(list(why = paste("it combines several coefficients whose",
                            "priors are not all normal, and the prior",
                            "of their combination has no closed form",
                            "here; test them one at a time, or give",
                            "them normal priors")))
  }
  mu <- aff$a + sum(aff$b * vapply(dists, function(d) d$location, 1))
  sdq <- sqrt(sum((aff$b * vapply(dists, function(d) d$scale, 1))^2))
  list(density = stats::dnorm(0, mu, sdq))
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
#' hypotheses, with the between-block Monte Carlo error of each ratio.
#'
#' `point` says which rows were written as a POINT null. A hypothesis
#' written with no comparison at all is this package's own spelling for
#' "summarize this quantity", not a test; it parses `two.sided` and
#' gets no ratio and no complaint, because there is nothing there to
#' weigh evidence about.
#'
#' @noRd
er_evidence <- function(x, exs, dirs, draws, labels, point) {
  fit <- x$fit
  n <- nrow(draws)
  k_n <- length(exs)
  er <- rep(NA_real_, k_n)
  mcse <- rep(NA_real_, k_n)
  # resolved only when some row actually needs a denominator: a call
  # with only directional rows must not pay for, or fail on, a prior
  # specification it never reads
  ent <- NULL
  vo <- hyp_vals_only(fit)
  blocks <- er_blocks(x, n)
  reasons <- character(0)
  for (k in seq_len(k_n)) {
    v <- draws[, k]
    if (identical(dirs[k], "less") || identical(dirs[k], "greater")) {
      # brms's posterior odds of the claim; the prior never enters a
      # directional ratio, so nothing here can refuse
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
    if (!isTRUE(point[k])) next
    if (is.null(ent)) ent <- er_entries(fit)
    pd <- if (is.null(ent$entries)) list(why = ent$why) else {
      er_prior_at_zero(fit, exs[[k]], vo, ent$entries)
    }
    if (!is.null(pd$why)) {
      reasons <- c(reasons, paste0(labels[k], ": ", pd$why))
      next
    }
    if (!is.finite(pd$density) || pd$density <= 0) {
      reasons <- c(reasons,
                   paste0(labels[k], ": its prior density at the tested ",
                          "point is zero, so the ratio is not defined"))
      next
    }
    er[k] <- er_kde_at(v) / pd$density
    if (!is.null(blocks)) {
      # the denominator is exact, so the spread of the ratio is the
      # spread of its numerator
      per <- vapply(blocks, function(ix) er_kde_at(v[ix]), 1) / pd$density
      per <- per[is.finite(per)]
      if (length(per) > 1L) mcse[k] <- stats::sd(per) / sqrt(length(per))
    }
  }
  if (length(reasons)) {
    np <- sum(point, na.rm = TRUE)
    warning("hypothesis(): no evidence ratio for ", length(reasons),
            " of ", np, " point hypothes", if (np == 1L) "is" else "es",
            ". A Savage-Dickey ratio divides the posterior density at ",
            "the tested point by the PRIOR density there, and\n  ",
            paste(reasons, collapse = "\n  "), call. = FALSE)
  }
  pp <- ifelse(is.infinite(er), 1, er / (1 + er))
  list(evid_ratio = er, post_prob = pp, mcse = mcse)
}
