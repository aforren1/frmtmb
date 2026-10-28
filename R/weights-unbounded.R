# Refusal of a weighted likelihood (or posterior) that has no maximum
# (or does not integrate) because every row owns a random effect.

# The dpars that concentrate each family's density on one point when
# they go to their limit together. Each family was measured to run away
# with an observation-level intercept and weights of 1.2 and to fit with
# weights of 0.9 (dev/fuzzfix-spike-families.R, dev/fuzzfix-fr1-families.R).
# exgaussian needs both its scales: with beta held, its density stays
# below 1 / beta. zero_inflated_asym_laplace is asym_laplace on its
# nonzero rows (dev/rel073-zial-probe.R: Laplace logLik 278.05 and
# 201.90 unweighted, quadrature bounded at unit weights and unbounded
# at 1.2, a group of 4 rows bounded).
unbounded_spike_dpar <- list(
  gaussian = "sigma", student = "sigma", lognormal = "sigma",
  Gamma = "shape", weibull = "shape", beta = "phi",
  inverse.gaussian = "shape", skew_normal = "sigma",
  shifted_lognormal = "sigma", huber = "sigma", von_mises = "kappa",
  frechet = "nu", gen_extreme_value = "sigma",
  exgaussian = c("sigma", "beta"), asym_laplace = "sigma",
  zero_inflated_asym_laplace = "sigma")

# The limit each spike dpar runs to: a scale shrinks to 0, a precision
# or shape grows without bound
unbounded_spike_toward_zero <- c("sigma", "beta")

# Which covariance structures can give a row its own effect at any value
# with a finite density is a property of the structure, read through
# cs_full_rank() (R/covstruct.R): a built-in entry states it, and a
# registered structure states it or has it read off its covariance.

#' Refuse a weighted likelihood that has no maximum, or a posterior that
#' does not integrate.
#'
#' When every row with a positive weight has a random effect of its own
#' in `mu`, each row's density can be centered on its own observation.
#' As the spike dpars go to their limit together (`sigma` to 0, `shape`
#' to infinity) each density at its center grows like 1 / width, and the
#' marginal likelihood, exact or Laplace, grows like width^(n - sum(w)):
#' a row contributes w log(1 / width) from its density and
#' -log(1 / width) from the curvature of its own effect. With
#' `sum(w) > n` it has no maximum, under ML, REML and `profile = TRUE`
#' alike. asym_laplace is worse under the Laplace approximation: its log
#' density is linear on each side of the mode, the data give a row's own
#' effect no curvature, and the approximation grows like
#' width^(-sum(w)) for any weights.
#'
#' The proof needs three things, and the check fails open without any
#' of them: (i) every spike dpar reaches its limit at every row at once,
#' which holds when its design has a column and spans the constant (an
#' intercept); (ii) each row's mean reaches its observation, which an
#' owned effect in a full-rank block gives; (iii) the effects' prior
#' density is positive and finite there. It also fails open for several
#' responses, `se()`, residual correlation and `mi()`.
#'
#' Under `cens()` the argument runs over the uncensored rows. A censored
#' row's probability is at most 1, and its own effect can bring it to 1,
#' so the exact marginal of such a row stays between a positive constant
#' and 1 while the uncensored rows run away as before: the excess is
#' theirs alone. Under the Laplace approximation a censored row's own
#' effect sits where its probability turns from 0 to 1, where the
#' curvature grows like 1 / width, so each such row takes back up to
#' half a power of the growth; the check subtracts 1/2 per censored row,
#' which fails open where that is too much. On 80 gaussian rows with 16
#' censored (dev/rel073-wcens-probe.R) an uncensored excess of 2 and 4
#' stopped at an interior estimate and 16 and 32 ran away.
#'
#' A prior on a spike dpar exempts the model only when it provably
#' bounds the growth (`prior_spike_verdict()`); a prior on anything else
#' never does. `route = "sample"` asks the posterior question instead:
#' the density must integrate, which a prior finite at the limit allows
#' only while sum(w) - n stays below 1.
#'
#' @param exact `TRUE` when the integral over the effects is exact
#'   (`quadrature = TRUE`, or full sampling), `FALSE` under Laplace.
#' @noRd
check_weights_unbounded <- function(spec, frame, prior, exact = FALSE,
                                    route = c("map", "sample")) {
  route <- match.arg(route)
  if (length(spec$responses) != 1L) return(invisible())
  r <- names(spec$responses)[1L]
  resp <- spec$responses[[r]]
  fam <- resp$family[["family"]] %||% ""
  spike <- unbounded_spike_dpar[[fam]]
  if (is.null(spike)) return(invisible())
  laplace_linear <- fam %in% c("asym_laplace", "zero_inflated_asym_laplace") &&
    !isTRUE(exact)
  if (!is.null(resp$aterms[["se"]]) ||
        length(frame[["autocor"]] %||% list()) ||
        length(frame[["mi_map"]] %||% list())) {
    return(invisible())
  }
  n <- frame[["n_obs"]]
  w <- frame[["aterm_values"]][[r]][["weights"]] %||% rep(1, n)
  lps <- lapply(spike, function(dp) frame[["linpreds"]][[linpred_key(r, dp)]])
  if (!all(vapply(lps, unbounded_spans_constant, NA))) return(invisible())
  lp <- frame[["linpreds"]][[linpred_key(r, "mu")]]
  Z <- lp[["Z"]]
  if (is.null(Z) || !is.null(lp[["nl_body"]]) || length(lp[["Z_mo"]])) {
    return(invisible())
  }
  owned_blocks <- Filter(cs_full_rank, frame[["re_blocks"]])
  cols <- unlist(lapply(owned_blocks, `[[`, "c_idx"))
  cols <- intersect(cols, seq_len(ncol(Z)))
  if (!length(cols)) return(invisible())
  Zc <- Matrix::drop0(methods::as(Z[, cols, drop = FALSE], "CsparseMatrix"))
  one <- which(diff(Zc@p) == 1L)
  owned <- unique(Zc@i[Zc@p[one] + 1L] + 1L)
  pos <- which(w > 0)
  if (!length(pos) || !all(pos %in% owned)) return(invisible())
  # the excess is the uncensored rows' (see above); a censored row needs
  # its own effect too, or its probability can fall to 0 as the
  # dispersion goes to its limit
  cen <- frame[["aterm_values"]][[r]][["cens"]]
  unc <- if (is.null(cen)) pos else pos[rep_len(cen, n)[pos] == 0]
  n_cens <- length(pos) - length(unc)
  # a zero of zero_inflated_asym_laplace is the point mass: log(zi),
  # whatever its own effect, so the row neither grows nor curves
  if (identical(fam, "zero_inflated_asym_laplace")) {
    yv <- frame[["y"]][[r]]
    unc <- unc[yv[unc] != 0]
  }
  if (!length(unc)) return(invisible())
  excess <- if (laplace_linear) sum(w[unc]) else sum(w[unc] - 1)
  if (!isTRUE(exact)) excess <- excess - n_cens / 2
  # rounding of weights normalized to mean 1 must not read as an excess
  if (excess <= 1e-8 * length(unc)) return(invisible())
  verdicts <- vapply(seq_along(spike), function(k) {
    prior_spike_verdict(prior, r, spike[k], lps[[k]], excess, route)
  }, "")
  # one bounded spike dpar keeps the densities from narrowing at all
  if (any(verdicts != "free")) return(invisible())
  terms <- unique(vapply(Filter(function(bk) {
    any(bk[["c_idx"]] %in% cols[one])
  }, owned_blocks), block_term_as_written, ""))
  frm_stop(weights_unbounded_message(terms, w[unc], excess, spike, fam,
                                     laplace_linear, route, n_cens),
           call. = FALSE)
}

#' A block's term as the formula spells it. block_written_label() falls
#' back to "tt + 0 | g [ar1]" for a covariance-structure wrapper, which
#' is not what the user wrote.
#'
#' @noRd
block_term_as_written <- function(bk) {
  wl <- block_written_label(bk)
  m <- regmatches(wl, regexec("^(.*) \\[([A-Za-z_0-9.]+)\\]$", wl))[[1]]
  if (length(m) != 3L) return(wl)
  # a dpar's block is labeled "sigma: 1 | g"; the prefix stays outside
  # the term, as the formula's sigma ~ (1 | g) has it
  pre <- regmatches(m[2], regexec("^([A-Za-z_.0-9]+): (.*)$", m[2]))[[1]]
  dp <- if (length(pre) == 3L) paste0(pre[2], ": ") else ""
  body <- if (length(pre) == 3L) pre[3] else m[2]
  # a plain bar term is stored as a us() block, and a Student-t one,
  # written with gr(dist = "student"), as us_t() or diag_t()
  if (m[3] %in% c("us_t", "diag_t")) {
    lr <- strsplit(body, " | ", fixed = TRUE)[[1]]
    if (length(lr) == 2L) {
      return(paste0(dp, "(", lr[1], if (m[3] == "diag_t") " || " else " | ",
                    "gr(", lr[2], ", dist = \"student\"))"))
    }
  }
  paste0(dp, if (identical(m[3], "us")) paste0("(", body, ")") else
    paste0(m[3], "(", body, ")"))
}

#' Whether a dpar's predictor can reach its limit at every row at once:
#' its fixed-effect design has a column and spans the constant. An
#' intercept-free `sigma ~ 0 + x` and an offset-only `sigma ~ 0 +
#' offset(o)` cannot, and their likelihood is bounded.
#'
#' @noRd
unbounded_spans_constant <- function(lp) {
  if (is.null(lp) || !is.null(lp[["constant"]]) ||
        !is.null(lp[["nl_body"]])) {
    return(FALSE)
  }
  X <- lp[["X"]]
  if (is.null(X) || !ncol(X)) return(FALSE)
  X <- as.matrix(X)
  one <- rep(1, nrow(X))
  res <- qr.resid(qr(X), one)
  sqrt(sum(res^2)) <= 1e-8 * sqrt(nrow(X))
}

#' What the prior on one spike dpar does to the growth: "bounds" when it
#' provably stops it (MAP) or makes the posterior integrate (sample),
#' "free" when it provably does not, "unknown" otherwise. No entry on
#' the dpar is "free": the objective is then the likelihood (MAP) or the
#' sampler's flat density on the link scale (sample).
#'
#' @noRd
prior_spike_verdict <- function(prior, r, dpar, lp, excess, route) {
  if (!length(prior)) return("free")
  if (!inherits(prior, "frmtmb_priorlist")) return("unknown")
  zero <- dpar %in% unbounded_spike_toward_zero
  link <- lp[["link"]][["name"]] %||% ""
  mine <- Filter(function(s) {
    identical(s$dpar %||% "", dpar) &&
      (!nzchar(s$resp %||% "") || identical(s$resp, r))
  }, prior)
  X <- as.matrix(lp[["X"]])
  has_icpt <- any(apply(X, 2L, function(v) v[1L] != 0 && all(v == v[1L])))
  if (!has_icpt) {
    return(prior_cellmeans_verdict(mine, X, zero, link, excess, route))
  }
  out <- character(0)
  for (s in mine) {
    if (!identical(s$class, "Intercept") || nzchar(s$coef %||% "")) next
    out <- c(out, prior_entry_verdict(s, zero, link, excess, route))
  }
  if (!length(out)) return("free")
  if (any(out == "bounds")) return("bounds")
  if (any(out == "unknown")) return("unknown")
  "free"
}

#' The verdict for a dispersion design that spans the constant without
#' an intercept column, as `sigma ~ 0 + fac` does through its cell
#' means. Every column then has to move to reach the limit at every
#' row, so a class `b` prior that bounds each column's coefficient on
#' the link scale bounds the model as an Intercept prior does. Some
#' columns bounded and some not is "unknown"; no `b` entry at all is
#' "free".
#'
#' @noRd
prior_cellmeans_verdict <- function(mine, X, zero, link, excess, route) {
  b <- Filter(function(s) identical(s$class, "b"), mine)
  if (!length(b)) return("free")
  cns <- colnames(X) %||% character(ncol(X))
  per_col <- vapply(cns, function(cn) {
    hit <- Filter(function(s) identical(s$coef %||% "", cn), b)
    if (!length(hit)) hit <- Filter(function(s) !nzchar(s$coef %||% ""), b)
    if (!length(hit)) return("none")
    # a coefficient's prior is on the link scale
    s <- hit[[1L]]
    s$natural <- FALSE
    prior_entry_verdict(s, zero, link, excess, route)
  }, "")
  if (all(per_col == "bounds")) "bounds" else "unknown"
}

#' The verdict of one prior entry on a spike dpar's intercept.
#'
#' On the natural scale, toward 0, the density goes like x^order and the
#' likelihood like x^(-excess): the MAP objective is bounded when order
#' exceeds the excess, and the posterior integrates when order - excess
#' is above -1. Toward infinity the likelihood grows like a power of x,
#' which any density with a faster than polynomial tail stops. On the
#' log link the likelihood grows like exp(excess |eta|) or slower, which
#' a normal tail stops and a student_t tail does not.
#'
#' @noRd
prior_entry_verdict <- function(s, zero, link, excess, route) {
  d <- s$dist
  kind <- d$kind %||% ""
  lb <- suppressWarnings(as.numeric(s$lb))
  ub <- suppressWarnings(as.numeric(s$ub))
  if (isTRUE(s$natural)) {
    if (zero && isTRUE(lb > 0)) return("bounds")
    if (!zero && isTRUE(is.finite(ub))) return("bounds")
    if (!zero) {
      return(if (kind %in% c("normal", "gamma", "exponential", "logistic",
                             "lognormal", "beta", "uniform")) {
        "bounds"
      } else "unknown")
    }
    # brms's (flat), frmtmb.sample's prior = "flat", is order 0 at 0
    order <- switch(kind,
      normal = , t = , logistic = , exponential = , flat = 0,
      uniform = if (isTRUE(d$lower > 0)) Inf else 0,
      gamma = d$shape - 1,
      beta = d$shape1 - 1,
      inv_gamma = , lognormal = Inf,
      NA_real_)
    if (is.na(order)) return("unknown")
    if (route == "map") {
      if (order > excess) "bounds" else if (order < excess) "free" else
        "unknown"
    } else {
      if (order - excess > -1) "bounds" else "free"
    }
  } else {
    if (!identical(link, "log")) return("unknown")
    if (zero && isTRUE(is.finite(lb))) return("bounds")
    if (!zero && isTRUE(is.finite(ub))) return("bounds")
    switch(kind,
      normal = "bounds",
      t = , flat = "free",
      uniform = "bounds",
      # a density on (0, Inf) never reaches eta = -Inf
      exponential = , gamma = , inv_gamma = , lognormal = , beta =
        if (zero) "bounds" else if (kind %in% c("lognormal", "beta"))
          "bounds" else if (kind == "inv_gamma") "free" else "unknown",
      logistic = if (1 / d$scale > excess) "bounds" else "unknown",
      "unknown")
  }
}

#' The refusal's text, as the model's own terms and numbers.
#'
#' @noRd
weights_unbounded_message <- function(terms, w, excess, spike, fam,
                                      laplace_linear, route, n_cens = 0L) {
  zero <- spike[1L] %in% unbounded_spike_toward_zero
  dp <- paste(spike, collapse = " and ")
  lim <- if (zero) {
    paste(dp, if (length(spike) > 1L) "shrink" else "shrinks", "to 0")
  } else {
    paste(dp, "grows")
  }
  head <- if (route == "sample") {
    "The posterior of this model is improper. "
  } else if (laplace_linear) {
    "The Laplace approximation of this model's likelihood has no maximum. "
  } else {
    "The weighted likelihood of this model has no maximum. "
  }
  owns <- paste0(paste0("`", terms, "`", collapse = ", "), " gives every ",
                 "row with a positive weight its own random effect in mu",
                 if (length(w) && any(w != 1)) {
                   paste0(", and the weights sum to ",
                          format(sum(w), digits = 6), " over ", length(w),
                          if (n_cens) " uncensored", " rows")
                 }, ". ",
                 if (n_cens) {
                   paste0("The ", n_cens, " censored rows do not stop ",
                          "it: a censored row's probability is at most 1, ",
                          "and its own effect can bring it to 1",
                          if (!laplace_linear) {
                            paste("; under the Laplace approximation",
                                  "each takes 1/2 off the power below")
                          }, ". ")
                 })
  why <- if (laplace_linear) {
    paste0("The ", fam, " log density is linear on each side of its ",
           "mode, so the data give a row's own effect no curvature, and as ",
           lim, " the approximation grows like the peak height to the ",
           "power sum(w) = ", format(excess, digits = 4), ", whatever the ",
           "weights. quadrature = TRUE integrates the effect exactly and ",
           "has a maximum when the weights sum to less than the number of ",
           "rows; so does the model without the observation-level term.")
  } else {
    paste0("A weight above 1 counts a row more than once, so the row's ",
           "own effect then stands for replicates with no residual ",
           "variance: as ", lim, ", each density narrows onto its own ",
           "observation, and the likelihood grows like the height of ",
           "that peak to the power sum(w) - n = ", format(excess, digits = 4),
           ", without bound. A usable fit needs weights that sum to the ",
           "number of rows or fewer, a random effect that is not ",
           "observation-level, or a fixed ", spike[1L], ".")
  }
  # a prior strong enough to bound the objective puts the estimate where
  # it and the excess balance: lognormal(0, 1) on sigma stopped at a log
  # sigma of -19.8 and sampled with Rhat 2.83
  pri <- if (laplace_linear) "" else {
    paste0(" A prior on ", spike[1L], " can make the objective bounded, ",
           "but not the estimate informative: the mode then sits where ",
           "the prior and the weight excess balance.")
  }
  paste0(head, owns, why, pri)
}
