# The joint covariance of a grid prediction, taken from core's seam.
#
# A simultaneous band, a delta-method derivative and an implicit-function
# feature all need the same object: the covariance of the WHOLE grid
# prediction, Sigma = A V A', where A maps the coefficient vector to the
# curve and V is the joint covariance of the fixed and random
# coefficients. A smooth's wiggly part is a random-effect block even when
# the smooth is a population term, so V must carry the b block.
# `vcov(fit, full = TRUE)` does not carry it under either of its
# branches, and `frm_linpred(se.fit = TRUE)` forms Sigma internally and
# returns only sqrt(diag(Sigma)).
#
# Until frmtmb 0.52.0 this package rebuilt A by unit perturbation, one
# predict() call per contributing coefficient, and read V out of
# `fit$cache$Vjoint`, which was an internal. Core now exports both
# halves: `frm_lp_basis()` returns A, its coefficient positions, the
# covariance at those positions and the variance that is NOT coefficient
# uncertainty, and `frm_joint_cov()` returns V on its own. Everything
# below is a consumer of the first of those.
#
# The check survives the change and is still not decoration:
# `diag(A V A') + extra_var` must equal `frm_linpred(se.fit = TRUE)^2`
# wherever core has an independent route to that number. What the check
# now verifies is that this package reads the seam correctly, rather
# than that a reconstruction reproduced it.
#
# The variance that is NOT coefficient uncertainty comes from the same
# call as a full covariance, `extra_cov` (frmtmb's `frm_lp_basis(extra_cov
# = TRUE)`): an exact `gp()`'s kriging residual, which every unseen row
# of one field shares, and a new grouping level's draw, which every row
# at that level shares. The grid's covariance is `A V A' + extra_cov`,
# and the band, the derivative and the feature are all built from it.
# Before that seam existed this package drew the band from `A V A'` and
# divided by a standard error that carried the residual, which made a
# band over an exact `gp()` past the observed positions at least 17
# percent too narrow (dev/diffcurve-findings.md).
#
# A DIFFERENCE CURVE is the same object twice, read in ONE call on the
# stacked grid `rbind(newdata, contrast)`: the off-diagonal block of
# `extra_cov` is the covariance between the two grids, so two grids at
# the same `gp()` position cancel their residual and two at different
# positions carry its difference, with no predicate deciding which. What
# the seam does NOT hand over is a second route to the difference's
# standard error: `frm_linpred(se.fit = TRUE)` returns a marginal
# standard error per row. So the check is run on each half and
# `cov_rel_error` is the worse of the two, and what it licenses is that
# both designs were read correctly, not that the difference's own
# standard error was verified against anything.

#' Run `expr`, holding back any `ps()` knot-span warning frmtmb raises
#' inside it.
#'
#' The seam is a CLASS, `frmtmb_ps_span_warning`, not the text: the
#' three exported functions here evaluate the curve at three, five or
#' one grid position per point the user asked about, so core's warning
#' would arrive once per internal evaluation and count stencil rows
#' rather than grid rows. Held back and returned, it can be surfaced
#' once under the name of the function the user actually called, or
#' turned into a refusal where a warning is the wrong answer.
#'
#' @noRd
sp_catch_span <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    expr,
    frmtmb_ps_span_warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = val, span = unique(msgs))
}

#' Refuse a feature search whose bracket leaves a `ps()` knot span.
#'
#' A warning is the right answer for a band, which the reader can see
#' bending to the intercept, and the wrong one for a feature: a peak or
#' a crossing located past the span is a peak or a crossing OF THE
#' DECAYING PARTIAL SUM, the search reports it as a number with a
#' standard error, and nothing in the output says which curve it belongs
#' to. One template, called from both places the bracket is checked.
#'
#' @noRd
sp_span_stop <- function(span) {
  if (!length(span)) return(invisible(NULL))
  frm_stop("frm_curve_feature(): the search bracket leaves a ps() term's ",
           "knot span, so a root located in it would be a root of the ",
           "decaying partial sum rather than of the fitted curve, and the ",
           "implicit-function standard error beside it would describe ",
           "neither. Narrow newdata to the span. ",
           paste(span, collapse = " "), call. = FALSE)
}

#' The `ps()` span messages for the grid the USER passed, rather than
#' for the stencil the caller widened it into.
#'
#' [frm_curve_deriv()] evaluates the design at `c(x - e, x, x + e)` and
#' [frm_curve_feature()]'s stationary-point scan at `c(x - e1, x + e1)`,
#' with `e` a millionth of the grid's range. Core is handed those points
#' and counts those rows, so a grid whose endpoint sits exactly on a
#' knot is outside the span by `e` and gets a warning, or in the feature
#' case a refusal, for a bracket the user chose entirely inside it.
#' Following the refusal's own advice ("Narrow newdata to the span")
#' reproduced the refusal.
#'
#' One `predict()` on the grid itself settles it. Callers only reach
#' here when the widened call already reported something, which for the
#' derivative stencil is exact (its point set CONTAINS the grid, so a
#' clean stencil implies a clean grid) and for the feature scan holds
#' whenever the spline argument is monotone in `var`, which is every
#' shape this package documents.
#'
#' @noRd
sp_span_on_grid <- function(fit, nd, dpar, resp, re_formula,
                            allow_new_levels) {
  sp_catch_span(sp_predict_eta(fit, nd, dpar, resp, re_formula,
                               allow_new_levels))$span
}

#' The same question over BOTH grids of a difference curve, with each
#' message saying which grid raised it.
#'
#' A difference is built from two frames, so every span statement this
#' package makes has two grids to make it about, and a reader who is
#' told "a value is outside the frozen knot span" without being told
#' which frame holds it has to guess. `contrast = NULL` is the ordinary
#' one-grid case and costs exactly what it did.
#'
#' The tag is pasted here rather than inside the `stop()` or `warning()`
#' that carries it, so the two doors keep one message template each.
#'
#' @noRd
sp_span_both <- function(span_a, span_b) {
  unique(c(span_a,
           if (length(span_b)) paste0("In `contrast`: ", span_b)))
}

#' The grid span messages over whichever grids the search actually
#' holds: one for a curve, both for a difference.
#'
#' @noRd
sp_grid_span <- function(sp, nd, ct) {
  sp_span_both(
    sp_span_on_grid(sp$fit, nd, sp$dpar, sp$resp, sp$re_formula,
                    sp$allow_new_levels),
    if (is.null(ct)) character(0) else
      sp_span_on_grid(sp$fit, ct, sp$dpar, sp$resp, sp$re_formula,
                      sp$allow_new_levels))
}

#' Does the grid hold every column but `var` at row 1's value?
#'
#' [frm_curve_feature()] evaluates nothing but row 1 with `var` moved:
#' the scan, the Newton steps and the five-point stencil are all
#' `row1[rep(1L, n), ]` with that one column overwritten. So when the
#' rest of the grid is pinned to row 1, the scan has already predicted
#' at every covariate combination the grid holds, and its own span
#' check covers the grid. Only a column that VARIES down the grid can
#' put a value outside a second `ps()` term's span in a row the search
#' never evaluates, which is the one case the second check is for.
#'
#' `duplicated()` rather than a comparison: it compares factors,
#' characters and matrix columns as themselves, it treats `NA` as equal
#' to `NA`, and "duplicated from row 2 on" is exactly "one distinct
#' value". A numeric tolerance would be wrong here, because a value a
#' hair past a knot is outside the span.
#'
#' @noRd
sp_grid_pinned <- function(nd, var) {
  if (nrow(nd) < 2L) return(TRUE)
  for (nm in setdiff(names(nd), var)) {
    if (!all(duplicated(nd[[nm]])[-1L])) return(FALSE)
  }
  TRUE
}

#' One prediction on the link scale, as a plain numeric vector.
#'
#' `allow_new_levels` has no default here or in any helper below: the
#' design, the estimate and the covariance check must all read the grid
#' under ONE setting, and a helper that silently fell back to `FALSE`
#' would refuse the unseen level on a route the caller never chose.
#'
#' @noRd
sp_predict_eta <- function(fit, newdata, dpar, resp, re_formula,
                           allow_new_levels) {
  as.numeric(frmtmb::frm_linpred(fit, newdata = newdata, type = "link",
                                 dpar = dpar, resp = resp,
                                 re_formula = re_formula,
                                 allow_new_levels = allow_new_levels))
}

#' The seam read at ONE grid: the design, its covariance and the
#' standard error the check compares.
#'
#' `E` is the covariance of the variance that is not coefficient
#' uncertainty, as core returns it (sparse when its blocks cover little
#' of the grid, all zero when nothing contributes), and `Sigma` the
#' grid's whole covariance `A V A' + E`, dense. A difference reads the
#' stacked grid through here once, so both halves come off the same seam
#' at the same coefficient rows.
#'
#' The standard error is `diag(A V A')` and NOT the `rowSums((A %*% V)
#' * A)` core writes it with, even though a difference discards both
#' halves' full covariances and keeps only its own.
#'
#' What that buys is one ulp, and it is worth being exact about how
#' little that is. Both spellings read the same `A`, `V` and
#' `extra_var` out of the same `frm_lp_basis()` call, so their only
#' independence is summation order. A real misreading of the seam is
#' caught either way: with `V`'s first two rows and columns swapped both
#' spellings refuse at 0.101 relative. What core's spelling costs is the
#' ability of `tol` to mean anything at zero. Measured: it reports 0
#' exactly, and `test-curve.R`'s "a tolerance no covariance could meet
#' refuses rather than returning" stops refusing (PASS=61 FAIL=1). That
#' is the whole reason for the triple product, and the extra work is
#' paid on a call whose cost is dominated by one joint-precision solve.
#'
#' @noRd
sp_one_basis <- function(fit, nd, dpar, resp, re_formula, allow_new_levels) {
  lbc <- sp_catch_span(frmtmb::frm_lp_basis(
    fit, newdata = nd, dpar = dpar, resp = resp, re_formula = re_formula,
    allow_new_levels = allow_new_levels, extra_cov = TRUE))
  lb <- lbc$value
  C <- as.matrix(lb$A)
  Sigma <- unname(C %*% lb$V %*% t(C))
  se <- unname(sqrt(pmax(diag(Sigma) + lb$extra_var, 0)))
  E <- lb$extra_cov
  Sigma <- sp_add_extra(Sigma, E)
  list(lb = lb, C = C, E = E, Sigma = Sigma, se = se, span = lbc$span)
}

#' `S + E`, without a dense copy of `E` when it has no nonzero entry,
#' which is every fit with no exact `gp()` off its positions and no
#' unseen level.
#'
#' @noRd
sp_add_extra <- function(S, E) {
  if (is.null(E)) return(S)
  if (inherits(E, "Matrix")) {
    if (!length(E@x) || !any(E@x != 0)) return(S)
    E <- as.matrix(E)
  }
  S + unname(E)
}

#' `sqrt(diag(A V A') + extra_var)` against `frm_linpred(se.fit = TRUE)`, or
#' a refusal.
#'
#' One template, called once for an ordinary curve and twice for a
#' difference; `side` names the grid at run time so that the two calls
#' still resolve to one line of source.
#'
#' @noRd
sp_cov_check <- function(fit, nd, se, dpar, resp, re_formula,
                         allow_new_levels, tol, side) {
  ref <- frmtmb::frm_linpred(fit, newdata = nd, type = "link", dpar = dpar,
                             resp = resp, re_formula = re_formula,
                             allow_new_levels = allow_new_levels,
                             se.fit = TRUE)
  se_ref <- as.numeric(ref$se.fit)
  rel <- max(abs(se / pmax(se_ref, .Machine$double.eps) - 1))
  if (!is.finite(rel) || rel > tol) {
    frm_stop("frm_curve(): the assembled covariance of ", side,
             " disagrees with frm_linpred(se.fit = TRUE) by ",
             format(rel, digits = 3),
             " relative, which is above the tolerance ", format(tol),
             ". Both come from frm_lp_basis(); a disagreement means this ",
             "package is reading the seam wrongly, and a fit where the two ",
             "disagree is one it must not report a band for", call. = FALSE)
  }
  rel
}

#' Everything the three exported functions share: the grid, the design,
#' the covariance, and the check that the covariance is the right one.
#'
#' `frm_lp_basis(extra_cov = TRUE)` (frmtmb's development version)
#' returns `A`, `V` and the variance that is not coefficient uncertainty
#' as a full covariance over the rows. For a linear predictor it is the
#' same object `frm_linpred(se.fit = TRUE)` reduces to a diagonal, so the
#' two are compared on every call. For a NONLINEAR body `A` is a
#' Jacobian, `frm_linpred(se.fit = TRUE)` refuses outright, and there is
#' nothing to compare against: the check is skipped and `rel` is `NA`,
#' which `print()` reports rather than hides.
#'
#' With `contrast`, the reported functional is `(A1 - A2) c` and its
#' covariance is `(A1 - A2) V (A1 - A2)' + E11 + E22 - E12 - E21`, with
#' `E` the extra covariance of the stacked grid `rbind(newdata,
#' contrast)`, read in ONE call so that its off-diagonal block is the
#' covariance between the two grids. Two grids at one exact `gp()`
#' position load the same kriging residual, and it cancels; two at
#' different positions load correlated ones, and their difference is
#' what is left; a grid at an unseen grouping level and one at the SAME
#' unseen level load one draw of its effect, and two different unseen
#' levels independent draws. The check compares each half against
#' `frm_linpred(se.fit = TRUE)`: there is no second route to the
#' difference's own standard error.
#'
#' @noRd
sp_curve_parts <- function(fit, newdata, dpar, resp, re_formula,
                           allow_new_levels, tol, contrast = NULL) {
  if (!inherits(fit, "frmtmb_fit")) {
    frm_stop("frm_curve(): `object` must be a frmtmb fit, the model a curve ",
             "is read off, not an object of class ", class(fit)[1L],
             call. = FALSE)
  }
  if (!is.data.frame(newdata) || !nrow(newdata)) {
    frm_stop("`newdata` must be a data frame with at least one row: it is ",
             "the grid the curve is evaluated on", call. = FALSE)
  }
  nl <- sp_is_nl(fit, dpar, resp)
  if (is.null(contrast)) {
    a <- sp_one_basis(fit, newdata, dpar, resp, re_formula,
                      allow_new_levels)
    out <- list(eta = a$lb$eta, C = a$C, V = a$lb$V, Sigma = a$Sigma,
                E = a$E, se = a$se, rel = NA_real_,
                n_predict = 0L, newdata = newdata, contrast = contrast,
                dpar = dpar, resp = resp, re_formula = re_formula,
                allow_new_levels = allow_new_levels, fit = fit,
                span = a$span, extra_var = a$lb$extra_var)
    # A nonlinear body is the case core refuses se.fit for, so there is
    # no second number to check against. Everything else is checked.
    if (!nl) {
      out$rel <- sp_cov_check(fit, newdata, a$se, dpar, resp, re_formula,
                              allow_new_levels, tol, "this grid")
      out$n_predict <- 1L
    }
    return(out)
  }
  n <- nrow(newdata)
  i1 <- seq_len(n)
  i2 <- n + i1
  st <- sp_one_basis(fit, sp_stack(newdata, contrast), dpar, resp,
                     re_formula, allow_new_levels)
  dif <- function(M) M[i1, i1] + M[i2, i2] - M[i1, i2] - M[i2, i1]
  C <- st$C[i1, , drop = FALSE] - st$C[i2, , drop = FALSE]
  E <- dif(st$E)
  Sigma <- sp_add_extra(unname(C %*% st$lb$V %*% t(C)), E)
  # the stacked call counts the rows of both grids together, so a span
  # message is re-asked of each grid to say which one holds the value
  span <- if (length(st$span)) {
    sp_span_both(
      sp_span_on_grid(fit, newdata, dpar, resp, re_formula,
                      allow_new_levels),
      sp_span_on_grid(fit, contrast, dpar, resp, re_formula,
                      allow_new_levels))
  } else {
    character(0)
  }
  out <- list(eta = st$lb$eta[i1] - st$lb$eta[i2], C = C, V = st$lb$V,
              Sigma = Sigma, E = E,
              se = unname(sqrt(pmax(diag(Sigma), 0))), rel = NA_real_,
              n_predict = 0L, newdata = newdata, contrast = contrast,
              dpar = dpar, resp = resp, re_formula = re_formula,
              allow_new_levels = allow_new_levels, fit = fit, span = span,
              extra_var = as.numeric(E[cbind(i1, i1)]))
  if (!nl) {
    out$rel <- max(
      sp_cov_check(fit, newdata, st$se[i1], dpar, resp, re_formula,
                   allow_new_levels, tol, "`newdata`"),
      sp_cov_check(fit, contrast, st$se[i2], dpar, resp, re_formula,
                   allow_new_levels, tol, "`contrast`"))
    out$n_predict <- 2L
  }
  out
}

#' The two grids of a difference as one frame, `newdata`'s rows first.
#'
#' Only the columns both carry: a column one grid lacks is one the
#' prediction does not read, or the grid that lacks it is refused by
#' core under its own name either way.
#'
#' @noRd
sp_stack <- function(a, b) {
  nm <- intersect(names(a), names(b))
  out <- rbind(a[nm], b[nm])
  rownames(out) <- NULL
  out
}

#' Is this linear predictor computed by a nonlinear body?
#'
#' @noRd
sp_is_nl <- function(fit, dpar, resp) {
  rn <- resp %||% names(fit$spec$responses)[1L]
  rspec <- fit$spec$responses[[rn]]
  dp <- dpar %||% if ("mu" %in% names(rspec$dpars)) "mu" else
    rspec$primary_dpars[1L]
  lp <- fit$frame[["linpreds"]][[paste0(rn, ".", dp)]]
  !is.null(lp) && !is.null(lp[["nl_body"]])
}

#' The max-deviation simultaneous critical value of Ruppert, Wand and
#' Carroll (2003, ch. 6).
#'
#' Draw the curve's own deviation process, standardize it by `div`, take
#' the largest absolute value over the grid, and report the `level`
#' quantile. `div` is a separate argument rather than `sqrt(diag(S))`
#' because the critical value is comparable between two packages only
#' when the divisor is: gratia standardizes a smooth-only deviation by
#' the FULL predictor's standard error, and the critical value that
#' comes back is 8 percent smaller than the self-standardized one on the
#' same fit. The BAND is right either way, because the same divisor
#' calibrates it and scales it.
#'
#' @section A divisor of exactly zero:
#' A grid point whose standard error is exactly zero has a
#' DETERMINISTIC deviation, not a small one. `S` is positive
#' semi-definite, so a zero diagonal entry forces that whole row to zero,
#' and `(L z)` is exactly 0 there: the standardized deviation is `0 / 0`.
#' Such a point is covered with probability one, so it cannot be the
#' argmax and it leaves the maximization rather than turning it into
#' `NaN` and killing `quantile()` one line later.
#'
#' This is not a corner case. Past a [frmtmb::ps()] term's knot span the
#' basis is exactly zero, so its DERIVATIVE design is exactly zero and
#' every such row of a [frm_curve_deriv()] grid has a standard error of
#' exactly zero. Before this guard, `frm_curve_deriv()` on its own
#' defaults raised the span warning and then died in `quantile.default`
#' with a message naming neither the span nor the function.
#'
#' Dropping nothing is bit-identical to not having the guard, so a grid
#' with no zero divisor is unaffected.
#'
#' @noRd
sp_sim_crit <- function(S, div, nsim, level, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- nrow(S)
  keep <- is.finite(div) & div > 0
  if (!any(keep)) {
    frm_stop("simultaneous = TRUE: every point on this grid has a standard ",
             "error of exactly zero, so the deviation process is degenerate ",
             "and there is no maximum to take a quantile of. Past a ps() ",
             "term's knot span the basis is exactly zero and so is the ",
             "derivative design, which is the usual way to arrive here. Use ",
             "simultaneous = FALSE, or move the grid inside the span",
             call. = FALSE)
  }
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = m)
  z <- matrix(stats::rnorm(m * nsim), m, nsim)
  mx <- apply(abs((L[keep, , drop = FALSE] %*% z) / div[keep]), 2L, max)
  crit <- unname(stats::quantile(mx, level, type = 8))
  # The standard error of a sample quantile, sqrt(p(1-p)/n) / f(q). A
  # critical value reported without it invites a comparison across
  # packages that simulation noise alone would fail.
  d <- stats::density(mx)
  f <- stats::approx(d$x, d$y, xout = crit)$y
  list(crit = crit,
       mcse = if (isTRUE(is.finite(f) && f > 0)) {
         sqrt(level * (1 - level) / nsim) / f
       } else NA_real_)
}

#' The covariance of the curve's derivative that is not coefficient
#' uncertainty, from frmtmb's `frm_extra_cov_deriv()`, which forms it
#' from its sources rather than by differencing `extra_cov`; for a
#' difference, from one call on the stacked grid. `NULL` when nothing
#' contributes, and for a nonlinear body, whose route refuses both
#' sources.
#'
#' @noRd
sp_extra_deriv <- function(sp, nd, ct, var, order, e) {
  if (sp_is_nl(sp$fit, sp$dpar, sp$resp)) return(NULL)
  one <- function(g) {
    frmtmb::frm_extra_cov_deriv(sp$fit, g, var = var, order = order,
                                eps = e, dpar = sp$dpar, resp = sp$resp,
                                re_formula = sp$re_formula,
                                allow_new_levels = sp$allow_new_levels)
  }
  E <- if (is.null(ct)) {
    one(nd)
  } else {
    n <- nrow(nd)
    i1 <- seq_len(n)
    i2 <- n + i1
    E2 <- one(sp_stack(nd, ct))
    if (inherits(E2, "Matrix") && !length(E2@x)) return(NULL)
    E2[i1, i1] + E2[i2, i2] - E2[i1, i2] - E2[i2, i1]
  }
  # an empty sparse matrix is the no-source case, and costs nothing
  if (inherits(E, "Matrix")) {
    if (!length(E@x) || !any(E@x != 0)) return(NULL)
    E <- as.matrix(E)
  }
  if (!any(E != 0)) NULL else unname(E)
}
