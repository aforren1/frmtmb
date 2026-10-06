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
# A DIFFERENCE CURVE is the same object twice. `contrast = ` builds the
# seam at a second grid, subtracts the two designs and reports
# `(A1 - A2) V (A1 - A2)'`. What the seam does NOT hand over is a second
# route to that number: `frm_linpred(se.fit = TRUE)` returns a marginal
# standard error per row and never the covariance BETWEEN the two grids,
# which is the whole content of a difference. So the check is run on
# each half and `cov_rel_error` is the worse of the two, and what it
# licenses is that both designs were read correctly, not that the
# difference's own standard error was verified against anything.
#
# `extra_var` meets the same limit from the other side. An exact `gp()`
# at an unseen position contributes a kriging residual per ROW, and the
# seam returns no covariance between the two grids, so `var(g1 - g2)`
# has no public route. It does not need one when the two grids load the
# same residual, which they do whenever their non-fixed design columns
# agree, and then the term is exactly zero rather than unknown. That is
# the ordinary case and it is now computed; the refusal is what is left.
# `sp_same_latent()` carries the test and the measurements behind it.

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
#' Split out because a difference curve needs it twice and the halves
#' must come off the same seam at the same coefficient rows before the
#' subtraction means anything.
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
    allow_new_levels = allow_new_levels))
  lb <- lbc$value
  C <- as.matrix(lb$A)
  Sigma <- unname(C %*% lb$V %*% t(C))
  # a row that moves along a direction the fit's Hessian lost has no
  # standard error; `V` leaves those directions out, so its finite value
  # there would be too small (frmtmb >= 0.68.0 marks the rows)
  lost <- lb$se_nonest %||% rep(FALSE, nrow(C))
  se <- unname(sqrt(pmax(diag(Sigma) + lb$extra_var, 0)))
  se[lost] <- NaN
  list(lb = lb, C = C, Sigma = Sigma, se = se, lost = lost,
       span = lbc$span)
}

#' Do the two grids of a difference load the SAME latent draw?
#'
#' Variance that is not coefficient uncertainty, which for a curve is an
#' exact `gp()` kriging residual, arrives from the seam as one
#' number per ROW with no covariance between the two grids, so a
#' difference cannot form `var(g1 - g2)` in general. It does not have to
#' when `g1` and `g2` are the SAME random variable: the term is then
#' exactly zero and the difference is `(A1 - A2) V (A1 - A2)'` with
#' nothing left over.
#'
#' The test is on the DESIGN and not on the numbers. A latent block
#' reaches a prediction only through its own columns of `A`, so two
#' grids whose non-fixed columns are bit-identical load one functional
#' of one field, and whatever that functional leaves unexplained is one
#' residual rather than two. Equality of `extra_var` alone would NOT do:
#' two different levels of one grouping block have identical marginal
#' variances by construction and are different draws. Measured on
#' `y ~ fac + s(x, k = 6) + (1 | g)` at `re_formula = NULL`, contrasting
#' level 1 against level 2: `extra_var` identical, non-fixed design
#' columns not identical, so this returns `FALSE` where the numbers
#' agree.
#'
#' Sensitivity, on `y ~ fac + gp(x)` at n = 90 with the grid off the
#' observed positions: a contrast across `fac` at one `x` gives `TRUE`,
#' and moving the second grid's `x` by 1e-10, or mirroring it about the
#' range midpoint, gives `FALSE`.
#'
#' It is STRICTER than the mathematics needs, and the refusal says so.
#' Only the block that carries `extra_var` has to match for the residual
#' to cancel; this asks it of every non-fixed column. So a contrast
#' across `fac` on `y ~ fac + s(x, by = fac) + gp(x)` is refused even
#' though the `gp()` columns are bit-identical, because the by-factor
#' smooth's own columns differ, which is what a by-factor smooth is for.
#' Narrowing to the right block needs the block identity, and the only
#' public route to it is the `b.<block>.<level>` shape of
#' `frm_joint_cov()$labels`, which is name parsing this package does not
#' do. It fails closed, so the cost is a refused answer and not a wrong
#' one, and the core seam that returns a matrix-valued `extra_cov`
#' removes the need for the predicate entirely.
#'
#' What it assumes: for an exact `gp()` the kriging weights fix the
#' position uniquely unless every observed position is equidistant from
#' the two grids' positions. In one dimension that needs every
#' observation at one point, which makes the kernel matrix singular; in
#' two it needs them collinear with the two positions mirrored across
#' that line. That design was built: `gp(u, v)` with 60 observations all
#' at `v = 0` and grids at `v = +0.75` and `v = -0.75`. This returns
#' TRUE there and the omitted standard error is 0.888, so the escape is
#' real and large. It is blocked by a SECOND condition rather than by
#' this one: a `gp()` carries one length scale per dimension, so
#' observations on a line do not identify the scale across it, the fit
#' comes back singular and `frm_curve()` refuses at the covariance
#' check. Moving three observations off the line makes the fit healthy
#' and makes this predicate return `FALSE`. A fit that pinned the second
#' length scale with a prior or a bound would be non-singular and would
#' take the escape, so this is a second line of defence and not an
#' airtight one.
#'
#' Everything else fails closed: a component vector that does not line
#' up with the design, or no non-fixed column at all, returns `FALSE`
#' and the call refuses.
#'
#' @noRd
sp_same_latent <- function(fit, a, b) {
  nm <- frmtmb::frm_joint_cov(fit)$names[a$lb$coef_pos]
  if (length(nm) != ncol(a$C) || anyNA(nm)) return(FALSE)
  keep <- !(nm %in% c("beta", "betad"))
  if (!any(keep)) return(FALSE)
  identical(a$C[, keep, drop = FALSE], b$C[, keep, drop = FALSE]) &&
    identical(a$lb$extra_var, b$lb$extra_var)
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
                         allow_new_levels, tol, side,
                         lost = rep(FALSE, length(se))) {
  # frm_linpred() says the same thing about the same rows; the caller
  # says it once, under the name of the function the user called
  ref <- withCallingHandlers(
    frmtmb::frm_linpred(fit, newdata = nd, type = "link", dpar = dpar,
                        resp = resp, re_formula = re_formula,
                        allow_new_levels = allow_new_levels,
                        se.fit = TRUE),
    frmtmb_se_lost_prediction = function(w) invokeRestart("muffleWarning"))
  se_ref <- as.numeric(ref$se.fit)
  # rows without a standard error are NaN on both routes and carry
  # nothing to compare
  if (all(lost)) return(NA_real_)
  rel <- max(abs(se[!lost] / pmax(se_ref[!lost], .Machine$double.eps) - 1))
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

#' Refuse a route that would drop a new level's variance.
#'
#' `frm_lp_basis()` returns a new grouping level's variance one number
#' per ROW, as `extra_var`, and no covariance between rows. One unseen
#' level is ONE draw shared by every row, so the grid's covariance is
#' `C V C' + Z S Z'` and the seam hands over only the diagonal of the
#' second part. The pointwise band adds that diagonal and is right. The
#' simultaneous band draws from `C V C'` and divides by a standard error
#' that includes it, so its critical value fell BELOW `qnorm(0.975)`
#' (0.646 on `(1 | subject)`, coverage 0.352); the derivative and the
#' feature standard errors are built from `C V C'` alone and came out
#' at 0.40 to 0.50 and about a quarter of the right values
#' (dev/reviews/2026-09-29-splinecurve.md, Finding A). Before this
#' package took `allow_new_levels` every one of these calls refused, so
#' refusing keeps them where they were until the seam returns the
#' covariance.
#'
#' Only under `allow_new_levels = TRUE`. Under `FALSE` the one source
#' of `extra_var` is an exact `gp()` off the observed positions, which is
#' older than this argument and of a recorded size of 7e-07; it is filed
#' for the consolidating session rather than refused here.
#'
#' `what` names the route at run time so the three callers share one
#' template.
#'
#' @noRd
sp_new_level_stop <- function(parts, fn, what) {
  if (!isTRUE(parts$allow_new_levels) || !any(parts$extra_var != 0)) {
    return(invisible(NULL))
  }
  frm_stop(fn, ": with allow_new_levels = TRUE this grid carries a new ",
           "level's marginal variance, from ",
           sp_new_level_terms(parts$fit, parts$newdata), ". One unseen ",
           "level is one draw shared by every row, and frm_lp_basis() ",
           "returns its variance one number per row with no covariance ",
           "between rows, so ", what, " would omit it and come out too ",
           "narrow. The pointwise frm_curve(simultaneous = FALSE) band ",
           "carries it. Use re_formula = NA to drop the term and read the ",
           "population curve, or name a level the fit saw. When no row is ",
           "at an unseen level, as with an exact gp() read between its ",
           "observed positions, pass allow_new_levels = FALSE",
           call. = FALSE)
}

#' Which grouping factor of the grid holds a level the fit did not see.
#'
#' Read through the public `ngrps()` and `ranef()` rather than the
#' fitted frame. A grouping expression that is not a column, such as
#' `g1:g2`, cannot be matched to the grid by name and falls through to
#' the generic phrase; so does an exact `gp()`, the other source of
#' `extra_var`. Wrapped in `tryCatch()` because a refusal must not fail
#' while it names the thing it refuses.
#'
#' @noRd
sp_new_level_terms <- function(fit, nd) {
  found <- tryCatch({
    re <- frmtmb::ranef(fit)
    out <- character(0)
    for (g in names(frmtmb::ngrps(fit))) {
      if (!g %in% names(nd)) {
        # core fills an absent grouping column with NA under TRUE, so
        # every row is at a new level of it
        out <- c(out, paste0("`", g, "`, absent from the grid"))
        next
      }
      new <- setdiff(unique(as.character(nd[[g]])), rownames(re[[g]]))
      if (length(new)) {
        out <- c(out, paste0("`", g, "` at ",
                             paste(new, collapse = ", ")))
      }
    }
    out
  }, error = function(e) character(0))
  if (length(found)) paste(found, collapse = " and ") else {
    "a grouping term or an exact gp() this grid reaches off the fit"
  }
}

#' Everything the three exported functions share: the grid, the design,
#' the covariance, and the check that the covariance is the right one.
#'
#' `frm_lp_basis()` (frmtmb >= 0.52.0) returns `A`, `V` and the variance
#' that is not coefficient uncertainty. For a linear predictor it is the
#' same object `frm_linpred(se.fit = TRUE)` reduces to a diagonal, so the
#' two are compared on every call. For a NONLINEAR body `A` is a
#' Jacobian, `frm_linpred(se.fit = TRUE)` refuses outright, and there is
#' nothing to compare against: the check is skipped and `rel` is `NA`,
#' which `print()` reports rather than hides.
#'
#' With `contrast`, the reported functional is `(A1 - A2) c` and its
#' covariance is `(A1 - A2) V (A1 - A2)'`. Three things the one-grid
#' path takes for granted have to be established first, and they are the
#' whole difference between "the covariance is already assembled" and a
#' difference curve:
#'
#' 1. The two designs must sit at the SAME rows of `V`. `coef_pos` says
#'    where each one sits, so the subtraction is checked rather than
#'    assumed.
#' 2. `extra_var` is a per-row variance with no covariance between the
#'    grids. Where the two grids load the SAME latent draw the cross
#'    term is exactly zero and there is nothing to fetch, which
#'    `sp_same_latent()` decides; otherwise the difference is refused.
#' 3. The check compares each half against `frm_linpred(se.fit = TRUE)`.
#'    There is no second route to the difference's own standard error.
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
  a <- sp_one_basis(fit, newdata, dpar, resp, re_formula, allow_new_levels)
  nl <- sp_is_nl(fit, dpar, resp)
  out <- list(eta = a$lb$eta, C = a$C, V = a$lb$V, Sigma = a$Sigma,
              se = a$se, lost = a$lost, rel = NA_real_, n_predict = 0L,
              newdata = newdata, contrast = contrast, dpar = dpar,
              resp = resp, re_formula = re_formula,
              allow_new_levels = allow_new_levels, fit = fit, span = a$span,
              extra_var = a$lb$extra_var)
  if (is.null(contrast)) {
    # A nonlinear body is the case core refuses se.fit for, so there is
    # no second number to check against. Everything else is checked.
    if (!nl) {
      out$rel <- sp_cov_check(fit, newdata, a$se, dpar, resp, re_formula,
                              allow_new_levels, tol, "this grid", a$lost)
      out$n_predict <- 1L
    }
    return(out)
  }
  b <- sp_one_basis(fit, contrast, dpar, resp, re_formula, allow_new_levels)
  if (!identical(a$lb$coef_pos, b$lb$coef_pos)) {
    frm_stop("frm_curve(contrast = ): the two grids load on different ",
             "coefficients (", length(a$lb$coef_pos), " and ",
             length(b$lb$coef_pos), " of them), so subtracting their ",
             "designs would pair columns that belong to different ",
             "parameters. Both grids must reach the same linear predictor ",
             "under the same re_formula", call. = FALSE)
  }
  # Variance that is not coefficient uncertainty arrives per row with no
  # covariance between the grids, so a difference can only report it
  # when the two rows carry the SAME residual and it cancels exactly.
  # That is the ordinary case, a contrast across a factor at one gp()
  # position, and refusing it would refuse an answer that is right.
  has_extra <- any(a$lb$extra_var != 0) || any(b$lb$extra_var != 0)
  # An unseen level loads NO column of the design, so sp_same_latent()
  # cannot tell two unseen levels apart: both show all-zero columns and
  # the same marginal variance, and the predicate would call two
  # independent draws one draw and cancel a variance that belongs in the
  # answer. Refused before that predicate is ever asked.
  if (has_extra && isTRUE(allow_new_levels)) {
    frm_stop("frm_curve(contrast = , allow_new_levels = TRUE): a grid here ",
             "carries variance that is not coefficient uncertainty, a new ",
             "grouping level's marginal variance or an exact gp() kriging ",
             "residual. An unseen level loads no column of the design, so ",
             "nothing the seam returns says whether the two grids name the ",
             "SAME unseen level, whose draw cancels, or two different ones, ",
             "whose draws add. Drop the grouping term with re_formula = NA ",
             "to difference the population curves, or use ",
             "allow_new_levels = FALSE when no row is at an unseen level",
             call. = FALSE)
  }
  if (has_extra && !sp_same_latent(fit, a, b)) {
    frm_stop("frm_curve(contrast = ): this prediction carries variance that ",
             "is not coefficient uncertainty, which for a curve is an exact ",
             "gp() kriging residual. A difference can only report it when ",
             "both grids load the same one, and that is decided on the ",
             "whole latent design: EVERY column of it outside the fixed ",
             "effects has to match, not only the gp() block's. Here they do ",
             "not, so the gp() positions may well be identical and some ",
             "other random-effect or smooth term is what differs. ",
             "frm_lp_basis() returns that variance per row and no covariance ",
             "between the grids, so there is no cross term to fall back on. ",
             "Hold every latent term equal between the grids and contrast a ",
             "fixed effect, or read the two curves separately",
             call. = FALSE)
  }
  out$C <- a$C - b$C
  out$Sigma <- unname(out$C %*% out$V %*% t(out$C))
  out$se <- unname(sqrt(pmax(diag(out$Sigma), 0)))
  # conservative: a difference of two rows that each move along a lost
  # direction may itself be determined, but the seam marks rows, not
  # differences, so either half takes the row
  out$lost <- a$lost | b$lost
  out$se[out$lost] <- NaN
  out$eta <- a$lb$eta - b$lb$eta
  out$span <- sp_span_both(a$span, b$span)
  if (!nl) {
    r <- c(
      sp_cov_check(fit, newdata, a$se, dpar, resp, re_formula,
                   allow_new_levels, tol, "`newdata`", a$lost),
      sp_cov_check(fit, contrast, b$se, dpar, resp, re_formula,
                   allow_new_levels, tol, "`contrast`", b$lost))
    out$rel <- if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)
    out$n_predict <- 2L
  }
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
