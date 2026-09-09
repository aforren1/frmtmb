#' A fitted curve on a grid, with pointwise and simultaneous bands
#'
#' Evaluates a fitted linear predictor on a grid of covariate values and
#' returns it with two intervals: the usual pointwise interval, and a
#' SIMULTANEOUS band that covers the whole curve at once.
#'
#' A pointwise interval is the wrong tool for the question a curve
#' usually raises. "Is the velocity above zero at 300 ms" is pointwise;
#' "does this curve have the shape I claim" is a statement about every
#' point at once, and a 95 percent pointwise band covers the whole curve
#' far less than 95 percent of the time. The simultaneous band is the
#' max-deviation simulation of Ruppert, Wand and Carroll (2003, ch. 6):
#' draw the curve's own deviation process from its joint covariance,
#' standardize each draw by the pointwise standard error, take the
#' largest absolute value over the grid, and use the `level` quantile of
#' those maxima in place of `qnorm(0.975)`. It is the construction
#' `gratia::confint(type = "simultaneous")` uses on an mgcv fit, and this
#' package's simulation reproduces gratia's critical value inside its
#' Monte Carlo error.
#'
#' @section What the covariance is, and how it is checked:
#' A penalized smooth's wiggly part is a random-effect block in the
#' fitted objective even when the smooth is a population term, so the
#' covariance of a curve needs the joint covariance of the fixed AND
#' random coefficients. frmtmb exports no route to it: `vcov(full =
#' TRUE)` returns the outer parameter vector, which excludes `b` under
#' both of its branches, and `predict(se.fit = TRUE)` forms the grid
#' covariance internally and returns only its diagonal.
#'
#' So this function rebuilds it. The linear predictor is LINEAR in the
#' coefficients, so the difference between a prediction and the same
#' prediction with one coefficient raised by one is that coefficient's
#' design column, exactly. The joint covariance comes from the fit's own
#' joint precision matrix.
#'
#' Neither piece was handed over by an exported function, so neither is
#' trusted. Every call recomputes `sqrt(diag(Sigma))` and compares it
#' with `predict(se.fit = TRUE)`, and refuses when the two disagree by
#' more than `tol`. The measured agreement is in the `"check"` attribute
#' and is reported by `print()`. On the package's own test models it is
#' at the tenth significant figure or better.
#'
#' @section The route to the covariance:
#' Both halves come from frmtmb's own exported seam,
#' [frmtmb::frm_lp_basis()], which returns the design `A` over the
#' coefficient vector, the joint covariance `V` at exactly the rows `A`'s
#' columns sit at, and the variance that is NOT coefficient uncertainty
#' (a new grouping level's marginal variance, an exact `gp()`'s kriging
#' variance) as a separate element. `Sigma` is `A V A'`.
#'
#' Up to frmtmb 0.51.0 there was no such seam. This package rebuilt `A`
#' by unit perturbation, one `predict()` call per contributing
#' coefficient, and read `V` out of `fit$cache$Vjoint`, which was an
#' internal with no precedent to point at. Both are gone: the
#' reconstruction and the reach were replaced by one call, this package
#' now requires frmtmb (>= 0.52.0), and what the covariance check
#' verifies has changed from "the reconstruction reproduced core's
#' number" to "the seam is being read correctly".
#'
#' The check itself stays. Every call recomputes `sqrt(diag(Sigma))` and
#' compares it with `predict(se.fit = TRUE)`, and refuses when the two
#' disagree by more than `tol`. The measured agreement is in the
#' `"check"` attribute and is reported by `print()`.
#'
#' The one case with nothing to check against is a nonlinear (`nl =
#' TRUE`) body: `predict(se.fit = TRUE)` is refused there, so
#' `frm_lp_basis()` is the only route to the number and `cov_rel_error`
#' is `NA`. `print()` says so rather than reporting a check that never
#' ran.
#'
#' @section Cost:
#' What this call costs is dominated by ONE thing: the single
#' `predict(se.fit = TRUE)` check call, inside which core inverts the
#' fit's joint precision matrix over EVERY coefficient, including the
#' ones this curve does not touch. Measured at `re.form = NA` on a
#' 20-point grid, one process each:
#'
#' \itemize{
#'   \item `s(x, k = 10)`, 8 random coefficients: 0.29 s.
#'   \item `s(t, k = 8) + (1 + t | subject)`, 1000 subjects and 2006
#'     random coefficients: 0.98 s.
#'   \item the same over 4000 subjects, 8006 random coefficients:
#'     6.87 s.
#' }
#'
#' The design rebuild that used to sit beside those figures, and that
#' was a tenth of them at every size, is gone: `frm_lp_basis()` returns
#' the design core already had, so the `predict()` call count no longer
#' depends on the number of coefficients at all. The joint-precision
#' solve is now the whole cost, it is paid once because core memoizes
#' it, and it grows with the total number of coefficients in the fit
#' rather than with the grid.
#'
#' @section A difference curve:
#' `contrast` is a second grid of the same height. The curve returned is
#' then the DIFFERENCE of the two linear predictors, row by row, and its
#' covariance is `(A1 - A2) V (A1 - A2)'`, so both bands describe the
#' difference and the simultaneous one answers "is this difference
#' anywhere other than zero" over the whole grid at once. It is the
#' quantity `gratia::difference_smooths(group_means = TRUE)` reports for
#' a factor-by smooth, and on the same mgcv fit the two agree. Name that
#' argument when you compare: gratia's default, `group_means = FALSE`,
#' zeroes the intercept and the parametric group columns and reports the
#' smooth-only difference, which is a different quantity. On the fixture
#' `test-difference.R` uses it is 0.52 away.
#'
#' A difference is NOT the difference of two calls to this function. The
#' two curves share coefficients, so their covariance is what the
#' difference is made of, and adding two standard errors in quadrature
#' would ignore it.
#'
#' What the difference path cannot do, and refuses by name:
#'
#' \itemize{
#'   \item `transform = TRUE`. A difference of linear predictors is not
#'     a difference of responses under any link but the identity, so
#'     there is nothing to map it through.
#'   \item Two grids that load on different coefficients.
#'   \item Two grids that load DIFFERENT draws of a latent field whose
#'     variance is not coefficient uncertainty. See the next section.
#' }
#'
#' The covariance check also means less here, and `print()` says so.
#' `predict(se.fit = TRUE)` returns a marginal standard error per row
#' and never the covariance between the grids, so the check runs on each
#' half and `cov_rel_error` is the worse of the two: what it licenses is
#' that both designs were read correctly.
#'
#' @section An exact `gp()` under a difference:
#' An exact `gp()` evaluated off the observed positions carries
#' a kriging residual that is not coefficient uncertainty.
#' [frmtmb::frm_lp_basis()] returns its variance one number per ROW and
#' returns no covariance BETWEEN two grids, so `var(g1 - g2)` has no
#' public route.
#'
#' It needs none in the ordinary case. A contrast taken across a factor
#' at ONE `gp()` position leaves both grids loading the same residual,
#' so it cancels exactly and the difference is `(A1 - A2) V (A1 - A2)'`
#' with nothing left over. That case is computed rather than refused.
#'
#' Sameness is decided on the design and not on the numbers: the two
#' grids must agree bit for bit on EVERY column of `A` outside the fixed
#' effects. Equality of the variances would not be enough, because two
#' levels of one grouping block have identical marginal variances by
#' construction and are different draws.
#'
#' That test is stricter than the mathematics needs, and it is worth
#' knowing where the extra strictness bites. Only the block carrying the
#' kriging residual has to match for the residual to cancel, but the
#' test asks it of every latent column, so a contrast across `fac` on
#' `y ~ fac + s(x, by = fac) + gp(x)` is REFUSED even though the `gp()`
#' columns are identical: the by-factor smooth's own columns differ,
#' which is what a by-factor smooth is for. It fails closed, so the cost
#' is an answer you do not get rather than one you should not trust.
#'
#' The rest is refused because the SEAM cannot supply it, not because
#' the mathematics is missing. The conditional cross-covariance is
#' `k(x1, x2) - Xr1 K Xr2'`, and core forms every piece of it while
#' predicting, but reduces the result to one variance per row before the
#' seam returns. Hold every latent term equal between the grids and
#' contrast a fixed effect, or read the two curves separately.
#'
#' @param object A `frmtmb_fit` from [frmtmb::frm()].
#' @param newdata The grid, as a data frame. Every variable the linear
#'   predictor reads must be a column, held at the value the curve is
#'   wanted at.
#' @param contrast A second grid with the same number of rows, or `NULL`
#'   for an ordinary curve. With it the curve is `newdata` minus
#'   `contrast`, row by row.
#' @param dpar Distributional parameter to read the curve off. `NULL`,
#'   the default, is the location parameter `mu`.
#' @param resp Response name, for a multivariate fit.
#' @param re.form `NA` (the default) evaluates the population curve, the
#'   convention `mgcv` and `gratia` plot. `NULL` keeps every random
#'   effect, so the grid must carry the grouping columns and the curve is
#'   that group's own.
#' @param level Coverage of both intervals.
#' @param simultaneous Compute the simultaneous band. `FALSE` returns the
#'   pointwise interval alone and skips the simulation.
#' @param nsim Draws in the max-deviation simulation. The default 10000
#'   puts the Monte Carlo error of the critical value near 0.013; 200000
#'   puts it near 0.003.
#' @param transform Return the curve and both bands through the link
#'   inverse. The bands are transformed end to end rather than rebuilt,
#'   which keeps their coverage under any monotone link.
#' @param seed Seed for the simulation, for a reproducible band.
#' @param tol Largest relative disagreement with `predict(se.fit = TRUE)`
#'   the assembled covariance may show before the call refuses.
#'
#' @section Past a `ps()` knot span:
#' A [frmtmb::ps()] basis is a partition of unity only between its
#' frozen outer knots. Past them it is a partial sum that decays to
#' zero, so a curve drawn there bends smoothly to whatever the rest of
#' the body gives, which is exactly the shape a reader does not
#' question. `predict(newdata = )` says so, and so does
#' [frmtmb::frm_lp_basis()], the seam this function reads. It is
#' surfaced again here, ONCE per `ps()` term per call and carrying the
#' span, so that the sentence names the function you called: this one
#' reads the seam on the grid, but [frm_curve_deriv()] reads it on a
#' three-point difference stencil and [frm_curve_feature()] on a
#' five-point one, and core counts the rows it was handed.
#'
#' [frm_curve_feature()] REFUSES instead of warning. A band past the
#' span is visibly wrong on the page; a peak located past it is a
#' number with a standard error beside it and nothing to give it away.
#'
#' @return A data frame of class `frmtmb_curve`: the columns of
#'   `newdata`, then `.estimate`, `.se`, `.crit`, `.lower_ci`,
#'   `.upper_ci`, and when `simultaneous = TRUE` also `.crit_sim`,
#'   `.lower_sim` and `.upper_sim`. The grid covariance is the `"Sigma"`
#'   attribute, the fit is the `"fit"` attribute, and `"check"` carries
#'   the covariance agreement and the `predict()` call count.
#'
#' @references
#' Ruppert, D., Wand, M. P. and Carroll, R. J. (2003) *Semiparametric
#' Regression*. Cambridge University Press, ch. 6.
#'
#' @seealso [frm_curve_deriv()] for the derivative of the same curve,
#'   [frm_curve_feature()] for the location of a peak or a crossing.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200)))
#' dd$y <- 2 * sin(pi * dd$x) + rnorm(200, 0, 0.4)
#' fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
#'                    family = stats::gaussian(), data = dd)
#' cv <- frm_curve(fit, newdata = data.frame(x = seq(0, 1, length.out = 25)),
#'                 nsim = 2000)
#' head(cv[, c("x", ".estimate", ".se", ".lower_ci", ".lower_sim")])
#' @export
frm_curve <- function(object, newdata, contrast = NULL, dpar = NULL,
                      resp = NULL, re.form = NA, level = 0.95,
                      simultaneous = TRUE, nsim = 10000L,
                      transform = FALSE, seed = NULL, tol = 1e-6) {
  sp_check_level(level)
  sp_check_flag(simultaneous, "simultaneous")
  sp_check_flag(transform, "transform")
  sp_check_contrast(newdata, contrast)
  if (!is.null(contrast) && isTRUE(transform)) {
    stop("frm_curve(contrast = , transform = TRUE): a difference of two ",
         "linear predictors is not the difference of two responses ",
         "under any link but the identity, so there is no inverse to ",
         "return it through. Leave transform = FALSE", call. = FALSE)
  }
  sp_rp_gate(object)
  parts <- sp_curve_parts(object, newdata, dpar, resp, re.form, tol,
                          contrast)
  # re-raised under this function's own name rather than let out of the
  # seam as it stands: the user called frm_curve(), not frm_lp_basis(),
  # and the sibling functions hand the seam a stencil rather than the
  # grid, so core's row count needs a caller to explain it
  for (msg in parts$span) {
    warning(warningCondition(paste0(
      "frm_curve(): this grid leaves a ps() term's knot span, so the ",
      "band below is drawn around a decaying partial sum rather than ",
      "around the fitted curve. ", msg),
      class = "frmtmb_ps_span_warning"))
  }
  sp_assemble(parts, parts$eta, parts$se, parts$Sigma, level, simultaneous,
              nsim, transform, seed, newdata,
              what = if (is.null(contrast)) "value" else "difference")
}

#' Build the returned data frame from an estimate, a standard error and
#' a covariance. Shared by the curve and its derivatives, which differ
#' only in which linear functional of the coefficients they report.
#'
#' @noRd
sp_assemble <- function(parts, est, se, Sigma, level, simultaneous, nsim,
                        transform, seed, newdata, what) {
  crit <- stats::qnorm(1 - (1 - level) / 2)
  out <- newdata
  out[[".estimate"]] <- est
  out[[".se"]] <- se
  out[[".crit"]] <- crit
  out[[".lower_ci"]] <- est - crit * se
  out[[".upper_ci"]] <- est + crit * se
  sim <- NULL
  if (isTRUE(simultaneous)) {
    sp_check_count(nsim, "nsim")
    if (nrow(newdata) < 2L) {
      stop("frm_curve(simultaneous = TRUE) needs a grid of at least two ",
           "points: a band over one point is the pointwise interval",
           call. = FALSE)
    }
    sim <- sp_sim_crit(Sigma, se, nsim, level, seed)
    out[[".crit_sim"]] <- sim$crit
    out[[".lower_sim"]] <- est - sim$crit * se
    out[[".upper_sim"]] <- est + sim$crit * se
  }
  if (isTRUE(transform)) {
    linkinv <- sp_linkinv(parts)
    for (nm in c(".estimate", ".lower_ci", ".upper_ci", ".lower_sim",
                 ".upper_sim")) {
      if (!is.null(out[[nm]])) out[[nm]] <- linkinv(out[[nm]])
    }
    # a transformed band is no longer symmetric about the estimate, so a
    # standard error on the transformed scale would be read as if it
    # were; dropping it says outright that it is not there
    out[[".se"]] <- NULL
  }
  structure(out,
            class = c("frmtmb_curve", "data.frame"),
            Sigma = Sigma, fit = parts$fit, what = what,
            level = level,
            spec = list(newdata = newdata, contrast = parts$contrast,
                        dpar = parts$dpar, resp = parts$resp,
                        re.form = parts$re.form),
            check = list(cov_rel_error = parts$rel,
                         n_predict = parts$n_predict,
                         crit_mcse = if (is.null(sim)) NA_real_ else sim$mcse),
            row.names = seq_len(nrow(out)))
}

#' The link inverse of the linear predictor the curve came off.
#'
#' Read through `predict()` rather than off the family object: a dpar
#' carries its own link, and a nonlinear parameter carries none at all.
#' Two predictions on one row, one on each scale, identify the map.
#'
#' @noRd
sp_linkinv <- function(parts) {
  fit <- parts$fit
  nd <- parts$newdata[1L, , drop = FALSE]
  lk <- sp_predict_eta(fit, nd, parts$dpar, parts$resp, parts$re.form)
  rs <- try(as.numeric(stats::predict(fit, newdata = nd, type = "response",
                                      dpar = parts$dpar, resp = parts$resp,
                                      re.form = parts$re.form)),
            silent = TRUE)
  if (inherits(rs, "try-error") || length(rs) != 1L) {
    stop("frm_curve(transform = TRUE): this linear predictor has no ",
         "response scale to transform onto. predict(type = \"response\") ",
         "refuses it, so the curve stays on the link scale",
         call. = FALSE)
  }
  fam <- stats::family(fit)
  lnk <- fam[["links"]][[parts$dpar %||% "mu"]]
  if (!is.null(lnk) && is.function(lnk[["linkinv"]])) return(lnk[["linkinv"]])
  if (isTRUE(all.equal(lk, rs))) return(identity)
  stop("frm_curve(transform = TRUE): the link inverse of dpar '",
       parts$dpar %||% "mu", "' is not reachable from the family object, ",
       "so the band cannot be transformed. Leave transform = FALSE and ",
       "transform the columns yourself", call. = FALSE)
}

#' @export
print.frmtmb_curve <- function(x, ...) {
  ck <- attr(x, "check")
  cat("<frmtmb curve> ", attr(x, "what"), ", ", nrow(x),
      " grid points, level ", attr(x, "level"), "\n", sep = "")
  if (!is.null(x[[".crit_sim"]])) {
    cat("  critical value: pointwise ", format(x[[".crit"]][1L], digits = 5),
        ", simultaneous ", format(x[[".crit_sim"]][1L], digits = 5),
        " (mcse ", format(ck$crit_mcse, digits = 2), ")\n", sep = "")
  } else {
    cat("  critical value: pointwise ", format(x[[".crit"]][1L], digits = 5),
        "\n", sep = "")
  }
  # a SUBSET of a curve keeps its rows and loses its attributes, so
  # this reads a zero-length value rather than a missing one
  if (!length(ck$cov_rel_error) || is.na(ck$cov_rel_error)) {
    cat("  covariance NOT checked: predict(se.fit = TRUE) is refused",
        " for a nonlinear predictor, so there is no second route to",
        " compare against\n", sep = "")
  } else if (!is.null(attr(x, "spec")[["contrast"]])) {
    # the number is the worse of the two grids, and it says the two
    # DESIGNS were read correctly. predict(se.fit = TRUE) has no
    # covariance between the grids to offer, so the difference's own
    # standard error has no second route and print() must not imply one
    cat("  each grid checked against predict(se.fit = TRUE) to ",
        format(ck$cov_rel_error, digits = 3),
        " relative; the difference itself has no second route\n",
        sep = "")
  } else {
    cat("  covariance checked against predict(se.fit = TRUE) to ",
        format(ck$cov_rel_error, digits = 3), " relative\n",
        sep = "")
  }
  n <- min(6L, nrow(x))
  print(as.data.frame(x)[seq_len(n), , drop = FALSE])
  if (nrow(x) > n) cat("  ... ", nrow(x) - n, " more rows\n", sep = "")
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @noRd
sp_check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between 0 and 1",
         call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
sp_check_flag <- function(x, nm) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", nm, "` must be TRUE or FALSE", call. = FALSE)
  }
  invisible(NULL)
}

#' The second grid of a difference curve, checked before anything is
#' predicted.
#'
#' Equal height is the contract, not recycling: `A1 - A2` is row by row,
#' and a one-row contrast silently recycled against a fifty-row grid is
#' a different quantity from the one the argument name promises. A user
#' who wants a fixed reference profile repeats the row themselves, which
#' says so in their own code.
#'
#' @noRd
sp_check_contrast <- function(newdata, contrast) {
  if (is.null(contrast)) return(invisible(NULL))
  if (!is.data.frame(contrast) || !nrow(contrast)) {
    stop("`contrast` must be a data frame with at least one row: it is ",
         "the second grid the curve is differenced against",
         call. = FALSE)
  }
  if (!is.data.frame(newdata) || nrow(contrast) != nrow(newdata)) {
    stop("`contrast` must have the same number of rows as `newdata`, ",
         "because the difference is taken row by row. It has ",
         nrow(contrast), " against ",
         if (is.data.frame(newdata)) nrow(newdata) else "none",
         ". Repeat the row yourself to difference against one profile",
         call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
sp_check_count <- function(x, nm, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min ||
      x != round(x)) {
    stop("`", nm, "` must be a single whole number of at least ", min,
         call. = FALSE)
  }
  invisible(NULL)
}
