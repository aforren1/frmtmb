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
#' @section The one internal this reaches into:
#' The covariance is read from `fit$cache$Vjoint`, and that is a read
#' into frmtmb's INTERNALS rather than a sanctioned seam. The function
#' that writes it, `get_joint_cov()`, is `@noRd`; neither `fit$cache`
#' nor the `list(V =, names =)` shape of the memo appears in
#' `?frmtmb::`frmtmb-extension-api``. Unlike `fit$obj` and
#' `fit$estimates`, which other extensions already read, this one has
#' no precedent to point at. It is stated here rather than buried in a
#' development note because a user is entitled to know which of a
#' package's dependencies are contractual and which are not.
#'
#' It is made anyway, because every alternative was worse. Computing
#' the covariance here instead means a second `sdreport()` per call, a
#' dense Schur complement over coefficients the curve never touches
#' (114 s and 2.1 GB at 8000 random coefficients, against 8.2 s and
#' 1.2 GB for the cache read), and a covariance that goes round
#' `autoscale_sdreport()` and is therefore WRONG on an autoscaled fit.
#' A public-but-wrong route was traded for a private-but-correct one.
#'
#' What happens if core changes it, in full:
#'
#' \itemize{
#'   \item **The name or the shape changes.** The read returns `NULL`,
#'     the sparse fallback runs, and the answer is the same one about
#'     16 times slower. No wrong number.
#'   \item **The meaning of `V` changes without the name changing.**
#'     The covariance check catches it: every call compares
#'     `sqrt(diag(Sigma))` against `predict(se.fit = TRUE)` and refuses
#'     above `tol`. No wrong number.
#'   \item **It is absent.** It is absent on every FRESH fit, because
#'     it is a memo rather than a slot. This is not a dependency on the
#'     cache being warm: the `predict(se.fit = TRUE)` check runs first
#'     and warms it, which is why that call is ordered ahead of the
#'     covariance and why the ordering is enforced by an argument
#'     rather than by a comment.
#' }
#'
#' So the reach cannot produce a wrong answer; it can only become slow,
#' or refuse. The standing ask is an exported accessor,
#' `dev/spline-seam-proposal.md` Part 1a, which would remove it.
#'
#' @section Cost:
#' What this call costs is dominated by ONE thing, and it is not the
#' `predict()` call count that the `"check"` attribute reports. It is
#' the single `predict(se.fit = TRUE)` call, inside which core inverts
#' the fit's joint precision matrix over EVERY coefficient, including
#' the ones this curve does not touch. Measured at `re.form = NA` on a
#' 20-point grid, one process each:
#'
#' \itemize{
#'   \item `s(x, k = 10)`, 8 random coefficients: design rebuild 0.01 s,
#'     `predict(se.fit = TRUE)` 0.29 s.
#'   \item `s(t, k = 8) + (1 + t | subject)`, 1000 subjects and 2006
#'     random coefficients: 0.07 s against 0.98 s.
#'   \item the same over 4000 subjects, 8006 random coefficients:
#'     0.28 s against 6.87 s.
#' }
#'
#' The design rebuild is a tenth of the cost at every size, and the
#' covariance itself is FREE: it is read from the object
#' `predict(se.fit = TRUE)` has already cached, not recomputed. Size a
#' job from the joint-precision solve, which grows with the total number
#' of coefficients in the fit, and not from the grid or the call count.
#'
#' The design is rebuilt with one `predict()` call per contributing
#' coefficient, plus one probe per block of `24` that contributes
#' nothing, so the count does not grow with the number of LEVELS of a
#' grouping factor. On the model `vignette("curve-inference")` fits,
#' `v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)` over 20
#' subjects, that is 32 calls against 110 random coefficients; a
#' factor-smooth model with NO population smooth needs far fewer,
#' because at `re.form = NA` the `fs` term contributes nothing and the
#' population curve is a constant. `tests/testthat/test-curve.R` pins
#' both. The call count is in the `"check"` attribute.
#'
#' @param object A `frmtmb_fit` from [frmtmb::frm()].
#' @param newdata The grid, as a data frame. Every variable the linear
#'   predictor reads must be a column, held at the value the curve is
#'   wanted at.
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
frm_curve <- function(object, newdata, dpar = NULL, resp = NULL,
                      re.form = NA, level = 0.95, simultaneous = TRUE,
                      nsim = 10000L, transform = FALSE, seed = NULL,
                      tol = 1e-6) {
  sp_check_level(level)
  sp_check_flag(simultaneous, "simultaneous")
  sp_check_flag(transform, "transform")
  sp_rp_gate(object)
  parts <- sp_curve_parts(object, newdata, dpar, resp, re.form, tol)
  sp_assemble(parts, parts$eta, parts$se, parts$Sigma, level, simultaneous,
              nsim, transform, seed, newdata, what = "value")
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
            spec = list(newdata = newdata, dpar = parts$dpar,
                        resp = parts$resp, re.form = parts$re.form),
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
  cat("  covariance checked against predict(se.fit = TRUE) to ",
      format(ck$cov_rel_error, digits = 3), " relative, in ", ck$n_predict,
      " predict() calls\n", sep = "")
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

#' @noRd
sp_check_count <- function(x, nm, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min ||
      x != round(x)) {
    stop("`", nm, "` must be a single whole number of at least ", min,
         call. = FALSE)
  }
  invisible(NULL)
}
