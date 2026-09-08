#' Coherence and phase from a fitted coupling model, with intervals
#'
#' Reads the `coh` and `phase` distributional parameters off a
#' [cross_wishart()] fit and returns them with confidence intervals.
#'
#' @param fit A [frmtmb::frm()] fit whose family is [cross_wishart()].
#' @param newdata Optional data frame of predictor values. Defaults to
#'   the data the model was fitted to.
#' @param re.form `NULL` keeps the random effects, so the answer is per
#'   group; `NA` drops them, so the answer is the population one. Passed
#'   through to [stats::predict()].
#' @param level Confidence level.
#' @param allow_new_levels Passed through to [stats::predict()].
#'
#' @return A data frame with one row per row of `newdata` and columns
#'   `.estimate`, `.se`, `.lower` and `.upper`. `frm_coherence()` adds
#'   `.eta`, the value on the logit scale the interval was built on.
#'
#' @section Why the interval is built on the link scale:
#' Coherence lives on `(0, 1)` and a Wald interval on that scale walks
#' off the end whenever the estimate is near a boundary, which is
#' exactly where few segments put it. The interval here is a Wald
#' interval on the logit scale pushed through `plogis()`, so it cannot
#' leave `(0, 1)` and it is asymmetric in the direction the sampling
#' distribution actually is.
#'
#' Measured coverage of that interval, 40 subjects with a coherence
#' random effect, 150 replicates: 0.912 at 4 segments per subject and
#' 0.893 at 16, when every distributional parameter carries a random
#' effect. With the two channel powers left as free per-subject
#' parameters instead, coverage falls to 0.622 and 0.820, because two
#' free numbers per subject that gain no information as subjects are
#' added starve the variance component. `dev/xspec-findings.md` in the
#' repository has the table. **Put `(1 | id)` on every dpar, not only on
#' `coh`.**
#'
#' @section Phase is an angle:
#' `frm_phase()` returns radians and does no wrapping: the linear
#' predictor is unbounded and a phase of `-3.1` and one of `3.18` are
#' the same angle reached from two sides. When a group's phase is near
#' the wrap point its interval is wide for a reason that is a
#' parameterization artifact rather than data, and the fit is worth
#' re-centering with an offset. Nothing here detects that for you.
#'
#' @examples
#' set.seed(5)
#' src <- rnorm(4096)
#' xs <- frm_cross_spectrum(src + rnorm(4096), 0.9 * src + rnorm(4096),
#'                          segments = 16)
#' xs <- xs[xs$freq < 0.1, ]
#' fit <- frmtmb::frm(
#'   frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
#'              pow2 ~ 1, coh ~ 1, phase ~ 1),
#'   family = cross_wishart(), data = xs)
#' frm_coherence(fit, newdata = xs[1, ])
#' frm_phase(fit, newdata = xs[1, ])
#' @export
frm_coherence <- function(fit, newdata = NULL, re.form = NULL,
                          level = 0.95, allow_new_levels = FALSE) {
  p <- cp_link_se(fit, "coh", newdata, re.form, level, allow_new_levels)
  data.frame(.estimate = stats::plogis(p$fit), .se = p$se.fit,
             .lower = stats::plogis(p$fit - p$z * p$se.fit),
             .upper = stats::plogis(p$fit + p$z * p$se.fit),
             .eta = p$fit)
}

#' @rdname frm_coherence
#' @export
frm_phase <- function(fit, newdata = NULL, re.form = NULL,
                      level = 0.95, allow_new_levels = FALSE) {
  p <- cp_link_se(fit, "phase", newdata, re.form, level, allow_new_levels)
  data.frame(.estimate = p$fit, .se = p$se.fit,
             .lower = p$fit - p$z * p$se.fit,
             .upper = p$fit + p$z * p$se.fit)
}

#' The one route both extractors take.
#'
#' Everything is read through predict(se.fit = TRUE), which is a
#' documented seam, rather than by reassembling a design from the fit's
#' internals. That is why re.form and allow_new_levels work here at all:
#' they are predict()'s, passed through.
#'
#' @noRd
cp_link_se <- function(fit, dpar, newdata, re.form, level,
                       allow_new_levels) {
  cp_require_family(fit, if (identical(dpar, "coh")) "frm_coherence()"
                    else "frm_phase()")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
        level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between 0 and 1.",
         call. = FALSE)
  }
  p <- stats::predict(fit, newdata = newdata, type = "link", dpar = dpar,
                      re.form = re.form, se.fit = TRUE,
                      allow_new_levels = allow_new_levels)
  if (!is.list(p) || is.null(p$se.fit)) {
    stop("predict() returned no standard error for `", dpar,
         "`, so no interval can be formed. This happens when the fit did ",
         "not produce a positive definite Hessian; check summary(fit).",
         call. = FALSE)
  }
  # A zero or non-finite standard error would otherwise be pushed through
  # plogis() and returned as an interval whose lower, estimate and upper
  # are the same number. A 95 percent interval of width zero is a wrong
  # answer rather than a wide one, so it is refused here at the point the
  # number is consumed.
  bad <- !is.finite(p$se.fit) | p$se.fit <= 0
  if (any(bad)) {
    stop("the standard error of `", dpar, "` is ",
         if (all(!is.finite(p$se.fit[bad]))) "not finite" else
           "zero or negative",
         " in ", sum(bad), " of ", length(bad),
         " rows, so no interval exists for them. This is what a coherence ",
         "estimated at 1 to machine precision looks like: the linear ",
         "predictor has run off to where the fit carries no curvature. ",
         "Refit with more segments per unit, or with a prior on `coh`.",
         call. = FALSE)
  }
  list(fit = as.numeric(p$fit), se.fit = as.numeric(p$se.fit),
       z = stats::qnorm(1 - (1 - level) / 2))
}

#' @noRd
cp_require_family <- function(fit, what) {
  fam <- tryCatch(stats::family(fit)$family, error = function(e) NULL)
  if (!identical(fam, "cross_wishart")) {
    stop(what, " reads the coherence and phase parameters of a ",
         "cross_wishart() fit; this fit's family is ",
         if (is.null(fam)) "not readable" else paste0("`", fam, "`"), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Draw whole cross-spectral matrices from a fitted coupling model
#'
#' [stats::simulate()] refuses on a [cross_wishart()] fit, because a
#' draw is a Hermitian matrix and the response slot carries one of its
#' four numbers. This returns all four, so a simulated data frame can be
#' refitted or compared against the observed one.
#'
#' @param fit A [cross_wishart()] fit.
#' @param nsim Number of replicate data frames.
#' @param seed Optional seed, set with [set.seed()] before drawing.
#' @param newdata Optional data frame of predictor values.
#' @param re.form Passed to [stats::predict()]; `NULL` keeps the random
#'   effects.
#'
#' @return A list of `nsim` data frames, each with the columns
#'   [frm_cross_spectrum()] produces: `w11`, `w22`, `w12r`, `w12i`, `n`.
#'
#' @section How a draw is made:
#' Each row's fitted spectral matrix is factored as `L L^H` and `n`
#' independent complex normal vectors are drawn through `L`; the
#' returned matrix is their outer-product sum. That is the definition of
#' the complex Wishart rather than a transformation of it, so a
#' simulated row and the density are one statement of the model. Checked
#' against the density's own first moment: over 20000 draws the mean of
#' the simulated matrix divided by `n` matches the fitted matrix to
#' 0.0011 at `n = 16`.
#'
#' @examples
#' set.seed(6)
#' src <- rnorm(2048)
#' xs <- frm_cross_spectrum(src + rnorm(2048), 0.9 * src + rnorm(2048),
#'                          segments = 16)
#' xs <- xs[xs$freq < 0.1, ]
#' fit <- frmtmb::frm(
#'   frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
#'              pow2 ~ 1, coh ~ 1, phase ~ 1),
#'   family = cross_wishart(), data = xs)
#' str(frm_cross_simulate(fit, nsim = 2, seed = 1)[[1]])
#' @export
frm_cross_simulate <- function(fit, nsim = 1L, seed = NULL,
                               newdata = NULL, re.form = NULL) {
  cp_require_family(fit, "frm_cross_simulate()")
  nsim <- cp_count(nsim, "nsim")
  if (!is.null(seed)) set.seed(seed)
  get <- function(dp) as.numeric(stats::predict(fit, newdata = newdata,
                                                type = "link", dpar = dp,
                                                re.form = re.form))
  s11 <- exp(get("mu")); s22 <- exp(get("pow2"))
  ch <- stats::plogis(get("coh")); ph <- get("phase")
  dat <- if (is.null(newdata)) fit$frame[["data"]] else newdata
  n <- as.numeric(if (!is.null(dat) && !is.null(dat$n)) dat$n
                  else stats::model.frame(fit)$n)
  if (length(n) != length(s11)) {
    stop("the degrees of freedom could not be recovered for the rows ",
         "being simulated; pass `newdata` carrying its own `n` column.",
         call. = FALSE)
  }
  ## L = [[a, 0], [c, b]] with a = sqrt(S11), |c|^2 = C S22, b^2 = (1-C) S22.
  ## Conjugated so that the drawn cross term has argument +phase, which is
  ## the convention frm_cross_spectrum() writes and frm_phase() reads.
  a <- sqrt(s11); cm <- sqrt(ch * s22); b <- sqrt((1 - ch) * s22)
  lapply(seq_len(nsim), function(k) {
    m <- length(s11)
    w11 <- numeric(m); w22 <- numeric(m)
    w12r <- numeric(m); w12i <- numeric(m)
    for (i in seq_len(m)) {
      ni <- n[i]
      z1 <- complex(real = stats::rnorm(ni, 0, sqrt(0.5)),
                    imaginary = stats::rnorm(ni, 0, sqrt(0.5)))
      z2 <- complex(real = stats::rnorm(ni, 0, sqrt(0.5)),
                    imaginary = stats::rnorm(ni, 0, sqrt(0.5)))
      d1 <- a[i] * z1
      d2 <- cm[i] * complex(modulus = 1, argument = -ph[i]) * z1 + b[i] * z2
      cr <- sum(d1 * Conj(d2))
      w11[i] <- sum(Mod(d1)^2); w22[i] <- sum(Mod(d2)^2)
      w12r[i] <- Re(cr); w12i[i] <- Im(cr)
    }
    data.frame(w11 = w11, w22 = w22, w12r = w12r, w12i = w12i, n = n)
  })
}
