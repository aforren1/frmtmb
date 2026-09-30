# Parametric bootstrap over refits (the bootMer analog).

#' Parametric bootstrap
#'
#' Simulates `nsim` response vectors from the fitted model (by default
#' with new group effects each draw; see `re_formula`), refits the model
#' to each through [refit()] (warm-started, no re-parsing), and collects
#' `FUN` of every refit. Draws whose refit fails are kept as `NA` rows;
#' draws whose optimizer does not report convergence are kept but
#' flagged.
#'
#' There is no standard `bootstrap` generic to implement
#' (`boot::boot` and `lme4::bootMer` are plain functions), hence the
#' `frm_` prefix.
#'
#' @param fit A `frmtmb_fit` for a univariate model.
#' @param FUN Function of a `frmtmb_fit` returning a numeric vector.
#'   Default: `fixef(f, flatten = TRUE)`, the fixed effects named as
#'   `vcov()` names its rows, so a bootstrap standard error lines up
#'   with the Wald one by name.
#' @param nsim Number of bootstrap draws.
#' @param seed Optional seed.
#' @param re_formula Which group-level terms stay at their fitted
#'   values; the rest are redrawn in every replicate, from the fitted
#'   covariance. It is read exactly as [predict()] and
#'   [simulate.frmtmb_fit()] read it. The default `NA` (or `~0`, `~1`)
#'   redraws every group-level effect, in the manner of lme4's
#'   `bootMer(use.u = FALSE)`; `NULL` holds all of them and redraws the
#'   observation noise alone; a one-sided formula holds the terms it
#'   names and redraws the others.
#'
#'   A smooth, `gp()` or `hsgp()` term is a term of the formula, not a
#'   group effect, and it is never redrawn, under any `re_formula`: a
#'   smooth indexed by a grouping factor (`s(g, bs = "re")`,
#'   `s(x, g, bs = "fs")`, a `t2()` with an `re` margin) included. Every
#'   replicate holds it at its fitted coefficients and the refit
#'   estimates it again, so the spread is the uncertainty of the fitted
#'   curve.
#'
#' @section Ordinal thresholds:
#' Every refit reuses the assembled design, so the number of thresholds
#' is the fitted model's and `FUN` receives the same parameters in the
#' same order in every replicate. A replicate whose response never
#' reaches the top category keeps that threshold, which the replicate's
#' own data no longer identify: its estimate is large, and the
#' replicate's log-likelihood is the one a fit with `thres(K)` pinned on
#' that replicate gives, to optimizer tolerance. With `thres(gr = )` the
#' per-level counts are kept the same way.
#'
#' `simulate()` followed by `frm()`, written out by hand, is not a refit:
#' `frm()` counts the thresholds of whatever data it is given, so a
#' replicate that never reaches the top category fits one threshold
#' fewer, at the same maximized log-likelihood. Write `thres(K)` there to
#' pin it.
#'
#' @return A `frmtmb_boot` object: `t0` (FUN at the original fit), `t`
#'   (`nsim` x `length(t0)` matrix), and `converged`. `confint()` gives
#'   percentile intervals.
#' @examples
#' set.seed(3)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' bs <- frm_bootstrap(fit, nsim = 20, seed = 1)
#' bs
#' confint(bs)
#' @export
frm_bootstrap <- function(fit, FUN = function(f) fixef(f, flatten = TRUE),
                          nsim = 500, seed = NULL, re_formula = NA) {
  boot_refits(fit, FUN, nsim, seed, re_formula)
}

#' The body of `frm_bootstrap()`, with one setting the public function
#' does not offer: `re_plan`, a `sim_re_plan()` list saying which blocks
#' each simulation redraws. `conditional_effects()` holds the group-level
#' terms it reads at observed levels and redraws the others, a set it
#' reads off its grid rather than off a formula.
#'
#' @noRd
boot_refits <- function(fit, FUN, nsim, seed, re_formula, re_plan = NULL) {
  # nsim is passed straight to simulate() as a length, and a bad one
  # used to surface as "invalid 'length' argument" from inside the
  # bootstrap loop
  check_count(nsim, "nsim", min = 1L)
  if (!is.null(seed)) set.seed(seed)
  t0 <- FUN(fit)
  if (!is.numeric(t0)) {
    frm_stop("FUN must return a numeric vector", call. = FALSE)
  }
  # refit() replaces the response of the FITTED rows, so the na.exclude
  # padding simulate() adds has to come back off
  sims <- na_unpad(fit, sim_fit_draws(fit, nsim = nsim,
                                     re_formula = re_formula,
                                     re_plan = re_plan))
  # An ordinal draw arrives as an ordered factor carrying the response's
  # own levels, but the fit stores the 1..K codes and refit() takes
  # newresp as given: handed a factor, as.vector() turns it into text and
  # every refit dies inside the objective. Match what the fit holds.
  if (length(fit$frame[["y"]]) == 1L && !is.factor(fit$frame[["y"]][[1L]]) &&
      any(vapply(sims, is.factor, TRUE))) {
    sims[] <- lapply(sims, function(v) if (is.factor(v)) as.integer(v) else v)
  }
  # nsim refits run off this fit's control; a verbose original would
  # otherwise narrate every draw
  if (!is.null(fit$control)) fit$control$verbose <- FALSE
  tmat <- matrix(NA_real_, nsim, length(t0))
  if (!is.null(names(t0))) colnames(tmat) <- names(t0)
  ok <- logical(nsim)
  for (s in seq_len(nsim)) {
    ft <- tryCatch(suppressWarnings(refit(fit, sims[[s]])),
                   error = function(e) NULL)
    if (is.null(ft)) next
    ok[s] <- ft$opt$convergence == 0
    tmat[s, ] <- FUN(ft)
  }
  structure(list(t0 = t0, t = tmat, nsim = nsim, converged = ok),
            class = "frmtmb_boot")
}

#' @export
print.frmtmb_boot <- function(x, ...) {
  frm_check_dots(...)
  cat("Parametric bootstrap:", x$nsim, "refits,",
      sum(!x$converged), "failed or not converged\n\n")
  bmean <- colMeans(x$t, na.rm = TRUE)
  tab <- data.frame(
    estimate = x$t0,
    bias = bmean - x$t0,
    se = apply(x$t, 2, stats::sd, na.rm = TRUE),
    lwr = apply(x$t, 2, stats::quantile, 0.025, na.rm = TRUE),
    upr = apply(x$t, 2, stats::quantile, 0.975, na.rm = TRUE)
  )
  rownames(tab) <- names(x$t0) %||% paste0("t", seq_along(x$t0))
  print(signif(tab, 5))
  invisible(x)
}

#' @export
confint.frmtmb_boot <- function(object, parm = NULL, level = 0.95, ...) {
  frm_check_dots(...)
  a <- (1 - level) / 2
  ci <- cbind(
    lwr = apply(object$t, 2, stats::quantile, a, na.rm = TRUE),
    upr = apply(object$t, 2, stats::quantile, 1 - a, na.rm = TRUE),
    est = object$t0
  )
  rownames(ci) <- names(object$t0) %||% paste0("t", seq_along(object$t0))
  if (!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
