# Internal predictor standardization (the lme4 >= 1.1.37 autoscale
# idea). Instead of rescaling accessors post hoc, exploit cheap
# re-taping: fit a standardized copy of the frame, map that optimum
# back to the original parameterization exactly (ML is invariant under
# affine reparameterization of X columns), and use it to warm-start
# the ordinary unscaled fit. Everything downstream - sdreport, vcov,
# predict, profile - sees only the original-scale fit.

#' Which columns qualify: parametric X columns (1..n_param_cols) that
#' are numeric with more than two distinct values. That excludes the
#' intercept, factor dummies and other binary contrasts, smooth basis
#' columns (appended after n_param_cols), and the zero-placeholder
#' `mo()`/`mi()` columns. Centering needs an intercept column to absorb
#' the shift, so an intercept-free linpred is scaled without centering.
#' Returns NULL when nothing qualifies anywhere.
#'
#' @noRd
autoscale_plan <- function(frame) {
  plan <- list()
  for (key in names(frame[["linpreds"]])) {
    lp <- frame[["linpreds"]][[key]]
    if (!is.null(lp[["nl_body"]])) next     # no design matrix of its own
    if (!is.null(lp[["constant"]])) next    # mapped to a fixed value
    if (is.null(lp[["X"]]) || lp[["n_param_cols"]] == 0L) next
    icpt <- match("(Intercept)",
                  colnames(lp[["X"]])[seq_len(lp[["n_param_cols"]])])
    cols <- integer(0)
    center <- numeric(0)
    scale <- numeric(0)
    for (j in seq_len(lp[["n_param_cols"]])) {
      if (!is.na(icpt) && j == icpt) next
      xj <- lp[["X"]][, j]
      if (length(unique(xj)) <= 2L) next
      s <- stats::sd(xj)
      if (!is.finite(s) || s == 0) next
      cols <- c(cols, j)
      center <- c(center, if (is.na(icpt)) 0 else mean(xj))
      scale <- c(scale, s)
    }
    if (!length(cols)) next
    plan[[key]] <- list(cols = cols, center = center, scale = scale,
                        icpt = icpt, par = lp[["par"]], idx = lp[["idx"]],
                        z = autoscale_plan_z(frame, key, lp, cols, scale))
  }
  if (!length(plan)) NULL else plan
}

#' The random slopes on a planned column, scaled with it.
#'
#' Scaling the X column alone leaves a random slope on the same column
#' as badly scaled as before. Its log standard deviation starts at 0,
#' which on a column spread 1e-6 is a slope variance of 1e-12 in the
#' predictor's units: effectively zero, where the gradient of a log
#' standard deviation vanishes, and the fit stays there (lane
#' wt-predfix, punch round 1, `dev/predfix-p1-slope.R`). Scaling the Z
#' column by the same `s` is an exact reparameterization, `b * s` and
#' `log sd + log s`, and leaves the correlation parameters alone,
#' because it is a diagonal change of the block's coordinates.
#'
#' Only the blocks where that holds by construction: `us`, `diag`,
#' `us_t` and `diag_t`, whose first `dim` theta entries are per-component
#' log standard deviations (or log scales), with plain single-membership
#' components and an identity coefficient map. The Z columns are checked
#' against the X column before they are trusted.
#'
#' @noRd
autoscale_plan_z <- function(frame, key, lp, cols, scale) {
  Z <- lp[["Z"]]
  if (is.null(Z)) return(NULL)
  cn <- colnames(lp[["X"]])[cols]
  out <- list()
  for (bk in frame[["re_blocks"]]) {
    if (!bk[["covstruct"]] %in% c("us", "diag", "us_t", "diag_t")) next
    if (!identical(as.integer(bk[["b_idx"]]), as.integer(bk[["c_idx"]]))) {
      next
    }
    D <- bk[["dim"]]
    sdi <- covstruct_registry[[bk[["covstruct"]]]]$sd_idx(D)
    for (comp in bk[["components"]]) {
      if (!identical(comp$lp_key, key) || !is.null(comp$mm)) next
      for (k in seq_len(comp$dim)) {
        m <- match(comp$cnms[k], cn)
        if (is.na(m)) next
        pos <- comp$offset + k
        zc <- bk[["c_idx"]][(seq_len(bk[["n_levels"]]) - 1L) * D + pos]
        xj <- lp[["X"]][, cols[m]]
        zx <- as.numeric(Matrix::rowSums(Z[, zc, drop = FALSE]))
        if (max(abs(zx - xj)) > 1e-10 * max(abs(xj))) next
        out[[length(out) + 1L]] <- list(
          zcols = zc, theta = bk[["theta_idx"]][sdi[pos]], scale = scale[m])
      }
    }
  }
  if (length(out)) out
}

#' The column spread below which the default engages autoscale. It is
#' the boundary `diagnose()` and lme4's `check.scaleX` already call
#' badly scaled (three orders of magnitude from one). Measured on
#' 0.62.0 (`dev/predfix-scalescan.R`, `dev/predfix-scalecal.R`): below
#' it the plain fit can stall with the coefficient near zero, because
#' the coefficient's absolute gradient is under `grad_tol` before the
#' optimizer has moved, and it reports convergence 0. Up to 91
#' log-likelihood units were lost at a spread of 1e-6. At a spread of
#' 1e-3 the plain fit fell short by at most 4.6e-4 units (n = 20), and
#' on the LARGE side it was right in 450 of 450 fits up to 1e10, so a
#' large spread does not engage it.
#'
#' @noRd
autoscale_small_sd <- 1e-3

#' The spread below which a column carrying a random slope (a `pl$z`
#' entry) engages the default. A random slope stalls at a much larger
#' spread than a fixed effect: its log standard deviation starts at 0,
#' so the slope variance starts at `sd(x)^2` in the response's units,
#' near zero once the column is small. Measured with
#' `autoscale = FALSE` on `y ~ x + (1 + x | g)`, gaussian and poisson,
#' two slope sds, 736 fits (`dev/predfix-p1-slopecal*.R`): short by more
#' than 1e-3 units at 0.03 and below, by up to 14.4 units; at most
#' 2.2e-4 at 0.06 and 1.8e-6 from 0.1 up. `autoscale = TRUE` was within
#' 1.3e-6 everywhere.
#'
#' @noRd
autoscale_slope_sd <- 5e-2

#' The plan the fit runs with, or `NULL` for a plain fit.
#' `autoscale = TRUE` plans every qualifying column and `FALSE` none.
#' The default, `NULL`, plans them only when a column carrying a random
#' slope is spread below `autoscale_slope_sd`, or a column is spread
#' below `autoscale_small_sd` AND its coefficient is an outer parameter:
#' under `REML = TRUE` or `profile = TRUE` the `mu` coefficients are
#' integrated by the inner Newton solver, which measured nearly
#' scale-free: at spread 1e-6 the plain fit was within 9.6e-7 units of
#' autoscale on all 48 such fits (within 1e-9 on 45), against up to 202
#' units short under ML. Such a column alone leaves the fit bit for bit
#' where it was.
#'
#' @noRd
autoscale_decide <- function(frame, control, REML) {
  ask <- control[["autoscale"]]
  if (isFALSE(ask)) return(NULL)
  plan <- autoscale_plan(frame)
  if (is.null(plan) || isTRUE(ask)) return(plan)
  inner_beta <- isTRUE(REML) || isTRUE(control[["profile"]])
  small <- vapply(plan, function(pl) {
    x_small <- !(inner_beta && identical(pl$par, "beta")) &&
      any(pl$scale < autoscale_small_sd)
    # a slope's theta is outer under REML and profile too
    z_small <- any(vapply(pl$z, function(zz) zz$scale < autoscale_slope_sd,
                          NA))
    x_small || z_small
  }, NA)
  if (any(small)) plan
}

#' Standardized copy of the frame: only the planned X columns change.
#'
#' @noRd
autoscale_frame <- function(frame, plan) {
  for (key in names(plan)) {
    pl <- plan[[key]]
    X <- frame[["linpreds"]][[key]]$X
    Xs <- sweep(X[, pl$cols, drop = FALSE], 2, pl$center)
    X[, pl$cols] <- sweep(Xs, 2, pl$scale, "/")
    frame[["linpreds"]][[key]]$X <- X
    for (zz in pl$z) {
      Z <- frame[["linpreds"]][[key]]$Z
      Z[, zz$zcols] <- Z[, zz$zcols, drop = FALSE] / zz$scale
      frame[["linpreds"]][[key]]$Z <- Z
    }
  }
  frame
}

#' Map a parameter list between the parameterizations. Column j scaled
#' as `(x - c_j)/s_j` gives scaled coefficient `b_j * s_j`, and the
#' linpred intercept absorbs `sum_j b_j * c_j` (`c_j` is zero without an
#' intercept, so absorption vanishes there). A random slope scaled with
#' its column (`pl$z`) moves its effects by the same factor and its log
#' standard deviation by `log s_j`. Every other `theta`, `b` and extra is
#' untouched: those column changes leave them identical at the optimum.
#'
#' @noRd
autoscale_map <- function(pars, plan, to = c("original", "scaled")) {
  to <- match.arg(to)
  for (pl in plan) {
    for (zz in pl$z) {
      sg <- if (to == "scaled") 1 else -1
      if (!is.null(pars[["theta"]])) {
        pars[["theta"]][zz$theta] <- pars[["theta"]][zz$theta] +
          sg * log(zz$scale)
      }
      if (length(pars[["b"]])) {
        pars[["b"]][zz$zcols] <- pars[["b"]][zz$zcols] * zz$scale^sg
      }
    }
    p <- pars[[pl$par]]
    j <- pl$idx[pl$cols]
    if (to == "scaled") {
      shift <- sum(p[j] * pl$center)
      p[j] <- p[j] * pl$scale
      if (!is.na(pl$icpt)) {
        p[pl$idx[pl$icpt]] <- p[pl$idx[pl$icpt]] + shift
      }
    } else {
      p[j] <- p[j] / pl$scale
      if (!is.na(pl$icpt)) {
        p[pl$idx[pl$icpt]] <- p[pl$idx[pl$icpt]] - sum(p[j] * pl$center)
      }
    }
    pars[[pl$par]] <- p
  }
  pars
}

#' Stage one of `frmtmb_control(autoscale = TRUE)`: fit the standardized
#' frame and return its optimum mapped back to the original
#' parameterization, as the warm-start template for the real fit. User
#' starts are translated into the scaled parameterization so they mean
#' the same model.
#'
#' @noRd
autoscale_prefit <- function(spec, frame, bform, cl, REML, start,
                             control, lower, upper, prior, quadrature,
                             plan) {
  ctl <- control
  ctl$autoscale <- FALSE
  # the caller times and reports this whole stage as one "autoscale
  # pre-fit" line; a nested copy of every stage line would only
  # interleave with the fit that is actually reported
  ctl$verbose <- FALSE
  if (!is.null(start)) {
    # through make_start(), so a partially NAMED start leaves the
    # entries it does not address at the defaults the pre-fit would
    # otherwise have used, rather than at zero
    full <- make_start(frame, start)
    full <- autoscale_map(full, plan, "scaled")
    # theta and b too: a random slope scaled with its column moves them
    for (nm in intersect(names(start), c("beta", "betad", "theta", "b"))) {
      start[[nm]] <- unname(full[[nm]])
    }
  }
  # The pre-fit only supplies a starting point, and the fit that follows
  # re-derives every check on the model the user wrote, so its warnings
  # would be duplicates at best. At worst they are false: a bound on the
  # Intercept binds on the CENTERED intercept here and not in the model,
  # and the pre-fit then warned "Large maximum absolute gradient" on a
  # fit that converged cleanly. What the outcome means for the reported
  # fit is decided by autoscale_prefit_verdict(), so an error is
  # returned rather than raised.
  sfit <- tryCatch(withCallingHandlers(
    fit_assembled(spec, autoscale_frame(frame, plan), bform, cl,
                  REML = REML, start = start, control = ctl,
                  se = FALSE, lower = lower, upper = upper,
                  prior = prior, quadrature = quadrature),
    warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) e)
  if (inherits(sfit, "error")) return(list(error = sfit))
  list(template = autoscale_map(sfit$estimates, plan, "original"),
       converged = autoscale_prefit_converged(sfit),
       message = sfit$opt$message)
}

#' Whether a pre-fit reached a point a reported fit may start from.
#'
#' A pre-fit that did not converge IN ITS COEFFICIENTS poisons the fit
#' it seeds. On a completely separated bernoulli fit the pre-fit ran its
#' coefficients to 1e4 to 4.5e5 with code 1 and a gradient of 1e-250;
#' the reported fit then started there, stopped at once with code 0 and
#' said nothing, where autoscale = FALSE warned (the wt-predfix
#' reviewer, punch round 2).
#'
#' A pre-fit whose only trouble is a variance component on its boundary
#' is different: a random slope with no slope variance stops with code 1
#' and a correlation parameter running off, and it is still the right
#' start. Rejecting it on the code alone cost 25.1 log-likelihood units
#' at a spread of 1e-6. So the test is on the coefficient block of the
#' outer Hessian (`beta`, `betad`): it must be finite and positive
#' definite with a reciprocal condition number above 1e-4. Measured on
#' standardized pre-fits (`dev/predfix-p2-prefit.R`): that ratio is at
#' most 3.4e-7, or an eigenvalue is not positive, on all 19 separated
#' fits that returned, and at least 0.041 on all 21 healthy and
#' boundary-variance fits. The block is the whole Hessian when the
#' pre-fit converged with code 0, which is the common case.
#'
#' @noRd
autoscale_prefit_converged <- function(sfit) {
  op <- sfit$opt
  if (!is.finite(op$objective) || !all(is.finite(op$par))) return(FALSE)
  if (!length(op$par)) return(identical(as.integer(op$convergence), 0L))
  H <- tryCatch(stats::optimHess(op$par, sfit$obj$fn, sfit$obj$gr),
                error = function(e) NULL)
  if (is.null(H) || !all(is.finite(H))) return(FALSE)
  pd <- function(M) {
    ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    all(ev > 0) && min(ev) / max(ev) > 1e-4
  }
  if (identical(as.integer(op$convergence), 0L) && pd(H)) return(TRUE)
  cb <- names(op$par) %in% c("beta", "betad")
  any(cb) && pd(H[cb, cb, drop = FALSE])
}

#' A finished fit that is a verified optimum: code 0, and an outer
#' Hessian, taken in the fit's natural units (`par_units`, as
#' autoscale_sdreport() takes it), that is finite and positive definite
#' with a reciprocal condition number above 1e-8. A separated or
#' otherwise flat fit fails the Hessian test.
#'
#' @noRd
autoscale_fit_sound <- function(fit) {
  op <- fit$opt
  if (!identical(as.integer(op$convergence), 0L)) return(FALSE)
  if (!length(op$par)) return(TRUE)
  u <- fit$par_units %||% rep(1, length(op$par))
  H <- tryCatch(stats::optimHess(op$par / u,
                                 function(q) fit$obj$fn(q * u),
                                 function(q) fit$obj$gr(q * u) * u),
                error = function(e) NULL)
  if (is.null(H) || !all(is.finite(H))) return(FALSE)
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  all(ev > 0) && min(ev) / max(ev) > 1e-8
}

#' Whether the reported fit may start from the pre-fit: `TRUE`, `FALSE`
#' (fit as `autoscale = FALSE` does), or `"choose"`. The default
#' (`autoscale = NULL`) was not asked for, so a pre-fit that errors is
#' dropped and the fit runs exactly as `autoscale = FALSE` runs it, and
#' one that did not converge is `"choose"`, settled by
#' autoscale_choose(). `autoscale = TRUE` was asked for, so its error
#' stands, naming `autoscale = FALSE` as a remedy, and a pre-fit that
#' did not converge still seeds the fit, with a warning.
#'
#' @noRd
autoscale_prefit_verdict <- function(pre, control, vb = 0L) {
  asked <- isTRUE(control[["autoscale"]])
  if (!is.null(pre$error)) {
    msg <- conditionMessage(pre$error)
    if (asked) {
      frm_stop(msg, ". This happened in the autoscale pre-fit, the ",
               "standardized copy of the model; ",
               "frmtmb_control(autoscale = FALSE) fits the model as ",
               "written", call. = FALSE,
               class = c("frmtmb_autoscale_prefit_error",
                         "frmtmb_fit_error"))
    }
    if (vb) vb_say("autoscale pre-fit failed (", msg, "); fitting without it")
    return(FALSE)
  }
  if (isTRUE(pre$converged)) return(TRUE)
  if (asked) {
    frm_warning("The autoscale pre-fit did not converge (", pre$message,
                "), and the fit starts from where it stopped. ",
                "frmtmb_control(autoscale = FALSE) fits the model ",
                "without it", call. = FALSE)
    return(TRUE)
  }
  if (vb) {
    vb_say("autoscale pre-fit did not converge (", pre$message,
           "); fitting both from it and without it")
  }
  "choose"
}

#' Under the default, a pre-fit that did not converge cleanly may still
#' be the right start or may be poison, and the pre-fit alone cannot
#' tell which. A completely separated fit and a skew-normal `alpha` on a
#' 1e-6 column both leave a flat coefficient direction; seeding the
#' first hid the separation, and refusing the second lost up to 14.46
#' log-likelihood units that seeding recovers
#' (`dev/predfix-p2-skewcmp.R`). So the fit is run both ways and the
#' seeded one is reported only if its objective is lower AND either it
#' warns at least as often as the plain one or it is a verified optimum
#' (autoscale_fit_sound()). The second clause is for the skew-normal
#' case: there the plain fit's only warning is about its own stall,
#' 9.6 to 14.5 units short, and the seeded fit is a clean optimum. A
#' separated fit is never sound, so it keeps the plain fit and its
#' warning. The reported fit's own warnings are replayed, and the plain
#' fit's messages (where a start came from) are shown once. Only a
#' questionable pre-fit pays for the two extra fits.
#'
#' @noRd
autoscale_choose <- function(run, control, template, vb = 0L) {
  capture <- function(expr) {
    w <- list()
    m <- list()
    v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
      warning = function(x) {
        w[[length(w) + 1L]] <<- x
        invokeRestart("muffleWarning")
      },
      message = function(x) {
        m[[length(m) + 1L]] <<- x
        invokeRestart("muffleMessage")
      })
    list(value = v, warnings = w, messages = m,
         error = inherits(v, "error"))
  }
  plain <- capture(run(utils::modifyList(control, list(autoscale = FALSE)),
                       NULL))
  seeded <- capture(run(control, template))
  obj <- function(r) if (r$error) Inf else r$value$opt$objective
  tol <- 1e-8 * max(1, abs(obj(plain)[is.finite(obj(plain))]))
  pick_seeded <- !seeded$error &&
    (plain$error || (obj(seeded) < obj(plain) - tol &&
                       (length(seeded$warnings) >= length(plain$warnings) ||
                          (!length(seeded$warnings) &&
                             autoscale_fit_sound(seeded$value)))))
  # replayed as the condition objects they were, class and call kept
  for (x in plain$messages) frm_message(x)
  chosen <- if (pick_seeded) seeded else plain
  if (vb) {
    vb_say("kept the fit ", if (pick_seeded) "from the pre-fit" else
             "without the pre-fit")
  }
  for (x in chosen$warnings) frm_warning(x)
  if (chosen$error) frm_stop(chosen$value)
  chosen$value
}

#' sdreport on the unscaled optimum needs two repairs when non-unit
#' autoscale units are in play: optimHess steps are absolute (far larger
#' than a 1e-6-magnitude coefficient, so the default outer Hessian is
#' NaN or garbage), and no solver can invert a Hessian whose scaled-in
#' condition number exceeds double precision. Both problems vanish in
#' the coordinates `q = par / unit`, the scaled problem, so compute the
#' Hessian there, hand its exact back-map to sdreport (keeping the
#' joint precision and `diag.cov.random` finite), and overwrite
#' `cov.fixed` with the exactly transformed well-conditioned inverse.
#'
#' @noRd
autoscale_sdreport <- function(fit, jp = needs_jp(fit)) {
  u <- fit$par_units
  if (is.null(u) || all(u == 1)) {
    return(sdr_name_joint(RTMB::sdreport(fit$obj, getJointPrecision = jp),
                          fit$obj))
  }
  obj <- fit$obj
  q0 <- fit$opt$par / u
  Hq <- stats::optimHess(q0, function(q) obj$fn(q * u),
                         function(q) obj$gr(q * u) * u)
  sdr <- RTMB::sdreport(obj, par.fixed = fit$opt$par,
                        hessian.fixed = Hq / outer(u, u),
                        getJointPrecision = jp)
  Vq <- try(solve(Hq), silent = TRUE)
  if (!inherits(Vq, "try-error")) {
    V <- Vq * outer(u, u)
    dimnames(V) <- dimnames(sdr$cov.fixed)
    sdr$cov.fixed <- V
    sdr$pdHess <- !inherits(try(chol(Hq), silent = TRUE), "try-error")
  }
  sdr_name_joint(sdr, obj)
}

#' Name the joint precision when sdreport left it unnamed.
#'
#' TMB names the joint precision only when there is an outer parameter;
#' with none it returns the inner Hessian bare. That is a REML or
#' `profile = TRUE` fit whose family has no free dispersion (poisson,
#' bernoulli) and no random effect: every parameter is a `beta` and
#' inner. Every reader here finds its blocks by row name, so the bare
#' matrix killed `vcov()` with a `dimnames` error. With no outer
#' parameter the rows are exactly the inner parameters in order.
#'
#' @noRd
sdr_name_joint <- function(sdr, obj) {
  Q <- sdr$jointPrecision
  if (is.null(Q) || !is.null(rownames(Q))) return(sdr)
  nm <- names(obj$env$par)[obj$env$random]
  if (length(nm) == nrow(Q) && !length(sdr$par.fixed)) {
    dimnames(sdr$jointPrecision) <- list(nm, nm)
  }
  sdr
}

#' Natural magnitude per outer parameter of the warm-started unscaled
#' fit, aligned with `names(obj$par)`: `1/s_j` for the coefficient of a
#' column with sample SD `s_j`, 1 elsewhere. Multiplying the gradient by
#' these units reproduces the well-conditioned scaled-fit gradient.
#'
#' @noRd
autoscale_units <- function(frame, plan, par_names) {
  unit <- lapply(frame[["par_template"]], function(v) rep(1, length(v)))
  for (pl in plan) {
    unit[[pl$par]][pl$idx[pl$cols]] <- 1 / pl$scale
  }
  out <- numeric(0)
  for (cp in unique(par_names)) {
    u <- unit[[cp]]
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      u <- u[-frame[["betad_fixed_idx"]]]
    }
    out <- c(out, u)
  }
  stopifnot(length(out) == length(par_names))
  out
}
