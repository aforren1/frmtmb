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
    if (!is.null(lp[["equate"]])) next      # its target is planned
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
#' the coefficient's absolute gradient is under `grad_tol`, the
#' convergence check's trip-wire, before the optimizer has moved, and it
#' reports convergence 0. That is the one gap the trip-wire leaves, and
#' this default is what closes it (`grad_verdict()`). Up to 91
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
  # would be duplicates. They used to be worse: a bound on the Intercept
  # binds on the CENTERED intercept here and not in the model, and the
  # pre-fit then warned "Large maximum absolute gradient" on a fit that
  # converged cleanly. grad_verdict() no longer counts a gradient
  # component a bound holds in place, so that case is gone at its source
  # and the muffling stays only for the duplication. What the outcome
  # means for the reported
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
  # in the pre-fit's own optimizer units, as autoscale_fit_sound() reads
  # a fit: a smooth's null-space column has a small spread in brms's
  # basis (smooth_fx_units()), and in raw units its coefficient alone
  # put the block's condition ratio under the cutoff
  u <- sfit$par_units %||% rep(1, length(op$par))
  H <- tryCatch(stats::optimHess(op$par / u,
                                 function(q) sfit$obj$fn(q * u),
                                 function(q) sfit$obj$gr(q * u) * u),
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
  # the fit-time check (se_check()) built this Hessian already; it is
  # reused only at the point it was taken. Otherwise it is built here,
  # at the cost sdreport() would pay for it, because sdr_rescue() needs
  # its noise to tell a row of noise from curvature
  hc <- if (is.environment(fit$cache)) fit$cache$hessian_fixed
  fresh <- function() {
    h <- tryCatch(fit_outer_hessian(fit), error = function(e) NULL)
    if (!is.null(h) && is.environment(fit$cache)) {
      fit$cache$hessian_fixed <- h
    }
    h
  }
  if (is.null(u) || all(u == 1)) {
    at <- sdr_outer_point(fit$obj)
    if (is.null(hc) || !is.null(hc$Hq) || !identical(hc$at, at)) {
      hc <- if (length(at)) fresh()
    }
    H <- if (!is.null(hc) && is.null(hc$Hq) && identical(hc$at, at)) hc$H
    sdr <- RTMB::sdreport(fit$obj, hessian.fixed = H,
                          getJointPrecision = jp)
    return(sdr_rescue(fit, sdr_name_joint(sdr, fit$obj), H))
  }
  obj <- fit$obj
  q0 <- fit$opt$par / u
  if (is.null(hc$Hq) || !identical(hc$at, fit$opt$par)) hc <- fresh()
  Hq <- if (!is.null(hc$Hq) && identical(hc$at, fit$opt$par)) hc$Hq else {
    stats::optimHess(q0, function(q) obj$fn(q * u),
                     function(q) obj$gr(q * u) * u)
  }
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
  sdr_rescue(fit, sdr_name_joint(sdr, obj), Hq / outer(u, u))
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

#' Optimizer units for the null-space coefficients of the smooths:
#' `1/s_j` for a null-space column with sample SD `s_j`, 1 elsewhere,
#' aligned with `names(obj$par)`; `NULL` when the model has none.
#'
#' brms's basis (`diagonal.penalty = TRUE`, see the frame) gives these
#' columns a small spread, 0.13 to 0.16 on gamSim data, so their
#' coefficients are large. nlminb is not scale invariant, and in that
#' geometry it ran a smoothing SD that belongs at zero far down its log
#' scale, -90 where the old basis stopped at -8.5, until the outer
#' Hessian was singular and every fixed-effect standard error NaN: 17
#' of 140 gamSim fits against 3 (dev/fixes-sx-conv3.R, -conv6.R). The
#' reported coefficient stays brms's; only the steps the optimizer
#' takes, and the gradient its convergence check reads, are in unit-SD
#' coordinates, as autoscale's are.
#'
#' @noRd
smooth_fx_units <- function(frame, par_names) {
  unit <- lapply(frame[["par_template"]], function(v) rep(1, length(v)))
  hit <- FALSE
  for (lp in frame[["linpreds"]]) {
    X <- lp[["X"]]
    idx <- lp[["idx"]]
    if (is.null(X) || length(idx) != ncol(X)) next
    for (si in lp[["smooths"]] %||% list()) {
      for (j in si[["xf_idx"]] %||% integer(0)) {
        s <- stats::sd(as.numeric(X[, j]))
        if (!is.finite(s) || s <= 0) next
        unit[[lp[["par"]]]][idx[j]] <- 1 / s
        hit <- TRUE
      }
    }
  }
  if (!hit) return(NULL)
  out <- numeric(0)
  for (cp in unique(par_names)) {
    u <- unit[[cp]]
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      u <- u[-frame[["betad_fixed_idx"]]]
    }
    out <- c(out, u)
  }
  if (length(out) != length(par_names)) return(NULL)
  out
}

#' The unit of each linear predictor that is on the response's scale,
#' as a named vector over `frame$linpreds` keys; empty when none is.
#'
#' The `mu` of a location-scale family (`predictor_unit_families`) with
#' the identity link is in the units of the response, and its `sigma` is
#' the spread of what that `mu` leaves unexplained, in the same units.
#' The unit is that spread, measured before the fit as the sd of the
#' response (`log y` for the lognormal families) with `mu`'s fixed
#' design taken out (predictor_resid_sd()). It scales with the response
#' by construction. Every other predictor (a log or logit link, a
#' nonlinear parameter, a family with no `sigma`) is in units the
#' response's scale does not reach, and keeps 1.
#'
#' A mixture's components and an extension's family are left out. Their
#' likelihoods have several local maxima, and which one a fit reaches
#' follows its start and its path: with every `mu` with a `sigma` given
#' a unit, 68 fits of the eight suites ended at a lower maximum and 38
#' at a higher one, among them a `mixture(gaussian, gaussian)` that ran
#' a component's sigma to 0 (111 log-likelihood units lower) and
#' frmtmb.latent's hidden Markov models, which reached the optimum
#' `hmm_starts()` is tested to escape from (dev/optrap-findings.md).
#'
#' Without it a fit was not invariant to the units of the response. On
#' `y ~ x + (1 | g)` (20 groups of 5) with `y` in units 1e3 times
#' larger, 20 of 40 fits stopped short of lme4's maximum, by up to 0.67
#' log-likelihood units at optimizer code 0, and 4 of 40 (up to 0.82)
#' with `y` 1e3 times smaller (dev/optrap-findings.md). The group sd
#' started at 1 in the response's units, 1e-3 or 1e3 times the residual
#' sd, and nlminb stepped the mu coefficients in absolute units. Either
#' cure alone left 22 or 24 of the 40 short at 1e3 (dev/optrap-scale.R
#' with `OPTRAP_ABLATE`); both together left none: the sds start at the
#' unit (`unit_theta_start()`) and the mu coefficients are stepped in it
#' (`predictor_par_units()`).
#'
#' The unit is a power of two, so a residual spread between 0.71 and
#' 1.41 gives a unit of exactly 1 and the fit is the one 0.69.0 made,
#' and two data sets whose residuals spread alike get the same unit. It
#' does not make the fit of `y * 2^k` bitwise the scaled fit of `y`:
#' the log-scale parameters shift by `k log 2` and the objective by
#' `n k log 2`, neither exactly (13 of 30 fits followed the same path
#' to 1 ulp, 17 another; dev/optrap-rev-pow2.R).
#'
#' @noRd
predictor_units <- function(frame, fits = NULL) {
  out <- numeric(0)
  for (key in names(frame[["linpreds"]])) {
    u <- predictor_unit_of(frame, frame[["linpreds"]][[key]], fits[[key]])
    if (is.finite(u) && u > 0 && u != 1) out[[key]] <- u
  }
  out
}

#' The location-scale families whose `sigma` starts on the scale of
#' their identity-link `mu` (predictor_units()).
#'
#' huber() stays in although its fits stop short as often at unit 1 as
#' with a unit: its `k` is in units of sigma, so the unit does not
#' change it, and its Laplace objective jumps where a residual crosses
#' `k` (the inner Hessian is a step function), so a fit can stop on
#' either side of a jump. Without the unit, at y x 1e3 and 1e6 the fits
#' are 17 and 20 of 20 short, by up to 12.4; with it 11 and 10, by up to
#' 0.36 (dev/optrap-findings.md, punch round 1).
#'
#' @noRd
predictor_unit_families <- c("gaussian", "student", "lognormal",
                             "shifted_lognormal", "hurdle_lognormal",
                             "skew_normal", "exgaussian", "asym_laplace",
                             "zero_inflated_asym_laplace", "huber")

#' One linear predictor's unit for predictor_units(), or 1.
#'
#' `fit` is the least-squares fit of predictor_unit_input()'s response
#' on this predictor (predictor_ls()), when make_start() has already
#' made it for the starting values (start_ls_plan()), so that the
#' design is solved once.
#'
#' @noRd
predictor_unit_of <- function(frame, lp, fit = NULL) {
  inp <- predictor_unit_input(frame, lp)
  if (is.null(inp)) return(1)
  if (is.null(fit)) {
    fit <- tryCatch(predictor_ls(lp, inp$z, inp$ok, inp$w),
                    error = function(e) NULL)
  }
  raw <- tryCatch(predictor_resid_sd(lp, inp$z, inp$ok, fit = fit),
                  error = function(e) NA)
  if (!is.finite(raw) || raw <= 0) return(1)
  u <- 2^round(log2(raw))
  if (is.finite(u) && u > 0) u else 1
}

#' The response a predictor's unit is measured on, for
#' predictor_unit_of(): `z` (the response on the scale of `mu`), the
#' rows `ok` that set it and their weights `w` (`NULL` without
#' weights); `NULL` when the predictor gets no unit.
#'
#' @noRd
predictor_unit_input <- function(frame, lp) {
  if (!is.null(lp[["constant"]]) || !is.null(lp[["equate"]]) ||
        !is.null(lp[["nl_body"]])) return(NULL)
  if (!identical(lp[["link"]]$name %||% "", "identity")) return(NULL)
  # Not a mixture component's mu<k>: a unit there fixed the review's
  # mixture scale case (0 short at every scale) but broke four files of
  # the mixture suite. A random intercept started at the residual
  # spread, which holds the components' separation, put test-v18's
  # group-level mixture in a multimodal inner problem (false
  # convergence, 111 lower); the rest were the precision of the stop
  # (dev/optrap-findings.md, punch round 1)
  if (!identical(lp[["dpar"]], "mu")) return(NULL)
  fam <- frame[["spec"]]$responses[[lp[["resp"]]]]$family
  fname <- fam[["family"]]
  if (!is.character(fname) || length(fname) != 1L ||
        !fname %in% predictor_unit_families) return(NULL)
  # The unit is the scale of an estimated sigma. A sigma held at a
  # constant (se(), or `sigma = 2` in the formula) gives the residual
  # scale itself, and the student fuzz specs with se() that stopped
  # without converging stopped elsewhere with a unit, 64.7 higher on one
  # (test-fuzz.R, dev/optrap-findings.md)
  sig <- Filter(function(l) {
    identical(l[["resp"]], lp[["resp"]]) && identical(l[["dpar"]], "sigma")
  }, frame[["linpreds"]])
  if (length(sig) != 1L || !is.null(sig[[1L]][["constant"]]) ||
        !is.null(sig[[1L]][["equate"]])) return(NULL)
  y <- frame[["y"]][[lp[["resp"]]]]
  if (!is.numeric(y) || !is.null(dim(y))) return(NULL)
  z <- predictor_scale_response(fname, y)
  # an imputed value is not data (response_observed())
  ok <- is.finite(z) & response_observed(frame, lp[["resp"]])
  # a row of weight 0 is not data, so it does not set the unit either:
  # zero weights and the subset they leave fit the same model, and with
  # the unit of the full response the two took different paths
  # (test-open-issues.R, lme4#880)
  w <- frame[["aterm_values"]][[lp[["resp"]]]][["weights"]]
  if (is.numeric(w) && length(w) == length(y)) {
    ok <- ok & w > 0
  } else {
    w <- NULL
  }
  list(z = z, ok = ok, w = w)
}

#' The response of a family with an identity-link `mu` on the scale that
#' `mu` is in, `NA` on a row that does not inform `mu`: log y for the
#' lognormal families, with the shift the shifted lognormal's own start
#' uses, only the positive rows of a hurdle, and only the nonzero rows
#' of a zero-inflated family, whose zeros the inflation explains.
#'
#' @noRd
predictor_scale_response <- function(fname, y) {
  switch(fname,
    lognormal = log(y),
    shifted_lognormal = log(y - min(y) / 2),
    hurdle_lognormal = ifelse(y > 0, log(y), NA_real_),
    zero_inflated_asym_laplace = ifelse(y != 0, y, NA_real_),
    y)
}

#' The spread of `z[ok]` with the predictor's fixed design (and offset)
#' taken out by least squares, for predictor_unit_of().
#'
#' The unit is the spread of what the model leaves to sigma and the
#' random effects, not of the response: `y = 1 + k x + b_g + e` has
#' `sd(y)` near `k`, and a unit that large left fits at lme4's maximum
#' at "false convergence (8)" with a warning (k = 1e4: 33 of 40;
#' dev/optrap-rev-r2b.R, the review of lane optrap, B1). Grouping
#' structure is not taken out; its variance stays in the residual,
#' which is the scale the group sds start at.
#'
#' The spread is the MAD of the residuals (`stats::mad()`, the sd for
#' normal residuals), so a few outliers do not set it: on student fits
#' of t(1) data the residuals' sd gave a unit of 64 against a scale
#' near 1, and one seed stopped 29.8 and 39.6 short (the review's
#' dev/optrap-rev2-mad.R). Where half the residuals or more are equal,
#' the MAD is 0 and the sd is taken.
#'
#' A design with more than `cells` entries (`n * p`) is fitted on
#' evenly spaced rows up to that size: a spread from a few 1e5 rows is
#' within a percent, far inside the power of two it is rounded to.
#' Skipping the least squares there gave the raw spread of `y`, the rule
#' of B1, and a fit 101000 by 50 then ran to its evaluation limit about
#' 8.5e3 log-likelihood units short (dev/optrap-rev2-b1.R wide).
#'
#' A residual below 1e-10 of the response's magnitude is the rounding
#' of an exact fit, not a scale, and gives `NA` (no unit): test-
#' edgecases.R's constant response left residuals near 1e-16, a unit of
#' 2^-53, and a fit that stopped at code 0 where 0.69.0's ran on, both
#' toward a sigma of 0 (dev/optrap-findings.md, punch round 1). So does
#' a design whose rows do not match the response's.
#'
#' With weights the least squares is weighted (predictor_ls()), the
#' estimate of `mu`'s coefficients that the likelihood weights ask for,
#' and the spread is of its residuals.
#'
#' @noRd
predictor_resid_sd <- function(lp, z, ok, w = NULL,
                               cells = getOption("frmtmb.unit_ls_cells",
                                                 5e6),
                               fit = predictor_ls(lp, z, ok, w, cells)) {
  if (is.null(fit)) return(NA_real_)
  z <- fit[["resid"]]
  mag <- fit[["mag"]]
  s <- stats::mad(z)
  if (!is.finite(s) || !(s > 1e-10 * mag)) s <- stats::sd(z)
  if (!is.finite(s) || !(s > 1e-10 * mag)) NA_real_ else s
}

#' The least-squares fit of `z[ok]` (less the predictor's offset) on the
#' predictor's fixed design, for the unit (predictor_resid_sd()) and the
#' starting values (start_ls_plan()), so that one solve serves both.
#'
#' A list: `rows`, the indices of the rows used; `keep`, the positions
#' within `rows` the design was solved on (all of them up to `cells`
#' design entries, evenly spaced ones above); `coef`, the coefficients
#' (`NA` for an aliased column; `NULL` for a predictor without a design);
#' `resid`, the residuals on the rows kept; `mag`, the largest `|z|`
#' over `rows`. `NULL` when fewer than two rows are usable or the design
#' does not match the response's rows.
#'
#' `w`, the likelihood weights, gives the weighted fit. Equal weights
#' give the unweighted fit, so that a weight of 1 on every row changes
#' nothing.
#'
#' @noRd
predictor_ls <- function(lp, z, ok, w = NULL,
                         cells = getOption("frmtmb.unit_ls_cells", 5e6)) {
  rows <- which(ok)
  if (length(rows) < 2L) return(NULL)
  z <- z[rows]
  mag <- max(abs(z))
  X <- if (is.null(lp[["constant"]])) lp[["X"]]
  off <- lp[["offset"]]
  if (is.numeric(off) && length(off) == length(ok) && all(is.finite(off))) {
    z <- z - off[rows]
  }
  out <- list(rows = rows, keep = seq_along(rows), coef = NULL, resid = z,
              mag = mag)
  if (!design_ok(X)) return(out)
  if (nrow(X) != length(ok)) return(NULL)
  p <- ncol(X)
  m <- max(2L, floor(cells / p))
  keep <- if (length(rows) > m) {
    unique(round(seq(1, length(rows), length.out = m)))
  } else {
    seq_along(rows)
  }
  Xk <- as.matrix(X[rows[keep], , drop = FALSE])
  zk <- z[keep]
  wk <- if (is.numeric(w) && length(w) == length(ok)) w[rows[keep]]
  if (is.null(wk) || all(wk == wk[1L])) {
    q <- base::qr(Xk)
    out$resid <- as.numeric(base::qr.resid(q, zk))
    out$coef <- as.numeric(base::qr.coef(q, zk))
  } else {
    sw <- sqrt(wk)
    q <- base::qr(Xk * sw)
    out$coef <- as.numeric(base::qr.coef(q, zk * sw))
    b <- out$coef
    b[is.na(b)] <- 0
    out$resid <- as.numeric(zk - Xk %*% b)
  }
  out$keep <- keep
  out
}

#' The largest condition number among the fixed designs of a frame, each
#' with its columns scaled to unit length and its exactly aliased columns
#' left out (the QR's rank), estimated from the R factor (base::rcond()),
#' on evenly spaced rows up to `cells` as in predictor_resid_sd(). For
#' newton_polish(): raw polynomials on x in [1, 2] give 3e5 at degree 5
#' and 4e6 at degree 6, where nlminb stopped short in silence; sleepstudy
#' gives 3.4. 1 for a frame with no design.
#'
#' @noRd
predictor_design_kappa <- function(frame,
                                   cells = getOption("frmtmb.unit_ls_cells",
                                                     5e6)) {
  out <- 1
  for (lp in frame[["linpreds"]] %||% list()) {
    k <- lp_design_kappa(lp, cells)
    if (is.finite(k) && k > out) out <- k
  }
  out
}

#' One linear predictor's condition number for predictor_design_kappa(),
#' or 1 for a predictor without a design of two columns or more. The
#' least-squares start (start_ls_one()) reads the same number for `mu`,
#' so the start and newton_polish() take the same side of their bound:
#' an estimate of its own, from the start's QR, differed from this one
#' by 2.8e-7 at 1e5 and put one design on both sides (the review of lane
#' optstart, dev/optstart-rev-kappa.R).
#'
#' @noRd
lp_design_kappa <- function(lp,
                            cells = getOption("frmtmb.unit_ls_cells", 5e6)) {
  X <- if (is.null(lp[["constant"]])) lp[["X"]]
  if (!design_ok(X) || ncol(X) < 2L) return(1)
  tryCatch({
    p <- ncol(X)
    n <- nrow(X)
    m <- max(2L, floor(cells / p))
    rows <- if (n > m) unique(round(seq(1, n, length.out = m))) else
      seq_len(n)
    Xk <- as.matrix(X[rows, , drop = FALSE])
    nrm <- sqrt(colSums(Xk^2))
    Xk <- Xk[, nrm > 0, drop = FALSE] / rep(nrm[nrm > 0], each = nrow(Xk))
    q <- base::qr(Xk)
    if (q$rank < 2L) 1 else {
      R <- qr.R(q)[seq_len(q$rank), seq_len(q$rank), drop = FALSE]
      1 / base::rcond(R, triangular = TRUE)
    }
  }, error = function(e) 1)
}

#' Start every random-effect standard deviation on the scale of the
#' predictor its column enters (`predictor_units()`), so that a group
#' sd starts at the response's spread as lme4's relative-covariance
#' start does, whatever the response's units. A structure whose sds are
#' shared across columns is moved only when all of its columns agree;
#' an entry the frame maps to a constant is never moved.
#'
#' @noRd
unit_theta_start <- function(frame, tpl, units = predictor_units(frame)) {
  if (!length(units) || !length(tpl[["theta"]])) return(tpl)
  mapped <- frame[["map"]][["theta"]]
  for (bk in frame[["re_blocks"]]) {
    sdi <- tryCatch(as.integer(block_sd_idx(bk)),
                    error = function(e) integer(0))
    D <- bk[["dim"]] %||% 0L
    if (!length(sdi) || !D) next
    col_u <- rep(1, D)
    for (cp in bk[["components"]] %||% list()) {
      u <- units[cp[["lp_key"]] %||% ""]
      if (is.null(cp[["offset"]]) || is.null(cp[["dim"]]) ||
            is.na(u)) next
      cols <- cp[["offset"]] + seq_len(cp[["dim"]])
      col_u[cols[cols <= D]] <- u
    }
    shift <- if (length(sdi) == D) log(col_u) else
      if (all(col_u == col_u[1L])) rep(log(col_u[1L]), length(sdi)) else
        rep(0, length(sdi))
    pos <- bk[["theta_idx"]][sdi]
    keep <- shift != 0 & !is.na(pos)
    if (!is.null(mapped)) keep <- keep & !is.na(mapped[pos])
    tpl[["theta"]][pos[keep]] <- tpl[["theta"]][pos[keep]] + shift[keep]
  }
  tpl
}

#' The natural units of the coefficients of each predictor on the
#' response's scale (`predictor_units()`), aligned with
#' `names(obj$par)`; `NULL` when every unit is 1. A coefficient moves
#' the predictor by its column's spread per unit, so its natural unit is
#' the response's unit over that spread: `u / sd(x_j)` (each to a power
#' of two) for a numeric column with more than two values, `u` for the
#' intercept and a binary column, and for a `mo()` coefficient
#' `u / sd(codes * mult)`, the spread of the term at equal weights.
#' `scaled` marks coefficients that
#' already carry their column's unit (autoscale_units(),
#' smooth_fx_units()); those get `u` alone.
#'
#' Multiplied into `par_units`, so nlminb steps these coefficients in
#' their natural units and the standard-error check reads their Hessian
#' in them. The response's unit alone was too coarse for a large
#' column: on `ls ~ mo(income) * age` (age near 40, unit 16) the
#' interaction's coefficient, whose maximum is near 0.03, was stepped in
#' units of 16, and 6 of the 200 fits of dev/optima-mo-study.R ended in
#' the wrong sign's cone, up to 0.95 short (dev/optrap-findings.md).
#'
#' @noRd
predictor_par_units <- function(frame, par_names,
                                units = predictor_units(frame),
                                scaled = NULL) {
  if (!length(units)) return(NULL)
  unit <- lapply(frame[["par_template"]], function(v) rep(1, length(v)))
  for (key in names(units)) {
    lp <- frame[["linpreds"]][[key]]
    if (is.null(unit[[lp[["par"]]]])) next
    s <- predictor_col_spread(lp, frame)
    if (length(s) != length(lp[["idx"]])) s <- rep(1, length(lp[["idx"]]))
    unit[[lp[["par"]]]][lp[["idx"]]] <- units[[key]] / s
  }
  out <- numeric(0)
  for (cp in unique(par_names)) {
    u <- unit[[cp]]
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      u <- u[-frame[["betad_fixed_idx"]]]
    }
    out <- c(out, u)
  }
  if (length(out) != length(par_names) || all(out == 1)) return(NULL)
  # a coefficient already in its column's units keeps that and gains
  # the response's unit only
  if (length(scaled) == length(out)) {
    resp_only <- predictor_par_units_resp(frame, par_names, units)
    out[scaled] <- resp_only[scaled]
  }
  out
}

#' The response's unit alone per coefficient (no column spread), aligned
#' like predictor_par_units().
#'
#' @noRd
predictor_par_units_resp <- function(frame, par_names, units) {
  unit <- lapply(frame[["par_template"]], function(v) rep(1, length(v)))
  for (key in names(units)) {
    lp <- frame[["linpreds"]][[key]]
    if (is.null(unit[[lp[["par"]]]])) next
    unit[[lp[["par"]]]][lp[["idx"]]] <- units[[key]]
  }
  out <- numeric(0)
  for (cp in unique(par_names)) {
    u <- unit[[cp]]
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      u <- u[-frame[["betad_fixed_idx"]]]
    }
    out <- c(out, u)
  }
  out
}

#' The spread of each column of a linear predictor's design, for
#' predictor_par_units(): the sample sd of a numeric column with more
#' than two values, of `codes * mult` for a `mo()` placeholder, of the
#' observed values times `mult` for a `mi()` placeholder and of the
#' measured values times `mult` for an `me()` placeholder (`frame`
#' supplies those two), and 1 for the intercept, a binary column and
#' any other placeholder.
#'
#' A `mi()` or `me()` column is a zero placeholder in `X`, so it read as
#' spread 1, and its coefficient got the response's unit itself where its
#' natural unit is that over the spread of the variable: on brms_missings'
#' `bmi | mi() ~ age * mi(chl)` (sd of chl 45) `bmi_michl` was stepped
#' in units of 8 against a standard error of 0.034. The outer Hessian
#' the standard errors come from is taken by finite differences in those
#' units, and steps that large across the latent missing values gave a
#' Hessian that was not positive definite: every standard error NaN at
#' the optimum of 0.67.0, or the sigmas "flat" and the rest up to 35%
#' too large (lane remeasure, R1; dev/optstart-findings.md).
#'
#' @noRd
predictor_col_spread <- function(lp, frame = NULL) {
  X <- lp[["X"]]
  if (is.null(X)) return(numeric(0))
  k <- ncol(X)
  s <- rep(1, k)
  cn <- colnames(X)
  mos <- lp[["mo"]] %||% list()
  mo_cols <- vapply(mos, function(m) as.integer(m[["col"]]), 0L)
  mis <- lp[["mi"]] %||% list()
  mi_cols <- vapply(mis, function(m) as.integer(m[["col"]]), 0L)
  mes <- lp[["me"]] %||% list()
  me_cols <- vapply(mes, function(m) as.integer(m[["col"]]), 0L)
  for (j in seq_len(k)) {
    if (identical(cn[j], "(Intercept)")) next
    m <- match(j, mo_cols)
    i <- match(j, mi_cols)
    e <- match(j, me_cols)
    xj <- if (!is.na(m)) {
      tryCatch(suppressWarnings(as.numeric(mos[[m]][["codes"]]) *
                                  as.numeric(mos[[m]][["mult"]] %||% 1)),
               error = function(e) NA_real_)
    } else if (!is.na(i)) {
      placeholder_values(spread_mi_values(mis[[i]], frame), mis[[i]],
                         nrow(X))
    } else if (!is.na(e)) {
      placeholder_values(spread_me_values(mes[[e]], frame), mes[[e]],
                         nrow(X))
    } else {
      as.numeric(X[, j])
    }
    # a placeholder whose codes or multiplier do not line up with the
    # rows (a structure this function does not know) keeps spread 1,
    # the step 0.69.0 took, rather than a spread of recycled values
    if (length(xj) != nrow(X)) next
    # a mi() value is missing where it is imputed; the observed rows set
    # its spread
    if (!is.na(m) || is.na(i)) {
      if (anyNA(xj)) next
    } else {
      xj <- xj[!is.na(xj)]
    }
    if (length(unique(xj)) <= 2L) next
    sj <- stats::sd(xj)
    # a power of two, as the response's unit is, so that two data sets
    # whose columns spread alike step their coefficients alike
    if (is.finite(sj) && sj > 0) s[j] <- 2^round(log2(sj))
  }
  s
}

#' The values a `mi()` placeholder stands for, on the rows of its
#' predictor: the observed response `var` (missing where imputed), read
#' through `idxl` when the two responses do not share their rows.
#' `NULL` where `frame` does not have them.
#'
#' @noRd
spread_mi_values <- function(info, frame) {
  v <- frame[["y"]][[info[["var"]]]]
  if (!is.numeric(v) || !is.null(dim(v))) return(NULL)
  v[!response_observed(frame, info[["var"]])] <- NA
  idx <- info[["idxl"]]
  if (!is.null(idx)) {
    if (!is.numeric(idx) || anyNA(idx) || any(idx < 1) ||
          any(idx > length(v))) return(NULL)
    v <- v[idx]
  }
  as.numeric(v)
}

#' The measured values an `me()` placeholder stands for: the product of
#' its terms' observed values (`frame$me$terms[[k]]$Xn`). `NULL` where
#' `frame` does not have them.
#'
#' @noRd
spread_me_values <- function(info, frame) {
  terms <- frame[["me"]][["terms"]]
  ks <- info[["terms"]]
  if (is.null(terms) || !length(ks) || anyNA(ks)) return(NULL)
  out <- NULL
  for (k in ks) {
    xn <- terms[[k]][["Xn"]]
    if (!is.numeric(xn)) return(NULL)
    if (is.null(out)) {
      out <- as.numeric(xn)
    } else {
      if (length(xn) != length(out)) return(NULL)
      out <- out * as.numeric(xn)
    }
  }
  out
}

#' The values of a placeholder column times its multiplier, or `NA` when
#' they do not line up with the predictor's `n` rows.
#'
#' @noRd
placeholder_values <- function(v, info, n) {
  if (is.null(v) || length(v) != n) return(NA_real_)
  mult <- info[["mult"]]
  if (is.null(mult)) return(v)
  mult <- suppressWarnings(as.numeric(mult))
  if (length(mult) != n) return(NA_real_)
  v * mult
}

#' The inner (random-effect) solver's control for a fit whose response
#' units are `units` (`predictor_units()`).
#'
#' TMB's inner Newton solver returns without a step once the largest
#' absolute gradient in the random effects is below `grad.tol`, 1e-8 by
#' default, and it starts from the last solution. That gradient is per
#' unit of a random effect, so with random effects in units `u` it is
#' `1/u` as large for the same relative error: at `u = 2^20` a mode can
#' be left of the order of 1e3 units (in an sd near 1e6) from where it
#' belongs, an error of the order of 1e-5 in the Laplace objective.
#' Without this, with the response in units 1e6 times larger, 6 of 40
#' ri20 fits and 3 of 40 sleepstudy-like fits ended at a nonzero code,
#' up to 2.6e-6 short of lme4 (dev/optrap-scale.R with
#' `OPTRAP_ABLATE=inner`; dev/optrap-findings.md). The tolerance is
#' divided by the largest unit, which leaves it at TMB's default when no
#' unit is above 1.
#' `maxit` is TMB's default, which an `inner.control` list replaces.
#'
#' @noRd
predictor_inner_control <- function(units) {
  u <- max(c(1, units))
  if (u == 1) return(list(maxit = 1000))
  list(maxit = 1000, grad.tol = 1e-8 / u)
}
