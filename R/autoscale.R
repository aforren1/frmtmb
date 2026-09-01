# Internal predictor standardization (the lme4 >= 1.1.37 autoscale
# idea). Instead of rescaling accessors post hoc, exploit cheap
# re-taping: fit a standardized copy of the frame, map that optimum
# back to the original parameterization exactly (ML is invariant under
# affine reparameterization of X columns), and use it to warm-start
# the ordinary unscaled fit. Everything downstream - sdreport, vcov,
# predict, profile - sees only the original-scale fit.

# Which columns qualify: parametric X columns (1..n_param_cols) that
# are numeric with more than two distinct values. That excludes the
# intercept, factor dummies and other binary contrasts, smooth basis
# columns (appended after n_param_cols), and the zero-placeholder
# mo()/mi() columns. Centering needs an intercept column to absorb the
# shift, so an intercept-free linpred is scaled without centering.
# Returns NULL when nothing qualifies anywhere.
autoscale_plan <- function(frame) {
  plan <- list()
  for (key in names(frame$linpreds)) {
    lp <- frame$linpreds[[key]]
    if (!is.null(lp$nl_body)) next     # no design matrix of its own
    if (!is.null(lp$constant)) next    # mapped to a fixed value
    if (is.null(lp$X) || lp$n_param_cols == 0L) next
    icpt <- match("(Intercept)",
                  colnames(lp$X)[seq_len(lp$n_param_cols)])
    cols <- integer(0)
    center <- numeric(0)
    scale <- numeric(0)
    for (j in seq_len(lp$n_param_cols)) {
      if (!is.na(icpt) && j == icpt) next
      xj <- lp$X[, j]
      if (length(unique(xj)) <= 2L) next
      s <- stats::sd(xj)
      if (!is.finite(s) || s == 0) next
      cols <- c(cols, j)
      center <- c(center, if (is.na(icpt)) 0 else mean(xj))
      scale <- c(scale, s)
    }
    if (!length(cols)) next
    plan[[key]] <- list(cols = cols, center = center, scale = scale,
                        icpt = icpt, par = lp$par, idx = lp$idx)
  }
  if (!length(plan)) NULL else plan
}

# Standardized copy of the frame: only the planned X columns change.
autoscale_frame <- function(frame, plan) {
  for (key in names(plan)) {
    pl <- plan[[key]]
    X <- frame$linpreds[[key]]$X
    Xs <- sweep(X[, pl$cols, drop = FALSE], 2, pl$center)
    X[, pl$cols] <- sweep(Xs, 2, pl$scale, "/")
    frame$linpreds[[key]]$X <- X
  }
  frame
}

# Map a parameter list between the parameterizations. Column j scaled
# as (x - c_j)/s_j gives scaled coefficient b_j * s_j, and the linpred
# intercept absorbs sum_j b_j * c_j (c_j is zero without an intercept,
# so absorption vanishes there). theta, b, and extras are untouched:
# affine column changes leave them identical at the optimum.
autoscale_map <- function(pars, plan, to = c("original", "scaled")) {
  to <- match.arg(to)
  for (pl in plan) {
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

# Stage one of frmtmb_control(autoscale = TRUE): fit the standardized
# frame and return its optimum mapped back to the original
# parameterization, as the warm-start template for the real fit. User
# starts are translated into the scaled parameterization so they mean
# the same model.
autoscale_prefit <- function(spec, frame, bform, cl, REML, start,
                             control, lower, upper, priors, quadrature,
                             plan) {
  ctl <- control
  ctl$autoscale <- FALSE
  if (!is.null(start)) {
    full <- frame$par_template
    for (nm in intersect(names(start), names(full))) {
      full[[nm]][] <- start[[nm]]
    }
    full <- autoscale_map(full, plan, "scaled")
    for (nm in intersect(names(start), c("beta", "betad"))) {
      start[[nm]] <- unname(full[[nm]])
    }
  }
  sfit <- fit_assembled(spec, autoscale_frame(frame, plan), bform, cl,
                        REML = REML, start = start, control = ctl,
                        se = FALSE, lower = lower, upper = upper,
                        priors = priors, quadrature = quadrature)
  autoscale_map(sfit$estimates, plan, "original")
}

# sdreport on the unscaled optimum needs two repairs when non-unit
# autoscale units are in play: optimHess steps are absolute (far larger
# than a 1e-6-magnitude coefficient, so the default outer Hessian is
# NaN or garbage), and no solver can invert a Hessian whose scaled-in
# condition number exceeds double precision. Both problems vanish in
# the coordinates q = par / unit - the scaled problem - so compute the
# Hessian there, hand its exact back-map to sdreport (keeping the
# joint precision and diag.cov.random finite), and overwrite cov.fixed
# with the exactly transformed well-conditioned inverse.
autoscale_sdreport <- function(fit, jp = needs_jp(fit)) {
  u <- fit$par_units
  if (is.null(u) || all(u == 1)) {
    return(RTMB::sdreport(fit$obj, getJointPrecision = jp))
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
  sdr
}

# Natural magnitude per outer parameter of the warm-started unscaled
# fit, aligned with names(obj$par): 1/s_j for the coefficient of a
# column with sample SD s_j, 1 elsewhere. Multiplying the gradient by
# these units reproduces the well-conditioned scaled-fit gradient.
autoscale_units <- function(frame, plan, par_names) {
  unit <- lapply(frame$par_template, function(v) rep(1, length(v)))
  for (pl in plan) {
    unit[[pl$par]][pl$idx[pl$cols]] <- 1 / pl$scale
  }
  out <- numeric(0)
  for (cp in unique(par_names)) {
    u <- unit[[cp]]
    if (cp == "betad" && length(frame$betad_fixed_idx)) {
      u <- u[-frame$betad_fixed_idx]
    }
    out <- c(out, u)
  }
  stopifnot(length(out) == length(par_names))
  out
}
