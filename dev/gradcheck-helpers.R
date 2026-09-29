# Scoring helpers for the gradient-check lane. Nothing here changes the
# package: every candidate criterion is computed from a finished fit, so
# the reference build and the lane build can be scored by the same code.

gc_libpaths <- function(which = c("lane", "base")) {
  which <- match.arg(which)
  user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
  first <- if (identical(which, "lane")) {
    "C:/Users/adf44/source/r/wt-gradcheck-lib"
  } else {
    "C:/Users/adf44/source/r/rellib-r3"
  }
  .libPaths(c(first, user))
  invisible(.libPaths())
}

# The box the optimizer actually ran under, on the outer parameter
# vector. `fit$lower`/`fit$upper` hold only the caller's bounds; a bound
# spelled on a prior (`set_prior(ub = )`) is merged into the box AFTER
# those slots are filled, so reading the slots alone misses exactly the
# bounds the arcov lane reported false alarms on.
gc_bounds <- function(fit) {
  shim <- list(frame = fit$frame, REML = fit$REML, control = fit$control)
  lower <- fit$lower
  upper <- fit$upper
  if (!is.null(fit$prior)) {
    ri <- frmtmb:::resolve_prior_input(
      list(frame = fit$frame, spec = fit$spec), fit$prior)
    if (length(ri$lower)) {
      lower <- unlist(utils::modifyList(as.list(ri$lower),
                                        as.list(lower %||% c())))
    }
    if (length(ri$upper)) {
      upper <- unlist(utils::modifyList(as.list(ri$upper),
                                        as.list(upper %||% c())))
    }
  }
  frmtmb:::resolve_bounds(shim, lower, upper)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# One row of evidence per fit: the criterion the package uses today and
# every candidate, on the same gradient evaluation.
gc_probe <- function(fit, hess = TRUE) {
  p <- fit$opt$par
  nm <- frmtmb:::outer_par_names(fit)
  u <- fit$par_units %||% rep(1, length(p))
  g <- tryCatch(drop(fit$obj$gr(p)) * u, error = function(e) NULL)
  if (is.null(g) || !length(p)) {
    return(list(np = length(p), gmax = NA_real_, gmax_par = NA_character_,
                gproj = NA_real_, nactive = 0L, decr = NA_real_,
                gscaled = NA_real_, grel = NA_real_,
                obj = fit$opt$objective, conv = fit$opt$convergence))
  }
  bd <- tryCatch(gc_bounds(fit), error = function(e) NULL)
  lo <- if (is.null(bd)) rep(-Inf, length(p)) else bd$lower
  hi <- if (is.null(bd)) rep(Inf, length(p)) else bd$upper
  # KKT: a parameter pinned at a bound whose gradient pushes it further
  # out cannot move, so its gradient component is not evidence about
  # convergence. The tolerance on "at the bound" is relative to the
  # bound's own magnitude, so nothing here is an absolute constant.
  atol <- 1e-8 * pmax(1, abs(ifelse(is.finite(lo), lo, 0)),
                      abs(ifelse(is.finite(hi), hi, 0)))
  at_lo <- is.finite(lo) & (p - lo) <= atol
  at_hi <- is.finite(hi) & (hi - p) <= atol
  # nll is minimized, so at a lower bound the feasible direction is +
  # and the component is inactive when the gradient is positive
  active <- (at_lo & g > 0) | (at_hi & g < 0)
  gp <- g
  gp[active] <- 0
  free <- which(!active)
  decr <- NA_real_
  gsc <- NA_real_
  if (hess && length(free)) {
    H <- tryCatch(
      stats::optimHess(p, function(q) fit$obj$fn(q),
                       function(q) drop(fit$obj$gr(q))),
      error = function(e) NULL)
    if (!is.null(H)) {
      H <- (H + t(H)) / 2
      # the gradient was measured in natural units, so the Hessian must
      # be too: H_natural = D H D with D = diag(par_units)
      H <- diag(u, nrow = length(u)) %*% H %*% diag(u, nrow = length(u))
      Hf <- H[free, free, drop = FALSE]
      ch <- tryCatch(chol(Hf), error = function(e) NULL)
      if (!is.null(ch)) {
        z <- backsolve(ch, gp[free], transpose = TRUE)
        decr <- 0.5 * sum(z^2)
      }
      d <- diag(Hf)
      if (all(is.finite(d)) && all(d > 0)) {
        gsc <- max(abs(gp[free]) / sqrt(d))
      }
    }
  }
  list(np = length(p),
       gmax = max(abs(g)), gmax_par = nm[which.max(abs(g))],
       gproj = if (length(free)) max(abs(gp)) else 0,
       nactive = sum(active),
       decr = decr, gscaled = gsc,
       grel = max(abs(g)) / max(1, abs(fit$opt$objective)),
       obj = fit$opt$objective, conv = fit$opt$convergence)
}

# Did the fit warn the way the package does today?
gc_catch <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                class = "gc_error")),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = w)
}

gc_warned_grad <- function(w) {
  any(grepl("Large maximum absolute gradient", w, fixed = TRUE))
}

gc_row <- function(label, seed, fit, warns, extra = list()) {
  pr <- gc_probe(fit)
  out <- c(list(label = label, seed = seed,
                warned = gc_warned_grad(warns)), pr, extra)
  out
}

gc_print <- function(row) {
  cat(sprintf(
    "%-34s seed %-5s np %3d conv %d obj %14.6f gmax %10.3e (%s)\n",
    row$label, as.character(row$seed), row$np, row$conv, row$obj,
    row$gmax, row$gmax_par))
  cat(sprintf(
    paste0("%-34s warned %-5s active %2d gproj %10.3e decr %10.3e",
           " gsc %10.3e grel %10.3e\n"),
    "", row$warned, row$nactive, row$gproj, row$decr, row$gscaled,
    row$grel))
}
