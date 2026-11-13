# Restarts and stop checks for fits that can end at a lower maximum with
# nothing in the fit to mark it: a prior in conflict with the data, and a
# nonlinear body on a ridge that runs to a limit.

#' How many of its scales a coefficient's estimate must lie from the
#' location of its prior before prior_conflict_escape() restarts there.
#'
#' A normal likelihood with a free scale and a normal prior on a
#' coefficient is bimodal only when the estimate lies far beyond the
#' prior's scale: one mode near the data's value with a small residual
#' scale, one near the prior's location with the residual scale inflated.
#'
#' @noRd
prior_conflict_z <- 4

#' Restart the optimizer with each fixed coefficient that is in conflict
#' with its prior at that prior's location, and keep the better optimum.
#'
#' `y ~ x` at y about 1e6 with `normal(0, 1000)` on the slope has two
#' posterior modes: the slope near the data's 5e4 with sigma near 1e3,
#' and the slope near 4400 with sigma near 4.6e4, 470 higher in the
#' penalized log-likelihood. The default start is the least-squares
#' solution, so the fit took the lower mode in silence, where the body
#' spelling `a + b * x`, whose start is the prior's location, took the
#' higher one (dev/nlstop-prior.R, seeds 6201 to 6203). The same holds for
#' a prior on the centered intercept: `normal(0, 1e4)` on it stopped
#' 3613.6 below the mode near 0 (review of lane nlstop, MAJOR 2). A
#' coefficient is in conflict when its estimate, or the centered
#' intercept a prior reads (`offset`), lies more than `prior_conflict_z`
#' scales of a normal, Student t or logistic prior from its location. The
#' restart begins at the tape's own start with those coefficients moved,
#' a centered intercept moved so that the centered value is the location.
#'
#' It runs only where the second mode can be: the likelihood has a free
#' dispersion or scale in the outer vector (free_dispersion()). A
#' logistic GLMM of 20000 rows with five slopes in conflict paid 2.1 times
#' its fit for a restart that never won, and refit() 2.5 times (review,
#' MAJOR 1). The caller runs it on a cold start only (frm()'s own call,
#' no warm template), so refits, influence(), frm_allfit(), anova()'s ML
#' refit and the autoscale pre-fit keep the mode they are handed.
#' Nonlinear coefficients already start at their prior's location
#' (prior_nl_starts()), and a caller's `start` is kept as given.
#'
#' Returns the optimum to keep, with `prior_escape` holding the two
#' objective values when the restart won by more than nl_probe_gain().
#'
#' @noRd
prior_conflict_escape <- function(obj, opt, frame, entries, control,
                                  bounds, par_units = NULL, tally = NULL) {
  if (!length(entries) || is.function(control$optimizer) ||
        !identical(control$optimizer %||% "nlminb", "nlminb") ||
        optctrl_stops_earlier(control$optCtrl) ||
        !is.finite(opt$objective)) {
    return(opt)
  }
  pn <- names(obj$par)
  bpos <- which(pn == "beta")
  if (!length(bpos) || !free_dispersion(frame) ||
        length(bpos) != length(frame[["par_template"]][["beta"]])) {
    return(opt)
  }
  nl_pos <- nl_beta_positions(frame)
  # the centered value a prior with an offset reads, at outer vector p
  off_at <- function(e, p) {
    o <- e$offset
    (e[["shift"]] %||% 0) +
      if (is.null(o)) 0 else sum(p[bpos[o$idx]] * o$w)
  }
  hits <- list()
  for (e in entries) {
    o <- e$offset
    if (!identical(e$comp, "beta") ||
          !identical(e$scale %||% "internal", "internal") ||
          !is.null(e$link) ||
          (!is.null(o) && (!identical(o$comp, "beta") ||
                             length(o$w) != length(o$idx)))) {
      next
    }
    d <- e$dist
    if (!(d$kind %||% "") %in% c("normal", "t", "logistic")) next
    m <- prior_dist_location(d)
    s <- d$scale
    if (is.na(m) || !is.numeric(s) || length(s) != 1L || !is.finite(s) ||
          s <= 0) {
      next
    }
    for (i in setdiff(e$idx, nl_pos)) {
      v <- opt$par[[bpos[i]]] + off_at(e, opt$par)
      if (isTRUE(abs(v - m) > prior_conflict_z * s)) {
        hits[[length(hits) + 1L]] <- list(j = bpos[i], loc = m, e = e)
      }
    }
  }
  if (!length(hits)) return(opt)
  q <- obj$par
  # plain coefficients first, then each centered value from the moved q
  plain <- vapply(hits, function(h) is.null(h$e$offset), NA)
  for (h in hits[plain]) q[h$j] <- h$loc
  for (h in hits[!plain]) q[h$j] <- h$loc - off_at(h$e, q)
  q <- pmin(pmax(q, bounds$lower), bounds$upper)
  saved <- obj_state_save(obj)
  o <- tryCatch(optimize_obj(obj, control, bounds, par_units, start_par = q,
                             tally = tally),
                error = keep_missing(NULL))
  if (is.null(o) || !is.finite(o$objective) ||
        !(opt$objective - o$objective > nl_probe_gain(opt$objective))) {
    obj_settle(obj, opt, opt, saved)
    return(opt)
  }
  obj_settle(obj, opt, o, saved)
  at <- vapply(hits, `[[`, 0L, "j")
  nms <- par_template_names(frame[["par_template"]][["beta"]], "beta")
  nms <- ifelse(plain, nms[match(at, bpos)], "the centered Intercept")
  o$prior_escape <- list(from = opt$objective, to = o$objective,
                         names = nms)
  o
}

#' Does a distributional parameter that is a dispersion or a scale have a
#' coefficient in the outer vector? A probability (zi, hu, zoi, coi, a
#' mixture weight theta<k>) is not one: it does not widen the likelihood
#' around a coefficient the prior pulls away from the data, and
#' zero_inflated_poisson() took the restart on its zi alone (short check
#' of the review of lane nlstop, S2). A parameter fixed by `dpar = value`
#' has none.
#'
#' @noRd
free_dispersion <- function(frame) {
  fixed <- frame[["betad_fixed_idx"]] %||% integer(0)
  prob <- c("zi", "hu", "zoi", "coi", "theta")
  for (lp in frame[["linpreds"]]) {
    if (identical(lp[["par"]], "betad") &&
          is.null(lp[["constant"]]) &&
          !sub("[0-9]+$", "", lp[["dpar"]] %||% "") %in% prob &&
          length(setdiff(lp[["idx"]], fixed))) {
      return(TRUE)
    }
  }
  FALSE
}

#' The message of a fit that prior_conflict_escape() moved to another
#' mode, on the user's own call.
#'
#' @noRd
prior_conflict_message <- function(opt) {
  pe <- opt$prior_escape
  if (is.null(pe)) return(invisible(NULL))
  ll <- nl_ll_pair(-pe$from, -pe$to)
  frm_message(
    "The fit from the default start reached penalized log-likelihood ",
    ll[1L], ", with ", paste(unique(pe$names), collapse = ", "),
    " more than ", prior_conflict_z, " prior scales from the prior's ",
    "location. A start at that location reached ", ll[2L], ", and frm() ",
    "returns that fit: the prior and the data disagree, and the posterior ",
    "has more than one mode. Give `start` to choose your own",
    class = "frmtmb_prior_multistart")
  invisible(NULL)
}

#' The factor by which nl_limit_message() scales the coefficients along
#' the log-likelihood's flattest direction in their log magnitudes.
#'
#' A ridge that runs to a limit of a body is nearly flat over any factor:
#' on the bernoulli Michaelis-Menten row of the review's fresh set (seed
#' 55, dev/nlstop-mm55.R) Vm and K scaled together by 0.5 to 1e4 moved the
#' log-likelihood by at most 1.7e-5. A coefficient that is merely small,
#' b2 near 0 in b1 * exp(b2 * x), is moved by this factor to a value the
#' data do read.
#'
#' @noRd
nl_limit_factor <- 1000

#' Warn when a nonlinear fit that nothing else has warned about is on a
#' ridge along which its coefficients are not determined.
#'
#' A body such as `Vm * xp / (K + xp)` approaches `(Vm / K) * xp` as Vm and
#' K grow together, so the likelihood has no maximum there, only a limit,
#' and the optimizer stops anywhere along the ridge with a small gradient
#' and finite standard errors. On the review's fresh rows (bernoulli, seed
#' 55) the fit returned Vm = 1.3e4 and K = 1.8e5 at code 0 with no warning,
#' 1.7e-5 below the linear model's log-likelihood. Such ridges are
#' straight in the log magnitudes of the coefficients, (1, 1) there and (1,
#' -1) for `b1 * (1 - exp(-b2 * x))`, so the check takes the Hessian of the
#' objective over the body's nonzero coefficients in those coordinates
#' (2 k gradients), and moves along its flattest direction by
#' `nl_limit_factor` in the coefficient that moves most, both ways (2
#' objective evaluations). If every coefficient the move carries is at
#' least one conditional standard error from 0, and a move that makes
#' some coefficient larger lowers the log-likelihood by no more than
#' `grad_tol`, the data do not
#' determine the coefficients, and it says so. On mm 55 Vm and K times
#' 1000 raised it by 1.7e-5, and divided by 1000, where the body bends
#' within the data, lowered it by 1.8e-2 (dev/nlstop-limit.R).
#'
#' Returns the warning's text, or NULL.
#'
#' @noRd
nl_limit_message <- function(fit) {
  obj <- fit$obj
  frame <- fit$frame
  if (is.null(obj) || is.null(fit$opt) || !is.finite(fit$opt$objective) ||
        isTRUE(fit$quadrature) || !is.null(fit$importance) ||
        !is.null(fit$beta_whiten) || isTRUE(fit$control$profile)) {
    return(NULL)
  }
  pn <- names(obj$par)
  bpos <- which(pn == "beta")
  if (length(bpos) != length(frame[["par_template"]][["beta"]])) return(NULL)
  nl_pos <- nl_beta_positions(frame)
  p0 <- fit$opt$par
  j <- bpos[nl_pos]
  j <- j[is.finite(p0[j]) & p0[j] != 0]
  if (!length(j)) return(NULL)
  tol <- fit$control$grad_tol %||% 1e-3
  saved <- obj_state_save(obj)
  on.exit({
    obj_state_restore(obj, saved)
    try(obj$gr(p0), silent = TRUE)
  }, add = TRUE)
  g0 <- tryCatch(drop(obj$gr(p0))[j], error = keep_missing(NULL))
  if (is.null(g0) || !all(is.finite(g0))) return(NULL)
  H <- tryCatch({
    H <- matrix(NA_real_, length(j), length(j))
    for (i in seq_along(j)) {
      h <- 1e-4 * abs(p0[j[i]])
      up <- p0
      dn <- p0
      up[j[i]] <- up[j[i]] + h
      dn[j[i]] <- dn[j[i]] - h
      H[, i] <- (drop(obj$gr(up))[j] - drop(obj$gr(dn))[j]) / (2 * h)
    }
    H
  }, error = keep_missing(NULL))
  if (is.null(H) || !all(is.finite(H))) return(NULL)
  # the Hessian in u = log|theta|: d theta = theta du
  D <- p0[j]
  Hu <- (H + t(H)) / 2 * outer(D, D) + diag(D * g0, length(j))
  e <- eigen(Hu, symmetric = TRUE)
  v <- e$vectors[, length(j)]
  big <- abs(v) >= 0.1 * max(abs(v))
  # each coefficient the move carries must be one the objective reads
  # where it is, at least one conditional standard error from 0: a prior
  # that holds `a` near 0 in `a * xs`, xs about 1e-20 (test-nl-start.R),
  # lets a 1000-fold move change nothing, and that is the prior's answer,
  # not a ridge
  if (!all(abs(D[big]) * sqrt(abs(diag(H))[big]) >= 1)) return(NULL)
  lf <- log(nl_limit_factor)
  # a move whose quadratic model is far above grad_tol is no ridge
  if (!isTRUE(0.5 * e$values[length(j)] * (lf / max(abs(v)))^2 <=
                1e3 * tol)) {
    return(NULL)
  }
  t <- lf / max(abs(v))
  f <- function(s) {
    q <- p0
    q[j] <- p0[j] * exp(s * t * v)
    tryCatch(as.numeric(se_quiet(obj$fn(q))), error = function(e) NA_real_)
  }
  side <- c(1, -1)
  gain <- fit$opt$objective - c(f(1), f(-1))
  # a side counts where the log-likelihood does not fall by more than
  # grad_tol and some coefficient grows: moving a coefficient near 0
  # toward 0 changes nothing either, and is no ridge
  grows <- vapply(side, function(s) any(s * v[big] > 0), NA)
  ok <- is.finite(gain) & gain >= -tol & grows
  if (!any(ok)) return(NULL)
  k <- which(ok)[which.max(gain[ok])]
  rises <- gain[k] > tol
  nm <- par_template_names(frame[["par_template"]][["beta"]],
                           "beta")[match(j, bpos)]
  fac <- exp(side[k] * t * v[big])
  say <- paste0(nm[big], " x ", signif(fac, 2), collapse = ", ")
  paste0("The data do not determine the coefficients ",
         paste(nm[big], collapse = ", "), " of the nonlinear body: ",
         "moved together (", say, ") they ",
         if (rises) {
           paste0("raise the log-likelihood by ", signif(gain[k], 2),
                  ", so the fit is not at a maximum")
         } else {
           paste0("change the log-likelihood by ", signif(abs(gain[k]), 2),
                  ", within grad_tol (", format(tol), ")")
         },
         ". The fit is on a ridge that can run to a limit of the body, ",
         "where it reduces to a simpler form (Vm * x / (K + x), as Vm and ",
         "K grow together, reaches (Vm / K) * x). The estimates are one ",
         "point on it, and their standard errors describe no maximum. Fit ",
         "the limit's form, or put a prior on one of the coefficients")
}
