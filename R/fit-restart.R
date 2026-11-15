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
#' `nl_limit_factor` in the coefficient that moves most, both ways, each
#' move read as it lands and after one Newton step across the stiff
#' directions (4 objective evaluations and 2 gradients). If every
#' coefficient the move carries is at least one conditional standard error
#' from 0, and a move that makes some coefficient larger lowers the
#' log-likelihood by no more than `grad_tol`, the data do not determine the
#' coefficients, and it says so. On mm 55 Vm and K times 1000 raised it by
#' 1.7e-5, and divided by 1000, where the body bends within the data,
#' lowered it by 1.8e-2 (dev/nlstop-limit.R).
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
  val <- function(x) {
    tryCatch(as.numeric(se_quiet(obj$fn(x))), error = function(e) NA_real_)
  }
  # the objective after a move by t along w in the log magnitudes, or
  # after one Newton step across w from there, whichever is lower: the
  # flattest direction of a ridge that bends leaves it after a 1000-fold
  # move (on a * (exp(b * x) - 1) the move lost 0.0200 where the ridge
  # itself, (1000, 1/1000), gains 0.0033; review of lane nlstop, MINOR 1),
  # and a step on the curvature at the fit comes back to it
  f <- function(w, t) {
    q <- p0
    q[j] <- p0[j] * exp(t * w)
    f1 <- val(q)
    if (length(j) < 2L) return(f1)
    N <- qr.Q(qr(matrix(w)), complete = TRUE)[, -1L, drop = FALSE]
    A <- crossprod(N, Hu %*% N)
    if (!all(eigen(A, symmetric = TRUE, only.values = TRUE)$values > 0)) {
      return(f1)
    }
    gq <- tryCatch(drop(se_quiet(obj$gr(q)))[j], error = function(e) NULL)
    if (is.null(gq) || !all(is.finite(gq))) return(f1)
    q2 <- q
    q2[j] <- q[j] * exp(-drop(N %*% solve(A, crossprod(N, q[j] * gq))))
    f2 <- val(q2)
    if (!is.finite(f1) || (is.finite(f2) && f2 < f1)) f2 else f1
  }
  w <- v / max(abs(v))
  t <- log(nl_limit_factor)
  # a move whose quadratic model is far above grad_tol is no ridge
  if (!isTRUE(0.5 * sum(w * (Hu %*% w)) * t^2 <= 1e3 * tol)) return(NULL)
  side <- c(1, -1)
  gain <- fit$opt$objective - c(f(w, t), f(w, -t))
  # a side counts where the log-likelihood does not fall by more than
  # grad_tol and some coefficient grows: moving a coefficient near 0
  # toward 0 changes nothing either, and is no ridge
  grows <- vapply(side, function(s) any(s * w[big] > 0), NA)
  ok <- is.finite(gain) & gain >= -tol & grows
  if (!any(ok)) return(NULL)
  k <- which(ok)[which.max(gain[ok])]
  best <- list(gain = gain[k], fac = exp(side[k] * t * w[big]))
  rises <- best$gain > tol
  nm <- par_template_names(frame[["par_template"]][["beta"]],
                           "beta")[match(j, bpos)]
  say <- paste0(nm[big], " x ", signif(best$fac, 2), collapse = ", ")
  paste0("The data do not determine the coefficients ",
         paste(nm[big], collapse = ", "), " of the nonlinear body: ",
         "moved together (", say, ") they ",
         if (rises) {
           paste0("raise the log-likelihood by ", signif(best$gain, 2),
                  ", so the fit is not at a maximum")
         } else {
           paste0("change the log-likelihood by ", signif(abs(best$gain), 2),
                  ", within grad_tol (", format(tol), ")")
         },
         ". The fit is on a ridge that can run to a limit of the body, ",
         "where it reduces to a simpler form (Vm * x / (K + x), as Vm and ",
         "K grow together, reaches (Vm / K) * x). The estimates are one ",
         "point on it, and their standard errors describe no maximum. Fit ",
         "the limit's form, or put a prior on one of the coefficients")
}

#' The denominators of a body: the second operand of each `/`, and the
#' base of each `^` with a negative literal power.
#'
#' @noRd
nl_denoms <- function(e) {
  if (!is.call(e)) return(list())
  out <- unlist(lapply(as.list(e)[-1L], nl_denoms), recursive = FALSE)
  if (length(e) != 3L) return(out)
  f <- e[[1L]]
  neg <- function(p) {
    (is.numeric(p) && length(p) == 1L && p < 0) ||
      (is.call(p) && identical(p[[1L]], as.name("-")) && length(p) == 2L &&
         is.numeric(p[[2L]]) && p[[2L]] > 0)
  }
  if (identical(f, as.name("/"))) out <- c(out, list(e[[3L]]))
  if (identical(f, as.name("^")) && neg(e[[3L]])) out <- c(out, list(e[[2L]]))
  out
}

#' Where a body has poles in one coefficient, the rest held at `tpl`.
#'
#' A denominator affine in the coefficient `c` is 0 on row i at one value,
#' its wall `c_i`: `Vm * xp / (K + xp)` has a pole at `K = -xp_i` for each
#' row. The likelihood cannot be continued across a wall, so the walls cut
#' the line of `c` into pieces, and an optimizer stays in the piece it
#' starts in. Returns, for each such denominator with two or more distinct
#' walls, the walls' range `lo`, `hi`, the distance `dist` from `c` to the
#' nearest wall and the median gap between adjacent walls; `NULL` where no
#' denominator reads `c` affinely. A denominator that is not affine in `c`
#' (`EC50^h + x^h`, `1 + exp(b * x)`) is not read.
#'
#' @noRd
nl_pole_walls <- function(frame, lp, tpl, par, at, dens) {
  c0 <- tpl[[par]][at]
  if (!is.finite(c0)) return(NULL)
  h <- 1e-3 * max(1, abs(c0))
  at_c <- function(v) {
    t <- tpl
    t[[par]][at] <- v
    tryCatch(suppressWarnings(nl_body_at(frame, lp, t, exprs = dens)),
             error = function(e) NULL)
  }
  d0 <- at_c(c0)
  dp <- at_c(c0 + h)
  dm <- at_c(c0 - h)
  if (is.null(d0) || is.null(dp) || is.null(dm)) return(NULL)
  out <- list()
  for (k in seq_along(dens)) {
    a <- d0[[k]]
    b <- dp[[k]]
    m <- dm[[k]]
    if (!length(a) || length(b) != length(a) || length(m) != length(a) ||
          !all(is.finite(c(a, b, m)))) {
      next
    }
    s <- (b - m) / (2 * h)
    sc <- max(abs(b - m))
    if (!(sc > 0) || max(abs(b - 2 * a + m)) > 1e-8 * sc) next
    reads <- abs(s) > 1e-12 * max(abs(s))
    w <- sort(unique(c0 - a[reads] / s[reads]))
    if (length(w) < 2L) next
    gaps <- diff(w)
    gaps <- gaps[gaps > 1e-9 * (w[length(w)] - w[1L])]
    if (!length(gaps)) next
    out[[length(out) + 1L]] <- list(lo = w[1L], hi = w[length(w)],
                                    dist = min(abs(c0 - w)),
                                    gap = stats::median(gaps),
                                    mixed = any(a[reads] > 0) &&
                                      any(a[reads] < 0))
  }
  if (length(out)) out
}

#' Does a band of walls hold the fit: does the denominator take both signs
#' over the rows, so that a pole lies between rows of the data, or is the
#' fit nearer a wall than the walls' median gap?
#'
#' On the fresh rows the fits held by a band (binomial seeds 11, 33 and 44)
#' sit 0.32 to 0.53 median gaps from a wall, and the fits clear of it 68 to
#' 669 (gaussian, poisson, gamma, cumulative; dev/nlstop2-polegeo.R). A fit
#' between rows' poles can sit further from a wall than one gap: the
#' pieces between the walls of K + xp on those rows are up to 10.4 median
#' gaps wide (dev/nlstop2-gaps.R). Whether `c` is between the first and
#' last wall does not decide it: the walls of `1 + b * x` with x of both
#' signs lie on both sides of every b, the fit sits in the wide piece
#' between the two groups, 7.5 to 20 gaps from a wall, and its denominator
#' is positive on every row.
#'
#' @noRd
nl_pole_pressed <- function(wl) {
  isTRUE(wl$mixed) || wl$dist < wl$gap
}

#' Optimize a nonlinear fit again from beyond the poles of its bodies.
#'
#' On the review's fresh rows (binomial `Vm * xp / (K + xp)`, xp in 0.1
#' to 5) the default fit of seed 11 stopped at code 0 with no warning at K
#' -0.1825, between the walls at -0.18 and -0.183 of rows inside the data,
#' 4.32 below the maximum at K -6.0152; seed 44 stopped at K -0.1103, 0.004
#' beside the wall of the smallest xp, 0.92 below the maximum at K
#' -6.4531 (dev/reviews/nlstop-review-pole.R). No restart the fit makes
#' crosses a wall: the sign probe starts K at -1, inside the band of walls
#' from -5 to -0.1. Where a band of walls of some coefficient holds the fit
#' (nl_pole_walls(), nl_pole_pressed()), the optimizer runs
#' from each side of the band the fit is not beyond, half the band's width
#' past its last wall, the rest of the model at the fit. A fit that is
#' clear of every band pays only the evaluations of its denominators.
#'
#' A run is kept if it ends clear of the band it started beside and higher
#' than the fit by more than nl_probe_gain(); then frm() fits again from
#' its point by `run()`, as the sign probe does, and that fit must end
#' clear of the band too. A run that ends against a wall is not kept:
#' bernoulli seed 55 rises toward the wall of its largest xp, where the
#' likelihood has a supremum and no maximum.
#'
#' Not run under `frmtmb_control(nl_starts = FALSE)`, which turns off this
#' restart and the sign probe.
#'
#' Returns `NULL`, or a list as nl_sign_probe() returns it.
#'
#' @noRd
nl_pole_escape <- function(run, frame, fit) {
  obj <- fit$obj
  ctl <- fit$control
  if (isFALSE(ctl$nl_starts) ||
        is.null(obj) || is.null(fit$opt) || !is.finite(fit$opt$objective) ||
        isTRUE(fit$quadrature) || !is.null(fit$importance) ||
        !is.null(fit$beta_whiten) || isTRUE(ctl$profile) ||
        is.function(ctl$optimizer) ||
        !identical(ctl$optimizer %||% "nlminb", "nlminb") ||
        optctrl_stops_earlier(ctl$optCtrl)) {
    return(NULL)
  }
  pn <- names(obj$par)
  est <- fit$estimates
  starts <- list()
  checks <- list()
  for (lp in frame[["linpreds"]]) {
    if (is.null(lp[["nl_body"]])) next
    dens <- nl_denoms(lp[["nl_body"]])
    if (!length(dens)) next
    bc <- nl_body_coefs(frame, lp)
    for (par in names(bc$coefs)) {
      pos <- which(pn == par)
      if (length(pos) != length(frame[["par_template"]][[par]])) return(NULL)
      for (at in bc$coefs[[par]]) {
        wls <- nl_pole_walls(frame, lp, est, par, at, dens)
        c0 <- est[[par]][at]
        for (wl in wls) {
          if (!nl_pole_pressed(wl)) next
          half <- (wl$hi - wl$lo) / 2
          to <- c(if (c0 >= wl$lo) wl$lo - half, if (c0 <= wl$hi) wl$hi + half)
          for (v in to) {
            starts[[length(starts) + 1L]] <- list(j = pos[at], v = v)
            checks[[length(checks) + 1L]] <- list(lp = lp, par = par,
                                                  at = at, dens = dens)
          }
        }
      }
    }
  }
  if (!length(starts)) return(NULL)
  pu <- fit$par_units
  bounds <- fit_outer_box(fit)
  saved <- obj_state_save(obj)
  ob_fit <- obj$env$.frmtmb_opt_basis
  # a basis from nl_basis_descent() holds at the first fit's stop only
  if (isTRUE(ob_fit$nl)) {
    obj$env$.frmtmb_opt_basis <- opt_basis(frame, obj, ctl$optimizer)
  }
  tpl_at <- function(par_vec) {
    st <- list()
    for (cp in unique(pn)) {
      v <- est[[cp]]
      pos <- seq_along(v)
      if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
        pos <- setdiff(pos, frame[["betad_fixed_idx"]])
      }
      v[pos] <- unname(par_vec[pn == cp])
      st[[cp]] <- unname(v)
    }
    st
  }
  # an end against a wall of the band it started beside
  pressed_at <- function(ck, end) {
    wl <- nl_pole_walls(frame, ck$lp, end, ck$par, ck$at, ck$dens)
    any(vapply(wl, nl_pole_pressed, NA))
  }
  best <- fit$opt
  best_ck <- NULL
  for (i in seq_along(starts)) {
    p <- fit$opt$par
    p[starts[[i]]$j] <- starts[[i]]$v
    p <- pmin(pmax(p, bounds$lower), bounds$upper)
    o <- tryCatch(suppressWarnings(optimize_obj(obj, ctl, bounds, pu,
                                                start_par = p)),
                  error = function(e) NULL)
    # the run's code is not read: the refit below reports its own, and a
    # run along a ridge to a limit of the body ends at code 7 or 8 where
    # the refit from its point ends at 0 (bernoulli seed 55 from K = -2.5)
    if (is.null(o) || !is.finite(o$objective) ||
          !(best$objective - o$objective > nl_probe_gain(best$objective)) ||
          pressed_at(checks[[i]], tpl_at(o$par))) {
      next
    }
    best <- o
    best_ck <- checks[[i]]
  }
  obj$env$.frmtmb_opt_basis <- ob_fit
  obj_settle(obj, fit$opt, fit$opt, saved)
  if (is.null(best_ck)) return(NULL)
  r <- nl_refit_held(run, tpl_at(best$par))
  if (is.null(r$fit) || pressed_at(best_ck, r$fit$estimates)) return(NULL)
  ll0 <- as.numeric(logLik_penalized(fit))
  ll1 <- suppressWarnings(as.numeric(logLik_penalized(r$fit)))
  if (!is.finite(ll1) || !(ll1 - ll0 > nl_probe_gain(ll0))) return(NULL)
  list(fit = r$fit, held = r$held, ll0 = ll0, ll1 = ll1)
}

#' The message of a default fit that nl_pole_escape() moved past the poles
#' of a body.
#'
#' @noRd
nl_pole_message <- function(pe) {
  ll <- nl_ll_pair(pe$ll0, pe$ll1)
  llw <- if (fit_penalized(pe$fit)) "penalized log-likelihood" else
    "log-likelihood"
  frm_message(
    "The fit of the nonlinear model reached ", llw, " ", ll[1L], " with a ",
    "pole of the body (a denominator at 0) among or beside the rows of the ",
    "data, where the optimizer cannot cross it. A start beyond the poles ",
    "reached ", ll[2L], ", optimizer return code ",
    format(pe$fit$opt$convergence), ", and frm() returns that fit. Give ",
    "`start` to choose your own",
    class = "frmtmb_nl_pole_start")
  invisible(NULL)
}

#' The warning on a fit from the caller's start that nl_pole_escape() beat.
#'
#' @noRd
nl_pole_warning <- function(pe, fit) {
  ll <- nl_ll_pair(pe$ll0, pe$ll1)
  llw <- if (fit_penalized(fit)) "penalized log-likelihood" else
    "log-likelihood"
  frm_warning(
    "The fit from the given start reached ", llw, " ", ll[1L], " with a ",
    "pole of the body (a denominator at 0) among or beside the rows of the ",
    "data, where the optimizer cannot cross it. A start beyond the poles ",
    "reached ", ll[2L], ", so the fit is not at the maximum. Give a `start` ",
    "on the other side of the poles",
    call. = FALSE, class = "frmtmb_nl_pole")
  invisible(NULL)
}

#' A fit from `st` by `run()`, its warnings held: `list(fit, held)`, with
#' `fit` `NULL` where it stops.
#'
#' @noRd
nl_refit_held <- function(run, st) {
  held <- new.env(parent = emptyenv())
  held$w <- list()
  r <- tryCatch(
    withCallingHandlers(suppressMessages(run(st, FALSE)),
                        warning = function(w) {
                          held$w[[length(held$w) + 1L]] <- w
                          invokeRestart("muffleWarning")
                        }),
    error = function(e) NULL)
  list(fit = r, held = held$w)
}
