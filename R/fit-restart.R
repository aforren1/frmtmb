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
#' nearest wall, the median gap between adjacent walls, and the gaps that
#' nl_pole_holes() reads as holes in the data; `NULL` where no
#' denominator reads `c` affinely. A denominator that is not affine in `c`
#' is not read. The two such that the bodies of the fresh rows have,
#' `EC50^h + x^h` with `x >= 0` and `1 + exp(b * x)`, are positive for
#' every real coefficient at which the body is defined, so they have no
#' wall to read.
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
                                      any(a[reads] < 0),
                                    holes = nl_pole_holes(w))
  }
  if (length(out)) out
}

#' The sizes of nl_pole_holes()'s rule: the half width of the window of
#' gaps that sets a gap's local scale, the walls a hole needs on each
#' side, the odds against a hole on rows of a smooth density, and the
#' most holes one band starts in.
#'
#' Measured by dev/nlband-holes2.R on 400 draws of each design of 20 to
#' 1000 rows (100 and 20 draws at 1e4 and 1e5): with no hole in the data,
#' a start ran on at most 3.2 percent of uniform, normal, lognormal,
#' exponential and t3 draws up to 1000 rows (lognormal at 1000), and on at
#' most 10 percent at 1e4 and 1e5 rows (exponential at 1e5, 2 of 20);
#' never more than 2 starts on one draw. A hole of half width 0.3 in 400
#' uniform rows on (0, 4), the band data of the review of lane nlstop2
#' (F4), was found on 400 of 400 draws, one of half width 0.1 on 398. A
#' window of 21 gaps with 5 walls a side started on up to 75 percent of
#' the draws without a hole: the median of few gaps is noisy, and the gaps
#' of a tail grow.
#'
#' @noRd
nl_hole_window <- 20L
nl_hole_side <- 10L
nl_hole_odds <- 1000
nl_hole_max <- 2L

#' The gaps between adjacent walls `w` (sorted) that are holes in the
#' data: wider than the rows around them make likely, with rows on both
#' sides.
#'
#' A pole the optimizer cannot reach can sit in a gap of the data: `a /
#' (x + b)` on y = 1 / (x - 2) + noise, with no rows within 0.3 of 2,
#' stopped at b 0.39 with every denominator positive, 592 below the
#' maximum at b -2.00 (review of lane nlstop2, F4). Each of the n - 1 gaps
#' between rows' walls is a piece of the line with a maximum of its own,
#' so a start in each is out of reach on dense data. The maximum can be in
#' a narrow gap: y = 1 / (x - 2) + noise on 400 rows of U(0, 4), with no
#' hole, has its pole between two adjacent rows, and the fit stops 1575 to
#' 2100 short in silence on 0.83.0 and here (review of lane nlband, item
#' 3). This rule does not reach that; it starts only where a gap is wide
#' enough to be a hole in the data. A gap is a
#' hole when it is wider than log2(nl_hole_odds * m) times the running
#' median of the 2 * nl_hole_window + 1 gaps around it, m gaps in all:
#' on rows of a smooth density a gap is wider than t local means with
#' probability exp(-t), so m gaps hold one wider than log(1000 m) means,
#' log2(1000 m) medians, with odds of about 1 in 1000 if the local median
#' were exact (it is a median of 41 gaps, so the rate measured at
#' nl_hole_window is higher). A hole needs `nl_hole_side` walls or more on
#' each side. Returns a matrix with columns `lo` and `hi`, the widest
#' `nl_hole_max` holes first by that ratio, or `NULL`.
#'
#' @noRd
nl_pole_holes <- function(w) {
  w <- w[c(TRUE, diff(w) > 1e-9 * (w[length(w)] - w[1L]))]
  d <- diff(w)
  m <- length(d)
  if (m <= 2L * nl_hole_side) return(NULL)
  win <- 2L * nl_hole_window + 1L
  loc <- if (m <= win) {
    rep(stats::median(d), m)
  } else {
    stats::runmed(d, win, endrule = "constant")
  }
  q <- d / loc
  i <- seq_len(m)
  hit <- which(i > nl_hole_side & i <= m - nl_hole_side &
                 q > log2(nl_hole_odds * m))
  if (!length(hit)) return(NULL)
  hit <- hit[order(-q[hit], hit)][seq_len(min(length(hit), nl_hole_max))]
  cbind(lo = w[hit], hi = w[hit + 1L])
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
#' clear of every band, on data with no hole, pays only the evaluations of
#' its denominators.
#'
#' A run is kept if it ends clear of the band it started beside and higher
#' than the fit by more than nl_probe_gain(); then frm() fits again from
#' its point by `run()`, as the sign probe does, and that fit must end
#' clear of the band too. A run that ends against a wall is not kept:
#' bernoulli seed 55 rises toward the wall of its largest xp, where the
#' likelihood has a supremum and no maximum.
#'
#' Whether or not a band holds the fit, the optimizer also runs from the
#' middle of each hole in the data that the fit is not in
#' (nl_pole_holes()): a fit clear of every wall cannot reach a maximum
#' with a pole between rows. It runs twice, the rest of the model at the
#' fit and at the default start (nl_cold_par(); once where the two are
#' the same). Such a run, and its refit, must end no nearer a wall than
#' the walls' median gap, a proxy for a supremum on a wall that also
#' refuses a true maximum with its pole that near a row (a pole 0.003
#' below the first row above the hole, review of lane nlband, item 4).
#' The refit must have a standard error of the hole's coefficient below
#' the span of the walls (nl_hole_se()), which refuses a run that rises
#' without end. Every run that beats the fit is a candidate, the highest
#' first; a candidate whose refit is refused passes to the next.
#'
#' Not run under `frmtmb_control(nl_starts = FALSE)`, which turns off this
#' restart and the sign probe.
#'
#' Returns `NULL`, or a list as nl_sign_probe() returns it, with `where`
#' the start that won: its `kind` ("past" the band or in a "hole"), the
#' coefficient's `name`, and for a hole its walls `lo` and `hi`; and
#' `refused`, the first hole refit that the standard error refused, for
#' nl_pole_refused_message(). Where every candidate is refused and one by
#' its standard error, `fit` is `NULL` and `refused` is set.
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
  cold <- NULL
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
          ck <- list(lp = lp, par = par, at = at, dens = dens, kind = "past")
          if (nl_pole_pressed(wl)) {
            half <- (wl$hi - wl$lo) / 2
            to <- c(if (c0 >= wl$lo) wl$lo - half,
                    if (c0 <= wl$hi) wl$hi + half)
            for (v in to) {
              starts[[length(starts) + 1L]] <- list(j = pos[at], v = v)
              checks[[length(checks) + 1L]] <- ck
            }
          }
          hs <- wl$holes
          # the rest of the model at the fit, and at the default start:
          # from a fit on a ridge to a limit of the body (Vm 3e8 in Vm * x
          # / (K + x)) the run from the hole ended at log-likelihood
          # -5.5e9, where a start at Vm 0 to 100 reached the hole's
          # maximum (dev/nlband-dbg2.R), and from a student fit with log
          # nu at its limit, 19.7, the kept fit was 16 to 44 below the
          # hole's maximum (review of lane nlband, item 2)
          if (NROW(hs) && is.null(cold)) cold <- nl_cold_par(frame, pn)
          from <- list(NULL)
          q <- cold
          if (length(q)) q[pos[at]] <- fit$opt$par[pos[at]]
          if (length(q) && !identical(unname(q), unname(fit$opt$par))) {
            from <- c(from, list(q))
          }
          for (h in seq_len(NROW(hs))) {
            lo <- unname(hs[h, "lo"])
            hi <- unname(hs[h, "hi"])
            if (c0 > lo && c0 < hi) next
            ck$kind <- "hole"
            ck$j <- pos[at]
            ck$lo <- lo
            ck$hi <- hi
            ck$span <- wl$hi - wl$lo
            for (q in from) {
              starts[[length(starts) + 1L]] <- list(j = pos[at],
                                                    v = (lo + hi) / 2,
                                                    from = q)
              checks[[length(checks) + 1L]] <- ck
            }
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
  # an end against a wall of the band it started beside; from a hole, the
  # denominator takes both signs by design, so only the distance counts
  pressed_at <- function(ck, end) {
    wl <- nl_pole_walls(frame, ck$lp, end, ck$par, ck$at, ck$dens)
    held <- if (identical(ck$kind, "hole")) {
      function(x) x$dist < x$gap
    } else {
      nl_pole_pressed
    }
    any(vapply(wl, held, NA))
  }
  # every run that beats the fit, best first: a refit that is refused
  # passes the turn to the next (short check of lane nlband, item 2)
  runs <- list()
  for (i in seq_along(starts)) {
    p <- starts[[i]]$from %||% fit$opt$par
    p[starts[[i]]$j] <- starts[[i]]$v
    p <- pmin(pmax(p, bounds$lower), bounds$upper)
    o <- tryCatch(suppressWarnings(optimize_obj(obj, ctl, bounds, pu,
                                                start_par = p)),
                  error = function(e) NULL)
    # the run's code is not read: the refit below reports its own, and a
    # run along a ridge to a limit of the body ends at code 7 or 8 where
    # the refit from its point ends at 0 (bernoulli seed 55 from K = -2.5)
    if (is.null(o) || !is.finite(o$objective) ||
          !(fit$opt$objective - o$objective >
              nl_probe_gain(fit$opt$objective)) ||
          pressed_at(checks[[i]], tpl_at(o$par))) {
      next
    }
    runs[[length(runs) + 1L]] <- list(o = o, ck = checks[[i]])
  }
  obj$env$.frmtmb_opt_basis <- ob_fit
  obj_settle(obj, fit$opt, fit$opt, saved)
  if (!length(runs)) return(NULL)
  runs <- runs[order(vapply(runs, function(x) x$o$objective, 0))]
  ll0 <- as.numeric(logLik_penalized(fit))
  name_of <- function(ck) {
    par_template_names(frame[["par_template"]][[ck$par]], ck$par)[ck$at]
  }
  refused <- NULL
  for (x in runs) {
    ck <- x$ck
    r <- nl_refit_held(run, tpl_at(x$o$par))
    if (is.null(r$fit) || pressed_at(ck, r$fit$estimates)) next
    ll1 <- suppressWarnings(as.numeric(logLik_penalized(r$fit)))
    if (!is.finite(ll1) || !(ll1 - ll0 > nl_probe_gain(ll0))) next
    if (identical(ck$kind, "hole")) {
      s <- nl_hole_se(r$fit, ck)
      if (!isTRUE(s < ck$span)) {
        if (is.null(refused)) {
          refused <- list(ll0 = ll0, ll1 = ll1, se = s, span = ck$span,
                          name = name_of(ck), lo = ck$lo, hi = ck$hi)
        }
        next
      }
    }
    return(list(fit = r$fit, held = r$held, ll0 = ll0, ll1 = ll1,
                where = list(kind = ck$kind, name = name_of(ck),
                             lo = ck$lo, hi = ck$hi),
                refused = refused))
  }
  if (is.null(refused)) NULL else list(fit = NULL, refused = refused)
}

#' The outer parameter vector, named as `pn`, at the cold start
#' make_start() builds with no `start` and no prior: the start a default
#' fit takes, with each intercept at the family's initializer. `NULL`
#' where it cannot be built or does not match `pn`.
#'
#' @noRd
nl_cold_par <- function(frame, pn) {
  st <- tryCatch(
    suppressWarnings(suppressMessages(make_start(frame, NULL))),
    error = function(e) NULL)
  if (is.null(st)) return(NULL)
  q <- rep(NA_real_, length(pn))
  for (cp in unique(pn)) {
    v <- st[[cp]]
    pos <- seq_along(v)
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      pos <- setdiff(pos, frame[["betad_fixed_idx"]])
    }
    if (length(pos) != sum(pn == cp)) return(NULL)
    q[pn == cp] <- unname(v[pos])
  }
  if (anyNA(q)) return(NULL)
  stats::setNames(q, pn)
}

#' The standard error of the hole's coefficient in a hole fit, from the
#' standard-error check's own analysis of the fit (se_fit_analysis()), or
#' NaN. nl_pole_escape() keeps the fit only where it is below the span of
#' the walls, so that the data place the pole among the rows' poles.
#'
#' A run from a hole can rise without end. On bernoulli data that the
#' hole separates (y = 1 for x > 2 on the band data), the run and its
#' refit ended at log-likelihood -5e-268 with a = 1221, code 0 and no
#' warning, and the escape returned it (review of lane nlband, item 1).
#' Neither the limit check nor the standard-error check marks that fit:
#' nl_limit_message() returns NULL with and without its gate, and the
#' analysis keeps every standard error, 8.6e133 for a and 1.4e131 for b
#' (dev/nlband-dbg3.R). That is a supremum, not a maximum. The bound is
#' the span of the walls, about 4 on the band data, not half the hole's
#' width: at half the width the rule refused 5 of 101 maxima in a sweep
#' of noisy data, each 0.056 to 1.462 above the fit returned, with
#' standard errors up to 1.4 half widths (short check of lane nlband,
#' item 2); the separated fits sit 1e130 spans above it.
#'
#' @noRd
nl_hole_se <- function(fit, ck) {
  an <- tryCatch(se_fit_analysis(fit, force = TRUE),
                 error = function(e) NULL)
  V <- an$shown
  if (!is.matrix(V) || NROW(V) < ck$j) return(NaN)
  sqrt(V[ck$j, ck$j])
}

#' The sentence of nl_pole_message() and nl_pole_warning() that says
#' where the start that won began.
#'
#' @noRd
nl_pole_from <- function(pe) {
  w <- pe$where
  if (!identical(w$kind, "hole")) return("A start beyond the poles")
  paste0("A start with ", w$name, " between ", signif(w$lo, 4), " and ",
         signif(w$hi, 4), ", a gap that the rows' poles leave,")
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
    "The fit of the nonlinear model reached ", llw, " ", ll[1L],
    nl_pole_where(pe), " ", nl_pole_from(pe), " reached ", ll[2L],
    ", optimizer return code ", format(pe$fit$opt$convergence),
    ", and frm() returns that fit. Give `start` to choose your own",
    class = "frmtmb_nl_pole_start")
  invisible(NULL)
}

#' Where the fit that nl_pole_escape() beat sits against the poles.
#'
#' @noRd
nl_pole_where <- function(pe) {
  if (identical(pe$where$kind, "hole")) {
    return(paste0(". The body has a pole (a denominator at 0) at one ",
                  "value of ", pe$where$name, " for each row of the data, ",
                  "and the optimizer cannot cross one."))
  }
  paste0(" with a pole of the body (a denominator at 0) among or beside ",
         "the rows of the data, where the optimizer cannot cross it.")
}

#' The warning on a fit from the caller's start that nl_pole_escape() beat.
#'
#' @noRd
nl_pole_warning <- function(pe, fit) {
  ll <- nl_ll_pair(pe$ll0, pe$ll1)
  llw <- if (fit_penalized(fit)) "penalized log-likelihood" else
    "log-likelihood"
  frm_warning(
    "The fit from the given start reached ", llw, " ", ll[1L],
    nl_pole_where(pe), " ", nl_pole_from(pe), " reached ", ll[2L],
    ", so the fit is not at the maximum. Give a `start` ",
    if (identical(pe$where$kind, "hole")) "in that gap" else
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
  # a boundary message describes the refit itself, which may be the fit
  # returned: student() on the band data, seed 2, kept nu at its limit
  # with no word of it (short check of lane nlband, item 1)
  r <- tryCatch(
    withCallingHandlers(run(st, FALSE),
                        warning = function(w) {
                          held$w[[length(held$w) + 1L]] <- w
                          invokeRestart("muffleWarning")
                        },
                        message = function(m) {
                          if (inherits(m, "frmtmb_boundary_fit")) {
                            held$w[[length(held$w) + 1L]] <- m
                          }
                          invokeRestart("muffleMessage")
                        }),
    error = function(e) NULL)
  list(fit = r, held = held$w)
}

#' The message that nl_pole_escape() refused a higher hole fit, because
#' the data do not place its pole: the standard error of the coefficient
#' is not below the span of the walls.
#'
#' @noRd
nl_pole_refused_message <- function(pe) {
  rf <- pe$refused
  if (is.null(rf)) return(invisible(NULL))
  ll <- nl_ll_pair(rf$ll0, rf$ll1)
  frm_message(
    "A start with ", rf$name, " between ", signif(rf$lo, 4), " and ",
    signif(rf$hi, 4), ", a gap that the rows' poles leave, reached ",
    "log-likelihood ", ll[2L], " against ", ll[1L], ", but frm() does not ",
    "return that fit: the standard error of ", rf$name, " there, ",
    signif(rf$se, 3), ", is not below the span of the poles, ",
    signif(rf$span, 3), ", so the data do not place the pole (on data ",
    "the gap separates, the likelihood rises without end). Give `start` ",
    "in that gap to see the fit",
    class = "frmtmb_nl_pole_refused")
  invisible(NULL)
}
