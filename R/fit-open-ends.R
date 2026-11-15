# Box ends where a prior's density is 0 or infinite. frm() fits in the
# closed box first, as it always has, so every fit whose maximum is
# inside the box is the same fit. Only a run that stops on such an end
# is refitted in an open chart, where the end is at an infinite
# coordinate (backlog 5005).

#' Mark the box ends, by position in the outer vector, where a prior's
#' density is 0 or infinite (box_end_density()): `bounds$open_lower`
#' and `bounds$open_upper`, logical, the infinite ones also in
#' `bounds$inf_lower` and `bounds$inf_upper`, and `bounds$end_label_lower`
#' and `bounds$end_label_upper`, one label per position (`""` where
#' none). `open` is the resolver's
#' output, and `shim` is what resolve_bounds() reads the names against.
#'
#' Callers that build a box of their own pass no flags, and their box
#' stays closed: the trust-region rounds of `quad_settle()`, the rounds
#' of a nested quadrature fit (`quad_nest_fit()`) and `mo_search()`.
#' Each starts from a converged point, so none is a stop on such an end
#' that 0.82.0 did not have (backlog, filed by lane priors2).
#'
#' @noRd
open_end_flags <- function(bounds, open, shim) {
  np <- length(bounds$lower)
  at <- function(nm) {
    if (!length(nm) || !np) return(rep(NA_integer_, np))
    v <- stats::setNames(seq_along(nm), nm)
    pos <- tryCatch(resolve_bounds(shim, v, NULL)$lower,
                    error = keep_missing(NULL))
    k <- rep(NA_integer_, np)
    if (length(pos) != np) return(k)
    ok <- is.finite(pos)
    k[ok] <- as.integer(pos[ok])
    k
  }
  flag <- function(nm) !is.na(at(nm))
  bounds$open_lower <- flag(open$open_lower) & is.finite(bounds$lower)
  bounds$open_upper <- flag(open$open_upper) & is.finite(bounds$upper)
  bounds$inf_lower <- flag(open$inf_lower) & bounds$open_lower
  bounds$inf_upper <- flag(open$inf_upper) & bounds$open_upper
  label <- function(lab) {
    lab <- lab %||% character(0)
    k <- at(names(lab))
    ifelse(is.na(k), "", unname(lab)[k])
  }
  bounds$end_label_lower <- label(open$end_labels$lower)
  bounds$end_label_upper <- label(open$end_labels$upper)
  bounds
}

#' The positions of `par` that sit on an open box end (open_end_flags()),
#' judged as grad_bound_active() judges a bound: within
#' `sqrt(.Machine$double.eps)` of the bound's own magnitude.
#'
#' @noRd
open_end_at <- function(par, bounds, which = c("open", "inf")) {
  which <- match.arg(which)
  np <- length(par)
  lo <- rep_len(bounds$lower %||% -Inf, np)
  hi <- rep_len(bounds$upper %||% Inf, np)
  fl <- rep_len(bounds[[paste0(which, "_lower")]] %||% FALSE, np)
  fu <- rep_len(bounds[[paste0(which, "_upper")]] %||% FALSE, np)
  slack <- sqrt(.Machine$double.eps) *
    pmax(1, abs(ifelse(is.finite(lo), lo, 0)),
         abs(ifelse(is.finite(hi), hi, 0)))
  p <- as.numeric(par)
  (fl & is.finite(lo) & !(p - lo > slack)) |
    (fu & is.finite(hi) & !(hi - p > slack))
}

#' The optimizer's coordinates `cs` (opt_coords()) with each open box
#' end moved to an infinite coordinate: `q = lb + exp(u)` for an open
#' lower end, `q = ub - exp(u)` for an open upper end, and
#' `q = lb + (ub - lb) plogis(u)` where both ends are open. A closed end
#' beside an open one stays a box, at `log(ub - lb)` in `u`. The map is
#' elementwise in `q`, so `grad(g, u)` is the gradient in `q` times
#' `dq/du`. The objective is not changed, so there is no Jacobian:
#' brms's optimizing default. The start is moved inside where it is on
#' an open end. `cs$open` is TRUE when any coordinate is mapped. A
#' caller's own optimizer keeps the box it is handed.
#'
#' @noRd
open_end_coords <- function(cs, bounds, optimizer) {
  np <- length(cs$lower)
  ol <- rep_len(bounds$open_lower %||% FALSE, np)
  ou <- rep_len(bounds$open_upper %||% FALSE, np)
  if (!np || !any(ol | ou) || user_optimizer(optimizer)) return(cs)
  lo <- rep_len(cs$lower, np)
  hi <- rep_len(cs$upper, np)
  # the box is on q, the coordinate the open chart maps; a box the
  # design's coordinates moved off its parameter is not one to open
  ol <- ol & is.finite(lo)
  ou <- ou & is.finite(hi)
  if (!any(ol | ou)) return(cs)
  both <- ol & ou
  only_lo <- ol & !ou
  only_hi <- ou & !ol
  w <- hi - lo
  t_q <- function(u) {
    q <- u
    q[only_lo] <- lo[only_lo] + exp(u[only_lo])
    q[only_hi] <- hi[only_hi] - exp(u[only_hi])
    q[both] <- lo[both] + w[both] * stats::plogis(u[both])
    q
  }
  # a start on an open end, or outside it, moves a hundredth of the
  # box's width inside (a unit where the box is open on one side)
  t_u <- function(q) {
    u <- q
    span <- ifelse(is.finite(w), w, pmax(1, abs(ifelse(only_lo, lo, hi))))
    inset <- 0.01 * span
    d_lo <- q - lo
    d_hi <- hi - q
    d_lo[!(d_lo > 0)] <- inset[!(d_lo > 0)]
    d_hi[!(d_hi > 0)] <- inset[!(d_hi > 0)]
    u[only_lo] <- log(pmin(d_lo, ifelse(is.finite(w), w, Inf))[only_lo])
    u[only_hi] <- log(pmin(d_hi, ifelse(is.finite(w), w, Inf))[only_hi])
    pb <- pmin(pmax(d_lo / w, inset / w), 1 - inset / w)
    u[both] <- stats::qlogis(pb[both])
    u
  }
  dq_du <- function(u) {
    d <- rep(1, length(u))
    d[only_lo] <- exp(u[only_lo])
    d[only_hi] <- -exp(u[only_hi])
    pu <- stats::plogis(u[both])
    d[both] <- w[both] * pu * (1 - pu)
    d
  }
  to_q0 <- cs$to_q
  to_p0 <- cs$to_p
  grad0 <- cs$grad
  lower <- lo
  upper <- hi
  lower[ol | ou] <- -Inf
  upper[ol | ou] <- Inf
  # the closed end beside an open one
  cl <- (only_lo & is.finite(hi)) | (only_hi & is.finite(lo))
  upper[cl] <- log(w[cl])
  units <- cs$units
  if (!is.null(units)) {
    units <- rep_len(units, np)
    units[ol | ou] <- 1
  }
  out <- cs
  out$open <- TRUE
  out$to_q <- function(p) t_u(to_q0(p))
  out$to_p <- function(u) to_p0(t_q(u))
  out$grad <- function(g, q) grad0(g) * dq_du(q)
  out$lower <- lower
  out$upper <- upper
  out$units <- units
  out
}

#' nlminb's "NA/NaN function evaluation" from the box run, held where a
#' box end is open (open_end_flags()). The box run steps onto such an
#' end on its way, and when open_end_rescue() keeps its refit at the
#' mode, with code 0, the warning is about a point the fit left: A, B,
#' the pushed fits and case H all raised it (short review of lane
#' priors2). `$catch` is the handler, and `$release(opt)` raises the
#' held warnings again where the box run's own result is returned.
#'
#' @noRd
open_end_hold <- function(bounds, np, optimizer) {
  env <- new.env(parent = emptyenv())
  env$w <- list()
  on <- np > 0L && !user_optimizer(optimizer) &&
    any(rep_len(bounds$open_lower %||% FALSE, np) |
          rep_len(bounds$open_upper %||% FALSE, np))
  txt <- "NA/NaN function evaluation"
  txt <- unique(c(txt, gettext(txt, domain = "stats"),
                  gettext(txt, domain = "R-stats")))
  env$catch <- function(w) {
    if (on && conditionMessage(w) %in% txt) {
      env$w[[length(env$w) + 1L]] <- w
      invokeRestart("muffleWarning")
    }
  }
  env$release <- function(opt) {
    if (!isTRUE(opt$open_refit)) for (w in env$w) frm_warning(w)
    invisible(NULL)
  }
  env
}

#' A box run that failed: with an open box end (open_end_flags()) a
#' stand-in stop at the start, which open_end_rescue() refits in the
#' open chart, carrying the error it raises if that refit fails too;
#' without one, the error itself. On 0.82.0 a `gamma(0.5, 1)` with
#' `lb = 0` on a poisson("identity") slope failed with "NA/NaN gradient
#' evaluation" where its maximum is at 1.3833 (review of lane priors2,
#' F1).
#'
#' @noRd
open_end_failed <- function(par, e, bounds, optimizer) {
  np <- length(par)
  ends <- any(rep_len(bounds$open_lower %||% FALSE, np) |
                rep_len(bounds$open_upper %||% FALSE, np))
  if (!ends || user_optimizer(optimizer)) frm_stop(e)
  list(par = par, objective = NaN, convergence = 1L,
       message = conditionMessage(e), rescue_error = e)
}

#' A fit that stopped on an open box end (open_end_at()), or at a
#' non-finite objective beside one, refitted in the open chart
#' (open_end_coords()) and then once more in the box `cs` from there.
#'
#' The box run comes first, so every fit whose maximum is inside the box
#' is the fit 0.82.0 made, to the bit (dev/priors2-survey.R, the review's
#' 15 fits). Only a stop on the end is moved: under `beta(2, 2)` with
#' `lb = 0` a slope the data put below 0 stopped at 0 with the objective
#' Inf, and the open chart reaches the mode. The refit is kept where it
#' is finite and no higher, or where the box stop was not finite
#' (`opt$open_refit`). A fit that ends beside an infinite-density end
#' goes to open_end_probe(), which alone decides a refusal.
#'
#' @noRd
open_end_rescue <- function(obj, opt, optimizer, fn, gr, cs, bounds,
                            control, tally, verbose = 0L) {
  hit <- open_end_at(opt$par, bounds)
  ends <- any(rep_len(bounds$open_lower %||% FALSE, length(opt$par)) |
                rep_len(bounds$open_upper %||% FALSE, length(opt$par)))
  # a fit that ends inside is checked at its infinite-density ends; that
  # evaluation is not the optimizer's, so it stays out of the count
  # se_check_at_fit() weighs the Hessian against
  fn_probe <- function(p) obj$fn(p)
  if (!any(hit) && !(ends && !isTRUE(is.finite(opt$objective)))) {
    if (user_optimizer(optimizer)) return(opt)
    return(open_end_probe(obj, fn_probe, opt, bounds, opt$objective))
  }
  co <- open_end_coords(cs, bounds, optimizer)
  if (!isTRUE(co$open)) {
    if (!is.null(opt$rescue_error)) frm_stop(opt$rescue_error)
    return(opt)
  }
  saved <- obj_state_save(obj)
  fn_u <- function(u) fn(co$to_p(u))
  gr_u <- function(u) co$grad(gr(co$to_p(u)), u)
  o <- tryCatch(run_optimizer(optimizer, co$to_q(opt$par), fn_u, gr_u,
                              co$lower, co$upper, control$optCtrl,
                              co$units, tally),
                error = keep_missing(NULL))
  if (!is.null(o)) {
    o$par <- stats::setNames(as.numeric(co$to_p(o$par)), names(opt$par))
    if (all(is.finite(o$par)) && isTRUE(is.finite(o$objective))) {
      o <- open_end_finish(obj, o, optimizer, fn, gr, cs, control, tally,
                           verbose)
    }
  }
  # -Inf is the refit running to an infinite-density end, which the
  # refusal below reads
  keep <- !is.null(o) && all(is.finite(o$par)) &&
    (isTRUE(is.finite(o$objective)) || isTRUE(o$objective == -Inf)) &&
    (!isTRUE(is.finite(opt$objective)) ||
       isTRUE(o$objective <= opt$objective))
  if (verbose) {
    vb_say("the box run stopped on a bound where the prior's density is ",
           "0 or infinite; a run where that bound is out of reach ",
           if (keep) "is kept" else "is dropped")
  }
  if (!keep) {
    if (!is.null(opt$rescue_error)) frm_stop(opt$rescue_error)
    obj_settle(obj, opt, opt, saved)
    o <- opt
  }
  o$open_refit <- keep
  open_end_probe(obj, fn, o, bounds, c(opt$objective, o$objective))
}

#' Whether the negative log posterior falls without bound toward an
#' infinite-density end that the fit `o` stopped beside (open_end_at()).
#'
#' Position alone proves nothing: poisson("identity") with `gamma(0.5,
#' 1)`, `lb = 0` and an exposure of 1e7 per row has its maximum at b =
#' 4.854e-10, inside the slack of the end, and a first version refused
#' it (short review of lane priors2, case H). So the objective is
#' evaluated along that coordinate at points that approach the end, two
#' decades at a time, from the slack down to the nearest double. The
#' fit is refused (`o$open_unbounded`, read by check_map_unbounded())
#' only where it falls at every step and ends below the best objective
#' found (`seen`) by a margin; the density's own pole then outruns any
#' likelihood that is finite there. Otherwise the best point is kept,
#' and where that is a probe point, nearer the end than the optimizer
#' stopped, it is returned with code 1 and a warning
#' (`o$open_end_stop`).
#'
#' A fit that ends away from an infinite-density end is checked at the
#' nearest double past it; where the objective there is `-Inf` or below
#' the fit's by the same margin, the fit is kept and marked a local mode
#' (`o$open_local`, a warning). The state of `obj` is left as the fit
#' left it where no end is near.
#'
#' @noRd
open_end_probe <- function(obj, fn, o, bounds, seen) {
  np <- length(o$par)
  lo <- rep_len(bounds$lower %||% -Inf, np)
  hi <- rep_len(bounds$upper %||% Inf, np)
  il <- rep_len(bounds$inf_lower %||% FALSE, np)
  iu <- rep_len(bounds$inf_upper %||% FALSE, np)
  if (!any(il | iu)) return(o)
  slack <- sqrt(.Machine$double.eps) *
    pmax(1, abs(ifelse(is.finite(lo), lo, 0)),
         abs(ifelse(is.finite(hi), hi, 0)))
  p0 <- as.numeric(o$par)
  near_lo <- il & is.finite(lo) & !(p0 - lo > slack)
  near_hi <- iu & is.finite(hi) & !(hi - p0 > slack)
  far_lo <- il & is.finite(lo) & !near_lo
  far_hi <- iu & is.finite(hi) & !near_hi
  if (!any(near_lo | near_hi | far_lo | far_hi)) return(o)
  labs <- function(side) {
    lab <- rep_len(bounds[[paste0("end_label_", side)]] %||% "", np)
    ifelse(nzchar(lab), lab, names(o$par) %||% "")
  }
  lab_lo <- labs("lower")
  lab_hi <- labs("upper")
  saved <- obj_state_save(obj)
  best <- o
  moved <- FALSE
  unb <- character(0)
  stops <- character(0)
  f_at <- function(p) {
    v <- tryCatch(suppressWarnings(as.numeric(fn(p))),
                  error = function(e) NaN)
    if (length(v) == 1L) v else NaN
  }
  for (i in which(near_lo | near_hi)) {
    side_lo <- near_lo[i] &&
      (!near_hi[i] || p0[i] - lo[i] <= hi[i] - p0[i])
    e <- if (side_lo) lo[i] else hi[i]
    sgn <- if (side_lo) 1 else -1
    q <- open_end_steps(e, sgn, slack[i])
    p <- as.numeric(best$par)
    f <- vapply(q, function(qk) {
      p[i] <- qk
      f_at(p)
    }, 0)
    ok <- !anyNA(f) && !any(f == Inf)
    fin <- is.finite(f)
    # once -Inf, the rest must be -Inf too: the pole reached
    tail_ok <- ok && (all(fin) || all(!fin[seq(which(!fin)[1L], length(f))]))
    falls <- tail_ok && all(diff(f[fin]) < 0)
    ref <- c(seen, best$objective, f[1L])
    ref <- if (any(is.finite(ref))) min(ref[is.finite(ref)]) else NaN
    margin <- max(1, sqrt(.Machine$double.eps) * abs(ref))
    if (falls && is.finite(ref) && f[length(f)] < ref - margin) {
      unb <- c(unb, if (side_lo) lab_lo[i] else lab_hi[i])
      next
    }
    j <- if (any(!is.na(f))) which.min(f) else NA_integer_
    if (!is.na(j) && (!isTRUE(is.finite(best$objective)) ||
                        f[j] < best$objective)) {
      best$par[i] <- q[j]
      best$objective <- f[j]
      moved <- TRUE
      stops <- c(stops, if (side_lo) lab_lo[i] else lab_hi[i])
    }
  }
  if (moved) {
    best$convergence <- 1L
    best$message <- "stopped at the end of the prior's support"
    best$open_end_stop <- stops
  }
  # a fit inside, away from an infinite-density end: below it at the
  # nearest double past the end, the posterior is unbounded there and
  # the fit a local mode. y ~ x + z, gaussian, gamma(0.5, 1) with lb = 0
  # on x ended at code 0 with nothing said, where 0.82.0's box run went
  # to the pole at -Inf (short review of lane priors2, P3)
  local <- character(0)
  f_fit <- best$objective
  if (isTRUE(is.finite(f_fit))) {
    margin <- max(1, sqrt(.Machine$double.eps) * abs(f_fit))
    for (k in 1:2) {
      far <- if (k == 1L) far_lo else far_hi
      for (i in which(far)) {
        e <- if (k == 1L) lo[i] else hi[i]
        p <- as.numeric(best$par)
        p[i] <- open_end_nearest(e, if (k == 1L) 1 else -1)
        f <- f_at(p)
        if (isTRUE(f == -Inf) || isTRUE(f < f_fit - margin)) {
          local <- c(local, if (k == 1L) lab_lo[i] else lab_hi[i])
        }
      }
    }
  }
  if (!any(near_lo | near_hi)) {
    # the state the fit left, untouched, as where no end is checked
    obj_state_restore(obj, saved)
    if (length(local)) o$open_local <- local
    return(o)
  }
  obj_settle(obj, best, best, saved)
  if (length(unb)) best$open_unbounded <- unb
  if (length(local)) best$open_local <- local
  best
}

#' The nearest double past the end `e` in the direction `sgn`.
#'
#' @noRd
open_end_nearest <- function(e, sgn) {
  u <- if (e == 0) 2^-1074 else 2^(floor(log2(abs(e))) - 53)
  while (e + sgn * u == e) u <- 2 * u
  e + sgn * u
}

#' Points from `e + sgn * slack` toward the end `e`, two decades apart,
#' then the nearest double past `e`: the smallest distance a parameter
#' can be from that end.
#'
#' @noRd
open_end_steps <- function(e, sgn, slack) {
  d <- slack * 100^-(0:400)
  d <- d[d > 0]
  q <- e + sgn * d
  q <- q[q != e]
  unique(c(q, open_end_nearest(e, sgn)))
}

#' A last run in the box's own coordinates `cs` (opt_coords()) from the
#' point the open chart reached, kept where it converges no higher, so
#' that the reported point is one of the box's stops. At an end where
#' the density is zero the objective is infinite, so the run does not
#' step onto it.
#'
#' @noRd
open_end_finish <- function(obj, opt, optimizer, fn, gr, cs, control,
                            tally, verbose = 0L) {
  saved <- obj_state_save(obj)
  mapped <- !is.null(cs$basis)
  fn_q <- if (!mapped) fn else function(q) fn(cs$to_p(q))
  gr_q <- if (!mapped) gr else function(q) cs$grad(gr(cs$to_p(q)))
  o <- tryCatch(run_optimizer(optimizer, cs$to_q(opt$par), fn_q, gr_q,
                              cs$lower, cs$upper, control$optCtrl,
                              cs$units, tally),
                error = keep_missing(NULL))
  if (!is.null(o)) o$par <- cs$to_p(o$par)
  keep <- !is.null(o) && all(is.finite(o$par)) &&
    isTRUE(is.finite(o$objective)) && isTRUE(o$objective <= opt$objective) &&
    (isTRUE(o$convergence == 0) || !isTRUE(opt$convergence == 0))
  if (verbose) {
    vb_say("the run in the open chart of a bound stopped at code ",
           opt$convergence, "; a run in the box from there ",
           if (keep) "is kept" else "is dropped")
  }
  if (!keep) {
    obj_settle(obj, opt, opt, saved)
    return(opt)
  }
  o$par <- stats::setNames(as.numeric(o$par), names(opt$par))
  o
}

#' Refuse a MAP fit where the log posterior was shown to rise without
#' bound toward a box end at which the prior's density is infinite
#' (open_end_probe()): there is no maximum to report. The refusal
#' follows the fit, because whether the likelihood stays finite at the
#' end is a property of the data and the family (review of lane
#' priors2, F1). frmtmb.sample's `frm_sample()` samples such a
#' posterior with the bound written class-wide, as brms takes it. A fit
#' the probe moved to the end without that proof is kept, with a
#' warning.
#'
#' @noRd
check_map_unbounded <- function(opt, bounds) {
  l <- opt$open_local
  if (length(l)) {
    frm_warning("The posterior is unbounded toward the end of a prior's ",
                "support where its density is infinite: ",
                paste(l, collapse = "; "), ". The log posterior at the ",
                "nearest value past that end, the other parameters held ",
                "at the fit, is above the fit's, so the point reported is ",
                "a local mode and no maximum exists. Move the bound inside ",
                "the support, or use a density that is finite there",
                call. = FALSE, class = "frmtmb_prior_unbounded_local")
  }
  s <- opt$open_end_stop
  if (length(s)) {
    frm_warning("The fit stopped at the end of the support of a prior ",
                "whose density is infinite there: ",
                paste(s, collapse = "; "), ". The log posterior was not ",
                "shown to rise without bound toward that end, so the ",
                "point reported is the best one found, not a proven ",
                "maximum (code 1)", call. = FALSE)
  }
  u <- opt$open_unbounded
  if (!length(u)) return(invisible(opt))
  frm_stop("A prior whose density is infinite at a bound has no maximum ",
           "here: ", paste(u, collapse = "; "), ". Toward that end of ",
           "the density's support the log posterior rises at every ",
           "point, down to the nearest double, and passes the best fit ",
           "found, so a MAP fit has no estimate to report. Move ",
           "the bound inside the support, use a density that is finite ",
           "there, or sample the posterior with frmtmb.sample's ",
           "frm_sample(), which takes the bound on the class row, not on ",
           "one coefficient", call. = FALSE)
}
