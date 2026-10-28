# frmtmb_frame -> negative log-likelihood closure for RTMB::MakeADFun.
#
# Discipline (see SPEC.md section 2): everything observation-length is
# vectorized; loops run only over model structure (blocks, linear
# predictors, responses), which is resolved at tape time; nothing branches
# on parameter values. Data referenced by the closure is baked into the
# tape as constants.

#' Both tails of a family's CDF on the log scale at `q` (data), one
#' value per row: `lF = log F(q)` and `lS = log S(q)`, from the family's
#' `cdf_log` slot, with what else the pair carries for `log_window()`.
#'
#' @noRd
fam_log_pair <- function(fam, q, dpv, av, extra, n) {
  full <- function(v) if (is.null(v) || length(v) == n) v else v + numeric(n)
  pf <- fam[["cdf_log"]]
  lp <- if (length(formals(pf)) >= 4L) {
    pf(q, dpv, av, extra)
  } else {
    pf(q, dpv, av)
  }
  # a point-mass family's base part, or a window function of the
  # family's own and the per-row values it reads
  out <- lp
  for (nm in setdiff(names(lp), c("base", "in0", "lmass"))) {
    out[[nm]] <- full(lp[[nm]])
  }
  if (!is.null(lp[["base"]])) {
    out$base <- list(lF = full(lp[["base"]][["lF"]]),
                     lS = full(lp[["base"]][["lS"]]))
  }
  if (!is.null(lp[["in0"]])) out$in0 <- rep_len(lp[["in0"]], n)
  out
}

#' `log P(lo < Y <= hi)` on the rows `rows`, from the log pairs `pl` at
#' the lower edge and `ph` at the upper edge. `lo_abs` and `hi_abs`
#' (data, one per row) mark an edge that is absent, where the window
#' is open on that side.
#'
#' With both edges present the mass is a difference of two tails, and
#' it is taken in the tail the lower edge leaves more of: as
#' `S(lo) - S(hi)` where `S(lo) >= 1/2` and as `F(hi) - F(lo)` below
#' that, each on the log scale. A log near zero keeps every digit of
#' its complement, so either form keeps the digits the other loses: the
#' probability-scale `F(hi) - F(lo)` had both terms round toward one
#' above the median, and a hurdle's `(0, y2]` was `(hu + tiny) - hu`
#' (dev/reviews/2026-10-20-cdfs.md, finding 3).
#'
#' A zero-inflated or hurdle pair also carries its base family's pair
#' (`base`), the log weight of the base part (`lw`) and the rows whose
#' `q` is at or above zero (`in0`, data). A window whose lower edge is
#' at or above zero holds none of the point mass, and its mass is the
#' base's window times the weight. Taken from the mixture's own tails
#' it would be a difference of two numbers that share the point mass:
#' on a hurdle `(0, y2]` the log survivors are `l1m` and
#' `l1m + lS_b(y2)`, and `lS_b(y2)` came back from their difference to
#' the doubles' absolute accuracy, not its own.
#'
#' @noRd
log_window <- function(pl, ph, lo_abs, hi_abs, rows) {
  "[<-" <- RTMB::ADoverload("[<-")
  la <- lo_abs[rows]
  ha <- hi_abs[rows]
  # a window open on the same side on every row is one tail, read with
  # no sub-assignment: every row of a trunc(lb = ) normalizer is
  pick <- function(v) if (length(rows) == length(v)) v else v[rows]
  if (all(!la & ha)) return(pick(pl[["lS"]]))
  if (all(la & !ha)) return(pick(ph[["lF"]]))
  out <- numeric(length(rows))
  i <- which(!la & !ha)
  if (length(i)) {
    r <- rows[i]
    out[i] <- log_window_two(pl, ph, r)
    if (!is.null(pl[["lw"]]) && !is.null(ph[["lw"]])) {
      j <- which(pl[["in0"]][r])
      if (length(j)) {
        rj <- r[j]
        out[i[j]] <- pl[["lw"]][rj] +
          log_window_two(pl[["base"]], ph[["base"]], rj)
      }
    }
  }
  i <- which(la & !ha)
  if (length(i)) out[i] <- ph[["lF"]][rows[i]]
  i <- which(!la & ha)
  if (length(i)) out[i] <- pl[["lS"]][rows[i]]
  out
}

#' `log P(lo < Y <= hi)` on rows `r` with both edges present, in the
#' tail the lower edge leaves more of (`log_window()`).
#'
#' A family whose window is not a difference of its two tails anywhere
#' supplies `lmass` with its pair, which this calls instead: a von Mises
#' arc on the far side of the mean has a mass far below either tail at
#' its ends, which both run on to the other side of the circle, so
#' neither difference holds a digit of it (`vm_lmass()`).
#'
#' @noRd
log_window_two <- function(pl, ph, r) {
  if (is.function(pl[["lmass"]])) return(pl[["lmass"]](pl, ph, r))
  lSl <- pl[["lS"]][r]
  lFh <- ph[["lF"]][r]
  w <- 0.5 * (1 + sign_pos(lSl - log(0.5)))
  viaS <- lSl + log1m_exp_ad(ph[["lS"]][r] - lSl)
  viaF <- lFh + log1m_exp_ad(pl[["lF"]][r] - lFh)
  w * viaS + (1 - w) * viaF
}

#' The edges of the `trunc()` window, one per row (data), with the rows
#' where an edge is absent marked.
#'
#' brms 2.23.0's conventions: the window is `lb <= Y <= ub`, so on a
#' discrete response its lower edge enters the CDF as `F(lb - 1)`, and
#' `lb <= 0` on a count is no truncation at all, which brms accepts. An
#' infinite bound is absent: it is evaluated at a point inside the
#' support and not read, so `lb = -Inf` costs no `Inf * 0` in a
#' derivative.
#'
#' @noRd
trunc_edges <- function(fam, av, n) {
  disc <- identical(fam[["type"]], "discrete")
  lb <- av[["trunc_lb"]]
  ub <- av[["trunc_ub"]]
  lo <- if (is.null(lb)) rep(-Inf, n) else rep_len(as.numeric(lb), n)
  if (disc) lo <- lo - 1
  hi <- if (is.null(ub)) rep(Inf, n) else rep_len(as.numeric(ub), n)
  list(lo = lo, hi = hi,
       # an NA bound is not an absent one: it stays NA, as the row does
       lo_abs = !is.na(lo) & (lo == -Inf | (disc & lo < 0)),
       hi_abs = !is.na(hi) & hi == Inf)
}

#' The per-row log-density of one response: the family's `lpdf` with
#' `cens()` and `trunc()` folded in, in that composition order.
#'
#' Factored out of the objective below so that `log_lik()` (R/loo.R) can
#' evaluate the SAME density numerically at each posterior draw. Two
#' copies would drift: an addition term that reshapes a row's density -
#' censoring replaces it by a CDF difference, truncation renormalizes it
#' by the window mass - has to be applied identically, or the LOO
#' estimate describes a different model than the one that was fitted.
#' Case weights stay OUTSIDE: the objective sums `w * ll`, and
#' `log_lik()` multiplies the returned vector the same way, which is
#' also brms's convention (`log_lik_weight()`).
#'
#' `yobs` is what the density is evaluated at - an `OBS()`-tagged
#' response on the tape, or an `mi()` response's observed-or-latent
#' values - while `yraw` is the recorded response the censoring CDFs
#' read.
#'
#' For a family that declares `cdf_log`, every censored or truncated
#' quantity is the mass of a window `(lo, hi]`, `log_window()` on the
#' log pairs of `fam_log_pair()` (a family without one is read by
#' `row_lpdf_prob()`):
#'
#' \tabular{ll}{
#'   left at `y` \tab `(lb', y]` \cr
#'   right at `y` \tab `(y, ub]` \cr
#'   interval `(y, y2]` \tab `(y, y2]` \cr
#'   the `trunc()` normalizer \tab `(lb', ub]`
#' }
#'
#' with `lb'` the lower edge `trunc_edges()` gives and an absent edge
#' open. The bounds are read as brms 2.23.0 reads them, a discrete
#' response included: right censoring at `y` is `P(Y > y)`, an interval
#' is `P(y < Y <= y2)`, and only the truncation bound `lb` is
#' inclusive, as `F(lb - 1)`. The censored rows read the same edges the
#' normalizer does, so a censored row under `trunc()` observes its
#' event inside the window. brms does not: it divides the event over
#' the whole line by the window's mass, which is not a probability
#' (brms-34 in dev/upstream-bugs.md), so the two differ for right
#' censoring with an upper bound and left censoring with a lower one.
#'
#' @noRd
row_lpdf <- function(fam, yobs, yraw, dpv, av, extra) {
  "[<-" <- RTMB::ADoverload("[<-")
  ll <- if (is.null(fam[["extra_pars"]])) {
    fam[["lpdf"]](yobs, dpv, av)
  } else {
    fam[["lpdf"]](yobs, dpv, av, extra)
  }
  # The backstop for the silent-wrong-answer case that
  # `required_aterms` closes by declaration, kept here because it needs
  # no declaration: an absent addition term reaches the density as NULL,
  # NULL in arithmetic gives numeric(0), and `sum(numeric(0))` is 0. The
  # fit then RETURNS, with a log-likelihood of zero and estimates that
  # are whatever the optimizer did on a flat surface. This is a LENGTH
  # check, not a value check, so it resolves while the tape is being
  # built and leaves no branch on it.
  if (length(ll) == 0L && length(yobs) > 0L) {
    frm_stop("The '", fam[["family"]], "' density returned no values for ",
             length(yobs), " observations. The usual cause is an addition ",
             "term the density reads that the model does not supply: an ",
             "absent one is NULL, and NULL in arithmetic is a zero-length ",
             "result, not an error. Declare it with ",
             "frmtmb_family(required_aterms =) to get this refusal by name",
             call. = FALSE,
             package = frm_family_package(fam))
  }
  cen <- av[["cens"]]
  trunc_on <- !is.null(av[["trunc_lb"]]) || !is.null(av[["trunc_ub"]])
  if (is.null(cen) && !trunc_on) return(ll)
  n <- length(yraw)
  if (length(ll) != n) ll <- ll + numeric(n)
  te <- trunc_edges(fam, av, n)
  trunc_on <- trunc_on && !all(te$lo_abs & te$hi_abs)
  if (is.null(fam[["cdf_log"]])) {
    return(row_lpdf_prob(fam, ll, yraw, dpv, av, extra, te, trunc_on))
  }
  i_l <- i_r <- i_i <- integer(0)
  if (!is.null(cen)) {
    i_r <- which(cen == 1)
    i_l <- which(cen == -1)
    i_i <- which(cen == 2)
  }
  none <- rep(FALSE, n)
  # Which tail each edge's pair must supply. An edge read as the LOWER
  # end of a window needs S there, as the upper end F, and both where
  # the window's other end is present too.
  need <- list(y = c(F = FALSE, S = FALSE), y2 = c(F = FALSE, S = FALSE),
               lo = c(F = FALSE, S = FALSE), hi = c(F = FALSE, S = FALSE))
  mark <- function(need, lo, hi, lo_abs, hi_abs, rows) {
    if (!length(rows)) return(need)
    la <- lo_abs[rows]
    ha <- hi_abs[rows]
    both <- any(!la & !ha)
    if (both) {
      need[[lo]][] <- TRUE
      need[[hi]][] <- TRUE
    }
    if (any(la & !ha)) need[[hi]][["F"]] <- TRUE
    if (any(!la & ha)) need[[lo]][["S"]] <- TRUE
    need
  }
  need <- mark(need, "lo", "y", te$lo_abs, none, i_l)
  need <- mark(need, "y", "hi", none, te$hi_abs, i_r)
  need <- mark(need, "y", "y2", none, none, i_i)
  if (trunc_on) {
    need <- mark(need, "lo", "hi", te$lo_abs, te$hi_abs, seq_len(n))
  }
  # An absent edge is evaluated at the recorded response, a point inside
  # the support, and never read
  at <- list(y = yraw, y2 = av[["cens_y2"]],
             lo = ifelse(te$lo_abs, yraw, te$lo),
             hi = ifelse(te$hi_abs, yraw, te$hi))
  pairs <- lapply(names(need), function(k) {
    if (!any(need[[k]])) return(NULL)
    fam_log_pair(fam, rep_len(at[[k]], n), dpv, av, extra, n)
  })
  names(pairs) <- names(need)
  if (length(i_l)) {
    ll[i_l] <- log_window(pairs$lo, pairs$y, te$lo_abs, none, i_l)
  }
  if (length(i_r)) {
    ll[i_r] <- log_window(pairs$y, pairs$hi, none, te$hi_abs, i_r)
  }
  if (length(i_i)) {
    # an interval reaching outside [lb, ub] is a data contradiction (the
    # response could not have been observed there at all), so it is not
    # cut to the window
    ll[i_i] <- log_window(pairs$y, pairs$y2, none, none, i_i)
  }
  if (trunc_on) {
    ll <- ll - log_window(pairs$lo, pairs$hi, te$lo_abs, te$hi_abs,
                          seq_len(n))
  }
  ll
}

#' A truncation bound evaluated at the censored rows. An absent bound is
#' the scalar 0 / 1 the untruncated case wants, and a present one is a
#' full-length CDF vector, so the two need different subsetting.
#'
#' @noRd
bound_rows <- function(v, i) if (length(v) == 1L) v else v[i]

#' `row_lpdf()` for a family with no `cdf_log`: its one-tail slots read
#' on the probability scale, as core read them through 0.72.0.
#'
#' A family that supplies `lcdf` and `lccdf` and nothing else is held to
#' the call pattern it was written for: the frmtmb.eam Wiener families
#' read which rows a call is for from the censoring codes in `aterms`
#' and fill the others with a neutral value, and their `lcdf` is a
#' boundary's defective distribution function on a left- or
#' interval-censored row, which is not the complement of `lccdf`. So
#' each slot is called where it was: `lccdf` at the recorded responses
#' for the right-censored rows, `lcdf` there for the left- and
#' interval-censored ones and at `cens_y2`, and `lcdf` at the bounds for
#' the window `F(ub) - F(lb')`. Two things are new, and both are brms's
#' likelihood: no censoring bound on a count is shifted, and
#' `trunc(lb = )` with no upper bound divides by the log survivor
#' `lccdf(lb')`, which `log(1 - F(lb'))` was not once `F(lb')` rounded to
#' one (a weibull fit stopped at a NaN gradient there). That call is
#' made without the censoring codes, because the window holds for every
#' row and a slot that reads them would fill the uncensored rows with
#' its neutral value. An absent bound (`trunc_edges()`) is evaluated at
#' the response and not read.
#'
#' @noRd
row_lpdf_prob <- function(fam, ll, yraw, dpv, av, extra, te, trunc_on) {
  "[<-" <- RTMB::ADoverload("[<-")
  n <- length(yraw)
  disc <- identical(fam[["type"]], "discrete")
  # a bound as the slot always saw it, the scalar a constant bound is,
  # unless some rows of it are absent
  edge <- function(given, at, abs) {
    if (any(abs)) ifelse(abs, yraw, at) else given
  }
  absent_to <- function(v, abs, value) {
    if (!any(abs)) return(v)
    v <- v + numeric(n)
    v[which(abs)] <- value
    v
  }
  lo_on <- trunc_on && !all(te$lo_abs)
  hi_on <- trunc_on && !all(te$hi_abs)
  Flb <- 0
  Fub <- 1
  if (hi_on) {
    Fub <- absent_to(fam_lcdf(fam, edge(av[["trunc_ub"]], te$hi, te$hi_abs),
                              dpv, av, extra), te$hi_abs, 1)
  }
  if (lo_on) {
    ql <- edge(av[["trunc_lb"]] - as.numeric(disc), te$lo, te$lo_abs)
    Flb <- absent_to(fam_lcdf(fam, ql, dpv, av, extra), te$lo_abs, 0)
  }
  cen <- av[["cens"]]
  if (!is.null(cen)) {
    # censoring codes are data, so grouped sub-assignment replaces the
    # density on censored rows without parameter branching
    i_r <- which(cen == 1)
    i_l <- which(cen == -1)
    i_i <- which(cen == 2)
    lS <- if (has_lccdf(fam) && length(i_r)) {
      fam_lccdf(fam, yraw, dpv, av, extra)
    }
    need_F <- length(i_l) || length(i_i) || (length(i_r) && is.null(lS))
    Fv <- if (need_F) fam_lcdf(fam, yraw, dpv, av, extra)
    if (length(i_r)) {
      ll[i_r] <- if (is.null(lS)) {
        log(bound_rows(Fub, i_r) - Fv[i_r])
      } else if (!hi_on) {
        lS[i_r]
      } else {
        # P(y < Y <= ub) = S(y) - S(ub), as a log difference, so that the
        # exact log survivor is never round-tripped through a probability
        lSu <- absent_to(fam_lccdf(fam, edge(av[["trunc_ub"]], te$hi,
                                             te$hi_abs), dpv, av, extra),
                         te$hi_abs, -Inf)
        lS[i_r] + log1p(-exp(bound_rows(lSu, i_r) - lS[i_r]))
      }
    }
    if (length(i_l)) ll[i_l] <- log(Fv[i_l] - bound_rows(Flb, i_l))
    if (length(i_i)) {
      F2 <- fam_lcdf(fam, av[["cens_y2"]], dpv, av, extra)
      ll[i_i] <- log(F2[i_i] - Fv[i_i])
    }
  }
  if (trunc_on) {
    if (!hi_on && has_lccdf(fam)) {
      av_w <- av
      av_w[["cens"]] <- NULL
      av_w[["cens_y2"]] <- NULL
      ll <- ll - absent_to(fam_lccdf(fam, ql, dpv, av_w, extra), te$lo_abs, 0)
    } else {
      ll <- ll - log(Fub - Flb)
    }
  }
  ll
}

#' The simplex of a `mo()` term from its `D - 1` free coordinates `z`.
#'
#' `z` is a point of the plane that stereographic projection maps onto
#' the unit sphere in `D` dimensions, and the simplex is the squared
#' coordinates of that sphere point, `w = u^2`. The projection's pole is
#' `-(1, ..., 1) / sqrt(D)`, so `z = 0` is the barycenter, the start a
#' softmax at zero also gave.
#'
#' Why not the softmax, `w = exp(c(0, z)) / sum(.)`: a weight of 0 is
#' at an infinite coordinate there, and the gradient toward a vertex is
#' of the order of the weight that remains. A fit whose simplex runs to
#' a vertex stops on that plateau with code 0, wherever the likelihood
#' still rises: on brms_monotonic's `ls ~ mo(income) * age`, 84 of 200
#' data sets stopped more than 1e-6 below the exact maximum, 53 by more
#' than 0.01 and up to 1.95 (dev/optima-findings.md). Here every face
#' of the simplex is at a finite point, `u_j = 0`, where `w_j = u_j^2`
#' folds: the likelihood is quadratic in `u_j` there, with curvature
#' twice the derivative along the face's normal, so a face the
#' likelihood should leave is a saddle the optimizer leaves, and a face
#' that is the maximum is an ordinary minimum of the objective with a
#' finite Hessian. The one point at infinity is the pole, whose image
#' is the barycenter; `mo_search()` brings a fit that drifted toward it
#' back to the chart near zero.
#'
#' `chart = "softmax"` is the sampler's: frmtmb.sample sets it on its
#' own copy of the frame (`mo_chart_frame()`). A sampler needs a density
#' on the coordinates, and the folds `u_j = 0` that serve the optimizer
#' are zeros of this chart's Jacobian, so a uniform density on the
#' simplex is walls of zero density between the chart's `2^D` sheets,
#' which chains do not cross. The softmax is one sheet with no fold.
#'
#' @noRd
mo_simplex <- function(z, chart = NULL) {
  if (identical(chart, "softmax")) {
    "c" <- RTMB::ADoverload("c")
    e <- exp(c(0, z))
    return(e / sum(e))
  }
  D <- length(z) + 1L
  P <- -1 / sqrt(D)
  r2 <- sum(z * z)
  u <- (2 * as.vector(mo_sphere_basis(D) %*% z) + (r2 - 1) * P) / (r2 + 1)
  w <- u * u
  # one division keeps the weights summing to 1 to rounding, as the
  # softmax's did; the top category reads D * sum(w)
  w / sum(w)
}

#' An orthonormal basis of the plane orthogonal to `(1, ..., 1)`, the
#' `D x (D - 1)` matrix `mo_simplex()` projects through.
#'
#' @noRd
mo_sphere_basis <- function(D) {
  H <- stats::contr.helmert(D)
  sweep(H, 2L, sqrt(colSums(H^2)), "/")
}

#' The coordinates `mo_simplex()` maps to the simplex `w`, taken on the
#' sheet with every `u_j = sqrt(w_j) >= 0`. That sheet lies inside the
#' unit ball of the coordinates, so these are the canonical coordinates
#' of `w`: every other preimage differs by the sign of some `u_j`, and
#' lies farther out.
#'
#' `chart = "softmax"` gives the log ratios to the first weight, with
#' every weight held at `floor` or above first, because a face is at an
#' infinite softmax coordinate and a sampler's start must be finite.
#'
#' @noRd
mo_coords <- function(w, chart = NULL, floor = 1e-8) {
  if (identical(chart, "softmax")) {
    w <- pmax(w, floor)
    return(log(w[-1L] / w[1L]))
  }
  D <- length(w)
  u <- sqrt(pmax(w, 0) / sum(pmax(w, 0)))
  as.numeric(crossprod(mo_sphere_basis(D), u)) / (1 + sum(u) / sqrt(D))
}

#' A copy of `frame` whose simplexes read their coordinates in `chart`
#' (see `mo_simplex()`): every `mo()` simplex, and every simplex a
#' response's family declares (`frmtmb_family(simplex_pars =)`), which
#' the family reads from the response's addition term `<name>_chart`.
#' The objective built from it is the same model in other coordinates;
#' `simplex_frame_terms()` lists where the coordinates sit.
#'
#' @noRd
mo_chart_frame <- function(frame, chart = "softmax") {
  for (k in seq_along(frame[["linpreds"]])) {
    mos <- frame[["linpreds"]][[k]][["mo"]]
    if (!length(mos)) next
    for (j in seq_along(mos)) mos[[j]][["chart"]] <- chart
    frame[["linpreds"]][[k]][["mo"]] <- mos
  }
  for (tm in family_simplex_terms(frame)) {
    frame[["aterm_values"]][[tm$resp]][[tm$chart_key]] <- chart
  }
  frame
}

#' Every simplex of a fit's frame, in the shape `mo_frame_terms()` gives
#' a `mo()` term: the `mo()` simplexes, then each simplex a response's
#' family declares, under its brms name, `sbhaz[k]` for a `cox()`
#' baseline (`sbhaz_<resp>[k]` in a multivariate model). `kind` is the
#' declaring family's name. The sampler reads every one in softmax
#' coordinates under brms's default `dirichlet(1)`.
#'
#' @noRd
simplex_frame_terms <- function(fit) {
  out <- mo_frame_terms(fit)
  for (tm in family_simplex_terms(fit$frame)) {
    out[[length(out) + 1L]] <- list(
      zeta = tm$zeta, D = tm$D,
      names = paste0(tm$name, if (tm$multi) paste0("_", tm$resp), "[",
                     seq_len(tm$D), "]"),
      kind = tm$family)
  }
  out
}

#' The simplexes the responses' families declare
#' (`frmtmb_family(simplex_pars =)`), one entry each: the response, the
#' brms name, the template component `zeta` that holds the `D - 1`
#' coordinates, `D`, the addition term that carries the chart, and the
#' family's name.
#'
#' @noRd
family_simplex_terms <- function(frame, spec = frame[["spec"]]) {
  rs <- spec$responses %||% list()
  multi <- length(rs) > 1L
  out <- list()
  for (resp in rs) {
    sp <- resp$family[["simplex_pars"]]
    if (is.null(sp)) next
    r <- resp$resp_name
    for (nm in names(sp)) {
      zeta <- extra_tpl_name(frame, if (multi) r, sp[[nm]])
      out[[length(out) + 1L]] <- list(
        resp = r, name = nm, zeta = zeta,
        D = length(frame[["par_template"]][[zeta]]) + 1L,
        chart_key = paste0(nm, "_chart"), family = resp$family[["family"]],
        multi = multi)
    }
  }
  out
}

#' Every `mo()` term of a frame: its simplex component, `D`, and the
#' names brms gives its simplex, `simo<prefix>_<label>1[k]`, after
#' brms's `rename_sp()`, where the prefix is that of the term's
#' predictor (`_sigma` for a simplex in sigma's formula, empty in mu's).
#'
#' @noRd
mo_frame_terms <- function(fit) {
  out <- list()
  for (lp in fit$frame[["linpreds"]]) {
    pre <- brms_lp_prefix(fit, lp)
    for (mi in lp[["mo"]] %||% list()) {
      lab <- brms_rename(mi[["label"]])
      nm <- paste0("simo", if (nzchar(pre)) paste0("_", pre), "_", lab,
                   "1[", seq_len(mi[["D"]]), "]")
      out[[length(out) + 1L]] <- list(zeta = mi[["zeta"]], D = mi[["D"]],
                                      names = nm)
    }
  }
  out
}

#' The values `D * cumsum(w)[code + 1]` of one `mo()` term at the
#' coordinates `z`, with no multiplier: the transform a population-level
#' coefficient and the group-level ones that vary it multiply.
#'
#' @noRd
mo_transform <- function(mi, z) {
  mo_min(mi)   # refuses a frame that frmtmb 0.75.0 or earlier coded
  "c" <- RTMB::ADoverload("c")
  w <- mo_simplex(z, mi[["chart"]])
  cz0 <- c(0, cumsum(w))
  mi$D * cz0[mi$codes + 1L]
}

#' What the group-level `mo()` columns add to `Z b` (AD-safe).
#'
#' `lp$Z` holds such a column as its placeholder value (1, or the
#' multiplier of an interaction), so `Z b` counts the column's effect
#' once per row; the effect is the coefficient times the transform, so
#' this adds it times the transform less one. Kept apart from `Z b`, so
#' the tape of a model without such a term does not change.
#'
#' @noRd
mo_re_zterm <- function(lp, pars, bvec) {
  out <- 0
  for (zm in lp[["Z_mo"]]) {
    mi <- lp[["mo"]][[zm$mo]]
    m <- mo_transform(mi, pars[[mi$zeta]])
    out <- out + as.vector(lp[["Z"]][, zm$cols, drop = FALSE] %*%
                             bvec[zm$cols]) * (m - 1)
  }
  out
}

#' The part of one linear predictor that the random-effect coefficients
#' do not enter: `X beta`, the offset, and the `mo()`, `mi()` and `me()`
#' terms.
#' `zterm` is the block contribution `Z b`, added in the position it has
#' always occupied so that the fitted Laplace tape is unchanged to the
#' last bit.
#'
#' Split out of the loop below because the importance correction
#' (R/importance.R) needs exactly this piece and nothing else: it is the
#' same for every importance draw, so it is formed once at length `n`
#' and expanded, while only the `Z u` contribution varies from draw to
#' draw. One copy of the arithmetic means the corrected objective cannot
#' drift away from the objective it corrects.
#'
#' @noRd
lp_eta_fixed <- function(lp, pars, n, mivals, yfall, zterm = NULL,
                         mevals = NULL) {
  # The overloads have to be re-established here, and this is the same
  # trio frmtmb_ad_overload() installs (R/ad-env.R). This code used to
  # sit inside the objective closure, which sets them once at the top;
  # a plain base `c(0, advector)` silently drops the advector class,
  # the mo() simplex is built exactly that way, and the tape then ends
  # in a bare complex vector instead of a value. Spelled out rather
  # than wrapped because a wrapper rewrites the body at load time,
  # which is worth it for an extension author's function and not for
  # one line of package-internal arithmetic.
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  "diag<-" <- RTMB::ADoverload("diag<-")
  # as.vector, not drop: it collapses both Matrix and advector results
  eta <- if (ncol(lp[["X"]])) {
    as.vector(lp[["X"]] %*% pars[[lp[["par"]]]][lp[["idx"]]])
  } else {
    # threshold-only ordinal model; the design's rows, because a
    # subset() response has fewer than the frame
    rep(0, nrow(lp[["X"]]))
  }
  if (!is.null(zterm)) {
    eta <- eta + zterm
  }
  if (!is.null(lp[["offset"]])) {
    eta <- eta + lp[["offset"]]
  }
  # monotonic terms: scale coefficient (in beta, zero X column)
  # times D * cumulative simplex at the observed category
  for (mi in lp[["mo"]] %||% list()) {
    mo_min(mi)   # refuses a frame that frmtmb 0.75.0 or earlier coded
    zeta <- mo_simplex(pars[[mi$zeta]], mi[["chart"]])
    cz0 <- c(0, cumsum(zeta))
    term <- mi$D * cz0[mi$codes + 1L]
    if (!is.null(mi$mult)) term <- term * mi$mult
    eta <- eta + pars[[lp[["par"]]]][lp[["idx"]][mi$col]] * term
  }
  # mi(x) terms: coefficient times the observed-or-latent values
  for (mt in lp[["mi"]] %||% list()) {
    xv <- mivals[[mt$var]] %||% yfall[[mt$var]]
    # mi(x, idx = ): the rows of x that this predictor's rows read
    if (!is.null(mt$idxl)) xv <- xv[mt$idxl]
    if (!is.null(mt$mult)) xv <- xv * mt$mult
    eta <- eta + pars[[lp[["par"]]]][lp[["idx"]][mt$col]] * xv
  }
  # me() terms: coefficient times the product of latent values
  for (mt in lp[["me"]] %||% list()) {
    eta <- eta + pars[[lp[["par"]]]][lp[["idx"]][mt$col]] *
      me_col_value(mt, mevals)
  }
  eta
}

#' Turns an assembled `frmtmb_frame` into the negative log-likelihood
#' closure that `RTMB::MakeADFun()` tapes. The returned function takes the
#' parameter list and gives back one value, with the design matrices, the
#' response, and the addition terms baked in as constants. It is a closure
#' and not compiled code because RTMB tapes plain R: the model structure
#' is resolved once here, outside the function the tape records, so the
#' taped code loops over data only through vectorized operations.
#'
#' @noRd
build_objective <- function(frame) {
  require_model_packages(frame)
  require_structure_packages(frame)
  lps <- frame[["linpreds"]]
  blocks <- frame[["re_blocks"]]
  block_fns <- lapply(blocks, function(bk) {
    cs_entry(bk[["covstruct"]])$nll
  })
  # Non-centered blocks (frm_sample(reparameterize = TRUE) only, set on
  # a private copy of the frame): the sampled vector for these is z, and
  # the coefficients b = L(theta) z are computed on the tape. Empty on
  # every fitted objective, where the Laplace approximation makes the
  # parameterization of the integrated variable irrelevant.
  ncp_idx <- frame[["ncp_blocks"]] %||% integer(0)
  is_ncp <- seq_along(blocks) %in% ncp_idx
  # asked once, not per evaluation, and asked of the BLOCKS: a frame
  # that predates `has_expand` must not be read as "nothing to expand"
  needs_expand <- frame_needs_expand(frame)
  spec <- frame[["spec"]]
  resps <- spec$responses
  rescor <- isTRUE(spec$rescor)
  # the response carrying the one nu of a Student-t rescor model
  rescor_nu <- spec[["rescor_nu"]]
  y <- frame[["y"]]
  atv <- frame[["aterm_values"]]
  n <- frame[["n_obs"]]
  acs <- frame[["autocor"]] %||% list()
  # the registry entries of the residual terms, read once per tape
  ac_ents <- lapply(acs, function(ac) {
    if (autocor_is_seam(ac)) ac_block_entry(ac)
  })

  extra_names <- frame[["extra_names"]] %||% character(0)
  me_fr <- frame[["me"]]

  # Cluster-robust scores (R/sandwich.R) need the per-cluster pieces of
  # the objective as a function of one extra parameter each, so that
  # `obj$gr()` at a cluster indicator IS that cluster's score. The hook
  # is a per-response row -> cluster index; it is NULL for every
  # ordinary fit and is only ever set on a private copy of the frame,
  # never on the one stored in the fit. vcov_cluster() refuses every
  # likelihood whose contribution is not the per-row `sum(w * ll)`
  # below (autocor, hmm, mixtures, rescor), so this is the only place a
  # cluster weight has to enter.
  clw_idx <- frame[["cluster_w"]]

  function(pars) {
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    nll <- 0

    # `bfull` is the coefficient vector every predictor below reads. It
    # IS pars$b on the centered route, so nothing about that tape moves.
    bfull <- pars[["b"]]
    for (i in ncp_idx) {
      bk <- blocks[[i]]
      bfull[bk[["b_idx"]]] <- ncp_scale_b(bk, pars[["b"]][bk[["b_idx"]]],
                                     pars[["theta"]][bk[["theta_idx"]]])
    }

    for (i in seq_along(blocks)) {
      bk <- blocks[[i]]
      if (is_ncp[i]) {
        # the block's whole density: the transform took the entire
        # factor out, so what is left is standard normal
        nll <- nll - sum(RTMB::dnorm(pars[["b"]][bk[["b_idx"]]], 0, 1,
                                     log = TRUE))
        next
      }
      # Student-t blocks of one grouping factor share a latent scale
      # (student_share_sets()) and are one density, taken at the first
      sh <- bk[["t_share"]]
      if (!is.null(sh)) {
        if (i == sh[1L]) {
          nll <- nll - student_shared_lpdf(blocks[sh], bfull,
                                           pars[["theta"]])
        }
        next
      }
      nll <- nll - block_fns[[i]](bfull[bk[["b_idx"]]],
                                  pars[["theta"]][bk[["theta_idx"]]], bk)
    }

    extra <- NULL
    if (length(extra_names)) {
      extra <- lapply(stats::setNames(extra_names, extra_names),
                      function(nm) pars[[nm]])
    }
    # each response's density reads its own family extras under the
    # family's names; identical to `extra` outside a multivariate frame
    extra_r <- lapply(stats::setNames(nm = names(resps)), function(r) {
      resp_extras(frame, extra, r)
    })

    # coefficient-space vector for the Z products (rr blocks expand
    # their factors through the loadings, esicar blocks center)
    bvec <- if (needs_expand) {
      expand_b(frame, bfull, pars[["theta"]])
    } else {
      bfull
    }

    # mi(): observed-or-latent value vectors per imputation response
    mivals <- list()
    for (vn in names(frame[["mi_map"]] %||% list())) {
      mm_ <- frame[["mi_map"]][[vn]]
      xv <- y[[vn]]
      xv[mm_$rows] <- pars[["miss"]][mm_$idx]
      mivals[[vn]] <- xv
      if (!is.null(mm_$se)) {
        # measurement model: observed values scatter around the latent
        # truth with known SD (brms me())
        nll <- nll - sum(RTMB::dnorm(y[[vn]][mm_$obs], xv[mm_$obs],
                                     mm_$se[mm_$obs], log = TRUE))
      }
    }

    # me(): the latent values and their measurement and latent
    # densities (R/me.R)
    mevals <- NULL
    if (!is.null(me_fr)) {
      mel <- me_latent(me_fr, pars)
      nll <- nll - mel$ll
      mevals <- mel$values
    }

    dparv <- list()
    for (lp in lps) {
      if (!is.null(lp[["nl_body"]])) {
        # nonlinear predictor: arbitrary R code over the nlpar values and
        # raw data columns, evaluated straight onto the tape. The linear
        # predictors were put in dependency order at parse time, so
        # every parameter this body names already has a value; only the
        # ones it names are in scope, which keeps a data column that
        # shares a name with another dpar meaning the column - what the
        # same body means in brms. `nl_dpar_refs` are the names that had
        # no column behind them and read another parameter's per-row
        # value instead (a variance function of the fitted mean).
        ev <- c(dparv[[lp[["resp"]]]][c(lp[["nl_pars"]], lp[["nl_dpar_refs"]])],
                lp[["data_list"]],
                # ps() blocks: one closure per term, built from the
                # CURRENT coefficients. It has to be rebuilt here, per
                # call, rather than captured in `nl_env`, which parse
                # time fixed before any coefficient existed.
                ps_env(lp, pars, bvec))
        # the body is taped once, so the handler costs nothing per
        # gradient; it exists because a body name that resolved to an
        # environment object instead of a column fails here, far from
        # the cause
        eta <- tryCatch(eval(lp[["nl_body"]], ev,
                             ad_overload_env(lp[["nl_env"]], lp[["nl_body"]])),
                        error = function(e) nl_body_error(e, lp))
        dparv[[lp[["resp"]]]][[lp[["dpar"]]]] <- lp[["link"]]$linkinv(eta)
        dparv[[lp[["resp"]]]][[paste0(".eta_", lp[["dpar"]])]] <- eta
        next
      }
      zterm <- if (!is.null(lp[["Z"]])) as.vector(lp[["Z"]] %*% bvec)
      if (length(lp[["Z_mo"]])) zterm <- zterm + mo_re_zterm(lp, pars, bvec)
      eta <- lp_eta_fixed(lp, pars, n, mivals, y, zterm = zterm,
                          mevals = mevals)
      # cs(x) terms: n x (K-1) threshold-specific offsets, consumed by
      # the sequential ordinal lpdfs through dpars$.cs (an ordinal
      # mixture's component k through `.cs_mu<k>`, cs_slot())
      if (length(lp[["cs"]] %||% list())) {
        CS <- 0
        for (ct in lp[["cs"]]) {
          bcs <- pars[[ct$par]]
          CS <- CS + RTMB::matrix(ct$vals, length(ct$vals), 1) %*%
            RTMB::matrix(bcs, 1, length(bcs))
        }
        dparv[[lp[["resp"]]]][[cs_slot(lp[["dpar"]])]] <- CS
      }
      dparv[[lp[["resp"]]]][[lp[["dpar"]]]] <- lp[["link"]]$linkinv(eta)
      # The linear predictor rides along beside the inverse-linked value
      # under a reserved `.eta_` name. An inverse link saturates -
      # plogis(40) is exactly 1, exp(-800) is exactly 0 - so a density
      # that recomputes 1 - mu or log(mu) from the dpar reads -Inf or
      # NaN with a useless gradient in a region the linear predictor
      # itself describes perfectly well. The families that can use the
      # eta scale read it from here through the public accessors
      # (`dpar_log()`, `dpar_log_complement()`, `dpar_complement()`),
      # and undoing the link to recover it would defeat the purpose.
      # The entry is RESERVED rather than API: what it holds is on the
      # link scale, so its meaning depends on the dpar's link, which is
      # why every accessor takes that link and why a density must not
      # read the entry itself. Only the taped objective supplies it:
      # the numeric post-fit paths (fitted(), simulate(), the CDFs)
      # pass the dpar values alone, and there nothing is differentiated,
      # so the saturation is harmless.
      dparv[[lp[["resp"]]]][[paste0(".eta_", lp[["dpar"]])]] <- eta
    }

    # brms's cov = FALSE ARMA: mu moves by a regression on the earlier
    # residuals of its group, on the response scale (brms adds it after
    # the inverse link), and every density below, rowwise or rescor,
    # then reads the shifted mu. The residual is taken against the
    # observed-or-imputed response, brms's Yl under mi().
    for (r in names(acs)) {
      ac <- acs[[r]]
      if (!autocor_is_cond(ac)) next
      mu_r <- dparv[[r]]$mu
      if (length(mu_r) == 1L) mu_r <- mu_r + numeric(n)
      yl <- mivals[[r]] %||% y[[r]]
      dparv[[r]]$mu <- mu_r +
        autocor_cond_shift(yl - mu_r, pars[["thetaac"]][ac[["theta_idx"]]],
                           ac)
      # the link-scale value no longer describes mu, and a density
      # reading it through the public accessors would miss the shift
      dparv[[r]][[".eta_mu"]] <- NULL
    }

    if (rescor) {
      # gaussian joint likelihood: standardized residuals against a
      # constant correlation matrix keeps this vectorized even with
      # observation-varying sigmas
      K <- length(resps)
      zvec <- NULL
      logsig <- 0
      for (r in names(resps)) {
        zr <- (y[[r]] - dparv[[r]]$mu) / dparv[[r]]$sigma
        zvec <- if (is.null(zvec)) zr else c(zvec, zr)
        logsig <- logsig + sum(log(dparv[[r]]$sigma))
      }
      Zstd <- RTMB::matrix(zvec, n, K)
      C <- us_chol_cor(pars[["thetar"]], K)
      if (is.null(rescor_nu)) {
        nll <- nll - sum(RTMB::dmvnorm(Zstd, 0, C, log = TRUE)) + logsig
      } else {
        # brms's multi_student_t(nu, Mu, Sigma) with Sigma = D C D: the
        # scale of row i is sigma_i, so the standardized rows are
        # multivariate t with scale C, and logsig is the same Jacobian
        nll <- nll - mvt_std_loglik(Zstd, C, dparv[[rescor_nu]]$nu[1]) +
          logsig
      }
    } else {
      for (r in names(resps)) {
        fam <- resps[[r]]$family
        w <- atv[[r]]$weights %||% 1
        if (!is.null(clw_idx)) w <- w * pars[["clw"]][clw_idx[[r]]]
        if (!is.null(acs[[r]]) && !autocor_is_cond(acs[[r]])) {
          # R-side residual correlation: the response's density is a
          # joint (multivariate normal / t) one per group, not a
          # product over rows, so it replaces fam$lpdf entirely. Every
          # aterm that would reshape a per-row contribution was refused
          # at frame assembly, so nothing below this point applies.
          ac <- acs[[r]]
          if (autocor_is_seam(ac)) {
            # a registry term (fcor(), sar()): the density of its map
            mu_r <- dparv[[r]]$mu
            if (length(mu_r) == 1L) mu_r <- mu_r + numeric(n)
            nu <- if (isTRUE(ac[["student"]])) dparv[[r]]$nu[1] else NULL
            nll <- nll - autocor_seam_loglik(
              ac_ents[[r]], y[[r]], mu_r, dparv[[r]]$sigma,
              pars[["thetaac"]][ac[["theta_idx"]]], ac, nu)
            next
          }
          tac <- pars[["thetaac"]][ac[["theta_idx"]]]
          R <- autocor_cor(tac, ac)
          # D R D with D the marginal sds; under ar(), ma() and arma()
          # sigma is the innovation sd, as in brms
          sg <- autocor_marg_sigma(dparv[[r]]$sigma, tac, ac)
          z <- (y[[r]] - dparv[[r]]$mu) / sg
          lsig <- if (length(sg) == 1L) n * log(sg) else sum(log(sg))
          nu <- if (isTRUE(ac[["student"]])) dparv[[r]]$nu[1] else NULL
          nll <- nll - autocor_loglik(z, R, ac, lsig, nu)
          next
        }
        # A structure that supplies `loglik` REPLACES the rowwise
        # density; one that declares only capabilities (lca) leaves the
        # rowwise path below alone, which is why the branch is on the
        # slot and not on the structure.
        stll <- fam_structure(fam)[["loglik"]]
        if (!is.null(stll)) {
          # Such a likelihood does not factorize over rows: a hidden
          # state sequence is summed out by a forward recursion per
          # SEQUENCE, a latent-class mixture sums a group's
          # per-observation densities BEFORE the logsumexp over
          # classes. Either way there is no per-row factor left, so none
          # of the aterm machinery below applies - the frame refused
          # every term that would have needed it. The dispatch is on
          # DATA (the family object and the frame block), so it resolves
          # while the tape is being built and puts no branch on it.
          nll <- nll - stll(y[[r]], dparv[[r]], atv[[r]], w,
                            frame_block_of(frame, r), extra_r[[r]])
          next
        }
        # OBS() drives simulation/OSA machinery, but registers data under
        # the DEPARSED ARGUMENT EXPRESSION: in this loop every response
        # would collide on "y[[r]]" and silently swap data. Univariate
        # fits only; it also has no matrix method.
        yobs <- if (!is.null(mivals[[r]])) {
          # mi() response: the density is evaluated at the
          # observed-or-latent values (the latent entries' Gaussian
          # contribution IS this term)
          mivals[[r]]
        } else if (length(resps) == 1L && !is.matrix(y[[r]])) {
          # a stable symbol gives OSA machinery its observation.name
          .frm_obs <- y[[r]]
          RTMB::OBS(.frm_obs)
        } else {
          y[[r]]
        }
        ll <- row_lpdf(fam, yobs, y[[r]], dparv[[r]], atv[[r]], extra_r[[r]])
        nll <- nll - sum(w * ll)
      }
    }
    nll
  }
}
