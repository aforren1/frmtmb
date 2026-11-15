# Refusal of a weighted likelihood (or posterior) that has no maximum
# (or does not integrate) because every row owns a random effect.

# The dpars that concentrate each family's density on one point when
# they go to their limit together. Each family was measured to run away
# with an observation-level intercept and weights of 1.2 and to fit with
# weights of 0.9 (dev/fuzzfix-spike-families.R, dev/fuzzfix-fr1-families.R).
# exgaussian needs both its scales: with beta held, its density stays
# below 1 / beta. zero_inflated_asym_laplace is asym_laplace on its
# nonzero rows (dev/rel073-zial-probe.R: Laplace logLik 278.05 and
# 201.90 unweighted, quadrature bounded at unit weights and unbounded
# at 1.2, a group of 4 rows bounded).
unbounded_spike_dpar <- list(
  gaussian = "sigma", student = "sigma", lognormal = "sigma",
  Gamma = "shape", weibull = "shape", beta = "phi",
  inverse.gaussian = "shape", skew_normal = "sigma",
  shifted_lognormal = "sigma", huber = "sigma", von_mises = "kappa",
  frechet = "nu", gen_extreme_value = "sigma",
  exgaussian = c("sigma", "beta"), asym_laplace = "sigma",
  zero_inflated_asym_laplace = "sigma",
  # lognormal, Gamma and Beta on the rows off the point mass
  # (unbounded_point_mass); dev/wcheck2-4968.R profiles each
  hurdle_lognormal = "sigma",
  hurdle_gamma = "shape",
  zero_inflated_beta = "phi",
  zero_one_inflated_beta = "phi")

# The responses at which a family puts a point mass that its density
# reads off a probability dpar alone (hu, zi, zoi and coi): such a row's
# density neither depends on mu nor on the spike dpar, so its own effect
# keeps only its prior, and the row neither grows nor curves as the
# spike dpar goes to its limit. The excess is the other rows'.
unbounded_point_mass <- list(
  zero_inflated_asym_laplace = 0, hurdle_lognormal = 0, hurdle_gamma = 0,
  zero_inflated_beta = 0, zero_one_inflated_beta = c(0, 1))

# The limit each spike dpar runs to: a scale shrinks to 0, a precision
# or shape grows without bound
unbounded_spike_toward_zero <- c("sigma", "beta")

# How fast each density narrows as its spike dpar k grows: the width
# goes like k^-rate. A shape or precision that sets a variance (Gamma,
# Beta, inverse.gaussian, von_mises) gives 1/2; weibull's and frechet's
# shape set the scale of log y, which gives 1. A scale shrinking to 0
# is its own width (rate 1). The likelihood then grows like
# k^(rate * excess), which is what a prior's tail is weighed against.
unbounded_spike_rate <- list(
  Gamma = c(shape = 0.5), beta = c(phi = 0.5),
  hurdle_gamma = c(shape = 0.5), zero_inflated_beta = c(phi = 0.5),
  zero_one_inflated_beta = c(phi = 0.5),
  inverse.gaussian = c(shape = 0.5), von_mises = c(kappa = 0.5),
  weibull = c(shape = 1), frechet = c(nu = 1))

# Which covariance structures can give a row its own effect at any value
# with a finite density is a property of the structure, read through
# cs_full_rank() (R/covstruct.R): a built-in entry states it, and a
# registered structure states it or has it read off its covariance.

#' Refuse a weighted likelihood that has no maximum, or a posterior that
#' does not integrate.
#'
#' When every row with a positive weight has a random effect of its own
#' in `mu`, each row's density can be centered on its own observation.
#' As the spike dpars go to their limit together (`sigma` to 0, `shape`
#' to infinity) each density at its center grows like 1 / width, and the
#' marginal likelihood, exact or Laplace, grows like width^(n - sum(w)):
#' a row contributes w log(1 / width) from its density and
#' -log(1 / width) from the curvature of its own effect. With
#' `sum(w) > n` it has no maximum, under ML, REML and `profile = TRUE`
#' alike. asym_laplace is worse under the Laplace approximation: its log
#' density is linear on each side of the mode, the data give a row's own
#' effect no curvature, and the approximation grows like
#' width^(-sum(w)) for any weights.
#'
#' The proof needs three things, and the check fails open without any
#' of them: (i) the spike dpars reach their limit together on the rows
#' in question, which holds on every row when each design has a column
#' and spans the constant (an intercept), and on the rows of one cell
#' when the family has one spike dpar whose design has a column that is
#' constant on that cell and 0 elsewhere (`sigma ~ fac`); (ii) each
#' row's mean reaches its observation, which an owned effect in a
#' full-rank block gives; (iii) the effects' prior density is positive
#' and finite there.
#' It also fails open for `se()`, autocorrelation, `rescor` and `mi()`.
#'
#' Each response of a multivariate model is checked on its own rows. Its
#' density is one factor of the likelihood, and the other responses'
#' factors stay finite and positive while its own runs away, whatever
#' block its effects share with theirs.
#'
#' Under `cens()` the argument runs over the uncensored rows. A censored
#' row's probability is at most 1, and its own effect can bring it to 1,
#' so the exact marginal of such a row stays between a positive constant
#' and 1 while the uncensored rows run away as before: the excess is
#' theirs alone. Under the Laplace approximation a censored row costs
#' nothing when its effect's prior mean already puts its probability at
#' 1, and at most half a power of the growth otherwise (its own effect
#' then sits where the probability turns from 0 to 1, where the
#' curvature grows like 1 / width). When the mean's design spans the
#' constant, there is no `trunc()`, and no prior bounds the constant on
#' that side (`unbounded_shift_sides()`), a constant mean far enough up
#' (down) puts every right (left) censored row's probability at 1, so
#' only the other side's rows and the interval rows cost 1/2 each; the
#' check takes the cheaper side. Otherwise every censored row costs 1/2,
#' which fails open where that is too much. On 80 gaussian rows with 16
#' right censored (dev/wcheck-cens-path.R) the Laplace objective along
#' such a mean grew with a slope of exactly the uncensored excess, 0.5
#' to 7.5, where the 0.81.0 allowance of 1/2 per row let all of them
#' through.
#'
#' A hurdle or zero-inflated family puts some responses in a point mass
#' whose probability is a dpar of its own (`unbounded_point_mass`). Such
#' a row's density reads neither mu nor the spike dpar, so its own effect
#' keeps only its prior: the row neither grows nor costs, and the excess
#' is the other rows'.
#'
#' A prior on a spike dpar exempts the model only when it provably
#' bounds the growth (`prior_spike_verdict()`); a prior on anything else
#' never does. `route = "sample"` asks the posterior question instead:
#' the density must integrate, which a prior finite at the limit allows
#' only while sum(w) - n stays below 1.
#'
#' @param exact `TRUE` when the integral over the effects is exact
#'   (`quadrature = TRUE`, or full sampling), `FALSE` under Laplace.
#' @noRd
check_weights_unbounded <- function(spec, frame, prior, exact = FALSE,
                                    route = c("map", "sample")) {
  route <- match.arg(route)
  # the joint density of a rescor model is not a product over the
  # responses, and frame() refuses weights() with rescor = TRUE
  if (isTRUE(spec$rescor) ||
        length(frame[["autocor"]] %||% list()) ||
        length(frame[["mi_map"]] %||% list())) {
    return(invisible())
  }
  multi <- length(spec$responses) > 1L
  for (r in names(spec$responses)) {
    hit <- weights_unbounded_response(spec, frame, prior, r, exact, route)
    if (!is.null(hit)) {
      # a prior the fit would refuse is the first thing to fix, and the
      # refusal below would hide its error (frm() reaches this check
      # before it resolves the prior); frm_sample() resolved its own
      if (route == "map" && length(prior)) {
        resolve_prior_input(list(frame = frame, spec = spec), prior)
      }
      hit$resp <- if (multi) r
      frm_stop(do.call(weights_unbounded_message, hit), call. = FALSE)
    }
  }
  invisible()
}

#' The weights check on one response: `NULL` when no unbounded direction
#' is proved, or the arguments of its refusal's message.
#'
#' @noRd
weights_unbounded_response <- function(spec, frame, prior, r, exact,
                                       route) {
  resp <- spec$responses[[r]]
  fam <- resp$family[["family"]] %||% ""
  spike <- unbounded_spike_dpar[[fam]]
  if (is.null(spike)) return(NULL)
  laplace_linear <- fam %in% c("asym_laplace", "zero_inflated_asym_laplace") &&
    !isTRUE(exact)
  if (!is.null(resp$aterms[["se"]])) return(NULL)
  lps <- lapply(spike, function(dp) frame[["linpreds"]][[linpred_key(r, dp)]])
  if (any(vapply(lps, is.null, NA))) return(NULL)
  lp <- frame[["linpreds"]][[linpred_key(r, "mu")]]
  Z <- lp[["Z"]]
  if (is.null(Z) || !is.null(lp[["nl_body"]]) || length(lp[["Z_mo"]])) {
    return(NULL)
  }
  # the response's own rows, which subset() can make fewer than the data's
  n <- nrow(Z)
  av <- frame[["aterm_values"]][[r]]
  w <- av[["weights"]] %||% rep(1, n)
  owned_blocks <- Filter(cs_full_rank, frame[["re_blocks"]])
  cols <- unlist(lapply(owned_blocks, `[[`, "c_idx"))
  cols <- intersect(cols, seq_len(ncol(Z)))
  if (!length(cols)) return(NULL)
  Zc <- Matrix::drop0(methods::as(Z[, cols, drop = FALSE], "CsparseMatrix"))
  one <- which(diff(Zc@p) == 1L)
  owned <- unique(Zc@i[Zc@p[one] + 1L] + 1L)
  pos <- which(w > 0)
  if (!length(pos) || !all(pos %in% owned)) return(NULL)
  # the excess is the uncensored rows' (see above); a censored row needs
  # its own effect too, or its probability can fall to 0 as the
  # dispersion goes to its limit
  cen <- av[["cens"]]
  cen <- if (is.null(cen)) rep(0, n) else rep_len(cen, n)
  unc <- pos[cen[pos] == 0]
  # a row at a point mass reads log(zi) or log(hu), whatever its own
  # effect, so it neither grows nor curves
  mass <- unbounded_point_mass[[fam]]
  at_mass <- rep(FALSE, n)
  if (!is.null(mass)) {
    at_mass <- frame[["y"]][[r]] %in% mass
    unc <- unc[!at_mass[unc]]
  }
  # a constant mean can put one side's censored rows at probability 1
  # (see above); trunc() bounds the mean's reach, so there it cannot
  shift <- if (is.null(av[["trunc_lb"]]) && is.null(av[["trunc_ub"]])) {
    unbounded_shift_sides(prior, r, lp)
  } else c(up = FALSE, down = FALSE)
  # membership by mask: a cell's rows are read in time proportional to
  # their number, never to n
  is_pos <- w > 0
  is_unc <- logical(n)
  is_unc[unc] <- TRUE
  excess_of <- function(rows) {
    u <- rows[is_unc[rows]]
    e <- if (laplace_linear) sum(w[u]) else sum(w[u] - 1)
    cr <- rows[cen[rows] != 0]
    if (isTRUE(exact) || !length(cr)) return(list(e = e, u = u, cost = 0))
    free_side <- max(0, if (shift[["up"]]) sum(cen[cr] == 1),
                     if (shift[["down"]]) sum(cen[cr] == -1))
    cost <- length(cr) - free_side
    list(e = e - cost / 2, u = u, cost = cost)
  }
  multi <- length(spec$responses) > 1L
  hit <- function(rows, ex, coef = NULL) {
    terms <- unique(vapply(Filter(function(bk) {
      any(bk[["c_idx"]] %in% cols[one])
    }, owned_blocks), function(bk) {
      if (multi) block_term_in_response(bk, r, resp$resp_name %||% r) else
        block_term_as_written(bk)
    }, ""))
    list(terms = terms, w = w[ex$u], excess = ex$e, spike = spike,
         fam = fam, laplace_linear = laplace_linear, route = route,
         n_cens = sum(cen[rows] != 0), cost = ex$cost, coef = coef,
         n_mass = sum(at_mass[rows] & cen[rows] == 0), mass = mass)
  }
  rate <- unbounded_spike_rate[[fam]][spike] %||% rep(1, length(spike))
  rate[is.na(rate)] <- 1
  # rounding of weights normalized to mean 1 must not read as an excess
  tol <- function(ex) ex$e <= 1e-8 * max(1, length(ex$u))
  # every row at once, through each spike dpar's constant
  if (all(vapply(lps, unbounded_spans_constant, NA))) {
    ex <- excess_of(pos)
    if (!tol(ex)) {
      verdicts <- vapply(seq_along(spike), function(k) {
        prior_spike_verdict(prior, r, spike[k], lps[[k]], ex$e * rate[k],
                            route)
      }, "")
      # one bounded spike dpar keeps the densities from narrowing at all
      if (all(verdicts == "free")) return(hit(pos, ex))
    }
  }
  # the rows of one cell, through a coefficient that moves them alone;
  # with two spike dpars both would have to reach their limit on the
  # same rows, which the check does not try to prove
  if (length(spike) == 1L && unbounded_spike_free_lp(lps[[1L]])) {
    cells <- unbounded_spike_cells(lps[[1L]])
    # the design once, not once per cell: per cell it cost O(L n p)
    Xs <- if (length(cells) && length(prior)) as.matrix(lps[[1L]][["X"]])
    icpt_s <- if (!is.null(Xs)) unbounded_has_icpt(Xs)
    for (cell in cells) {
      rows <- cell$rows[is_pos[cell$rows]]
      ex <- excess_of(rows)
      if (tol(ex)) next
      v <- prior_column_verdict(prior, r, spike, lps[[1L]], cell$col,
                                ex$e * rate[1L], route, X = Xs,
                                icpt = icpt_s)
      if (identical(v, "free")) return(hit(rows, ex, coef = cell$name))
    }
    # no unit can gain while every uncensored weight is at most 1; an
    # unweighted model stops here
    if (laplace_linear || any(w[is_unc] > 1)) {
      pr <- unbounded_cell_pair(cells, lps[[1L]], prior, r, spike, is_pos,
                                is_unc, w, excess_of)
      if (!is.null(pr)) {
        out <- hit(pr$rows, pr$ex)
        out$pair <- pr$pair
        return(out)
      }
    }
  }
  NULL
}

#' Two units of a scale's design whose rows can move apart with the
#' centered intercept held, which a prior holds: the scale shrinks on
#' one unit's rows and grows on another's while that prior stays
#' constant.
#'
#' frmtmb puts a class `Intercept` prior on `Ic = b0 + means_X'b`, the
#' mean of the predictor over the design's rows. With `Ic` held, unit
#' j's log scale can move by `s_j` and unit k's by `s_k` when `n_j s_j +
#' n_k s_k = 0` (`n` counting the design's rows of each unit, as
#' `means_X` does). The units are of two kinds:
#'
#' - When the parametric design is saturated in its rows (its distinct
#'   rows number `rank(X)`, as for one factor in any coding, treatment,
#'   `contr.sum` or Helmert, and for a full interaction of factors), each
#'   distinct row is a unit. Any move of the units is then a direction in
#'   the coefficients, `d = solve(D, s)` over the distinct rows `D`.
#' - Otherwise, the rows that one column moves alone (a treatment
#'   column, `unbounded_spike_cells()`), and, when those cells do not
#'   overlap, the rows that none of them reaches (the reference level: b0
#'   moves and every cell's column takes it back). A numeric column
#'   beside a factor, or two factors without their interaction, give this
#'   kind; their other directions are not searched.
#'
#' As unit j's scale shrinks its rows gain `E_j` per unit of `-log
#' sigma` (the excess, less the censored rows' cost); as unit k's scale
#' grows, each uncensored row off a point mass loses its weight per unit,
#' for every scale in the table (its density at a fixed point falls like
#' 1 / sigma, and its own effect's curvature goes to 0), and a censored
#' or point-mass row loses nothing. So the objective grows without bound
#' along that ray when `E_j / n_j > W_k / n_k`, with `W_k` the weight of
#' unit k's uncensored rows off the mass. For a saturated design whose
#' columns are all free this is also the only way: over moves of the
#' units with `sum n_u s_u = 0`, the growth is linear and its maximum
#' puts all its mass on one shrinking unit and one growing unit.
#'
#' Decided only for a scale spike (sigma) on the log link, the one loss
#' rate proved here, and only for a pair whose moving columns all have a
#' free prior (none, flat or student_t, with no bound); the intercept's
#' column is free, since its prior is on the held `Ic`. A column's own
#' prior would cost a power of t that the linear growth must beat, and a
#' held column restricts the moves to a linear program that this does
#' not solve. Otherwise `NULL` (not decided). A prior that does not bound
#' `Ic` is the cell rule's case, decided before this one.
#'
#' @noRd
unbounded_cell_pair <- function(cells, lp, prior, r, spike, is_pos, is_unc,
                                w, excess_of) {
  if (!(spike %in% unbounded_spike_toward_zero)) return(NULL)
  if (!identical(lp[["link"]][["name"]] %||% "", "log")) return(NULL)
  X <- as.matrix(lp[["X"]])
  if (isFALSE(lp[["center"]]) || !unbounded_has_icpt(X)) return(NULL)
  # means_X holds the parametric columns only (lp_center_offset()); the
  # rest are held
  np <- min(lp[["n_param_cols"]] %||% ncol(X), ncol(X))
  Xp <- X[, seq_len(np), drop = FALSE]
  n <- nrow(Xp)
  cns <- (colnames(X) %||% character(ncol(X)))[seq_len(np)]
  icpt <- apply(Xp, 2L, function(v) v[1L] != 0 && all(v == v[1L]))
  if (!any(icpt)) return(NULL)
  # whether the class b prior of one column is free both ways: an
  # exponential or gamma tail costs linearly too, which is not weighed
  # here
  mine_b <- if (length(prior)) {
    if (!inherits(prior, "frmtmb_priorlist")) return(NULL)
    Filter(function(s) {
      identical(s$dpar %||% "", spike) && prior_resp_matches(s, r) &&
        identical(s$class, "b")
    }, prior)
  } else list()
  col_free <- vapply(seq_len(np), function(j) {
    if (icpt[j]) return(TRUE)
    hit <- Filter(function(s) identical(s$coef %||% "", cns[j]), mine_b)
    if (!length(hit)) hit <- Filter(function(s) !nzchar(s$coef %||% ""),
                                    mine_b)
    if (!length(hit)) return(TRUE)
    s <- hit[[1L]]
    kind <- if (is.null(s$dist)) "flat" else s$dist$kind %||% ""
    kind %in% c("t", "flat") &&
      !isTRUE(is.finite(suppressWarnings(as.numeric(s$lb)))) &&
      !isTRUE(is.finite(suppressWarnings(as.numeric(s$ub))))
  }, NA)
  cells <- Filter(function(cl) cl$col <= np, cells)
  cell_cols <- vapply(cells, `[[`, 0L, "col")
  cell_rows <- unlist(lapply(cells, `[[`, "rows"))
  disjoint <- !anyDuplicated(cell_rows)
  pat <- unbounded_row_patterns(Xp)
  if (!is.null(pat)) {
    rows_of <- split(seq_len(n), pat$id)
    # a column's direction for a pair: d = solve(D, s), read only when a
    # column is not free
    moves_of <- function(j, k) {
      if (all(col_free)) return(TRUE)
      if (nrow(pat$D) != np) return(FALSE)
      s <- numeric(nrow(pat$D))
      s[j] <- -1 / length(rows_of[[j]])
      s[k] <- 1 / length(rows_of[[k]])
      d <- solve(pat$D, s)
      mv <- abs(d) > 1e-9 * max(abs(d))
      all(col_free[mv])
    }
  } else {
    if (!length(cells)) return(NULL)
    rows_of <- lapply(cells, `[[`, "rows")
    moves <- as.list(cell_cols)
    ref <- if (disjoint) which(!seq_len(n) %in% cell_rows) else integer()
    if (length(ref)) {
      # b0 moves the reference rows, and every cell's column takes it back
      rows_of[[length(rows_of) + 1L]] <- ref
      moves[[length(moves) + 1L]] <- cell_cols
    }
    moves_of <- function(j, k) {
      all(col_free[union(moves[[j]], moves[[k]])]) &&
        !anyDuplicated(c(rows_of[[j]], rows_of[[k]]))
    }
  }
  size <- vapply(rows_of, length, 0L)
  exs <- lapply(rows_of, function(rr) excess_of(rr[is_pos[rr]]))
  gain <- vapply(exs, `[[`, 0, "e") / size
  loss <- vapply(rows_of, function(rr) sum(w[rr[is_unc[rr]]]), 0) / size
  tol <- function(g, j) g > 1e-8 * max(1, gain[j])
  if (!tol(max(gain) - min(loss), which.max(gain))) return(NULL)
  best <- NULL
  lo <- order(loss)
  for (j in order(gain, decreasing = TRUE)) {
    if (!is.null(best) && gain[j] - loss[lo[1L]] <= best$g) break
    for (k in lo) {
      g <- gain[j] - loss[k]
      if (!tol(g, j)) break
      if (j == k || !moves_of(j, k)) next
      if (is.null(best) || g > best$g) best <- list(g = g, j = j, k = k)
      break
    }
  }
  if (is.null(best)) return(NULL)
  j <- best$j
  k <- best$k
  # how the message names a unit: a treatment column's cell, the
  # reference level, or the values of its design row
  name_of <- function(u, down) {
    rr <- rows_of[[u]]
    if (!is.null(pat)) {
      row <- pat$D[u, ]
      nz <- which(row != 0 & !icpt)
      cl <- match(nz, cell_cols)
      # the pattern is a column's cell only when it holds all its rows,
      # and the reference level only when it holds every row no cell
      # reaches
      if (length(nz) == 1L && !is.na(cl) &&
            length(rr) == length(cells[[cl]]$rows)) {
        u_cell <- cells[[cl]]$name
      } else if (!length(nz) && disjoint && length(cells) &&
                   length(rr) == n - length(cell_rows)) {
        u_cell <- NA_character_
      } else {
        return(paste0(" rows where ",
                      paste0("`", cns[nz], "` is ", format(row[nz]),
                             collapse = " and "),
                      " in the design of ", spike))
      }
    } else {
      u_cell <- if (u <= length(cells)) cells[[u]]$name else NA_character_
    }
    if (is.na(u_cell)) {
      return(paste0(" rows that none of ",
                    paste0("`", cns[cell_cols], "`", collapse = ", "),
                    " of ", spike, " reaches"))
    }
    if (down) {
      paste0(" rows that the coefficient `", u_cell, "` of ", spike,
             " moves alone")
    } else paste0(" rows that `", u_cell, "` moves alone")
  }
  list(rows = rows_of[[j]][is_pos[rows_of[[j]]]], ex = exs[[j]],
       pair = list(down = name_of(j, TRUE), up = name_of(k, FALSE),
                   gain = gain[j], loss = loss[k], n_up = size[k],
                   n_down = size[j]))
}

#' The distinct rows of a design and each row's pattern, when the design
#' is saturated in them: as many distinct rows as its rank. Then every
#' pattern can move alone. `NULL` otherwise, and as soon as the distinct
#' rows outnumber the columns. The rows are coded column by column, so
#' the cost is n times the columns, with no row turned into a string.
#'
#' @noRd
unbounded_row_patterns <- function(Xp) {
  p <- ncol(Xp)
  id <- rep(1L, nrow(Xp))
  for (j in seq_len(p)) {
    v <- Xp[, j]
    code <- match(v, unique(v))
    key <- (id - 1) * max(code) + code
    id <- match(key, unique(key))
    if (max(id) > p) return(NULL)
  }
  U <- max(id)
  if (U < 2L) return(NULL)
  D <- Xp[match(seq_len(U), id), , drop = FALSE]
  if (qr(D)$rank != U) return(NULL)
  list(id = id, D = D)
}

#' Whether a dpar's coefficients are free to move at all: no constant
#' value and no non-linear body.
#'
#' @noRd
unbounded_spike_free_lp <- function(lp) {
  !is.null(lp) && is.null(lp[["constant"]]) && is.null(lp[["nl_body"]]) &&
    !is.null(lp[["X"]]) && ncol(lp[["X"]]) > 0L
}

#' The columns of a dispersion design that move the rows of one cell
#' alone: a column with one nonzero value, on some rows but not all. As
#' its coefficient goes to its limit, those rows' dpar goes there and
#' every other row's stays where it is (`sigma ~ fac`, `sigma ~ 0 +
#' fac`).
#'
#' @noRd
unbounded_spike_cells <- function(lp) {
  X <- as.matrix(lp[["X"]])
  cns <- colnames(X) %||% character(ncol(X))
  out <- list()
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    nz <- which(v != 0)
    if (!length(nz) || length(nz) == nrow(X) || any(v[nz] != v[nz[1L]])) {
      next
    }
    out[[length(out) + 1L]] <- list(col = j, name = cns[j], rows = nz)
  }
  out
}

#' The verdict of the prior on one coefficient of a spike dpar, for the
#' rows its column moves alone. A design with an intercept names the
#' coefficient's prior class `b`; so does a design without one. No
#' entry is "free".
#'
#' frmtmb evaluates a class `Intercept` prior at the centered intercept
#' `b0 + means_X'b` (R/priors.R, "Where an intercept prior lands"), so
#' the column's coefficient going to its limit with `b0` held moves the
#' centered intercept too, and holding the centered intercept instead
#' moves every other row. When that prior is not free, whether either
#' ray is unbounded depends on the other cells' weights (with two cells,
#' on whether the mean excess of one exceeds the mean weight of the
#' other), which is not decided here: "unknown".
#'
#' @param X,icpt The design as a dense matrix and whether it has an
#'   intercept, which a caller looping over columns reads once.
#' @noRd
prior_column_verdict <- function(prior, r, dpar, lp, j, excess, route,
                                 X = as.matrix(lp[["X"]]),
                                 icpt = unbounded_has_icpt(X)) {
  if (!length(prior)) return("free")
  if (!inherits(prior, "frmtmb_priorlist")) return("unknown")
  zero <- dpar %in% unbounded_spike_toward_zero
  link <- lp[["link"]][["name"]] %||% ""
  mine <- Filter(function(s) {
    identical(s$dpar %||% "", dpar) && prior_resp_matches(s, r)
  }, prior)
  if (!isFALSE(lp[["center"]]) && icpt) {
    for (s in mine) {
      if (!identical(s$class, "Intercept") || nzchar(s$coef %||% "")) next
      if (!identical(prior_entry_verdict(s, zero, link, excess, route),
                     "free")) {
        return("unknown")
      }
    }
  }
  cn <- (colnames(X) %||% character(ncol(X)))[j]
  prior_b_column_verdict(Filter(function(s) identical(s$class, "b"), mine),
                         cn, zero, link, excess, route)
}

#' The verdict of the class `b` entries that reach one column: the
#' column's own coefficient entry, else a class-wide one; none is
#' "free". A coefficient's prior is on the link scale.
#'
#' @noRd
prior_b_column_verdict <- function(b, cn, zero, link, excess, route) {
  hit <- Filter(function(s) identical(s$coef %||% "", cn), b)
  if (!length(hit)) hit <- Filter(function(s) !nzchar(s$coef %||% ""), b)
  if (!length(hit)) return("free")
  s <- hit[[1L]]
  s$natural <- FALSE
  prior_entry_verdict(s, zero, link, excess, route)
}

#' Whether a prior entry addresses response `r`: no `resp`, the raw
#' name, or brms's spelling of it (`make_stan_names()`: `y_a` is "ya"),
#' which frmtmb accepts too.
#'
#' @noRd
prior_resp_matches <- function(s, r) {
  rs <- s$resp %||% ""
  !nzchar(rs) || rs %in% c(r, brms_stan_name(r))
}

#' Whether a design has an intercept column: constant and nonzero.
#'
#' @noRd
unbounded_has_icpt <- function(X) {
  any(apply(X, 2L, function(v) v[1L] != 0 && all(v == v[1L])))
}

#' The directions in which a constant mean can move without bound: up
#' (to put every right-censored row's probability at 1) and down (every
#' left-censored row's). It needs a mean design that spans the constant,
#' and a prior on whatever moves (the intercept, or every column of a
#' design without one) whose support is unbounded on that side. Any
#' full-support density costs only a constant along the way; a bound or
#' a density on a half-line or an interval stops it.
#'
#' @noRd
unbounded_shift_sides <- function(prior, r, lp) {
  none <- c(up = FALSE, down = FALSE)
  if (!unbounded_spans_constant(lp)) return(none)
  if (!length(prior)) return(c(up = TRUE, down = TRUE))
  if (!inherits(prior, "frmtmb_priorlist")) return(none)
  X <- as.matrix(lp[["X"]])
  mine <- Filter(function(s) {
    (s$dpar %||% "") %in% c("", "mu") && prior_resp_matches(s, r)
  }, prior)
  moved <- if (!isFALSE(lp[["center"]]) && unbounded_has_icpt(X)) {
    Filter(function(s) {
      identical(s$class, "Intercept") && !nzchar(s$coef %||% "")
    }, mine)
  } else {
    # the constant moves every column; a column without an entry of its
    # own is reached by a class-wide one
    Filter(function(s) identical(s$class, "b"), mine)
  }
  sides <- c(up = TRUE, down = TRUE)
  for (s in moved) sides <- sides & prior_support_sides(s)
  sides
}

#' Whether one prior entry's support reaches +Inf (up) and -Inf (down)
#' on the scale its parameter lives on. An unknown density reaches
#' neither.
#'
#' @noRd
prior_support_sides <- function(s) {
  # an entry with bounds and no density is flat between them
  kind <- if (is.null(s$dist)) "flat" else s$dist$kind %||% ""
  lb <- suppressWarnings(as.numeric(s$lb))
  ub <- suppressWarnings(as.numeric(s$ub))
  out <- switch(kind,
    normal = , t = , logistic = , flat = c(up = TRUE, down = TRUE),
    exponential = , gamma = , inv_gamma = , lognormal =
      c(up = TRUE, down = FALSE),
    c(up = FALSE, down = FALSE))
  if (isTRUE(is.finite(ub))) out[["up"]] <- FALSE
  if (isTRUE(is.finite(lb))) out[["down"]] <- FALSE
  out
}

#' A block's term as the formula spells it. block_written_label() falls
#' back to "tt + 0 | g [ar1]" for a covariance-structure wrapper, which
#' is not what the user wrote.
#'
#' @noRd
block_term_as_written <- function(bk) {
  wl <- block_written_label(bk)
  m <- regmatches(wl, regexec("^(.*) \\[([A-Za-z_0-9.]+)\\]$", wl))[[1]]
  if (length(m) != 3L) return(wl)
  # a dpar's block is labeled "sigma: 1 | g"; the prefix stays outside
  # the term, as the formula's sigma ~ (1 | g) has it
  pre <- regmatches(m[2], regexec("^([A-Za-z_.0-9]+): (.*)$", m[2]))[[1]]
  dp <- if (length(pre) == 3L) paste0(pre[2], ": ") else ""
  body <- if (length(pre) == 3L) pre[3] else m[2]
  # a plain bar term is stored as a us() block, and a Student-t one,
  # written with gr(dist = "student"), as us_t() or diag_t()
  if (m[3] %in% c("us_t", "diag_t")) {
    lr <- strsplit(body, " | ", fixed = TRUE)[[1]]
    if (length(lr) == 2L) {
      return(paste0(dp, "(", lr[1], if (m[3] == "diag_t") " || " else " | ",
                    "gr(", lr[2], ", dist = \"student\"))"))
    }
  }
  paste0(dp, if (identical(m[3], "us")) paste0("(", body, ")") else
    paste0(m[3], "(", body, ")"))
}

#' A block's term as one response of a multivariate model writes it. The
#' frame labels the block "y2 1 | obs", and a block shared across
#' responses through an ID joins every response's label; the response's
#' own component, less its name, is what its formula holds.
#'
#' @noRd
block_term_in_response <- function(bk, r, resp_name) {
  pre <- paste0(resp_name, " ")
  cps <- bk[["components"]] %||% list()
  key <- linpred_key(r, "mu")
  cp <- Filter(function(cp) identical(cp[["lp_key"]], key), cps)
  if (length(cps) > 1L && length(cp) && identical(bk[["covstruct"]], "us")) {
    lab <- cp[[1L]][["label"]] %||% ""
    if (startsWith(lab, pre)) lab <- substring(lab, nchar(pre) + 1L)
    return(paste0("(", lab, ")"))
  }
  # "(y2 1 | obs)", "ar1(y2 tt + 0 | g)"
  sub(paste0("(", pre), "(", block_term_as_written(bk), fixed = TRUE)
}

#' Whether a dpar's predictor can reach its limit at every row at once:
#' its fixed-effect design has a column and spans the constant. An
#' intercept-free `sigma ~ 0 + x` and an offset-only `sigma ~ 0 +
#' offset(o)` cannot, and their likelihood is bounded.
#'
#' @noRd
unbounded_spans_constant <- function(lp) {
  if (is.null(lp) || !is.null(lp[["constant"]]) ||
        !is.null(lp[["nl_body"]])) {
    return(FALSE)
  }
  X <- lp[["X"]]
  if (is.null(X) || !ncol(X)) return(FALSE)
  X <- as.matrix(X)
  one <- rep(1, nrow(X))
  res <- qr.resid(qr(X), one)
  sqrt(sum(res^2)) <= 1e-8 * sqrt(nrow(X))
}

#' What the prior on one spike dpar does to the growth: "bounds" when it
#' provably stops it (MAP) or makes the posterior integrate (sample),
#' "free" when it provably does not, "unknown" otherwise. No entry on
#' the dpar is "free": the objective is then the likelihood (MAP) or the
#' sampler's flat density on the link scale (sample).
#'
#' @noRd
prior_spike_verdict <- function(prior, r, dpar, lp, excess, route) {
  if (!length(prior)) return("free")
  if (!inherits(prior, "frmtmb_priorlist")) return("unknown")
  zero <- dpar %in% unbounded_spike_toward_zero
  link <- lp[["link"]][["name"]] %||% ""
  mine <- Filter(function(s) {
    identical(s$dpar %||% "", dpar) && prior_resp_matches(s, r)
  }, prior)
  X <- as.matrix(lp[["X"]])
  if (!unbounded_has_icpt(X)) {
    return(prior_cellmeans_verdict(mine, X, zero, link, excess, route))
  }
  # bf(center = FALSE) or `0 + Intercept`: the intercept is a class b
  # coefficient, and only its own column moves
  if (isFALSE(lp[["center"]])) {
    ic <- which(apply(X, 2L, function(v) v[1L] != 0 && all(v == v[1L])))
    cn <- (colnames(X) %||% character(ncol(X)))[ic[1L]]
    return(prior_b_column_verdict(
      Filter(function(s) identical(s$class, "b"), mine), cn, zero, link,
      excess, route))
  }
  out <- character(0)
  for (s in mine) {
    if (!identical(s$class, "Intercept") || nzchar(s$coef %||% "")) next
    out <- c(out, prior_entry_verdict(s, zero, link, excess, route))
  }
  if (!length(out)) return("free")
  if (any(out == "bounds")) return("bounds")
  if (any(out == "unknown")) return("unknown")
  "free"
}

#' The verdict for a dispersion design that spans the constant without
#' an intercept column, as `sigma ~ 0 + fac` does through its cell
#' means. Every column then has to move to reach the limit at every
#' row, so a class `b` prior that bounds each column's coefficient on
#' the link scale bounds the model as an Intercept prior does. Every
#' column free, or no `b` entry at all, is "free"; some columns bounded
#' and some not is "unknown" here, and a free column whose own rows run
#' away is refused through `unbounded_spike_cells()`.
#'
#' @noRd
prior_cellmeans_verdict <- function(mine, X, zero, link, excess, route) {
  b <- Filter(function(s) identical(s$class, "b"), mine)
  if (!length(b)) return("free")
  cns <- colnames(X) %||% character(ncol(X))
  per_col <- vapply(cns, function(cn) {
    hit <- Filter(function(s) identical(s$coef %||% "", cn), b)
    if (!length(hit)) hit <- Filter(function(s) !nzchar(s$coef %||% ""), b)
    if (!length(hit)) return("none")
    # a coefficient's prior is on the link scale
    s <- hit[[1L]]
    s$natural <- FALSE
    prior_entry_verdict(s, zero, link, excess, route)
  }, "")
  if (all(per_col == "bounds")) return("bounds")
  # every column free moves the constant freely too
  if (all(per_col %in% c("free", "none"))) return("free")
  "unknown"
}

#' The power j of x in the log-Jacobian `log |dx / d eta|` that frmtmb
#' adds to a natural-scale prior, as x goes to 0 (`zero`) or to
#' infinity, for a link with inverse x(eta). log: dx/deta = x, j = 1 at
#' both ends. identity: 1, j = 0. softplus, x = log(1 + exp(eta)):
#' dx/deta = 1 - exp(-x), which is x near 0 (j = 1) and 1 at infinity
#' (j = 0). squareplus, x = (eta + sqrt(eta^2 + 4)) / 2: dx/deta =
#' x / sqrt(eta^2 + 4), where x is about 1/|eta| near 0 (j = 2) and
#' about eta at infinity (j = 0). logm1 (frechet's and student's nu),
#' x = 1 + exp(eta): dx/deta = x - 1, j = 1 at infinity; x never goes
#' to 0 there. Any other link: NA, and the verdict is "unknown".
#'
#' This is the MAP objective's power only. A posterior's integrability
#' does not depend on the link: with the Jacobian, the integral over eta
#' is the integral of p(x) L(x) dx, so the sample route reads j = 1
#' under every link (`prior_entry_verdict()`).
#'
#' @noRd
natural_jacobian_power <- function(link, zero) {
  switch(link %||% "",
    log = 1,
    identity = 0,
    softplus = if (zero) 1 else 0,
    squareplus = if (zero) 2 else 0,
    logm1 = if (zero) NA_real_ else 1,
    NA_real_)
}

#' The verdict of one prior entry on a spike dpar's intercept.
#'
#' `excess` is the power the likelihood grows by in the dpar's own
#' units (the weight excess times `unbounded_spike_rate`): like
#' x^(-excess) as a scale x goes to 0, like x^excess as a shape x grows.
#'
#' A prior on the natural scale (`class = "sigma"`, `"shape"`) is
#' evaluated with the inverse link's log-Jacobian (R/priors.R,
#' `natural_log_jacobian()`), so frmtmb's MAP objective, on the link
#' scale, carries a further x^j whose power depends on the link and the
#' limit (`natural_jacobian_power()`; 1 for log). Toward 0 a density
#' like x^order gives the objective x^(order + j - excess): bounded
#' when order + j > excess. Toward infinity a density with a faster
#' than polynomial tail stops it, and a power tail x^-p (student_t: p =
#' nu + 1; inv_gamma: p = shape + 1) gives x^(excess + j - p): bounded
#' when p - j > excess. The posterior integrates when the same margins
#' are positive with j = 1, whatever the link: the integral over eta is
#' the integral of p(x) L(x) dx. At equality the MAP objective has a
#' finite limit and no maximum is proved either way ("unknown"), and the
#' posterior does not integrate ("free"): the integrand in x tends to a
#' positive constant times 1/x.
#'
#' A prior on a link-scale coefficient (class `Intercept` or `b` of a
#' dpar with a formula) carries no Jacobian. On the log link the
#' likelihood grows like exp(excess |eta|): a normal or bounded density
#' stops it, a student_t, inv_gamma or lognormal tail does not, and an
#' exponential tail exp(-r |eta|) stops it when r exceeds the excess.
#'
#' @noRd
prior_entry_verdict <- function(s, zero, link, excess, route) {
  d <- s$dist
  kind <- d$kind %||% ""
  lb <- suppressWarnings(as.numeric(s$lb))
  ub <- suppressWarnings(as.numeric(s$ub))
  # the objective's power at the limit: below 0 it is bounded, above 0
  # it runs away, and at 0 the next order decides, which is not read
  margin <- function(m) {
    if (!is.finite(m)) return("unknown")
    if (m > 0) "bounds" else if (m < 0 || route == "sample") "free" else
      "unknown"
  }
  if (isTRUE(s$natural)) {
    # the MAP objective's link power; the posterior's is 1 for any link
    jac <- if (identical(route, "sample")) 1 else
      natural_jacobian_power(link, zero)
    if (zero && isTRUE(lb > 0)) return("bounds")
    if (!zero && isTRUE(is.finite(ub))) return("bounds")
    if (!zero) {
      if (kind %in% c("normal", "gamma", "exponential", "logistic",
                      "lognormal", "beta", "uniform")) {
        return("bounds")
      }
      p <- switch(kind, t = d$df + 1, inv_gamma = d$shape + 1, flat = 0,
                  NA_real_)
      if (is.na(p) || is.na(jac)) return("unknown")
      return(margin(p - jac - excess))
    }
    # brms's (flat), frmtmb.sample's prior = "flat", is order 0 at 0
    order <- switch(kind,
      normal = , t = , logistic = , exponential = , flat = 0,
      uniform = if (isTRUE(d$lower > 0)) Inf else 0,
      gamma = d$shape - 1,
      beta = d$shape1 - 1,
      inv_gamma = , lognormal = Inf,
      NA_real_)
    if (is.na(order)) return("unknown")
    if (is.infinite(order)) return("bounds")
    if (is.na(jac)) return("unknown")
    margin(order + jac - excess)
  } else {
    if (!identical(link, "log")) return("unknown")
    if (zero && isTRUE(is.finite(lb))) return("bounds")
    if (!zero && isTRUE(is.finite(ub))) return("bounds")
    switch(kind,
      normal = "bounds",
      t = , flat = "free",
      uniform = "bounds",
      # a density on (0, Inf) never reaches eta = -Inf; toward +Inf the
      # log-normal tail exp(-(log eta)^2 / 2) is slower than exp(eta)
      exponential = , gamma = , inv_gamma = , lognormal = , beta =
        if (zero || kind == "beta") "bounds" else
          if (kind %in% c("inv_gamma", "lognormal")) "free" else
            margin(d$rate - excess),
      logistic = margin(1 / d$scale - excess),
      "unknown")
  }
}

#' The refusal's text, as the model's own terms and numbers.
#'
#' @param cost The censored rows that each take 1/2 off the power under
#'   the Laplace approximation.
#' @param coef The dispersion coefficient whose own rows run away, when
#'   the refusal is about one cell rather than every row.
#' @param resp The response, in a multivariate model.
#' @param pair From `unbounded_cell_pair()`: how the shrinking and the
#'   growing rows are named, with their gain, loss and row counts.
#' @param n_mass,mass The uncensored rows at the family's point mass, and
#'   the values it sits at (`unbounded_point_mass`).
#' @noRd
weights_unbounded_message <- function(terms, w, excess, spike, fam,
                                      laplace_linear, route, n_cens = 0L,
                                      cost = n_cens, coef = NULL,
                                      resp = NULL, pair = NULL,
                                      n_mass = 0L, mass = NULL) {
  zero <- spike[1L] %in% unbounded_spike_toward_zero
  dp <- paste(spike, collapse = " and ")
  lim <- if (zero) {
    paste(dp, if (length(spike) > 1L) "shrink" else "shrinks", "to 0")
  } else {
    paste(dp, "grows")
  }
  head <- if (route == "sample") {
    "The posterior of this model is improper. "
  } else if (laplace_linear) {
    "The Laplace approximation of this model's likelihood has no maximum. "
  } else {
    "The weighted likelihood of this model has no maximum. "
  }
  where <- if (is.null(resp)) "mu" else paste0("the mu of response ", resp)
  fmt <- function(x) format(x, digits = 4)
  rows <- if (!is.null(coef)) {
    paste0(" rows that the coefficient `", coef, "` of ", spike[1L],
           " moves alone")
  } else if (!is.null(pair)) pair$down else " rows"
  held <- if (!is.null(pair)) {
    paste0("A prior holds the Intercept of ", spike[1L], ", which frmtmb ",
           "places at the predictor means, but ", spike[1L], " can shrink ",
           "on these rows while it grows on the ", pair$n_up, pair$up, ", ",
           "and keep that intercept fixed. Each row there loses its weight ",
           "as ", spike[1L], " grows, ", fmt(pair$loss), " a row on ",
           "average, which is less than the excess a row here, ",
           fmt(pair$gain), ". ")
  }
  pw <- if (is.null(pair)) excess else excess - pair$n_down * pair$loss
  less <- paste0(if (cost && n_cens) paste0(" - ", cost, "/2"),
                 if (!is.null(pair)) {
                   paste0(" - ", pair$n_down, " x ", fmt(pair$loss))
                 })
  owns <- paste0(paste0("`", terms, "`", collapse = ", "), " gives every ",
                 "row with a positive weight its own random effect in ",
                 where,
                 if (!is.null(coef) || !is.null(pair) ||
                       (length(w) && any(w != 1))) {
                   paste0(", and the weights sum to ",
                          format(sum(w), digits = 6), " over ",
                          length(w), if (n_cens) " uncensored", rows)
                 }, ". ",
                 if (n_cens) {
                   paste0("The ", n_cens, " censored rows do not stop ",
                          "it: a censored row's probability is at most 1, ",
                          "and its own effect can bring it to 1",
                          if (cost) {
                            paste0("; under the Laplace approximation the ",
                                   cost, " that a constant mean cannot ",
                                   "bring to 1 with the rest take 1/2 ",
                                   "each off the power below")
                          }, ". ")
                 },
                 if (n_mass) {
                   paste0("The ", n_mass, " rows at ",
                          paste(mass, collapse = " or "), " sit in the ",
                          fam, " point mass, whose probability reads ",
                          "neither mu nor ", spike[1L], ", so they ",
                          "neither add to the growth nor stop it. ")
                 }, held)
  why <- if (laplace_linear) {
    paste0("The ", fam, " log density is linear on each side of its ",
           "mode, so the data give a row's own effect no curvature, and as ",
           lim, " the approximation grows like the peak height to the ",
           "power sum(w)", less, " = ", fmt(pw),
           ", whatever the weights. quadrature = TRUE integrates the ",
           "effect exactly and has a maximum when the weights sum to less ",
           "than the number of rows; so does the model without the ",
           "observation-level term.")
  } else {
    paste0("A weight above 1 counts a row more than once, so the row's ",
           "own effect then stands for replicates with no residual ",
           "variance: as ", lim, ", each density narrows onto its own ",
           "observation, and the likelihood grows like the height of ",
           "that peak to the power sum(w) - n", less, " = ",
           fmt(pw), ", without bound. A usable fit ",
           "needs weights that sum to the number of rows or fewer, a ",
           "random effect that is not observation-level, or a fixed ",
           spike[1L], ".")
  }
  # a prior strong enough to bound the objective puts the estimate where
  # it and the excess balance: lognormal(0, 1) on sigma stopped at a log
  # sigma of -19.8 and sampled with Rhat 2.83
  pri <- if (laplace_linear) "" else {
    paste0(" A prior on ", if (!is.null(coef)) {
             paste0("the coefficient `", coef, "`")
           } else if (!is.null(pair)) {
             paste0("the coefficients of ", spike[1L])
           } else spike[1L],
           " can make the objective bounded, but not the estimate ",
           "informative: the mode then sits where the prior and the ",
           "weight excess balance.")
  }
  paste0(head, owns, why, pri)
}
