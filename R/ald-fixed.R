# asym_laplace() and zero_inflated_asym_laplace() fits without random
# effects: the Hessian standard error refused by name, and the proof that
# the fit stops at the quantile regression's optimum.

#' The reason `se_lost` gives a parameter whose Hessian standard error is
#' refused by name (ald_se_lost()).
#'
#' @noRd
ald_se_reason <- "kink"

#' The routes to intervals that the refusal names, in this order. A
#' family with curvature that estimates the same quantile goes first
#' when there is one.
#'
#' frmtmb's own routes read the same working likelihood, so they are
#' named with that caveat: on gaussian data at quantile 0.9 (400 rows),
#' frm_bootstrap()'s slope standard error was 0.106 to 0.125 against a
#' Monte Carlo sd of rq()'s slope of 0.298 (lane squant,
#' dev/squant-boot-ald.R on rellib-r22; dev/alapse-findings.md for the
#' lane build).
#'
#' @noRd
ald_se_routes <- c(
  paste0("quantreg::rq() with se = \"nid\", or a bootstrap of the rows, ",
         "whose intervals do not assume the asymmetric Laplace"),
  paste0("frm_bootstrap() or frmtmb.sample::frm_sample() only where the ",
         "data are asymmetric Laplace: both read this likelihood, and on ",
         "other data the bootstrap's standard errors were less than half ",
         "the true spread"))

#' The responses of a fit whose Hessian standard errors are refused by
#' name: asym_laplace() or zero_inflated_asym_laplace() responses of a
#' fit with no random effect, no latent value and no integration. With
#' random effects the outer Hessian is that of a Laplace or quadrature
#' marginal, which the refusal of check_ald_laplace() and backlog item
#' 5449 are about.
#'
#' @noRd
ald_fixed_resps <- function(fit) {
  if (!length(fit$opt$par) || length(fit$obj$env$random) ||
        length(fit$frame[["re_blocks"]]) || !is.null(fit$importance) ||
        isTRUE(fit$quadrature)) {
    return(character(0))
  }
  fams <- vapply(fit$spec$responses, function(r) {
    r$family[["family"]] %||% ""
  }, "")
  rs <- names(fams)[fams %in% ald_laplace_families]
  # the proof is about the density of an exact row; a response whose
  # every row is censored reads only tail probabilities, which are
  # curved, and the usual analysis judges it (backlog item 5490)
  rs[vapply(rs, function(r) {
    av <- fit$frame[["aterm_values"]][[r]]
    y <- as.numeric(fit$frame[["y"]][[r]])
    ok <- is.finite(y) & (av[["cens"]] %||% 0) == 0 &
      (av[["weights"]] %||% 1) > 0
    if (fams[[r]] == "zero_inflated_asym_laplace") ok <- ok & y != 0
    any(ok)
  }, NA)]
}

#' The outer positions of a linear predictor's coefficients, NA for one
#' the fit holds fixed.
#'
#' @noRd
ald_lp_pos <- function(fit, lp, om) {
  cp <- lp[["par"]]
  idx <- lp[["idx"]]
  if (!length(idx) || !(cp %in% c("beta", "betad"))) return(integer(0))
  if (cp == "betad") {
    fx <- fit$frame[["betad_fixed_idx"]]
    keep <- setdiff(seq_along(fit$frame[["par_template"]][["betad"]]), fx)
    idx <- match(idx, keep)
  }
  which(om$comp == cp)[idx]
}

#' The outer parameters whose Hessian standard error is refused by name,
#' with reason `ald_se_reason`: every coefficient of an asym_laplace()
#' response but those of `zi`, whose rows the asymmetric Laplace does
#' not read. Empty when the fit has no such response.
#'
#' The log density is linear in mu on each side of its mode, so its
#' second derivative in mu is 0 wherever it is defined: the Hessian of
#' the location coefficients is exactly 0 (with the identity link) and
#' carries no information about their uncertainty, which depends on the
#' density of the data at the quantile (the sparsity), a quantity the
#' likelihood never reads. Another link adds `rho'(u) f''(eta)`, which is
#' not that either. The standard errors of sigma and the quantile need
#' the inverse of the whole matrix, so they go too. Measured on 0.82.0
#' (dev/selimit2-asym.R): every fit lost every parameter as "nocurv",
#' or kept standard errors of exactly 0 before that repair.
#'
#' @noRd
ald_se_lost <- function(fit) {
  none <- stats::setNames(character(0), character(0))
  rs <- ald_fixed_resps(fit)
  if (!length(rs)) return(none)
  om <- outer_par_map(fit)
  if (length(om$names) != length(fit$opt$par)) return(none)
  pos <- integer(0)
  for (lp in fit$frame[["linpreds"]]) {
    if (!(lp[["resp"]] %in% rs) || identical(lp[["dpar"]], "zi")) next
    pos <- c(pos, ald_lp_pos(fit, lp, om))
  }
  pos <- sort(unique(pos[!is.na(pos)]))
  stats::setNames(rep(ald_se_reason, length(pos)), om$names[pos])
}

#' The standard-error analysis of a fit with refused parameters, in the
#' form cov_from_hessian() returns: those parameters NaN in `V` and
#' `shown`, and lost with reason `ald_se_reason`. The other parameters
#' (a `zi` part, another response) keep the analysis of their own block
#' when the exact Hessian holds no cross term between the two sets, as a
#' zero-inflation probability that reads only the zero rows does not;
#' otherwise they are lost with them. NULL when nothing is refused.
#' Cached, since every reader of the covariance asks.
#'
#' @noRd
ald_se_analysis <- function(fit) {
  lost <- ald_se_lost(fit)
  if (!length(lost)) return(NULL)
  cached <- if (is.environment(fit$cache)) fit$cache$ald_se_an
  if (!is.null(cached) && identical(cached$at, fit$opt$par)) {
    return(cached$an)
  }
  p <- fit$opt$par
  np <- length(p)
  nm <- outer_par_names(fit)
  L <- match(names(lost), nm)
  K <- setdiff(seq_len(np), L)
  an <- NULL
  if (length(K)) {
    obj <- fit$obj
    saved <- obj_state_save(obj)
    H <- tryCatch(obj_he_cached(obj, p), error = function(e) NULL)
    obj_state_restore(obj, saved)
    if (is.matrix(H) && all(dim(H) == np) && all(is.finite(H))) {
      Hs <- (H + t(H)) / 2
      cross <- Hs[L, K, drop = FALSE]
      if (max(abs(cross)) <= 64 * .Machine$double.eps * max(abs(Hs), 1)) {
        # the refused block replaced by a unit one, which the analysis of
        # the rest then never reads
        H2 <- Hs
        H2[L, ] <- 0
        H2[, L] <- 0
        H2[cbind(L, L)] <- 1
        an <- tryCatch(cov_from_hessian(fit, H2, abs(H - t(H)) / 2,
                                        exact = TRUE),
                       error = keep_missing(NULL))
      }
    }
    if (is.null(an)) {
      lost <- stats::setNames(rep(ald_se_reason, np), nm)
      L <- seq_len(np)
    }
  }
  V <- if (is.null(an)) matrix(NaN, np, np) else an$V
  shown <- if (is.null(an)) V else an$shown
  V[L, ] <- NaN
  V[, L] <- NaN
  shown[L, ] <- NaN
  shown[, L] <- NaN
  other <- if (is.null(an)) NULL else an$lost[!(names(an$lost) %in% nm[L])]
  null <- if (!is.null(an) && NCOL(an$null)) {
    N <- an$null
    N[L, ] <- 0
    N[, colSums(N != 0) > 0, drop = FALSE]
  }
  out <- list(V = V, shown = shown, lost = c(other, lost),
              null = if (NCOL(null)) null, flat = an$flat, tier = 3L,
              root = NULL)
  if (is.environment(fit$cache)) {
    fit$cache$ald_se_an <- list(at = fit$opt$par, an = out)
  }
  out
}

#' `sdr` with the refusal applied: whatever the sdreport and the repairs
#' before this made of the refused parameters, every reader of the
#' covariance (summary(), vcov(), confint(), predictions and the
#' adapters that read vcov()) gets NaN for them, the same reason, and
#' the other parameters' analysis, also where the stop was short of the
#' maximum (sdr_rescue_hessian() keeps sdreport()'s inverse there).
#'
#' @noRd
ald_sdr_refuse <- function(fit, sdr) {
  V0 <- sdr$cov.fixed
  if (!is.matrix(V0)) return(sdr)
  an <- ald_se_analysis(fit)
  if (is.null(an) || !all(dim(an$V) == dim(V0))) return(sdr)
  dn <- dimnames(V0)
  sdr$cov_fixed_raw <- sdr$cov_fixed_raw %||% V0
  sdr$cov.fixed <- an$shown
  dimnames(sdr$cov.fixed) <- dn
  sdr$cov_fixed_prop <- an$V
  dimnames(sdr$cov_fixed_prop) <- dn
  sdr$cov_fixed_root <- NULL
  sdr$se_lost <- an$lost
  sdr$se_null <- an$null
  sdr$se_null_flat <- an$flat
  sdr
}

#' The families a refusal names, as written.
#'
#' @noRd
ald_family_names <- function(fit) {
  rs <- ald_fixed_resps(fit)
  unique(vapply(fit$spec$responses[rs], function(r) {
    paste0(r$family[["family"]], "()")
  }, ""))
}

#' Why the refused parameters have no standard error, and where
#' intervals come from instead: the clause summary() and the fit-time
#' warning share.
#'
#' @noRd
ald_se_why <- function(fit, one) {
  fams <- paste(ald_family_names(fit), collapse = " and ")
  paste0(if (one) "it has" else "they have", " no Hessian standard ",
         "error, by design: the ", fams, " log density is linear in mu ",
         "on each side of its mode, so its curvature in mu is 0 almost ",
         "everywhere and the Hessian carries no information about the ",
         "uncertainty of the location coefficients, while the scale's and ",
         "the quantile's standard errors need the inverse of that same ",
         "Hessian. For intervals use ",
         paste(ald_se_routes, collapse = "; or "))
}

#' The fit-time refusal, one condition of class `frmtmb_ald_no_se`,
#' under `frmtmb_control(check_se = )`: a warning by default, an error
#' under "stop", nothing under "ignore".
#'
#' @noRd
ald_se_report <- function(fit, lost, act) {
  if (!length(lost) || identical(act, "ignore")) return(invisible(NULL))
  np <- length(fit$opt$par)
  cert <- fit$opt$ald_cert
  msg <- paste0("Standard errors are not available for ", length(lost),
                " of ", np, " parameters. ",
                paste(names(lost), collapse = ", "), ": ",
                ald_se_why(fit, length(lost) == 1L),
                if (isTRUE(cert$ok)) {
                  paste0(". The estimates are the exact quantile ",
                         "regression solution at the fitted scale and ",
                         "quantile, and ", if (isTRUE(cert$global)) {
                           "the highest maximum of the likelihood"
                         } else "a local maximum of the likelihood")
                },
                "; summary() and diagnose() list them")
  if (identical(act, "stop")) {
    frm_stop(msg, call. = FALSE, class = "frmtmb_ald_no_se")
  }
  frm_warning(msg, call. = FALSE, class = "frmtmb_ald_no_se")
}

#' The clause a convergence warning of a refused fit ends with, in place
#' of conv_headroom_clause(), which would promise the standard errors of
#' the stop: what the proof of ald_polish() found instead. NULL for any
#' other fit.
#'
#' @noRd
ald_conv_clause <- function(fit) {
  if (!length(ald_se_lost(fit))) return(NULL)
  why <- fit$opt$ald_cert$why %||% "it was not attempted"
  paste0(". The ", paste(ald_family_names(fit), collapse = " and "),
         " log density has kinks, so a Newton step cannot measure how far ",
         "the stop is from the maximum, and the proof that it is the ",
         "quantile regression's optimum failed: ", why)
}

# The proof, and the move to the exact solution --------------------------

#' The largest Newton gain over the smooth parameters at which the polish
#' stops: rounding, two orders below what check_convergence() accepts at
#' a maximum (conv_at_max_tol).
#'
#' @noRd
ald_gain_tol <- 1e-10

#' The rounds of ald_polish(): each solves the quantile regression at the
#' current scale and quantile, then takes one Newton step in them.
#' Inside one interval of the quantile regression's solution the steps
#' are Newton's on a smooth function, so a handful settles them.
#'
#' @noRd
ald_polish_rounds <- 30L

#' How far outside its interval a basic row's subgradient may lie and
#' still count as inside: rounding of the k by k solve, on a quantity of
#' order 1.
#'
#' @noRd
ald_margin_tol <- 1e-8

#' The data of each asym_laplace() response whose location the proof can
#' read, or a sentence saying why it cannot. It needs the location to be
#' the identity of a linear predictor (no nonlinear body, no mo()), every
#' row exact (no cens() or trunc()), no me() or residual correlation term,
#' no prior with a density and no bound on the location's coefficients,
#' and the design of full column rank. The other parameters only need a
#' smooth objective.
#'
#' @noRd
ald_lp_blocks <- function(fit, rs, bounds) {
  fr <- fit$frame
  # (an mi() value the fit estimates is a random effect, which
  # ald_fixed_resps() has excluded; one no density reads leaves its row,
  # as the objective's weight of 0 does)
  if (length(fr[["me"]] %||% list()) || length(fr[["autocor"]] %||% list())) {
    return("the model has me() or a residual correlation term")
  }
  om <- outer_par_map(fit)
  np <- length(fit$opt$par)
  if (length(om$names) != np) return("its parameters could not be mapped")
  ent <- tryCatch(fit_prior_entries(fit), error = keep_missing(NULL))
  dens <- Filter(function(e) !prior_entry_constant(e), ent %||% list())
  reads <- unlist(lapply(dens, function(e) {
    c(e[["comp"]], e[["offset"]][["comp"]])
  }))
  if (is.null(ent) || "beta" %in% reads) {
    return("a prior density on the coefficients")
  }
  lower <- rep_len(bounds$lower %||% -Inf, np)
  upper <- rep_len(bounds$upper %||% Inf, np)
  out <- list()
  for (r in rs) {
    lp <- fr[["linpreds"]][[linpred_key(r, "mu")]]
    if (is.null(lp) || !is.null(lp[["nl_body"]]) ||
          length(lp[["mo"]] %||% list()) || !is.null(lp[["Z"]]) ||
          !identical(lp[["link"]][["name"]] %||% "", "identity") ||
          !identical(lp[["par"]], "beta")) {
      return(paste0("the location of response ", r, " is not linear in ",
                    "its coefficients on the identity link"))
    }
    av <- fr[["aterm_values"]][[r]]
    if (any((av[["cens"]] %||% 0) != 0) || !is.null(av[["trunc_lb"]]) ||
          !is.null(av[["trunc_ub"]])) {
      return("a censored or truncated row reads a curved tail probability")
    }
    pos <- ald_lp_pos(fit, lp, om)
    X <- as.matrix(lp[["X"]])
    if (!length(pos) || anyNA(pos) || length(pos) != ncol(X)) {
      return("a location coefficient is held fixed")
    }
    if (any(is.finite(lower[pos]) | is.finite(upper[pos]))) {
      return("a bound holds a location coefficient")
    }
    y <- as.numeric(fr[["y"]][[r]])
    w <- as.numeric(av[["weights"]] %||% rep(1, length(y)))
    if (length(y) != nrow(X) || length(w) != nrow(X)) {
      return("the rows of the response and its design differ")
    }
    off <- lp[["offset"]] %||% 0
    w[fr[["mi_map"]][[r]][["unread"]]] <- 0
    rows <- is.finite(y) & is.finite(w) & w > 0
    if (identical(fit$spec$responses[[r]]$family[["family"]],
                  "zero_inflated_asym_laplace")) {
      rows <- rows & y != 0
    }
    rows <- which(rows)
    Xr <- X[rows, , drop = FALSE]
    if (!length(rows) || qr(Xr)$rank < ncol(X)) {
      return(paste0("the location design of response ", r, " is not of ",
                    "full rank on its rows"))
    }
    out[[r]] <- list(resp = r, X = Xr, y = (y - rep_len(off, length(y)))[rows],
                     w = w[rows], rows = rows, pos = pos)
  }
  out
}

#' The outer parameter vector `par` as the fit's estimates, as fit_assembled()'s
#' est_at() writes them (no integrated component: the fit has none).
#'
#' @noRd
ald_est_at <- function(fit, par) {
  est <- fit$estimates
  tpl <- fit$frame[["par_template"]]
  pn <- names(fit$opt$par)
  for (cp in unique(pn)) {
    pos <- seq_along(tpl[[cp]])
    if (cp == "betad" && length(fit$frame[["betad_fixed_idx"]])) {
      pos <- setdiff(pos, fit$frame[["betad_fixed_idx"]])
    }
    est[[cp]][pos] <- unname(par[pn == cp])
  }
  est
}

#' The exact minimizer of `sum(w * rho_p(y - X b))` (the check function
#' at a quantile per row), from the basis of the k rows nearest the fit
#' at `b0`, or from the basis `h0`, by the simplex steps of an L1 problem
#' (Barrodale and Roberts): from a vertex, the edge along which a basic
#' row leaves zero, a weighted-median line search to the next vertex.
#'
#' A vertex `h` is optimal exactly when the subgradient condition holds:
#' with `psi_i = w_i (p_i - [r_i < 0])` off the basis, the `s` solving
#' `X_h' (w_h s) = -sum psi_i x_i` lies in `[p_i - 1, p_i]` row by row
#' (Koenker 2005, Theorem 2.1, with weights and a quantile per row). A
#' row off the basis with a zero residual takes `p_i`, one point of its
#' interval, so a pass is still a proof; at a degenerate vertex (tied
#' data) that choice can fail where another would pass, and the steps
#' can then cycle, which ald_vertex_ties() resolves. Returns `b`, the
#' basis `h`, `margin` (the least distance of an `s` inside its
#' interval; negative outside, 0 for a solution that is not unique) and
#' `ok`; NULL when no vertex was reached.
#'
#' @noRd
ald_vertex <- function(X, y, w, p, b0, h0 = NULL, iter = NULL) {
  n <- nrow(X)
  k <- ncol(X)
  if (n < k) return(NULL)
  h <- h0
  if (length(h) != k) {
    ord <- order(abs(y - drop(X %*% b0)) * w)
    m <- min(n, 4L * k)
    repeat {
      cand <- ord[seq_len(m)]
      q <- qr(t(X[cand, , drop = FALSE]))
      if (q$rank == k || m == n) break
      m <- min(n, 2L * m)
    }
    if (q$rank < k) return(NULL)
    h <- cand[q$pivot[seq_len(k)]]
  }
  iter <- iter %||% (50L + 10L * k)
  for (it in seq_len(iter)) {
    Xh <- X[h, , drop = FALSE]
    b <- tryCatch(solve(Xh, y[h]), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    r <- y - drop(X %*% b)
    r[h] <- 0
    psi <- w * (p - (r < 0))
    psi[h] <- 0
    g <- drop(crossprod(X, psi))
    u <- tryCatch(solve(t(Xh), -g), error = function(e) NULL)
    if (is.null(u)) return(NULL)
    s <- u / w[h]
    lo <- p[h] - 1
    hi <- p[h]
    margin <- min(hi - s, s - lo)
    if (margin >= -ald_margin_tol) {
      return(list(b = b, h = h, margin = margin, ok = TRUE, iter = it))
    }
    # the leaving row: the steepest descent of the objective along its
    # edge, per unit move of its residual
    dn_hi <- w[h] * (hi - s)
    dn_lo <- w[h] * (s - lo)
    slope <- pmin(dn_hi, dn_lo)
    j <- which.min(slope)
    # its residual leaves zero upward when s is above its interval
    delta <- if (dn_hi[j] < dn_lo[j]) -1 else 1
    e <- numeric(k)
    e[j] <- delta
    d <- tryCatch(solve(Xh, e), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    xd <- drop(X %*% d)
    xd[h] <- 0
    # rows whose residual crosses zero as t grows: a row at zero counts
    # as positive (psi above), so it crosses when it moves down
    cross <- (r > 0 & xd > 0) | (r < 0 & xd < 0) | (r == 0 & xd > 0)
    cross[h] <- FALSE
    ci <- which(cross)
    if (!length(ci)) return(NULL)
    t_i <- pmax(r[ci] / xd[ci], 0)
    o <- order(t_i)
    sl <- slope[j] + cumsum(w[ci][o] * abs(xd[ci][o]))
    hit <- which(sl >= 0)
    if (!length(hit)) return(NULL)
    h[j] <- ci[o][hit[1L]]
  }
  NULL
}

#' The size of the perturbation that breaks ties (ald_vertex_ties()),
#' relative to `1 + max|y|`: far above the rounding of a residual, far
#' below the spacing of rounded data.
#'
#' @noRd
ald_tie_eps <- 1e-9

#' ald_vertex() for tied data. Rows with equal `y` at equal `x` make
#' vertices degenerate: several rows off the basis have a zero residual,
#' the steps can pivot with a step of 0 and cycle, and the test that
#' gives each such row the subgradient `p_i` can fail at an optimal
#' vertex (on the review's integer-y fits, 17 of 18 kept their warnings,
#' 5 of them at rq()'s own optimum; dev/reviews/2026-11-17-alapse.md,
#' finding 2).
#'
#' The search runs on `y` perturbed by a fixed pattern of size
#' `ald_tie_eps` (a symbolic perturbation in the sense of the simplex
#' method's anti-cycling rules), where no vertex is degenerate. Its
#' basis then gives `b` on the data as they are, and the proof is made
#' on those data: every row whose residual is not zero keeps its sign
#' (else NULL), and each row at zero off the basis takes the end of its
#' interval that its perturbed residual's sign picks. Those are
#' subgradients, so when the basic rows' `s` falls inside their
#' intervals, 0 is in the subdifferential and `b` is optimal: the
#' perturbation only finds the basis, it is not in the proof.
#'
#' @noRd
ald_vertex_ties <- function(X, y, w, p, b0) {
  n <- length(y)
  scale <- 1 + max(abs(y))
  pat <- ((seq_len(n) * 0.6180339887498949) %% 1) - 0.5
  yp <- y + ald_tie_eps * scale * pat
  v <- ald_vertex(X, yp, w, p, b0)
  if (is.null(v)) return(NULL)
  h <- v$h
  Xh <- X[h, , drop = FALSE]
  b <- tryCatch(solve(Xh, y[h]), error = function(e) NULL)
  if (is.null(b)) return(NULL)
  r <- y - drop(X %*% b)
  rp <- yp - drop(X %*% v$b)
  r[h] <- 0
  rp[h] <- 0
  zero <- abs(r) <= 1e-12 * scale
  if (any(!zero & sign(r) != sign(rp))) return(NULL)
  neg <- ifelse(zero, rp < 0, r < 0)
  psi <- w * (p - neg)
  psi[h] <- 0
  g <- drop(crossprod(X, psi))
  u <- tryCatch(solve(t(Xh), -g), error = function(e) NULL)
  if (is.null(u)) return(NULL)
  s <- u / w[h]
  margin <- min(p[h] - s, s - (p[h] - 1))
  if (margin < -ald_margin_tol) return(NULL)
  list(b = b, h = h, margin = margin, ok = TRUE, iter = v$iter, ties = TRUE)
}

#' ald_vertex() from the last basis, then from the rows nearest `b0`,
#' then with the ties broken.
#'
#' @noRd
ald_vertex_try <- function(X, y, w, p, b0, h0 = NULL) {
  v <- ald_vertex(X, y, w, p, b0, h0)
  if (is.null(v) && length(h0)) v <- ald_vertex(X, y, w, p, b0)
  v %||% ald_vertex_ties(X, y, w, p, b0)
}

#' One local run of ald_polish(): from `par`, solve the location exactly
#' at the current scale and quantile, take a Newton step in the rest,
#' and repeat until the Newton gain is below `ald_gain_tol`. Returns
#' `par`, `gain`, `margin`, `rounds`, `bases` and `why` (NULL when the
#' run settled at a proved point).
#'
#' @noRd
ald_local <- function(fit, blocks, par, fn, inside, bases = list()) {
  obj <- fit$obj
  bpos <- unlist(lapply(blocks, `[[`, "pos"))
  tpos <- setdiff(seq_along(par), bpos)
  margin <- Inf
  gain <- 0
  why <- NULL
  settled <- FALSE
  at_rows <- function(v, rows) {
    v <- as.numeric(v)
    if (length(v) == 1L) rep(v, length(rows)) else v[rows]
  }
  for (rd in seq_len(ald_polish_rounds)) {
    fp <- fit
    fp$estimates <- ald_est_at(fit, par)
    dp <- tryCatch(eval_dpars(fp, b = NULL), error = keep_missing(NULL))
    if (is.null(dp)) {
      why <- "the distributional parameters could not be evaluated"
      break
    }
    margin <- Inf
    for (blk in blocks) {
      d <- dp[[blk[["resp"]]]]
      s <- at_rows(d[["sigma"]], blk[["rows"]])
      q <- at_rows(d[["quantile"]], blk[["rows"]])
      if (!all(is.finite(s) & s > 0) || !all(is.finite(q) & q > 0 & q < 1)) {
        why <- "a scale or quantile is at the end of its range"
        break
      }
      v <- ald_vertex_try(blk[["X"]], blk[["y"]], blk[["w"]] / s, q,
                          par[blk[["pos"]]], bases[[blk[["resp"]]]])
      if (is.null(v)) {
        why <- paste0("no optimal vertex of the quantile regression of ",
                      "response ", blk[["resp"]], " was reached")
        break
      }
      par[blk[["pos"]]] <- v$b
      bases[[blk[["resp"]]]] <- v$h
      margin <- min(margin, v$margin)
    }
    if (!is.null(why)) break
    if (!length(tpos)) {
      gain <- 0
      settled <- TRUE
      break
    }
    g <- tryCatch(drop(obj$gr(par))[tpos], error = function(e) NULL)
    H <- tryCatch(obj_he_cached(obj, par)[tpos, tpos, drop = FALSE],
                  error = function(e) NULL)
    R <- if (!is.null(g) && all(is.finite(g)) && is.matrix(H) &&
               all(is.finite(H))) {
      tryCatch(chol((H + t(H)) / 2), error = function(e) NULL)
    }
    if (is.null(R)) {
      why <- paste0("the curvature in the scale, the quantile or another ",
                    "smooth parameter is not positive definite")
      break
    }
    step <- backsolve(R, forwardsolve(t(R), g))
    gain <- sum(g * step) / 2
    # no step can gain on a point that is level to rounding, so a failed
    # line search ends the rounds as a small gain does
    settled <- gain <= ald_gain_tol
    if (settled) break
    fpar <- fn(par)
    a <- 1
    for (k in 1:30) {
      cand <- par
      cand[tpos] <- par[tpos] - a * step
      fc <- if (inside(cand)) fn(cand) else NA_real_
      if (is.finite(fc) && is.finite(fpar) && fc < fpar) break
      a <- a / 2
    }
    if (!(is.finite(fc) && is.finite(fpar) && fc < fpar)) {
      settled <- TRUE
      break
    }
    par <- cand
  }
  if (is.null(why) && !settled) {
    why <- paste0("the scale and the quantile did not settle in ",
                  ald_polish_rounds, " rounds")
  }
  list(par = par, gain = gain, margin = margin, rounds = rd, bases = bases,
       why = why)
}

#' Which global proof a fit admits, from the shape of its likelihood:
#' `"convex"`, `"search"` or a sentence saying why neither.
#'
#' With a fixed quantile and a constant sigma, the negative
#' log-likelihood of the asymmetric Laplace part is
#' `W log sigma + sum w_i rho_p(y_i - x_i b) / sigma`; in `t = 1 / sigma`
#' and `a = b / sigma` it is `-W log t + sum w_i rho_p(t y_i - x_i a)`,
#' since `rho_p` is positively homogeneous. Both terms are convex, so
#' every local maximum is the global one (a `zi` part is a logistic
#' regression on the zeros, convex and apart). With the quantile
#' estimated too, the profile in it has several local maxima
#' (ald_quantile_search() then decides). With a scale or quantile model
#' (`sigma ~ z`, `quantile ~ z`) neither holds: with the quantile fixed
#' at 0.2 and `sigma ~ z`, 4 of 40 fits stopped up to 0.153 below a
#' multi-start maximum (dev/alapse-fixsig.R), and the review found
#' `sigma ~ z` and `quantile ~ z` fits with an estimated quantile 0.313
#' and 1.89 below.
#'
#' The constant's link may be any that maps its coefficient one-to-one
#' onto the parameter's range: the likelihood is then the same function
#' of sigma or the quantile. A box on log sigma keeps the concave case
#' (an interval in `t` is convex, and the local proof needs an interior
#' point); the search ranges over every quantile and profiles sigma
#' freely, so a bound on either leaves the fit local (the short review
#' of the fix round, dev/reviews/alapse-short.md, finding 1: a bound
#' was crossed). Returns the case with each constant's link, or the
#' sentence with attribute `advice`, the remedy the warning gives.
#'
#' @noRd
ald_global_case <- function(fit, blocks, bounds = NULL) {
  local_why <- function(x, advice) {
    structure(x, advice = advice, class = "ald_local_why")
  }
  alone <- paste0("A model with a fixed quantile and a constant sigma, or ",
                  "an estimated constant quantile, has a proved highest ",
                  "maximum")
  rs <- names(fit$spec$responses)
  if (length(blocks) != 1L || length(rs) != 1L) {
    return(local_why("the model has more than one response",
                     paste0("Each response fitted alone with a fixed ",
                            "quantile or a constant sigma has a proved ",
                            "highest maximum")))
  }
  if (isTRUE(tryCatch(fit_prior_density(fit), error = function(e) TRUE))) {
    return(local_why("a prior density makes the objective not convex",
                     paste0("Without the prior density, ", tolower(alone))))
  }
  fr <- fit$frame
  om <- outer_par_map(fit)
  r <- rs[[1L]]
  # a constant whose link maps its one coefficient one-to-one onto the
  # parameter's range: the likelihood is the same function of the
  # parameter, so the proofs in sigma and the quantile carry over
  intercept_only <- function(lp, links) {
    !is.null(lp) && is.null(lp[["Z"]]) && is.null(lp[["nl_body"]]) &&
      is.null(lp[["offset"]]) && identical(lp[["par"]], "betad") &&
      (lp[["link"]][["name"]] %||% "") %in% links &&
      is.function(lp[["link"]][["linkfun"]]) &&
      is.function(lp[["link"]][["linkinv"]]) &&
      NCOL(lp[["X"]]) == 1L && all(lp[["X"]] == lp[["X"]][1L]) &&
      lp[["X"]][1L] > 0
  }
  ls <- fr[["linpreds"]][[linpred_key(r, "sigma")]]
  lq <- fr[["linpreds"]][[linpred_key(r, "quantile")]]
  ps <- ald_lp_pos(fit, ls, om)
  pq <- ald_lp_pos(fit, lq, om)
  if (!intercept_only(ls, dpar_links_positive) || length(ps) != 1L ||
        anyNA(ps)) {
    return(local_why(paste0("sigma is not one estimated constant, so the ",
                            "likelihood is not proved concave"),
                     paste0("Fix the quantile and drop the scale model for ",
                            "a quantile regression (bf(y ~ x, quantile = ",
                            "0.5)), whose maximum is then the highest")))
  }
  fixed_q <- !length(pq[!is.na(pq)])
  if (!fixed_q && (!intercept_only(lq, dpar_links_unit) ||
                     length(pq) != 1L)) {
    return(local_why(paste0("the quantile has a model, so the likelihood ",
                            "has several local maxima"),
                     paste0("Fix the quantile for a quantile regression ",
                            "(bf(y ~ x, quantile = 0.5)), whose maximum is ",
                            "then the highest")))
  }
  zpos <- ald_lp_pos(fit, fr[["linpreds"]][[linpred_key(r, "zi")]], om)
  known <- c(blocks[[1L]][["pos"]], ps, if (!fixed_q) pq, zpos)
  if (!all(seq_along(fit$opt$par) %in% known)) {
    return(local_why("the model has parameters the proof does not read",
                     alone))
  }
  # the search ranges over every quantile and profiles sigma freely, so
  # a bound on either leaves it open; the concave case holds in a box
  np <- length(fit$opt$par)
  lower <- rep_len(bounds$lower %||% -Inf, np)
  upper <- rep_len(bounds$upper %||% Inf, np)
  boxed <- c(ps, if (!fixed_q) pq)
  if (!fixed_q && any(is.finite(lower[boxed]) | is.finite(upper[boxed]))) {
    return(local_why(paste0("a bound on sigma or the quantile limits the ",
                            "search over the quantile"),
                     paste0("Without the bound, or with the quantile fixed, ",
                            "the maximum is proved the highest")))
  }
  list(kind = if (fixed_q) "convex" else "search", sig = ps,
       sig_c = ls[["X"]][1L], sig_link = ls[["link"]],
       q = if (!fixed_q) pq, q_c = if (!fixed_q) lq[["X"]][1L],
       q_link = if (!fixed_q) lq[["link"]])
}

#' The highest maximum of the asymmetric Laplace part over an estimated
#' constant quantile, with a constant sigma, by branch and bound on the
#' quantile `p`.
#'
#' With sigma profiled the log-likelihood in `p` is
#' `W (log p + log(1 - p) - log(S(p) / W) - 1)`, with `S(p)` the
#' minimum of `sum w_i rho_p(y_i - x_i b)` over b: a minimum of
#' functions linear in `p`, so concave and piecewise linear, with a
#' breakpoint wherever the quantile regression's solution changes. On
#' an interval `[a, b]` a concave `S` lies above its chord
#' `L(p) = c0 + c1 p` (toward the ends of `(0, 1)` the chord through
#' 0, since `S >= 0`), so the profile is at most
#' `W (log p + log(1 - p) - log L(p) - 1)`, whose maximum on the
#' interval is at an end or at a root of `c1 p^2 + 2 c0 p - c0 = 0`.
#' The bound is exact on an interval with no breakpoint inside, where
#' `S` is its chord. An interval whose bound is within `tol` of the best
#' value seen is dropped; the others are split at the bound's maximum,
#' which on an interval with no breakpoint is the profile's own, until
#' none is left. Each `S(p)` is one exact linear program
#' (ald_vertex_try(), warm from the nearest basis). Returns `p`, `S`,
#' `b`, `value` of the best point, `proved` (no interval left above
#' it) and `evals`.
#'
#' @noRd
ald_quantile_search <- function(blk, p_fit, b_fit, tol) {
  X <- blk[["X"]]
  y <- blk[["y"]]
  w <- blk[["w"]]
  W <- sum(w)
  n <- length(y)
  # state the inner functions share
  acc <- new.env(parent = emptyenv())
  acc$pts <- list()
  acc$evals <- 0L
  pvec <- function() vapply(acc$pts, `[[`, 0, "p")
  at <- function(p) {
    nb <- if (length(acc$pts)) acc$pts[[which.min(abs(pvec() - p))]]
    v <- ald_vertex_try(X, y, w, rep(p, n), nb$b %||% b_fit, nb$h)
    acc$evals <- acc$evals + 1L
    if (is.null(v)) return(NULL)
    r <- y - drop(X %*% v$b)
    S <- sum(w * r * (p - (r < 0)))
    if (!(is.finite(S) && S > 0)) return(NULL)
    pt <- list(p = p, S = S, b = v$b, h = v$h,
               value = W * (log(p) + log1p(-p) - log(S / W) - 1))
    acc$pts[[length(acc$pts) + 1L]] <- pt
    pt
  }
  prof_bound <- function(p, c0, c1) {
    W * (log(p) + log1p(-p) - log(c0 + c1 * p) + log(W) - 1)
  }
  # the bound's maximum on [a, b], and where it is
  bound <- function(a, b, Sa, Sb) {
    if (a == 0) return(c(W * (log(b / Sb) + log(W) - 1), (a + b) / 2))
    if (b == 1) {
      return(c(W * (log((1 - a) / Sa) + log(W) - 1), (a + b) / 2))
    }
    c1 <- (Sb - Sa) / (b - a)
    c0 <- Sa - c1 * a
    cand <- c(a, b)
    disc <- c0 * (c0 + c1)
    if (c1 == 0) {
      cand <- c(cand, 0.5)
    } else if (disc >= 0) {
      cand <- c(cand, (-c0 + c(-1, 1) * sqrt(disc)) / c1)
    }
    cand <- cand[cand >= a & cand <= b]
    v <- prof_bound(cand, c0, c1)
    k <- which.max(v)
    at_p <- cand[k]
    # a split at an end would not shrink the interval
    if (at_p <= a + 1e-3 * (b - a) || at_p >= b - 1e-3 * (b - a)) {
      at_p <- (a + b) / 2
    }
    c(v[k], at_p)
  }
  grid <- sort(unique(c(seq_len(31L) / 32, p_fit)))
  for (p in grid) {
    if (is.null(at(p))) return(list(proved = FALSE, evals = acc$evals))
  }
  bounds_now <- function() {
    o <- order(pvec())
    ps <- c(0, pvec()[o], 1)
    Ss <- c(0, vapply(acc$pts, `[[`, 0, "S")[o], 0)
    m <- length(ps) - 1L
    out <- vapply(seq_len(m), function(j) {
      bound(ps[j], ps[j + 1L], Ss[j], Ss[j + 1L])
    }, numeric(2))
    list(ub = out[1L, ], split = out[2L, ], width = diff(ps))
  }
  for (it in seq_len(ald_search_evals)) {
    best <- max(vapply(acc$pts, `[[`, 0, "value"))
    bn <- bounds_now()
    open <- which(bn$ub > best + tol & bn$width > 1e-12)
    if (!length(open)) break
    j <- open[which.max(bn$ub[open])]
    if (is.null(at(bn$split[j]))) {
      return(list(proved = FALSE, evals = acc$evals))
    }
  }
  vals <- vapply(acc$pts, `[[`, 0, "value")
  top <- acc$pts[[which.max(vals)]]
  bn <- bounds_now()
  c(top, list(proved = all(bn$ub <= max(vals) + tol), evals = acc$evals))
}


#' The most linear programs ald_quantile_search() solves: on the review's
#' design B (100 rows) and the lane's designs the search ends well
#' below it (dev/alapse-findings.md).
#'
#' @noRd
ald_search_evals <- 600L

#' Move an asym_laplace() fit without random effects to the exact
#' quantile regression solution, prove that it is a maximum, and prove
#' it the highest where the likelihood's shape allows.
#'
#' The optimizer stops on the kinks of the likelihood (code 8, or a
#' gradient of 0.1 to 3 at code 0 or 4): measured against
#' quantreg::rq() at the fit's own scale and quantile, 43 of 179 fits on
#' 0.83.0 stopped more than 2e-5 below it, up to 0.042
#' (dev/alapse-rq.R, dev/alapse-findings.md). With the scale and the
#' quantile held, the location coefficients minimize a weighted check
#' function, a linear program, so they are found exactly
#' (ald_vertex_try()); the scale, the quantile and every other parameter
#' enter smoothly, so a Newton step on them, with the location solved
#' again, settles them (ald_local()).
#'
#' The local proof: the location is the exact optimum for the scale and
#' quantile held (the vertex condition), and the Newton gain over the
#' rest is within `min(conv_at_max_tol, grad_tol)` on a Hessian that is
#' positive definite. The kinks sit on hyperplanes in the location
#' alone, so the derivative in the other parameters is continuous across
#' them and this is the first-order condition of the whole problem, at a
#' minimum in the location; so the point is a local maximum of the
#' likelihood (strict when the vertex condition holds with a positive
#' margin; a margin of 0 is a solution that is not unique).
#'
#' The global proof (ald_global_case()): with a fixed quantile and a
#' constant sigma the likelihood is concave in a reparameterization, so
#' the local maximum is the global one; with an estimated constant
#' quantile ald_quantile_search() finds the highest maximum over the
#' quantile, and the local run restarts there when it is higher. Only a
#' fit proved global gets code 0, its message says so, and the
#' optimizer's own code and message are kept in `ald_cert`. A fit
#' proved only local keeps the optimizer's code and gets the message and
#' the warning of a local maximum (ald_local_message()). Returns the
#' fit's `opt` with `ald_cert` (`ok`, `global`, `why` when not proved,
#' `local` when not global, and the numbers); NULL for a fit this does
#' not apply to.
#'
#' @noRd
ald_polish <- function(fit, bounds) {
  rs <- ald_fixed_resps(fit)
  if (!length(rs)) return(NULL)
  opt <- fit$opt
  fail <- function(why) {
    opt$ald_cert <- list(ok = FALSE, why = why)
    opt
  }
  blocks <- ald_lp_blocks(fit, rs, bounds)
  if (is.character(blocks)) return(fail(blocks))
  obj <- fit$obj
  saved <- obj_state_save(obj)
  restore <- function() obj_state_restore(obj, saved)
  p0 <- opt$par
  fn <- function(q) {
    v <- tryCatch(se_quiet(obj$fn(q)), error = function(e) NA_real_)
    if (is.finite(v)) v else NA_real_
  }
  f0 <- fn(p0)
  if (!is.finite(f0)) {
    restore()
    return(fail("the objective is not finite at the stop"))
  }
  bpos <- unlist(lapply(blocks, `[[`, "pos"))
  tpos <- setdiff(seq_along(p0), bpos)
  # a step stays strictly inside the optimizer's box, where the
  # first-order condition is the gradient's
  lower <- rep_len(bounds$lower %||% -Inf, length(p0))[tpos]
  upper <- rep_len(bounds$upper %||% Inf, length(p0))[tpos]
  inside <- function(q) all(q[tpos] > lower & q[tpos] < upper)
  if (!inside(p0)) {
    restore()
    return(fail("a bound holds a scale or quantile coefficient"))
  }
  tol <- min(conv_at_max_tol, fit$control$grad_tol %||% 1e-3)
  proved <- function(run) {
    if (!is.null(run$why)) return(run$why)
    if (!(run$gain <= tol)) {
      return(paste0("a Newton step in the scale and the quantile would ",
                    "still gain ", format(run$gain, digits = 3)))
    }
    if (run$margin < -ald_margin_tol) {
      return("the vertex condition does not hold")
    }
    NULL
  }
  run <- ald_local(fit, blocks, p0, fn, inside)
  why <- proved(run)
  if (!is.null(why)) {
    restore()
    return(fail(why))
  }
  # the global proof
  case <- ald_global_case(fit, blocks, bounds)
  global <- FALSE
  search <- NULL
  local <- if (inherits(case, "ald_local_why")) as.vector(case)
  advice <- attr(case, "advice")
  quantile_advice <- paste0("Fix the quantile for a quantile ",
                            "regression (bf(y ~ x, quantile = 0.5)), ",
                            "whose maximum is then the highest")
  if (is.null(local) && identical(case$kind, "convex")) global <- TRUE
  if (is.null(local) && identical(case$kind, "search")) {
    blk <- blocks[[1L]]
    p_fit <- case$q_link$linkinv(case$q_c * run$par[case$q])
    search <- ald_quantile_search(blk, p_fit, run$par[blk[["pos"]]], tol)
    r_fit <- blk[["y"]] - drop(blk[["X"]] %*% run$par[blk[["pos"]]])
    S_fit <- sum(blk[["w"]] * r_fit * (p_fit - (r_fit < 0)))
    W <- sum(blk[["w"]])
    v_fit <- W * (log(p_fit) + log1p(-p_fit) - log(S_fit / W) - 1)
    if (!isTRUE(search$proved)) {
      local <- paste0("the search over the quantile did not settle in ",
                      ald_search_evals, " linear programs")
    } else if (search$value > v_fit + tol) {
      # a higher maximum at another quantile: start there
      start <- run$par
      start[case$q] <- case$q_link$linkfun(search$p) / case$q_c
      start[case$sig] <- case$sig_link$linkfun(search$S / W) / case$sig_c
      start[blk[["pos"]]] <- search$b
      run2 <- ald_local(fit, blocks, start, fn, inside)
      p2 <- case$q_link$linkinv(case$q_c * run2$par[case$q])
      r2 <- blk[["y"]] - drop(blk[["X"]] %*% run2$par[blk[["pos"]]])
      S2 <- sum(blk[["w"]] * r2 * (p2 - (r2 < 0)))
      v2 <- W * (log(p2) + log1p(-p2) - log(S2 / W) - 1)
      if (is.null(proved(run2)) && inside(run2$par) &&
            v2 >= search$value - tol &&
            isTRUE(fn(run2$par) < fn(run$par))) {
        search$moved_from <- p_fit
        run <- run2
        global <- TRUE
      } else {
        local <- if (min(search$p, 1 - search$p) < 0.01) {
          paste0("the likelihood is higher, ", format(search$value -
                   v_fit, digits = 3), " above the fit, at quantile ",
                 format(search$p, digits = 3), ", toward the end of its ",
                 "range, where it need not have a maximum")
        } else "the run from the highest quantile found did not settle"
      }
    } else {
      global <- TRUE
    }
  }
  par <- run$par
  f1 <- fn(par)
  if (!(is.finite(f1) && f1 <= f0 + 1e-12 * max(1, abs(f0)))) {
    restore()
    return(fail(paste0("the exact solution's objective is not below the ",
                       "stop's (a term the proof does not read)")))
  }
  new <- opt
  new$par[] <- par
  new$objective <- f1
  if (global) {
    new$convergence <- 0L
    new$message <- paste0("the exact quantile regression solution, the ",
                          "maximum of the likelihood (asym_laplace())")
  } else {
    # code 0 says the maximum is found, which only the global proof says
    new$convergence <- 1L
    new$message <- paste0("the exact quantile regression solution, a ",
                          "local maximum of the likelihood ",
                          "(asym_laplace())")
  }
  if (!global && is.null(advice)) advice <- quantile_advice
  new$ald_cert <- list(ok = TRUE, global = global, local = local,
                       advice = advice,
                       code = opt$convergence, message = opt$message,
                       gain = run$gain, margin = run$margin, gap = f0 - f1,
                       moved = max(abs(par - p0)), rounds = run$rounds,
                       search = if (!is.null(search)) {
                         search[c("p", "value", "evals", "proved",
                                  intersect("moved_from", names(search)))]
                       })
  obj_settle(obj, opt, new, saved)
  new
}

#' The warning of a fit proved only a local maximum, class
#' `frmtmb_ald_local`: why the proof is local, and the remedy.
#'
#' @noRd
ald_local_message <- function(fit) {
  cert <- fit$opt$ald_cert
  why <- cert$local %||% "its shape was not read"
  paste0("The ", paste(ald_family_names(fit), collapse = " and "),
         " fit is the exact quantile regression solution at the fitted ",
         "sigma and quantile, and a local maximum of the likelihood, but ",
         "it is not proved the highest: ", why, ". Such a likelihood can ",
         "have several local maxima (with an estimated quantile, simulated ",
         "fits stopped up to 2.65 below the highest). ",
         cert$advice %||% "Fix the quantile for a quantile regression")
}

# Every row censored on one side ----------------------------------------

#' An asym_laplace() fit whose every row is censored on the same side, in
#' the form separation_check() returns (`list(message, pars, named)`),
#' or NULL. Its maximum likelihood estimate does not exist (backlog item
#' 5490), as for a gaussian response, which the limit rule of the
#' standard-error check names.
#'
#' The proof reads the data alone. A right-censored row reads
#' `P(Y > c)`, below 1 at every finite parameter value and rising to 1 as
#' mu grows, whatever sigma and the quantile are; a left-censored row the
#' same as mu falls. With every row on one side and an intercept in mu on
#' the identity link, free that way, the likelihood rises toward 1 along
#' the intercept and never reaches it. With random effects in mu the
#' same holds for the marginal (each row's probability goes to 1 for
#' every value of the effects) and for its Laplace approximation, which
#' is below 1 everywhere, since each row's log probability is concave in
#' mu (the asymmetric Laplace density is log-concave) so the data only
#' add curvature, and goes to 1 along the same path. A prior density, a
#' block that is not gaussian, a nonlinear or another link, truncation
#' or an exact row leave the fit to the usual checks.
#'
#' The limit rule does not reach this stop itself: the optimizer runs
#' until every probability is 1 in double precision (on
#' dev/wcheck2-allcens.R's data the quantile at 3.4e-217 without the
#' group term, sigma at exp(-38.6) with it), where every direction it
#' probes is level, so it calls the parameters flat.
#'
#' @noRd
ald_cens_sep <- function(fit) {
  rs <- names(fit$spec$responses)
  if (length(rs) != 1L || !length(fit$opt$par) ||
        !identical(fit$spec$responses[[rs]]$family[["family"]],
                   "asym_laplace") ||
        isTRUE(tryCatch(fit_prior_density(fit), error = function(e) TRUE))) {
    return(NULL)
  }
  fr <- fit$frame
  key <- linpred_key(rs, "mu")
  for (bk in fr[["re_blocks"]] %||% list()) {
    keys <- vapply(bk[["components"]] %||% list(), function(cp) {
      cp[["lp_key"]] %||% ""
    }, "")
    if (!all(keys == key) ||
          !identical(bk[["covstruct_package"]] %||% "frmtmb", "frmtmb") ||
          !is.null(bk[["dist_nu"]]) || !is.null(bk[["t_share"]]) ||
          bk[["covstruct"]] %in% c("us_t", "diag_t")) {
      return(NULL)
    }
  }
  av <- fr[["aterm_values"]][[rs]]
  if (!is.null(av[["trunc_lb"]]) || !is.null(av[["trunc_ub"]])) return(NULL)
  y <- as.numeric(fr[["y"]][[rs]])
  w <- rep_len(as.numeric(av[["weights"]] %||% 1), length(y))
  w[fr[["mi_map"]][[rs]][["unread"]]] <- 0
  cens <- rep_len(av[["cens"]] %||% 0, length(y))[is.finite(y) & w > 0]
  side <- if (length(cens) && all(cens == 1)) 1 else
    if (length(cens) && all(cens == -1)) -1 else 0
  lp <- fr[["linpreds"]][[key]]
  if (side == 0 || is.null(lp) || !is.null(lp[["nl_body"]]) ||
        length(lp[["mo"]] %||% list()) ||
        !identical(lp[["link"]][["name"]] %||% "", "identity") ||
        !identical(lp[["par"]], "beta") || !NCOL(lp[["X"]])) {
    return(NULL)
  }
  X <- as.matrix(lp[["X"]])
  j <- which(apply(X, 2L, function(v) all(v == v[1L]) && v[1L] > 0))
  pos <- ald_lp_pos(fit, lp, outer_par_map(fit))[j]
  box <- fit_outer_box(fit)
  np <- length(fit$opt$par)
  end <- if (side > 0) rep_len(box$upper %||% Inf, np)[pos] else
    rep_len(box$lower %||% -Inf, np)[pos]
  if (!length(j) || anyNA(pos) || !any(is.infinite(end))) return(NULL)
  nm <- outer_par_names(fit)
  msg <- paste0(
    "Every row of the asym_laplace() response is ",
    if (side > 0) "right" else "left", " censored, so ",
    if (side > 0) "a larger" else "a smaller", " intercept of mu raises ",
    "every row's probability toward 1, which no finite value reaches: ",
    "the likelihood has no maximum, and the maximum likelihood ",
    "estimates do not exist. The reported values are where the ",
    "optimizer stopped, and no parameter has a standard error. Exact ",
    "rows, or a prior on the intercept (set_prior()), give the model an ",
    "estimate")
  list(message = msg, pars = nm, named = nm)
}
