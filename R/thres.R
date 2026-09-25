# brms's thres() addition term for the ordinal families.
#
# `y | thres(x, gr) ~ ...` does two things in brms. `x` sets the number
# of thresholds, for a response whose top categories are not all
# observed. `gr` gives every level of a factor a threshold vector of its
# own, merged into one parameter vector, with the number of thresholds
# allowed to differ by level ("grouped thresholds"). The linear
# predictor is shared.
#
# A count alone changes nothing but the length of the threshold vector,
# so it reaches the existing densities through `extra_pars`. Groups
# change the density: a row's thresholds are a slice of the merged
# vector, and the slice has a length of its own. Those densities are
# built here, once the response is in hand, by the ordinal families'
# `family_finalize()` slot, and the resolved layout rides on the family
# object as `fam[["thres"]]`, where prediction and simulation on new
# data read it back.

#' The group codes of a `thres(gr = )` factor, with the levels attached.
#'
#' brms reads the grouping variable through `factor()`, so a character
#' vector takes the sorted order and a factor keeps its own, and it
#' refuses a variable that is not factor-like. The levels are carried as
#' an attribute so that the family can name its thresholds; the tape
#' only ever sees the codes.
#'
#' @noRd
thres_group_codes <- function(v) {
  if (!(is.factor(v) || is.character(v) || is.logical(v))) {
    frm_stop("Variable 'gr' in 'thres' needs to be factor-like. ",
             "thres(gr = ) gives each level its own threshold vector, so ",
             "it takes a factor, a character or a logical variable, not ",
             arg_desc(v), call. = FALSE)
  }
  f <- factor(v)
  structure(as.numeric(f), thres_levels = levels(f))
}

#' The layout of the merged threshold vector: per group, its count and
#' where its slice starts and ends.
#'
#' @noRd
thres_layout <- function(nthres) {
  nthres <- as.integer(nthres)
  end <- cumsum(nthres)
  list(nthres = nthres, start = end - nthres + 1L, end = end,
       G = length(nthres), K1max = max(nthres))
}

#' Whether a family carries grouped thresholds.
#'
#' @noRd
thres_grouped <- function(fam) {
  isTRUE(fam[["thres"]][["grouped"]])
}

#' The group code of every row, as the density indexes it.
#'
#' @noRd
thres_row_groups <- function(aterms, n) {
  gi <- aterms[["thres_gr"]]
  if (is.null(gi) || length(gi) != n) {
    frm_stop("This ordinal model has grouped thresholds, thres(gr = ), and ",
             "the rows reached the density without their group. Supply ",
             "the grouping variable in newdata", call. = FALSE)
  }
  as.integer(gi)
}

#' The merged internal vector mapped to the thresholds, one slice per
#' group. cumulative holds each slice as (first threshold, log
#' increments); sratio, cratio and acat hold the thresholds themselves.
#' The loop runs over the parameter layout, not over rows, and works on
#' the tape and on plain doubles alike.
#'
#' @noRd
thres_tau <- function(raw, lay, ordered) {
  if (!ordered) return(raw)
  "[<-" <- RTMB::ADoverload("[<-")
  tau <- raw
  for (g in seq_len(lay$G)) {
    s <- lay$start[g]
    e <- lay$end[g]
    if (e > s) for (k in (s + 1L):e) tau[k] <- tau[k - 1L] + exp(raw[k])
  }
  tau
}

#' The number of thresholds on each group, from the data, following
#' brms's `extract_thres_names()`: the count `thres(x = )` gives, or the
#' highest category observed in the group less one. An ordered-factor
#' response counts the same way, because the levels no row takes are
#' dropped with the model frame, in brms as here (measured on brms
#' 2.23.0, dev/thres-findings.md).
#'
#' @noRd
thres_counts <- function(y, x, gi, G, labels) {
  if (!is.null(x)) {
    if (any(!is.finite(x)) || any(x != round(x)) || any(x < 1)) {
      frm_stop("Number of thresholds must be a positive integer. ",
               "thres(x = ) takes the number of thresholds, one less than ",
               "the number of response categories", call. = FALSE)
    }
  }
  if (is.null(gi)) {
    if (is.null(x)) return(max(y) - 1L)
    if (length(x) > 1L) {
      frm_stop("Number of thresholds needs to be a single value. Without ",
               "gr = , one threshold vector serves every row, so thres() ",
               "takes one number; a per-row count needs thres(x, gr = )",
               call. = FALSE)
    }
    return(as.integer(x))
  }
  if (is.null(x)) {
    out <- vapply(seq_len(G), function(g) max(y[gi == g]) - 1, 0)
  } else if (length(x) == 1L) {
    out <- rep(x, G)
  } else {
    out <- vapply(seq_len(G), function(g) {
      v <- unique(x[gi == g])
      if (length(v) > 1L) {
        frm_stop("Number of thresholds should be unique for each group. ",
                 "Level '", labels[g], "' of thres(gr = ) has rows with ",
                 paste(sort(v), collapse = ", "), call. = FALSE)
      }
      v
    }, 0)
  }
  if (any(out < 1)) {
    frm_stop("Could not extract the number of thresholds for level(s) ",
             paste0("'", labels[out < 1], "'", collapse = ", "),
             " of thres(gr = ): every row there is in the first category, ",
             "so the level holds no threshold to estimate. Use ordered ",
             "factors or positive integers as your ordinal response and ",
             "ensure that more than one response category is present in ",
             "each group, or set the count with thres(x = , gr = )",
             call. = FALSE)
  }
  as.integer(out)
}

#' The ordinal families' `family_finalize()` slot: resolve `thres()`
#' against the response and rebuild the pieces of the family that read
#' the thresholds.
#'
#' Without `thres()` the family comes back unchanged. A count alone only
#' lengthens the threshold vector. Groups replace the density, the
#' simulator, the start values and the map `variables()` reads the
#' thresholds through.
#'
#' @noRd
thres_finalizer <- function(family, ordered, link) {
  force(family)
  force(ordered)
  force(link)
  function(fam, y, aterms) {
    x <- aterms[["thres"]]
    gi <- aterms[["thres_gr"]]
    if (is.null(x) && is.null(gi)) return(fam)
    labels <- attr(gi, "thres_levels")
    G <- length(labels)
    if (!is.null(gi)) gi <- as.integer(gi)
    nthres <- thres_counts(y, x, gi, G, labels)
    g_all <- if (is.null(gi)) rep(1L, length(y)) else gi
    seen <- vapply(seq_along(nthres), function(g) max(y[g_all == g]), 0)
    short <- which(seen > nthres + 1)
    if (length(short)) {
      frm_stop("Number of thresholds is smaller than required by the ",
               "response. thres(x = ) counts thresholds, one less than ",
               "the number of categories, and the response reaches ",
               "category ", seen[short[1L]],
               if (!is.null(gi)) paste0(" in level '", labels[short[1L]],
                                        "' of thres(gr = )"),
               call. = FALSE)
    }
    lay <- thres_layout(nthres)
    # the thresholds above a group's highest observed category bound
    # categories nobody chose; see thres_fit_check()
    unident <- unlist(lapply(seq_along(nthres), function(g) {
      if (seen[g] > nthres[g]) integer(0) else
        (lay$start[g] + seen[g] - 1L):lay$end[g]
    }))
    fam[["thres"]] <- list(grouped = !is.null(gi), nthres = nthres,
                           levels = if (is.null(gi)) "" else labels,
                           groups = if (is.null(gi)) "" else
                             brms_rename(labels),
                           unident = unident)
    if (length(unident)) fam[["post"]][["fit_check"]] <- thres_fit_check
    fam[["extra_pars"]] <- function(y, aterms) {
      if (is.null(aterms[["thres_gr"]])) {
        return(ord_tau_init(y, ordered, link, K = nthres + 1L))
      }
      g_row <- as.integer(aterms[["thres_gr"]])
      raw <- unlist(lapply(seq_len(lay$G), function(g) {
        ord_tau_init(y[g_row == g], ordered, link,
                     K = lay$nthres[g] + 1L)$tau_raw
      }))
      list(tau_raw = raw)
    }
    if (is.null(gi)) return(fam)
    fam[["lpdf"]] <- thres_lpdf(family, lay, ordered, link)
    fam[["sim"]] <- thres_sim(family, lay, ordered, link)
    fam[["post"]][["ord_thresholds"]] <- function(raw) {
      unlist(lapply(seq_len(lay$G), function(g) {
        ord_tau_from_raw(raw[lay$start[g]:lay$end[g]], ordered)
      }))
    }
    fam
  }
}

#' The grouped-threshold log-density of one ordinal family.
#'
#' Row `i` reads the slice of its group, `nthres[g_i]` thresholds long,
#' so the categories it can fall in are `1..nthres[g_i] + 1`; a row
#' asked about a category beyond that gets log probability `-Inf`,
#' which makes the density a proper pmf on `1..max(nthres) + 1` for
#' every row and lets `fitted()` read category probabilities straight
#' out of it. The indices are data, built once per tape.
#'
#' @noRd
thres_lpdf <- function(family, lay, ordered, link) {
  Fcdf <- link$linkinv
  q <- link[["logit_eta"]]
  lg <- ord_log_cdf_pair(link)
  # the stopping and continuing log-probabilities as functions of the
  # `tau_j - eta` column, exactly as fam_sratio() and fam_cratio() form
  # them: log space where the link has a log-odds form
  hz <- switch(family,
    sratio = if (!is.null(lg)) {
      list(stop = lg$lF, go = lg$l1mF)
    } else {
      list(stop = function(x) log(Fcdf(x)),
           go = function(x) log(1 - Fcdf(x)))
    },
    cratio = if (!is.null(lg)) {
      list(stop = function(x) lg$l1mF(-x), go = function(x) lg$lF(-x))
    } else {
      list(stop = function(x) log(1 - Fcdf(-x)),
           go = function(x) log(Fcdf(-x)))
    },
    NULL)
  function(y, dpars, aterms, extra) {
    "[<-" <- RTMB::ADoverload("[<-")
    if (!is.null(osa_unwrap(y))) {
      frm_stop("residuals(type = \"osa\") is not available with grouped ",
               "thresholds, thres(gr = ): the one-step density selects ",
               "the category arithmetically over one shared set of ",
               "categories, and here the set differs by group. ",
               "dharma_residuals() is the simulation-based check, and it ",
               "draws each row from its own group's thresholds",
               call. = FALSE)
    }
    n <- length(y)
    eta <- dpars[["mu"]]
    if (length(eta) < n) eta <- eta + numeric(n)
    tau <- thres_tau(extra$tau_raw, lay, ordered)
    gi <- thres_row_groups(aterms, n)
    s <- lay$start[gi]
    nk <- lay$nthres[gi]
    out <- switch(family,
      cumulative = thres_lpdf_cumulative(y, eta, tau, s, nk, Fcdf, q),
      acat = thres_lpdf_acat(y, eta, tau, s, nk, gi, lay),
      thres_lpdf_seq(y, eta, tau, s, nk, lay$K1max, hz))
    off <- which(y > nk + 1)
    if (length(off)) out[off] <- -Inf
    out
  }
}

#' Grouped cumulative: `F(tau_y - eta) - F(tau_{y-1} - eta)` on the
#' row's own slice, in log space for a link with a log-odds form, like
#' `fam_cumulative()`.
#'
#' @noRd
thres_lpdf_cumulative <- function(y, eta, tau, s, nk, Fcdf, q) {
  "[<-" <- RTMB::ADoverload("[<-")
  iK <- as.numeric(y == nk + 1)
  i1 <- as.numeric(y == 1)
  if (is.null(q)) {
    up <- Fcdf(tau[s + pmin(y, nk) - 1L] - eta) * (1 - iK) + iK
    lo <- Fcdf(tau[s + pmax(y - 1, 1) - 1L] - eta) * (1 - i1)
    return(log(up - lo))
  }
  out <- i1 * log_inv_logit(q(tau[s] - eta)) +
    iK * log1m_inv_logit(q(tau[s + nk - 1L] - eta))
  # the interior categories only on the rows that have one, so no row
  # evaluates a threshold pair outside its own slice
  rm <- which(y > 1 & y <= nk)
  if (length(rm)) {
    a <- q(tau[s[rm] + y[rm] - 1L] - eta[rm])
    b <- q(tau[s[rm] + y[rm] - 2L] - eta[rm])
    out[rm] <- out[rm] + RTMB::logspace_sub(-b, -a) +
      log_inv_logit(a) + log_inv_logit(b)
  }
  out
}

#' Grouped sratio and cratio: one pass per threshold position up to
#' the largest group, the positions past a row's own slice masked out.
#'
#' @noRd
thres_lpdf_seq <- function(y, eta, tau, s, nk, K1max, hz) {
  out <- 0
  for (j in seq_len(K1max)) {
    Mj <- tau[s + pmin(j, nk) - 1L] - eta
    live <- j <= nk
    out <- out + as.numeric(live & y == j) * hz$stop(Mj) +
      as.numeric(live & j < y) * hz$go(Mj)
  }
  out
}

#' Grouped acat: `P(y = r)` proportional to
#' `exp((r - 1) eta - sum_{j < r} tau_j)` over the row's own `1..K`, the
#' normalizer folded in log space and frozen past the row's top
#' category.
#'
#' @noRd
thres_lpdf_acat <- function(y, eta, tau, s, nk, gi, lay) {
  "[<-" <- RTMB::ADoverload("[<-")
  # c(0, cumsum(slice)) per group, merged: group g starts at s_g + g - 1
  ct <- numeric(lay$end[lay$G] + lay$G)
  for (g in seq_len(lay$G)) {
    o <- lay$start[g] + g - 1L
    for (k in seq_len(lay$nthres[g])) {
      ct[o + k] <- ct[o + k - 1L] + tau[lay$start[g] + k - 1L]
    }
  }
  s2 <- s + gi - 1L
  num <- 0
  den <- NULL
  for (r in seq_len(lay$K1max + 1L)) {
    Er <- (r - 1) * eta - ct[s2 + pmin(r, nk + 1L) - 1L]
    num <- num + as.numeric(y == r) * Er
    if (is.null(den)) {
      den <- Er
    } else {
      live <- as.numeric(r <= nk + 1L)
      den <- live * RTMB::logspace_add(den, Er) + (1 - live) * den
    }
  }
  num - den
}

#' The grouped-threshold simulator: each group's category distribution
#' from `ord_cat_probs()` on its own slice, padded with zeros.
#'
#' @noRd
thres_sim <- function(family, lay, ordered, link) {
  function(dpars, aterms, n, extra) {
    P <- thres_cat_probs(family, rep(dpars[["mu"]], length.out = n),
                         extra$tau_raw, thres_row_groups(aterms, n), lay,
                         ordered, link)
    cp <- t(apply(P, 1L, cumsum))
    if (n == 1L) cp <- matrix(cp, 1L, ncol(P))
    pmin(1L + rowSums(cp < stats::runif(n)), ncol(P))
  }
}

#' `n x (max(nthres) + 1)` category probabilities under grouped
#' thresholds, in plain doubles. A category past a group's own count
#' has probability zero, which is also what brms's `posterior_epred()`
#' reports there.
#'
#' @noRd
thres_cat_probs <- function(family, eta, raw, gi, lay, ordered, link) {
  tau <- thres_tau(raw, lay, ordered)
  P <- matrix(0, length(eta), lay$K1max + 1L)
  for (g in seq_len(lay$G)) {
    rows <- which(gi == g)
    if (!length(rows)) next
    tg <- tau[lay$start[g]:lay$end[g]]
    P[rows, seq_len(lay$nthres[g] + 1L)] <-
      ord_cat_probs(family, eta[rows], tg, NULL, link)
  }
  P
}

#' The group codes of newdata's `thres(gr = )` variable, against the
#' levels the fit estimated thresholds for.
#'
#' @noRd
thres_newdata_codes <- function(fam, v) {
  lv <- fam[["thres"]][["levels"]]
  code <- match(as.character(v), lv)
  if (anyNA(code)) {
    frm_stop("newdata holds level(s) ",
             paste0("'", unique(as.character(v)[is.na(code)]), "'",
                    collapse = ", "),
             " of the thres(gr = ) variable, which the fit has no ",
             "thresholds for. Its levels are ",
             paste0("'", lv, "'", collapse = ", "), call. = FALSE)
  }
  as.numeric(code)
}

#' The number of categories of an ordinal family: the largest group's
#' under grouped thresholds, else one more than the thresholds.
#'
#' @noRd
thres_ncat <- function(fam, raw) {
  th <- fam[["thres"]]
  if (isTRUE(th[["grouped"]])) return(max(th[["nthres"]]) + 1L)
  length(raw) + 1L
}

#' The index labels brms puts in a threshold's name: `k` in
#' `b_Intercept[k]`, or `g,k` in `b_Intercept[g,k]` under grouped
#' thresholds.
#'
#' @noRd
thres_labels <- function(fam, n) {
  th <- fam[["thres"]]
  if (!isTRUE(th[["grouped"]])) return(as.character(seq_len(n)))
  paste0(rep(th[["groups"]], th[["nthres"]]), ",",
         unlist(lapply(th[["nthres"]], seq_len)))
}

#' Warn when `thres(x = )` names categories above the highest one
#' observed and no prior holds the thresholds that bound them.
#'
#' Nothing in the data places such a threshold: the likelihood keeps
#' rising as it moves off to infinity, so the maximum likelihood
#' estimate is on the boundary and its standard error is meaningless.
#' brms never meets this, because its default student_t prior on class
#' "Intercept" holds every threshold. A prior here does the same, and
#' then the warning has nothing to say.
#'
#' @noRd
thres_fit_check <- function(fit, resp) {
  fam <- fit$spec$responses[[resp]]$family
  th <- fam[["thres"]]
  ent <- if (!is.null(fit$prior)) {
    resolve_prior_input(list(frame = fit$frame, spec = fit$spec),
                        fit$prior)$entries
  }
  held <- unlist(lapply(ent, function(e) {
    if (identical(e$comp, "tau_raw")) e$idx
  }))
  free <- setdiff(th[["unident"]], held)
  if (!length(free)) return(invisible(NULL))
  lab <- thres_labels(fam, sum(th[["nthres"]]))[free]
  frm_warning("thres(): the threshold", if (length(free) > 1L) "s",
              " Intercept[", paste(lab, collapse = "], Intercept["), "] ",
              if (length(free) > 1L) "bound" else "bounds",
              " response categories that no row takes, so the data do ",
              "not place ", if (length(free) > 1L) "them" else "it",
              ": the maximum likelihood estimate runs off toward ",
              "infinity and its standard error means nothing. brms ",
              "holds such thresholds with its default student_t(3, 0, ",
              "2.5) prior. Set a prior on class \"Intercept\" to do the ",
              "same, or lower the count in thres(x = )", call. = FALSE)
  invisible(NULL)
}
