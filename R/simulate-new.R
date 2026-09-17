# De novo simulation: responses from a bf() specification and supplied
# parameters, no fitted model (the glmmTMB::simulate_new analog; power
# analysis). Parameters come either on the internal scale (the
# par_template components) or on the natural scale, in the
# hypothesis()/variables() vocabulary; priors turn the same call into
# prior-predictive simulation (brms's sample_prior = "only").

# ---------------------------------------------------------------------
# Natural-scale vocabulary
# ---------------------------------------------------------------------

#' One slot per natural-scale name the model accepts, under brms's
#' names, the ones [variables()] returns:
#'   kind `"coef"`  coefficient on the LINK scale, `b_x` (comp, idx)
#'   kind `"dpar"`  a distributional parameter nobody wrote a formula
#'                  for, on its natural scale, `sigma` (comp, idx, link)
#'   kind `"sd"`    group-level SD, `sd_g__Intercept` (block, positions)
#'   kind `"cor"`   group-level correlation (block, `j < k`)
#'   kind `"simplex"` one of a mixture's weights `theta1 ... thetaK`
#'                  (group, `k`, `K`, the `idx` of the log ratio it
#'                  writes, `NA` for the last, and `to_logits`); the
#'                  whole simplex is set together, see `apply_natural()`
#'
#' `fitlike` carries `spec`, `frame` and `bform`, which is what brms's
#' names are built from. The `legacy` attribute maps the bare names
#' frmtmb used before it took brms's (`x`, `sigma_Intercept`,
#' `sd_g__Intercept` with the old sanitizing) to the brms name, so a
#' bare name is refused with its brms spelling. Names two parameters
#' would share are recorded as ambiguous rather than resolved to one.
#'
#' @noRd
nat_slots <- function(frame, fitlike) {
  tpl <- frame[["par_template"]]
  slots <- list()
  dup <- character(0)
  legacy <- character(0)
  put <- function(nm, s, old = character(0)) {
    if (nm %in% names(slots)) dup <<- c(dup, nm) else slots[[nm]] <<- s
    old <- setdiff(old, nm)
    legacy[old] <<- nm
  }

  tab <- brms_coef_table(fitlike)
  n_beta <- length(tpl[["beta"]])
  keep_d <- setdiff(seq_along(tpl[["betad"]]), frame[["betad_fixed_idx"]])
  inv <- attr(tab, "linkinv")
  fun <- attr(tab, "linkfun")
  smp <- attr(tab, "simplex")
  for (gi in seq_along(smp)) {
    sg <- smp[[gi]]
    K <- length(sg$names)
    for (k in seq_len(K)) {
      put(sg$names[k],
          list(kind = "simplex", group = gi, k = k, K = K, comp = "betad",
               idx = if (k < K) keep_d[sg$pos[k] - n_beta] else NA_integer_,
               to_logits = sg$to_logits))
    }
  }
  in_smp <- unlist(lapply(smp, `[[`, "pos"))
  for (i in setdiff(seq_len(nrow(tab)), in_smp)) {
    comp <- if (i <= n_beta) "beta" else "betad"
    idx <- if (i <= n_beta) i else keep_d[i - n_beta]
    old <- gsub("[()]", "", tab$internal[i])
    if (tab$natural[i]) {
      put(tab$brms[i],
          list(kind = "dpar", comp = comp, idx = idx,
               link = list(linkinv = inv[[i]], linkfun = fun[[i]])),
          old)
    } else {
      put(tab$brms[i], list(kind = "coef", comp = comp, idx = idx), old)
    }
  }

  san <- function(s) gsub("[^[:alnum:]_.]", "", gsub("[()]", "", s))
  sds_nm <- brms_sds_names(fitlike)
  for (bi in seq_along(frame[["re_blocks"]])) {
    bk <- frame[["re_blocks"]][[bi]]
    reg <- covstruct_registry[[bk[["covstruct"]]]]
    if (is.null(reg$from_natural)) next
    d <- bk[["dim"]]
    if (identical(bk[["covstruct"]], "smooth")) {
      # a smooth basis has one smoothing SD, not one per basis column
      put(sds_nm[bi], list(kind = "sd", blk = bi, pos = seq_len(d)),
          paste0("sds_", san(bk[["group_name"]])))
      next
    }
    g <- brms_group_name(bk)
    tn <- brms_re_rnames(fitlike, bk)
    g0 <- san(bk[["group_name"]])
    tn0 <- san(bk[["cnms"]])
    shared <- length(reg$sd_idx(d)) == 1L && d > 1L
    for (j in seq_len(d)) {
      put(paste0("sd_", g, "__", tn[j]),
          list(kind = "sd", blk = bi, pos = if (shared) seq_len(d) else j),
          paste0("sd_", g0, "__", tn0[j]))
    }
    if (reg$npar(d) > length(reg$sd_idx(d)) && d > 1L) {
      for (j in seq_len(d - 1L)) {
        for (k in seq(j + 1L, d)) {
          put(paste0("cor_", g, "__", tn[j], "__", tn[k]),
              list(kind = "cor", blk = bi, j = j, k = k),
              paste0("cor_", g0, "__", tn0[j], "__", tn0[k]))
        }
      }
    }
  }
  # a legacy spelling that is also somebody's brms name is not legacy
  legacy <- legacy[!names(legacy) %in% names(slots)]
  structure(slots, ambiguous = unique(dup), legacy = legacy)
}

#' The brms name of the coefficient slot at (`comp`, `idx`), or `NULL`.
#'
#' @noRd
nat_coef_name <- function(slots, comp, idx) {
  for (nm in names(slots)) {
    s <- slots[[nm]]
    if (s$kind %in% c("coef", "dpar") && identical(s$comp, comp) &&
        as.integer(s$idx) == as.integer(idx)) {
      return(nm)
    }
  }
  NULL
}

#' `theta` indices a natural sd slot writes to.
#'
#' @noRd
nat_sd_theta_idx <- function(frame, slot) {
  bk <- frame[["re_blocks"]][[slot$blk]]
  reg <- covstruct_registry[[bk[["covstruct"]]]]
  unique(bk[["theta_idx"]][reg$sd_idx(bk[["dim"]])[slot$pos]])
}

#' Whether `newparams` speaks the internal parameterization: every name
#' is a `par_template` component (or `"b"`). A model whose coefficients
#' are literally named `beta`/`theta` must use the internal spelling.
#'
#' @noRd
is_internal_newparams <- function(np, tpl) {
  length(np) > 0L &&
    all(names(np) %in% c(names(tpl), "b"))
}

#' Refuse the natural path where it cannot express the whole parameter
#' vector, rather than leaving parameters at silent defaults.
#'
#' @noRd
check_natural_supported <- function(frame) {
  extra <- setdiff(names(frame[["par_template"]]),
                   c("beta", "betad", "b", "theta"))
  if (length(extra)) {
    stop("Natural-scale newparams does not cover the parameter ",
         "component(s) ", paste(extra, collapse = ", "),
         "; use the internal spelling", call. = FALSE)
  }
  for (bk in frame[["re_blocks"]]) {
    if (is.null(covstruct_registry[[bk[["covstruct"]]]]$from_natural)) {
      stop("Natural-scale newparams cannot set the '", bk[["covstruct"]],
           "' block (", bk[["term_label"]], "): no inverse map from ",
           "standard deviations and correlations. Structures with one: ",
           paste(sort(names(Filter(function(e) !is.null(e$from_natural),
                                   covstruct_registry))),
                 collapse = ", "),
           ". Use the internal spelling (theta).", call. = FALSE)
    }
  }
}

#' Write natural-scale values into the internal estimate list. Blocks
#' the names touch are rebuilt from their SDs and correlations as a
#' whole, so an unset correlation is 0 rather than whatever `theta`
#' happened to hold.
#'
#' @noRd
apply_natural <- function(est, frame, slots, np) {
  bad <- setdiff(names(np), names(slots))
  old <- intersect(bad, names(attr(slots, "legacy")))
  if (length(old)) {
    # brms's names only, the strings variables() returns: accepting the
    # bare spelling too would give one parameter two names here and one
    # everywhere else
    stop("newparams takes brms's parameter names, the ones variables() ",
         "returns. Rename: ",
         paste0(old, " -> ", attr(slots, "legacy")[old], collapse = ", "),
         call. = FALSE)
  }
  if (length(bad)) {
    stop("Unknown newparams name(s): ", paste(bad, collapse = ", "),
         ". Available: ", paste(names(slots), collapse = ", "),
         call. = FALSE)
  }
  amb <- intersect(names(np), attr(slots, "ambiguous"))
  if (length(amb)) {
    stop("Ambiguous newparams name(s): ", paste(amb, collapse = ", "),
         " (two parameters share the name); use the internal spelling",
         call. = FALSE)
  }
  for (nm in names(np)) {
    if (!is.numeric(np[[nm]]) || length(np[[nm]]) != 1L) {
      stop("newparams$", nm, " must be a single number", call. = FALSE)
    }
  }

  blk_ids <- unique(unlist(lapply(names(np), function(nm) {
    s <- slots[[nm]]
    if (s$kind %in% c("sd", "cor")) s$blk else NULL
  })))
  sds <- list()
  cors <- list()
  for (bi in blk_ids) {
    bk <- frame[["re_blocks"]][[bi]]
    V <- covstruct_registry[[bk[["covstruct"]]]]$vcov(
      est[["theta"]][bk[["theta_idx"]]], bk)
    s <- sqrt(diag(as.matrix(V)))
    sds[[as.character(bi)]] <- s
    cors[[as.character(bi)]] <- if (all(s > 0)) {
      stats::cov2cor(as.matrix(V))
    } else {
      diag(nrow = length(s))
    }
  }

  smp_given <- list()
  for (nm in names(np)) {
    s <- slots[[nm]]
    v <- as.numeric(np[[nm]])
    if (s$kind == "simplex") {
      key <- as.character(s$group)
      smp_given[[key]][[nm]] <- v
      next
    }
    if (s$kind == "coef") {
      est[[s$comp]][s$idx] <- v
    } else if (s$kind == "dpar") {
      est[[s$comp]][s$idx] <- s$link$linkfun(v)
    } else if (s$kind == "sd") {
      if (v <= 0) {
        stop("newparams$", nm, " must be positive", call. = FALSE)
      }
      sds[[as.character(s$blk)]][s$pos] <- v
    } else {
      if (abs(v) >= 1) {
        stop("newparams$", nm, " must lie strictly between -1 and 1",
             call. = FALSE)
      }
      C <- cors[[as.character(s$blk)]]
      C[s$j, s$k] <- v
      C[s$k, s$j] <- v
      cors[[as.character(s$blk)]] <- C
    }
  }

  for (bi in blk_ids) {
    bk <- frame[["re_blocks"]][[bi]]
    key <- as.character(bi)
    est[["theta"]][bk[["theta_idx"]]] <-
      covstruct_registry[[bk[["covstruct"]]]]$from_natural(sds[[key]],
                                                      cors[[key]], bk)
  }

  # brms's mixing weights are one simplex: every theta<k> of it, each
  # in (0, 1), summing to one. Any subset would leave the rest to a
  # guess, so a partial simplex is refused
  for (key in names(smp_given)) {
    given <- smp_given[[key]]
    s1 <- slots[[names(given)[1L]]]
    all_nm <- names(slots)[vapply(slots, function(x) {
      identical(x$kind, "simplex") && identical(x$group, s1$group)
    }, TRUE)]
    all_nm <- all_nm[order(vapply(slots[all_nm], `[[`, 1L, "k"))]
    if (!setequal(names(given), all_nm)) {
      stop("newparams sets a mixture's weights together: give ",
           paste(all_nm, collapse = ", "), ", which sum to one (",
           "missing: ", paste(setdiff(all_nm, names(given)),
                              collapse = ", "), ")", call. = FALSE)
    }
    p <- unlist(given[all_nm])
    if (any(p <= 0 | p >= 1)) {
      stop("newparams: the mixture weights ",
           paste(all_nm, collapse = ", "), " must each lie strictly ",
           "between 0 and 1", call. = FALSE)
    }
    if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
      stop("newparams: the mixture weights ",
           paste(all_nm, collapse = ", "), " must sum to one, not ",
           format(sum(p), digits = 15), call. = FALSE)
    }
    lg <- s1$to_logits(p)
    for (nm in all_nm) {
      sl <- slots[[nm]]
      if (!is.na(sl$idx)) est[["betad"]][sl$idx] <- lg[1L, sl$k]
    }
  }
  est
}

# ---------------------------------------------------------------------
# Prior draws
# ---------------------------------------------------------------------

#' One draw from a parsed `set_prior()` distribution, on the scale the
#' class defines it (class `"sd"` draws a standard deviation).
#'
#' @noRd
rprior_one <- function(dist) {
  switch(dist$kind,
    normal = stats::rnorm(1L, dist$location, dist$scale),
    t = dist$location + dist$scale * stats::rt(1L, dist$df),
    exponential = stats::rexp(1L, dist$rate),
    logistic = stats::rlogis(1L, dist$location, dist$scale),
    gamma = stats::rgamma(1L, shape = dist$shape, rate = dist$rate),
    # the inverse gamma is the reciprocal of a gamma whose RATE is the
    # inverse gamma's scale, which is the same reading its density uses
    inv_gamma = 1 / stats::rgamma(1L, shape = dist$shape,
                                  rate = dist$scale),
    beta = stats::rbeta(1L, dist$shape1, dist$shape2),
    stop("No sampler for prior kind '", dist$kind, "'", call. = FALSE)
  )
}

#' Rejection draw inside the entry's bounds (and inside `(0, Inf)` for
#' an sd), on the entry's own scale.
#'
#' @noRd
draw_prior_entry <- function(e, label, max_try = 1000L) {
  lb <- if (is.null(e$lb) || is.na(e$lb)) -Inf else e$lb
  ub <- if (is.null(e$ub) || is.na(e$ub)) Inf else e$ub
  if (identical(e$scale, "sd")) lb <- max(lb, 0)
  if (identical(e$scale, "ordthres")) {
    # the entry is a density on a whole ORDERED vector, and one draw
    # from the marginal would not be ordered. Pinning the thresholds is
    # what a prior-predictive ordinal simulation needs anyway
    stop("A prior on ", label, " is a density on the whole threshold ",
         "vector, and drawing from it one threshold at a time would ",
         "not produce an ordered one. Pin the thresholds with ",
         "newparams = list(tau_raw = ) instead", call. = FALSE)
  }
  # a "natural" entry's density is about the dpar itself, so a draw
  # outside the link's own support has no internal value to write and is
  # rejected the same way a non-positive sd is
  ok_support <- function(v) {
    if (identical(e$scale, "sd")) return(v > 0)
    if (!identical(e$scale, "natural")) return(TRUE)
    is.finite(suppressWarnings(as.numeric(e$link$linkfun(v))))
  }
  for (it in seq_len(max_try)) {
    v <- rprior_one(e$dist)
    if (v >= lb && v <= ub && ok_support(v)) {
      return(v)
    }
  }
  stop("The prior on ", label, " did not produce a draw inside [",
       lb, ", ", ub, "] in ", max_try, " tries", call. = FALSE)
}

#' The internal value one drawn prior value stands for, undoing the
#' entry's change of variables and its centering offset.
#'
#' @noRd
prior_draw_to_internal <- function(e, v, est) {
  x <- if (identical(e$scale, "sd")) {
    log(v)
  } else if (identical(e$scale, "natural")) {
    as.numeric(e$link$linkfun(v))
  } else {
    v
  }
  o <- e$offset
  if (is.null(o)) return(x)
  x - sum(est[[o$comp]][o$idx] * o$w)
}

#' What the `pars` table reports for one drawn entry.
#'
#' The draw is on the scale the prior was WRITTEN on, and for a centered
#' intercept that is the intercept at the predictor means. The table has
#' to round-trip through `newparams =`, whose `Intercept` is the
#' intercept at zero, so an offset entry reports the value that was
#' actually written rather than the value that was drawn.
#'
#' @noRd
prior_draw_report <- function(e, v, x) {
  if (is.null(e$offset)) return(v)
  if (identical(e$scale, "sd")) return(exp(x))
  if (identical(e$scale, "natural")) return(as.numeric(e$link$linkinv(x)))
  x
}

#' Label for a prior entry: the natural name where one exists, else the
#' internal spelling.
#'
#' @noRd
prior_entry_label <- function(frame, slots, e) {
  if (e$comp %in% c("beta", "betad")) {
    nm <- nat_coef_name(slots, e$comp, e$idx)
    if (!is.null(nm)) return(nm)
    nms <- names(frame[["par_template"]][[e$comp]])
    return(gsub("[()]", "", nms[e$idx]))
  }
  for (nm in names(slots)) {
    s <- slots[[nm]]
    if (identical(s$kind, "sd") &&
        e$idx %in% nat_sd_theta_idx(frame, s)) {
      return(nm)
    }
  }
  paste0(e$comp, "_", e$idx)
}

#' Draw every prior-covered parameter once and write it into `est`.
#'
#' @noRd
draw_prior_pars <- function(est, entries, labels) {
  vals <- numeric(length(entries))
  # a centered intercept's internal value is its draw minus `means_X'b`,
  # so the coefficients that offset it have to be in `est` first. Their
  # own entries carry no offset, which is what puts them in the first
  # pass; slopes with no prior at all are already in `est`
  ord <- c(which(vapply(entries, function(e) is.null(e$offset), TRUE)),
           which(vapply(entries, function(e) !is.null(e$offset), TRUE)))
  for (i in ord) {
    e <- entries[[i]]
    v <- draw_prior_entry(e, labels[i])
    x <- prior_draw_to_internal(e, v, est)
    est[[e$comp]][e$idx] <- x
    vals[i] <- prior_draw_report(e, v, x)
  }
  list(est = est, vals = vals)
}

# ---------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------

#' Every fixed coefficient and every random-effect SD must be pinned by
#' `newparams` or by a prior. Leaving one out means simulating from a
#' default (a zero effect, a unit SD) that nobody chose.
#'
#' @noRd
check_coverage <- function(frame, slots, np_internal, np_natural,
                           entries) {
  tpl <- frame[["par_template"]]
  missing_nm <- character(0)

  for (comp in intersect(c("beta", "betad"), names(tpl))) {
    if (comp %in% names(np_internal)) next
    nms <- names(tpl[[comp]])
    for (i in seq_along(nms)) {
      if (comp == "betad" && i %in% frame[["betad_fixed_idx"]]) next
      set_nat <- any(vapply(names(np_natural), function(nm) {
        s <- slots[[nm]]
        !is.null(s) && s$kind %in% c("coef", "dpar", "simplex") &&
          identical(s$comp, comp) && !is.na(s$idx) &&
          as.integer(s$idx) == i
      }, TRUE))
      set_pr <- any(vapply(entries, function(e) {
        identical(e$comp, comp) && as.integer(e$idx) == i
      }, TRUE))
      if (!set_nat && !set_pr) {
        # name it the way the user would set it, brms's name
        smp_nm <- names(slots)[vapply(slots, function(x) {
          identical(x$kind, "simplex") && identical(x$comp, comp) &&
            identical(as.integer(x$idx), as.integer(i))
        }, TRUE)]
        if (length(smp_nm)) {
          g <- slots[[smp_nm]]$group
          smp_nm <- paste(names(slots)[vapply(slots, function(x) {
            identical(x$kind, "simplex") && identical(x$group, g)
          }, TRUE)], collapse = ", ")
        }
        missing_nm <- c(missing_nm,
                        if (length(smp_nm)) smp_nm else
                          nat_coef_name(slots, comp, i) %||%
                          gsub("[()]", "", nms[i]))
      }
    }
  }

  if (!"theta" %in% names(np_internal) && length(frame[["re_blocks"]])) {
    prior_theta <- unlist(lapply(entries, function(e) {
      if (identical(e$comp, "theta")) as.integer(e$idx) else NULL
    }))
    nat_theta <- unlist(lapply(names(np_natural), function(nm) {
      s <- slots[[nm]]
      if (!is.null(s) && identical(s$kind, "sd")) {
        nat_sd_theta_idx(frame, s)
      }
    }))
    for (bk in frame[["re_blocks"]]) {
      reg <- covstruct_registry[[bk[["covstruct"]]]]
      # a structure with no natural map needs its whole theta segment
      need <- if (is.null(reg$from_natural)) {
        bk[["theta_idx"]]
      } else {
        bk[["theta_idx"]][reg$sd_idx(bk[["dim"]])]
      }
      gone <- setdiff(need, c(prior_theta, nat_theta))
      if (length(gone)) {
        for (i in gone) {
          nm <- prior_entry_label(frame, slots,
                                  list(comp = "theta", idx = i))
          missing_nm <- c(missing_nm, nm)
        }
      }
    }
  }

  if (length(missing_nm)) {
    stop("No value for ", paste(unique(missing_nm), collapse = ", "),
         ": every coefficient and random-effect SD needs a newparams ",
         "entry or a prior", call. = FALSE)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------
# frm_simulate()
# ---------------------------------------------------------------------

#' Simulate responses from a formula and parameters
#'
#' Builds the design from `formula` and `data` exactly as [frm()]
#' would, sets the parameters from `newparams` (and/or draws them from
#' `prior`), and simulates responses. Random effects are redrawn from
#' their covariance for every simulation unless `newparams$b` supplies
#' them.
#'
#' `data` must contain a response column with values that are valid for
#' the family (any dummy values do; they only anchor the design).
#'
#' Draws come back in the response's own type, exactly as
#' [simulate()]'s do: an ordered factor for an ordinal family, an
#' unordered one for a categorical family, and a matrix column for a
#' matrix response ([multinomial()], [mixture_mvn()]).
#'
#' The structured families draw here through the same implementation
#' [simulate()] uses (see its Structured draws section):
#' `mixture(groups = ~g)` takes one class per group, [mixture_mvn()]
#' uses its class covariances, a residual
#' correlation term (`ar()`, `ma()`, ...) contributes one correlated
#' residual per group, and a family from an extension package draws
#' through whatever its structure declares. The de novo frame carries
#' those structures, so
#' nothing is lost relative to a fit; `ar()` and friends need their
#' `thetaac` entry in the internal `newparams` spelling, since a
#' correlation parameter has no natural-scale name here.
#'
#' @section Two spellings for `newparams`:
#' *Natural scale* (recommended): brms's parameter names, the strings
#' [variables()] returns, one number each.
#' - coefficients as `b_Intercept`, `b_x`, `b_sigma_Intercept`, on the
#'   LINK scale;
#' - `sigma` (and any other distributional parameter nobody wrote a
#'   formula for: `shape`, `phi`, `zi`, ...) on the RESPONSE scale, so
#'   `sigma = 0.7` is a residual SD of 0.7;
#' - a mixture's weights with no theta formula as brms's simplex, all of
#'   `theta1 ... thetaK` together, each in (0, 1) and summing to one:
#'   `theta1 = 0.9, theta2 = 0.1`;
#' - `sd_<group>__<coef>` for group-level standard deviations,
#'   `cor_<group>__<c1>__<c2>` for their correlations, and
#'   `sds_<label>_1` for a smooth's smoothing SD (`sds_sx_1` for
#'   `s(x)`). Unset correlations are 0.
#'
#' A bare name (`x`, `Intercept`, `sigma_Intercept`) is refused with
#' its brms spelling. Priors are different: `set_prior(coef = "x")`
#' takes the bare coefficient, as in brms.
#'
#' *Internal scale*: named after the parameter components - `beta`,
#' `betad`, `theta`, and optionally `b` - each a full-length vector, on
#' the internal parameterization (`theta` holds log SDs and Cholesky
#' correlation parameters, `betad` holds dispersion dpars on their link
#' scale). [par_template()] returns that layout for a formula and data
#' without fitting anything, filled with the default values and carrying
#' the name of every entry; edit it and pass it back as `newparams`.
#'
#' The two spellings cannot be mixed: `newparams` is read as internal
#' when every name is a `par_template` component (or `b`), and as
#' natural otherwise. The natural spelling covers `us`, `diag`,
#' `homdiag`, `smooth`, `gr_cov` and `gr_prec` random-effect
#' structures; other structures (`ar1`, `cs`, `toep`, GPs, ...) have no
#' inverse map from SDs and correlations and need the internal
#' spelling.
#'
#' @section Prior-predictive simulation:
#' With `prior`, each of the `nsim` simulations draws a fresh
#' parameter vector from the [set_prior()] specification and simulates
#' from it - the analog of brms's `sample_prior = "only"` followed by
#' `posterior_predict()`. Draws are taken on the scale each class
#' defines: class `"b"`/`"Intercept"` on the link scale, class `"sd"`
#' on the natural SD scale (mapped to `theta` afterwards), class
#' `"theta"` on the internal scale. `lb`/`ub` truncate by rejection.
#' The drawn values come back as `attr(result, "pars")`, one row per
#' simulation, so a prior-predictive check can relate parameters to
#' outcomes.
#'
#' That table reports each parameter on the scale `newparams` names it
#' on, which for `b_Intercept` is the intercept at ZERO. A
#' `class = "Intercept"` prior is a density on the intercept at the MEAN
#' of the predictors (see the Where an intercept prior lands section of
#' [set_prior()]), so on a design with uncentered predictors the number
#' drawn and the number reported differ by `colMeans(X)` times the
#' slopes. Reporting the written value is what lets a `pars` row be
#' handed straight back as `newparams`. A prior on a distributional
#' parameter's own brms class is likewise drawn on that parameter and
#' reported there, not on its link scale.
#'
#' The column name says which scale a dpar's draw is on, as brms's
#' names do: `sigma` holds sigma, drawn under `set_prior(class =
#' "sigma")` on a model that gives sigma no formula, and
#' `b_sigma_Intercept` holds log sigma, drawn under `set_prior(class =
#' "Intercept", dpar = "sigma")` on a model that does.
#'
#' Parameters without a prior keep their `newparams` value. Whenever
#' `prior` are used, or `newparams` uses the natural spelling, every
#' fixed coefficient and every random-effect SD must be pinned by one
#' or the other; an unpinned parameter is an error rather than a silent
#' zero effect or unit SD.
#'
#' @param formula A `bf()` formula (with a family attached) or a plain
#'   formula plus `family`.
#' @param data A data frame of model data, including a dummy response
#'   column.
#' @param family Family, when `formula` does not carry one.
#' @param newparams Named list of parameters, in either spelling (see
#'   Details). [par_template()] discovers the internal spelling for a
#'   formula and data. Optional when `prior` pins everything.
#' @param prior A [set_prior()] specification to draw parameters from,
#'   once per simulation, or a `brmsprior` object brms built, which is
#'   translated row by row. The argument takes brms's spelling,
#'   `prior`; the `priors` of releases before 0.43 is gone rather than
#'   aliased.
#' @param nsim,seed As in [simulate()].
#' @param data2 Structural objects, as in [frm()].
#' @return A data frame with `nsim` columns of simulated responses,
#'   carrying the drawn parameters in `attr(., "pars")` when `prior`
#'   is used.
#' @examples
#' # power analysis: simulate from a design with chosen parameters
#' dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)), y = 0)
#' sims <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
#'                      newparams = list(b_Intercept = 1, b_x = 0.5,
#'                                       sigma = 0.7,
#'                                       sd_g__Intercept = 0.5),
#'                      nsim = 3, seed = 1)
#' head(sims)
#' # the same thing on the internal scale
#' par_template(bf(y ~ x + (1 | g)) + gaussian(), dd)   # the layout
#' sims2 <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
#'                       newparams = list(beta = c(1, 0.5),
#'                                        betad = log(0.7),
#'                                        theta = log(0.5)),
#'                       nsim = 3, seed = 1)
#' # prior-predictive draws
#' pp <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
#'                    prior = set_prior("normal(0, 1)", class = "b") +
#'                      set_prior("normal(0, 2)", class = "Intercept") +
#'                      set_prior("exponential(1)", class = "sd") +
#'                      set_prior("exponential(1)", class = "sigma"),
#'                    nsim = 4, seed = 1)
#' head(attr(pp, "pars"))
#' @export
frm_simulate <- function(formula, data, family = NULL, newparams = NULL,
                         prior = NULL, nsim = 1, seed = NULL,
                         data2 = list()) {
  prior <- as_priorlist(prior)
  # nsim reaches replicate() as a length; "invalid 'length' argument"
  # names neither this function nor the argument
  check_count(nsim, "nsim", min = 1L)
  if (!is.null(seed)) set.seed(seed)
  bform <- resolve_deferred_families(as_bform(formula, family), data)
  spec <- parse_spec(bform)
  frame <- assemble_frame(spec, data, data2 = data2)
  # An extension that installs its SIMULATOR at finalize (frmtmb.eam's
  # gddm()) was refused here for having none, while frm() on the same
  # model simulated fine, because the draw was taken from the family as
  # written rather than the one the likelihood would score.
  spec <- carry_finalized_responses(spec, frame)
  if (length(spec$responses) > 1L) {
    stop("frm_simulate() supports univariate models", call. = FALSE)
  }
  rspec <- spec$responses[[1L]]
  if (!sim_can(rspec$family)) {
    stop("frm_simulate(): family '", rspec$family[["family"]],
         "' has no simulator yet", sim_note(rspec$family), call. = FALSE)
  }
  if (is.null(newparams) && is.null(prior)) {
    stop("frm_simulate() needs newparams, prior, or both", call. = FALSE)
  }
  newparams <- newparams %||% list()
  if (length(newparams) && is.null(names(newparams))) {
    stop("newparams must be a named list", call. = FALSE)
  }

  est <- frame[["par_template"]]
  internal <- is_internal_newparams(newparams, est)
  np_internal <- if (internal) newparams else list()
  np_natural <- if (internal) list() else newparams
  slots <- nat_slots(frame, list(spec = spec, frame = frame,
                                 bform = bform))

  if (internal) {
    # a nonlinear parameter may be named after a template component and
    # still fit; see nl_start_collision_msg()
    nl_named <- unique(unlist(lapply(spec$responses,
                                     function(r) r$nlpars %||% character(0))))
    for (nm in names(np_internal)) {
      if (length(np_internal[[nm]]) != length(est[[nm]])) {
        stop("newparams$", nm, " must have length ", length(est[[nm]]),
             if (nm %in% nl_named) {
               paste0(". ", nl_start_collision_msg(nm, est, "newparams"))
             },
             call. = FALSE)
      }
      if (nm %in% nl_named) {
        warning(nl_start_collision_msg(nm, est, "newparams"), call. = FALSE)
      }
      est[[nm]][] <- np_internal[[nm]]
    }
  } else if (length(np_natural)) {
    check_natural_supported(frame)
    est <- apply_natural(est, frame, slots, np_natural)
  }

  entries <- list()
  labels <- character(0)
  if (!is.null(prior)) {
    shim0 <- list(spec = spec, frame = frame, estimates = est)
    entries <- resolve_prior_input(shim0, prior)$entries
    if (!length(entries)) {
      stop("The prior specification targets no parameter", call. = FALSE)
    }
    labels <- vapply(entries, function(e) {
      prior_entry_label(frame, slots, e)
    }, "")
  }
  if (length(entries) || length(np_natural)) {
    check_coverage(frame, slots, np_internal, np_natural, entries)
  }

  fixed_b <- "b" %in% names(np_internal)
  av <- frame[["aterm_values"]][[rspec$resp_name]]
  n <- frame[["n_obs"]]
  out <- vector("list", nsim)
  pars <- if (length(entries)) {
    matrix(NA_real_, nsim, length(entries),
           dimnames = list(NULL, labels))
  }
  for (s in seq_len(nsim)) {
    est_s <- est
    if (length(entries)) {
      dr <- draw_prior_pars(est_s, entries, labels)
      est_s <- dr$est
      pars[s, ] <- dr$vals
    }
    # a minimal fit-shaped object: eval_dpars, draw_b and the simulator
    # contract only touch spec / frame / estimates, so the shim carries
    # the group, sequence and residual-correlation structures the
    # structured simulators read straight off the assembled frame
    shim <- list(spec = spec, frame = frame, estimates = est_s)
    b_use <- if (is.null(est_s[["b"]])) {
      NULL
    } else if (fixed_b) {
      est_s[["b"]]
    } else {
      draw_b(shim)
    }
    dp <- with_cs_offsets(shim, rspec,
                          eval_dpars(shim, b = b_use))[[rspec$resp_name]]
    ex <- frame[["extra_names"]] %||% character(0)
    out[[s]] <- sim_draw(sim_context(
      shim, rspec, dp, aterms = av, n = n,
      extra = if (length(ex)) est_s[ex]))
    out[[s]] <- sim_restore_type(shim, rspec, out[[s]])
  }
  names(out) <- paste0("sim_", seq_len(nsim))
  out <- sim_as_data_frame(out)
  if (!is.null(pars)) attr(out, "pars") <- as.data.frame(pars)
  out
}
