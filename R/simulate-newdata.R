# What `simulate()` draws for the group-level terms, and how it draws at
# `newdata`.
#
# `re_formula` means in simulate() what it means in predict(): `NULL`
# keeps every group-level term, `NA`, `~0` and `~1` keep none, and a
# one-sided formula keeps the terms it names. What differs is what
# "not kept" does to a DRAW. predict() and fitted() report an average
# group, so a dropped term contributes nothing. A simulated response
# needs a group, so a dropped term is REDRAWN from its estimated
# distribution each replicate, which is lme4's unconditional simulation
# and what `re_formula = NA` has always done here. Through 0.62.0 any
# formula meant "keep everything", so `~1` gave the conditional draws
# where predict(re_formula = ~1) gave the population prediction.
#
# The same resolution (`re_keep_plan()`) decides both, so a formula
# that predict() refuses is refused here with the same words.

#' Which group-level coefficients one replicate redraws.
#'
#' Returns `list(redraw, cond, blocks, every)`. `redraw` holds the `b`
#' positions taken from a fresh `draw_b()`; `cond` holds one entry per
#' block where a formula keeps some columns and drops others, whose
#' dropped columns are drawn given the kept ones at their estimates;
#' `blocks` is every block touched; `every` is `TRUE` when no
#' group-level term is kept at all and none is a factor smooth, which
#' is when an unseen level at `newdata` is simply one more fresh level:
#' the newdata design gives an unseen level of a factor smooth the
#' population curve, not a fresh one, so there the level is left to
#' `allow_new_levels`, as in `predict()`.
#'
#' `NA` redraws exactly what predict(re_formula = NA) drops: every block
#' but a POPULATION smooth, `gp()` or `hsgp()` curve. Those are
#' penalized coefficients of a population-level term, not group
#' effects. Redrawing them from the smoothing prior, as every release
#' through 0.62.0 did, replaced the fitted curve with a random one: on
#' `y ~ s(x)` the draws spread with an sd of 2.13 around a curve whose
#' residual sd is 0.28 (dev/simnewdata-log/probe2-base.txt).
#'
#' `smooths = TRUE` restores that for `NA`, `~0` and `~1`: every block
#' is redrawn. It is `frm_bootstrap()`'s setting (user decision,
#' 2026-09-24), a whole-model parametric bootstrap in the manner of
#' lme4's `bootMer(use.u = FALSE)` with the smooths treated as random
#' effects; `simulate()` does not offer it.
#'
#' @noRd
sim_re_plan <- function(fit, re_formula, smooths = FALSE) {
  none <- list(redraw = integer(0), cond = list(), blocks = integer(0))
  blocks <- fit$frame[["re_blocks"]] %||% list()
  if (!length(blocks) || is.null(re_formula)) return(none)
  kp <- re_keep_plan(fit, re_formula, "simulate()")
  if (identical(kp$kind, "all")) return(none)
  if (kp$kind %in% c("asis", "none")) {
    # asis here is NA: check_re_form() has refused anything else
    ids <- if (smooths) seq_along(blocks) else sim_group_block_ids(fit)
    fs <- any(vapply(blocks[ids], function(bk) {
      bk[["covstruct"]] %in% c("smooth", "gp", "hsgp")
    }, NA))
    return(list(redraw = unlist(lapply(blocks[ids], `[[`, "b_idx")),
                cond = list(), blocks = ids, every = !fs))
  }
  # a partial formula: per block, which columns stay
  colkeep <- lapply(blocks, function(bk) rep(TRUE, bk[["dim"]]))
  for (k in seq_along(kp$named)) {
    nc <- kp$named[[k]]
    cp <- blocks[[nc$block]][["components"]][[nc$comp]]
    colkeep[[nc$block]][cp[["offset"]] + which(!kp$keep[[k]])] <- FALSE
  }
  redraw <- integer(0)
  cond <- list()
  ids <- integer(0)
  for (i in seq_along(blocks)) {
    ck <- colkeep[[i]]
    if (all(ck)) next
    bk <- blocks[[i]]
    ids <- c(ids, i)
    if (!any(ck)) {
      redraw <- c(redraw, bk[["b_idx"]])
      next
    }
    if (bk[["covstruct"]] %in% c("rr", "gr_cov", "gr_prec", "car", "spde") ||
          is_student_block(bk)) {
      frm_stop("simulate(): re_formula keeps some columns of the ",
               "group-level term ", bk[["term_label"]] %||% bk[["covstruct"]],
               " and drops others. A dropped column is redrawn given the ",
               "kept ones, which needs the conditional law of one set of ",
               "columns given another; this term's '", bk[["covstruct"]],
               "' structure", if (is_student_block(bk)) " (Student-t)",
               " has no such law here. Keep or drop the whole term",
               call. = FALSE)
    }
    cond[[length(cond) + 1L]] <- list(block = i, keep = ck)
  }
  list(redraw = redraw, cond = cond, blocks = ids)
}

#' The blocks `re_formula = NA` removes from a prediction: every block
#' except the population smooth, `gp()` and `hsgp()` blocks
#' (`lp_eta_design()`'s rule, from the same `smooth_group_block_ids()`).
#'
#' @noRd
sim_group_block_ids <- function(fit) {
  blocks <- fit$frame[["re_blocks"]]
  grp_sm <- unique(unlist(lapply(fit$frame[["linpreds"]],
                                 smooth_group_block_ids)))
  which(!vapply(seq_along(blocks), function(i) {
    blocks[[i]][["covstruct"]] %in% c("smooth", "gp", "hsgp") &&
      !i %in% grp_sm
  }, NA))
}

#' One replicate's `b`: the estimates, with the planned positions
#' redrawn. `draw_b()` runs only when something is redrawn whole, and
#' then it runs once over every block, so a plan that redraws every
#' block consumes the caller's stream exactly as 0.62.0's
#' `re_formula = NA` did and returns the same vector.
#'
#' @noRd
sim_draw_b <- function(fit, plan) {
  b <- fit$estimates[["b"]]
  if (length(plan$redraw)) {
    bd <- draw_b(fit)
    b[plan$redraw] <- bd[plan$redraw]
  }
  for (cd in plan$cond) b <- sim_b_conditional(fit, cd, b)
  b
}

#' Draw a block's dropped columns given its kept columns at their
#' estimates, level by level: `b_D | b_K ~ N(V_DK V_KK^-1 b_K,
#' V_DD - V_DK V_KK^-1 V_KD)`. Drawing them from their marginal instead
#' would discard the fitted correlation between an intercept that is
#' kept and a slope that is not.
#'
#' @noRd
sim_b_conditional <- function(fit, cd, b) {
  bk <- fit$frame[["re_blocks"]][[cd$block]]
  th <- fit$estimates[["theta"]]
  V <- covstruct_registry[[bk[["covstruct"]]]]$vcov(th[bk[["theta_idx"]]],
                                                    bk)
  V <- as.matrix(V)
  K <- which(cd$keep)
  R <- which(!cd$keep)
  D <- bk[["dim"]]
  B <- matrix(b[bk[["b_idx"]]], nrow = D)
  A <- V[R, K, drop = FALSE] %*% solve(V[K, K, drop = FALSE])
  Vc <- V[R, R, drop = FALSE] - A %*% V[K, R, drop = FALSE]
  Vc <- (Vc + t(Vc)) / 2
  M <- A %*% B[K, , drop = FALSE]
  for (l in seq_len(ncol(B))) B[R, l] <- M[, l] + mvn_draw_cov(Vc)
  b[bk[["b_idx"]]] <- as.vector(B)
  b
}

#' Refuse a `newdata` a structured FAMILY cannot draw at.
#'
#' A family whose simulator is a structure (`sim_ctx`) draws whole
#' sequences or groups out of the frame it was fitted on:
#' `mixture(groups = )` takes one class per fitted group, a hidden
#' Markov family (frmtmb.latent) walks the fitted state sequence, and a
#' learning family (frmtmb.learn) walks the fitted trial sequence.
#' `mixture_mvn()` draws row by row, but it goes through the same slot
#' and is refused with the rest, as `predict(newdata = )` and
#' `posterior_predict(newdata = )` refuse it. A residual correlation
#' term is not a family structure and is rebuilt on the new rows
#' instead (`autocor_for_newdata()`).
#'
#' @noRd
sim_newdata_refuse_structure <- function(fam) {
  if (is.null(fam_sim_ctx(fam))) return(invisible(NULL))
  frm_stop("simulate(newdata = ) is not supported for family '",
           fam[["family"]], "': its draw is structured, and the structure ",
           "(a group-level latent class, a hidden state sequence, a trial ",
           "sequence) indexes the rows the model was fitted on. Drop ",
           "newdata to simulate those rows", call. = FALSE,
           package = frm_family_package(fam))
}

#' Refuse a `newdata` that lacks a grouping column the draw reads.
#'
#' A redrawn group effect is shared by the rows of one level, and a kept
#' one enters at that level's estimate, so the column says which rows
#' are which. Without it every row would be its own level, and the
#' within-group correlation would vanish with nothing said. Under a
#' kept term `allow_new_levels = TRUE` asks for exactly that, as it
#' does in `predict()`; a redrawn term has no such reading.
#'
#' A grouping variable that is not a column of newdata but exists in the
#' formula environment is read from there, because `predict()` and
#' `frm_linpred()` read it from there (`fill_new_group_vars()`, which
#' follows `model.frame()`; measured in `edge-predict-base.txt`). One of
#' the wrong length is refused here by name; in `predict()` it dies on
#' "non-conformable arrays".
#'
#' @noRd
sim_newdata_group_cols <- function(fit, rspec, newdata, plan,
                                   allow_new_levels) {
  env <- rspec$formula_env
  gvars <- function(ids) {
    out <- character(0)
    for (bk in fit$frame[["re_blocks"]][ids]) {
      # an spde block's levels are mesh nodes, read by its own design
      if (bk[["covstruct"]] == "spde") next
      for (comp in bk[["components"]]) {
        bar <- comp[["bar"]]
        out <- c(out, if (!is.null(comp[["mm"]])) comp[["mm"]]$gvars else
          if (is.call(bar)) all.vars(bar[[3L]]))
      }
    }
    setdiff(unique(out), names(newdata))
  }
  from_env <- gvars(sim_group_block_ids(fit))
  from_env <- from_env[vapply(from_env, exists, NA, envir = env)]
  for (v in from_env) {
    len <- NROW(get(v, envir = env))
    if (len != nrow(newdata)) {
      frm_stop("simulate(newdata = ): the grouping factor `", v, "` is not ",
               "a column of newdata, and the `", v, "` the formula's ",
               "environment holds has ", len, " value(s) for newdata's ",
               nrow(newdata), " rows. Add the column to newdata",
               call. = FALSE)
    }
  }
  ids <- if (allow_new_levels) plan$blocks else {
    union(plan$blocks, sim_group_block_ids(fit))
  }
  need <- gvars(ids)
  need <- need[!vapply(need, exists, NA, envir = env)]
  if (!length(need)) return(invisible(NULL))
  frm_stop("simulate(newdata = ): newdata has no column ",
           paste0("`", need, "`", collapse = ", "), ", the grouping ",
           "factor of a group-level term. A kept term enters at each ",
           "level's estimate and a redrawn one shares one draw across the ",
           "rows of a level, so the column says which rows are which. ",
           "allow_new_levels = TRUE treats every row as an unseen level of ",
           "a kept term", call. = FALSE)
}

#' The parts of the newdata design that do not change between
#' replicates, built once per call.
#'
#' Only `b` differs from one replicate to the next, so each linear
#' predictor's design pieces (`lp_eta_design()`: the fixed-effect
#' product, the offset, the group-level and smooth parts, the
#' estimability mask) and an ordinal `cs()` term's offsets are built
#' here once. Rebuilding them per replicate cost 6.2 ms against 0.64 ms
#' on the fitted rows (500 rows, `sigma ~ x`; `timing-before.txt`). A
#' nonlinear predictor, whose value is not a linear function of `b`, is
#' left to `frm_linpred()` in each replicate.
#'
#' @noRd
sim_newdata_design <- function(fit, rspec, newdata, allow_new_levels) {
  rn <- rspec$resp_name
  est <- fit$estimates
  per <- list()
  for (dnm in names(rspec$dpars)) {
    lp <- fit$frame[["linpreds"]][[linpred_key(rn, dnm)]]
    if (is.null(lp) || !is.null(lp[["nl_body"]])) {
      per[[dnm]] <- list(lp = lp, ed = NULL)
      next
    }
    ed <- lp_eta_design(fit, lp, newdata, TRUE, allow_new_levels)
    # the fixed part, formed exactly as lp_eta_design() forms it
    xb <- drop(as.matrix(ed[["X"]] %*% est[[lp[["par"]]]][lp[["idx"]]]))
    per[[dnm]] <- list(lp = lp, ed = ed, xb = xb)
  }
  n <- nrow(newdata)
  CS <- NULL
  for (lp in fit$frame[["linpreds"]]) {
    if (!identical(lp[["resp"]], rn) || !length(lp[["cs"]] %||% list())) {
      next
    }
    for (ct in ord_cs_values(fit, lp, newdata, n)) {
      CS <- (CS %||% 0) + outer(ct$vals, est[[ct$par]])
    }
  }
  list(per = per, cs = CS, n = n)
}

#' Every dpar of one response at `newdata` for one replicate's `b`, on
#' the scale the simulator consumes (the link inverse,
#' `dpars_natural()`'s scale, so a mixture weight is not passed through
#' its reporting softmax), with an unseen level's drawn effect added on
#' the link scale. The sums are `lp_eta_design()`'s, in its order, so
#' the draws are bitwise those of `frm_linpred()` at the same `b`.
#'
#' @noRd
sim_newdata_dpars <- function(fit, rspec, newdata, allow_new_levels, dz,
                              off = NULL) {
  rn <- rspec$resp_name
  cvec <- coef_b(fit)
  dp <- list()
  for (dnm in names(rspec$dpars)) {
    pd <- dz$per[[dnm]]
    ed <- pd$ed
    if (is.null(ed)) {
      eta <- as.vector(frm_linpred(fit, newdata = newdata, dpar = dnm,
                                   resp = rn, type = "link",
                                   allow_new_levels = allow_new_levels))
    } else {
      eta <- pd$xb
      if (length(ed[["re_parts"]])) {
        eta <- eta + re_eta(ed[["re_parts"]], cvec, ed[["n"]])
      }
      if (length(ed[["sm_parts"]])) {
        eta <- eta + sm_eta(ed[["sm_parts"]], cvec)
      }
      if (!is.null(ed[["off"]])) eta <- eta + ed[["off"]]
      if (any(ed[["nonest"]])) eta[ed[["nonest"]]] <- NA_real_
      eta <- as.vector(eta)
    }
    if (!is.null(off[[dnm]])) eta <- eta + off[[dnm]]
    dp[[dnm]] <- as.vector(pd$lp[["link"]]$linkinv(eta))
  }
  if (!is.null(dz$cs)) dp[[".cs"]] <- dz$cs
  dp
}

#' A residual correlation term rebuilt on the rows of `newdata`.
#'
#' The correlation matrix is the fitted one, over the fitted time
#' levels; what changes is which rows form a group and at which time
#' each row sits. So the time levels are matched against the FITTED set,
#' read the way `autocor_time_index()` read them at fit time, and a
#' time the fit never saw is refused: a lag is counted in levels, and a
#' new level has no place in that count. With no time variable the
#' time is the row's position within its group, as at fit time. The
#' groups themselves may be new, because the residual draw has no
#' per-group parameter.
#'
#' @noRd
autocor_for_newdata <- function(fit, ac, rspec, newdata) {
  env <- rspec$formula_env
  n <- nrow(newdata)
  gv <- if (is.null(ac[["gr_expr"]])) NULL else {
    autocor_gr_value(ac, newdata, env)
  }
  gidx <- if (is.null(gv)) rep(1L, n) else as.integer(factor(gv))
  if (is.null(ac[["time_expr"]])) {
    tidx <- integer(n)
    for (g in split(seq_len(n), gidx)) tidx[g] <- seq_along(g)
  } else {
    tv <- eval(ac[["time_expr"]], newdata, env)
    tf <- eval(ac[["time_expr"]], fit$frame[["data_frame"]], env)
    if (is.factor(tf)) {
      tidx <- match(as.character(tv), levels(droplevels(tf)))
    } else if (is.numeric(tf)) {
      tidx <- match(as.numeric(tv), sort(unique(as.numeric(tf))))
    } else {
      tidx <- match(as.character(tv), sort(unique(as.character(tf))))
    }
  }
  bad <- is.na(tidx) | tidx > ac[["d"]]
  if (any(bad)) {
    frm_stop("simulate(newdata = ): the residual correlation term ",
             ac[["label"]], " counts lags in the fitted time levels (",
             ac[["d"]], " of them), and newdata puts ", sum(bad),
             " row(s) at a time outside that set",
             if (!is.null(ac[["time_expr"]])) {
               paste0(" (", deparse1(ac[["time_expr"]]), " = ",
                      paste(unique(format(tv[bad])), collapse = ", "), ")")
             } else " (more rows in a group than time levels)",
             call. = FALSE)
  }
  key <- paste(gidx, tidx, sep = "\r")
  if (anyDuplicated(key)) {
    frm_stop("simulate(newdata = ): the residual correlation term ",
             ac[["label"]], " needs each group's rows at distinct times, ",
             "and newdata has ", sum(duplicated(key)), " row(s) that ",
             "repeat a time within their group", call. = FALSE)
  }
  by_g <- split(seq_len(n), gidx)
  pkeys <- vapply(by_g, function(rows) {
    paste(sort(tidx[rows]), collapse = ",")
  }, "")
  ac[["patterns"]] <- unname(lapply(split(seq_along(by_g), pkeys),
                                    function(gs) {
    tset <- sort(tidx[by_g[[gs[1L]]]])
    k <- length(tset)
    rows <- vapply(gs, function(g) {
      r <- by_g[[g]]
      r[order(tidx[r])]
    }, integer(k))
    list(k = k, G = length(gs), rows = as.vector(rows),
         gather = as.vector(outer(tset, tset, function(a, b) {
           (b - 1L) * ac[["d"]] + a
         })))
  }))
  ac[["n_groups"]] <- length(by_g)
  ac
}

#' One replicate at `newdata` when some rows have a parameter that is
#' not finite: the family draws the rows that are, and the others are
#' `NA`, as `predict()` leaves them.
#'
#' A residual correlation draws a group's rows jointly, so its structure
#' was built on the rows finite at the estimates (`ok0`); a replicate
#' whose drawn effects lose a further row cannot be drawn on that
#' structure and is `NA` whole, as `predict()`'s structured branch
#' leaves it.
#'
#' @noRd
sim_newdata_draw_rows <- function(fit, rspec, dp, av, ac, ok, ok0) {
  n <- length(ok)
  ys <- NULL
  if (any(ok) && (is.null(ac) || identical(ok, ok0))) {
    ctx <- sim_context(fit, rspec, lapply(dp, subset_rows, keep = ok),
                       aterms = lapply(av, subset_rows, keep = ok),
                       n = sum(ok), extra = fit_extras(fit))
    ctx[["autocor"]] <- ac
    ys <- sim_draw(ctx)
  }
  if (is.matrix(ys)) {
    out <- matrix(NA_real_, n, ncol(ys), dimnames = list(NULL, colnames(ys)))
    out[ok, ] <- ys
    return(out)
  }
  out <- rep(NA_real_, n)
  if (!is.null(ys)) out[ok] <- ys
  out
}

#' @noRd
sim_newdata_warn_rows <- function(rows) {
  frm_warning("simulate(newdata = ): ", length(rows), " row(s) of newdata ",
              "have a distributional parameter that is not finite at the ",
              "estimates, usually from a missing covariate: ",
              paste(utils::head(rows, 10L), collapse = ", "),
              if (length(rows) > 10L) ", ...", ". Their draws are NA, as ",
              "predict() leaves them", call. = FALSE)
}

#' @noRd
sim_newdata_warn_draws <- function(n_lost, nsim) {
  frm_warning("simulate(newdata = ): a redrawn group effect made a ",
              "distributional parameter not finite in ", n_lost, " cell(s) ",
              "over ", nsim, " replicate(s). Those cells are NA; the other ",
              "cells of each replicate stand, except under a residual ",
              "correlation term, whose group is drawn whole", call. = FALSE)
}
