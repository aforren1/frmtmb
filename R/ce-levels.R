# Which grouping level a conditional-effects grid row reads, block by
# block, and how a NEW level is drawn: shared by the fit method's
# bootstrap band and the draws method (frmtmb.sample).

#' The variables a prediction reads OUTSIDE its grouping terms: fixed
#' effects, smooths (a factor-smooth's own factor included), `mo()`,
#' `mi()`, `me()`, `gp()`, the covariates of a nonlinear body and the
#' `by` variable of a `gr(g, by = f)` term. A placeholder level may move
#' a grouping column only when nothing else reads it.
#'
#' @noRd
ce_locked_vars <- function(fit) {
  lps <- unlist(lapply(fit$frame[["linpreds"]] %||% list(), function(lp) {
    c(ce_lp_vars(lp),
      unlist(lapply(lp[["smooths"]] %||% list(), function(si) {
        smooth_pred_vars(si$sm)
      })),
      unlist(lapply(lp[["gps"]] %||% list(), function(gi) {
        unlist(lapply(gi$exprs, all.vars))
      })),
      names(lp[["data_list"]]),
      if (!is.null(lp[["nl_body"]])) all.vars(lp[["nl_body"]]))
  }))
  # a by variable picks WHICH block a row reads, so moving it would
  # move the row into another by-level's term
  by <- unlist(lapply(fit$frame[["re_blocks"]] %||% list(), function(bk) {
    all.vars(bk[["by"]][["expr"]])
  }))
  unique(c(lps, by))
}

#' The group-level blocks a display's prediction reads: every one under
#' `re_formula = NULL`, the terms a formula keeps, none under `NA`.
#' Smooth, `gp()` and `hsgp()` blocks are never groups here.
#'
#' @noRd
ce_kept_blocks <- function(fit, re_formula) {
  blocks <- fit$frame[["re_blocks"]] %||% list()
  grp <- which(vapply(blocks, function(bk) {
    !bk[["covstruct"]] %in% c("smooth", "gp", "hsgp")
  }, NA))
  if (!inherits(re_formula, "formula")) {
    return(if (is.null(re_formula)) grp else integer(0))
  }
  kp <- re_keep_plan(fit, re_formula, "conditional_effects()")
  if (identical(kp$kind, "all")) return(grp)
  if (!identical(kp$kind, "partial")) return(integer(0))
  kept <- unlist(lapply(seq_along(kp$named), function(k) {
    if (any(kp$keep[[k]])) kp$named[[k]]$block
  }))
  intersect(grp, unique(kept))
}

#' Which group-level block each grid row reads at an OBSERVED level and
#' which at a NEW one, and how a draw of the new ones is placed.
#'
#' brms's rule, per row and per block: `conditional_effects()` predicts
#' with `allow_new_levels = TRUE`, so a block reads its fitted effects
#' where every one of its grouping variables is set, in that row, to a
#' level the fit saw for that block, and a new level everywhere else:
#' a variable left unset (`NA`), a level the fit never saw (`"99"`), or
#' a combination a nested or interaction term never saw. A variable set
#' by `effects` counts as set, as one set by `conditions` does. A new
#' level is drawn once per draw (per bootstrap replicate on a fit) for
#' each block and each distinct written level, and shared by every row
#' and panel that names it.
#'
#' Two kinds of term read a row differently. A `gr(g, by = f)` term is
#' one block per level of `f`, and a row reads only the block of its own
#' `f`, so its new level is drawn with that by-level's covariance, as
#' brms draws it. An `mm(g1, g2)` term reads one level per member: each
#' distinct member value the fit never saw (an unset member is `NA`) is
#' one new level, drawn once, so two unset members are one new group
#' with their weights added, which is how the Wald band reads them too.
#'
#' The draw goes through the ordinary design: the rows of one pattern
#' of new blocks are predicted together, at a placeholder level of each
#' new block whose coefficients the draw overwrites. A placeholder moves
#' only grouping columns that nothing else reads (`ce_locked_vars()`)
#' and that no observed block of those rows reads, so a nested
#' `(1 | g / h)` row with `g` set and `h` unset keeps `g`'s own level
#' and draws a new `g:h` within it. A pattern no placeholder can serve
#' is refused by name.
#'
#' `status` gives, per block, whether any row reads it observed and
#' whether any reads it new, which is what the fit's bootstrap
#' simulation is chosen from.
#'
#' @noRd
ce_level_plan <- function(fit, grids, base, re_formula) {
  all_blocks <- fit$frame[["re_blocks"]] %||% list()
  blocks <- list()
  for (i in ce_kept_blocks(fit, re_formula)) {
    bk <- all_blocks[[i]]
    mmv <- unlist(lapply(bk[["components"]], function(cp) {
      cp[["mm"]][["gvars"]]
    }))
    mm <- length(mmv) > 0L
    gv <- if (mm) {
      mmv
    } else {
      unlist(strsplit(bk[["group_name"]] %||% "", ":", fixed = TRUE))
    }
    gv <- gv[nzchar(gv)]
    # a grouping built by a call has no column per variable to place a
    # level in; such a block is predicted as the design reads it and
    # gets no drawn level, as before
    if (!length(gv) || !all(gv %in% names(base))) next
    lv <- as.character(bk[["levels"]])
    parts <- if (mm) {
      matrix(lv, ncol = 1L)
    } else {
      do.call(rbind, strsplit(lv, ":", fixed = TRUE))
    }
    if (is.null(parts) || (!mm && ncol(parts) != length(gv))) next
    rr <- identical(bk[["covstruct"]], "rr")
    blocks[[length(blocks) + 1L]] <- list(
      id = i, bk = bk, gv = gv, parts = parts, rr = rr, mm = mm,
      # the levels a member may name: a by-split mm() block holds one
      # by-level's, and a member at another by-level's is still observed
      known = bk[["by"]][["group_levels"]] %||% lv,
      d = if (rr) bk[["rank"]] else bk[["dim"]],
      draw = rr || !bk[["covstruct"]] %in%
        c("gr_cov", "gr_prec", "car", "spde"))
  }
  locked <- ce_locked_vars(fit)
  status <- lapply(blocks, function(b) c(observed = FALSE, new = FALSE))
  plans <- lapply(grids, function(g) {
    nd <- g$nd
    n <- nrow(nd)
    if (!length(blocks)) {
      return(list(n = n, parts = list(list(rows = seq_len(n), nd = nd,
                                           new = list()))))
    }
    # per row and block: "" observed, "new:<level>" new, "-" not read
    st <- matrix("", n, length(blocks))
    for (k in seq_along(blocks)) {
      b <- blocks[[k]]
      vals <- lapply(b$gv, function(v) as.character(nd[[v]]))
      key <- do.call(paste, c(vals, sep = ":"))
      miss <- Reduce(`|`, lapply(vals, is.na))
      if (b$mm) {
        seen <- lapply(vals, function(v) !is.na(v) & v %in% b$known)
        # the members this block serves: all of them, or for one
        # by-level's block of mm(by = ) the members routed to it
        inb <- ce_mm_in_block(b, nd, vals, seen)
        obs <- Reduce(`&`, Map(function(s, i) s | !i, seen, inb))
        some <- Reduce(`|`, Map(`&`, seen, inb))
        read <- Reduce(`|`, inb)
        if (!is.null(b$bk[["by"]])) {
          # the key names only this block's members, so that a member
          # routed elsewhere does not split this block's rows
          key <- do.call(paste, c(Map(function(v, i) {
            ifelse(i, v, "-")
          }, vals, inb), sep = ":"))
        }
      } else {
        obs <- !miss & key %in% b$bk[["levels"]]
        some <- obs
        read <- ce_by_reads(b$bk, nd)
      }
      st[, k] <- ifelse(!read, "-", ifelse(obs, "", paste0("new:", key)))
      status[[k]] <<- status[[k]] | c(any(read & some), any(read & !obs))
    }
    sig <- apply(st, 1L, paste, collapse = "\r")
    groups <- split(seq_len(n), factor(sig, levels = unique(sig)))
    parts <- lapply(groups, ce_plan_part, nd = nd, st = st,
                    blocks = blocks, locked = locked, base = base)
    list(n = n, parts = unname(parts))
  })
  list(grids = plans, blocks = blocks, status = status)
}

#' Which members of an `mm()` term each row routes to one block, as a
#' list of logical vectors, one per member.
#'
#' Without a `by` variable a block serves every member. With one, the
#' term is one block per by-level, and the design routes a member the
#' way `by_route_rows()` does: a member at an observed group reads the
#' block that holds that group, and a member at a new level reads the
#' block of its own by-value in the row. A by-value the grid cannot
#' evaluate leaves every member in every block, as before.
#'
#' @noRd
ce_mm_in_block <- function(b, nd, vals, seen) {
  by <- b$bk[["by"]]
  all_in <- lapply(vals, function(v) rep(TRUE, length(v)))
  if (is.null(by)) return(all_in)
  byv <- tryCatch(
    as.matrix(eval(by[["expr"]], nd, globalenv())),
    error = function(e) NULL)
  if (is.null(byv) || nrow(byv) != nrow(nd) || ncol(byv) != length(vals)) {
    return(all_in)
  }
  lv <- as.character(b$bk[["levels"]])
  lapply(seq_along(vals), function(k) {
    bk_by <- as.character(byv[, k])
    ifelse(seen[[k]], vals[[k]] %in% lv,
           !is.na(bk_by) & bk_by == by[["level"]])
  })
}

#' The label the design gives a row for a block: the grouping
#' expression evaluated on the row, as the newdata design evaluates it
#' (`g:h` is the interaction of two factors, `"1:1"`). A row whose
#' expression cannot be evaluated here falls back to the joined values.
#'
#' @noRd
ce_design_label <- function(b, ndp) {
  bar <- b$bk[["components"]][[1L]][["bar"]]
  lab <- if (!is.null(bar)) {
    tryCatch(as.character(group_values(bar[[3L]], ndp[1L, , drop = FALSE],
                               globalenv())),
             error = function(e) NULL)
  }
  if (length(lab) != 1L) {
    vals <- vapply(b$gv, function(v) as.character(ndp[[v]][1L]), "")
    lab <- if (anyNA(vals)) NA_character_ else paste(vals, collapse = ":")
  }
  lab
}

#' Which rows read a block: all of them, except that the block of one
#' level of a `gr(g, by = f)` term is read only where `f` has that
#' level.
#'
#' @noRd
ce_by_reads <- function(bk, nd) {
  by <- bk[["by"]]
  if (is.null(by)) return(rep(TRUE, nrow(nd)))
  v <- tryCatch(as.character(eval(by[["expr"]], nd, globalenv())),
                error = function(e) NULL)
  # a by value the grid cannot evaluate is left to the design, which
  # reads every block as before
  if (!length(v)) return(rep(TRUE, nrow(nd)))
  v <- rep_len(v, nrow(nd))
  !is.na(v) & v == by[["level"]]
}

#' The refusal for a new level no placeholder can serve.
#'
#' @noRd
ce_plan_refuse <- function(b, nd, r1) {
  at <- paste0(b$gv, " = ", vapply(b$gv, function(v) {
    as.character(nd[[v]][r1])
  }, ""), collapse = ", ")
  frm_stop("conditional_effects() cannot draw a new level of the ",
           "group-level term (", b$bk[["term_label"]], ") at the grid ",
           "rows with ", at, ": to place one it would move a column ",
           "that another term or an observed group reads, to a level ",
           "the fit saw with the values the grid holds. Set the ",
           "variables in conditions to levels the fit saw together, ",
           "or use re_formula = NA", call. = FALSE)
}

#' One pattern of new blocks: the rows, the placeholder levels and the
#' sub-grid they are predicted on.
#'
#' @noRd
ce_plan_part <- function(rows, nd, st, blocks, locked, base) {
  r1 <- rows[1L]
  newk <- which(startsWith(st[r1, ], "new:"))
  ndp <- nd[rows, , drop = FALSE]
  if (!length(newk)) return(list(rows = rows, nd = ndp, new = list()))
  # a column an observed block of these rows reads keeps its value, and
  # so does an observed member of a multi-membership term; a nested
  # term's shorter parent is placed after it, so it takes the level the
  # longer term chose
  obsk <- which(st[r1, ] == "")
  # an observed multi-membership block holds only its observed members:
  # with a by variable, a member routed to another by-level's block is
  # not read here and may still need a placeholder there
  held <- union(locked, unlist(lapply(blocks[obsk], function(b) {
    if (!b$mm) return(b$gv)
    vals <- vapply(b$gv, function(v) as.character(nd[[v]][r1]), "")
    b$gv[!is.na(vals) & vals %in% b$known]
  })))
  mmk <- newk[vapply(blocks[newk], `[[`, NA, "mm")]
  for (k in mmk) {
    b <- blocks[[k]]
    vals <- vapply(b$gv, function(v) as.character(nd[[v]][r1]), "")
    held <- union(held, b$gv[!is.na(vals) & vals %in% b$known])
  }
  given <- list()
  moved <- character(0)
  new <- list()
  regk <- setdiff(newk, mmk)
  ord <- regk[order(-vapply(blocks[regk], function(b) length(b$gv), 1L))]
  for (k in ord) {
    b <- blocks[[k]]
    ok <- rep(TRUE, nrow(b$parts))
    for (j in seq_along(b$gv)) {
      v <- b$gv[j]
      want <- if (!is.null(given[[v]])) {
        given[[v]]
      } else if (v %in% held) {
        as.character(nd[[v]][r1])
      }
      if (!is.null(want)) ok <- ok & !is.na(want) & b$parts[, j] == want
    }
    if (!any(ok)) {
      # no level of this block agrees with the columns the rows must
      # keep: crossed (1 | g) + (1 | h) + (1 | g:h) with g and h at
      # observed levels never seen together. The columns stay, and in
      # this part's copy of the fit the block's first level is renamed
      # to the label the rows carry (ce_plan_eval()), so the design
      # reads the drawn effects there and nothing else moves. A row with
      # a grouping variable unset has no label to rename to, and a
      # label of NA would read NA wherever the variable is also a
      # predictor (y ~ trt + (1 | trt:subj) with nothing set), so it
      # keeps 0.66.0's refusal
      if (anyNA(vapply(b$gv, function(v) as.character(nd[[v]][r1]), ""))) {
        ce_plan_refuse(b, nd, r1)
      }
      for (v in b$gv) {
        if (is.null(given[[v]])) given[[v]] <- as.character(nd[[v]][r1])
      }
      new[[length(new) + 1L]] <- list(
        block = k, key = st[r1, k], idx = b$bk[["b_idx"]][seq_len(b$d)],
        relabel = 1L)
      next
    }
    lv <- which(ok)[1L]
    for (j in seq_along(b$gv)) {
      v <- b$gv[j]
      if (is.null(given[[v]])) {
        given[[v]] <- b$parts[lv, j]
        if (!v %in% held) moved <- c(moved, v)
      }
    }
    new[[length(new) + 1L]] <- list(
      block = k, key = st[r1, k],
      idx = b$bk[["b_idx"]][(lv - 1L) * b$d + seq_len(b$d)])
  }
  # a multi-membership term: each distinct new member value gets its
  # own placeholder, one the row's observed members do not name, so the
  # members sharing a value share its draw and their weights add up.
  # With a by variable the term is one block per by-level, and a new
  # member is placed in the block of its own by-value, at one of that
  # block's levels: the other by-levels' blocks then read that level as
  # an observed group of another by-level, which they skip
  for (k in mmk) {
    b <- blocks[[k]]
    vals <- vapply(b$gv, function(v) as.character(nd[[v]][r1]), "")
    seen <- !is.na(vals) & vals %in% b$known
    inb <- vapply(ce_mm_in_block(b, nd[r1, , drop = FALSE], as.list(vals),
                                 as.list(seen)), `[`, NA, 1L)
    isnew <- !seen & inb
    if (!any(isnew)) next
    u <- unique(vals[isnew])
    mt <- match(vals, u)
    free <- setdiff(b$parts[, 1L], vals[seen])
    for (v in b$gv[isnew]) {
      if (is.factor(base[[v]])) free <- intersect(free, levels(base[[v]]))
    }
    if (any(b$gv[isnew] %in% held)) ce_plan_refuse(b, nd, r1)
    # a member a plain new term already placed (`(1 | g1)` beside
    # `mm(g1, g2)`) keeps that placeholder: the two terms' coefficients
    # are separate slots, so one column value serves both draws
    place <- rep(NA_character_, length(u))
    for (j in which(isnew)) {
      gj <- given[[b$gv[j]]]
      if (is.null(gj)) next
      if (!is.na(place[mt[j]]) && place[mt[j]] != gj) {
        ce_plan_refuse(b, nd, r1)
      }
      place[mt[j]] <- gj
    }
    pre <- place[!is.na(place)]
    if (!all(pre %in% free) || anyDuplicated(pre)) ce_plan_refuse(b, nd, r1)
    need <- which(is.na(place))
    free <- setdiff(free, pre)
    if (length(free) < length(need)) ce_plan_refuse(b, nd, r1)
    place[need] <- free[seq_along(need)]
    for (j in which(isnew)) {
      if (is.null(given[[b$gv[j]]])) {
        given[[b$gv[j]]] <- place[mt[j]]
        moved <- c(moved, b$gv[j])
      }
    }
    for (j in seq_along(u)) {
      lv <- match(place[j], b$parts[, 1L])
      new[[length(new) + 1L]] <- list(
        block = k, key = paste0("new:", u[j]),
        idx = b$bk[["b_idx"]][(lv - 1L) * b$d + seq_len(b$d)])
    }
  }
  for (v in unique(moved)) {
    col <- base[[v]]
    val <- given[[v]]
    ndp[[v]] <- rep(if (is.factor(col)) {
      factor(val, levels = levels(col))
    } else if (is.numeric(col)) {
      as.numeric(val)
    } else {
      val
    }, length.out = length(rows))
  }
  # a renamed level takes the label the design gives these rows, read
  # after every move so that it is the label the prediction will see
  for (i in seq_along(new)) {
    if (!is.null(new[[i]]$relabel)) {
      new[[i]]$label <- ce_design_label(blocks[[new[[i]]$block]], ndp)
    }
  }
  list(rows = rows, nd = ndp, new = new)
}

#' One grid's display values under one parameter vector, each new level
#' drawn from THAT vector's covariance. `cache`, an environment, shares a
#' draw between the rows and panels that name the same new level.
#'
#' @noRd
ce_plan_eval <- function(f, plan, gi, cache, categorical, resp, dpar,
                         re_form, anl) {
  gp <- plan$grids[[gi]]
  if (length(gp$parts) == 1L && !length(gp$parts[[1L]]$new)) {
    return(ce_boot_one(f, gp$parts[[1L]]$nd, categorical, resp, dpar,
                       re_form, anl))
  }
  out <- NULL
  for (p in gp$parts) {
    fp <- f
    if (length(p$new)) {
      b <- fp$estimates[["b"]]
      for (e in p$new) {
        ck <- paste0(e$block, "\r", e$key)
        if (is.null(cache[[ck]])) {
          cache[[ck]] <- ce_plan_draw(fp, plan$blocks[[e$block]])
        }
        b[e$idx] <- cache[[ck]]
        if (!is.null(e$relabel)) {
          # the copy's level is renamed to the rows' label; the copy
          # lives for this part only, so the fit itself is untouched
          id <- plan$blocks[[e$block]]$id
          lv <- as.character(fp$frame[["re_blocks"]][[id]][["levels"]])
          lv[e$relabel] <- e$label
          fp$frame[["re_blocks"]][[id]][["levels"]] <- lv
        }
      }
      fp$estimates[["b"]] <- b
    }
    v <- ce_boot_one(fp, p$nd, categorical, resp, dpar, re_form, anl)
    K <- length(v) %/% length(p$rows)
    if (is.null(out)) out <- matrix(NA_real_, gp$n, K)
    out[p$rows, ] <- matrix(v, length(p$rows), K)
  }
  as.vector(out)
}

#' One new level's effects for one block, from the covariance the
#' parameters of `f` imply; zero for a block whose levels ARE its
#' structure (`gr_cov`, `gr_prec`, `car`, `spde`), which has no marginal
#' law for an unseen level, as the Wald band assumes too.
#'
#' @noRd
ce_plan_draw <- function(f, b) {
  if (!isTRUE(b$draw)) return(numeric(b$d))
  if (b$rr) return(stats::rnorm(b$d))
  bk <- b$bk
  th <- f$estimates[["theta"]]
  S <- covstruct_registry[[bk[["covstruct"]]]]$vcov(th[bk[["theta_idx"]]],
                                                    bk)
  if (is_student_block(bk)) S <- S * student_var_factor(bk[["dist_nu"]])
  mvn_draw_cov(as.matrix(S))
}

#' Whether a plan draws any new level at all.
#'
#' @noRd
ce_plan_has_new <- function(plan) {
  any(vapply(plan$grids, function(gp) {
    any(vapply(gp$parts, function(p) length(p$new) > 0L, NA))
  }, NA))
}

#' The blocks a bootstrap simulation holds at their fitted effects: the
#' ones some grid row reads at an observed level. A block read observed
#' in one row and new in another would need both kinds of simulation at
#' once, which one bootstrap cannot give, so that call is refused by
#' name.
#'
#' @noRd
ce_plan_kept <- function(plan) {
  kept <- integer(0)
  for (k in seq_along(plan$blocks)) {
    s <- plan$status[[k]]
    if (s[["observed"]] && s[["new"]]) {
      frm_stop("band = \"boot\" cannot cover this call: the group-level ",
               "term (", plan$blocks[[k]]$bk[["term_label"]], ") is read ",
               "at an observed level in some grid rows and at a new level ",
               "in others (or, for a multi-membership term, through ",
               "different members of one row). The observed levels need ",
               "refits that keep that term's fitted effects and the new ",
               "ones refits that redraw them, and one bootstrap cannot do ",
               "both. Use band = ",
               "\"wald\", or call conditional_effects() once per kind of ",
               "row", call. = FALSE)
    }
    if (s[["observed"]]) kept <- c(kept, plan$blocks[[k]]$id)
  }
  kept
}

#' A comparable key of a level plan: which rows of which grid draw which
#' new level of which block. The draws of a bootstrap depend on it, so
#' reuse compares it.
#'
#' @noRd
ce_plan_key <- function(plan) {
  # a plan that draws nothing is the population's key: which group the
  # draws belong to is then said by the simulation (`kept`), not here
  if (is.null(plan) || !ce_plan_has_new(plan)) return(list())
  lapply(plan$grids, function(gp) {
    lapply(gp$parts, function(p) {
      list(rows = p$rows, new = lapply(p$new, function(e) {
        list(block = plan$blocks[[e$block]]$id, key = e$key, idx = e$idx)
      }))
    })
  })
}
