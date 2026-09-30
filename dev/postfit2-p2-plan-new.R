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
        seen <- lapply(vals, function(v) {
          !is.na(v) & v %in% b$bk[["levels"]]
        })
        obs <- Reduce(`&`, seen)
        some <- Reduce(`|`, seen)
      } else {
        obs <- !miss & key %in% b$bk[["levels"]]
        some <- obs
      }
      read <- ce_by_reads(b$bk, nd)
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
  held <- union(locked, unlist(lapply(blocks[obsk], `[[`, "gv")))
  mmk <- newk[vapply(blocks[newk], `[[`, NA, "mm")]
  for (k in mmk) {
    b <- blocks[[k]]
    vals <- vapply(b$gv, function(v) as.character(nd[[v]][r1]), "")
    held <- union(held, b$gv[!is.na(vals) & vals %in% b$bk[["levels"]]])
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
    if (!any(ok)) ce_plan_refuse(b, nd, r1)
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
  # members sharing a value share its draw and their weights add up
  for (k in mmk) {
    b <- blocks[[k]]
    vals <- vapply(b$gv, function(v) as.character(nd[[v]][r1]), "")
    isnew <- is.na(vals) | !vals %in% b$bk[["levels"]]
    u <- unique(vals[isnew])
    free <- setdiff(b$parts[, 1L], vals[!isnew])
    for (v in b$gv[isnew]) {
      if (is.factor(base[[v]])) free <- intersect(free, levels(base[[v]]))
    }
    if (any(b$gv[isnew] %in% c(held, names(given))) ||
          length(free) < length(u)) {
      ce_plan_refuse(b, nd, r1)
    }
    mt <- match(vals, u)
    for (j in which(isnew)) {
      given[[b$gv[j]]] <- free[mt[j]]
      moved <- c(moved, b$gv[j])
    }
    for (j in seq_along(u)) {
      lv <- match(free[j], b$parts[, 1L])
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
  list(rows = rows, nd = ndp, new = new)
}

