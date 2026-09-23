# Which group-level terms a `re_formula` keeps.
#
# brms reads a one-sided `re_formula` term by term: `update_re_terms()`
# keeps the model's group-level terms that the formula names and drops
# the rest, so `~ (1 | g)` on a fit with `(1 | g) + (1 | h)` keeps `g`
# alone, and `~ (1 | g)` on a fit with `(1 + x | g)` keeps the intercept
# of `g` and drops its slope. frmtmb read the argument as a two-way
# switch, `NA` against everything else, so a formula naming SOME terms
# returned the prediction with ALL of them and said nothing
# (dev/adefects-findings.md section 11, item 7), and a formula naming a
# grouping factor the fit does not have did the same (dev/test-backlog.md).
#
# The resolution here produces a VIEW of the fit: the same estimates,
# with the dropped terms removed from the prediction design. Every
# prediction path already reads the design from the frame, in sample
# through each predictor's `Z` and on new data through the block
# components, so removing a term there removes it from the estimate,
# from the delta-method standard error and from any draw, with one
# construction rather than one per method.

#' The group-level terms of a one-sided formula, as brms reads them.
#'
#' `(lhs | group)`, `(lhs || group)` and brms's `(lhs | id | group)`.
#' Anything outside the bars is ignored, as brms ignores it. A nested
#' group `a/b` expands to `a` and `a:b`, and `||` to one term per
#' column, because that is how the fitted components are stored.
#'
#' @noRd
re_formula_bars <- function(re_formula) {
  rhs <- re_formula[[length(re_formula)]]
  out <- list()
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    # deparse1, not as.character: a namespaced call's head is itself a
    # call, and as.character() would give a length-3 vector here
    op <- deparse1(e[[1L]])
    if (op %in% c("+", "(")) {
      for (a in as.list(e)[-1L]) walk(a)
    } else if (op %in% c("|", "||")) {
      lhs <- e[[2L]]
      # brms's `(1 | p | g)` parses as `(1 | p) | g`: the middle part
      # is the correlation id, which says nothing about which columns
      if (is.call(lhs) && identical(as.character(lhs[[1L]]), "|")) {
        lhs <- lhs[[2L]]
      }
      for (grp in re_expand_group(e[[3L]])) {
        out[[length(out) + 1L]] <<- list(lhs = lhs, group = grp,
                                         double = op == "||",
                                         text = paste0("(", deparse1(e),
                                                       ")"))
      }
    }
    invisible(NULL)
  }
  walk(rhs)
  out
}

#' `a/b` is `a` plus `a:b`; `gr(g, ...)` is `g`. Anything else is kept.
#'
#' @noRd
re_expand_group <- function(g) {
  if (is.call(g) && identical(as.character(g[[1L]]), "gr")) {
    return(list(g[[2L]]))
  }
  if (is.call(g) && identical(as.character(g[[1L]]), "/")) {
    outer <- re_expand_group(g[[2L]])
    inner <- g[[3L]]
    return(c(outer, list(call(":", outer[[length(outer)]], inner))))
  }
  list(g)
}

#' One comparison key per grouping expression: an interaction is the
#' same factor whichever order its parts are written in.
#'
#' @noRd
re_group_key <- function(g) {
  g <- re_expand_group(g)[[1L]]
  parts <- strsplit(deparse1(g), ":", fixed = TRUE)[[1L]]
  paste(sort(trimws(parts)), collapse = ":")
}

#' The design columns one re_formula term asks for, named as the fitted
#' components name theirs.
#'
#' @noRd
re_term_cnms <- function(fit, term) {
  lhs <- term$lhs
  tt <- stats::terms(stats::as.formula(call("~", lhs)))
  if (isTRUE(term$double)) {
    # `(1 + x || g)` is stored as `(1 | g)` and `(0 + x | g)`; the
    # columns are the union, which the per-column components carry
    labs <- attr(tt, "term.labels")
    return(c(if (attr(tt, "intercept") == 1L) "(Intercept)",
             unlist(lapply(labs, function(l) {
               re_term_cnms(fit, list(lhs = call("+", 0, str2lang(l)),
                                      double = FALSE))
             }))))
  }
  if (!length(attr(tt, "term.labels"))) {
    return(if (attr(tt, "intercept") == 1L) "(Intercept)" else character(0))
  }
  df <- fit$frame[["data_frame"]]
  mf <- tryCatch(stats::model.frame(tt, df, na.action = stats::na.pass),
                 error = function(e) NULL)
  if (is.null(mf)) return(NULL)
  colnames(stats::model.matrix(tt, mf))
}

#' The fitted components a re_formula can address, and the ones it
#' cannot name.
#'
#' A component is addressable when it carries a `(lhs | group)` bar.
#' Blocks that `re_formula = NA` drops but that have no bar (a
#' factor-smooth term, a `car()` or `spde()` field) cannot be named in
#' a formula, so a formula that keeps SOME terms cannot say whether they
#' stay; those are returned so the caller can refuse rather than guess.
#'
#' @noRd
re_fit_components <- function(fit) {
  blocks <- fit$frame[["re_blocks"]]
  named <- list()
  unnamed <- character(0)
  for (i in seq_along(blocks)) {
    bk <- blocks[[i]]
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp")) next
    for (j in seq_along(bk[["components"]])) {
      cp <- bk[["components"]][[j]]
      bar <- cp[["bar"]]
      if (is.null(bar) || !is.call(bar)) {
        unnamed <- c(unnamed, cp[["label"]] %||% bk[["covstruct"]])
        next
      }
      named[[length(named) + 1L]] <- list(
        block = i, comp = j, lp_key = cp[["lp_key"]],
        group = re_group_key(bar[[3L]]), cnms = cp[["cnms"]],
        label = deparse1(bar))
    }
  }
  for (lp in fit$frame[["linpreds"]]) {
    for (si in lp[["smooths"]] %||% list()) {
      if (!is.null(si[["group_var"]])) unnamed <- c(unnamed, si[["label"]])
    }
  }
  list(named = named, unnamed = unique(unnamed))
}

#' Resolve `re_formula` against a fit.
#'
#' Returns `list(fit, re_formula)`. `NULL` and `NA` pass through
#' unchanged. A one-sided formula with no group-level term, `~0` or
#' `~1`, is `NA`: brms reads both as "no group-level effects". A
#' formula that keeps every column of every term is `NULL`, and the fit
#' comes back untouched, so its answer is bitwise the full one. Any
#' other formula returns a VIEW of the fit with the dropped terms taken
#' out of the prediction design, and `re_formula = NULL`, so that the
#' callers downstream see an ordinary "keep everything" request on the
#' reduced design.
#'
#' A term that names no group-level term of the fit is REFUSED. brms
#' drops such a term without a word (`check_re_formula()`), which makes
#' a misspelled grouping factor read as "drop the group effects"; a
#' refusal names the term and lists the ones the fit has.
#'
#' @noRd
re_resolve <- function(fit, re_formula, what = "predict()") {
  if (is.null(re_formula) || !inherits(re_formula, "formula")) {
    return(list(fit = fit, re_formula = re_formula))
  }
  terms <- re_formula_bars(re_formula)
  if (!length(terms)) return(list(fit = fit, re_formula = NA))
  fc <- re_fit_components(fit)
  named <- fc$named
  keep <- lapply(named, function(nc) rep(FALSE, length(nc$cnms)))
  have <- vapply(named, function(nc) nc$label %||% "", "")
  for (tm in terms) {
    want <- re_term_cnms(fit, tm)
    gk <- re_group_key(tm$group)
    hit <- FALSE
    if (!is.null(want) && length(want)) {
      for (lk in unique(vapply(named, `[[`, "", "lp_key"))) {
        at <- which(vapply(named, function(nc) {
          nc$lp_key == lk && nc$group == gk
        }, NA))
        if (!length(at)) next
        got <- unlist(lapply(named[at], `[[`, "cnms"))
        # brms's rule (check_re_formula()): the term is found when its
        # columns are a subset of the fitted term's, so `(1 | g)` on a
        # `(1 + x | g)` fit keeps the intercept alone
        if (!all(want %in% got)) next
        for (k in at) {
          keep[[k]] <- keep[[k]] | named[[k]]$cnms %in% want
        }
        hit <- TRUE
      }
    }
    if (!hit) {
      frm_stop(what, ": the re_formula term ", tm$text, " matches no ",
               "group-level term of this fit. A term is kept when its ",
               "grouping factor is one the fit has and its columns are ",
               "among that term's columns. The fit's group-level terms ",
               "are: ", if (length(have)) {
                 paste0("(", unique(have), ")", collapse = ", ")
               } else "none",
               ". brms drops such a term silently, which turns a ",
               "misspelled grouping factor into a population-level ",
               "prediction; use re_formula = NA to ask for that",
               call. = FALSE)
    }
  }
  # every column of every named term kept: the full prediction, and the
  # content a formula cannot name stays in, as it does under NULL and
  # as brms keeps a smooth under any re_formula
  if (all(vapply(keep, all, NA))) {
    return(list(fit = fit, re_formula = NULL))
  }
  if (length(fc$unnamed)) {
    frm_stop(what, ": re_formula = ", deparse1(re_formula), " keeps some ",
             "group-level terms, and this fit also has group-level ",
             "content a formula cannot name: ",
             paste(fc$unnamed, collapse = ", "), ". re_formula = NA drops ",
             "it with every other group-level term and NULL keeps it, ",
             "but a partial formula cannot say which, so it is refused ",
             "rather than guessed", call. = FALSE)
  }
  list(fit = re_view(fit, named, keep), re_formula = NULL)
}

#' A fit whose prediction design holds only the kept group-level
#' columns.
#'
#' In sample, each predictor's `Z` has the dropped columns zeroed, at
#' every level. On new data, a component with nothing kept leaves its
#' block, so its grouping column is neither read nor required, and a
#' partly kept component carries `keep_cols`, which `pred_design()`
#' applies to the rebuilt design. The estimates, the covariance and the
#' cache are the fit's own: nothing that depends on the frame is cached.
#'
#' @noRd
re_view <- function(fit, named, keep) {
  frame <- fit$frame
  blocks <- frame[["re_blocks"]]
  zmask <- list()
  drop_comp <- list()
  for (k in seq_along(named)) {
    nc <- named[[k]]
    kp <- keep[[k]]
    if (all(kp)) next
    bk <- blocks[[nc$block]]
    cp <- bk[["components"]][[nc$comp]]
    D <- bk[["dim"]]
    # from c_idx rather than from levels: every block has the first and
    # a structured one may carry its levels in another shape
    nlev <- length(bk[["c_idx"]]) %/% D
    cols <- cp[["offset"]] + which(!kp)
    idx <- bk[["c_idx"]][as.vector(outer(cols, (seq_len(nlev) - 1L) * D,
                                         `+`))]
    zmask[[nc$lp_key]] <- c(zmask[[nc$lp_key]], idx)
    if (!any(kp)) {
      drop_comp[[length(drop_comp) + 1L]] <- c(nc$block, nc$comp)
    } else {
      blocks[[nc$block]][["components"]][[nc$comp]][["keep_cols"]] <- kp
    }
  }
  # remove from the highest component index down, so the positions the
  # loop above recorded stay valid
  if (length(drop_comp)) {
    dc <- do.call(rbind, drop_comp)
    dc <- dc[order(dc[, 1L], -dc[, 2L]), , drop = FALSE]
    for (r in seq_len(nrow(dc))) {
      blocks[[dc[r, 1L]]][["components"]][[dc[r, 2L]]] <- NULL
    }
  }
  frame[["re_blocks"]] <- blocks
  for (key in names(zmask)) {
    lp <- frame[["linpreds"]][[key]]
    Z <- lp[["Z"]]
    if (is.null(Z)) next
    m <- rep(1, ncol(Z))
    m[unique(zmask[[key]])] <- 0
    frame[["linpreds"]][[key]][["Z"]] <- Z %*% Matrix::Diagonal(x = m)
  }
  fit$frame <- frame
  fit
}
