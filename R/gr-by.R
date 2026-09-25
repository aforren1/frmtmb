# brms's gr(g, by = f): one set of group-level standard deviations and
# correlations per level of `f`, where every level of `g` belongs to
# exactly one level of `f`.
#
# The covariance of the term's coefficients is therefore block-diagonal
# across the levels of `f`: the levels of `g` inside one by-level share
# that by-level's covariance, and levels in different by-levels are
# independent. brms writes it with one `sd_1` column and one Cholesky
# factor `L_1[k]` per by-level, indexed per level of `g` through `Jby_1`
# (`scale_r_cor_by()`), which is the same density.
#
# frmtmb builds it as that block-diagonal: the term becomes one ordinary
# random-effect block per by-level, over the levels of `g` in that
# by-level. Each is the covariance structure the term asks for, so a
# by-split `ar1()` or `cs()` term is one such block per by-level, and the
# Laplace machinery, the priors, `VarCorr()`, `confint()`, the sampler
# and the simulator read them with no by branch. The precision stays
# block-sparse: splitting adds no coupling between blocks. What knows
# about the split is the naming (brms's `sd_g__Intercept:fa`), the views
# that brms keys by the grouping factor (`ranef()`, `coef()`,
# `ngrps()`), and prediction on new data, where an unseen level of `g`
# takes the covariance of its own by-level.

#' The structures a by-split term can carry: every structure whose
#' density factorizes over the levels of the grouping factor, so that
#' one block per by-level is the whole model.
#'
#' @noRd
gr_by_structures <- c("us", "diag", "homdiag", "cs", "homcs", "ar1",
                      "hetar1", "toep", "homtoep", "ou", "exp", "gau",
                      "mat", "us_t", "diag_t")

#' Parse `gr(g, by = f)`.
#'
#' `by` is an expression evaluated against the data, as brms evaluates
#' it (`eval2(by, data)`); `label` is its deparsed form, which is what
#' brms pastes in front of each by-level in a parameter name.
#'
#' @noRd
parse_gr_by <- function(by, gvar, cls, bar, has_known) {
  if (has_known) {
    frm_stop("gr(g, by = , cov = ) / gr(g, by = , prec = ) is not ",
             "supported: ", deparse1(bar), ". brms scales each level's ",
             "effects by the covariance of its by-level and only then ",
             "correlates the levels through the relationship matrix, so ",
             "levels in different by-levels stay correlated and the ",
             "covariance is neither block-diagonal over the by-levels nor ",
             "a Kronecker product. Use gr(g, cov = A) without by, or ",
             "gr(g, by = f) without the matrix", call. = FALSE)
  }
  if (calls_function(gvar, "mm")) {
    frm_stop("gr(mm(...), by = ) is not supported: ", deparse1(bar),
             ". Write the multi-membership term as brms writes it, ",
             "mm(g1, g2, by = cbind(f1, f2))", call. = FALSE)
  }
  if (!cls %in% gr_by_structures) {
    frm_stop("gr(by = ) is not supported for ", cls, "(): ",
             deparse1(bar), ". A by-split term is one block per by-level ",
             "with its own parameters, which covers the structures whose ",
             "density factorizes over the grouping levels (",
             paste(setdiff(gr_by_structures, c("us_t", "diag_t")),
                   collapse = ", "),
             ")",
             switch(cls,
                    rr = paste0("; rr() shares its loadings across the ",
                                "whole factor, so write one rr() term per ",
                                "subset of the data instead"),
                    equalto = paste0("; equalto() has no parameters to ",
                                     "split, so drop by"),
                    ""),
             call. = FALSE)
  }
  list(expr = by, label = deparse1(by), fn = "gr")
}

#' The data columns a random-effect term's by-variable reads, from
#' `gr(g, by = )` or `mm(g1, g2, by = )`.
#'
#' @noRd
by_term_vars <- function(rt) {
  unique(c(all.vars(rt$by$expr), all.vars(rt$mm$by$expr)))
}

#' brms's `rm_wsp()`: level labels lose their whitespace in names.
#'
#' @noRd
by_rm_wsp <- function(x) gsub("[ \t\r\n]+", "", x)

#' brms's `extract_levels()`: a factor's own levels, unused ones
#' included, and otherwise the sorted distinct values.
#'
#' @noRd
by_extract_levels <- function(x) {
  if (!is.factor(x)) x <- factor(x)
  levels(x)
}

#' How a refusal spells the by-argument of a term.
#'
#' @noRd
by_arg <- function(by) paste0(by$fn %||% "gr", "(by = ", by$label, ")")

#' The value of a by-variable on some rows, with a refusal that names the
#' term when it cannot be evaluated. A multi-membership term's is a
#' matrix with one column per member (`members`), as brms requires.
#'
#' @noRd
by_eval <- function(by, data, env, n, what, members = NULL) {
  v <- tryCatch(eval(by$expr, data, env), error = function(e) {
    frm_stop(by_arg(by), ": cannot evaluate the by-variable ", what, ": ",
             conditionMessage(e), call. = FALSE)
  })
  if (!is.null(members)) {
    v <- as.matrix(v)
    if (!identical(dim(v), c(as.integer(n), as.integer(members)))) {
      # brms's frame_re() refusal, verbatim
      frm_stop("Grouping structure 'mm' expects 'by' to be a matrix with as ",
               "many columns as grouping factors. ", by_arg(by), " gives ",
               nrow(v), " x ", ncol(v), " ", what, ", and ", n, " x ",
               members, " is needed; build it with cbind(f1, f2)",
               call. = FALSE)
    }
    return(v)
  }
  if (length(v) == 1L && n != 1L) v <- rep(v, n)
  if (length(v) != n || !is.null(dim(v))) {
    frm_stop(by_arg(by), " must give one value per row ", what, "; got ",
             length(v), " for ", n, " rows", call. = FALSE)
  }
  v
}

#' The by-level of each grouping level, from the row-level pairs,
#' refused with brms's message when a level carries two.
#'
#' `J` indexes the grouping levels row by row (a matrix with one column
#' per member for a multi-membership term, whose rows then give one
#' by-value per member), `byv` is the by-variable in the same shape.
#'
#' @noRd
by_level_map <- function(J, byv, n_levels, groups, by) {
  byc <- as.character(byv)
  ok <- !is.na(J)
  if (anyNA(byc[ok])) {
    frm_stop(by_arg(by), " has missing values on rows that carry a level ",
             "of ", paste(groups, collapse = ", "),
             "; each level needs its by-level", call. = FALSE)
  }
  pairs <- unique(data.frame(J = as.vector(J)[ok], by = byc[ok],
                             stringsAsFactors = FALSE))
  if (anyDuplicated(pairs$J)) {
    # brms's frame_re() refusal, verbatim, with the reason after it
    frm_stop("Some levels of ", paste0("'", groups, "'", collapse = ", "),
             " correspond to multiple levels of '", by$label, "'. A by-split ",
             "term gives each level of the grouping factor the covariance ",
             "of ONE by-level, so the by-variable must be constant within ",
             "each grouping level", call. = FALSE)
  }
  out <- rep(NA_character_, n_levels)
  out[pairs$J] <- pairs$by
  out
}

#' Split one assembled random-effect component into one component per
#' by-level.
#'
#' `lev_by` is the by-level of each of the component's levels
#' (`by_level_map()`), `all_by` the by-variable's levels as brms reads
#' them. A by-level with no grouping level in the data has no data for
#' its covariance, so it gets no block, and a message says so; brms
#' keeps it, informed by its prior alone.
#'
#' Each sub-component keeps the grouping factor, the design columns and
#' the covariance structure of the term, so everything that reads a
#' block reads it as it reads any other. It carries `by`, which the
#' naming, the grouped views and prediction on new data read. `key`
#' ties the sub-components of one term together, and `group_levels`
#' with `level_by` rebuild the term's full level set.
#'
#' @noRd
by_split_component <- function(cp, by, lev_by, all_by) {
  used <- all_by[all_by %in% lev_by]
  unused <- setdiff(all_by, used)
  if (length(unused)) {
    frm_message(by_arg(by), " in `", cp$label, "`: no level of ",
                cp$group_name, " falls in by-level ",
                paste0("'", unused, "'", collapse = ", "),
                ", so there is no data for its covariance and it is left ",
                "out of the model")
  }
  d <- cp$dim
  info <- list(expr = by$expr, label = by$label, fn = by$fn, levels = used,
               group_levels = cp$levels, level_by = lev_by, key = cp$label)
  lapply(used, function(b) {
    sel <- which(lev_by == b)
    cols <- as.vector(outer(seq_len(d), (sel - 1L) * d, "+"))
    sub <- cp
    sub$levels <- cp$levels[sel]
    sub$n_levels <- length(sel)
    sub$Zlocal <- cp$Zlocal[, cols, drop = FALSE]
    sub$label <- paste0(cp$label, " [", by$label, " = ", b, "]")
    if (!is.null(cp$id)) sub$id <- paste0(cp$id, "\rby = ", b)
    sub$by <- c(info, list(level = b,
                           name = paste0(by$label, by_rm_wsp(b))))
    sub
  })
}

#' Refuse two different by-variables on one grouping factor, with
#' brms's message: brms attaches the by-levels to the grouping factor's
#' levels, so a factor has one by-variable across the whole model.
#'
#' @noRd
check_gr_by_groups <- function(spec) {
  seen <- list()
  for (resp in spec$responses) {
    for (dp in resp$dpars) {
      for (z in dp[["re"]] %||% list()) {
        by <- z$by %||% z$mm$by
        if (is.null(by)) next
        g <- if (is.null(z$mm)) deparse1(z$bar[[3L]]) else z$mm$label
        seen[[g]] <- unique(c(seen[[g]], by$label))
        if (length(seen[[g]]) > 1L) {
          frm_stop("Each grouping factor can only be associated with one ",
                   "'by' variable. '", g, "' is split by ",
                   paste0("'", seen[[g]], "'", collapse = " and "),
                   "; give the second term a copy of the grouping column ",
                   "under another name", call. = FALSE)
        }
      }
    }
  }
  invisible(NULL)
}

#' The blocks as brms keys them by grouping factor: the sub-blocks of
#' one by-split term merged back into one block over the term's full
#' level set, in the grouping factor's level order. `c_idx` and `b_idx`
#' are gathered level by level, so a reader that reshapes them by `dim`
#' sees the term as if it had never been split. Every other block is
#' returned as it is.
#'
#' @noRd
by_merged_blocks <- function(blocks) {
  keys <- vapply(blocks, function(bk) bk[["by"]][["key"]] %||% "", "")
  out <- list()
  done <- character(0)
  for (i in seq_along(blocks)) {
    k <- keys[i]
    if (!nzchar(k)) {
      out[[length(out) + 1L]] <- blocks[[i]]
      next
    }
    if (k %in% done) next
    done <- c(done, k)
    subs <- blocks[keys == k]
    by <- subs[[1L]][["by"]]
    D <- subs[[1L]][["dim"]]
    lev <- by$group_levels
    c_idx <- integer(0)
    b_idx <- integer(0)
    for (l in seq_along(lev)) {
      s <- match(by$level_by[l], by$levels)
      bk <- subs[[s]]
      p <- match(lev[l], bk[["levels"]])
      pos <- (p - 1L) * D + seq_len(D)
      c_idx <- c(c_idx, bk[["c_idx"]][pos])
      b_idx <- c(b_idx, bk[["b_idx"]][pos])
    }
    m <- subs[[1L]]
    m$levels <- lev
    m$n_levels <- length(lev)
    m$c_idx <- c_idx
    m$b_idx <- b_idx
    m$term_label <- gsub(paste0(" [", by$label, " = ", by$level, "]"), "",
                         m$term_label, fixed = TRUE)
    m$by <- NULL
    m$by_split <- by
    out[[length(out) + 1L]] <- m
  }
  out
}

#' Which rows of new data a by-split sub-block serves.
#'
#' A level of the grouping factor the fit saw keeps its own fitted
#' effect whatever the new row says about the by-variable, as in brms,
#' where the effect is indexed by the level alone; it is served by the
#' sub-block that holds it and by no other. A level the fit did not see
#' is drawn from the covariance of the by-level the new row names, as
#' brms's `get_new_rdraws()` draws it, so it is new in that sub-block
#' and absent from the others.
#'
#' `gv` is the grouping level of each row as character and `j` its match
#' in this sub-block's levels. For a multi-membership term they are the
#' per-member list and the rows x members matrix, and the by-variable is
#' read one column per member. Returns `off`, the entries this sub-block
#' does not serve (their design row is zeroed), and `new`, the entries
#' that are an unseen level of it, both shaped like `j`.
#'
#' @noRd
by_route_rows <- function(bk, gv, j, newdata, env) {
  by <- bk[["by"]]
  members <- if (is.list(gv)) length(gv) else NULL
  gm <- if (is.null(members)) gv else do.call(cbind, lapply(gv, as.character))
  gm <- matrix(gm, NROW(j), NCOL(j))
  unseen <- is.na(j) & !(gm %in% by$group_levels & !is.na(gm))
  off <- is.na(j) & !unseen
  if (any(unseen)) {
    byc <- as.character(by_eval(by, newdata, env, NROW(j), "on newdata",
                                members))
    miss <- unseen & is.na(byc)
    if (any(miss)) {
      frm_stop("A level of ", bk[["group_name"]], " the fit did not see ",
               "takes the covariance of its by-level, and ", by_arg(by),
               " is missing on ", sum(miss), " such entries of newdata. ",
               "Supply it, or predict at the population level with ",
               "re_formula = NA", call. = FALSE)
    }
    bad <- unseen & !byc %in% by$levels
    if (any(bad)) {
      frm_stop("New level(s) of ", bk[["group_name"]], " fall in ",
               by$label, " = ", paste0("'", unique(byc[bad]), "'",
                                       collapse = ", "),
               ", which has no fitted covariance: every by-level the fit ",
               "estimated is ", paste0("'", by$levels, "'", collapse = ", "),
               call. = FALSE)
    }
    nv <- which(!is.na(gm) & unseen)
    un <- unique(gm[nv])
    by_level_map(match(gm[nv], un), byc[nv], length(un),
                 bk[["group_name"]], by)
    off <- off | (unseen & byc != by$level)
  }
  new <- is.na(j) & !off
  if (is.null(members)) {
    off <- as.vector(off)
    new <- as.vector(new)
  }
  list(off = off, new = new)
}
