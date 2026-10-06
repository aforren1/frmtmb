# brms's gp(..., by = ), read from brms 2.23.0 (R/formula-gp.R and
# brms:::data_gp()).
#
# A factor `by` splits the term into one GP per column of the by-factor's
# design matrix: one per level under cmc = TRUE (cell means), and an
# intercept GP plus one contrast GP per remaining level under cmc =
# FALSE. Each sub-GP sees only the rows its column is nonzero on, is
# scaled and centered over those rows alone, and carries its own sd and
# lengthscales; the column's value multiplies it. A numeric `by`
# multiplies a single GP over every row.
#
# frmtmb builds each sub-GP as its own gp or hsgp block, so the
# objective, the Laplace machinery, the sampler and the covariance
# summaries read them with no by branch, exactly as gr(g, by = f) is
# built (R/gr-by.R). What knows about the split is the naming (brms's
# `sdgp_gpxfa`), the multiplier a prediction reads off newdata, and the
# refusal of a by-level the fit never saw, which brms makes too.

#' frmtmb's label for a `gp()` term, which names its blocks.
#'
#' @noRd
gp_term_label <- function(ge, vnames) {
  paste0("gp(", paste(vnames, collapse = ", "),
         if (!is.null(ge$by)) paste0(", by = ", deparse1(ge$by)),
         if (!is.null(ge$k)) paste0(", k = ", ge$k),
         ")")
}

#' brms's `label` of a gp term in `frame_gp()`: `gp` with the renamed
#' covariates and the renamed `by` variable pasted on.
#'
#' @noRd
gp_brms_label <- function(ge, vnames) {
  paste0("gp", brms_rename(paste(vnames, collapse = "")),
         if (!is.null(ge$by)) brms_rename(deparse1(ge$by)))
}

#' The sub-GPs of one `gp()` term: their rows, the multiplier on those
#' rows, a label suffix, brms's `sfx1` (the `sdgp` and `zgp` name) and
#' `sfx2` (the `lscale` names), and what a prediction needs to rebuild
#' the multiplier on new data.
#'
#' brms's `sfx2` is `outer(sfx1, covars, paste0)` for a non-isotropic GP,
#' with the covariates as written and not renamed, and `sfx1` itself
#' when isotropic.
#'
#' @noRd
gp_sub_terms <- function(ge, byv, vnames, iso, n) {
  base <- gp_brms_label(ge, vnames)
  sfx2_of <- function(s1) if (iso) s1 else paste0(s1, vnames)
  if (is.null(ge$by)) {
    return(list(list(rows = seq_len(n), mult = rep(1, n), lab_sfx = "",
                     sfx1 = base, sfx2 = sfx2_of(base), by = NULL)))
  }
  blab <- deparse1(ge$by)
  if (!brms_like_factor(byv)) {
    if (!is.numeric(byv)) {
      frm_stop("gp(by = ", blab, "): a by variable is a factor, which ",
               "splits the term into one GP per level, or numeric, which ",
               "multiplies the GP; this one is ", class(byv)[1L],
               call. = FALSE)
    }
    return(list(list(rows = seq_len(n), mult = as.numeric(byv),
                     lab_sfx = "", sfx1 = base, sfx2 = sfx2_of(base),
                     by = list(type = "numeric", expr = ge$by,
                               label = blab))))
  }
  # brms drops unused levels from the data before it builds the term
  # (brm(drop_unused_levels = TRUE)), so a level with no row has no GP
  fac <- droplevels(as.factor(byv))
  lev <- levels(fac)
  C <- gp_by_contrasts(fac, isTRUE(ge$cmc), attr(byv, "contrasts"))
  cons <- by_rm_wsp(sub("^byval", "", brms_rename(colnames(C))))
  out <- list()
  for (j in seq_len(ncol(C))) {
    mult <- unname(C[as.integer(fac), j])
    rows <- which(mult != 0)
    if (!length(rows)) {
      # a contrast column that is zero on every row: only a contrast
      # matrix with a zero column gets here, the levels being in use
      frm_stop("gp(by = ", blab, "): the design column '", colnames(C)[j],
               "' of the by variable is zero on every row, so its GP has ",
               "nothing to be fitted to", call. = FALSE)
    }
    s1 <- paste0(base, cons[j])
    out[[j]] <- list(
      rows = rows, mult = mult[rows],
      lab_sfx = paste0(" [", blab, " = ",
                       if (isTRUE(ge$cmc)) lev[j] else cons[j], "]"),
      sfx1 = s1, sfx2 = sfx2_of(s1),
      by = list(type = "factor", expr = ge$by, label = blab, levels = lev,
                contr = unname(C), j = j))
  }
  out
}

#' The by-factor's design over its levels, one row per level: brms's
#' `model.matrix(~ 0 + byval)` under `cmc = TRUE` and `~ 1 + byval`
#' otherwise, with the factor's own contrasts, as `as.factor()` keeps
#' them in brms.
#'
#' @noRd
gp_by_contrasts <- function(fac, cmc, ctr = NULL) {
  lev <- levels(fac)
  u <- factor(lev, levels = lev, ordered = is.ordered(fac))
  # droplevels() drops the attribute, so the caller passes it
  if (is.matrix(ctr) && nrow(ctr) == length(lev)) stats::contrasts(u) <- ctr
  form <- if (cmc) ~ 0 + byval else ~ 1 + byval
  stats::model.matrix(form, data.frame(byval = u))
}

#' The multiplier one sub-GP applies to the rows of `newdata`: 1 when
#' the term has no `by`, the by value when it is numeric, and the
#' by-factor's design column otherwise. A by-level the fit never saw is
#' refused with brms's own message: brms's `validate_newdata()` stops on
#' it before any prediction is made, `allow_new_levels = TRUE` included
#' (measured on brms 2.23.0, `dev/gpby-brms-newlevel.R`).
#'
#' @noRd
gp_by_mult <- function(by, newdata, env, n) {
  if (is.null(by)) return(rep(1, n))
  v <- tryCatch(eval(by$expr, newdata, env), error = function(e) NULL)
  if (is.null(v) || length(v) != n) {
    frm_stop("newdata must hold the by variable '", by$label, "' of a ",
             "gp(by = ) term, as in brms, whose prediction reads it to ",
             "choose the GP each row adds", call. = FALSE)
  }
  if (anyNA(v)) {
    frm_stop("gp(by = ", by$label, "): the by variable has missing ",
             "values in newdata", call. = FALSE)
  }
  if (identical(by$type, "numeric")) {
    if (!is.numeric(v)) {
      frm_stop("gp(by = ", by$label, "): the fit read a numeric by ",
               "variable, and newdata gives ", class(v)[1L], call. = FALSE)
    }
    return(as.numeric(v))
  }
  k <- match(as.character(v), by$levels)
  if (anyNA(k)) {
    found <- levels(as.factor(v))
    frm_stop("New factor levels are not allowed.\nLevels allowed: ",
             paste0("'", by$levels, "'", collapse = ", "),
             "\nLevels found: ", paste0("'", found, "'", collapse = ", "),
             call. = FALSE)
  }
  by$contr[k, by$j]
}

#' The `gp()` terms of a fit as brms names them: one entry per written
#' term, in block order, holding its sub-GP blocks in level order, the
#' brms prefix of its linear predictor, and brms's `sdgp_` and `lscale_`
#' names. `lscale` is a `levels x dimensions` matrix whose column-major
#' order is brms's draw order (`as.vector(sfx2)` in `rename_gp()`).
#'
#' @noRd
gp_brms_terms <- function(fit) {
  blocks <- fit$frame[["re_blocks"]] %||% list()
  out <- list()
  for (bi in seq_along(blocks)) {
    bk <- blocks[[bi]]
    if (!bk[["covstruct"]] %in% c("gp", "hsgp")) next
    meta <- bk[["gp_brms"]]
    if (is.null(meta)) next
    key <- meta$term
    if (is.null(out[[key]])) {
      cp <- bk[["components"]][[1L]]
      lp <- fit$frame[["linpreds"]][[cp[["lp_key"]]]]
      pre <- if (is.null(lp)) "" else brms_lp_prefix(fit, lp)
      out[[key]] <- list(blocks = integer(0), pre = pre,
                         sdgp = character(0), lscale = list())
    }
    e <- out[[key]]
    e$blocks <- c(e$blocks, bi)
    e$sdgp <- c(e$sdgp, paste0("sdgp_", brms_usc(e$pre, meta$sfx1)))
    e$lscale[[length(e$lscale) + 1L]] <-
      paste0("lscale_", brms_usc(e$pre, meta$sfx2))
    out[[key]] <- e
  }
  lapply(unname(out), function(e) {
    e$lscale <- do.call(rbind, e$lscale)
    e
  })
}

#' brms's `sdgp_` and `lscale_` values of every `gp()` term at one
#' `theta`, in brms's order: a term's standard deviations, level by
#' level, then its length scales, levels innermost. The length scale is
#' on brms's scaled inputs, `gp(scale = TRUE)`: the exact form holds it
#' in data units and divides by the largest distance between two
#' positions of its own rows.
#'
#' @noRd
gp_brms_values <- function(fit, th) {
  blocks <- fit$frame[["re_blocks"]]
  out <- numeric(0)
  for (e in gp_brms_terms(fit)) {
    sd <- numeric(0)
    ls <- matrix(NA_real_, length(e$blocks), ncol(e$lscale))
    for (k in seq_along(e$blocks)) {
      bk <- blocks[[e$blocks[k]]]
      t0 <- th[bk[["theta_idx"]]]
      sd <- c(sd, exp(t0[1L]))
      ls[k, ] <- exp(t0[-1L]) / (bk[["gp_lscale_div"]] %||% 1)
    }
    out <- c(out, stats::setNames(sd, e$sdgp),
             stats::setNames(as.vector(ls), as.vector(e$lscale)))
  }
  out
}

#' The `theta` segment's labels with every `gp()` hyperparameter under
#' brms's name, `sdgp_<term>` and `lscale_<term>`, the rest unchanged.
#'
#' @noRd
gp_theta_labels <- function(fit, v) {
  blocks <- fit$frame[["re_blocks"]]
  for (e in gp_brms_terms(fit)) {
    for (k in seq_along(e$blocks)) {
      ti <- blocks[[e$blocks[k]]][["theta_idx"]]
      v[ti[1L]] <- e$sdgp[k]
      v[ti[-1L]] <- e$lscale[k, ]
    }
  }
  v
}

#' The map between a `gp()` block's internal `theta` and brms's
#' `sdgp_` and `lscale_` draws, for a sampler that stores draws under
#' brms's names and on brms's scales: `theta_pos` is the position in
#' `theta`, `linkinv` takes the internal value to brms's and `linkfun`
#' takes it back. The sd is `exp(theta)`; the length scale is
#' `exp(theta) / div`, with `div` the distance the block's inputs are
#' divided by on brms's side and not on this one.
#'
#' @noRd
gp_brms_natural <- function(fit) {
  blocks <- fit$frame[["re_blocks"]]
  out <- list(names = character(0), theta_pos = integer(0),
              linkinv = list(), linkfun = list())
  add <- function(nm, pos, div) {
    out$names <<- c(out$names, nm)
    out$theta_pos <<- c(out$theta_pos, pos)
    out$linkinv[[length(out$linkinv) + 1L]] <<- local({
      d <- div
      function(x) exp(x) / d
    })
    out$linkfun[[length(out$linkfun) + 1L]] <<- local({
      d <- div
      function(x) log(x * d)
    })
  }
  for (e in gp_brms_terms(fit)) {
    for (k in seq_along(e$blocks)) {
      bk <- blocks[[e$blocks[k]]]
      ti <- bk[["theta_idx"]]
      add(e$sdgp[k], ti[1L], 1)
      for (d in seq_along(ti[-1L])) {
        add(e$lscale[k, d], ti[1L + d], bk[["gp_lscale_div"]] %||% 1)
      }
    }
  }
  out
}
