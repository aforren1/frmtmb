# brms's parameter names, built from a fit's structure.
#
# brms names a population-level coefficient `b_<dpar>_<resp>_<coef>`,
# a group-level standard deviation `sd_<group>__<dpar>_<resp>_<coef>`
# and a group-level coefficient `r_<group>__<dpar>_<resp>[<level>,<coef>]`,
# dropping every prefix part that is `mu` or that a univariate model
# does not have (`brms:::combine_prefix()`, whose part order is
# `vars_prefix()`: dpar, resp, nlpar). Every piece passes through the
# renamers brms itself applies, ported below and read from brms 2.23.0
# in dev/brmsnames-log/brms-bodies6.txt to brms-bodies9.txt:
#
# - a design column through `rename()` (`get_model_matrix()`), so
#   `I(x^2)` is `IxE2` and a level `c-d` is `cMd`;
# - a response through `make_stan_names()`, so `y_a` is `ya`;
# - a group through `rename()`, which keeps the `:` of `g:h2` and turns
#   `mm(g1, g2)` into `mmg1g2`;
# - a level of an interaction group through `combine_groups()`, which
#   joins the parts with `_`, and a level in an `r_` name through
#   `rename_re_levels()`, which turns whitespace into `.`;
# - a smooth's unpenalized columns through `rename_sm()`,
#   `bs_<label>_<k>`.
#
# frmtmb's internal names order the prefix the other way round and keep
# the raw column names, so the brms spelling cannot be had by pasting a
# prefix on the internal one. It is built here, in one place, so that
# `variables()`, `hypothesis()`, `VarCorr()`, `frm_simulate()` and the
# names frmtmb.sample puts on draws cannot drift apart.

#' brms's `rename()`: its default character map, applied literally, and
#' its duplicate refusal.
#'
#' @noRd
brms_rename <- function(x, pattern = NULL, replacement = NULL,
                        fixed = TRUE, check_dup = FALSE) {
  pattern <- as.character(pattern)
  replacement <- as.character(replacement)
  if (!length(pattern) && !length(replacement)) {
    pattern <- c(" ", "(", ")", "[", "]", ",", "\"", "'", "?", "+", "-",
                 "*", "/", "^", "=", "$")
    replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ", "USD")
  }
  if (length(replacement) == 1L) {
    replacement <- rep(replacement, length(pattern))
  }
  stopifnot(length(pattern) == length(replacement))
  keep <- nzchar(pattern)
  pattern <- pattern[keep]
  replacement <- replacement[keep]
  out <- x
  for (i in seq_along(pattern)) {
    out <- gsub(pattern[i], replacement[i], out, fixed = fixed)
  }
  dup <- duplicated(out)
  if (check_dup && any(dup)) {
    dup <- x[out %in% out[dup]]
    frm_stop("Internal renaming led to duplicated names. Consider renaming ",
             "your variables to have different suffixes.\nOccured for: ",
             paste0("'", dup, "'", collapse = ", "), call. = FALSE)
  }
  out
}

#' brms's `make_stan_names()`, which is how a response name enters every
#' parameter name.
#'
#' @noRd
brms_stan_name <- function(x) {
  gsub("[._]", "", make.names(x, unique = TRUE))
}

#' The brms prefix of one linear predictor: dpar, then response, with
#' `mu` and a univariate model's response left out. A nonlinear
#' parameter is a dpar here.
#'
#' @noRd
brms_lp_prefix <- function(fit, lp) {
  mv <- length(fit$spec$responses) > 1L
  parts <- c(if (!identical(lp[["dpar"]], "mu")) lp[["dpar"]],
             if (mv) brms_stan_name(lp[["resp"]]))
  paste(parts, collapse = "_")
}

#' `prefix_name`, or `name` when the prefix is empty.
#'
#' @noRd
brms_usc <- function(prefix, name) {
  if (nzchar(prefix)) paste0(prefix, "_", name) else name
}

#' Whether the user wrote a formula for this distributional parameter.
#' brms reports a dpar nobody modeled as the parameter itself, `sigma`
#' on its natural scale, and one written as `sigma ~ 1` as the
#' coefficient `b_sigma_Intercept`.
#'
#' @noRd
brms_dpar_written <- function(fit, lp) {
  if (identical(lp[["par"]], "beta")) return(TRUE)
  bform <- fit$bform
  if (is.null(bform)) return(TRUE)
  forms <- if (inherits(bform, "frmtmb_mvformula")) bform$forms else
    list(bform)
  i <- match(lp[["resp"]], names(fit$spec$responses))
  if (is.na(i) || i > length(forms)) return(TRUE)
  f <- forms[[i]]
  lp[["dpar"]] %in% c(names(f$pforms), names(f$nlforms))
}

#' The family of the response a linear predictor belongs to.
#'
#' @noRd
brms_lp_family <- function(fit, lp) {
  rs <- fit$spec$responses
  r <- if (is.null(lp[["resp"]])) rs[[1L]] else rs[[lp[["resp"]]]]
  r[["family"]]
}

#' The mixing weights of a finite mixture, as brms reports them: the
#' simplex `theta1 ... thetaK`, where frmtmb estimates the `K - 1` log
#' ratios against the last component. The map is not elementwise, so it
#' is a pair of closures over the whole row: `to_simplex()` takes an
#' `n x (K - 1)` matrix of log ratios to the `n x K` simplex, and
#' `to_logits()` takes it back.
#'
#' @noRd
brms_simplex_maps <- function(K) {
  list(
    to_simplex = function(E) {
      E <- cbind(matrix(E, ncol = K - 1L), 0)
      top <- apply(E, 1L, max)
      W <- exp(E - top)
      W / rowSums(W)
    },
    to_logits = function(P) {
      P <- matrix(P, ncol = K)
      log(P[, seq_len(K - 1L), drop = FALSE]) - log(P[, K])
    }
  )
}

#' One row per estimated coefficient, in `estimated_coef_names()` order:
#' its internal name, brms's name, and whether brms reports it as the
#' distributional parameter itself on the natural scale, in which case
#' the `linkinv` attribute maps the internal value onto that scale and
#' `linkfun` back. `dpar` and `resp` name the predictor of each row.
#'
#' An intercept-only distributional parameter nobody wrote a formula for
#' is `sigma`, `shape`, `nu` or, on a multivariate model, `sigma_ya`, as
#' in brms. Any other coefficient is `b_<prefix>_<column>`, a
#' smooth's unpenalized column is `bs_<prefix>_<label>_<k>`, and a
#' monotonic term's scale is `bsp_<prefix>_mo<x>`. A name brms
#' would give twice across predictors takes brms's `__1` suffix
#' (`repair_stanfit()`); within one predictor the model is refused when
#' it is assembled, as brms refuses it.
#'
#' Two kinds of distributional parameter are not elementwise on their
#' natural scale, and a `linkinv` of the one coefficient would name a
#' value that is not the parameter:
#'
#' - a mixture's mixing weights. The `K - 1` log ratios are natural
#'   only together, as brms's simplex `theta1 ... thetaK`; the
#'   `simplex` attribute lists each such group with its rows, its `K`
#'   names and the maps of `brms_simplex_maps()`, and those rows carry
#'   no `linkinv`. A mixture with any theta formula keeps every theta as
#'   a coefficient.
#' - the dpars a family lists in `link_scale_dpars`, such as an `hmm()`
#'   transition logit, which one row of a softmax reads. brms has no
#'   such parameter, so they stay coefficients.
#'
#' @noRd
brms_coef_table <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  raw <- list(beta = names(tpl[["beta"]]), betad = names(tpl[["betad"]]))
  nm <- lapply(raw, function(v) {
    if (length(v)) paste0("b_", brms_rename(v)) else character(0)
  })
  nat <- lapply(raw, function(v) logical(length(v)))
  dpar <- lapply(raw, function(v) rep(NA_character_, length(v)))
  resp <- dpar
  inv <- lapply(raw, function(v) vector("list", length(v)))
  fun <- inv
  # mixing weights that qualify, per response: betad index and prefix
  thetas <- list()
  for (lp in fit$frame[["linpreds"]]) {
    cn <- colnames(lp[["X"]])
    par <- lp[["par"]]
    if (!length(cn) || !par %in% c("beta", "betad")) next
    pre <- brms_lp_prefix(fit, lp)
    dpar[[par]][lp[["idx"]]] <- lp[["dpar"]]
    resp[[par]][lp[["idx"]]] <- lp[["resp"]] %||% NA_character_
    fam <- brms_lp_family(fit, lp)
    scalar <- identical(par, "betad") && identical(cn, "(Intercept)") &&
      is.null(lp[["Z"]]) && is.null(lp[["constant"]]) &&
      is.null(lp[["nl_body"]]) && !brms_dpar_written(fit, lp) &&
      !lp[["dpar"]] %in% fam[["link_scale_dpars"]]
    if (scalar) {
      nm[[par]][lp[["idx"]]] <- pre
      nat[[par]][lp[["idx"]]] <- TRUE
      mix <- fam[["post"]][["dpar_response"]][["dpars"]]
      if (lp[["dpar"]] %in% mix && grepl("^theta[0-9]+$", lp[["dpar"]])) {
        key <- lp[["resp"]] %||% ".univariate"
        thetas[[key]][[lp[["dpar"]]]] <- list(idx = lp[["idx"]], pre = pre,
                                              lp = lp, mix = mix)
        next
      }
      inv[[par]][[lp[["idx"]]]] <- lp[["link"]]$linkinv
      fun[[par]][[lp[["idx"]]]] <- lp[["link"]]$linkfun
      next
    }
    out <- paste0("b_", brms_usc(pre, brms_rename(cn)))
    mo <- vapply(lp[["mo"]] %||% list(), function(e) as.integer(e[["col"]]),
                 1L)
    if (length(mo)) {
      out[mo] <- paste0("bsp_", brms_usc(pre, brms_rename(cn[mo])))
    }
    fx <- grepl("[.]fx[0-9]+$", cn)
    if (any(fx)) {
      lab <- brms_rename(sub("[.]fx[0-9]+$", "", cn[fx]))
      k <- sub("^.*[.]fx([0-9]+)$", "\\1", cn[fx])
      out[fx] <- paste0("bs", if (nzchar(pre)) paste0("_", pre), "_", lab,
                        "_", k)
    }
    nm[[par]][lp[["idx"]]] <- out
  }
  keep_d <- setdiff(seq_along(raw$betad), fit$frame[["betad_fixed_idx"]])
  n_beta <- length(raw$beta)
  groups <- list()
  for (th in thetas) {
    mix <- th[[1L]]$mix
    if (!all(mix %in% names(th))) {
      # a theta with a formula: the log ratios are coefficients, and so
      # are the ones without, since no simplex can be formed from them
      for (e in th) {
        nm$betad[e$idx] <- paste0("b_", brms_usc(e$pre, "Intercept"))
        nat$betad[e$idx] <- FALSE
      }
      next
    }
    K <- length(mix) + 1L
    lp1 <- th[[mix[1L]]]$lp
    names_k <- vapply(seq_len(K), function(k) {
      lp1[["dpar"]] <- paste0("theta", k)
      brms_lp_prefix(fit, lp1)
    }, "")
    groups[[length(groups) + 1L]] <- c(
      list(pos = n_beta + match(vapply(th[mix], function(e) {
        as.integer(e$idx)
      }, 1L), keep_d),
           names = names_k),
      brms_simplex_maps(K))
  }
  tab <- data.frame(
    internal = c(raw$beta, raw$betad[keep_d]),
    brms = make.unique(c(nm$beta, nm$betad[keep_d]), sep = "__"),
    natural = c(nat$beta, nat$betad[keep_d]),
    dpar = c(dpar$beta, dpar$betad[keep_d]),
    resp = c(resp$beta, resp$betad[keep_d]),
    stringsAsFactors = FALSE
  )
  attr(tab, "linkinv") <- c(inv$beta, inv$betad[keep_d])
  attr(tab, "linkfun") <- c(fun$beta, fun$betad[keep_d])
  attr(tab, "simplex") <- groups
  # where a linear predictor's own coefficients sit in this table: a
  # `beta` index is its own row, and a `betad` index has to skip the
  # rows a fixed betad dropped. coef() needs the map to name its
  # columns as brms does, and rebuilding it at the call site is what
  # would drift.
  attr(tab, "n_beta") <- n_beta
  attr(tab, "keep_d") <- keep_d
  tab
}

#' The table rows of one linear predictor's design columns, in design
#' column order, `NA` where the table has no row for it.
#'
#' @noRd
brms_lp_rows <- function(lp, tab) {
  idx <- lp[["idx"]]
  if (identical(lp[["par"]], "beta")) return(as.integer(idx))
  if (identical(lp[["par"]], "betad")) {
    return(attr(tab, "n_beta") + match(idx, attr(tab, "keep_d")))
  }
  rep(NA_integer_, length(idx))
}

#' brms's names for one linear predictor's design columns, the class
#' prefix (`b_`, `bs_`, `bsp_`, `bcs_`) dropped so that they read as
#' `fixef()` and `vcov()` name their rows.
#'
#' The design column name is kept where the table has no row for it,
#' which is a coefficient held fixed, and where the row is a NATURAL
#' one: brms calls that a spec_par and reports it on its own scale, so
#' putting `sigma` on a coefficient that is log sigma would name a
#' value that is not the one under the name. Only the rows `fixef()`
#' and `vcov()` actually carry are renamed, which is the pairing that
#' has to agree.
#'
#' @noRd
brms_lp_coef_names <- function(lp, tab) {
  cn <- colnames(lp[["X"]])
  r <- brms_lp_rows(lp, tab)
  out <- cn
  ok <- !is.na(r) & !tab$natural[r]
  out[ok] <- sub("^(b|bs|bsp|bcs)_", "", tab$brms[r[ok]])
  out
}

#' brms names of the estimated coefficients, parallel to
#' `estimated_coef_names(fit)`; see `brms_coef_table()`.
#'
#' @noRd
brms_coef_names <- function(fit) brms_coef_table(fit)$brms

#' brms's group name for a block: the grouping term with brms's
#' `rename()` applied, which keeps an interaction's `:`.
#'
#' @noRd
brms_group_name <- function(bk) {
  g <- bk[["group_name"]] %||% bk[["term_label"]]
  # reformulas expands g/h into h:g and brms into g:h, so without this
  # (1 | h:g) + (1 | g/h), two blocks in brms, would share one name
  if (isTRUE(bk[["from_slash"]])) g <- slash_brms_order(g)
  brms_rename(g)
}

#' A `:` grouping factor or level from a slash, in brms's order.
#'
#' @noRd
slash_brms_order <- function(x) {
  vapply(strsplit(x, ":", fixed = TRUE),
         function(p) paste(rev(p), collapse = ":"), "")
}

#' brms's levels of a block's grouping factor. `combine_groups()` joins
#' the parts of an interaction group with `_`, and `for_r = TRUE` adds
#' `rename_re_levels()`'s whitespace-to-`.` for the `r_` names.
#'
#' @noRd
brms_levels <- function(bk, for_r = FALSE) {
  lev <- bk[["levels"]]
  if (isTRUE(bk[["from_slash"]])) lev <- slash_brms_order(lev)
  if (grepl(":", bk[["group_name"]] %||% "", fixed = TRUE)) {
    lev <- gsub(":", "_", lev, fixed = TRUE)
  }
  if (for_r) lev <- gsub("[ \t\r\n]", ".", lev)
  lev
}

#' The prefix, the renamed coefficient and brms's coefficient name
#' (`get_rnames()`) of each dimension of one random-effect block, in the
#' block's own order. An `|ID|`-merged block spans several components,
#' and each coefficient takes the prefix of the component it came from.
#'
#' @noRd
brms_re_parts <- function(fit, bk) {
  cf <- brms_rename(bk[["cnms"]])
  pre <- character(length(cf))
  for (cp in bk[["components"]]) {
    lp <- fit$frame[["linpreds"]][[cp[["lp_key"]]]]
    if (is.null(lp)) next
    pos <- cp[["offset"]] + seq_len(cp[["dim"]])
    pre[pos] <- brms_lp_prefix(fit, lp)
    cf[pos] <- brms_rename(cp[["cnms"]])
  }
  list(prefix = pre, coef = cf,
       rnames = ifelse(nzchar(pre), paste0(pre, "_", cf), cf))
}

#' brms's coefficient names for one random-effect block; see
#' `brms_re_parts()`.
#'
#' @noRd
brms_re_rnames <- function(fit, bk) brms_re_parts(fit, bk)$rnames

#' Whether a block's sampled `b` segment IS brms's `r_` coefficient: the
#' deviation of each level, laid out level by coefficient. A reduced-rank
#' block samples factor scores, an ICAR field samples its uncentered
#' values, and a smooth, Gaussian-process or SPDE block samples basis
#' weights; brms has no `r_` for any of those.
#'
#' @noRd
brms_block_has_r <- function(bk) {
  !bk[["covstruct"]] %in% c("rr", "smooth", "gp", "hsgp", "car", "spde") &&
    !is.null(bk[["levels"]]) && !is.null(bk[["group_name"]])
}

#' Refuse the models whose brms names would merge two things, each with
#' brms's own message where brms has one:
#'
#' - two responses that `make_stan_names()` spells alike, `y_a` and
#'   `ya` (brms: "Cannot use the same response variable twice");
#' - two random-effect blocks that give one group the same coefficient
#'   (`frame_re()`): brms refuses the model, so there is no brms name for
#'   either copy, and two parameters under one name cannot both be read;
#' - an interaction group two of whose levels `combine_groups()` joins
#'   to one string, `1_2:3` and `1:2_3` both `1_2_3`. brms pastes before
#'   it builds the factor, so it fits the two as ONE level and pools
#'   their effects. Matching that would merge two groups the data keeps
#'   apart, and suffixing would give names brms never has, so the model
#'   is refused by name.
#'
#' @noRd
brms_check_re_dups <- function(fit) {
  rs <- names(fit$spec$responses)
  if (length(rs) > 1L && anyDuplicated(brms_stan_name(rs))) {
    frm_stop("Cannot use the same response variable twice in the same model. ",
             "brms spells a response without '_' and '.', so ",
             paste0("'", rs[brms_stan_name(rs) %in%
                              brms_stan_name(rs)[
                                duplicated(brms_stan_name(rs))]],
                    "'", collapse = " and "),
             " are one name; rename one of the columns", call. = FALSE)
  }
  seen <- character(0)
  for (bk in fit$frame[["re_blocks"]]) {
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "car", "spde")) {
      next
    }
    g <- brms_group_name(bk)
    if (!is.null(bk[["levels"]])) {
      lev <- brms_levels(bk)
      if (anyDuplicated(lev)) {
        same <- bk[["levels"]][lev %in% lev[duplicated(lev)]]
        frm_stop("The levels ", paste0("'", same, "'", collapse = " and "),
                 " of group '", g, "' are one level in brms, which joins the ",
                 "parts of an interaction group with '_', and brms would ",
                 "pool their group-level effects. Recode the factors so that ",
                 "the joined levels differ, for example without '_' in them",
                 call. = FALSE)
      }
    }
    parts <- brms_re_parts(fit, bk)
    key <- paste(g, parts$rnames, sep = "\r")
    hit <- key %in% seen | duplicated(key)
    if (any(hit)) {
      i <- which(hit)[1L]
      frm_stop("Duplicated group-level effects are not allowed.\nOccured ",
               "for effect '", parts$coef[i], "' of group '", g, "'. Give ",
               "one of the two terms its own grouping column, a copy of the ",
               "factor under another name", call. = FALSE)
    }
    seen <- c(seen, key)
  }
  invisible(NULL)
}

#' Labels for a sampled parameter vector in template order, in brms's
#' spelling wherever brms has a parameter with the same content, and in
#' the internal spelling otherwise.
#'
#' - `beta`, `betad`: brms's names from `brms_coef_table()`, including
#'   the natural-scale name of an unmodeled distributional parameter,
#'   whose draws frmtmb.sample stores on that scale.
#' - `b`: `r_<group>[<level>,<coef>]`, with brms's `__<prefix>` after the
#'   group for a non-`mu` or multivariate component, for every block
#'   `brms_block_has_r()` accepts. Any other block keeps `b[<i>]`, the
#'   name the `stanfit` itself carries.
#' - everything else: the internal name, `theta_1`, `thetaac_1`,
#'   `miss_1`, or the family's own extra-parameter names, with
#'   parentheses dropped. brms has no parameter with that content.
#'
#' `include_random = FALSE` drops `b` and `miss`, the layout of draws
#' taken with the random effects integrated out.
#'
#' @noRd
brms_par_labels <- function(fit, include_random = TRUE) {
  brms_check_re_dups(fit)
  tpl <- fit$frame[["par_template"]]
  lab <- brms_coef_names(fit)
  n_beta <- length(tpl[["beta"]])
  out <- character(0)
  for (cp in names(tpl)) {
    if (cp == "beta") {
      out <- c(out, lab[seq_len(n_beta)])
      next
    }
    if (cp == "betad") {
      out <- c(out, lab[-seq_len(n_beta)])
      next
    }
    if (cp %in% c("b", "miss") && !include_random) next
    if (cp == "b") {
      out <- c(out, brms_r_labels(fit))
      next
    }
    v <- names(tpl[[cp]])
    if (is.null(v)) v <- paste0(cp, "_", seq_along(tpl[[cp]]))
    out <- c(out, par_name_bare(v))
  }
  # brms's repair_stanfit(): a label given twice, such as the r_ level
  # names of `lvl 1` and `lvl.1`, takes a __1 suffix on the later copy
  make.unique(out, sep = "__")
}

#' The `b` segment's labels, see `brms_par_labels()`.
#'
#' @noRd
brms_r_labels <- function(fit) {
  n_b <- length(fit$frame[["par_template"]][["b"]])
  lab <- paste0("b[", seq_len(n_b), "]")
  for (bk in fit$frame[["re_blocks"]]) {
    if (!brms_block_has_r(bk)) next
    D <- bk[["dim"]]
    g <- brms_group_name(bk)
    parts <- brms_re_parts(fit, bk)
    lev <- brms_levels(bk, for_r = TRUE)
    k <- seq_along(bk[["b_idx"]]) - 1L
    li <- k %/% D + 1L
    ci <- k %% D + 1L
    pre <- parts$prefix[ci]
    lab[bk[["b_idx"]]] <- paste0("r_", g,
                                 ifelse(nzchar(pre), paste0("__", pre), ""),
                                 "[", lev[li], ",", parts$coef[ci], "]")
  }
  lab
}

#' brms's name for the smoothing standard deviation of each smooth
#' block, `sds_<prefix>_<label>_<k>` (`rename_sm()`), `NA` for any other
#' block. `k` counts the penalties of one smooth term, which frmtmb
#' holds as separate blocks under one label.
#'
#' @noRd
brms_sds_names <- function(fit) {
  blocks <- fit$frame[["re_blocks"]]
  out <- rep(NA_character_, length(blocks))
  seen <- character(0)
  for (bi in seq_along(blocks)) {
    bk <- blocks[[bi]]
    if (!identical(bk[["covstruct"]], "smooth")) next
    pre <- ""
    for (cp in bk[["components"]]) {
      lp <- fit$frame[["linpreds"]][[cp[["lp_key"]]]]
      if (!is.null(lp)) pre <- brms_lp_prefix(fit, lp)
    }
    # brms keeps a by-smooth's ':' in bs_sx:f2u_1 and drops it here,
    # sds_sxf2u_1
    base <- paste0("sds", if (nzchar(pre)) paste0("_", pre), "_",
                   gsub(":", "", brms_rename(bk[["group_name"]] %||%
                                               bk[["term_label"]]),
                        fixed = TRUE))
    seen <- c(seen, base)
    out[bi] <- paste0(base, "_", sum(seen == base))
  }
  out
}
