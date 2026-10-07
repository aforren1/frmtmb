# frmtmb_spec + data -> frmtmb_frame: numeric design matrices, random-effect
# blocks, smooth bases, and the flat parameter template with index
# bookkeeping.
#
# Parameter layout (glmmTMB convention):
#   beta   - fixed effects of every primary (location) linear predictor of
#            every response (REML integrates these), including smooth
#            fixed (null-space) columns
#   betad  - fixed effects of all other dpar linear predictors
#   b      - all random-effect modes; each block level-major; |ID|-merged
#            blocks interleave their component coefficients within a level
#   theta  - all covariance / smoothing parameters
#   thetar - residual-correlation parameters (rescor only)
#
# Every linear predictor's Z spans the FULL b vector (sparse), so Z's
# column indices are b indices; |ID| merging is only a matter of column
# placement. Constant dpars keep their intercept column in betad but are
# excluded from estimation through the MakeADFun `map`.

#' Refuse a group-level coefficient that two terms of one linear
#' predictor both carry for the same grouping factor.
#'
#' `(1 | g) + (x | g)` gives g two independent intercept blocks, and
#' the likelihood sees only the sum of their variances, so the split
#' between them is not identified. lme4 and glmmTMB fit it anyway; brms
#' refuses it, and brms is the tiebreaker. The comparison is on the
#' coefficient names the design actually has, as brms's is, so
#' `(1 | g) + (0 + x | g)`, which splits one block into two without
#' repeating a coefficient, is still accepted. An exact twin term is
#' collapsed to one before this, as brms collapses it
#' (`drop_twin_bar_terms()`), and grouping factors are compared as
#' brms writes them (`slash_nested_bars()`).
#'
#' The glmmTMB covariance structures brms does not have are held to the
#' same rule: `rr(0 + x | g) + diag(0 + x | g)`, `equalto(...) + us(...)`
#' and `ar1(...) + diag(...)` on one factor and coefficient are refused.
#' brms has no verdict on them, and on the base commit they were not
#' identified at 150 groups either, apart from the known-matrix design
#' brms's own rule already covers (dev/reviews/20260916-priorform.md,
#' F3).
#'
#' @noRd
refuse_duplicated_re <- function(cps) {
  seen_key <- character(0)
  seen_label <- character(0)
  seen_known <- logical(0)
  for (cp in cps) {
    known <- isTRUE(cp[["covstruct"]] %in% c("gr_cov", "gr_prec"))
    # the term as written, gr() and a covariance wrapper included, so
    # the two terms the message names can be told apart
    wr <- cp[["written"]] %||% cp[["bar"]]
    lab <- if (is.null(wr)) cp[["label"]] else {
      cs <- cp[["covstruct"]] %||% "us"
      if (cs %in% c("us", "gr_cov", "gr_prec", "us_t")) {
        paste0("(", deparse1(wr), ")")
      } else {
        paste0(cs, "(", deparse1(wr), ")")
      }
    }
    # brms compares grouping factors as the strings it writes, and it
    # writes g/h as g:h where reformulas writes h:g, so only a term that
    # came from a slash is read back in brms's order
    grp <- cp[["group_name"]]
    # brms names a multi-membership group by its members alone
    # (mmg1g2), so weights or scale do not make it a second factor
    if (!is.null(cp[["mm"]])) {
      grp <- paste0("mm", paste(cp[["mm"]][["gvars"]], collapse = ""))
    }
    if (isTRUE(cp[["from_slash"]])) {
      grp <- paste(rev(strsplit(grp, ":", fixed = TRUE)[[1L]]),
                   collapse = ":")
    }
    keys <- paste0(grp, "\r", cp[["cnms"]])
    hit <- match(keys, seen_key)
    if (any(!is.na(hit))) {
      j <- which(!is.na(hit))[1L]
      frm_stop("Duplicated group-level effects are not allowed: the ",
               "coefficient '", cp[["cnms"]][j], "' of group '",
               grp, "' appears in both ",
               seen_label[hit[j]], " and ", lab, ". Only the sum ",
               "of the two variances is identified; write the coefficient ",
               "in one term",
               # an animal model's genetic and permanent-environment terms
               # ARE identified, through the matrix, and brms writes the
               # second on a copy of the column; so does this package
               if (known || seen_known[hit[j]]) {
                 paste0(". A term with a known covariance, gr(cov = ) or ",
                        "gr(prec = ), beside an ordinary one on the same ",
                        "factor is written on a copy of the column, as brms ",
                        "writes it: d$", cp[["group_name"]], "2 <- d$",
                        cp[["group_name"]], " and (1 | ",
                        cp[["group_name"]], "2)")
               }, call. = FALSE)
    }
    seen_key <- c(seen_key, keys)
    seen_label <- c(seen_label, rep(lab, length(keys)))
    seen_known <- c(seen_known, rep(known, length(keys)))
  }
  invisible(NULL)
}

#' Sum of offset() terms of one linear predictor, evaluated against the
#' combined model frame (whose columns are named by deparsed expressions).
#'
#' @noRd
extract_offset <- function(tt, mf, env) {
  oi <- attr(tt, "offset")
  if (is.null(oi) || !length(oi)) return(NULL)
  vars <- attr(tt, "variables")
  off <- 0
  for (i in oi) {
    expr <- vars[[i + 1L]]
    v <- mf[[deparse1(expr)]]
    if (is.null(v)) v <- eval(expr, mf, env)
    off <- off + as.numeric(v)
  }
  off
}

#' The variables of the `offset()` terms of a formula's right-hand side.
#'
#' @noRd
offset_vars <- function(f) {
  walk <- function(e) {
    if (!is.call(e)) return(character(0))
    if (identical(e[[1L]], as.name("offset"))) return(all.vars(e))
    unlist(lapply(as.list(e)[-1L], walk))
  }
  if (is.null(f)) return(character(0))
  unique(as.character(walk(if (inherits(f, "formula")) f[[length(f)]]
                           else f)))
}

#' The variables one dpar needs in the combined model frame: parametric
#' terms, bar variables, and raw smooth variables (never the s() calls
#' themselves, which model.frame cannot evaluate).
#'
#' @noRd
dpar_frame_rhs <- function(dp) {
  if (!is.null(dp[["nl_body"]])) {
    parts <- lapply(dp[["datavars"]], as.name)
    if (!length(parts)) return(1)
    out <- NULL
    for (p in parts) out <- if (is.null(out)) p else call("+", out, p)
    return(out)
  }
  parts <- list(reformulas::RHSForm(dp[["fixed"]]))
  for (rt in dp[["re"]] %||% list()) {
    for (v in by_term_vars(rt)) parts <- c(parts, list(as.name(v)))
    if (is.null(rt$mm)) {
      parts <- c(parts, list(rt$bar[[2]], rt$bar[[3]]))
      next
    }
    # mm(g1, g2) is not a model-frame variable: the member factors are,
    # and so is every mmc() argument. The mmc()-stripped left side goes
    # in as an expression so data-dependent bases still freeze.
    parts <- c(parts, list(rt$mm$lhs), rt$mm$groups)
    for (mc in rt$mm$mmc) parts <- c(parts, mc$exprs)
    for (v in all.vars(rt$mm$weights_expr)) {
      parts <- c(parts, list(as.name(v)))
    }
  }
  for (sspec in dp[["smooth"]] %||% list()) {
    for (tm in sspec$term) parts <- c(parts, list(as.name(tm)))
    if (!is.null(sspec$by) && sspec$by != "NA") {
      parts <- c(parts, list(as.name(sspec$by)))
    }
  }
  for (ent in dp[["mo"]] %||% list()) {
    for (v in all.vars(ent$expr)) parts <- c(parts, list(as.name(v)))
    if (!is.null(ent$mult)) parts <- c(parts, list(ent$mult))
  }
  for (ent in dp[["miterms"]] %||% list()) {
    if (!is.null(ent$mult)) parts <- c(parts, list(ent$mult))
    if (!is.null(ent$idx)) parts <- c(parts, list(ent$idx))
  }
  for (v in me_frame_vars(dp)) parts <- c(parts, list(as.name(v)))
  for (cexpr in dp[["csterms"]] %||% list()) {
    for (v in all.vars(cexpr)) parts <- c(parts, list(as.name(v)))
  }
  for (ge in dp[["gpterms"]] %||% list()) {
    for (ex in c(ge$exprs, if (!is.null(ge$by)) list(ge$by))) {
      for (v in all.vars(ex)) parts <- c(parts, list(as.name(v)))
    }
  }
  # car()/spde() need their grouping variable in the frame; the
  # adjacency and mesh matrices are external data, like gr(cov = A)'s.
  # The RAW variables go in, not the expression - a call-valued gr
  # (gr = factor(node)) is evaluated against the frame at assembly, the
  # way every other special resolves its arguments, and model.frame
  # cannot be asked to carry the call itself.
  for (ce in c(dp[["carterms"]] %||% list(), dp[["spdeterms"]] %||% list())) {
    for (v in all.vars(ce$gr_expr)) parts <- c(parts, list(as.name(v)))
  }
  out <- NULL
  for (p in parts) out <- if (is.null(out)) p else call("+", out, p)
  out
}

#' mo()/mi() interaction multipliers scale a single coefficient, so they
#' have to be one numeric column; a factor or character multiplier would
#' need contrast expansion, which the simplex/latent machinery has no
#' column for. Reject those up front: as.numeric() on a character vector
#' yields all-NA and the failure only surfaces as "NA/NaN gradient
#' evaluation" from the optimizer. `[brms#1828]`
#'
#' @noRd
check_special_mult <- function(mult, expr, fn) {
  if (is.logical(mult)) return(as.numeric(mult))
  if (!is.numeric(mult) || is.factor(mult)) {
    frm_stop(fn, "() interactions support numeric multipliers only: ",
             deparse1(expr), " is ", class(mult)[1L],
             "; expand it to numeric indicator columns first", call. = FALSE)
  }
  as.numeric(mult)
}

#' brms enumerates its special terms in terms() order, so every mo()
#' main effect is numbered before any mo() interaction and simo_1 is not
#' always the first mo() written. Our parser keeps the written order and
#' emits mo(x):z ahead of the mo(x) that `mo(x) * z` implies, so the mo
#' list is re-sorted before the frame hands out simplexes; otherwise
#' zeta<j> and brms's simo_<j> would name different terms. mi() terms
#' are sorted by the same call, for the same reason one step later: it is
#' the order their coefficients are numbered and reported in.
#'
#' Sorting by interaction order alone reproduces terms(): order() is
#' stable, so terms of equal order keep the order they were written in,
#' which is what terms() does with keep.order = FALSE.
#'
#' @noRd
mo_terms_in_brms_order <- function(mo_terms) {
  if (length(mo_terms) < 2L) return(mo_terms)
  ord <- vapply(mo_terms, function(e) mo_interaction_order(e[["mult"]]),
                integer(1))
  mo_terms[order(ord)]
}

#' How many `:`-separated components a term has, counting the special
#' itself. A NULL multiplier is a main effect, order 1.
#'
#' @noRd
mo_interaction_order <- function(mult) {
  if (is.null(mult)) return(1L)
  count <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name(":"))) {
      count(e[[2L]]) + count(e[[3L]])
    } else 1L
  }
  1L + count(mult)
}

#' ar1()/hetar1() correlate two levels by their POSITION in the ordering
#' factor, never by the distance between the labels: with times 1..6 and
#' 10 present, cor(t6, t10) is fitted as rho, not rho^4, and rho itself
#' comes out biased. glmmTMB reads the levels the same way and our
#' agreement tests pin that reading down, so the likelihood stays as it
#' is and the user gets told instead. Silence is right when the labels
#' carry no numbers: position is then the only meaning available.
#' `[glmmTMB#1278]`
#'
#' @noRd
warn_ar1_level_gaps <- function(bar, mf, cs_name) {
  v <- all.vars(bar[[2]])
  if (length(v) != 1L || !is.factor(mf[[v]])) return(invisible(NULL))
  lv <- levels(mf[[v]])
  pos <- suppressWarnings(as.numeric(lv))
  if (anyNA(pos) || any(pos != trunc(pos))) return(invisible(NULL))
  gap <- which(abs(diff(pos)) != 1)
  if (!length(gap)) return(invisible(NULL))
  i <- gap[1L]
  frm_warning(cs_name, "(): the levels of '", v, "' are whole numbers but ",
              "not consecutive ('", lv[i], "' is followed by '", lv[i + 1L],
              "'), and ", cs_name, "() correlates levels by position, so ",
              "that gap counts as a single step. For irregularly spaced ",
              "positions use ou() over num_factor(): ou(num_factor(", v,
              ") + 0 | ...)", call. = FALSE)
  invisible(NULL)
}

#' The two values of a bernoulli response in the order that codes them
#' 0 and 1: brms's `data_response()` takes the levels of
#' `as.factor(y)`, and a numeric response with one value is read as the
#' 1 of 0 and 1 unless that value is 0. A logical is read as the number
#' it is, so an all-`TRUE` response is all successes, where brms's
#' `as.factor(TRUE)` has one level and codes it 0.
#'
#' @noRd
bernoulli_levels <- function(y) {
  if (is.logical(y)) y <- as.numeric(y)
  lv <- levels(as.factor(y))
  if (is.numeric(y) && length(lv) == 1L) {
    lv <- if (lv == "0") c("0", "1") else c("0", lv)
  }
  lv
}

#' A bernoulli response coded 0 and 1 against `levels`, brms's
#' `as.integer(as_factor(y, levels)) - 1`, refused in brms's words when
#' a value is not one of the two.
#'
#' @noRd
bernoulli_code <- function(y, levels, what = "Family 'bernoulli'") {
  if (is.logical(y)) y <- as.numeric(y)
  code <- as.integer(factor(y, levels = levels)) - 1L
  bad <- !is.na(y) & (is.na(code) | code > 1L)
  if (any(bad)) {
    frm_stop(what, " requires responses to contain only two different ",
             "values", if (length(levels) <= 2L) {
               paste0(", and the fitted model codes ",
                      paste0("'", levels, "' as ", seq_along(levels) - 1L,
                             collapse = " and "), ". The value(s) ",
                      paste0("'", unique(as.character(y[bad])), "'",
                             collapse = ", "), " are neither")
             }, call. = FALSE)
  }
  as.numeric(code)
}

#' A bernoulli response read from newdata, on the 0 and 1 codes the fit
#' stores, which is what every quantity that compares it with the draws
#' or the fitted values needs. Any other response is returned as it is.
#'
#' @noRd
response_codes_newdata <- function(rspec, y, what) {
  lv <- rspec$family[["bin_levels"]]
  if (is.null(lv) || is.null(y)) return(y)
  bernoulli_code(y, lv, what = paste0(what, ": family 'bernoulli'"))
}

#' Pull one response out of the combined model frame and coerce it to the
#' numeric form the objective needs. Ordinal factors become category
#' codes and keep their labels in a `y_levels` attribute, a bernoulli
#' response is coded 0 and 1 by level order as brms codes it and keeps
#' the two values in a `bin_levels` attribute, binomial two-level
#' factors become 0/1, and unsupported factor responses or non-finite
#' values are rejected here. A family that already carries `bin_levels`
#' (a refit of a fitted model) is coded against those.
#'
#' @noRd
extract_y <- function(resp, mf) {
  y <- mf[[deparse1(resp$resp_expr)]]
  if (is.null(y)) {
    y <- eval(resp$resp_expr, mf, resp$formula_env)
  }
  lv <- NULL
  bin_lv <- NULL
  if (identical(resp$family[["family"]], "bernoulli") && !is.matrix(y)) {
    uy <- if (is.numeric(y)) unique(y[!is.na(y)])
    if (is.null(resp$family[["bin_levels"]]) && length(uy) &&
          all(uy > 0 & uy < 1)) {
      # coded as brms codes them, but two values inside (0, 1) are
      # usually a proportion, which the coding turns into successes and
      # failures without a word (lme4#682); -0.5 and 0.5 are not
      frm_warning("Family 'bernoulli': the response takes only the ",
                  "values ", paste(sort(uy), collapse = " and "),
                  ", which lie strictly between 0 and 1. They are coded ",
                  "0 and 1 by their order, as brms codes them. If they ",
                  "are proportions, fit them with binomial() and ",
                  "trials(), or with Beta()", call. = FALSE)
    }
    bin_lv <- resp$family[["bin_levels"]] %||% bernoulli_levels(y)
    y <- bernoulli_code(y, bin_lv)
  } else if (identical(resp$family[["type"]], "categorical")) {
    # a nominal response: the level order fixes the reference category
    # and the dpar names, and the codes carry no meaning without the
    # labels, so simulate() and the probability columns can restore them
    known <- resp$family[["cat_levels"]]
    lv <- if (!is.null(known)) {
      # frm() already read the levels off this response (and said so, if
      # coercing a character vector); repeating the message here would
      # print it twice for one fit
      obs <- levels(factor(y))
      if (!all(obs %in% known)) {
        frm_stop("categorical(): the response holds values (",
                 paste(setdiff(obs, known), collapse = ", "),
                 ") that are not among the family's categories (",
                 paste(known, collapse = ", "), ")", call. = FALSE)
      }
      known
    } else {
      categorical_y_levels(y, deparse1(resp$resp_expr))
    }
    if (!is.null(lv)) y <- as.numeric(factor(y, levels = lv))
  } else if (is.factor(y) && identical(resp$family[["type"]], "ordinal")) {
    # An unordered factor's level order is alphabetical unless someone
    # set it, and that order IS the model here. brms 2.23.0 refuses it
    # (dev/famlink-brms-behavior-log.txt); this used to warn and fit.
    code0 <- ord_code0(resp$family)
    if (!is.ordered(y)) {
      frm_stop("Family '", resp$family[["family"]], "' requires either ",
               if (code0 == 0L) "non-negative" else "positive",
               " integers or ordered factors as responses. '",
               deparse1(resp$resp_expr), "' is an unordered factor, whose ",
               "level order (", paste(levels(y), collapse = " < "),
               ") is not a category order anyone stated. Use ",
               "factor(..., levels = ..., ordered = TRUE)", call. = FALSE,
               package = frm_family_package(resp$family))
    }
    # the codes carry no meaning without the labels, and simulate() has
    # to hand draws back in the response's own type
    lv <- levels(y)
    # category codes 1..K in level order; brms codes a hurdle family's
    # first level as the category 0
    y <- as.numeric(y) - 1 + code0
  } else if (is.factor(y)) {
    if (!identical(resp$family[["family"]], "binomial") || nlevels(y) != 2L) {
      frm_stop("Factor responses are only supported for binomial families ",
               "with 2 levels (ordinal families accept ordered factors)",
               call. = FALSE)
    }
    y <- as.numeric(y) - 1
  }
  if (is.matrix(y) && ncol(y) == 1L) {
    y <- drop(y)   # scale() and friends return n x 1 matrices (glmmTMB#937)
  }
  if (is.matrix(y)) {
    storage.mode(y) <- "double"
  } else {
    y <- as.numeric(as.vector(y))
  }
  if (any(!is.finite(y) & !is.na(y))) {
    frm_stop("Non-finite (Inf/NaN) values in the response are not allowed",
             call. = FALSE)
  }
  # attached last: the numeric coercions above drop attributes
  if (!is.null(lv)) attr(y, "y_levels") <- lv
  if (!is.null(bin_lv)) attr(y, "bin_levels") <- bin_lv
  y
}

#' The families that need a number of trials, as brms's `has_trials()`
#' lists them among the families frmtmb has. `multinomial` could read
#' its trials from the row sums, and did; brms requires the term and
#' checks it against the row sums, and so does frmtmb now.
#'
#' @noRd
trials_families <- c("binomial", "beta_binomial", "zero_inflated_binomial",
                     "zero_inflated_beta_binomial",
                     "multinomial")

#' Refuse `se()` on a family whose density does not read it.
#'
#' A capability test, not a name test. `se()` is the one core addition
#' term whose whole effect lives inside the density: the core neither
#' multiplies it in nor reshapes anything with it, it only hands
#' `aterms[["se"]]` to the family and maps out a now redundant `sigma`.
#' So the question is whether the density READS it, and the family
#' answers by declaring it, which is what `accepts_aterms` and
#' `required_aterms` are for. An undeclared family (`accepts_aterms =
#' NULL`, which otherwise means "every term") is refused, because "did
#' not say" cannot mean yes for a term that changes nothing unless it is
#' read: that is the same silent wrong answer the allow-list exists to
#' close.
#'
#' It reads the formula and the family only, so `assemble_frame()` also
#' runs it before the data check when the `se()` column is missing, as
#' brms refuses the term from the formula alone.
#'
#' @noRd
check_se_declared <- function(resp) {
  if (family_declares_aterm(resp$family, "se")) return(invisible(NULL))
  frm_stop(
    "se() carries a known standard deviation into the density, ",
    "so only a family that reads it can be given one, and '",
    resp$family[["family"]], "' does not declare that it does. ",
    "A family declares it with frmtmb_family(accepts_aterms = ",
    "c(..., \"se\")) or with required_aterms = \"se\", which ",
    "also refuses a model that leaves the term out. The density ",
    "then reads aterms[[\"se\"]] as the standard deviation, and ",
    "honors se(x, sigma = TRUE) by reading aterms[[\"se_sigma\"]] ",
    "and using sqrt(dpars[[\"sigma\"]]^2 + aterms[[\"se\"]]^2) ",
    "where it is TRUE. The built-in families that read it are ",
    "gaussian and student", call. = FALSE,
    package = frm_family_package(resp$family))
}

#' Refuse a category-specific `cs(x)` on the left of a group-level bar.
#'
#' `(cs(period) | subject)` reached `model.frame()` as a call to a
#' function nobody defines and died there with R's "could not find
#' function \"cs\"". brms refuses it on a family without category
#' specific effects ("Category specific effects are not supported for
#' this family"), and so does this, in the words the population-level
#' `cs(x)` is refused with. On the ordinal families that do take `cs()`
#' brms fits a group-level category-specific effect; this package has no
#' such term, and says so. The bar form `cs(x | g)`, a compound-symmetry
#' covariance, is a different thing and is not reached here.
#'
#' @noRd
check_cs_in_bar <- function(resp) {
  fam <- resp$family
  cs_ok <- identical(fam[["type"]], "ordinal")
  for (dp in resp$dpars) {
    for (re in dp[["re"]] %||% list()) {
      lhs <- re$bar[[2L]]
      # with a function named cs() on the search path (brms's, attached
      # after frmtmb), parse_linpred() protects the call as .frm_cs()
      # before it gets here, and both spellings are the same term
      if (!calls_function(lhs, "cs") && !calls_function(lhs, ".frm_cs")) {
        next
      }
      re$bar <- str2lang(gsub(".frm_cs(", "cs(", deparse1(re$bar),
                              fixed = TRUE))
      if (!cs_ok) {
        frm_stop("Category specific effects are not supported for this ",
                 "family. cs() needs an ordinal family; the group-level ",
                 "term (", deparse1(re$bar), ") has one", call. = FALSE)
      }
      frm_stop("A category-specific effect inside a group-level term, (",
               deparse1(re$bar), "), is not supported: cs() is a ",
               "population-level term here. Keep cs() outside the bar ",
               "and write the group-level part without it, e.g. ",
               "(1 + x | g)", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Refuse a binomial-type response with no `trials()` term, as brms does.
#'
#' The densities read `trials %||% 1`, so a response without the term
#' used to be fitted as one trial per row. brms refuses it, and a
#' count response of 0 to 10 written without its trials was fitted as a
#' binomial that could not produce most of its own data. A mixture is
#' checked through its components, as brms checks them.
#'
#' @noRd
check_trials_given <- function(resp, av) {
  fam <- resp$family
  fams <- c(fam[["family"]], fam[["component_families"]])
  if (!any(fams %in% trials_families) || !is.null(av[["trials"]])) {
    return(invisible(NULL))
  }
  frm_stop("Specifying 'trials' is required for this model. ",
           paste(intersect(fams, trials_families), collapse = ", "),
           " reads the number of trials from the formula: write ",
           resp$resp_name, " | trials(n) ~ ...",
           if (!identical(fams, "multinomial")) {
             ", or use bernoulli() for a response of zeros and ones"
           }, call. = FALSE)
}

#' Say when a response has only two outcomes and `bernoulli()` would do.
#'
#' brms MESSAGES this, with this sentence, for a binomial-type family
#' whose trials are all one and for an ordinal or categorical response
#' with two categories (dev/famlink-brms-behavior-log.txt). It is a
#' message and not a warning because the model is still correct; it is
#' only doing more work than it has to, and on an ordinal family it
#' estimates one threshold where bernoulli estimates an intercept.
#'
#' Called from `frm()` alone, which is brms's `brm()` and `standata()`:
#' the internal re-assemblies (influence(), the prior table, the
#' simulator) would otherwise say it once per call. A mixture is read
#' through its components, as brms's `has_trials()` reads it.
#'
#' @noRd
suggest_bernoulli <- function(spec, frame) {
  for (resp in spec$responses) {
    fam <- resp$family
    nm <- resp$resp_name
    fams <- c(fam[["family"]], fam[["component_families"]])
    two <- if (any(fams %in% trials_families)) {
      tr <- frame[["aterm_values"]][[nm]][["trials"]]
      length(tr) > 0L && max(tr) == 1
    } else if (isTRUE(fam[["type"]] %in% c("ordinal", "categorical"))) {
      lv <- frame[["y_levels"]][[nm]]
      K <- fam[["cat_K"]] %||% if (length(lv)) length(lv) else {
        max(frame[["y"]][[nm]], na.rm = TRUE)
      }
      identical(as.numeric(K), 2)
    } else {
      FALSE
    }
    if (isTRUE(two)) {
      frm_message("Only 2 levels detected so that family 'bernoulli' might ",
                  "be a more efficient choice.")
    }
  }
  invisible(NULL)
}

#' Data-dependent bases (poly, ns, scale) must be frozen at fit time: the
#' combined model frame's terms carry predvars for every variable; patch
#' them onto a sub-formula's terms by deparsed-variable match before any
#' newdata evaluation (glmmTMB#402 and the largest bug class in
#' lme4/glmmTMB prediction history).
#'
#' @noRd
patch_predvars <- function(tt, map) {
  if (is.null(map) || !length(map)) return(tt)
  vars <- attr(tt, "variables")
  pv <- vars
  changed <- FALSE
  for (i in seq_along(vars)[-1]) {
    key <- deparse1(vars[[i]])
    if (!is.null(map[[key]])) {
      pv[[i]] <- map[[key]]
      changed <- TRUE
    }
  }
  if (changed) attr(tt, "predvars") <- pv
  tt
}

#' The design columns one `cs()` term contributes, and what a `newdata`
#' row is recoded against.
#'
#' brms hands a `cs()` term to `model.matrix()` like any other
#' population-level term, so a factor or character predictor becomes
#' treatment-contrast dummies with one coefficient per dummy per
#' threshold: `cs(fc)` is the pair `fcb`, `fcc`, `cs(fch)` on a
#' character column is `fchb`, `fchc`, and a factor with numeric-looking
#' levels is `fnum2`, `fnum3` (`dev/csfactor-log/brms.txt`). frmtmb
#' called `as.numeric()` on the value instead, so a factor model was
#' fitted on the INTEGER CODES and a character column died in the
#' optimizer.
#'
#' The returned `terms` object carries `predvars`, so a data-dependent
#' basis inside `cs()` (`poly()`, `scale()`, `ns()`) is frozen at the
#' fit as it is for the ordinary design; `xlevels` and `contrasts` are
#' the FIT's, which is what a `newdata` factor is recoded against
#' rather than against its own levels.
#'
#' @noRd
cs_term_design <- function(cexpr, mf, env) {
  tt <- stats::terms(stats::as.formula(call("~", cexpr), env = env))
  mfc <- stats::model.frame(tt, mf, na.action = stats::na.pass)
  tt <- attr(mfc, "terms")
  X <- stats::model.matrix(tt, mfc)
  contr <- attr(X, "contrasts")
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  if (!ncol(X)) {
    frm_stop("cs(", deparse1(cexpr), ") contributes no design column",
             call. = FALSE)
  }
  list(X = X, terms = tt, xlevels = stats::.getXlevels(tt, mfc),
       contrasts = contr, colnames = colnames(X),
       label = deparse1(cexpr))
}

#' The same columns on `newdata`, through the fit's terms, levels and
#' contrasts. A column set that does not match the fit's is a refusal
#' rather than a silent recoding.
#'
#' @noRd
cs_term_newdata <- function(spec, newdata) {
  tt <- spec[["terms"]]
  newdata <- check_newdata_frame(tt, newdata, spec[["xlevels"]])
  mfc <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                            xlev = spec[["xlevels"]])
  X <- stats::model.matrix(tt, mfc, contrasts.arg = spec[["contrasts"]])
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  if (!identical(colnames(X), spec[["colnames"]])) {
    frm_stop("cs(", spec[["label"]], "): newdata gives the design ",
             "column(s) ", paste(colnames(X), collapse = ", "),
             " where the fit has ",
             paste(spec[["colnames"]], collapse = ", "), call. = FALSE)
  }
  X
}

#' The population-level columns a refused `cs()` column is aliased with,
#' named rather than left to the reader. They come from the least-squares
#' weights of the cs column on `cbind(1, X)`, which is exact here because
#' the column IS in that span, so a weight below a thousandth of the
#' largest one is rounding.
#'
#' @noRd
cs_alias_cols <- function(M, z, xnames) {
  if (!length(xnames)) return(character(0))
  co <- tryCatch(stats::setNames(qr.coef(qr(M), z)[-1L], xnames),
                 error = function(e) NULL)
  if (is.null(co) || all(is.na(co))) return(character(0))
  co[is.na(co)] <- 0
  names(co)[abs(co) > 1e-3 * max(abs(co))]
}

#' The same names as a phrase, capped so a wide design does not produce a
#' paragraph.
#'
#' @noRd
cs_alias_phrase <- function(alias, cap = 6L) {
  if (!length(alias)) return("the design as a whole")
  paste0(paste0("'", utils::head(alias, cap), "'", collapse = ", "),
         if (length(alias) > cap) ", ..." else "")
}

#' The edit that actually works, which depends on what the formula holds.
#' Advising "write `x + cs(x)` as `cs(x)`" on `poly(x, 2) + cs(x)` sends
#' the reader to delete a term that is not there and loses the quadratic;
#' the edit that works is `I(x^2) + cs(x)`, which reaches the identical
#' maximum with one parameter fewer (-317.979080807 either way,
#' `dev/csfactor-rev-log/rev-poly.log`). A SMOOTH has no such spelling,
#' because it is the smooth's own unpenalized linear column
#' (`s(x).fx1`) that is aliased, so there one of the two terms has to go.
#' Each worked example stays in ITS OWN branch: the polynomial rewrite
#' pasted into the general case advertised a `poly()` that
#' `x + z + cs(I(x + z))` does not contain.
#'
#' @noRd
cs_alias_advice <- function(term, xnames, alias) {
  if (term %in% (xnames %||% character(0))) {
    return(paste0("Write '", term, " + cs(", term, ")' as 'cs(", term,
                  ")' alone, which fits the same set of distributions, ",
                  "or drop the cs()."))
  }
  if (any(grepl("\\.fx[0-9]+$", alias))) {
    return(paste0("That column is a SMOOTH's unpenalized linear part, ",
                  "and no spelling of the formula keeps both it and the ",
                  "cs() term, so drop one of the two."))
  }
  if (any(grepl("^poly\\(", alias))) {
    return(paste0("That column is a POLYNOMIAL basis's linear part. Drop ",
                  "the cs(), or respell the basis without it: ",
                  "'I(x^2) + cs(x)' in place of 'poly(x, 2) + cs(x)' ",
                  "reaches the same maximum with one parameter fewer."))
  }
  paste0("No population-level term spells '", term, "' on its own, so ",
         "there is nothing to delete by that name: drop either the cs() ",
         "term or the population-level term(s) carrying the column(s) ",
         "named above.")
}

#' The code-indicator basis of one `mo()` term, the space every simplex
#' the objective could choose puts that term's column in.
#'
#' `mo()` contributes `b * D * cumsum(zeta)[codes + 1]`. Over the simplex
#' those vectors span the functions of `codes` that vanish at code 0,
#' which is exactly the treatment-dummy basis of the codes, times the
#' interaction multiplier when the term carries one (`mo(x):z`). `NULL`
#' when the term has no usable codes.
#'
#' @noRd
mo_indicator_basis <- function(mi) {
  codes <- mi[["codes"]]
  D <- mi[["D"]]
  if (!length(codes) || !length(D) || D < 1L) return(NULL)
  B <- outer(as.integer(codes), seq_len(D), "==") * 1
  mult <- mi[["mult"]]
  if (!is.null(mult)) B <- B * as.numeric(mult)
  colnames(B) <- paste0(mi[["label"]], ".", seq_len(D))
  B
}

#' Refuse a `cs()` model whose category-specific columns the data cannot
#' separate from the rest of the predictor.
#'
#' `cs()` gives a column one coefficient per threshold and the
#' thresholds are already one intercept per threshold, so a cs column
#' inside the span of the population-level design plus the constant adds
#' nothing: in `y ~ x + cs(x)` the substitution `b -> b + t`,
#' `bcs[k] -> bcs[k] - t` leaves every row's density unchanged, and the
#' same holds for a factor written on both sides. brms 2.23.0 builds
#' both blocks and samples the ridge (`dev/csfactor-log/brms.txt`, cases
#' B1 and B2); frmtmb refuses, because maximum likelihood has no prior
#' to hold the ridge and what it reported was an estimate with a
#' standard error of 2.7e5 (`dev/csfactor-log/before.txt`).
#'
#' The test is on RANKS, not column counts, because the stored `X`
#' carries zero placeholder columns for `mo()`, `mi()` and `me()` terms
#' that the objective fills later. That is also why the rank test ALONE
#' is blind to those terms, and why `mo()` gets the separate test below.
#' `mi()` and `me()` do not: their columns are observed-or-latent values
#' rather than a function of a data column, so there is no basis to test
#' them against at frame time. `?frm` says so.
#'
#' The rank is `qr()`'s, at its default `tol = 1e-7` on columns scaled to
#' unit norm, so the boundary is a condition number near 1e7: measured, a
#' `cs()` column correlated with a population-level one at `eps = 1e-6`
#' (condition number 2.2e6) fits and `eps = 1e-7` (2.3e7) is refused
#' (`dev/csfactor-rev-log/rev-rank-lane.log`).
#'
#' @noRd
check_cs_identified <- function(X, cs_info, lp_key, mo_info = list()) {
  if (length(cs_info) < 1L) return(invisible(NULL))
  Z <- matrix(unlist(lapply(cs_info, `[[`, "vals"), use.names = FALSE),
              nrow = length(cs_info[[1L]][["vals"]]))
  Xm <- if (is.null(X) || !ncol(X)) NULL else as.matrix(X)
  M <- if (is.null(Xm)) matrix(1, nrow(Z), 1L) else cbind(1, Xm)
  scal <- function(A) {
    s <- sqrt(colSums(A^2))
    s[!(s > 0)] <- 1                 # a zero column keeps its zeros
    sweep(A, 2L, s, "/")
  }
  rnk <- function(A) qr(scal(A))$rank
  base <- rnk(M)
  for (j in seq_len(ncol(Z))) {
    if (rnk(cbind(M, Z[, seq_len(j), drop = FALSE])) >= base + j) next
    col <- sub("^cs", "", cs_info[[j]][["label"]])
    term <- deparse1(cs_info[[j]][["expr"]] %||% str2lang(col))
    # the column names the culprit, the TERM names what to edit: a
    # factor's dummy is not something the formula mentions
    which_col <- if (identical(col, term)) "the column " else {
      paste0("its column '", col, "' ")
    }
    alone <- rnk(cbind(M, Z[, j, drop = FALSE])) < base + 1L
    frm_stop(
      "cs(", term, ") is not identified in '", lp_key, "': ", which_col,
      if (alone && diff(range(Z[, j])) == 0) {
        paste0("is constant over the rows, and an ordinal family ",
               "already has one threshold per category boundary, so ",
               "nothing separates a category-specific coefficient from ",
               "its threshold. Drop the term.")
      } else if (alone) {
        alias <- cs_alias_cols(M, Z[, j], colnames(Xm))
        paste0("is already spanned by the population-level part of the ",
               "predictor (", cs_alias_phrase(alias),
               "), so adding a constant to that coefficient and ",
               "subtracting the same constant from every threshold's ",
               "coefficient leaves the likelihood unchanged. ",
               cs_alias_advice(term, colnames(Xm), alias),
               " brms 2.23.0 builds both blocks and samples the ridge; ",
               "maximum likelihood has no prior to hold it.")
      } else {
        paste0("is a linear combination of the other cs() columns and ",
               "the population-level design, so its coefficients are ",
               "not separately estimable. Drop the redundant term.")
      },
      call. = FALSE)
  }
  # mo() terms. The stored column is a zero placeholder, so the rank
  # loop above cannot see it, but what the objective puts there is
  # D * cumsum(zeta)[codes + 1] for a simplex zeta, and every such vector
  # lies in the span of the code INDICATORS. When the cs() columns plus
  # the rest of the predictor already span that indicator space, no
  # simplex escapes: shifting the mo() coefficient is absorbed exactly.
  # Measured on `yo ~ mo(m) + cs(m)`, m an ordered factor with 4 levels:
  # 3 extra parameters bought 6.4e-09 of log likelihood and every one of
  # 9 standard errors came back NaN (dev/csfactor-rev-log/rev-gap-lane.log,
  # seed 1907). The second condition keeps the refusal ABOUT cs(): a mo()
  # term that the population-level design alone already absorbs is a
  # different and PRE-EXISTING problem, which nothing reports today -
  # `yo ~ mo(m) + m` fits with a negative covariance eigenvalue (-6441)
  # and 6 of 6 non-finite standard errors, in silence. Filed in
  # dev/test-backlog.md; refusing it here would be a change about mo()
  # rather than about cs().
  bz <- rnk(cbind(M, Z))
  for (mi in mo_info) {
    Dm <- mo_indicator_basis(mi)
    if (is.null(Dm)) next
    if (rnk(cbind(M, Z, Dm)) > bz) next
    if (rnk(cbind(M, Dm)) <= base) next
    mnm <- deparse1(mi[["expr"]])
    frm_stop(
      "mo(", mnm, ") and the cs() term(s) of '",
      lp_key, "' are not identified together: cs() already gives every ",
      "value of '", mnm, "' its own coefficient at ",
      "each category boundary, so the monotonic term adds no direction ",
      "the cs() coefficients cannot absorb, whatever its simplex. Keep ",
      "one of the two: mo() alone when the effect is monotone and the ",
      "same at every boundary, cs() alone when it is not. An INTERACTION ",
      "of the two survives, because it is not a function of '", mnm,
      "' alone: write 'z + mo(", mnm, "):z + cs(", mnm, ")' rather than ",
      "'mo(", mnm, ") * z + cs(", mnm, ")', whose main effect is the part ",
      "that is refused.", call. = FALSE)
  }
  invisible(NULL)
}

#' sparse.model.matrix names the columns of matrix-valued terms (poly,
#' ns, scale) by bare index instead of the term-prefixed dense names.
#' Names are load-bearing (frozen param_colnames, coefficient labels), so
#' take them from a one-row dense header over the same frame.
#'
#' This is the sparse counterpart of `stats::model.matrix()` used by
#' `assemble_frame()` when `sparse_x = TRUE`. It takes the terms of one
#' linear predictor, the combined model frame, and the contrasts, and
#' returns a sparse design matrix whose column names match the dense
#' path exactly. It errors if the two paths disagree on the column
#' count, because every later stage indexes coefficients by those names.
#'
#' @noRd
sparse_mm <- function(tt, mf, contrasts.arg = NULL) {
  X <- Matrix::sparse.model.matrix(tt, mf, contrasts.arg = contrasts.arg)
  # keeping the frame's terms attribute on the one-row slice routes
  # model.matrix through deparse-name column matching, never
  # re-evaluating data-dependent bases (poly, scale) on a single row
  mf1 <- mf[1L, , drop = FALSE]
  attr(mf1, "terms") <- attr(mf, "terms")
  hdr <- stats::model.matrix(tt, mf1, contrasts.arg = contrasts.arg)
  if (ncol(hdr) != ncol(X)) {
    frm_stop("Internal error: sparse and dense fixed-effect designs disagree ",
             "on columns; refit without frmtmb_control(sparse_x = TRUE)",
             call. = FALSE)
  }
  colnames(X) <- colnames(hdr)
  X
}

#' Cheap rank screen for a sparse design: |diag(R)| of the sparse QR
#' collapses on aliased columns. The threshold is 100x looser than dense
#' qr()'s 1e-7 default, so anything the dense path would drop gets
#' flagged; a flagged design is re-checked densely, which decides the
#' exact drop set (identical to the dense path by construction).
#'
#' @noRd
sparse_maybe_deficient <- function(X) {
  d <- tryCatch(
    abs(Matrix::diag(Matrix::qrR(Matrix::qr(X), backPermute = FALSE))),
    error = function(e) NULL
  )
  is.null(d) || !all(is.finite(d)) || min(d) <= 0 ||
    min(d) < 1e-5 * max(d)
}

#' A random-effect bar whose left-hand side has no intercept, with one
#' added, or `NULL` when it has one already. `cmc = FALSE` builds the
#' term's design this way and then drops the intercept column, which
#' leaves a factor's treatment contrasts (brms:::validate_terms()).
#'
#' @noRd
bar_with_intercept <- function(bar) {
  tl <- stats::terms(stats::as.formula(call("~", bar[[2L]])))
  if (attr(tl, "intercept") != 0L) return(NULL)
  labs <- attr(tl, "term.labels")
  bar[[2L]] <- str2lang(paste(c("1", labs), collapse = " + "))
  bar
}

#' Whether `cmc = FALSE` changes the columns of a random-effect bar's
#' left-hand side on this model frame: only a factor main effect
#' without an intercept is coded differently, and numeric slopes or an
#' `f:x` without `f` keep their columns. `TRUE` when the columns cannot
#' be built here, so the caller keeps its rewrite and its checks.
#'
#' @noRd
cmc_changes_columns <- function(bar, mf, env) {
  cc <- cmc_bar_columns(bar, mf, env)
  is.null(cc) || !identical(cc$cell, cc$contrast)
}

#' The columns of a random-effect bar's left-hand side under cell-mean
#' coding (`cell`, the default) and under `cmc = FALSE` (`contrast`), or
#' `NULL` when they cannot be built on this model frame. The refusal
#' quotes both, because what `cmc = FALSE` does depends on the factor's
#' contrasts: treatment contrasts leave a level out, and an ordered
#' factor's polynomial contrasts leave no level out but are not levels.
#'
#' @noRd
cmc_bar_columns <- function(bar, mf, env) {
  cols <- function(icpt) {
    tt <- stats::terms(stats::as.formula(call("~", bar[[2L]]), env = env))
    if (icpt) attr(tt, "intercept") <- 1L
    fr <- stats::model.frame(tt, mf, na.action = stats::na.pass)
    setdiff(colnames(stats::model.matrix(tt, fr)),
            if (icpt) "(Intercept)")
  }
  tryCatch(list(cell = cols(FALSE), contrast = cols(TRUE)),
           error = function(e) NULL)
}

#' A column list for a message, cut after `max` names.
#'
#' @noRd
cmc_col_list <- function(x, max = 6L) {
  if (length(x) > max) x <- c(x[seq_len(max)], "...")
  paste(x, collapse = ", ")
}

# Operators that reformulas expands structurally inside a grouping
# expression; everything else on the right of a bar is an ordinary call
# whose value has to exist as a single model-frame column.
grp_structural_ops <- c(":", "/", "*", "+", "%in%")

#' A grouping factor written as a call - `(1 | factor(x))`,
#' `(1 | interaction(a, b))` - is one variable of the combined model frame,
#' stored under its deparsed name. reformulas instead re-evaluates the
#' expression inside the frame, where the call's own arguments (x, a, b)
#' are not columns, and dies with an error raised several frames down.
#' Point the bar at the existing column by name instead; the original bar
#' is what the fit keeps, so prediction still evaluates the expression
#' against newdata. `[lme4#464, #156]`
#'
#' @noRd
resolve_group_calls <- function(bars, fr, env) {
  sub_call <- function(e) {
    if (!is.call(e)) return(e)
    if (as.character(e[[1L]])[1L] %in% grp_structural_ops) {
      for (i in seq_along(e)[-1]) e[[i]] <- sub_call(e[[i]])
      return(e)
    }
    key <- deparse1(e)
    if (!key %in% names(fr)) {
      # defensive: a call nested inside ':' need not be a frame variable
      fr[[key]] <<- eval(e, fr, env)
    }
    as.name(key)
  }
  bars <- lapply(bars, function(b) {
    b[[3L]] <- sub_call(b[[3L]])
    b
  })
  list(bars = bars, fr = fr)
}

#' The pooled level set of a multi-membership term.
#'
#' brms pools the members into ONE grouping factor
#' (`frame_re()`: `unique(ulapply(groups, extract_levels))`), so the
#' levels are each member variable's own levels concatenated in the
#' order the variables were written, deduplicated. A factor contributes
#' `levels()`, anything else the sorted unique values, and a level
#' present in only one member still gets a coefficient.
#'
#' @noRd
mm_pooled_levels <- function(gvals) {
  unique(unlist(lapply(gvals, function(v) {
    if (is.factor(v)) levels(v) else levels(factor(v))
  }), use.names = FALSE))
}

#' Read the member columns of one multi-membership term out of a data
#' frame (the combined model frame at fit time, `newdata` at prediction
#' time).
#'
#' @noRd
mm_member_values <- function(mmspec, data, env) {
  lapply(mmspec$groups, function(g) {
    v <- data[[deparse1(g)]]
    if (is.null(v)) v <- eval(g, data, env)
    if (is.null(v)) {
      frm_stop("mm(): membership variable '", deparse1(g),
               "' is not in the data", call. = FALSE)
    }
    v
  })
}

#' Membership indices and row weights.
#'
#' Returns an `n x J` integer matrix of level indices (`NA` for a level
#' outside `levels`) and an `n x J` weight matrix. The default weights
#' are `1/J` on every row and are NOT rescaled, which is what brms's
#' `data_gr_local()` does: `scale =` only ever touches a supplied
#' weight matrix.
#'
#' @noRd
mm_index_weights <- function(mmspec, data, env, levels) {
  gvals <- mm_member_values(mmspec, data, env)
  n <- length(gvals[[1L]])
  ng <- length(gvals)
  J <- matrix(NA_integer_, n, ng)
  for (k in seq_len(ng)) {
    J[, k] <- match(as.character(gvals[[k]]), levels)
  }
  if (is.null(mmspec$weights_expr)) {
    W <- matrix(1 / ng, n, ng)
  } else {
    W <- eval(mmspec$weights_expr, data, env)
    W <- as.matrix(W)
    if (!identical(dim(W), c(n, ng))) {
      frm_stop("mm(weights = ", deparse1(mmspec$weights_expr),
               "): expected a matrix with one row per observation and one ",
               "column per membership variable (", n, " x ", ng, "), got ",
               nrow(W), " x ", ncol(W),
               ". Build it with cbind(w1, w2)", call. = FALSE)
    }
    storage.mode(W) <- "double"
    if (any(!is.finite(W))) {
      frm_stop("mm(weights = ", deparse1(mmspec$weights_expr),
               "): the weights must all be finite", call. = FALSE)
    }
    if (isTRUE(mmspec$scale)) {
      if (any(W < 0)) {
        frm_stop("mm(scale = TRUE) cannot scale negative weights; pass ",
                 "scale = FALSE to use them as they are", call. = FALSE)
      }
      rs <- rowSums(W)
      if (any(rs == 0)) {
        frm_stop("mm(scale = TRUE): row(s) of the weight matrix sum to ",
                 "zero, so the scaled weights are undefined (first at row ",
                 which(rs == 0)[1L], ")", call. = FALSE)
      }
      W <- W / rs
    }
  }
  list(J = J, W = W, n = n, n_members = ng)
}

#' Per-member design matrices of one multi-membership term.
#'
#' Every member shares the ordinary columns of the bar's left side, and
#' each `mmc()` term contributes ONE column whose values are member
#' specific: member `k` uses `mmc()`'s `k`-th argument. So the returned
#' list holds `J` matrices of identical column count and names, one per
#' member, which is the same encoding brms writes into `Z_..._k`.
#'
#' @noRd
mm_member_designs <- function(mmspec, data, env, n_members,
                              predvar_map = NULL, xlev = NULL,
                              use_model_frame = FALSE) {
  tt <- stats::terms(stats::as.formula(call("~", mmspec$lhs), env = env))
  tt <- patch_predvars(tt, predvar_map)
  # cmc = FALSE: treatment contrasts without the intercept column, as
  # for any other group-level term (brms:::validate_terms())
  cmc_drop <- isFALSE(mmspec[["cmc"]]) && attr(tt, "intercept") == 0L
  if (cmc_drop) attr(tt, "intercept") <- 1L
  Xp <- if (use_model_frame) {
    data <- check_newdata_frame(tt, data, xlev)
    mf2 <- stats::model.frame(tt, data, na.action = stats::na.pass,
                              xlev = xlev)
    stats::model.matrix(tt, mf2)
  } else {
    stats::model.matrix(tt, data)
  }
  if (cmc_drop) Xp <- Xp[, colnames(Xp) != "(Intercept)", drop = FALSE]
  n <- nrow(Xp)
  mmc_vals <- lapply(mmspec$mmc, function(mc) {
    cols <- lapply(mc$exprs, function(ex) {
      v <- data[[deparse1(ex)]]
      if (is.null(v)) v <- eval(ex, data, env)
      if (is.factor(v) || is.character(v)) {
        frm_stop("mmc() requires numeric variables; '", deparse1(ex),
                 "' is a ", if (is.factor(v)) "factor" else "character",
                 " column", call. = FALSE)
      }
      as.numeric(v)
    })
    matrix(unlist(cols, use.names = FALSE), nrow = n)
  })
  cnms <- c(colnames(Xp),
            vapply(mmspec$mmc, `[[`, "", "label"))
  if (!length(cnms)) {
    frm_stop("A multi-membership term needs at least one coefficient: ",
             deparse1(mmspec$lhs), " | ", mmspec$label,
             " has an empty design", call. = FALSE)
  }
  designs <- lapply(seq_len(n_members), function(k) {
    out <- Xp
    for (m in seq_along(mmc_vals)) {
      out <- cbind(out, mmc_vals[[m]][, k])
    }
    colnames(out) <- cnms
    out
  })
  list(designs = designs, cnms = cnms)
}

#' Build the local Z of one multi-membership component: `n` rows by
#' `dim * n_levels` columns, level-major within a level exactly as the
#' single-membership blocks are, so phase 3 places it with no special
#' case. Row `i` puts `w_ik * x_ik` on member `k`'s level block, summed
#' over members, which is why a degenerate `mm(g, g)` reproduces
#' `(1 | g)` bit for bit.
#'
#' @noRd
mm_local_Z <- function(J, W, designs, n_levels) {
  n <- nrow(J)
  D <- ncol(designs[[1L]])
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (k in seq_len(ncol(J))) {
    jk <- J[, k]
    for (cc in seq_len(D)) {
      val <- W[, k] * designs[[k]][, cc]
      keep <- which(!is.na(jk) & val != 0)
      if (!length(keep)) next
      ii <- c(ii, keep)
      jj <- c(jj, (jk[keep] - 1L) * D + cc)
      xx <- c(xx, val[keep])
    }
  }
  # duplicated (i, j) pairs are SUMMED by sparseMatrix(), which is what
  # makes a row that names the same level twice add its two weights
  Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                       dims = c(n, D * n_levels))
}

# Internal censoring codes, shared with brms: -1 left, 0 observed,
# 1 right, 2 interval. The names are the accepted spellings.
cens_code_map <- c(none = 0, left = -1, right = 1, interval = 2)

#' brms lets cens() carry spelled-out codes as well as numbers, and a
#' plain as.numeric() would turn those into NAs with only a coercion
#' warning, so decode before the numeric conversion the other addition
#' terms use. Prefix matching mirrors brms prepare_cens(); unlike brms
#' it is case-insensitive, because "Right" reads as a typo, not garbage.
#'
#' @noRd
decode_cens <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (!is.character(v)) return(as.numeric(v))
  key <- tolower(trimws(v))
  # the numeric codes the error message advertises can arrive as
  # strings (a character column, or a factor built from one); prefix
  # matching would reject every one of them, so read those directly.
  # Out-of-range numbers fall through to the -1/0/1/2 check, which
  # names the offending value
  num <- suppressWarnings(as.numeric(key))
  idx <- vapply(key, function(k) {
    hit <- if (nzchar(k)) which(startsWith(names(cens_code_map), k)) else
      integer(0)
    if (length(hit) == 1L) hit else NA_integer_
  }, integer(1L), USE.NAMES = FALSE)
  bad <- unique(v[is.na(idx) & is.na(num)])
  if (length(bad)) {
    frm_stop("cens() cannot decode: ",
             paste0("\"", bad, "\"", collapse = ", "),
             "; use \"none\", \"left\", \"right\", or \"interval\" ",
             "(any unambiguous prefix), or the codes 0, -1, 1, 2",
             call. = FALSE)
  }
  out <- unname(cens_code_map[idx])
  ifelse(is.na(num), out, num)
}

#' Normalize and check the `data2` argument of [frm()].
#'
#' `data2` must be a named list; anything else is a user mistake that
#' would otherwise surface far away, inside one of the structural
#' lookups. Takes the user's value and returns a named list, with `NULL`
#' becoming an empty one. Errors on any other shape.
#'
#' @noRd
validate_data2 <- function(data2) {
  if (is.null(data2)) return(list())
  nms <- names(data2)
  if (!is.list(data2) || (length(data2) && is.null(nms)) ||
      any(!nzchar(nms))) {
    frm_stop("data2 must be a named list, e.g. data2 = list(W = W)",
             call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    frm_stop("data2 has duplicate names: ",
             paste(unique(nms[duplicated(nms)]), collapse = ", "),
             call. = FALSE)
  }
  data2
}

#' Resolve a structural object named by a special term.
#'
#' Takes the unevaluated argument of a term such as `car(M)`,
#' `gr(cov =)` or `spde(fem)`, together with `data2`, `data` and the
#' formula environment, and returns the object it names.
#'
#' Structural objects (adjacency matrices, precisions, covariance
#' matrices, FEM triples) resolve from `data2` before anything else, so
#' a fit that names them there is self-contained: `saveRDS()` carries
#' them on the fit and a later refit never reaches back into the calling
#' environment that built them, which by then may be gone.
#'
#' brms accepts a bare name in `data2` and nothing more. This keeps that
#' rule and adds one: a compound expression is evaluated with `data2` in
#' front of the data mask, so `car(solve(P))` finds `P` in `data2` too.
#' The historical data-then-formula-env evaluation stays as the
#' fallback, so models written before `data2` existed keep working.
#'
#' @noRd
lookup_structural <- function(expr, data2, data, env, what) {
  e2 <- NULL
  if (length(data2)) {
    if (is.symbol(expr)) {
      nm <- as.character(expr)
      if (nm %in% names(data2)) return(data2[[nm]])
    }
    # the wrapper list separates "evaluated to NULL" from "failed"
    val <- tryCatch(
      list(v = eval(expr, data2, list2env(as.list(data), parent = env))),
      error = function(e) {
        e2 <<- e
        NULL
      }
    )
    if (!is.null(val)) return(val$v)
  }
  tryCatch(eval(expr, data, env), error = function(e) {
    # the data2-mask attempt saw the widest scope, so when both paths
    # fail its error names the real cause (the fallback just repeats
    # "not found" for objects that only exist in data2)
    frm_stop(structural_lookup_msg(expr, data2, what,
                                   if (is.null(e2)) e else e2),
             call. = FALSE)
  })
}

#' Do two resolved structural matrices describe the same thing?
#'
#' Compared on the RESOLVED objects rather than on the deparsed
#' expressions: `|ID|`-linked terms live in different `bf()` formulas,
#' and two formula environments can bind the same name to different
#' matrices. Both sides are already reordered onto the block's levels,
#' so a plain elementwise comparison answers the question. The tolerance
#' is relative to the larger of the two, because a relationship matrix
#' assembled twice (a pedigree recomputed per formula) can differ in the
#' last bits without describing a different structure.
#'
#' @noRd
same_structural_matrix <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  if (!identical(dim(a), dim(b))) return(FALSE)
  scale <- max(1, suppressWarnings(max(abs(a))))
  d <- suppressWarnings(max(abs(a - b)))
  isTRUE(is.finite(d) && d <= 1e-10 * scale)
}

#' Message for a structural object that a lookup could not find.
#'
#' Takes the unevaluated expression, the `data2` list, the term name and
#' the caught condition, and returns one string that names what was
#' looked for, what `data2` holds, and the underlying error. It exists
#' so that every structural term reports a miss the same way.
#'
#' @noRd
structural_lookup_msg <- function(expr, data2, what, e) {
  txt <- deparse1(expr)
  held <- if (length(data2)) {
    paste0(" data2 holds: ", paste(names(data2), collapse = ", "), ".")
  } else {
    " data2 is empty."
  }
  if (is.symbol(expr)) {
    paste0(what, ": cannot find '", txt, "'. Pass structural objects ",
           "in data2, e.g. frm(..., data2 = list(", txt, " = ", txt,
           ")); a fit whose matrices come from data2 also survives ",
           "saveRDS() and refits in a new session.", held)
  } else {
    paste0(what, ": cannot evaluate '", txt, "' (", conditionMessage(e),
           "). Put the objects it needs in data2, e.g. ",
           "frm(..., data2 = list(...)).", held)
  }
}

#' Is this object something `model.frame()` could never hold as a
#' column?
#'
#' The set is decided by `model.frame()` itself, not by taste. Its C
#' code accepts a vector type - logical, integer, double, complex,
#' character, raw - and everything built on one, so a factor, a `Date`,
#' a `POSIXct`, a `difftime` and a matrix are all legal columns. Every
#' other type is refused with `invalid type (...) for variable 'x'`:
#' lists (a data.frame and a `POSIXlt` are lists), functions,
#' environments, and language objects (a formula is a call). Those are
#' exactly the objects that cannot be a column reference, so a body name
#' that resolves to one is a reference to the object.
#'
#' A vector of any type is deliberately not in the set. A numeric or
#' character vector in the formula environment is a legal model-frame
#' variable, and `model.frame()` already resolves it through the
#' formula environment when `data` has no such column; taking it
#' lexically instead would change a fit that works today. A matrix is
#' not in the set either: matrix columns are a real feature here
#' (smooths and functional terms), so a matrix has to keep reaching the
#' frame.
#'
#' `NULL` is not in the set. `get0()` returns `NULL` for a name that is
#' absent from the environment too, and that name has to stay in the
#' frame so the failure is `model.frame()`'s "object 'x' not found"
#' rather than a coercion error later.
#'
#' @noRd
nl_lexical_only <- function(obj) {
  if (is.null(obj)) return(FALSE)
  is.function(obj) || is.list(obj) || is.environment(obj) ||
    is.language(obj)
}

#' Drop nonlinear-body names that are objects, not data.
#'
#' `parse_one_response()` collects every symbol of the nonlinear body
#' that is not a nonlinear parameter and asks the combined model frame
#' for it. That is right for `pk_ode(exp(lka), time, dose)`, where only
#' the arguments are symbols, and wrong as soon as a helper takes
#' another object as an argument - `solve_pk(pk_dyn, ...)` would ask
#' `model.frame()` for a column named `pk_dyn`, and
#' `solve_pk(..., events = doses)` for a column named `doses`. The body
#' is evaluated in its own formula environment anyway, so a name that
#' is not a column of `data` and resolves there to something
#' `model.frame()` could never accept as a column of any data - see
#' `nl_lexical_only()` - is left to resolve lexically.
#'
#' A column always wins over a same-named object in the environment, so
#' a variable called `t`, `c` or `data` is unaffected, and a list column
#' named in a body still fails at `model.frame()` the way an unusable
#' column should.
#'
#' The names left to resolve lexically are recorded, because they are
#' also the suspects when the body later fails. A misspelled column that
#' shares a name with a base function (`t`, `c`, `df`) is no longer
#' caught by `model.frame()`, so the symptom moves to a coercion error
#' deep in the objective; `nl_body_error()` turns that back into a
#' message that names the candidates.
#'
#' @noRd
drop_nl_lexical_datavars <- function(spec, data) {
  dn <- names(data)
  for (i in seq_along(spec$responses)) {
    for (nm in names(spec$responses[[i]]$dpars)) {
      dp <- spec$responses[[i]]$dpars[[nm]]
      if (is.null(dp[["nl_body"]]) || !length(dp[["datavars"]])) next
      keep <- vapply(dp[["datavars"]], function(v) {
        if (v %in% dn) return(TRUE)
        obj <- tryCatch(get0(v, envir = dp[["nl_env"]], ifnotfound = NULL),
                        error = function(e) NULL)
        !nl_lexical_only(obj)
      }, TRUE)
      spec$responses[[i]]$dpars[[nm]]$datavars <- dp[["datavars"]][keep]
      spec$responses[[i]]$dpars[[nm]]$nl_lexical <- dp[["datavars"]][!keep]
    }
  }
  spec
}

#' Decide which names in a nonlinear body are references to another
#' distributional parameter's value.
#'
#' `nl_dpar()` marks every body name that matches one of the family's
#' dpars as a CANDIDATE reference and leaves it in `datavars` as well,
#' because the formula alone cannot tell the two apart. A column of the
#' data wins: that keeps every body brms accepts meaning exactly what it
#' means there (in brms a body name is a column or a nonlinear
#' parameter, never a dpar). Only a name with no column behind it reads
#' the other parameter's per-row value, which is what a variance
#' function of the model's own mean needs -
#' `nlf(sigma ~ ls + th * log(abs(mu)))`, nlme's
#' `varPower(form = ~ fitted(.))`.
#'
#' Dropping the resolved names from `datavars` also keeps them out of
#' the combined model frame, which would otherwise ask `model.frame()`
#' for a column called `mu` and fail before the body ever runs.
#'
#' @noRd
resolve_nl_dpar_refs <- function(spec, data) {
  dn <- names(data)
  for (i in seq_along(spec$responses)) {
    for (nm in names(spec$responses[[i]]$dpars)) {
      dp <- spec$responses[[i]]$dpars[[nm]]
      refs <- dp[["nl_dpar_refs"]] %||% character(0)
      if (!length(refs)) next
      refs <- setdiff(refs, dn)
      spec$responses[[i]]$dpars[[nm]]$nl_dpar_refs <- refs
      spec$responses[[i]]$dpars[[nm]]$datavars <-
        setdiff(dp[["datavars"]], refs)
    }
  }
  spec
}

#' Refuse a nonlinear body that names its OWN parameter.
#'
#' `nl_dpar()` excludes a parameter's own name from `nl_dpar_refs` -
#' a body cannot read the value it is itself computing - so the name
#' stays in `datavars` and the combined model frame is asked for a
#' column of that name. With no such column that surfaced as R's own
#' "object 'mu' not found" from `eval(predvars, data, env)`, naming
#' neither the parameter nor the reason.
#'
#' A real column of that name still wins, exactly as it does for a dpar
#' reference: `bf(I ~ mu - chi * logw, chi ~ 1, nl = TRUE)` on data
#' carrying a `mu` column is a body reading that column, and is left
#' alone. The refusal is therefore made HERE, where the data is known,
#' rather than at parse time. `parse_one_response()` catches the other
#' half of the collision - a name given both a formula and a body - and
#' needs no data to do it.
#'
#' @noRd
check_nl_self_reference <- function(spec, data) {
  dn <- names(data)
  for (resp in spec$responses) {
    for (nm in names(resp$dpars)) {
      dp <- resp$dpars[[nm]]
      if (is.null(dp[["nl_body"]])) next
      if (!nm %in% (dp[["datavars"]] %||% character(0))) next
      if (nm %in% dn) next
      reserved <- resp$family[["dpars"]] %||% character(0)
      frm_stop("The body of '", nm, "' refers to '", nm, "' itself, and ",
               "`data` has no column of that name. A nonlinear body is ",
               "computed FROM its parameters, so it cannot read the ",
               "parameter it computes.",
               if (nm %in% reserved) {
                 paste0(" '", nm, "' is a distributional parameter of family '",
                        resp$family[["family"]], "': a nonlinear parameter ",
                        "cannot be named after one, so rename it. This ",
                        "family reserves: ",
                        paste(reserved, collapse = ", "), ".")
               },
               " Declare the parameter under another name with ",
               "bf(..., a ~ 1, nl = TRUE), or add the column to `data`",
               call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Re-raise a nonlinear-body failure with the lexical names attached.
#'
#' @noRd
nl_body_error <- function(e, lp) {
  lex <- lp[["nl_lexical"]] %||% character(0)
  extra <- if (length(lex)) {
    paste0(" These names in the body are not columns of `data` and were ",
           "resolved outside the data, in the formula environment: ",
           paste(lex, collapse = ", "),
           ". A misspelled column name that happens to match an object ",
           "there - a function such as t, c or df - fails exactly this ",
           "way.")
  } else {
    ""
  }
  frm_stop("The nonlinear formula body could not be evaluated: ",
           conditionMessage(e), extra, call. = FALSE)
}

#' The language objects a structured family needs in the model frame:
#' a grouping column, a time column, a sequence id. Read off the family
#' alone, before any data is seen.
#'
#' @noRd
structure_frame_vars <- function(fam) {
  st <- fam_structure(fam)
  fv <- st[["frame_vars"]]
  if (is.null(fv)) return(list())
  out <- fv(fam)
  if (is.null(out)) list() else Filter(Negate(is.null), as.list(out))
}

#' Model-frame columns that are not predictors: responses, and the
#' grouping variables of random-effect, structured-family, `car()` and
#' `spde()` terms.
#'
#' @noRd
nonpredictor_frame_vars <- function(spec) {
  out <- character(0)
  for (resp in spec$responses) {
    out <- c(out, deparse1(resp$resp_expr))
    # both spellings: a bare name is its own column, and a compound
    # expression reaches the frame under its deparsed form
    for (ex in structure_frame_vars(resp$family)) {
      out <- c(out, all.vars(ex), deparse1(ex))
    }
    out <- c(out, all.vars(resp$autocor$time_expr),
             all.vars(resp$autocor$gr_expr))
    for (dp in resp$dpars) {
      for (rt in dp[["re"]] %||% list()) {
        out <- c(out, if (is.null(rt$mm)) deparse1(rt$bar[[3L]]) else
                        rt$mm$gvars, by_term_vars(rt))
      }
      for (ce in c(dp[["carterms"]] %||% list(),
           dp[["spdeterms"]] %||% list())) {
        out <- c(out, all.vars(ce$gr_expr))
      }
    }
  }
  unique(out)
}

#' Report date and time columns, which reach the design as a number.
#'
#' `model.matrix()` reduces a `Date`, a `POSIXct` or a `difftime` to its
#' underlying number, so the model that is fitted is the right one, but
#' it is expressed in an origin and a unit the user did not choose. The
#' coefficient is per day (`Date`) or per second (`POSIXct`), and the
#' intercept is the fitted value at 1970-01-01, tens of thousands of
#' units away from any modern data. That extrapolated intercept is what
#' breaks the fit: a `Date` predictor of about 18000 makes the objective
#' badly conditioned, and the same model on days-since-the-first-day
#' converges where the raw column reports false convergence.
#'
#' Only predictors are reported. The combined model frame also holds the
#' responses and the grouping factors, and neither has the problem the
#' message describes: a response is converted by an explicit
#' `as.numeric()` and a location shift of it is absorbed by the
#' intercept, and a grouping variable is used for its distinct levels,
#' where a `Date` behaves exactly like a factor. Advice about
#' coefficients and intercepts on either of those would be wrong.
#'
#' It is a `message()`, not a warning: the coercion is deliberate and
#' correct in meaning, and `suppressMessages()` silences it for a caller
#' who has already centered the column or wants the epoch origin.
#'
#' @srrstats {G2.5} `Date`, `POSIXct` and `difftime` predictors are
#'   accepted, and the class they are silently reduced to is surfaced
#'   rather than left to be discovered: frame assembly names each such
#'   column, states the unit and the origin its coefficient will be in,
#'   and points at centering. `vignette("inputs")` documents the same in
#'   its "Predictor classes" table.
#' @srrstats {G2.9} This is the package's diagnostic for a conversion
#'   that loses information. `Date`, `POSIXct` and `difftime` are the
#'   only predictor classes `model.matrix()` reduces to a bare number,
#'   dropping the class, the calendar and the unit; the number that
#'   survives is measured from an origin (1970-01-01) that the user did
#'   not choose and that is far outside modern data. Rather than let
#'   that be discovered through a non-converging fit, assembly emits one
#'   message per fit naming each affected column, the unit and origin
#'   its coefficient and intercept will be expressed in, and the
#'   centering that avoids the loss of conditioning. No other accepted
#'   class is reduced this way: factors keep their contrasts, matrix
#'   terms keep their frozen basis, and a response is converted by an
#'   explicit documented `as.numeric()`.
#' @srrstats {G2.4d} Conversion to `factor` happens in two places, both
#'   deliberate, both documented and both tested. A grouping variable
#'   that is not already a factor is converted so that its level set can
#'   be read, which is why an integer, a character and a factor spelling
#'   of the same grouping variable give the same fit
#'   (`tests/testthat/test-edgecases.R`); and the `cluster` argument of
#'   [vcov_cluster()] is converted and dropped of empty levels. No
#'   predictor is converted to build a design matrix: factor structure
#'   there belongs to the user's data and is resolved by
#'   `stats::model.frame()` and `stats::model.matrix()` with the stored
#'   contrasts, and creating one silently would change the model.
#' @noRd
report_datetime_columns <- function(mf, exclude = character(0)) {
  kind <- function(v) {
    if (inherits(v, "Date")) "Date, days since 1970-01-01"
    else if (inherits(v, "POSIXt")) "POSIXct, seconds since 1970-01-01"
    else if (inherits(v, "difftime")) {
      paste0("difftime, ", attr(v, "units") %||% "unknown units")
    } else NA_character_
  }
  mf <- mf[setdiff(names(mf), exclude)]
  if (!length(mf)) return(invisible(NULL))
  hits <- vapply(mf, kind, "")
  hits <- hits[!is.na(hits)]
  if (!length(hits)) return(invisible(NULL))
  frm_message(
    "Date/time column", if (length(hits) > 1L) "s" else "",
    " used as ", if (length(hits) > 1L) "numbers" else "a number", ": ",
    paste0(names(hits), " (", hits, ")", collapse = ", "),
    ". The coefficient is per unit of that origin and the intercept is ",
    "the value at it, which is far outside the data and can stop the ",
    "optimizer converging. Center the column, for example ",
    "as.numeric(x - min(x)), to put the intercept back in range."
  )
  invisible(NULL)
}

#' Column permutation taking the smooth's own basis to the
#' `smooth2random()` layout, for smooths that have no `trans.U`.
#'
#' `smooth2random(type = 2)` reports the split of a t2 basis as
#' `pen.ind` (which penalty each original column belongs to, 0 = the
#' unpenalized null space) instead of a `trans.U` rotation, so the
#' mapping is a scaling by `trans.D` followed by a stable sort that puts
#' penalty 1, penalty 2, ... first and the null-space columns last. The
#' identity is verified here on a few rows rather than asserted, so a
#' smooth class that reaches this path with some other meaning of
#' `pen.ind` returns NULL and prediction refuses instead of quietly
#' returning wrong numbers.
#'
#' @noRd
smooth_pen_order <- function(sm, re2) {
  if (!is.null(re2$trans.U)) return(NULL)
  pen <- re2$pen.ind
  if (is.null(pen) || is.null(re2$trans.D) ||
      length(pen) != ncol(sm$X) || length(re2$trans.D) != ncol(sm$X)) {
    return(NULL)
  }
  ord <- order(ifelse(pen == 0L, max(pen) + 1L, pen), seq_along(pen))
  rows <- seq_len(min(nrow(sm$X), 50L))
  Z <- cbind(do.call(cbind, re2$rand), re2$Xf)[rows, , drop = FALSE]
  if (ncol(Z) != ncol(sm$X)) return(NULL)
  rec <- sweep(sm$X[rows, , drop = FALSE], 2, re2$trans.D, `*`)
  if (max(abs(rec[, ord, drop = FALSE] - Z)) > 1e-8 * max(1, max(abs(Z)))) {
    return(NULL)
  }
  ord
}

#' Move the intercept column of a `0 + Intercept` design to where brms
#' puts it: after the columns of the `pos` terms written before
#' `Intercept`. model.matrix() orders columns by term, and the rewritten
#' formula keeps the other terms in their written order, so `assign`
#' says which columns those are. See rsv_intercept_fixed().
#'
#' @noRd
rsv_intercept_order <- function(X, pos) {
  asg <- attr(X, "assign")
  contr <- attr(X, "contrasts")
  j <- match("(Intercept)", colnames(X))
  if (pos < 1L || is.na(j) || is.null(asg)) return(X)
  before <- setdiff(which(asg >= 1L & asg <= pos), j)
  ord <- c(before, j, setdiff(seq_len(ncol(X)), c(before, j)))
  out <- X[, ord, drop = FALSE]
  attr(out, "assign") <- asg[ord]
  attr(out, "contrasts") <- contr
  out
}

#' brms's deprecated lower-case `intercept` in a formula without an
#' intercept is a data column of ones (`data_rsv_intercept()`), filled
#' in here, and in newdata by `pred_design()`. A column the data already
#' carry must be all ones, as brms requires.
#'
#' @noRd
rsv_lower_fill <- function(spec, data) {
  lower <- any(vapply(spec$responses, function(r) {
    any(vapply(r$dpars, function(dp) isTRUE(dp[["rsv_lower"]]), NA))
  }, NA))
  if (!lower || is.environment(data) || !length(data)) return(data)
  x <- data[["intercept"]]
  if (!is.null(x) && (!(is.numeric(x) || is.logical(x)) || anyNA(x) ||
                        any(x != 1))) {
    frm_stop("Variable name 'intercept' is reserved in models without a ",
             "population-level intercept.", call. = FALSE)
  }
  data[["intercept"]] <- rep(1, NROW(data[[1L]]))
  data
}

#' Refuse a `0 + Intercept` model whose data carry an `Intercept` (or
#' `intercept`) column that is not all ones.
#'
#' brms fills both names with ones in such a model and refuses data that
#' say otherwise ("Variable name 'Intercept' is reserved in models
#' without a population-level intercept", `data_rsv_intercept()`), so a
#' column the user meant as a covariate is never silently replaced by
#' the intercept. frmtmb never reads the column; the refusal is kept so
#' that a model brms refuses is not quietly given a meaning here.
#'
#' @noRd
check_rsv_intercept_data <- function(spec, data) {
  rsv <- any(vapply(spec$responses, function(r) {
    any(vapply(r$dpars, function(dp) !is.null(dp[["rsv_intercept"]]), NA))
  }, NA))
  if (!rsv || is.environment(data)) return(invisible(NULL))
  for (v in intersect(c("Intercept", "intercept"), names(data))) {
    x <- data[[v]]
    if (!(is.numeric(x) || is.logical(x)) || anyNA(x) || any(x != 1)) {
      frm_stop("`data` has a column `", v, "` that is not all ones, and ",
               "the model writes 0 + Intercept, where `", v, "` is the ",
               "reserved name of the intercept's column of ones; brms ",
               "refuses it too. Rename the column", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Whether `v` is a variable of `data`: a column of a data frame or a
#' named list, or an object of an environment passed as data.
#'
#' @noRd
in_data_var <- function(v, data) {
  if (is.environment(data)) exists(v, envir = data) else v %in% names(data)
}

#' Refuse, by name, an addition-term variable that is not in `data`, as
#' brms's `validate_data()` refuses it ("can neither be found in 'data'
#' nor in 'data2'"). A value from the formula environment would be read
#' again by every prediction on newdata and every refit, and after
#' saveRDS() and readRDS() in another session it is missing or changed;
#' on the response side, as `trials(k)`, the change is silent. Function
#' calls stay allowed: only the free variables are checked.
#'
#' @noRd
check_aterm_data_vars <- function(resp, data) {
  for (nm_at in names(resp$aterms)) {
    a <- resp$aterms[[nm_at]]
    if (is.numeric(a) || is.logical(a) || is.character(a)) next
    miss <- Filter(function(v) !in_data_var(v, data), all.vars(a))
    if (!length(miss)) next
    v <- miss[[1L]]
    env <- resp$formula_env %||% globalenv()
    fn_hint <- if (is.function(tryCatch(get(v, envir = env),
                                        error = function(e) NULL))) {
      paste0(" (R finds only the function ", v, "() from the formula)")
    } else ""
    frm_stop("Addition term ", aterm_label(nm_at, a), " of response '",
             resp$resp_name, "' reads `", v, "`, which is not a column of ",
             "`data`", fn_hint, ". An addition term reads its variables ",
             "from the data alone, as brms does: a value from the ",
             "formula environment would be read again, possibly changed, ",
             "by predictions and refits. Put `", v, "` in the data",
             call. = FALSE)
  }
  invisible(NULL)
}

#' The variables of an expression that belong in the model frame when
#' the expression itself cannot be a frame column: a column of `data`,
#' or an object of the formula environment with one value per row. `a`
#' is the expression or the names of its variables.
#'
#' @noRd
expr_frame_vars <- function(a, data, env) {
  vars <- if (is.character(a)) a else all.vars(a)
  n <- if (is.data.frame(data)) nrow(data)
  env <- env %||% globalenv()
  Filter(function(v) {
    if (if (is.environment(data)) exists(v, envir = data)
        else v %in% names(data)) {
      return(TRUE)
    }
    x <- tryCatch(get(v, envir = env), error = function(e) NULL)
    if (is.function(x)) {
      # weights(t * 2) without a column t found base::t(), and the
      # arithmetic then failed without naming the variable
      frm_stop("The model uses `", v, "`, which is not a column of ",
               "`data`; R finds only the function ", v, "() from the ",
               "formula. Add the column to `data` or correct the name",
               call. = FALSE)
    }
    # an object that is not found goes to the frame so that the frame's
    # own check names it; a constant such as `k` in weights(wt * k) is
    # read when the term is evaluated, since it has no row to drop
    is.null(x) || (!is.null(n) && (is.atomic(x) || is.factor(x)) &&
                     NROW(x) == n)
  }, vars)
}

#' Refuse, by name, the variables stats::model.frame() would refuse in
#' its own words: a name that neither `data` nor the formula
#' environment holds, a list used as a variable, and an object from the
#' environment whose length is not the number of rows.
#'
#' The lookup is model.frame()'s: `data` first, then `env` and its
#' parents. An environment `data` is where model.frame() evaluates, so
#' there the parents of `data` are searched and `env` is not. Only a
#' whole term is tested for type and length, because a list or a scalar
#' inside a call such as `I(x / k)` is legitimate.
#'
#' @noRd
check_frame_variables <- function(rhs, data, env) {
  data_env <- is.environment(data)
  in_data <- function(v) {
    if (data_env) exists(v, envir = data) else v %in% names(data)
  }
  if ("." %in% all.vars(rhs)) {
    # frm() expands `.` against a data frame before the formula is
    # parsed (expand_dot_bform()), so one that reaches here had no
    # columns to stand for
    frm_stop("The formula has a `.`, which stands for the columns of ",
             "`data` that the formula does not otherwise use, and this ",
             "`data` has no columns to expand it into. Pass a data frame ",
             "as `data`, or write the predictors out, e.g. y ~ x1 + x2",
             call. = FALSE)
  }
  for (v in all.vars(rhs)) {
    if (!in_data(v) && (data_env || !exists(v, envir = env))) {
      # brms's reserved name reaches here only where it is not reserved,
      # so say where it is
      hint <- if (v %in% c("Intercept", "intercept")) {
        paste0(". `Intercept` is the reserved name of the intercept only ",
               "in a population-level formula without one, as in ",
               "y ~ 0 + Intercept + x; elsewhere it is an ordinary ",
               "variable. A group-level intercept is (1 | g)")
      } else ""
      frm_stop("The model uses `", v, "`, which is not a column of `data` ",
               "and not an object that R finds from the formula. Add the ",
               "column to `data` or correct the name", hint, call. = FALSE)
    }
  }
  leaves <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L) {
      c(leaves(e[[2L]]), leaves(e[[3L]]))
    } else if (is.name(e)) {
      as.character(e)
    }
  }
  n <- if (is.data.frame(data)) nrow(data)
  for (v in unique(leaves(rhs))) {
    from_data <- in_data(v)
    x <- if (!from_data) {
      get(v, envir = env)
    } else if (data_env) {
      get(v, envir = data)
    } else {
      data[[v]]
    }
    if (is.list(x)) {
      where <- if (from_data) "a column of `data`" else {
        "found from the formula"
      }
      frm_stop("`", v, "` is a list, ", where, ", and a model variable ",
               "must be a vector, a factor or a matrix", call. = FALSE)
    }
    if (!from_data && !is.null(n) && (is.atomic(x) || is.factor(x)) &&
          NROW(x) != n) {
      frm_stop("`", v, "` is not a column of `data`; the object R finds ",
               "from the formula has ", NROW(x), " rows and `data` has ", n,
               ". Add the column to `data` or correct the name",
               call. = FALSE)
    }
  }
  invisible(NULL)
}

#' The raw variables the model reads through a transform, such as `x` in
#' `poly(x, 2)` or `z` in `log(abs(z) + 1)`, on the rows of the model
#' frame. The model frame keeps the transformed columns only; brms's
#' keeps the variables, and `conditional_effects()` varies a variable,
#' not a basis matrix. A variable read only inside `offset()` is left
#' out: the offset has its own handling in the display. `NULL` when the
#' model frame already holds every variable, or when its rows cannot be
#' matched to the data's by name.
#'
#' @noRd
frame_raw_vars <- function(mf, data, rhs, env = NULL) {
  if (is.null(rhs) || !is.data.frame(data)) return(NULL)
  off_only <- character(0)
  other <- character(0)
  walk <- function(e, in_off) {
    if (is.name(e)) {
      nm <- as.character(e)
      if (in_off) off_only <<- c(off_only, nm) else other <<- c(other, nm)
    } else if (is.call(e)) {
      off <- in_off || identical(e[[1L]], as.name("offset"))
      for (a in as.list(e)[-1L]) walk(a, off)
    }
  }
  walk(rhs, FALSE)
  v <- setdiff(unique(other), names(mf))
  rows <- match(rownames(mf), rownames(data))
  if (anyNA(rows)) return(NULL)
  # a predictor keeps R's lookup in the formula environment (user
  # decision of 2026-09-30), so a variable found there, one value per
  # row of the data, is a raw variable as a column is; a function or a
  # constant of another length is not
  from_env <- list()
  if (!is.null(env)) {
    for (w in setdiff(v, names(data))) {
      x <- tryCatch(get(w, envir = env), error = function(e) NULL)
      if ((is.numeric(x) || is.factor(x) || is.character(x) ||
             is.logical(x)) && is.null(dim(x)) && length(x) == nrow(data)) {
        from_env[[w]] <- x[rows]
      }
    }
  }
  v <- intersect(v, names(data))
  if (!length(v) && !length(from_env)) return(NULL)
  out <- data[rows, v, drop = FALSE]
  for (w in names(from_env)) out[[w]] <- from_env[[w]]
  rownames(out) <- rownames(mf)
  out
}

#' Turn a parsed spec plus data into the numeric `frmtmb_frame` the
#' objective is built from. This is the second and last stage of the
#' formula-to-design-matrix pipeline: `parse_spec()` reads the formulas,
#' `assemble_frame()` reads the data.
#'
#' It builds one combined model frame over every response, every dpar
#' variable and every addition-term variable, so `na.action` drops rows
#' from all of them together. It then evaluates the responses and the
#' addition terms, builds the fixed-effect design of each linear
#' predictor (dropping aliased columns the way `lm()` does), and turns
#' every random-effect, smooth, Gaussian-process, spatial, monotonic,
#' missing-data and category-specific term into a component. Components
#' sharing an `|ID|` key merge into one block, blocks are laid out in
#' the flat `b` and `theta` vectors, and each linear predictor gets a
#' sparse `Z` over the whole `b` vector.
#'
#' The returned object carries the designs (`linpreds`), the
#' random-effect blocks (`re_blocks`), the response and addition-term
#' values, the `par_template` with the MakeADFun `map`, and everything
#' prediction needs to reproduce the designs on newdata: terms,
#' xlevels, contrasts, frozen predvars, and the model frame itself.
#' Structural matrices (adjacency, precision, covariance, FEM triples)
#' resolve through `lookup_structural()`: data2 first, then data, then
#' the formula environment.
#'
#' `check_response = FALSE` is for the prior table and its validator
#' alone. brms's `get_prior()` and `validate_prior()` never read the
#' response's values, because no prior slot depends on them: they answer
#' for a binomial-type response whatever its trials say, and for a
#' `Beta()` response outside (0, 1) (`tests.priors.R`, "auxiliary
#' parameters", draws `rnorm()`). So the table skips the trials refusal
#' and `valid_y()`, except on an ordinal family, whose threshold count
#' is read off the response. A fit never passes it.
#'
#' `thres_pin` is for a refit INSIDE the package that reassembles the
#' frame from a subset or a replacement of the fitted data
#' (`influence()`). It carries the fitted model's ordinal threshold
#' counts and category labels, from `thres_pin_of_fit()`. A factor
#' response is recoded against the fitted categories first
#' (`thres_pin_recode()`), because `drop.unused.levels = TRUE` renumbers
#' them when a subset leaves a level empty; the count is then written
#' into the addition-term values as if `thres(x = )` had been in the
#' formula (`thres_pin_apply()`), so that a subset whose top category is
#' absent refits the FITTED model rather than a shorter one. A count the
#' user wrote wins over the pinned count, but not over the recoding.
#'
#' @noRd
assemble_frame <- function(spec, data, na.action = stats::na.omit,
                           sparse_x = FALSE, data2 = list(),
                           check_response = TRUE, thres_pin = NULL,
                           drop_unused_levels = TRUE) {
  # `data = NULL` is not "no data": model.frame() falls back to the
  # formula environment and reports the first variable it cannot find
  # there ("object 'y' not found"), which sends the reader looking for a
  # typo in the formula rather than for the missing argument. A
  # data.frame, a tibble, a data.table and a plain named list all reach
  # model.frame() unchanged and are all supported.
  if (is.null(data)) {
    frm_stop("`data` is NULL: frm() needs the data frame holding the model ",
             "variables, e.g. frm(bf(y ~ x) + gaussian(), data = d)",
             call. = FALSE)
  }
  if (!is.data.frame(data) && !is.environment(data) &&
        !(is.list(data) && !is.null(names(data)))) {
    frm_stop("`data` must be a data frame holding the model variables, ",
             "not ", arg_desc(data), call. = FALSE)
  }
  data2 <- validate_data2(data2)
  data <- rsv_lower_fill(spec, data)
  check_nl_self_reference(spec, data)
  spec <- resolve_nl_dpar_refs(spec, data)
  spec <- drop_nl_lexical_datavars(spec, data)
  # One combined model frame holds every response, every variable of every
  # dpar of every response, and the aterm variables, so na.omit keeps rows
  # aligned. Responses go on the RHS of the frame formula (a multivariate
  # LHS is not a valid model.frame response) and are extracted by name.
  rhs_comb <- NULL
  # what each response adds, for subset()'s NA rule (R/subset.R)
  rhs_resp <- list()
  add_part <- function(part) {
    rhs_comb <<- if (is.null(rhs_comb)) part else call("+", rhs_comb, part)
    rhs_resp[[cur_resp]] <<- if (is.null(rhs_resp[[cur_resp]])) {
      part
    } else call("+", rhs_resp[[cur_resp]], part)
  }
  subset_check_spec(spec, na.action)
  for (resp in spec$responses) {
    cur_resp <- resp$resp_name
    add_part(resp$resp_expr)
    for (dp in resp$dpars) {
      add_part(dpar_frame_rhs(dp))
      # an offset(log(time)) column alone leaves `time` out of the
      # frame, and the conditional-effects and emmeans grids are built
      # from the frame's columns, where brms's grid holds `time` itself
      for (v in expr_frame_vars(offset_vars(dp[["fixed"]]), data,
                                resp$formula_env)) {
        add_part(as.name(v))
      }
    }
    # the time and grouping variables of a residual correlation term
    # live on the RESPONSE (the term is not part of any dpar's design),
    # so they are added here rather than in dpar_frame_rhs()
    for (v in c(all.vars(resp$autocor$time_expr),
                all.vars(resp$autocor$gr_expr))) {
      add_part(as.name(v))
    }
    # a structured family's own frame variables (a class grouping, the
    # time and grouping variables that define an HMM's sequences)
    # belong to the RESPONSE and to no dpar's design
    for (ex in structure_frame_vars(resp$family)) add_part(ex)
    for (nm_at in names(resp$aterms)) {
      # literal constants (e.g. trials(10)) are not frame variables, and
      # interval bounds (cens_y2) may be NA on non-interval rows, so they
      # stay out of the na.omit frame (brms#1070)
      a <- resp$aterms[[nm_at]]
      if (nm_at == "cens_y2" || is.numeric(a) || is.logical(a)) next
      # brms evaluates an addition term's expression on the data rows
      # (get_ad_values()), so `weights(wt * 2)` is `wt` times 2. Only
      # its variables enter the frame, never the expression: in a
      # formula `wt * 2` is an interaction and `s / 2` a nesting, and a
      # summary such as `min(y)` has one value where a frame column
      # needs one per row. A variable that is not a column of `data` is
      # refused below (check_aterm_data_vars())
      for (v in all.vars(a)) {
        if (in_data_var(v, data)) add_part(as.name(v))
      }
    }
  }
  # Whether a family takes se() is a fact about the formula, and brms
  # refuses it from the formula alone, so a missing se() column must not
  # hide that refusal behind "not a column of data". Only then is it
  # raised here: with the column present the order is the one below, in
  # which a structured family's own refusal of se() speaks first
  # (frmtmb.latent's hmm() says why it cannot take the term).
  for (resp in spec$responses) {
    se_expr <- resp$aterms[["se"]]
    if (!is.null(se_expr) &&
          !all(vapply(all.vars(se_expr), function(v) {
            v %in% names(data) ||
              exists(v, envir = resp$formula_env %||% globalenv())
          }, NA))) {
      check_se_declared(resp)
    }
    check_cs_in_bar(resp)
    check_aterm_data_vars(resp, data)
  }
  env <- spec$responses[[1]]$formula_env
  fr_formula <- stats::as.formula(call("~", rhs_comb), env = env)
  check_frame_variables(rhs_comb, data, env)
  check_rsv_intercept_data(spec, data)
  # x | mi() responses may carry NAs (they become latent parameters);
  # rows are dropped only for NAs in every OTHER variable. A structured
  # family that declares `keep_na` reads the NAs itself and takes the
  # same exemption: an hmm() NA is a time point the chain passes through
  # without emitting, and an lca(na.rm = FALSE) NA masks one item out of
  # that subject's likelihood. Either way the row must survive, because
  # dropping it changes the estimand rather than the sample.
  mi_cols <- vapply(
    Filter(function(r) {
      isTRUE(r$aterms[["mi"]]) ||
        isTRUE(fam_structure(r$family)[["keep_na"]])
    }, spec$responses),
    function(r) deparse1(r$resp_expr), ""
  )
  has_sub <- spec_has_subset(spec)
  if (length(mi_cols) || has_sub) {
    mf <- stats::model.frame(fr_formula, data = data,
                             drop.unused.levels = drop_unused_levels,
                             na.action = stats::na.pass)
    bad <- if (has_sub) {
      # brms's na_omit(): an NA is harmless on a row that every
      # response using the variable leaves out
      subset_na_rows(spec, mf, lapply(rhs_resp, frame_rhs_columns), mi_cols)
    } else {
      Reduce(`|`, lapply(setdiff(names(mf), mi_cols), function(cn) {
        v <- mf[[cn]]
        if (is.matrix(v)) rowSums(is.na(v)) > 0 else is.na(v)
      }), rep(FALSE, nrow(mf)))
    }
    if (any(bad)) {
      dropped <- which(bad)
      mf <- mf[!bad, , drop = FALSE]
      attr(mf, "na.action") <- structure(
        dropped,
        class = if (identical(na.action, stats::na.exclude)) "exclude"
                else "omit"
      )
    }
  } else {
    mf <- stats::model.frame(fr_formula, data = data,
                             drop.unused.levels = drop_unused_levels,
                             na.action = na.action)
  }
  # Dropping rows changes the estimand and the n every later standard
  # error is built on, so the loss is reported instead of being inferred
  # from nobs(). One message per fit; suppressMessages() silences it.
  n_dropped <- length(attr(mf, "na.action"))
  if (n_dropped > 0L) {
    frm_message(n_dropped, if (n_dropped == 1L) " row" else " rows",
                " removed because of missing values (na.action)")
  }
  n <- nrow(mf)
  if (n == 0L) {
    # zero rows in and zero rows left are different faults, and the
    # generic "after removing NAs" wording sends the second one hunting
    # for missing values that were never there
    frm_stop(if (n_dropped > 0L) {
               "No complete observations after removing NAs"
             } else {
               "`data` has no rows; nothing to fit"
             }, call. = FALSE)
  }
  # subset(): the rows each response uses. A univariate model has one
  # set of rows, so the frame itself is cut to them and every later
  # stage sees an ordinary frame.
  sub_rows <- if (has_sub) subset_rows_of(spec, mf) else list()
  # brms's nobs(): the rows of the data, before any subset() cut
  n_data <- nrow(mf)
  uni_rows <- NULL
  if (has_sub && length(spec$responses) == 1L) {
    uni_rows <- sub_rows[[1L]]
    mf <- frame_rows(mf, uni_rows)
    n <- nrow(mf)
    sub_rows <- list()
  }
  if (!has_sub && anyNA(mf[setdiff(names(mf), mi_cols)])) {
    frm_stop("NA values remain in the model variables after applying ",
             "na.action; use na.omit (default) or na.exclude", call. = FALSE)
  }
  report_datetime_columns(mf, exclude = nonpredictor_frame_vars(spec))
  # freeze data-dependent bases: map deparsed variable -> predvar call
  predvar_map <- local({
    tt_all <- attr(mf, "terms")
    vars <- as.list(attr(tt_all, "variables"))[-1]
    pv <- as.list(attr(tt_all, "predvars") %||%
                    attr(tt_all, "variables"))[-1]
    stats::setNames(pv, vapply(vars, deparse1, ""))
  })

  y <- list()
  y_levels <- list()
  aterm_values <- list()
  extras <- list()
  extra_map <- list()  # per response: family extra name -> template name
  mi_map <- list()   # per mi() response: missing rows + miss indices
  n_miss <- 0L
  miss_init <- numeric(0)
  blocks <- list()   # per response: the structured family's data block
  autocor <- list()  # per response: R-side residual correlation block
  n_thetaac <- 0L
  index_vals <- list()   # per index() response: its rows' index values
  mf_all <- mf
  n_all <- n
  for (resp in spec$responses) {
    # a subset() response sees only its own rows from here on
    mf <- frame_rows(mf_all, sub_rows[[resp$resp_name]])
    n <- nrow(mf)
    if (!is.null(resp$aterms[["index"]])) {
      index_vals[[resp$resp_name]] <- index_eval(resp, mf)
    }
    # A name that is both a nonlinear parameter and a data column is
    # ambiguous, and the nonlinear body resolves it to the PARAMETER,
    # silently ignoring the column - the fit runs and reports numbers
    # for a model the user did not write. Refuse rather than document a
    # precedence nobody would remember. [brms#391, #734]
    if (length(resp$nlpars)) {
      clash <- intersect(resp$nlpars, names(data))
      if (length(clash)) {
        frm_stop("Nonlinear parameter(s) ",
                 paste0("'", clash, "'", collapse = ", "),
                 " also name columns of the data. The nonlinear formula ",
                 "would use the parameter and ignore the column; rename ",
                 "one of them", call. = FALSE)
      }
    }
    yv0 <- extract_y(resp, mf)
    lv0 <- attr(yv0, "y_levels")
    attr(yv0, "y_levels") <- NULL   # nothing on the tape carries labels
    bin_lv <- attr(yv0, "bin_levels")
    attr(yv0, "bin_levels") <- NULL
    if (!is.null(bin_lv) && is.null(resp$family[["bin_levels"]])) {
      # brms stores the coding with the fit (frame$basis$resp_levels), so
      # a refit on a subset holding one of the two values, and a response
      # read from newdata, are coded as the fit was; the family is what
      # the fit's spec carries to both
      resp$family[["bin_levels"]] <- bin_lv
      spec$responses[[resp$resp_name]] <- resp
    }
    # FIRST, before anything reads the codes: a refit inside the package
    # codes its response against the FITTED model's categories, because
    # `drop.unused.levels = TRUE` above renumbers them when a subset
    # leaves a level empty
    rc <- thres_pin_recode(thres_pin, resp, yv0, lv0)
    y_levels[[resp$resp_name]] <- rc$levels
    y[[resp$resp_name]] <- rc$y
    at_names <- setdiff(names(resp$aterms),
                        c("cens_y2", "se_sigma", "mi", row_aterms))
    av <- stats::setNames(lapply(at_names, function(nm_at) {
      a <- resp$aterms[[nm_at]]
      v <- mf[[deparse1(a)]]
      if (is.null(v)) {
        v <- eval(a, mf, resp$formula_env)
        # brms recycles a single value, as in trunc(lb = min(y) - 1); a
        # literal such as trunc(lb = -5) stays one value, as it always has
        if (is.call(a) && length(all.vars(a))) {
          if (length(v) == 1L) v <- rep(v, n)
          if (anyNA(v)) {
            # the variables are complete by here, so the expression made
            # the NA; brms refuses it at standata(), and a row dropped
            # for it would change the sample in silence
            frm_stop("Addition term ", aterm_label(nm_at, a), " of ",
                     "response '", resp$resp_name, "' is NA on ",
                     sum(is.na(v)), " of ", n, " rows, where its ",
                     "variables are not. Give every row a value",
                     call. = FALSE)
          }
          if (!is.matrix(v) && length(v) != n) {
            frm_stop("Addition term ", aterm_label(nm_at, a), " of ",
                     "response '", resp$resp_name, "' has ", length(v),
                     " values where the data have ", n, " rows; it takes ",
                     "one value per row, or a single value", call. = FALSE)
          }
        }
      }
      if (nm_at == "cens") return(decode_cens(v))
      if (nm_at == "thres_gr") return(thres_group_codes(v))
      # a registered term brings its own coercion, which is the point of
      # registering one: the spelling a literature uses (a factor, a
      # two-level character) becomes the numbers the density indexes
      reg_at <- registered_aterm_of(nm_at)
      if (is.null(reg_at)) return(as.numeric(v))
      v <- reg_at$coerce(v)
      if (!is.numeric(v)) {
        frm_stop("The coercion registered for `", reg_at$name,
                 "()` returned ", arg_desc(v), "; an addition term's value ",
                 "is baked into the tape as data and must be numeric",
                 call. = FALSE)
      }
      as.numeric(v)
    }), at_names)
    if (!is.null(resp$aterms[["se_sigma"]])) {
      av[["se_sigma"]] <- resp$aterms[["se_sigma"]]   # logical flag, not data
    }
    # before every guard and before family_finalize(), which is what
    # reads the count: an in-package refit's threshold count comes from
    # the fitted model, never from the refit's own response
    av <- thres_pin_apply(thres_pin, resp, av)
    if (check_response) check_trials_given(resp, av)
    # Before EVERY other guard, including the structured one, because
    # each of them is handed `av`: a declared term that is absent leaves
    # a hole in it, and a hole reads as NULL rather than as an error.
    # The failure this replaces is silent - NULL in the density's
    # arithmetic gives numeric(0), the log-likelihood sums over nothing,
    # and the fit RETURNS, with a log-likelihood of zero.
    have_at <- unique(c(names(av), names(resp$aterms)))
    # Each group is a set of spellings the density reads ANY of, and a
    # plain character vector makes every entry a group of one, so the
    # conjunction the argument has always meant is unchanged.
    miss_at <- Filter(
      function(g) !any(g %in% have_at),
      required_aterm_groups(resp$family[["required_aterms"]]))
    # Plain requirements first, so a message carrying both reads
    # "`vreal1`, one of `dec` or `vint1`" rather than trailing a bare
    # name off the end of a choice. A message of plain ones only is
    # unaffected, which is what every family declaring the old spelling
    # gets.
    miss_at <- miss_at[order(lengths(miss_at) > 1L)]
    if (length(miss_at)) {
      needs <- vapply(miss_at, function(g) {
        if (length(g) == 1L) paste0("`", g, "`")
        else paste0("one of ", paste0("`", g, "`", collapse = " or "))
      }, "")
      # The example writes the FIRST alternative of each unmet group.
      # Which others there are is in the sentence above it; a formula
      # has to pick one spelling to be a formula at all.
      spell <- vapply(miss_at, function(g) aterm_spelling(g[[1L]]), "")
      frm_stop(resp$family[["family"]], ": the density needs ",
               paste(needs, collapse = ", "),
               ", which nothing on this response supplies. Write the ",
               "addition term: ", resp$resp_name, " | ",
               paste(spell, collapse = " + "),
               " ~ ...", call. = FALSE,
               package = frm_family_package(resp$family))
    }
    # The same declaration read the other way round: the check above
    # asks whether the datum arrived, this one whether two spellings of
    # it arrived together. An allow-list cannot catch that, because both
    # spellings are legitimately on it.
    check_exclusive_aterms_supplied(resp, av)
    # before the generic aterm guards, so a structured response is
    # refused for the shape of its likelihood rather than for its
    # family's missing CDF further down
    st_ <- fam_structure(resp$family)
    if (!is.null(st_[["check_spec"]])) st_[["check_spec"]](resp, spec, av)
    if (any(c("thres", "thres_gr") %in% names(resp$aterms)) &&
          !identical(resp$family[["type"]], "ordinal")) {
      # brms's own sentence: the term is refused for any family without
      # thresholds, and a custom family that declares no allow-list
      # would otherwise read nothing and fit as if it were absent
      frm_stop("thres() is not a valid addition term for family '",
               resp$family[["family"]], "': it sets the number of ",
               "thresholds of an ordinal family, and this family has ",
               "none. The ordinal families are cumulative(), sratio(), ",
               "cratio() and acat()", call. = FALSE)
    }
    if (isTRUE(resp$aterms[["mi"]])) {
      if (!resp$family[["family"]] %in% c("gaussian", "student")) {
        frm_stop("mi() responses need a gaussian or student model",
                 call. = FALSE)
      }
      if (any(c("cens", "trunc_lb", "trunc_ub", "se") %in%
                names(resp$aterms))) {
        frm_stop("mi() cannot be combined with cens(), trunc(), or se() ",
                 "on the same response", call. = FALSE)
      }
      if (spec$rescor) {
        frm_stop("mi() cannot be combined with rescor = TRUE", call. = FALSE)
      }
      yv <- y[[resp$resp_name]]
      if (is.matrix(yv)) {
        frm_stop("mi() responses must be numeric vectors", call. = FALSE)
      }
      if (!is.null(av[["mi_sd"]])) {
        # measurement error (brms me()): every true value is latent;
        # observed values get a N(latent, sd) term in the objective
        if (any(av[["mi_sd"]] <= 0)) {
          frm_stop("mi(sd): measurement SDs must be positive", call. = FALSE)
        }
        obs <- which(!is.na(yv))
        if (!length(obs)) {
          frm_stop("mi(sd): the response has no observed values",
                   call. = FALSE)
        }
        rows <- seq_along(yv)
        mi_map[[resp$resp_name]] <- list(
          rows = rows, idx = n_miss + rows,
          se = av[["mi_sd"]], obs = obs
        )
        n_miss <- n_miss + length(rows)
        miss_init <- c(miss_init,
                       ifelse(is.na(yv), mean(yv[obs]), yv))
        yv[is.na(yv)] <- 0
        y[[resp$resp_name]] <- yv
      } else {
        rows <- which(is.na(yv))
        if (length(rows)) {
          mi_map[[resp$resp_name]] <- list(rows = rows,
                                           idx = n_miss + seq_along(rows))
          n_miss <- n_miss + length(rows)
          miss_init <- c(miss_init,
                         rep(mean(yv[-rows]), length(rows)))
          yv[rows] <- 0
          y[[resp$resp_name]] <- yv
        }
      }
    }
    if (!is.null(resp$aterms[["cens_y2"]])) {
      v <- as.numeric(eval(resp$aterms[["cens_y2"]], data, resp$formula_env))
      nd_ <- if (is.data.frame(data)) nrow(data) else if (is.list(data)) {
        NROW(data[[1L]])
      }
      if (!is.null(nd_) && length(v) != nd_) {
        # brms: "Argument 'y2' needs to have length equal to the number
        # of data rows"; it is not recycled, unlike the censoring code
        frm_stop("cens(", deparse1(resp$aterms[["cens"]]), ", ",
                 deparse1(resp$aterms[["cens_y2"]]), "): the interval ",
                 "upper bound has ", length(v), " value(s) where the data ",
                 "have ", nd_, " rows. It takes one value per row, as ",
                 "brms requires", call. = FALSE)
      }
      if (!is.null(attr(mf, "na.action"))) {
        v <- v[-attr(mf, "na.action")]
      }
      # the rows this response reads, under subset()
      rows_ <- sub_rows[[resp$resp_name]] %||% uni_rows
      if (!is.null(rows_)) v <- v[rows_]
      av[["cens_y2"]] <- v
    }
    for (vn in grep("^vint", names(av), value = TRUE)) {
      if (any(av[[vn]] != round(av[[vn]]))) {
        frm_stop(vn, " values must be integers (use vreal() for reals)",
                 call. = FALSE)
      }
    }
    if (!is.null(av[["rate"]]) && (anyNA(av[["rate"]]) ||
                                     any(av[["rate"]] <= 0))) {
      # brms's own check, "Rate denomiators should be positive": the
      # exposure enters as log(denom)
      frm_stop("rate(", deparse1(resp$aterms[["rate"]]), "): rate ",
               "denominators should be positive", call. = FALSE)
    }
    if (!is.null(av[["weights"]])) {
      if (any(av[["weights"]] < 0)) {
        frm_stop("weights() must be non-negative", call. = FALSE)
      }
      if (spec$rescor) {
        frm_stop("weights() cannot be combined with rescor = TRUE",
                 call. = FALSE)
      }
    }
    # the joint-gaussian rescor likelihood has no censoring,
    # truncation, or known-se machinery; without this guard those
    # terms were silently dropped (wrong likelihood, no warning)
    if (spec$rescor &&
        (!is.null(av[["cens"]]) || !is.null(av[["trunc_lb"]]) ||
         !is.null(av[["trunc_ub"]]) || !is.null(av[["se"]]))) {
      frm_stop("cens()/trunc()/se() cannot be combined with rescor = TRUE",
               call. = FALSE)
    }
    if (!is.null(av[["cens"]]) || !is.null(av[["trunc_lb"]]) ||
        !is.null(av[["trunc_ub"]])) {
      # A family that declares only `lccdf` scores a RIGHT-censored row
      # exactly and has nothing to say about a left bound, so it is
      # admitted for right censoring alone. Everything else here still
      # needs `lcdf`.
      right_only <- is.null(av[["trunc_lb"]]) && is.null(av[["trunc_ub"]]) &&
        !is.null(av[["cens"]]) && all(av[["cens"]] %in% c(0, 1))
      if (is.null(resp$family[["lcdf"]]) &&
          !(right_only && !is.null(resp$family[["lccdf"]]))) {
        frm_stop("cens()/trunc() need a family with a CDF (currently: ",
                 "gaussian, lognormal, poisson, exponential, weibull, ",
                 "inverse.gaussian, cox). The list is not closed: a family ",
                 "supplies one through the lcdf argument of ",
                 "frmtmb_family(), and a family that only ever sees RIGHT ",
                 "censoring may supply the log survivor function through ",
                 "lccdf instead", call. = FALSE,
                 package = frm_family_package(resp$family))
      }
      if (!is.null(av[["cens"]]) && !all(av[["cens"]] %in% c(-1, 0, 1, 2))) {
        frm_stop("cens() codes must be -1 (left), 0 (observed), 1 (right), ",
                 "or 2 (interval), or the matching names \"left\", \"none\", ",
                 "\"right\", \"interval\"; got: ",
                 paste(unique(av[["cens"]][!av[["cens"]] %in% c(-1, 0, 1, 2)]),
                       collapse = ", "), call. = FALSE)
      }
      # The discrete censoring convention (see row_lpdf() in
      # R/objective.R, where it is the arithmetic). A bound on a count
      # NAMES a value the response can take and is INCLUDED: right
      # censoring at k is Y >= k, an interval is k <= Y <= k2, and left
      # censoring at k is Y <= k. That is the rule trunc(lb = ) already
      # follows, so one number means one thing on a response however it
      # is bounded, and it is what "5 or more" means where a count was
      # recorded that way.
      #
      # It reads a lower edge as F(edge - 1), so it assumes the support
      # is the unit integer lattice. A family whose support is not that
      # would be shifted onto a point it has no mass at, silently, so a
      # non-integer edge is refused here rather than rounded there.
      if (!is.null(av[["cens"]]) &&
          identical(resp$family[["type"]], "discrete")) {
        yv <- y[[resp$resp_name]]
        i_c <- which(av[["cens"]] != 0)
        edges <- c(yv[i_c],
                   if (!is.null(av[["cens_y2"]])) {
                     av[["cens_y2"]][av[["cens"]] == 2]
                   })
        edges <- edges[!is.na(edges)]
        if (length(edges) && any(edges != round(edges))) {
          frm_stop("cens() on a discrete family reads every bound as a ",
                   "value the response can take: right censoring at k is ",
                   "Y >= k, an interval is k <= Y <= k2, and a lower bound ",
                   "enters the CDF as F(k - 1). That step assumes the ",
                   "support is the integers, so a censoring bound must be ",
                   "one; got ",
                   paste(utils::head(unique(edges[edges != round(edges)]), 3),
                         collapse = ", "), call. = FALSE)
        }
      }
      if (!is.null(av[["cens"]]) && any(av[["cens"]] == 2)) {
        i2 <- av[["cens"]] == 2
        if (is.null(av[["cens_y2"]])) {
          frm_stop("Interval censoring (code 2) needs upper bounds: ",
                   "cens(c, y2)", call. = FALSE)
        }
        if (anyNA(av[["cens_y2"]][i2])) {
          frm_stop("cens() upper bounds must not be NA on interval-censored ",
                   "rows", call. = FALSE)
        }
        yv <- y[[resp$resp_name]]
        # A discrete interval includes both ends, so one whose ends
        # agree is the exact observation P(Y = y) = F(y) - F(y - 1) and
        # is a legitimate way to write a row that happens to be known.
        # On a continuous response that same interval is an event of
        # probability zero, so it stays refused.
        if (identical(resp$family[["type"]], "discrete")) {
          if (any(av[["cens_y2"]][i2] < yv[i2])) {
            frm_stop("Interval upper bounds (y2) must be at least the lower ",
                     "bounds (the response). Both ends of a discrete ",
                     "interval are included, so y2 == y is the exact ",
                     "observation and is allowed", call. = FALSE)
          }
        } else if (any(av[["cens_y2"]][i2] <= yv[i2])) {
          frm_stop("Interval upper bounds (y2) must exceed the lower bounds ",
                   "(the response)", call. = FALSE)
        }
        # NA bounds on non-interval rows are legal and unused; make them
        # harmless for the taped CDF evaluation
        av[["cens_y2"]][!i2 | is.na(av[["cens_y2"]])] <-
          yv[!i2 | is.na(av[["cens_y2"]])]
      }
      # RIGHT and INTERVAL censoring of a count is where the inclusive
      # reading DIVERGES from brms, which emits P(Y > y) for both. Left
      # censoring and continuous censoring agree exactly, so neither
      # says anything here, and neither does truncation alone. A model
      # ported from brms otherwise changes its answer with nothing at
      # the call site to say so: on 200 poisson draws at lambda = 4
      # right censored at 6 the two readings differ by 21.8 to 33.3 log
      # units, median 29.4 over 40 seeds, about 0.66 per censored row.
      # The 20.8 the 0.54.0 NEWS gives for the same setup is below
      # every one of those 40 draws, so it is not the figure to quote.
      # Last of the CENSORING guards, so a bad bound does not spend the
      # one notice a session gets. A refusal raised after this point
      # still spends it: se() on the same response, or anything the
      # fit refuses later, such as a prior on a class the model has
      # not got. Measured both ways.
      if (!is.null(av[["cens"]]) && any(av[["cens"]] %in% c(1, 2)) &&
          identical(resp$family[["type"]], "discrete")) {
        notify_once(
          "cens_discrete_inclusive",
          "cens() on a discrete response reads every bound as ",
          "INCLUSIVE: right censoring at k is P(Y >= k), an interval ",
          "is P(k <= Y <= k2), and a lower edge enters the CDF as ",
          "F(k - 1). brms reads a RIGHT or INTERVAL bound on a count ",
          "as exclusive, P(Y > y), so subtract one from those bounds ",
          "to reproduce a brms fit. Left censoring and continuous ",
          "censoring agree in both packages. The argument is in ",
          "vignette(\"brms-migration\"). This notice is given once ",
          "per session; options(frmtmb.notices = FALSE) turns it off.")
      }
    }
    if (!is.null(av[["se"]])) {
      check_se_declared(resp)
      if (any(av[["se"]] <= 0)) {
        frm_stop("se() values must be positive", call. = FALSE)
      }
    }
    # cbind(successes, failures) reaches here already rewritten to
    # successes + trials(successes + failures); a fractional failure
    # count would otherwise surface as the generic "response must be
    # integer counts in [0, trials]" message, which names neither
    # column. [glmmTMB#1319]
    if (isTRUE(resp$cbind_resp)) {
      fails <- av[["trials"]] - y[[resp$resp_name]]
      if (any(fails < 0) || any(fails != round(fails))) {
        frm_stop("cbind(successes, failures): the failure column must hold ",
                 "non-negative integer counts", call. = FALSE)
      }
    }
    # glm/glmer compatibility: a proportion response with trials()
    # becomes integer counts before validation
    if (resp$family[["family"]] %in% c("binomial", "beta_binomial") &&
        !is.null(av[["trials"]])) {
      yv <- y[[resp$resp_name]]
      if (is.numeric(yv) && !is.matrix(yv) && any(yv != round(yv))) {
        yc <- yv * av[["trials"]]
        if (max(abs(yc - round(yc))) > 1e-6) {
          frm_stop(resp$family[["family"]], ": a proportion response times ",
                   "trials() must give integer counts", call. = FALSE,
                   package = frm_family_package(resp$family))
        }
        y[[resp$resp_name]] <- round(yc)
      }
    }
    # A structured family assembles its own data block here, once, with
    # the response coerced and the random-effect blocks not yet built.
    # `block[["y"]]` is how a family that keeps NA rows substitutes the
    # masked placeholder every later stage sees.
    if (!is.null(st_[["frame_block"]])) {
      blk <- st_[["frame_block"]](resp, spec, av, mf,
                                  y[[resp$resp_name]], n)
      if (!is.null(blk[["y"]])) y[[resp$resp_name]] <- blk[["y"]]
      blocks[[resp$resp_name]] <- blk
    }
    # The reserved names a factorization slot rests on, checked where a
    # refusal can still name the family and the response rather than
    # surfacing as a length mismatch inside the correction.
    check_structure_block(st_, blocks[[resp$resp_name]], resp$family,
                          resp$resp_name, n)
    if (!is.null(resp$autocor)) {
      ac <- check_autocor_response(resp, spec, av, y[[resp$resp_name]])
      ac <- autocor_block(ac, resp$resp_name, mf, resp$formula_env, n)
      ac[["theta_idx"]] <- n_thetaac + seq_len(ac[["npar"]])
      ac[["block_label"]] <- if (length(spec$responses) > 1L) {
        paste0(resp$resp_name, " ", ac[["label"]])
      } else ac[["label"]]
      ac[["student"]] <- identical(resp$family[["family"]], "student")
      n_thetaac <- n_thetaac + ac[["npar"]]
      autocor[[resp$resp_name]] <- ac
    }
    # the prior table skips the response check (check_response = FALSE),
    # except on an ordinal family, whose threshold rows are counted from
    # the response: brms's extract_nthres() refuses there too, "Could
    # not extract the number of thresholds"
    if (!is.null(resp$family[["valid_y"]]) &&
          (check_response || identical(resp$family[["type"]], "ordinal"))) {
      resp$family[["valid_y"]](y[[resp$resp_name]], av)
    }
    # The allow-list, LAST of the addition-term guards and after
    # valid_y, so that a family refusing a term in its own words still
    # gets to say them: this one has only the declaration to go on.
    # What it catches is what nothing else does - a term the density
    # never reads and the core never acts on, parsed, stored on the
    # fit, and silently ignored.
    check_accepted_aterms(resp, av)
    # A family that is not fully determined until the response is in
    # hand - a link bounded above by min(y), a default only the data can
    # supply - gets its one chance here, after the response is coerced
    # and checked and before any link function has run. The order is
    # documented in ?frmtmb_family, because an extension that derives a
    # link this way is betting on it. Whatever comes back becomes the
    # family every later stage reads, so the per-dpar link copies the
    # parse made from the ORIGINAL family are refreshed with it.
    if (!is.null(resp$family[["family_finalize"]])) {
      fam_fin <- resp$family[["family_finalize"]](resp$family,
                                             y[[resp$resp_name]], av)
      if (!inherits(fam_fin, "frmtmb_family")) {
        frm_stop(resp$family[["family"]], ": family_finalize() must return a ",
                 "family object, not ", arg_desc(fam_fin),
                 ". Modify the family it is given and return it",
                 call. = FALSE)
      }
      fam_fin$links <- Map(function(lk, dp) get_link(lk, dpar = dp),
                           fam_fin$links, names(fam_fin$links))
      resp$family <- fam_fin
      for (i_dp in seq_along(resp$dpars)) {
        lk_dp <- fam_fin$links[[resp$dpars[[i_dp]]$name]]
        if (!is.null(lk_dp)) resp$dpars[[i_dp]]$link <- lk_dp
      }
      spec$responses[[resp$resp_name]] <- resp
    }
    # an ordered factor names its categories, and simulate() hands draws
    # back as that factor, so thres(x = ) may not ask for more categories
    # than it has levels; brms returns bare codes there instead
    th_ <- resp$family[["thres"]]
    lv_ <- y_levels[[resp$resp_name]]
    # a hurdle family's category 0 is one level more
    ncat_ <- if (!is.null(th_)) {
      max(th_[["nthres"]]) + 2L - ord_code0(resp$family)
    }
    if (!is.null(th_) && !is.null(lv_) && ncat_ > length(lv_)) {
      frm_stop("thres(x = ", max(th_[["nthres"]]), ") asks for ",
               ncat_, " categories, and the response ",
               "is an ordered factor with ", length(lv_), " levels in the ",
               "data: a level no row takes is dropped with the model ",
               "frame, as brms drops it. frmtmb returns simulated ",
               "responses as that factor, so it cannot hold more ",
               "categories than levels. Code the response as integers ",
               ord_code0(resp$family), "..", max(th_[["nthres"]]) + 1L,
               " to fit the categories nobody chose, or lower the count",
               call. = FALSE)
    }
    # Family-level DATA a likelihood needs but no addition term supplies
    # (the Cox baseline's spline bases). It is a function of the
    # validated response, so it is built here, once, and rides with the
    # addition-term values the objective already bakes into the tape.
    if (!is.null(resp$family[["aterm_data"]])) {
      av <- c(av, resp$family[["aterm_data"]](y[[resp$resp_name]], av))
    }
    if (!is.null(resp$family[["extra_pars"]])) {
      ex_r <- resp$family[["extra_pars"]](y[[resp$resp_name]], av)
      if (length(spec$responses) > 1L) {
        # Every response's extras live in one parameter list, so each
        # response's block is namespaced by the response and handed
        # back to its own density under the family's own names.
        check_mv_extra_family(resp$family)
        tpl_nm <- mv_extra_name(resp$resp_name, names(ex_r))
        extra_map[[resp$resp_name]] <- stats::setNames(tpl_nm, names(ex_r))
        names(ex_r) <- tpl_nm
      }
      extras <- c(extras, ex_r)
    }
    aterm_values[[resp$resp_name]] <- av
  }
  mf <- mf_all
  n <- n_all

  # me() latent values: after every mi() response, so the mi() slots
  # of `miss` keep the positions they have always had
  me_fr <- me_build_frame(spec, mf, env, n, n_miss)
  if (!is.null(me_fr)) {
    n_miss <- n_miss + me_fr$n_latent
    miss_init <- c(miss_init, me_fr$latent_init)
  }

  ## Phase 1: per-linpred design matrices and random-effect components.
  linpreds <- list()
  components <- list()
  beta_names <- character(0)
  betad_names <- character(0)
  betad_fixed_idx <- integer(0)

  for (resp in spec$responses) {
    # a subset() response's designs are built on its own rows
    mf <- frame_rows(mf_all, sub_rows[[resp$resp_name]])
    n <- nrow(mf)
    for (dp in resp$dpars) {
      lp_key <- linpred_key(resp$resp_name, dp[["name"]])
      is_primary <- dp[["name"]] %in% resp$primary_dpars
      par_name <- if (is_primary) "beta" else "betad"
      dp_prefix <- if (identical(dp[["name"]], "mu")) "" else
        paste0(dp[["name"]], ": ")
      if (length(spec$responses) > 1) {
        dp_prefix <- paste0(resp$resp_name, " ", dp_prefix)
      }
      if (!is.null(dp[["equate"]])) {
        # filled in after this response's other dpars (below), keeping
        # this dpar's place in the order of the predictors
        linpreds[[lp_key]] <- list()
        next
      }

      if (!is.null(dp[["nl_body"]])) {
        # a nonlinear dpar has no design of its own: it is evaluated
        # from the nonlinear-parameter values and the raw data columns
        # inside the objective, after the parameters it names
        data_list <- lapply(stats::setNames(dp[["datavars"]], dp[["datavars"]]),
                            function(v) {
                              val <- mf[[v]]
                              if (is.null(val)) {
                                frm_stop("Variable '", v, "' from the ",
                                         "nonlinear formula not found",
                                         call. = FALSE)
                              }
                              val
                            })
        # ps() terms: the frozen knots and the eigensplit of the
        # second-difference penalty. The null space joins the fixed
        # coefficients here; the range space becomes a component below,
        # deliberately WITHOUT a comp_id, because a nonlinear body reads
        # the block through a closure rather than through `Z b` and
        # there is no Z column for it to occupy.
        ps_terms <- dp[["ps_terms"]] %||% list()
        if (length(ps_terms) && length(spec$responses) > 1L) {
          frm_stop("ps() is not supported in a multivariate model yet: its ",
                   "block reaches the body through the per-call evaluation ",
                   "frame, and which response's frame that is has not been ",
                   "measured", call. = FALSE)
        }
        for (ti in seq_along(ps_terms)) {
          pt <- ps_build(ps_terms[[ti]], mf, resp$nlpars %||% character(0),
                         resp$formula_env)
          cn <- paste0(dp_prefix, pt[["label"]], ".fx",
                       seq_len(pt[["n_fixed"]]))
          # The null space of the penalty is a FIXED coefficient of this
          # predictor, so it belongs with the location coefficients. An
          # nl `mu` is not "primary" (its own design is empty and core
          # sends it to `betad`), and following that rule here would
          # file a curve's unpenalized directions among the
          # distributional parameters.
          pt[["par"]] <- if (is_primary ||
              dp[["name"]] %in% (resp$family[["primary_dpars"]] %||% "mu")) {
            "beta"
          } else {
            "betad"
          }
          if (pt[["par"]] == "beta") {
            pt[["beta_idx"]] <- length(beta_names) + seq_len(pt[["n_fixed"]])
            beta_names <- c(beta_names, cn)
          } else {
            pt[["beta_idx"]] <- length(betad_names) + seq_len(pt[["n_fixed"]])
            betad_names <- c(betad_names, cn)
          }
          components[[length(components) + 1L]] <- list(
            lp_key = lp_key, dpar = dp[["name"]], resp = resp$resp_name,
            covstruct = "smooth", id = NULL,
            dim = pt[["n_pen"]], n_levels = 1L,
            levels = NULL,
            cnms = paste0(pt[["label"]], ".", seq_len(pt[["n_pen"]])),
            bar = NULL,
            Zlocal = methods::as(Matrix::Matrix(0, n, pt[["n_pen"]],
                                                sparse = TRUE),
                                 "CsparseMatrix"),
            group_name = pt[["label"]],
            label = paste0(dp_prefix, pt[["label"]])
          )
          pt[["comp_id"]] <- length(components)
          ps_terms[[ti]] <- pt
        }
        linpreds[[lp_key]] <- list(
          resp = resp$resp_name, dpar = dp[["name"]], X = NULL,
          n_param_cols = 0L, Z = NULL, par = par_name, idx = integer(0),
          offset = NULL, link = dp[["link"]], terms = NULL, xlevels = NULL,
          contrasts = NULL, smooths = list(), comp_ids = integer(0),
          constant = NULL, nl_body = dp[["nl_body"]], data_list = data_list,
          nl_env = dp[["nl_env"]], nl_lexical = dp[["nl_lexical"]],
          nl_pars = dp[["nl_pars"]] %||% character(0),
          nl_dpar_refs = dp[["nl_dpar_refs"]] %||% character(0),
          ps_terms = ps_terms
        )
        next
      }

      tt <- stats::terms(dp[["fixed"]])
      # bf(cmc = FALSE): R codes a factor by its cell means when the
      # formula has no intercept, and brms's cmc = FALSE keeps the
      # treatment contrasts and removes only the intercept column
      # (brms:::validate_terms()). The terms keep the intercept, so a
      # prediction builds the same columns and selects param_colnames.
      cmc_drop <- isFALSE(dp[["cmc"]]) && attr(tt, "intercept") == 0L
      if (cmc_drop) attr(tt, "intercept") <- 1L
      X <- if (sparse_x) sparse_mm(tt, mf) else stats::model.matrix(tt, mf)
      contr <- attr(X, "contrasts")   # subsetting X below drops the attr
      if (cmc_drop) X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      if (isTRUE(resp$family[["drop_intercept"]]) && is_primary) {
        # thresholds replace the intercept; a threshold-only model
        # (y ~ 1) leaves X with zero columns, which is fine
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      }
      if (!is.null(dp[["rsv_intercept"]])) {
        X <- rsv_intercept_order(X, dp[["rsv_intercept"]])
      }
      # rank-deficient designs: drop aliased columns like lm() (lme4#144).
      # Sparse X densifies a copy only when the cheap screen flags a
      # possible deficiency, so the dropped-column set never differs
      # from the dense path.
      alias_null <- NULL
      dropped_colnames <- NULL
      if (ncol(X) > 1L) {
        Xq <- if (!sparse_x) X
              else if (sparse_maybe_deficient(X)) as.matrix(X)
              else NULL
        if (!is.null(Xq)) {
          Xq <- as.matrix(Xq)
          qrX <- qr(Xq)
          if (qrX$rank < ncol(X)) {
            dropped <- colnames(X)[qrX$pivot[(qrX$rank + 1L):ncol(X)]]
            frm_message("Fixed-effect design of '", lp_key,
                        "' is rank deficient; dropping column(s): ",
                        paste(dropped, collapse = ", "))
            # Directions the data could not identify: null(X) is the
            # orthogonal complement of the row space, so the trailing
            # columns of the complete Q of t(X) span it. Frozen here so
            # prediction can refuse rows that load on them. [lme4#303]
            alias_null <- qr.Q(qr(t(Xq)), complete = TRUE)[
              , (qrX$rank + 1L):ncol(Xq), drop = FALSE]
            rownames(alias_null) <- colnames(X)
            dropped_colnames <- dropped
            X <- X[, setdiff(colnames(X), dropped), drop = FALSE]
          }
        }
      }
      param_colnames <- colnames(X)
      n_param_cols <- ncol(X)
      off <- extract_offset(tt, mf, resp$formula_env)
      xlev <- stats::.getXlevels(stats::terms(mf), mf)
      comp_ids <- integer(0)

      if (length(dp[["re"]])) {
        bars <- lapply(dp[["re"]], `[[`, "bar")
        # mm() bars never reach mkReTrms: their grouping expression is a
        # call reformulas cannot evaluate, and their design row loads
        # several levels at once. They are built below instead, into a
        # component of exactly the same shape, so everything downstream
        # (blocks, theta, ranef, VarCorr, simulate) is unchanged.
        is_mm <- vapply(dp[["re"]], function(z) !is.null(z$mm), TRUE)
        plain_k <- which(!is_mm)
        rt_pos <- integer(length(bars))
        rt_pos[plain_k] <- seq_along(plain_k)
        rt <- NULL
        fassign <- integer(0)
        cmc_re <- logical(length(bars))
        if (length(plain_k)) {
          # bars keeps the user's expressions (labels, prediction); only
          # the copy handed to reformulas is name-resolved
          rt_bars <- bars[plain_k]
          # cmc = FALSE reaches the group-level terms too, as in brms:
          # the term gets an intercept here, which goes again below
          if (isFALSE(dp[["cmc"]])) {
            for (j in seq_along(rt_bars)) {
              wi <- bar_with_intercept(rt_bars[[j]])
              # a term whose columns cmc = FALSE leaves as they are
              # (numeric slopes, f:x without f) needs no rewrite, and
              # is no reason to refuse a structure
              if (!is.null(wi) &&
                    !cmc_changes_columns(rt_bars[[j]], mf,
                                         resp$formula_env)) {
                wi <- NULL
              }
              cs_j <- dp[["re"]][[plain_k[j]]]$covstruct
              # these read one coefficient per level of the factor
              # (a time, a position), which contrasts are not, treatment
              # and polynomial alike; brms has none of them, so there is
              # no brms reading
              if (!is.null(wi) &&
                    cs_j %in% c("ar1", "hetar1", "cs", "homcs", "toep",
                                "homtoep", "ou", "exp", "gau", "mat")) {
                cc <- cmc_bar_columns(rt_bars[[j]], mf, resp$formula_env)
                frm_stop("cmc = FALSE does not apply to ", cs_j, "(",
                         deparse1(bars[[plain_k[j]]]), "): the structure ",
                         "reads one coefficient per level",
                         if (!is.null(cc)) {
                           paste0(" (", cmc_col_list(cc$cell), ")")
                         },
                         ", and cmc = FALSE would code the term by its ",
                         "contrasts",
                         if (!is.null(cc)) {
                           paste0(" (", cmc_col_list(cc$contrast), ")")
                         },
                         ", which are not coefficients of single levels. ",
                         "Drop cmc = FALSE from this formula", call. = FALSE)
              }
              if (!is.null(wi)) {
                rt_bars[[j]] <- wi
                cmc_re[plain_k[j]] <- TRUE
              }
            }
          }
          grp <- resolve_group_calls(rt_bars, mf, resp$formula_env)
          rt <- reformulas::mkReTrms(grp$bars, fr = grp$fr,
                                     reorder.terms = FALSE)
          fassign <- attr(rt$flist, "assign")
        }
        dup_cps <- list()
        for (k in seq_along(bars)) {
          cs_name <- dp[["re"]][[k]]$covstruct
          if (is_mm[k]) {
            mms <- dp[["re"]][[k]]$mm
            # cmc = FALSE rides on the spec, which the component keeps,
            # so a prediction builds the members' designs the same way
            if (isFALSE(dp[["cmc"]])) mms[["cmc"]] <- FALSE
            gvals <- mm_member_values(mms, mf, resp$formula_env)
            levs <- mm_pooled_levels(gvals)
            iw <- mm_index_weights(mms, mf, resp$formula_env, levs)
            md <- mm_member_designs(mms, mf, resp$formula_env,
                                    iw$n_members)
            Zk <- mm_local_Z(iw$J, iw$W, md$designs, length(levs))
            components[[length(components) + 1L]] <- list(
              lp_key = lp_key, dpar = dp[["name"]], resp = resp$resp_name,
              covstruct = cs_name, id = NULL, rank = NULL,
              dim = length(md$cnms), n_levels = length(levs),
              levels = levs, cnms = md$cnms,
              bar = bars[[k]], Zlocal = methods::as(Zk, "CsparseMatrix"),
              mm = mms,
              written = dp[["re"]][[k]]$written,
              group_name = mms$label,
              label = paste0(dp_prefix, deparse1(bars[[k]]))
            )
            comp_ids <- c(comp_ids, length(components))
            dup_cps[[length(dup_cps) + 1L]] <- components[[length(components)]]
            if (!is.null(mms$by)) {
              # mm(g1, g2, by = cbind(f1, f2)): brms maps each POOLED
              # level to one by-level, read across every member column
              cp <- components[[length(components)]]
              byv <- by_eval(mms$by, mf, resp$formula_env, n, "in the data",
                             members = iw$n_members)
              lev_by <- by_level_map(iw$J, byv, length(levs), mms$gvars,
                                     mms$by)
              subs <- by_split_component(cp, mms$by, lev_by,
                                         by_extract_levels(byv))
              components[[length(components)]] <- subs[[1L]]
              for (s in subs[-1L]) {
                components[[length(components) + 1L]] <- s
                comp_ids <- c(comp_ids, length(components))
              }
            }
            next
          }
          kk <- rt_pos[k]
          zrows <- rt$Gp[kk] + seq_len(rt$Gp[kk + 1L] - rt$Gp[kk])
          if (cmc_re[k]) {
            # the intercept added for cmc = FALSE goes, and the treatment
            # contrasts stay; Zt is level-major within a term
            keep <- rt$cnms[[kk]] != "(Intercept)"
            zrows <- zrows[rep(keep, length(zrows) %/% length(keep))]
            rt$cnms[[kk]] <- rt$cnms[[kk]][keep]
          }
          d_k <- length(rt$cnms[[kk]])
          len_k <- length(zrows)
          dist_cs <- c("ou", "exp", "gau", "mat")
          if (cs_name %in% c("ar1", "hetar1", "cs", "homcs", "toep",
                             "homtoep", "rr", dist_cs) && d_k < 2L) {
            frm_stop(cs_name, "() needs at least 2 terms per level",
                     call. = FALSE)
          }
          if (cs_name %in% c("ar1", "hetar1", dist_cs) &&
              "(Intercept)" %in% rt$cnms[[kk]]) {
            frm_stop(cs_name, "() requires a factor without intercept on ",
                     "the left of the bar, e.g. ", cs_name,
                     "(times + 0 | g)", call. = FALSE)
          }
          if (cs_name %in% c("ar1", "hetar1")) {
            warn_ar1_level_gaps(bars[[k]], mf, cs_name)
          }
          fac <- rt$flist[[fassign[kk]]]
          aux_A <- NULL
          aux_D <- NULL
          aux_kron <- NULL
          if (cs_name %in% dist_cs) {
            v <- all.vars(bars[[k]][[2]])
            if (length(v) != 1L || !is.factor(mf[[v]])) {
              frm_stop(cs_name, "() needs a single factor built with ",
                       "num_factor(): ", cs_name, "(pos + 0 | g)",
                       call. = FALSE)
            }
            coords <- parse_num_levels(levels(mf[[v]]))
            if (is.matrix(coords)) {
              aux_D <- as.matrix(stats::dist(coords))
            } else {
              aux_D <- abs(outer(coords, coords, "-"))
            }
          }
          aux_Q <- NULL
          aux_Qk <- NULL
          if (cs_name == "gr_prec") {
            Q <- lookup_structural(dp[["re"]][[k]]$cov_expr, data2, data,
                                   resp$formula_env, "gr(prec = )")
            if (is.null(rownames(Q)) ||
                !all(levels(fac) %in% rownames(Q))) {
              frm_stop("gr(prec=): prec needs dimnames covering all ",
                       "grouping levels", call. = FALSE)
            }
            lv <- levels(fac)
            aux_Q <- methods::as(Matrix::Matrix(Q[lv, lv], sparse = TRUE),
                                 "generalMatrix")
            if (d_k > 1L) aux_Qk <- kron_prec_parts(aux_Q, d_k)
          }
          if (cs_name == "gr_cov") {
            A <- lookup_structural(dp[["re"]][[k]]$cov_expr, data2, data,
                                   resp$formula_env, "gr(cov = )")
            if (!is.matrix(A) || nrow(A) != ncol(A)) {
              frm_stop("gr(cov=): cov must be a square matrix", call. = FALSE)
            }
            lv <- levels(fac)
            if (is.null(rownames(A)) || !all(lv %in% rownames(A))) {
              frm_stop("gr(cov=): cov needs dimnames covering all grouping ",
                       "levels", call. = FALSE)
            }
            aux_A <- unname(A[lv, lv])
            if (d_k > 1L) {
              aux_kron <- kron_cov_index(d_k, len_k %/% d_k)
            }
          }
          if (cs_name == "equalto") {
            V <- lookup_structural(dp[["re"]][[k]]$cov_expr, data2, data,
                                   resp$formula_env, "equalto()")
            if (!is.matrix(V) || nrow(V) != d_k || ncol(V) != d_k) {
              frm_stop("equalto(): V must be a ", d_k, " x ", d_k,
                       " matrix", call. = FALSE)
            }
            aux_A <- unname(V)
          }
          Zk <- Matrix::t(rt$Zt[zrows, , drop = FALSE])
          components[[length(components) + 1L]] <- list(
            lp_key = lp_key, dpar = dp[["name"]], resp = resp$resp_name,
            covstruct = cs_name, id = dp[["re"]][[k]]$id,
            rank = dp[["re"]][[k]]$rank,
            # gr(dist = "student"): the FIXED degrees of freedom of the
            # t latent, NULL on every gaussian block
            dist_nu = dp[["re"]][[k]]$dist_nu,
            written = dp[["re"]][[k]]$written,
            from_slash = dp[["re"]][[k]]$from_slash,
            dim = d_k, n_levels = len_k %/% d_k,
            levels = levels(fac), cnms = rt$cnms[[kk]],
            bar = bars[[k]], Zlocal = Zk, aux_A = aux_A,
            aux_D = aux_D, aux_kron = aux_kron, aux_Q = aux_Q,
            aux_Qk = aux_Qk,
            group_name = names(rt$flist)[fassign[kk]],
            label = paste0(dp_prefix, deparse1(bars[[k]]))
          )
          if (cmc_re[k]) {
            # cmc = FALSE: a prediction rebuilds the term with the
            # intercept and drops it, as the fit did
            components[[length(components)]][["cmc_intercept"]] <- TRUE
          }
          comp_ids <- c(comp_ids, length(components))
          # gr(g, by = f): the term as written is what the duplicate
          # check reads, and one block per by-level is what is fitted
          dup_cps[[length(dup_cps) + 1L]] <- components[[length(components)]]
          by_k <- dp[["re"]][[k]]$by
          if (!is.null(by_k)) {
            cp <- components[[length(components)]]
            byv <- by_eval(by_k, mf, resp$formula_env, n, "in the data")
            lev_by <- by_level_map(as.integer(fac), byv, length(cp$levels),
                                   cp$group_name, by_k)
            subs <- by_split_component(cp, by_k, lev_by,
                                       by_extract_levels(byv))
            components[[length(components)]] <- subs[[1L]]
            for (s in subs[-1L]) {
              components[[length(components) + 1L]] <- s
              comp_ids <- c(comp_ids, length(components))
            }
          }
        }
        refuse_duplicated_re(dup_cps)
      }

      # Smooths: fixed (null-space) part into X, wiggly part as an
      # iid-Gaussian component whose variance is the inverse smoothing
      # parameter.
      sm_info <- list()
      for (sspec in dp[["smooth"]] %||% list()) {
        # modCon = 3 ("set fit and predict constraint to fit constraint").
        # A t2 smooth carries two identifiability constraints: `C`, absorbed
        # into `X`/`S` for fitting, and `Cp`, absorbed into `Xp`/`Sp` and
        # honored by PredictMat. Left at the default the two differ, so
        # PredictMat does not rebuild the fitted basis and newdata
        # prediction is impossible. gamm4 keeps both and maps the fitted
        # coefficients into the predict parameterization afterwards (its
        # `G$P`, "important for t2 smooths, where fit constraint is not
        # good for component wise prediction s.e.s"); we have no component
        # wise s.e. to protect, so dropping `Cp` is simpler. This does not
        # change the fit: for t2 `X` and `S` are bit-identical either way
        # (`modCon >= 3` only sets `sm$Cp <- NULL`), and for s() `Cp` is
        # already NULL so the argument is inert.
        #
        # diagonal.penalty = TRUE is brms's (data_sm(), frame_basis_sm()
        # since brms 2.8.7). It reparameterizes a single-penalty smooth
        # so the penalty is the identity on its range, which fixes the
        # scale of the null-space column: without it frmtmb's `sx_1` was
        # brms's coefficient times a data-dependent factor (5.0 to 7.7
        # on gamSim data, of either sign), so a prior or a hypothesis on
        # it meant another parameter. The random part, and so the fit,
        # does not change. mgcv ignores the flag for a smooth with
        # several penalties (t2()), as it does in brms's call.
        scl <- mgcv::smoothCon(sspec, data = mf, absorb.cons = TRUE,
                               modCon = 3, diagonal.penalty = TRUE)
        for (sm in scl) {
          re2 <- mgcv::smooth2random(sm, names(mf), type = 2)
          if (isTRUE(re2$fixed)) {
            frm_stop("Fixed (fx = TRUE) smooths are not supported: ", sm$label,
                     call. = FALSE)
          }
          nr <- integer(0)
          sm_comp_ids <- integer(0)
          for (r in seq_along(re2$rand)) {
            Xr <- re2$rand[[r]]
            k_r <- ncol(Xr)
            components[[length(components) + 1L]] <- list(
              lp_key = lp_key, dpar = dp[["name"]], resp = resp$resp_name,
              covstruct = "smooth", id = NULL,
              dim = k_r, n_levels = 1L,
              levels = NULL, cnms = paste0(sm$label, ".", seq_len(k_r)),
              bar = NULL, Zlocal = methods::as(Xr, "CsparseMatrix"),
              group_name = sm$label,
              label = paste0(dp_prefix, sm$label)
            )
            sm_comp_ids <- c(sm_comp_ids, length(components))
            nr <- c(nr, k_r)
          }
          comp_ids <- c(comp_ids, sm_comp_ids)
          nf <- if (is.null(re2$Xf)) 0L else ncol(re2$Xf)
          xf_idx <- integer(0)
          if (nf > 0L) {
            Xf <- re2$Xf
            colnames(Xf) <- paste0(sm$label, ".fx", seq_len(nf))
            xf_idx <- ncol(X) + seq_len(nf)
            X <- cbind(X, Xf)
          }
          sm_info[[length(sm_info) + 1L]] <- list(
            sm = sm, U = re2$trans.U, D = re2$trans.D,
            ord = smooth_pen_order(sm, re2),
            nr = nr, nf = nf, xf_idx = xf_idx,
            comp_ids = sm_comp_ids, block_ids = NULL,
            # the grouping factor this basis is indexed by, if any: what
            # separates a per-level curve (which re_formula = NA drops)
            # from a population smooth (which it keeps). Read off the
            # smooth object here, where the model frame is still around
            # to say which of its terms are factors.
            group_var = smooth_group_var(sm, mf),
            # the fitted levels of that factor: fs smooths carry them as
            # sm$flev, but a factor bs = "re" smooth does not, and an
            # unseen level would otherwise die inside PredictMat with a
            # non-conformable-arguments error instead of the named
            # new-levels refusal
            group_levels = local({
              gv <- smooth_group_var(sm, mf)
              if (is.null(gv) || is.null(mf[[gv]])) NULL else
                levels(as.factor(mf[[gv]]))
            }),
            label = sm$label,
            # the written term this basis came from: a factor `by`
            # splits one term into a basis per level, and
            # conditional_smooths() draws the term, not the level
            term = attr(sspec, "frm_term") %||% sm$label
          )
        }
      }

      # Gaussian-process terms: gp(x1, ...) is exact (a dense SE-kernel
      # block over the unique coordinate rows); gp(..., k=) is the
      # Hilbert-space approximation (tensor-product sine basis in Z,
      # spectral-density prior SDs). D up to 3; iso shares one
      # lengthscale, brms's default; iso = FALSE gives one per
      # dimension. A factor `by` is one independent GP per contrast
      # column (one per level under brms's cmc = TRUE), each over its own
      # rows with its own scaling and its own sd and lengthscales, as
      # brms:::data_gp() builds them; a numeric `by` multiplies one GP.
      gp_info <- list()
      for (ge in dp[["gpterms"]] %||% list()) {
        Xc <- do.call(cbind, lapply(ge$exprs, function(ex) {
          as.numeric(eval(ex, mf, resp$formula_env))
        }))
        Dg <- ncol(Xc)
        iso <- isTRUE(ge$iso) || Dg == 1L
        vnames <- vapply(ge$exprs, deparse1, "")
        lab0 <- gp_term_label(ge, vnames)
        cvec <- ge$c
        if (!is.null(ge$k)) {
          if (length(cvec) == 1L) cvec <- rep(cvec, Dg)
          if (length(cvec) != Dg) {
            frm_stop("gp(): c = must be length 1 or the number of ",
                     "variables (", Dg, ")", call. = FALSE)
          }
        }
        byv <- if (!is.null(ge$by)) {
          eval(ge$by, mf, resp$formula_env)
        }
        subs <- gp_sub_terms(ge, byv, vnames, iso, nrow(Xc))
        for (sb in subs) {
          rows <- sb$rows
          Xs <- Xc[rows, , drop = FALSE]
          lab <- paste0(lab0, sb$lab_sfx)
          # brms's gr = TRUE keeps one latent value per distinct
          # position; gr = FALSE one per observation
          uq <- if (isTRUE(ge$gr)) {
            Xs[!duplicated(pos_rowkey(Xs)), , drop = FALSE]
          } else {
            Xs
          }
          # brms's scale = TRUE divides by the largest pairwise distance
          # over those rows (brms:::.data_gp), level by level
          dmax <- if (isTRUE(ge$scale)) gp_max_dist(uq) else 1
          if (isTRUE(ge$scale) && !isTRUE(dmax > 0)) {
            frm_stop("gp(", paste(vnames, collapse = ", "), "): the ",
                     "coordinates have no spread",
                     if (nzchar(sb$lab_sfx)) paste0(" in", sb$lab_sfx),
                     ", so they cannot be scaled. brms says: Could not ",
                     "scale GP covariates. Please set 'scale' to FALSE ",
                     "in 'gp'", call. = FALSE)
          }
          # which written term this sub-GP belongs to: brms names and
          # orders a term's parameters together, levels innermost
          brms_meta <- list(sfx1 = sb$sfx1, sfx2 = sb$sfx2,
                            term = paste0(lp_key, "\r", lab0))
          if (is.null(ge$k)) {
            if (isTRUE(ge$gr)) {
              posdf <- as.data.frame(uq)
              posdf <- posdf[do.call(order, posdf), , drop = FALSE]
              pos <- unname(as.matrix(posdf))
              jpos <- match(pos_rowkey(Xs), pos_rowkey(pos))
            } else {
              pos <- unname(Xs)
              jpos <- seq_len(nrow(Xs))
            }
            npos <- nrow(pos)
            if (npos > 500L) {
              frm_stop("gp() without k= builds a dense ", npos,
                       "-point covariance; use k= for the Hilbert-space ",
                       "approximation", call. = FALSE)
            }
            Zg <- Matrix::sparseMatrix(i = rows, j = jpos, x = sb$mult,
                                       dims = c(nrow(Xc), npos))
            components[[length(components) + 1L]] <- list(
              lp_key = lp_key, dpar = dp[["name"]], resp = resp$resp_name,
              covstruct = "gp", id = NULL,
              dim = npos, n_levels = 1L,
              levels = NULL, cnms = paste0(lab, ".", seq_len(npos)),
              bar = NULL, Zlocal = methods::as(Zg, "CsparseMatrix"),
              aux_D2 = lapply(seq_len(Dg), function(j) {
                outer(pos[, j], pos[, j], "-")^2
              }),
              gp_D = Dg, gp_iso = iso, gp_vars = vnames,
              # the exact form keeps data units inside; brms reports
              # its lengthscale on the scaled inputs
              gp_lscale_div = dmax, gp_brms = brms_meta,
              group_name = lab,
              label = paste0(dp_prefix, lab)
            )
            gp_info[[length(gp_info) + 1L]] <- list(
              exprs = ge$exprs, type = "exact", positions = pos,
              prior_x = uq / dmax,
              by = sb$by, comp_id = length(components), block_id = NULL,
              label = lab
            )
          } else {
            # brms's input convention (brms:::.data_gp): rescale by the
            # largest pairwise distance, center on that scale, then take
            # a shared boundary L_j = c_j * max(1, range of the whole
            # centered matrix). The same gp(x, k, c) call is then the
            # same approximation here and in brms. The rows are the
            # distinct positions under gr = TRUE, because brms collapses
            # duplicates before computing the scale, so ties would
            # otherwise shift the center.
            ctr <- colMeans(uq / dmax)
            Lb <- gp_choose_L(sweep(uq / dmax, 2, ctr), cvec)
            m <- ge$k
            if (m^Dg > 1000) {
              frm_stop("gp(): k = ", m, " over ", Dg, " dimensions gives ",
                       m^Dg, " basis columns (cap 1000); lower k=",
                       call. = FALSE)
            }
            idx <- as.matrix(do.call(expand.grid,
                                     rep(list(seq_len(m)), Dg)))
            omega <- sweep(idx * pi, 2, 2 * Lb, "/")
            Phi <- matrix(0, nrow(Xc), nrow(omega))
            Phi[rows, ] <- sb$mult *
              hsgp_basis(sweep(Xs / dmax, 2, ctr), omega, Lb)
            M_b <- nrow(omega)
            components[[length(components) + 1L]] <- list(
              lp_key = lp_key, dpar = dp[["name"]], resp = resp$resp_name,
              covstruct = "hsgp", id = NULL,
              dim = M_b, n_levels = 1L,
              levels = NULL, cnms = paste0(lab, ".", seq_len(M_b)),
              bar = NULL, Zlocal = methods::as(Phi, "CsparseMatrix"),
              aux_omega = omega,
              gp_D = Dg, gp_iso = iso, gp_vars = vnames,
              gp_dmax = dmax, gp_lscale_div = 1, gp_brms = brms_meta,
              group_name = lab,
              label = paste0(dp_prefix, lab)
            )
            gp_info[[length(gp_info) + 1L]] <- list(
              exprs = ge$exprs, type = "hsgp", center = ctr, L = Lb,
              dmax = dmax, omega = omega, by = sb$by, prior_x = uq / dmax,
              comp_id = length(components), block_id = NULL, label = lab
            )
          }
          comp_ids <- c(comp_ids, length(components))
        }
      }

      # Spatial GMRF terms: car(M, gr = g) over an adjacency matrix and
      # spde(fem, gr = node) over a mesh. Both are one intercept per
      # location, so the block looks like (1 | g) with a structured
      # precision; a synthetic bar keeps the prediction, ranef() and
      # VarCorr() paths unchanged.
      for (ce in c(dp[["carterms"]] %||% list(),
           dp[["spdeterms"]] %||% list())) {
        is_car <- !is.null(ce$M_expr)
        fn <- if (is_car) "car" else "spde"
        gv <- eval(ce$gr_expr, mf, resp$formula_env)
        if (anyNA(gv)) {
          frm_stop(fn, "(): the grouping variable '", deparse1(ce$gr_expr),
                   "' has missing values", call. = FALSE)
        }
        aux_car <- NULL
        aux_spde <- NULL
        if (is_car) {
          # the adjacency matrix names its locations, so the block's
          # level order is the data's and the matrix is permuted to it
          locs <- if (is.factor(gv)) levels(droplevels(gv)) else {
            sort(unique(as.character(gv)))
          }
          j_loc <- match(as.character(gv), locs)
          M <- lookup_structural(ce$M_expr, data2, data,
                                 resp$formula_env, "car()")
          W <- car_adjacency(M, locs)
          aux_car <- car_aux(W, ce$type, ce$con_sd)
        } else {
          # the mesh names nothing, so the block's levels ARE the mesh
          # rows and the data is permuted to them: every node gets a
          # column, observed or not (an unobserved one keeps its prior)
          fem <- lookup_structural(ce$fem_expr, data2, data,
                                   resp$formula_env, "spde()")
          aux_spde <- spde_matrices(fem)
          n_node <- nrow(aux_spde[["M0"]])
          j_loc <- spde_node_index(gv, n_node, ce$gr_expr)
          locs <- as.character(seq_len(n_node))
        }
        Zc <- Matrix::sparseMatrix(i = seq_along(j_loc), j = j_loc, x = 1,
                                   dims = c(length(j_loc), length(locs)))
        components[[length(components) + 1L]] <- list(
          lp_key = lp_key, dpar = dp[["name"]], resp = resp$resp_name,
          covstruct = fn, id = NULL,
          dim = 1L, n_levels = length(locs),
          levels = locs, cnms = "(Intercept)",
          bar = call("|", 1, ce$gr_expr),
          Zlocal = methods::as(Zc, "CsparseMatrix"),
          aux_car = aux_car, aux_spde = aux_spde,
          car_type = ce$type,
          group_name = deparse1(ce$gr_expr),
          label = paste0(dp_prefix, ce$label)
        )
        comp_ids <- c(comp_ids, length(components))
      }

      # Monotonic terms: one scale coefficient in beta (a zero column in
      # the stored X keeps the bookkeeping - names, idx, vcov - while
      # the objective and the numeric prediction paths supply the
      # simplex-weighted values); the simplex parameters join `extras`.
      mo_info <- list()
      # brms gives every special TERM its own simplex and enumerates the
      # terms in terms() order, which puts every main effect ahead of
      # every interaction. The parser emits them in written order and
      # puts mo(x):z ahead of the mo(x) that mo(x) * z implies, so the
      # list is sorted here: the j-th entry below has to be brms's j-th
      # special term, or no prior or parameter map could line the two
      # up. It is the POSITION that carries that, not the zeta number,
      # which starts past whatever the family already put in `extras`.
      for (ent in mo_terms_in_brms_order(dp[["mo"]] %||% list())) {
        mexpr <- ent$expr
        v <- eval(mexpr, mf, resp$formula_env)
        if (is.factor(v)) {
          if (!is.ordered(v)) {
            frm_stop("mo(): factor variables must be ordered factors",
                     call. = FALSE)
          }
          codes <- as.integer(v) - 1L
          D_mo <- nlevels(v) - 1L
          mo_levels <- levels(v)
        } else {
          if (any(v < 0) || any(v != round(v))) {
            frm_stop("mo(): variable must be an ordered factor or ",
                     "non-negative integers", call. = FALSE)
          }
          codes <- as.integer(v)
          D_mo <- max(codes)
          mo_levels <- NULL
        }
        if (D_mo < 2L) {
          frm_stop("mo() needs at least 3 ordered categories", call. = FALSE)
        }
        vkey <- deparse1(mexpr)
        # one simplex per term occurrence, never shared between a main
        # effect and an interaction on the same variable: the two terms
        # describe different shapes and brms fits them separately
        zname <- paste0("zeta", length(extras) + 1L)
        extras[[zname]] <- numeric(D_mo - 1L)
        mult <- NULL
        if (!is.null(ent$mult)) {
          mult <- check_special_mult(eval(ent$mult, mf, resp$formula_env),
                                     ent$mult, "mo")
        }
        lab <- paste0("mo", vkey,
                      if (!is.null(ent$mult)) {
                        paste0(":", deparse1(ent$mult))
                      } else "")
        X <- cbind(X, matrix(0, nrow(X), 1,
                             dimnames = list(NULL, lab)))
        mo_info[[length(mo_info) + 1L]] <- list(
          expr = mexpr, codes = codes, D = D_mo, levels = mo_levels,
          zeta = zname, col = ncol(X), label = lab,
          mult = mult, mult_expr = ent$mult
        )
      }

      # mi(x) predictor terms: one coefficient (zero placeholder column,
      # as with mo); the values are observed-or-latent, supplied by the
      # objective and the numeric prediction paths
      mi_info <- list()
      # brms enumerates mi() terms in terms() order too, so `mi(x) * z`
      # gives the main effect a column before its interaction and
      # `variables()` reads bsp_<resp>_mi<x> before bsp_<resp>_mi<x>:z,
      # which is brms's order (punch round 1, minor 4)
      for (ent in mo_terms_in_brms_order(dp[["miterms"]] %||% list())) {
        vn <- deparse1(ent$expr)
        tgt <- spec$responses[[vn]]
        if (is.null(tgt) || !isTRUE(tgt$aterms[["mi"]])) {
          frm_stop("mi(", vn, ") needs a matching imputation model: ",
                   "add bf(", vn, " | mi() ~ ...)", call. = FALSE)
        }
        if (identical(vn, resp$resp_name)) {
          frm_stop("mi(", vn, ") cannot appear in its own model",
                   call. = FALSE)
        }
        mult <- NULL
        if (!is.null(ent$mult)) {
          mult <- check_special_mult(eval(ent$mult, mf, resp$formula_env),
                                     ent$mult, "mi")
        }
        # brms's name for mi(x, idx = g) is its rename(), "mixidxEQg"
        lab <- paste0("mi", vn,
                      if (!is.null(ent$idx)) {
                        paste0("idxEQ", deparse1(ent$idx))
                      },
                      if (!is.null(ent$mult)) {
                        paste0(":", deparse1(ent$mult))
                      } else "")
        X <- cbind(X, matrix(0, nrow(X), 1, dimnames = list(NULL, lab)))
        mi_info[[length(mi_info) + 1L]] <- list(
          var = vn, col = ncol(X), label = lab,
          mult = mult, mult_expr = ent$mult,
          # brms's idxl: the row of `vn` each row here reads, when the
          # two responses do not share their rows (R/subset.R)
          idx_expr = ent$idx,
          idxl = mi_idx_rows(ent, vn, resp, tgt, mf, index_vals, sub_rows)
        )
      }

      # me() terms: zero placeholder columns, filled from the latent
      # values by the objective and the post-fit paths (R/me.R)
      mec <- me_lp_columns(dp, me_fr, X, mf, resp$formula_env)
      X <- mec$X
      me_info <- mec$info

      # Category-specific ordinal effects cs(x): K-1 coefficients per
      # design COLUMN (extras), entering the threshold-specific
      # predictors. A factor or character term is several columns, as it
      # is in brms, so the term's model matrix is kept beside the
      # columns for the newdata paths to recode against.
      cs_info <- list()
      cs_mm <- list()
      if (length(dp[["csterms"]] %||% list())) {
        cs_fam <- cs_target_family(resp$family, dp[["name"]])
        ord_mix <- !is.null(resp$family[["mix"]][["ord"]])
        ordinal <- identical(resp$family[["type"]], "ordinal")
        if (is.null(cs_fam) &&
            (ordinal || !grepl("^mu[0-9]*$", dp[["name"]]))) {
          # brms's sentence first, on any family (sigma ~ cs(w) on a
          # gaussian too). The base build put these offsets in the one
          # slot the densities read, so cs() in disc's formula moved
          # the thresholds as if written in mu's
          # (dev/ordmix-base-behavior.R)
          frm_stop("Category specific effects are only supported for the ",
                   "main parameter 'mu'. cs() moves the thresholds of an ",
                   "ordinal family's latent predictor",
                   if (ord_mix) {
                     ", which in an ordinal mixture is mu1, mu2, ..."
                   },
                   ", and `", dp[["name"]], "` is not that predictor",
                   call. = FALSE)
        }
        if (!ordinal) {
          frm_stop("Category specific effects are not supported for this ",
                   "family. cs() needs an ordinal family: cumulative(), ",
                   "sratio(), cratio(), acat() or hurdle_cumulative()",
                   call. = FALSE)
        }
        if (thres_grouped(resp$family)) {
          # brms 2.23.0 refuses the pair in the same words
          frm_stop("Cannot use category specific effects in models with ",
                   "multiple thresholds. cs() gives each threshold ",
                   "position one coefficient, and with thres(gr = ) the ",
                   "positions differ by group", call. = FALSE)
        }
        # the threshold count, not max(y): thres(x = ) may name
        # categories above the highest one observed. A multivariate
        # frame holds this response's thresholds under its own name,
        # and a structure other than flexible holds fewer parameters
        # than thresholds, so its count is the family's
        tau_nm <- extra_map[[resp$resp_name]][["tau_raw"]] %||% "tau_raw"
        K_cs <- (resp$family[["thres"]][["nthres"]] %||%
                   length(extras[[tau_nm]])) + 1L
        for (ti in seq_along(dp[["csterms"]])) {
          cexpr <- dp[["csterms"]][[ti]]
          cd <- cs_term_design(cexpr, mf, resp$formula_env)
          cs_mm[[ti]] <- cd[setdiff(names(cd), "X")]
          for (j in seq_len(ncol(cd[["X"]]))) {
            csname <- paste0("bcs", length(extras) + 1L)
            extras[[csname]] <- numeric(K_cs - 1L)
            cs_info[[length(cs_info) + 1L]] <- list(
              vals = as.numeric(cd[["X"]][, j]), par = csname,
              label = paste0("cs", cd[["colnames"]][j]),
              expr = cexpr, tid = ti, col = j
            )
          }
        }
        check_cs_identified(X, cs_info, lp_key, mo_info)
        # brms 2.23.0 fits cs() on the two families with ordered
        # thresholds and warns in these words (brms:::check_cs(), once
        # per predictor; dev/rel068-cs-brms.R). The reason is that a
        # row's thresholds can cross, which the second sentence says. It
        # comes after the refusals, so that a model refused for another
        # reason does not also warn
        if (cs_fam[["family"]] %in% ord_ordered_families) {
          frm_warning("Category specific effects for this family should ",
                      "be considered experimental and may have convergence ",
                      "issues. cs() moves each threshold of a row by its ",
                      "own amount, so under ", cs_fam[["family"]], "() a ",
                      "row's thresholds can cross; such a row has no ",
                      "category distribution: its density and fitted() ",
                      "probabilities are NaN, and simulate() gives NA",
                      call. = FALSE)
        }
      }

      # `paste()` recycles to the LONGEST argument, so a design with no
      # columns at all - `b ~ 0 + (1 | g)`, a parameter that is purely a
      # random effect - used to come back from these two lines as the
      # single name "b_". The parameter template then carried a
      # coefficient no linear predictor indexed (`idx` stayed empty), it
      # entered no likelihood, and the outer Hessian was singular in
      # exactly that direction: `vcov()` and every standard error came
      # back NaN with "the model is probably overparameterized".
      cn <- colnames(X)
      # brms's get_model_matrix() renames the columns with check_dup, so
      # `y ~ Intercept + x` is refused there: `(Intercept)` and the
      # covariate `Intercept` would share the parameter name b_Intercept
      if (length(cn)) brms_rename(cn, check_dup = TRUE)
      if (length(cn)) {
        if (!identical(dp[["name"]], "mu")) {
          cn <- paste(dp[["name"]], cn, sep = "_")
        }
        # the shared nu of a Student-t rescor model belongs to no one
        # response, so it is named as in a univariate model
        if (length(spec$responses) > 1 && !isTRUE(dp[["shared"]])) {
          cn <- paste(resp$resp_name, cn, sep = "_")
        }
      } else {
        cn <- character(0)
      }
      if (par_name == "beta") {
        idx <- length(beta_names) + seq_len(ncol(X))
        beta_names <- c(beta_names, cn)
      } else {
        idx <- length(betad_names) + seq_len(ncol(X))
        betad_names <- c(betad_names, cn)
        if (!is.null(dp[["constant"]])) {
          betad_fixed_idx <- c(betad_fixed_idx, idx)
        }
      }

      linpreds[[lp_key]] <- list(
        resp = resp$resp_name,
        dpar = dp[["name"]],
        X = X,
        n_param_cols = n_param_cols,
        param_colnames = param_colnames,
        alias_null = alias_null,
        dropped_colnames = dropped_colnames,
        Z = NULL,               # filled in phase 3
        par = par_name,
        idx = idx,
        offset = if (!is.null(off)) as.numeric(off),
        link = dp[["link"]],
        terms = tt,
        xlevels = xlev,
        contrasts = contr,
        smooths = sm_info,
        gps = gp_info,
        mo = mo_info,
        mi = mi_info,
        me = me_info,
        cs = cs_info,
        cs_mm = cs_mm,
        comp_ids = comp_ids,
        constant = dp[["constant"]],
        # FALSE: the intercept is class "b" and not centered (brms's
        # `0 + Intercept` and `center = FALSE`); see rsv_intercept_fixed()
        center = !isFALSE(dp[["center"]]),
        shared = isTRUE(dp[["shared"]]),
        # newdata gets brms's column of ones too; see rsv_lower_fill()
        rsv_lower = isTRUE(dp[["rsv_lower"]])
      )
    }
    # bf(sigma1 = "sigma2"): the equated dpar is its target's predictor
    # under its own name, so it reads the target's coefficients and adds
    # none. The target may come later in the family's order, hence the
    # second pass.
    for (dp in resp$dpars) {
      if (is.null(dp[["equate"]])) next
      lp <- linpreds[[linpred_key(resp$resp_name, dp[["equate"]])]]
      lp[["dpar"]] <- dp[["name"]]
      lp[["link"]] <- dp[["link"]]
      lp[["equate"]] <- dp[["equate"]]
      linpreds[[linpred_key(resp$resp_name, dp[["name"]])]] <- lp
    }
  }
  mf <- mf_all
  n <- n_all

  ## Phase 2: components -> blocks. Components sharing an |ID| key merge
  ## into one block - unstructured by default, or one gr_cov/gr_prec
  ## Kronecker block when every linked term carries the same known
  ## covariance; everything else gets its own block.
  comp_group <- integer(length(components))
  id_keys <- vapply(components, function(cp) cp$id %||% "", "")
  group_defs <- list()
  for (ci in seq_along(components)) {
    key <- id_keys[ci]
    if (key != "" && key %in% names(group_defs)) {
      group_defs[[key]] <- c(group_defs[[key]], ci)
    } else if (key != "") {
      group_defs[[key]] <- ci
    } else {
      group_defs[[paste0(".solo", ci)]] <- ci
    }
  }

  re_blocks <- list()
  comp_block <- integer(length(components))   # component -> block index
  comp_offset <- integer(length(components))  # coef offset within level
  n_b <- 0L      # parameter space (rr blocks hold factors here)
  n_c <- 0L      # coefficient space (the Z columns)
  n_theta <- 0L
  has_rr <- FALSE
  has_esicar <- FALSE

  for (gd in group_defs) {
    cps <- components[gd]
    mrg_kron <- NULL   # rebuilt at the MERGED dimension, not cps[[1]]'s
    mrg_Qk <- NULL
    if (length(cps) > 1L) {
      cs_set <- unique(vapply(cps, `[[`, "", "covstruct"))
      if ("rr" %in% cs_set) {
        frm_stop("rr() terms cannot share an |ID| key", call. = FALSE)
      }
      lv <- cps[[1]]$levels
      for (cp in cps) {
        if (!identical(cp$levels, lv)) {
          frm_stop("|ID|-linked terms must share identical grouping-factor ",
                   "levels (", cps[[1]]$label, " vs ", cp$label, ")",
                   if (length(sub_rows)) {
                     paste0(". A response with subset() has the levels ",
                            "its own rows carry, and brms's levels of the ",
                            "whole data are not followed here")
                   }, call. = FALSE)
        }
      }
      D <- sum(vapply(cps, `[[`, 0L, "dim"))
      n_levels <- cps[[1]]$n_levels
      # A merged group of gr(cov =) / gr(prec =) terms is ONE Kronecker
      # block of the total merged dimension: b ~ N(0, A (x) Sigma) with
      # Sigma unstructured across the merged coefficients. That is the
      # same joint density as writing the traits long with a single
      # (0 + trait | gr(id, cov = A)) term, so the two spellings agree.
      # check_id_covstructs() already refused a key that mixes
      # structures; what is left to verify is that the cov = expressions
      # RESOLVE to the same matrix, which per-formula environments can
      # make them not do.
      if (length(cs_set) == 1L && cs_set %in% c("gr_cov", "gr_prec")) {
        cs_name <- cs_set
        akey <- if (cs_name == "gr_cov") "aux_A" else "aux_Q"
        for (cp in cps[-1]) {
          if (!same_structural_matrix(cps[[1]][[akey]], cp[[akey]])) {
            frm_stop("|ID|-linked ",
                     if (cs_name == "gr_cov") "gr(cov = )" else "gr(prec = )",
                     " terms must resolve to the same matrix (",
                     cps[[1]]$label, " vs ", cp$label,
                     "). The terms merge into one Kronecker block, which ",
                     "carries a single relationship matrix; put it in ",
                     "data2 so every formula resolves the same object.",
                     call. = FALSE)
          }
        }
        if (cs_name == "gr_cov") {
          mrg_kron <- kron_cov_index(D, n_levels)
        } else {
          mrg_Qk <- kron_prec_parts(cps[[1]]$aux_Q, D)
        }
      } else if (all(cs_set == "us")) {
        cs_name <- "us"
      } else {
        # Unreachable through parse_spec(), which refuses a mixed key
        # up front. Kept because the failure mode it guards is silent:
        # falling through to "us" here would drop a relationship matrix
        # into a density that never reads it.
        frm_stop("|ID|-linked terms mix covariance structures (",
                 paste(cs_set, collapse = ", "),
                 "), which a single merged block cannot carry",
                 call. = FALSE)
      }
      cnms <- unlist(lapply(cps, function(cp) {
        paste0(cp$lp_key, ":", cp$cnms)
      }))
      label <- paste0(paste(vapply(cps, `[[`, "", "label"),
                            collapse = " + "), " [ID]")
      rank_k <- NULL
    } else {
      cp <- cps[[1]]
      D <- cp$dim
      cs_name <- cp$covstruct
      cnms <- cp$cnms
      label <- cp$label
      n_levels <- cp$n_levels
      rank_k <- cp$rank
    }
    if (cs_name == "rr") {
      if (is.null(rank_k) || rank_k > D) {
        frm_stop("rr(): the rank d must not exceed the term dimension (",
                 D, ")", call. = FALSE)
      }
      has_rr <- TRUE
      npar_k <- rr_npar(D, rank_k)
      nb_k <- rank_k * n_levels
    } else if (cs_name %in% c("gp", "hsgp")) {
      # parameter count depends on the gp dimension count, not the
      # block dimension (positions / basis size)
      npar_k <- gp_npar(cps[[1]]$gp_D, cps[[1]]$gp_iso)
      nb_k <- D * n_levels
    } else if (cs_name == "car") {
      # the CAR type decides whether there is a mixing parameter, and
      # esicar's coefficients are its parameters CENTERED, so the frame
      # has to say that expand_b() must run
      npar_k <- car_npar(cps[[1]]$car_type)
      has_esicar <- has_esicar || identical(cps[[1]]$car_type, "esicar")
      nb_k <- D * n_levels
    } else if (cs_name == "spde") {
      npar_k <- spde_npar()
      nb_k <- D * n_levels
    } else {
      npar_k <- covstruct_registry[[cs_name]]$npar(D)
      nb_k <- D * n_levels
    }
    blk_i <- length(re_blocks) + 1L
    ofs <- 0L
    for (k in seq_along(gd)) {
      comp_block[gd[k]] <- blk_i
      comp_offset[gd[k]] <- ofs
      ofs <- ofs + cps[[k]]$dim
    }
    re_blocks[[blk_i]] <- list(
      covstruct = cs_name,
      dim = D,
      rank = rank_k,
      dist_nu = cps[[1]]$dist_nu,
      # gr(g, by = f): which by-level this block is (R/gr-by.R)
      by = cps[[1]]$by,
      n_levels = n_levels,
      b_idx = n_b + seq_len(nb_k),
      c_idx = n_c + seq_len(D * n_levels),
      theta_idx = n_theta + seq_len(npar_k),
      levels = cps[[1]]$levels,
      aux_A = cps[[1]]$aux_A,
      aux_D = cps[[1]]$aux_D,
      aux_D2 = cps[[1]]$aux_D2,
      aux_kron = mrg_kron %||% cps[[1]]$aux_kron,
      aux_Q = cps[[1]]$aux_Q,
      aux_Qk = mrg_Qk %||% cps[[1]]$aux_Qk,
      aux_car = cps[[1]]$aux_car,
      aux_spde = cps[[1]]$aux_spde,
      car_type = cps[[1]]$car_type,
      aux_omega = cps[[1]]$aux_omega,
      gp_D = cps[[1]]$gp_D,
      gp_iso = cps[[1]]$gp_iso,
      gp_vars = cps[[1]]$gp_vars,
      gp_dmax = cps[[1]]$gp_dmax,
      gp_lscale_div = cps[[1]]$gp_lscale_div,
      gp_brms = cps[[1]]$gp_brms,
      cnms = cnms,
      group_name = cps[[1]]$group_name,
      # reformulas writes g/h's nested factor as h:g and brms as g:h; the
      # frame keeps reformulas's factor, and brms's names read this flag
      from_slash = isTRUE(cps[[1]]$from_slash),
      term_label = label,
      dpar = cps[[1]]$dpar,
      components = lapply(seq_along(gd), function(k) {
        out <- list(lp_key = cps[[k]]$lp_key, offset = comp_offset[gd[k]],
                    dim = cps[[k]]$dim, bar = cps[[k]]$bar,
                    mm = cps[[k]]$mm, by = cps[[k]]$by,
                    cnms = cps[[k]]$cnms, label = cps[[k]]$label)
        # cmc = FALSE on this term, which a prediction has to repeat
        if (isTRUE(cps[[k]][["cmc_intercept"]])) {
          out[["cmc_intercept"]] <- TRUE
        }
        out
      })
    )
    n_b <- n_b + nb_k
    n_c <- n_c + D * n_levels
    n_theta <- n_theta + npar_k
  }

  ## Phase 3: per-linpred Z over the full b vector; resolve smooth block
  ## ids; assemble the parameter template.
  for (key in names(linpreds)) {
    lp <- linpreds[[key]]
    if (length(lp[["comp_ids"]])) {
      ii <- integer(0); jj <- integer(0); xx <- numeric(0)
      for (ci in lp[["comp_ids"]]) {
        cp <- components[[ci]]
        bk <- re_blocks[[comp_block[ci]]]
        Tk <- methods::as(cp$Zlocal, "TsparseMatrix")
        loc_col <- Tk@j + 1L
        lev <- (loc_col - 1L) %/% cp$dim + 1L
        cf <- (loc_col - 1L) %% cp$dim + 1L
        glob <- bk[["c_idx"]][(lev - 1L) * bk[["dim"]] + comp_offset[ci] + cf]
        ii <- c(ii, Tk@i + 1L)
        jj <- c(jj, glob)
        xx <- c(xx, Tk@x)
      }
      # a subset() response has its own number of rows
      lp[["Z"]] <- Matrix::sparseMatrix(
        i = ii, j = jj, x = xx,
        dims = c(nrow(components[[lp[["comp_ids"]][1L]]]$Zlocal), n_c))
    }
    if (length(lp[["smooths"]])) {
      lp[["smooths"]] <- lapply(lp[["smooths"]], function(si) {
        si$block_ids <- comp_block[si$comp_ids]
        si
      })
    }
    if (length(lp[["ps_terms"]] %||% list())) {
      lp[["ps_terms"]] <- lapply(lp[["ps_terms"]], function(pt) {
        bk <- re_blocks[[comp_block[pt[["comp_id"]]]]]
        pt[["block_id"]] <- comp_block[pt[["comp_id"]]]
        # coefficient space, not parameter space: an rr block elsewhere
        # in the same model makes the two differ, and the closure reads
        # the expanded vector. `b_idx` is kept as well because the joint
        # covariance's rows are the PARAMETER vector.
        pt[["c_idx"]] <- bk[["c_idx"]]
        pt[["b_idx"]] <- bk[["b_idx"]]
        pt
      })
    }
    if (length(lp[["gps"]])) {
      lp[["gps"]] <- lapply(lp[["gps"]], function(gi) {
        gi$block_id <- comp_block[gi$comp_id]
        gi
      })
    }
    lp[["comp_ids"]] <- NULL
    linpreds[[key]] <- lp
  }

  par_template <- list(beta = stats::setNames(numeric(length(beta_names)),
                                              beta_names))
  if (length(betad_names)) {
    par_template[["betad"]] <- stats::setNames(numeric(length(betad_names)),
                                          betad_names)
  }
  if (n_b) par_template[["b"]] <- numeric(n_b)
  if (n_theta) {
    th0 <- numeric(n_theta)
    for (bk in re_blocks) {
      th0[bk[["theta_idx"]]] <- if (bk[["covstruct"]] == "rr") {
        rr_start(bk[["dim"]], bk[["rank"]])
      } else if (bk[["covstruct"]] == "hsgp") {
        hsgp_start(bk[["gp_D"]], bk[["gp_iso"]])
      } else if (bk[["covstruct"]] == "gp") {
        gp_start(bk[["gp_D"]], bk[["gp_iso"]])
      } else if (bk[["covstruct"]] == "car") {
        car_start(bk[["car_type"]])
      } else if (bk[["covstruct"]] == "spde") {
        spde_start()
      } else {
        covstruct_registry[[bk[["covstruct"]]]]$start(bk[["dim"]])
      }
    }
    par_template[["theta"]] <- th0
  }
  if (n_thetaac) {
    # R-side residual correlation parameters. They are covariance
    # parameters, so they stay OUTER under REML exactly as theta does:
    # REML integrates the mu fixed effects and nothing else.
    thac0 <- numeric(n_thetaac)
    for (ac in autocor) thac0[ac[["theta_idx"]]] <- autocor_start(ac)
    par_template[["thetaac"]] <- thac0
  }
  if (spec$rescor) {
    K <- length(spec$responses)
    par_template[["thetar"]] <- numeric(K * (K - 1L) / 2L)
  }
  if (n_miss) par_template[["miss"]] <- miss_init
  if (!is.null(me_fr)) {
    par_template[["meanme"]] <- me_fr$meanme
    par_template[["logsdme"]] <- me_fr$logsdme
    if (me_fr$n_thetame) {
      par_template[["thetame"]] <- numeric(me_fr$n_thetame)
    }
  }
  for (nm in names(extras)) {
    if (nm %in% names(par_template)) {
      frm_stop("Extra-parameter name collides with the template: ", nm,
               call. = FALSE)
    }
    par_template[[nm]] <- extras[[nm]]
  }

  # Constant dpars: fix their betad entries at link(constant) via map.
  map <- list()
  if (length(betad_fixed_idx)) {
    mp <- seq_along(par_template[["betad"]])
    mp[betad_fixed_idx] <- NA
    map$betad <- factor(mp)
    for (lp in linpreds) {
      if (!is.null(lp[["constant"]])) {
        par_template[["betad"]][lp[["idx"]]] <-
          lp[["link"]]$linkfun(lp[["constant"]])
      }
    }
  }

  # A refusal that depends on the DESIGN rather than on the response
  # cannot be stated from valid_y(): the predictors do not exist yet
  # when that runs. The frame handed over is the assembled one minus
  # the parameter template, which is built below this point.
  frame_so_far <- list(spec = spec, n_obs = n, y = y, y_levels = y_levels,
                       aterm_values = aterm_values, linpreds = linpreds,
                       re_blocks = re_blocks, blocks = blocks,
                       autocor = autocor, data_frame = mf,
                       # the supported way for a check outside this
                       # package to reach a linear predictor, so the key
                       # format of `linpreds` stays ours to change
                       linpred = function(resp_name, dpar) {
                         linpreds[[linpred_key(resp_name, dpar)]]
                       })

  # Checks contributed from another package run first: a feature that
  # adds its own syntax refuses its own data problems before a family
  # gets to speak about the response.
  run_frame_checks(spec, frame_so_far)

  for (resp_ in spec$responses) {
    cf_ <- fam_structure(resp_$family)[["check_frame"]]
    if (!is.null(cf_)) cf_(spec, frame_so_far)
  }

  out <- structure(
    list(spec = spec, n_obs = n, y = y, y_levels = y_levels,
         aterm_values = aterm_values,
         linpreds = linpreds, re_blocks = re_blocks,
         n_c = n_c, has_rr = has_rr, has_expand = has_rr || has_esicar,
         mi_map = mi_map, me = me_fr, blocks = blocks,
         autocor = autocor,
         par_template = par_template, map = map,
         betad_fixed_idx = betad_fixed_idx,
         extra_names = names(extras),
         extra_map = if (length(extra_map)) extra_map,
         predvar_map = predvar_map,
         sparse_x = isTRUE(sparse_x),
         # a refit inside the package keeps the fit's choice
         drop_unused_levels = drop_unused_levels,
         data_frame = mf,
         raw_vars = frame_raw_vars(mf, data, rhs_comb, env),
         na_action = attr(mf, "na.action"),
         # subset(): the rows of `data_frame` each such response uses,
         # and the rows before a univariate model's cut, for nobs()
         subset_rows = if (length(sub_rows)) sub_rows,
         nobs_data = if (!is.null(uni_rows)) n_data),
    class = "frmtmb_frame"
  )
  # brms refuses a group-level effect given twice (frame_re()), and the
  # two copies would share every name this package reports
  brms_check_re_dups(list(frame = out, spec = spec))
  out
}

#' @export
print.frmtmb_frame <- function(x, ...) {
  frm_check_dots(...)
  cat("<frmtmb frame> ", x$n_obs, " observations, ",
      length(x$spec$responses), " response(s)\n", sep = "")
  for (lp in x$linpreds) {
    cat("  ", lp[["resp"]], ".", lp[["dpar"]], ": X[", nrow(lp[["X"]]), " x ",
        ncol(lp[["X"]]), "]", sep = "")
    if (!is.null(lp[["Z"]])) cat(", Z[",
                                 nrow(lp[["Z"]]), " x ", ncol(lp[["Z"]]), "]",
                            sep = "")
    cat(" -> ", lp[["par"]], "[", min(lp[["idx"]]), ":", max(lp[["idx"]]),
      "]", sep = "")
    if (!is.null(lp[["constant"]])) cat("  (fixed at ",
                                        lp[["constant"]], ")", sep = "")
    cat("\n")
  }
  for (bk in x$re_blocks) {
    cat("  RE block: ", bk[["term_label"]], " [",
                           bk[["covstruct"]], "] dim=", bk[["dim"]],
        " levels=", bk[["n_levels"]], "\n", sep = "")
  }
  for (ac in x$autocor %||% list()) {
    cat("  R-side: ", ac[["block_label"]], " [",
                                     ac[["struct"]], "] times=", ac[["d"]],
        " groups=", ac[["n_groups"]], " patterns=", length(ac[["patterns"]]),
        "\n", sep = "")
  }
  cat("  parameters:",
      paste0(names(x$par_template), "(",
             lengths(x$par_template), ")", collapse = ", "), "\n")
  invisible(x)
}
