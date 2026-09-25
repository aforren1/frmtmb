# Confidence intervals, convergence diagnostics, model comparison.

#' One-line label for a fit including dpar formulas, so two models that
#' share a primary formula (e.g. plain vs distributional) stay
#' distinguishable in anova tables.
#'
#' @noRd
model_label <- function(fit) {
  one <- function(bform) {
    parts <- deparse1(bform$formula)
    for (nm in names(bform$pforms)) {
      parts <- c(parts, deparse1(bform$pforms[[nm]]))
    }
    for (nm in names(bform$pfix)) {
      parts <- c(parts, paste(nm, "=", bform$pfix[[nm]]))
    }
    paste(parts, collapse = ", ")
  }
  bf0 <- fit$bform
  if (inherits(bf0, "frmtmb_mvformula")) {
    paste(vapply(bf0$forms, one, ""), collapse = " + ")
  } else {
    one(bf0)
  }
}

#' Names of the outer (optimized) parameters, in obj$par order, with the
#' template component each one came from: template component order,
#' minus `random` components, minus mapped entries. The component is
#' what lets a natural-scale alias (a `theta` or `thetaac` position) be
#' turned into a position in this vector.
#'
#' @noRd
outer_par_map <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  # mirror the MakeADFun random= construction in fit_assembled: b and
  # the mi() latent component are always inner, beta under REML or
  # control profile = TRUE
  random <- c("b", "miss")
  if (fit$REML || isTRUE(fit$control$profile)) {
    random <- c(random, "beta")
  }
  nm <- character(0)
  comp <- character(0)
  for (cp in names(tpl)) {
    if (cp %in% random) next
    v <- names(tpl[[cp]])
    if (is.null(v)) v <- paste0(cp, "_", seq_along(tpl[[cp]]))
    if (cp == "betad" && length(fit$frame[["betad_fixed_idx"]])) {
      v <- v[-fit$frame[["betad_fixed_idx"]]]
    }
    nm <- c(nm, v)
    comp <- c(comp, rep(cp, length(v)))
  }
  list(names = nm, comp = comp)
}

#' @noRd
outer_par_names <- function(fit) outer_par_map(fit)$names

## Addressing a parameter by name. Three vocabularies meet here: the
## internal names of outer_par_map() (`tarsus_(Intercept)`, `theta_1`),
## the parenthesis-free spelling hypothesis() and variables() use
## (`tarsus_Intercept`), and the natural-scale summaries
## (`sd_dam__tarsus.muIntercept`). Every parm-style argument accepts all
## three; output rows keep the internal name.

#' The comparison form of a parameter name: parentheses dropped, which
#' is the only difference between the internal spelling and the one
#' hypothesis() expressions must use (a name in an R expression cannot
#' carry parentheses).
#'
#' @noRd
par_name_bare <- function(x) gsub("[()]", "", x)

#' Positions of `x` in `nm`, exact first and then with parentheses
#' dropped on both sides. `NA` where nothing matches. Parentheses only
#' ever wrap `Intercept`, so a stripped match is unique in practice; the
#' ambiguity check guards the case rather than assuming it.
#'
#' @noRd
match_par_name <- function(x, nm) {
  idx <- match(x, nm)
  if (!anyNA(idx)) return(idx)
  bare <- par_name_bare(nm)
  for (k in which(is.na(idx))) {
    hit <- which(bare == par_name_bare(x[k]))
    if (length(hit) > 1L) {
      frm_stop("Parameter name '", x[k], "' is ambiguous once parentheses ",
               "are dropped: it matches ", paste(nm[hit], collapse = ", "),
               ". Write the full internal name", call. = FALSE)
    }
    if (length(hit) == 1L) idx[k] <- hit
  }
  idx
}

#' Natural-scale names that address exactly ONE internal parameter, as
#' positions in `outer_par_map()`. A standard deviation is one log-sd
#' entry of `theta` (`sd_idx` from the covariance-structure registry).
#' A correlation is one internal parameter only when the block has a
#' single non-sd `theta` entry - true of a 2x2 `us` block, `cs` and
#' `ar1`, false of a wider `us` whose Cholesky terms mix - and an
#' autocorrelation parameter only when its block has one `thetaac`
#' entry. Anything else is a function of several internal parameters
#' and is refused rather than aliased to one of them.
#'
#' @noRd
par_alias_index <- function(fit) {
  map <- outer_par_map(fit)
  th_pos <- which(map$comp == "theta")
  ac_pos <- which(map$comp == "thetaac")
  nms <- character(0)
  pos <- integer(0)
  add <- function(nm, p) {
    if (nm %in% nms || is.na(p) || length(p) != 1L) return(invisible())
    nms <<- c(nms, nm)
    pos <<- c(pos, p)
  }
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    # the same blocks hyp_env_vals names, for the same reason
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "car", "spde")) next
    reg <- covstruct_registry[[bk[["covstruct"]]]]
    if (is.null(reg) || is.null(reg$sd_idx)) next
    si <- tryCatch(as.integer(reg$sd_idx(bk[["dim"]])),
                   error = function(e) integer(0))
    if (!length(si)) next
    g <- brms_group_name(bk)
    tn <- brms_re_rnames(fit, bk)
    at <- function(i) {
      ti <- bk[["theta_idx"]][i]
      if (is.na(ti) || ti < 1L || ti > length(th_pos)) NA_integer_ else
        th_pos[ti]
    }
    for (j in seq_along(tn)) {
      add(paste0("sd_", g, "__", tn[j]), at(si[min(j, length(si))]))
    }
    if (length(tn) > 1L) {
      rest <- setdiff(seq_len(reg$npar(bk[["dim"]])), si)
      if (length(rest) == 1L) {
        for (j in seq_len(length(tn) - 1L)) {
          for (k in seq(j + 1L, length(tn))) {
            add(paste0("cor_", g, "__", tn[j], "__", tn[k]), at(rest))
          }
        }
      }
    }
  }
  for (ac in fit$frame[["autocor"]] %||% list()) {
    if (length(ac[["theta_idx"]]) != 1L) next
    nat <- autocor_natural(fit$estimates[["thetaac"]][ac[["theta_idx"]]], ac)
    if (length(nat) != 1L) next
    p <- ac[["theta_idx"]][1L]
    add(names(nat)[1L],
        if (p >= 1L && p <= length(ac_pos)) ac_pos[p] else NA_integer_)
  }
  stats::setNames(pos, nms)
}

#' Bare nonlinear-parameter names, as positions in `outer_par_map()`.
#' An `nl` parameter declared `la ~ 1` has exactly one coefficient,
#' `la_(Intercept)`, so the bare `la` names it without ambiguity. That
#' is the one-to-one rule the `sd_`/`ar[1]` aliases already follow. A
#' parameter
#' with a design matrix wider than one column names several, and is
#' reported in `ambiguous` for the caller to refuse by name.
#'
#' Only nonlinear parameters, deliberately: a distributional parameter's
#' bare name (`sigma`) is already the NATURAL-scale summary that
#' [hypothesis()] reports, and aliasing it to the internal
#' `sigma_(Intercept)` would make the same word mean a value and its
#' link transform in two neighbouring arguments.
#'
#' @noRd
nlpar_bare_alias <- function(fit) {
  spec <- fit$frame[["spec"]]
  nm <- outer_par_names(fit)
  tpl <- fit$frame[["par_template"]]
  mv <- length(spec$responses) > 1L
  pos <- integer(0)
  pnm <- character(0)
  amb <- list()
  for (resp in spec$responses) {
    for (np in resp$nlpars %||% character(0)) {
      bare <- if (mv) paste0(resp$resp_name, "_", np) else np
      # never displace a real parameter name, and never claim a bare
      # name twice (two responses can declare the same nlpar)
      if (bare %in% nm || bare %in% pnm || bare %in% names(amb)) next
      for (lp in fit$frame[["linpreds"]]) {
        if (!identical(lp[["resp"]], resp$resp_name) ||
            !identical(lp[["dpar"]], np)) {
          next
        }
        full <- names(tpl[[lp[["par"]]]])[lp[["idx"]]]
        hit <- stats::na.omit(match(full, nm))
        if (length(hit) == 1L) {
          pnm <- c(pnm, bare)
          pos <- c(pos, hit)
        } else if (length(hit) > 1L) {
          amb[[bare]] <- nm[hit]
        }
        break
      }
    }
  }
  list(pos = stats::setNames(pos, pnm), ambiguous = amb)
}

#' Fold the bare nonlinear-parameter aliases into a partly resolved
#' index vector, refusing a bare name that stands for several
#' coefficients rather than picking one of them.
#'
#' @noRd
apply_nlpar_alias <- function(fit, x, idx) {
  if (!anyNA(idx)) return(idx)
  al <- nlpar_bare_alias(fit)
  hit <- match(x, names(al$pos))
  took <- which(is.na(idx) & !is.na(hit))
  idx[took] <- al$pos[hit[took]]
  bad <- intersect(x[is.na(idx)], names(al$ambiguous))
  if (length(bad)) {
    frm_stop("Nonlinear parameter '", bad[1L], "' has more than one ",
             "coefficient, so the bare name does not identify one of them. ",
             "Name the coefficient in full: ",
             paste(al$ambiguous[[bad[1L]]], collapse = ", "), call. = FALSE)
  }
  idx
}

#' Resolve a `parm`-style argument to positions in `outer_par_map()`,
#' accepting the internal names, their parenthesis-free spelling, and
#' the one-to-one natural-scale aliases. An alias is reported, because
#' what gets profiled and returned is the internal parameter it names,
#' on the internal scale.
#'
#' @noRd
resolve_par_index <- function(fit, parm, what) {
  map <- outer_par_map(fit)
  nm <- map$names
  if (is.numeric(parm)) return(as.integer(parm))
  idx <- match_par_name(parm, nm)
  # a bare nlpar is a SPELLING of one internal parameter, like dropping
  # the parentheses, so it resolves silently; the natural-scale aliases
  # below name a different scale and say so
  idx <- apply_nlpar_alias(fit, parm, idx)
  if (anyNA(idx)) {
    # brms's `b_` spelling, which variables() and hypothesis() now use,
    # is a spelling of one coefficient too, so it resolves silently
    # a natural-scale name such as sigma is a transform of its
    # coefficient, not a spelling of it, so it does not resolve here
    tab <- brms_coef_table(fit)
    bn <- match(parm, ifelse(tab$natural, NA_character_, tab$brms))
    took <- which(is.na(idx) & !is.na(bn))
    if (length(took)) {
      idx[took] <- match_par_name(estimated_coef_names(fit)[bn[took]], nm)
    }
  }
  if (anyNA(idx)) {
    alias <- par_alias_index(fit)
    hit <- match(parm, names(alias))
    took <- which(is.na(idx) & !is.na(hit))
    idx[took] <- alias[hit[took]]
    if (length(took)) {
      frm_message(what, "(): ",
                  paste0("'", parm[took], "' is ", nm[idx[took]],
                         collapse = ", "),
                  ". The result is on that parameter's internal ",
                  "(unconstrained) scale, not the natural one; ",
                  "confint_varcorr() and hypothesis() report the natural ",
                  "scale.")
    }
  }
  if (anyNA(idx)) {
    bad <- parm[is.na(idx)]
    known <- variables(fit)
    if (any(bad %in% known)) {
      b <- bad[bad %in% known][1L]
      frm_stop("'", b, "' is a natural-scale summary rather than a fitted ",
               "parameter, and it does not stand for a single internal one ",
               "here (a correlation of a wider us() block mixes several, ",
               "and a response-scale summary such as sigma is a transform ",
               "of one). Use hypothesis(fit, \"", b, " = 0\", method = ",
               "'profile'), which profiles the combination itself, read ",
               "confint_varcorr() for natural-scale variance components, ",
               "or name the internal parameter: confint(fit) lists them",
               call. = FALSE)
    }
    frm_stop("Unknown parameter(s) in ", what, "(parm =): ",
             paste(bad, collapse = ", "), ". Available: ",
             paste(nm, collapse = ", "),
             ". Parentheses may be dropped, intercept-only nonlinear ",
             "parameters may be named bare, and the one-to-one natural-scale ",
             "names of variables() (sd_<group>__<term>, and a correlation ",
             "with a single internal parameter) are accepted as aliases",
             call. = FALSE)
  }
  idx
}

#' Confidence intervals for frmtmb fits
#'
#' Covariance parameters (`theta_*`) are reported on their internal
#' (unconstrained) scale.
#'
#' @param object A `frmtmb_fit`.
#' @param parm Parameter names (see `rownames` of the Wald result) or
#'   indices. Required for `"profile"` and `"uniroot"`; defaults to all
#'   parameters for `"wald"` and `"boot"`. Three spellings address the
#'   same parameter: the internal name (`tarsus_(Intercept)`,
#'   `theta_1`), that name without its parentheses
#'   (`tarsus_Intercept`, the spelling [hypothesis()] and [variables()]
#'   use), and a natural-scale name that stands for exactly one
#'   internal parameter, `sd_<group>__<term>` and a correlation whose
#'   block has a single internal correlation parameter. The bare name of
#'   an intercept-only nonlinear parameter (`la` for `la_(Intercept)`)
#'   is a fourth spelling of the same internal parameter, so it resolves
#'   silently; a nonlinear parameter with several coefficients is
#'   refused rather than resolved to one of them. An alias is
#'   reported and the row keeps the internal name, because the interval
#'   is on the internal scale: see [confint_varcorr()] and
#'   [hypothesis()] for natural-scale intervals. A natural-scale
#'   summary that is a function of several internal parameters (the
#'   correlations of a wider `us` block, an ICC) has no single
#'   parameter to work on and is refused, naming
#'   `hypothesis(method = "profile")`.
#' @param level Confidence level.
#' @param method `"wald"` (fast, from the sdreport covariance),
#'   `"profile"` (likelihood profile via [TMB::tmbprofile()]),
#'   `"uniroot"` (likelihood-root search via [TMB::tmbroot()]), or
#'   `"boot"` (parametric-bootstrap percentile intervals through
#'   [frm_bootstrap()], the `lme4::confint(method = "boot")` analog;
#'   like the other methods it works on the internal parameter scale).
#'   `"Wald"` is accepted as an alias for `"wald"`.
#' @param nsim,seed Bootstrap draws and seed for `method = "boot"`.
#' @param vcov `method = "wald"` only: a covariance matrix over the
#'   whole outer parameter vector to use in place of the model-based
#'   one - [vcov_cluster()] with `full = TRUE`, or a function of the
#'   fit returning such a matrix. A matrix that carries reference
#'   degrees of freedom (as `vcov_cluster()`'s does, `G - 1`) switches
#'   the interval from a normal to a `t` quantile.
#' @param ... Passed to the TMB profiling functions, or to
#'   [frm_bootstrap()] for `method = "boot"` (e.g. `re_formula`).
#' @return A matrix with columns `lwr`, `upr`, `est`.
#'
#' @srrstats {RE4.3} Confidence intervals on the model coefficients are
#'   returned by `confint()`, by four methods: Wald from the sdreport
#'   covariance, likelihood profile, likelihood-root search, and
#'   parametric bootstrap. `confint_varcorr()` gives natural-scale
#'   intervals for standard deviations and correlations. Row names match
#'   those of `vcov(full = TRUE)`, which the test suite asserts.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # Wald intervals for every parameter, covariance ones included
#' confint(fit)
#'
#' # the likelihood profile does not assume a quadratic log-likelihood,
#' # so it is the one to trust for a variance component
#' confint(fit, parm = "theta_1", method = "profile")
#'
#' # confint_varcorr() puts the same information on the SD scale
#' confint_varcorr(fit)
#' \donttest{
#' # a parametric bootstrap, the lme4 confint(method = "boot") analog
#' confint(fit, parm = "x", method = "boot", nsim = 50, seed = 1)
#' }
#' @export
confint.frmtmb_fit <- function(object, parm = NULL, level = 0.95,
                               method = c("wald", "Wald", "profile",
                                          "uniroot", "boot"),
                               nsim = 500, seed = NULL, vcov = NULL,
                               ...) {
  method <- frm_match_arg(method)
  if (method == "Wald") method <- "wald"
  # The dots are forwarded on three of the four methods and swallowed on
  # the fourth, so they are checked against whatever THIS method really
  # passes them to. Before this, confint(method = "wald", re_formula =)
  # changed nothing and said nothing.
  frm_check_dots(..., .allow = switch(method,
    boot = names(formals(frm_bootstrap)),
    profile = names(formals(TMB::tmbprofile)),
    uniroot = names(formals(TMB::tmbroot)),
    NULL))
  # A length-2 level makes a length-2 quantile, which then RECYCLES
  # against the parameter vector: rows 1 and 3 of one table came back at
  # 90% and rows 2 and 4 at 95%, with nothing in the output recording
  # it. A level outside (0, 1) produced a table of NaN or of Inf.
  check_probability(level, "level")
  check_count(nsim, "nsim", min = 1L)
  if (!is.null(parm) && !is.character(parm)) {
    frm_stop("`parm` must be a character vector of parameter names, or NULL ",
             "for all of them, not ", arg_desc(parm), call. = FALSE)
  }
  if (!is.null(vcov) && method != "wald") {
    frm_stop("confint(vcov = ) applies to method = 'wald' only: ",
             "method = '", method, "' does not go through a covariance ",
             "matrix", call. = FALSE)
  }
  nm <- outer_par_names(object)
  est <- object$opt$par
  a <- (1 - level) / 2

  idx <- if (is.null(parm)) {
    seq_along(nm)
  } else {
    resolve_par_index(object, parm, "confint")
  }

  if (method == "wald") {
    q <- stats::qnorm(a)
    if (is.null(vcov)) {
      sdr <- sdr_of(object)
      se <- sqrt(diag(sdr$cov.fixed))
      # profiled betas are absent from BOTH opt$par and par.fixed in
      # current RTMB; the fallback stays as a defensive alignment only
      if (length(est) != length(se)) est <- sdr$par.fixed
    } else {
      rv <- resolve_vcov_arg(object, vcov, "confint")
      se <- sqrt(pmax(0, diag(rv$V)))
      if (!is.null(rv$df)) q <- stats::qt(a, rv$df)
    }
    ci <- cbind(lwr = est + q * se,
                upr = est - q * se,
                est = est)
    rownames(ci) <- nm
    return(ci[idx, , drop = FALSE])
  }

  if (method == "boot") {
    # refits share the fit's control, so opt$par lines up with nm
    # (profile = TRUE excludes beta from both)
    bs <- frm_bootstrap(object, FUN = function(f) f$opt$par,
                        nsim = nsim, seed = seed, ...)
    ci <- cbind(lwr = apply(bs$t, 2, stats::quantile, a, na.rm = TRUE),
                upr = apply(bs$t, 2, stats::quantile, 1 - a,
                            na.rm = TRUE),
                est = est)
    rownames(ci) <- nm
    return(ci[idx, , drop = FALSE])
  }

  if (isTRUE(object$control$profile)) {
    frm_stop("confint(method = '", method, "') needs a fit without ",
             "frmtmb_control(profile = TRUE)", call. = FALSE)
  }
  if (is.null(parm)) {
    frm_stop("`parm` is required for method = '", method, "'", call. = FALSE)
  }
  ci <- matrix(NA_real_, length(idx), 3,
               dimnames = list(nm[idx], c("lwr", "upr", "est")))
  # An importance-corrected fit profiles its own FROZEN tape, whose
  # proposal sits at the estimate; a bound far from there can come from
  # a region the proposal no longer covers, and the effective sample
  # sizes the fit reports describe the anchor only. Rebuilt once for
  # the whole call, then read at each bound.
  prop <- if (!is.null(object$importance)) imp_frozen_proposal(object)
  for (k in seq_along(idx)) {
    i <- idx[k]
    if (method == "profile") {
      pr <- TMB::tmbprofile(object$obj, name = i, trace = FALSE, ...)
      ci[k, 1:2] <- unname(stats::confint(pr, level = level))
    } else {
      r <- TMB::tmbroot(object$obj, name = i,
                        target = 0.5 * stats::qchisq(level, df = 1), ...)
      ci[k, 1:2] <- unname(r)
    }
    ci[k, 3] <- est[i]
    if (!is.null(prop)) {
      imp_profile_ess_warn(object, prop, nm[i], i, ci[k, 1:2])
    }
  }
  ci
}

#' Transformed-scale Wald rows for the covariance parameters of one fit:
#' one row per SD/range (log scale) and per correlation (Fisher-z
#' scale), with the delta-method se on that scale. These scales are
#' where a normal approximation is defensible, which makes the rows the
#' right currency both for confint_varcorr's intervals and for Rubin
#' pooling across imputations in frm_multiple.
#'
#' The correlation the Fisher-z transform is clamped at inside the
#' finite-difference jacobian. A correlation this close to +/-1 is on
#' the boundary of the parameter space: `atanh` runs away there and the
#' clamp flattens the jacobian row, so those rows report NA bounds
#' rather than a zero-width interval at the clamp value. One constant so
#' the clamp and the boundary test cannot drift apart.
#'
#' @noRd
varcorr_cor_clamp <- 0.9999

#' @noRd
varcorr_trans_rows <- function(fit) {
  sdr <- sdr_of(fit)
  Vfull <- sdr$cov.fixed
  th_pos <- which(rownames(Vfull) == "theta")
  th <- fit$estimates[["theta"]]
  rows <- list()
  add <- function(term, type, est_t, se_t, bk) {
    rows[[length(rows) + 1L]] <<- data.frame(
      block = bk[["term_label"]], term = term, type = type,
      est_t = est_t, se_t = se_t
    )
  }
  for (bk in fit$frame[["re_blocks"]]) {
    Vth <- Vfull[th_pos[bk[["theta_idx"]]], th_pos[bk[["theta_idx"]]],
                 drop = FALSE]
    t0 <- th[bk[["theta_idx"]]]
    if (bk[["covstruct"]] == "smooth") {
      add("sd(wiggle)", "sd", t0[1], sqrt(Vth[1, 1]), bk)
      next
    }
    if (bk[["covstruct"]] %in% c("gp", "hsgp")) {
      se_t <- sqrt(diag(Vth))
      add("sd(gp)", "sd", t0[1], se_t[1], bk)
      # hsgp estimates the lengthscale on brms's rescaled inputs, but the
      # reported range belongs in data units. The scale factor is a data
      # constant, so the shift on the log scale is exact and the se rides
      # through unchanged. The exact gp keeps the raw scale (dmax NULL).
      log_dmax <- log(bk[["gp_dmax"]] %||% 1)
      # iso: one shared range; otherwise one per dimension
      nr <- length(t0) - 1L
      for (j in seq_len(nr)) {
        term_j <- if (nr == 1L) "range(gp)" else {
          paste0("range(gp, ", bk[["gp_vars"]][j], ")")
        }
        add(term_j, "range", t0[1 + j] + log_dmax, se_t[1 + j], bk)
      }
      next
    }
    if (bk[["covstruct"]] == "car") {
      se_t <- sqrt(diag(Vth))
      add("sd(car)", "sd", t0[1], se_t[1], bk)
      if (length(t0) > 1L) {
        # brms's names for the two mixing parameters, both on (0, 1)
        nm <- if (identical(bk[["car_type"]], "bym2")) "rhocar" else "car"
        add(nm, "prop", t0[2], se_t[2], bk)
      }
      next
    }
    if (bk[["covstruct"]] == "spde") {
      # sigma and range are analytic functions of (log tau, log kappa):
      # log sigma = -log tau - log kappa - log(4 pi) / 2 and
      # log range = log(8) / 2 - log kappa, so the delta method is one
      # exact linear map with no differencing
      g_sd <- c(-1, -1)
      add("sd(spde)", "sd", log(spde_sd(t0)),
          sqrt(max(drop(g_sd %*% Vth %*% g_sd), 0)), bk)
      add("range(spde)", "range", log(spde_range(t0)),
          sqrt(Vth[2, 2]), bk)
      next
    }
    if (bk[["covstruct"]] == "equalto") next   # nothing estimated
    # g(theta): log-sds then atanh-correlations, via the block's vcov.
    # The clamp keeps the CENTRAL DIFFERENCES finite; it must not reach
    # the reported estimate, which is why est0 below is computed
    # separately from the covariance itself.
    gfun <- function(tt) {
      V <- covstruct_registry[[bk[["covstruct"]]]]$vcov(tt, bk)
      sds <- sqrt(diag(V))
      out <- log(sds)
      if (nrow(V) > 1) {
        C <- stats::cov2cor(V)
        out <- c(out, atanh(pmin(pmax(C[lower.tri(C)],
                                      -varcorr_cor_clamp),
                                 varcorr_cor_clamp)))
      }
      out
    }
    g0 <- gfun(t0)
    # numeric jacobian, central differences
    J <- vapply(seq_along(t0), function(i) {
      h <- 1e-5 * max(abs(t0[i]), 1)
      tp <- t0; tp[i] <- tp[i] + h
      tm <- t0; tm[i] <- tm[i] - h
      (gfun(tp) - gfun(tm)) / (2 * h)
    }, numeric(length(g0)))
    J <- matrix(J, nrow = length(g0))
    se_g <- sqrt(pmax(diag(J %*% Vth %*% t(J)), 0))

    d <- bk[["dim"]]
    # Transformed-scale estimates read off the covariance without the
    # clamp. A component ON the boundary of its parameter space - a
    # standard deviation collapsed to zero, a correlation at +/-1 - maps
    # to +/-Inf here, which is the honest answer: the estimate
    # back-transforms to 0 or to +/-1, and the Wald interval around an
    # infinite point does not exist. Without this the clamp inside gfun
    # flattened the jacobian row to zero and the block reported
    # lwr == est == upr == 0.9999, a zero-width interval AT the clamp.
    V0 <- covstruct_registry[[bk[["covstruct"]]]]$vcov(t0, bk)
    sds0 <- sqrt(diag(V0))
    est0 <- log(sds0)
    bd <- sds0 <= 0
    if (d > 1L) {
      C0 <- stats::cov2cor(V0)
      r0 <- C0[lower.tri(C0)]
      est0 <- c(est0, atanh(r0))
      # the same threshold the clamp uses, on purpose: a correlation the
      # clamp touches has a zero jacobian row and would otherwise report
      # a zero-width interval AT the clamp
      bd <- c(bd, abs(r0) >= varcorr_cor_clamp)
    }
    # A non-finite se (an inverted Hessian that did not invert) is no
    # more usable than a boundary estimate; both become NA bounds.
    se_g[bd | !is.finite(est0) | !is.finite(se_g)] <- NA_real_

    n_sd <- length(g0) - if (d > 1) d * (d - 1) / 2 else 0
    for (i in seq_len(n_sd)) {
      add(bk[["cnms"]][min(i,
                           length(bk[["cnms"]]))], "sd", est0[i], se_g[i], bk)
    }
    if (d > 1) {
      pairs <- which(lower.tri(diag(d)), arr.ind = TRUE)
      for (k in seq_len(nrow(pairs))) {
        add(paste0("cor(", bk[["cnms"]][pairs[k, 2]], ",",
                   bk[["cnms"]][pairs[k, 1]], ")"),
            "cor", est0[n_sd + k], se_g[n_sd + k], bk)
      }
    }
  }
  acr <- autocor_trans_rows(fit)
  if (!length(rows) && is.null(acr)) return(NULL)
  out <- if (length(rows)) do.call(rbind, rows) else NULL
  out <- if (is.null(out)) acr else if (is.null(acr)) out else {
    rbind(out, acr)
  }
  rownames(out) <- NULL
  out
}

#' Back-transform to the natural scale (elementwise over types): log for
#' scales, Fisher-z for correlations, logit for the CAR mixing
#' proportions (which live on (0, 1), as they do in brms), and the
#' identity for the "raw" components - the AR/MA coefficients of a
#' higher-order residual process, which are bounded by stationarity as a
#' SET but not one at a time.
#'
#' @noRd
varcorr_untrans <- function(type, v) {
  ifelse(type == "raw", v,
         ifelse(type == "cor", tanh(v),
                ifelse(type == "prop", 1 / (1 + exp(-v)), exp(v))))
}

#' Natural-scale confidence intervals for covariance parameters
#'
#' Wald intervals for random-effect standard deviations (on the log
#' scale, back-transformed) and correlations (on the Fisher-z scale,
#' back-transformed), delta-method-propagated from the internal `theta`
#' covariance. One row per SD and per correlation of every block.
#'
#' A component sitting on the boundary of its parameter space - a
#' standard deviation collapsed to zero, or a correlation at `+/-1` -
#' has no interval on these scales, because the transform is infinite
#' there. Those rows report the estimate with `NA` bounds and warn,
#' rather than a zero-width interval at an arbitrary clamp. A
#' bootstrap ([hypothesis()] with `method = "boot"`) or a likelihood
#' profile of the underlying `theta` is the alternative.
#'
#' @param fit A `frmtmb_fit`.
#' @param level Confidence level.
#' @return A data frame with columns `block`, `term`, `type`,
#'   `estimate`, `lwr`, `upr`. Boundary components carry `NA` bounds.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
#' u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
#' dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
#' fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
#'
#' # one row per SD and per correlation, on the scale they are read on
#' confint_varcorr(fit)
#'
#' # confint() reports the same parameters on their internal scale, so
#' # the bounds there are log-SDs and Fisher-z correlations
#' confint(fit)[grep("^theta", rownames(confint(fit))), ]
#'
#' # a fit with no random effects has no covariance parameters
#' confint_varcorr(frm(bf(y ~ x) + gaussian(), data = dd))
#' @export
confint_varcorr <- function(fit, level = 0.95) {
  tr <- varcorr_trans_rows(fit)
  if (is.null(tr)) return(NULL)
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(
    block = tr$block, term = tr$term, type = tr$type,
    estimate = varcorr_untrans(tr$type, tr$est_t),
    lwr = varcorr_untrans(tr$type, tr$est_t - z * tr$se_t),
    upr = varcorr_untrans(tr$type, tr$est_t + z * tr$se_t)
  )
  label <- function(i) paste0(tr$block[i], " ", tr$term[i])
  bad <- which(is.na(tr$se_t))
  if (length(bad)) {
    frm_warning("No interval for ", length(bad), " component",
                if (length(bad) > 1L) "s" else "",
                " on the boundary of the parameter space (a standard ",
                "deviation at zero, or a correlation at +/-1). The interval ",
                "is a Wald interval on the log / Fisher-z scale, which is ",
                "infinite there, so it is reported as NA rather than as a ",
                "zero-width interval at an arbitrary clamp. The estimates ",
                "stand; for an interval use hypothesis(method = \"boot\") ",
                "or confint(method = \"profile\") on the theta parameter. ",
                "Affected: ",
                paste(vapply(bad, label, ""), collapse = "; "), call. = FALSE)
  }
  # An se above 10 on the LOG scale spans more than 17 orders of
  # magnitude each way: the component is not identified by the data, and
  # the numbers are noise dressed as an interval. Reporting them is
  # still better than dropping them (the point estimate is usually
  # fine), but saying nothing is not. The bounded types are exempt: a
  # correlation and a mixing proportion back-transform into (-1, 1) and
  # (0, 1), so a huge se there widens the interval to the whole
  # parameter space and stops, which reads as the non-identification it
  # is.
  wide <- which(!is.na(tr$se_t) & tr$type %in% c("sd", "range") &
                  tr$se_t > 10)
  if (length(wide)) {
    frm_warning("Uninformative interval for ", length(wide), " component",
                if (length(wide) > 1L) "s" else "",
                ": the standard error on the log scale exceeds 10, so the ",
                "reported bounds span many orders of magnitude and the ",
                "data do not identify the component. diagnose() reports ",
                "the boundary and curvature checks. Affected: ",
                paste(vapply(wide, label, ""), collapse = "; "),
                call. = FALSE)
  }
  out
}

#' Effective degrees of freedom of the smooth blocks: for an iid wiggly
#' block, `edf = k - tr(posterior cov)/prior variance` (the ridge
#' identity).
#'
#' @noRd
smooth_edf <- function(fit) {
  blocks <- Filter(function(bk) bk[["covstruct"]] == "smooth",
                   fit$frame[["re_blocks"]])
  if (!length(blocks)) return(NULL)
  sdr <- sdr_of(fit)
  dcr <- sdr$diag.cov.random
  if (is.null(dcr)) return(NULL)
  # par.random holds the `random` components in template order; b entries
  # are the ones named "b"
  b_pos <- which(names(sdr$par.random) == "b")
  th <- fit$estimates[["theta"]]
  out <- vapply(blocks, function(bk) {
    prior_var <- exp(th[bk[["theta_idx"]]])^2
    k <- bk[["dim"]]
    # +1 null-space columns live in beta; conventionally reported as the
    # penalized-part edf
    k - sum(dcr[b_pos[bk[["b_idx"]]]]) / prior_var
  }, numeric(1))
  stats::setNames(out, vapply(blocks, `[[`, "", "term_label"))
}

# Families whose linear predictor lives on a bounded probability scale,
# so an unbounded coefficient is evidence of separated data rather than
# of a large effect.
separation_families <- c("binomial", "bernoulli", "beta_binomial",
                         "zero_inflated_binomial",
                         "zero_inflated_beta_binomial")

#' Complete (or quasi-complete) separation: the maximum likelihood sits
#' at infinity, so the optimizer stops wherever its tolerances bite and
#' reports a huge coefficient with a standard error to match. lme4 and
#' glmmTMB both flag the pair rather than either half, because a
#' genuinely large effect on a well-populated cell keeps a small se.
#' `[glmmTMB diagnose()]`
#'
#' @noRd
diagnose_separation <- function(fit, ps) {
  rows <- list()
  for (lp in fit$frame[["linpreds"]]) {
    fam <- fit$spec$responses[[lp[["resp"]]]]$family
    if (!fam[["family"]] %in% separation_families) next
    if (!lp[["dpar"]] %in% (fam[["primary_dpars"]] %||% "mu")) next
    if (is.null(lp[["X"]]) || !ncol(lp[["X"]])) next
    est <- ps$est[[lp[["par"]]]][lp[["idx"]]]
    se <- ps$se[[lp[["par"]]]][lp[["idx"]]]
    hit <- which(abs(est) > 10 & (!is.finite(se) | se > 10))
    for (i in hit) {
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = paste0(coef_block_key(fit, lp), ": ",
                           colnames(lp[["X"]])[i]),
        estimate = est[i], std.error = se[i]
      )
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' A DISTRIBUTIONAL parameter whose maximum likelihood sits outside its
#' own range.
#'
#' Same evidence as `diagnose_separation()` and a different cause. An
#' estimate far out on the link with a standard error bigger than
#' itself means the likelihood was still rising when the optimizer's
#' tolerances bit, so what is reported is where it stopped rather than
#' an estimate. `student()`'s `nu` does it on any data with no heavy
#' tails: the student-t reaches the gaussian only in the limit
#' `nu -> Inf`, so the maximum is never attained and `nu` comes back at
#' whatever the last step reached.
#'
#' Measured over 84 student fits (seven error laws, `n` of 60 and 200,
#' six replicates each): `se` of `log(nu - 1)` ran from 4878 to 1.2e4
#' on the 35 fits whose `nu` ran off, and from 0.35 to 22 on the 49
#' where it did not. `|log(nu - 1)|` stayed below 5.1 on every one of
#' those, so both conditions together fired on none of them.
#'
#' The condition reads `abs(est)` and not `est`, because the same
#' parameter runs off both ways. On Cauchy errors `nu` goes to its
#' floor instead: 8 of 12 fits at `n` of 60 and 200 returned
#' `log(nu - 1)` near -19 with a standard error of 4000 to 8500, which
#' is the same defect and the same evidence. The message names both
#' directions, because `gaussian()` is the right answer to only one.
#'
#' `negbinomial()`'s `shape` does it too, on data with no
#' overdispersion: 4 of 4 Poisson-generated fits returned a log shape
#' near 18 with a non-finite standard error, and the check named them.
#'
#' Primary dpars are left to `diagnose_separation()`. A runaway `mu` is
#' separation in a binomial and collinearity elsewhere, and this
#' message describes neither.
#'
#' @noRd
diagnose_unbounded_dpar <- function(fit, ps) {
  rows <- list()
  for (lp in fit$frame[["linpreds"]]) {
    fam <- fit$spec$responses[[lp[["resp"]]]]$family
    if (lp[["dpar"]] %in% (fam[["primary_dpars"]] %||% "mu")) next
    if (is.null(lp[["X"]]) || !ncol(lp[["X"]])) next
    est <- ps$est[[lp[["par"]]]][lp[["idx"]]]
    se <- ps$se[[lp[["par"]]]][lp[["idx"]]]
    # the pair, not either half: a dpar legitimately far out on its
    # link (a log sd of a column scaled by 1e6) keeps a small se
    hit <- which(abs(est) > 10 & (!is.finite(se) | se > abs(est)))
    lk <- fam[["links"]][[lp[["dpar"]]]]
    for (i in hit) {
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = paste0(coef_block_key(fit, lp), ": ",
                           colnames(lp[["X"]])[i]),
        estimate = est[i], std.error = se[i],
        # the natural scale is what makes the message readable: a nu
        # of 2e9 says "gaussian" where a log(nu - 1) of 21.5 does not
        value = tryCatch(as.numeric(lk$linkinv(est[i])),
                         error = function(e) NA_real_)
      )
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Predictor columns whose spread is many orders of magnitude away from
#' one: the objective's curvature then spans the same range, and the
#' optimizer's convergence tolerances are absolute. `autoscale = TRUE`
#' fixes it without touching the model. `[glmmTMB diagnose()]`
#'
#' @noRd
diagnose_predictor_scale <- function(fit, tol = 3) {
  rows <- list()
  for (lp in fit$frame[["linpreds"]]) {
    if (is.null(lp[["X"]]) || !ncol(lp[["X"]])) next
    X <- as.matrix(lp[["X"]])
    for (j in seq_len(ncol(X))) {
      if (identical(colnames(X)[j], "(Intercept)")) next
      s <- stats::sd(X[, j])
      if (!is.finite(s) || s <= 0) next
      if (abs(log10(s)) <= tol) next
      rows[[length(rows) + 1L]] <- data.frame(
        column = paste0(coef_block_key(fit, lp), ": ", colnames(X)[j]),
        sd = s
      )
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Outer parameters the likelihood does not depend on at all: zero
#' gradient AND an empty Hessian row.
#'
#' A singular Hessian has two very different causes and the same
#' symptom. Two parameters that trade off against each other are
#' over-parameterization, which is what the standard-error warning has
#' always said. A parameter the likelihood is FLAT in is something else:
#' a nonlinear term that has left its own support, where the remedy is a
#' starting value, not a smaller model. A Gaussian peak
#' `exp(lamp) * exp(-0.5 * ((w - pk) / exp(lsig))^2)` with `pk` at the
#' default 0 and `w` running from 1 to 45 underflows to zero everywhere,
#' and `lamp`, `pk`, `lsig` and `pk`'s group standard deviation then
#' have no gradient at all. The old report named every parameter as
#' having a bad standard error, which names none of them.
#'
#' Measured, not assumed. `fit$obj$he()` is unavailable on a model with
#' random effects ("Hessian not yet implemented"), so the row is read
#' from the gradient: perturb one parameter and see whether the WHOLE
#' gradient vector moves. A direction the likelihood is flat in moves it
#' by nothing. On the peak fit above that separates the four flat
#' parameters from the five identified ones by ten orders of magnitude
#' (row norms 3e-30 to 3e-10 against 20 to 700) in 0.08 s for nine
#' gradient evaluations.
#'
#' Only parameters with a vanishing gradient are candidates, so a
#' healthy fit pays for nothing, and the caller gates the whole check on
#' a covariance that already failed.
#'
#' @noRd
diagnose_flat <- function(fit, gr = NULL, tol = 1e-8, gtol = 1e-6) {
  p <- fit$opt$par
  if (!length(p)) return(character(0))
  if (is.null(gr)) {
    gr <- tryCatch(drop(fit$obj$gr(p)), error = function(e) NULL)
  }
  if (is.null(gr) || length(gr) != length(p) || anyNA(gr)) {
    return(character(0))
  }
  cand <- which(abs(gr) <= gtol)
  if (!length(cand)) return(character(0))
  flat <- logical(length(cand))
  for (k in seq_along(cand)) {
    j <- cand[k]
    # a step scaled to the parameter, so a coefficient at 1e4 is
    # perturbed meaningfully and one at zero is perturbed at all
    h <- max(1e-4, 1e-4 * abs(p[j]))
    pj <- p
    pj[j] <- pj[j] + h
    g1 <- tryCatch(drop(fit$obj$gr(pj)), error = function(e) NULL)
    if (is.null(g1) || length(g1) != length(gr) || anyNA(g1)) next
    flat[k] <- max(abs(g1 - gr)) / h <= tol
  }
  if (!any(flat)) return(character(0))
  outer_par_names(fit)[cand[flat]]
}

#' `diagnose_flat()` once per fit. vcov(), summary(), confint() and
#' diagnose() all reach the same verdict on the same degenerate object.
#'
#' @noRd
flat_pars <- function(fit) {
  cache <- fit$cache
  if (is.environment(cache) && !is.null(cache$flat_pars)) {
    return(cache$flat_pars)
  }
  out <- tryCatch(diagnose_flat(fit), error = function(e) character(0))
  if (is.environment(cache)) cache$flat_pars <- out
  out
}

#' Does this fit carry nonlinear parameters at all?
#'
#' `diagnose_flat()` measures flatness and nothing else, so it finds
#' flat directions in models that have no nonlinear term: a saturated
#' `mo()` simplex is one (its softmax sits in a corner, and moving the
#' last `zeta` leaves the objective unchanged to fifteen digits). The
#' nonlinear EXPLANATION must therefore be gated on the model the
#' explanation is about, or the report asserts a cause it never
#' measured - which is the defect the flat-direction check was added to
#' remove, relocated to another model class.
#'
#' @noRd
fit_has_nlpars <- function(fit) {
  resp <- fit$frame[["spec"]]$responses %||% fit$spec$responses %||% list()
  any(vapply(resp, function(r) length(r$nlpars %||% character(0)) > 0L, TRUE))
}

#' A flat parameter's name, with the block its standard deviation
#' belongs to when the name is a bare `theta_k`.
#'
#' `theta_3` in a flat set is the standard deviation OF a flat
#' parameter's random effect, and a bare index leaves the reader to work
#' that out. `log_sd_theta_index()` already carries the labels.
#'
#' @noRd
flat_par_display <- function(fit, nms) {
  sd_i <- tryCatch(log_sd_theta_index(fit), error = function(e) integer(0))
  if (!length(sd_i)) return(nms)
  key <- paste0("theta_", as.integer(sd_i))
  vapply(nms, function(n) {
    j <- match(n, key)
    if (is.na(j)) n else paste0(n, " (sd of ", names(sd_i)[j], ")")
  }, "", USE.NAMES = FALSE)
}

#' The sentence that names the flat directions, or "" when there are
#' none. Appended to whichever warning reports the failed covariance.
#'
#' @noRd
flat_par_note <- function(fit) {
  fl <- flat_pars(fit)
  if (!length(fl)) return("")
  one <- length(fl) == 1L
  paste0(". The likelihood is FLAT in ", length(fl),
         if (one) " direction: " else " directions: ",
         paste(flat_par_display(fit, fl), collapse = ", "),
         " - zero gradient and an empty Hessian row, so ",
         if (one) "that parameter is" else "those parameters are",
         " not identified AT THIS POINT rather than over-parameterized",
         if (fit_has_nlpars(fit)) {
           paste0(". A nonlinear term that has left its own support does ",
                  "this; give it a starting value that puts it back (see ",
                  "par_template())")
         } else {
           paste0(": moving ", if (one) "it" else "them",
                  " does not change the likelihood at all")
         })
}

#' lme4's isSingular: a variance component sitting on the boundary of its
#' parameter space - a standard deviation at zero, or a correlation at
#' +/-1. The verdict is read off the estimates alone, so it stands even
#' when the Hessian is positive definite and every gradient is tiny (a
#' collapsed component is a well-behaved optimum of a model the data
#' cannot support). `[lme4 test-isSingular.R, #660]`
#'
#' @noRd
diagnose_singular <- function(fit, tol = 1e-4) {
  if (!length(fit$frame[["re_blocks"]])) return(NULL)
  vc <- tryCatch(as.data.frame(varcorr_matrices(fit)),
                 error = function(e) NULL)
  if (is.null(vc) || !nrow(vc)) return(NULL)
  is_cor <- !is.na(vc$var2)
  bad <- ifelse(is_cor, abs(vc$sdcor) > 1 - tol, vc$sdcor < tol)
  bad[is.na(bad)] <- FALSE
  if (!any(bad)) return(NULL)
  out <- data.frame(
    block = vc$grp[bad],
    term = ifelse(is_cor[bad],
                  paste0("cor(", vc$var1[bad], ",", vc$var2[bad], ")"),
                  paste0("sd(", vc$var1[bad], ")")),
    value = vc$sdcor[bad]
  )
  rownames(out) <- NULL
  out
}

#' The theta components that are LOG STANDARD DEVIATIONS, named the way
#' confint_varcorr() names their rows.
#'
#' The `|theta| > 8` near-singularity heuristic only reads as a boundary
#' fit on a log sd: `e^-8` is a variance no data supports. The other
#' components live on their own scales, where the same magnitude is
#' ordinary - a CAR/BYM2 mixing proportion and an AR(1) phi are logit-
#' and arctan-like, so a rho legitimately at the boundary sits at
#' `|theta| >> 8`, and the SPDE's (log tau, log kappa) are a precision
#' and an inverse range, neither of which is a standard deviation.
#' Reading those as singular fits is a false alarm on a converged model.
#' Structures whose registry declares no sd (rr, equalto, spde)
#' contribute nothing.
#'
#' @noRd
log_sd_theta_index <- function(fit) {
  idx <- integer(0)
  nms <- character(0)
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    # blocks whose confint() rows are hand-written carry their own name
    shared <- switch(bk[["covstruct"]],
                     smooth = "sd(wiggle)", gp = "sd(gp)",
                     hsgp = "sd(gp)", car = "sd(car)",
                     spde = NA_character_, equalto = NA_character_,
                     NULL)
    if (!is.null(shared)) {
      if (is.na(shared)) next
      idx <- c(idx, bk[["theta_idx"]][1L])
      nms <- c(nms, paste0(bk[["term_label"]], " ", shared))
      next
    }
    reg <- covstruct_registry[[bk[["covstruct"]]]]
    si <- if (is.null(reg)) integer(0) else {
      tryCatch(as.integer(reg$sd_idx(bk[["dim"]])),
               error = function(e) integer(0))
    }
    for (i in seq_along(si)) {
      idx <- c(idx, bk[["theta_idx"]][si[i]])
      # the generic confint() path labels sd rows by column name, and
      # falls back to the first when one sd is shared across columns
      nms <- c(nms, paste0(bk[["term_label"]], " ",
                           bk[["cnms"]][min(i, length(bk[["cnms"]]))]))
    }
  }
  th_n <- length(fit$estimates[["theta"]] %||% numeric(0))
  keep <- !is.na(idx) & idx >= 1L & idx <= th_n
  stats::setNames(idx[keep], nms[keep])
}

#' Convergence diagnostics for a frmtmb fit
#'
#' Reports the optimizer's own verdict plus six checks that a converged
#' fit can still fail: non-finite standard errors, flat directions,
#' complete separation in a binomial-type fit, a distributional
#' parameter whose maximum likelihood is outside its own range,
#' predictor columns scaled far from one, and variance components on
#' the boundary of their parameter space (lme4's `isSingular()`, read
#' off the estimates rather than the Hessian).
#'
#' NON-FINITE TRIALS are points where the optimizer found the objective
#' undefined and stepped back, for example a line search that crosses
#' below zero on the `1/mu^2` link. They do not make a fit wrong, so they
#' are reported as a count and do not count against "No convergence
#' problems detected". A large count on a model with no restricted link
#' can point to a density that is undefined where it should not be.
#'
#' A DISTRIBUTIONAL PARAMETER AT THE END OF ITS LINK is one whose
#' estimate is far out on the link scale AND whose standard error is
#' larger than the estimate itself. The likelihood was still rising
#' where the optimizer stopped, so the number reported is the stopping
#' point and not an estimate. `student()`'s `nu` does this on any data
#' with no heavy tails, because a student-t reaches the gaussian only
#' as `nu` goes to infinity: the maximum is never attained, and two
#' runs of the same model can report degrees of freedom orders of
#' magnitude apart while agreeing on every coefficient. Refit with
#' `gaussian()`, or hold `nu` somewhere finite with [set_prior()].
#' The same parameter runs the other way, down to one, on data whose
#' tails are heavier than any identified `nu` can hold. The check
#' names that too, with a large NEGATIVE estimate and a natural-scale
#' value of one, and there `gaussian()` is the wrong answer: the data
#' is the message, and a prior is the way to hold `nu` finite.
#'
#' A FLAT DIRECTION is an outer parameter the likelihood does not depend
#' on: zero gradient and an empty Hessian row. It separates the two
#' causes of `NaN` standard errors. Parameters that trade off against
#' each other are over-parameterization, and the model is too big.
#' Parameters the likelihood is flat in are unidentified AT THIS POINT,
#' and the remedy is a starting value: a nonlinear term evaluated
#' outside its own support (a bump whose centre starts far from the
#' data) is flat in several of its parameters at once. One unusable
#' direction makes EVERY standard error `NaN`, so `bad_se` names the
#' whole vector and `flat` names the cause. The check is measured by
#' perturbing each candidate and seeing whether the gradient moves, and
#' runs only when the covariance has already failed.
#'
#' @param fit A `frmtmb_fit`.
#' @param quiet If `TRUE`, return the diagnostics without printing.
#' @return Invisibly, a list of diagnostics.
#'
#' @srrstats {RE2.4b} Perfect collinearity between the predictors and the
#'   response is reported as complete separation: a binomial-type fit
#'   whose coefficients diverge because a predictor perfectly predicts
#'   the response is flagged by name, with the offending estimate, when
#'   [diagnose()] is called on the fit. The test suite checks both that a
#'   separating design is flagged and that a well-behaved binomial fit is
#'   not. For a continuous response an exact linear relationship is not a
#'   degenerate fit but a zero-residual one, so it is left to the
#'   dispersion estimate rather than reported as collinearity.
#' @srrstats {RE4.7} Convergence statistics are available from the model
#'   object. `fit$opt$convergence` and `fit$opt$message` carry the
#'   optimizer's verdict, and `diagnose()` returns the maximum absolute
#'   gradient, the worst-offending parameter, the positive-definiteness
#'   of the Hessian, non-finite standard errors, the smallest eigenvalue
#'   of the covariance, the flat directions, boundary (singular)
#'   variance components, separation, distributional parameters at the
#'   end of their link, and predictor scaling. `frm_allfit()` refits across
#'   optimizers as a further convergence check, and `check_laplace()`
#'   audits the approximation itself.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' diagnose(fit)
#'
#' # a random effect the data cannot support collapses to the boundary,
#' # which is a valid fit but a warning about the model
#' dd$h <- factor(rep(1:10, each = 10))
#' fit_s <- frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), data = dd)
#' d <- diagnose(fit_s, quiet = TRUE)
#' d$singular
#'
#' # a predictor scaled far from one slows the optimizer down; the
#' # remedy is frmtmb_control(autoscale = TRUE)
#' dd$xbig <- dd$x * 1e5
#' diagnose(frm(bf(xbig ~ 1) + gaussian(), data = dd), quiet = TRUE)$scale
#' @export
diagnose <- function(fit, quiet = FALSE) {
  if (!inherits(fit, "frmtmb_fit")) {
    frm_stop("diagnose() needs a model fitted by frm(), not ",
             arg_desc(fit), call. = FALSE)
  }
  check_flag(quiet, "quiet")
  nm <- outer_par_names(fit)
  # a degenerate fit (no free outer parameters) has no gradient, no
  # covariance and no theta to report on
  degenerate <- !length(fit$opt$par)
  gr <- if (degenerate) numeric(0) else drop(fit$obj$gr(fit$opt$par))
  V <- sdr_of(fit)$cov.fixed
  # on the pathological fits diagnose() exists for, cov.fixed can carry
  # negative diagonal entries; the resulting NaN SEs are the finding
  # (reported through bad_se), not a warning to relay
  se <- if (!length(V)) numeric(0) else suppressWarnings(sqrt(diag(V)))
  ev <- if (!length(V)) NULL else {
    tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
             error = function(e) NULL)
  }
  # theta is absent from fits with no random effects; abs(NULL) is an
  # error, not an empty result
  th <- fit$estimates[["theta"]] %||% numeric(0)
  # only the log-sd components; see log_sd_theta_index()
  sd_i <- log_sd_theta_index(fit)
  ps <- tryCatch(suppressWarnings(par_est_se(fit)),
                 error = function(e) NULL)
  out <- list(
    convergence = fit$opt$convergence,
    message = fit$opt$message,
    # NULL for an optimizer that does not count them (optim, a custom one)
    nonfinite_trials = fit$opt$nonfinite_trials,
    max_grad =if (length(gr)) max(abs(gr)) else NA_real_,
    worst_grad = if (length(gr)) nm[which.max(abs(gr))] else NA_character_,
    pdHess = isTRUE(sdr_of(fit)$pdHess),
    bad_se = nm[!is.finite(se)],
    # gated on a covariance that already failed: on a healthy fit every
    # parameter would be a candidate (a converged gradient is tiny
    # everywhere) and the check would cost one gradient evaluation per
    # parameter for a finding that cannot arise
    flat = if (!isTRUE(sdr_of(fit)$pdHess) || any(!is.finite(se))) {
      flat_pars(fit)
    } else character(0),
    min_cov_eigenvalue = if (!is.null(ev) && length(ev)) min(ev),
    extreme_theta = sd_i[abs(th[sd_i]) > 8],
    separation = if (!is.null(ps)) diagnose_separation(fit, ps),
    unbounded_dpar = if (!is.null(ps)) diagnose_unbounded_dpar(fit, ps),
    predictor_scale = diagnose_predictor_scale(fit),
    singular = diagnose_singular(fit)
  )
  if (!quiet) {
    cat("Optimizer convergence code:", out$convergence,
        if (!is.null(out$message)) paste0("(", out$message, ")"), "\n")
    if (degenerate) {
      cat("No free parameters: the model is degenerate and the ",
          "likelihood was evaluated once\n", sep = "")
    } else {
      cat("Max |gradient|:", format(out$max_grad, digits = 4),
          "at", out$worst_grad, "\n")
    }
    cat("Hessian positive definite:", out$pdHess, "\n")
    # the optimizer no longer warns for these (see nlminb_trial_fn()),
    # so this line is where a user reading diagnostics sees them
    if (isTRUE(out$nonfinite_trials > 0L)) {
      one <- out$nonfinite_trials == 1L
      cat("Non-finite objective at ", out$nonfinite_trials, " trial ",
          if (one) "point" else "points",
          "; the optimizer stepped back from ", if (one) "it" else "them",
          "\n", sep = "")
    }
    if (length(out$bad_se)) {
      cat("Non-finite standard errors:",
          paste(out$bad_se, collapse = ", "), "\n")
    }
    if (length(out$flat)) {
      # named BEFORE the covariance verdict is read as
      # over-parameterization: these are the parameters the likelihood
      # does not depend on, and the others inherit their NaN
      one <- length(out$flat) == 1L
      cat("Flat directions (zero gradient, empty Hessian row): ",
          paste(flat_par_display(fit, out$flat), collapse = ", "),
          "\n  The likelihood does not depend on ",
          if (one) "this parameter" else "these parameters",
          " at this point, so the model is UNIDENTIFIED HERE rather ",
          "than overparameterized, and every other standard error is ",
          "NaN because the Hessian cannot be inverted. ",
          # the CAUSE is gated on the model: the check measures flatness
          # and nothing else, and finds it in models with no nonlinear
          # term at all (a saturated mo() simplex, say). Asserting a bump
          # there is a false lead, and a false lead is worse than the
          # vague message this block replaced.
          if (fit_has_nlpars(fit)) {
            paste0("A nonlinear term evaluated outside its own support ",
                   "does this - a bump whose centre starts far from the ",
                   "data is flat in its centre, amplitude and width at ",
                   "once. Give those parameters a starting value that ",
                   "puts the term back on the data: par_template() names ",
                   "them, and start = takes the names straight back.\n")
          } else {
            paste0("Moving ", if (one) "it" else "them", " does not ",
                   "change the likelihood at all, so no amount of ",
                   "optimization will pin ", if (one) "it" else "them",
                   " down: what the fit reports for ",
                   if (one) "that parameter" else "those parameters",
                   " is wherever the optimizer stopped. par_template() ",
                   if (one) "names it." else "names them.", "\n")
          }, sep = "")
    }
    if (length(out$extreme_theta)) {
      cat("Extreme covariance parameters (|log sd| > 8): ",
          paste(paste0(names(out$extreme_theta), " (log sd ",
                       format(th[out$extreme_theta], digits = 3), ")"),
                collapse = "; "),
          "\n  The fit is near-singular; consider simplifying the ",
          "random effects (e.g. diag() instead of a correlated term)\n",
          sep = "")
    }
    if (!is.null(out$singular)) {
      cat("Singular fit: ",
          paste(paste0(out$singular$block, " ", out$singular$term,
                       " = ", format(out$singular$value, digits = 3)),
                collapse = "; "),
          "\n  A variance component is on the boundary of its ",
          "parameter space. The fit is valid but the random-effect ",
          "structure is more complex than the data support; drop the ",
          "collapsed term or use diag() instead of a correlated ",
          "block.\n", sep = "")
    }
    if (!is.null(out$separation)) {
      cat("Likely complete separation: ",
          paste(paste0(out$separation$parameter, " = ",
                       format(out$separation$estimate, digits = 3),
                       " (se ", format(out$separation$std.error,
                                       digits = 3), ")"),
                collapse = "; "),
          "\n  Coefficients this large on the link scale with standard ",
          "errors to match mean the maximum likelihood is at infinity: ",
          "some combination of the predictors separates the outcome ",
          "perfectly. Drop or pool the offending predictor, or add a ",
          "prior (see set_prior()).\n", sep = "")
    }
    if (!is.null(out$unbounded_dpar)) {
      cat("Distributional parameter at the end of its link: ",
          paste(paste0(out$unbounded_dpar$parameter, " = ",
                       format(out$unbounded_dpar$estimate, digits = 3),
                       " (se ", format(out$unbounded_dpar$std.error,
                                       digits = 3),
                       ", natural scale ",
                       format(out$unbounded_dpar$value, digits = 3), ")"),
                collapse = "; "),
          "\n  The likelihood was still rising when the optimizer ",
          "stopped, so this is where it stopped and not an estimate; ",
          "its standard error is larger than the estimate itself. Read ",
          "the natural-scale value as the end of a range and not as a ",
          "number. student()'s nu runs UP on data with no heavy tails, ",
          "and the fit is then the gaussian one, so refit with ",
          "gaussian() to say so; it runs DOWN to one on tails heavier ",
          "than any identified nu can hold, and there the data is the ",
          "message. Either way a prior holds the parameter somewhere ",
          "finite (see set_prior()).\n", sep = "")
    }
    if (!is.null(out$predictor_scale)) {
      cat("Badly scaled predictors: ",
          paste(paste0(out$predictor_scale$column, " (sd ",
                       format(out$predictor_scale$sd, digits = 3), ")"),
                collapse = "; "),
          if (!is.null(fit$par_units)) {
            # this fit already ran the standardized pre-fit, by default
            # or on request, so "refit with autoscale = TRUE" repeats it
            paste0("\n  The fit was already standardized internally ",
                   "(frmtmb_control(autoscale = )); rescaling the column ",
                   "is still the cleaner model.\n")
          } else {
            paste0("\n  Rescale the column, or refit with ",
                   "frmtmb_control(autoscale = TRUE).\n")
          }, sep = "")
    }
    clean <- out$convergence == 0 && out$pdHess && !length(out$bad_se) &&
      is.null(out$singular) && is.null(out$separation) &&
      is.null(out$unbounded_dpar) && is.null(out$predictor_scale) &&
      (degenerate || out$max_grad < 1e-3)
    if (clean) cat("No convergence problems detected\n")
  }
  invisible(out)
}

#' Largest absolute entry of a matrix, and 0 for an empty one.
#'
#' @noRd
maxabs <- function(M) if (!length(M)) 0 else max(abs(M))

#' Residual of B after projecting onto the column space of A.
#'
#' @noRd
proj_resid <- function(A, B) if (!ncol(A)) B else qr.resid(qr(A), B)

#' Test whether two design matrices span the same column space, up to a
#' relative tolerance. Used to decide when two REML likelihoods are
#' comparable.
#'
#' @noRd
same_column_space <- function(A, B, tol = 1e-8) {
  if (nrow(A) != nrow(B)) return(FALSE)
  if (!ncol(A) && !ncol(B)) return(TRUE)
  s <- max(1, maxabs(A), maxabs(B))
  if (!ncol(A) || !ncol(B)) return(FALSE)
  maxabs(proj_resid(A, B)) <= tol * s &&
    maxabs(proj_resid(B, A)) <= tol * s
}

#' Designs REML integrates out: the primary-dpar linear predictors, whose
#' coefficients live in the `beta` template component (dpar formulas keep
#' their coefficients in `betad` and stay outer).
#'
#' @noRd
reml_designs <- function(fit) {
  parts <- list()
  for (lp in fit$frame[["linpreds"]]) {
    if (!identical(lp[["par"]], "beta")) next
    X <- if (is.null(lp[["X"]])) {
      matrix(numeric(0), fit$frame[["n_obs"]], 0L)
    } else {
      as.matrix(lp[["X"]])
    }
    parts[[linpred_key(lp[["resp"]], lp[["dpar"]])]] <- X
  }
  parts[order(names(parts))]
}

#' A REML likelihood carries a `-1/2 log|X' V^-1 X|` term, so it is a
#' likelihood for a DIFFERENT quantity - the error contrasts - once X
#' changes, and differencing two of them is meaningless. It is perfectly
#' meaningful when the error contrasts are the same, which is exactly
#' when the fixed-effect designs span the same column space; that is the
#' usual REML comparison of variance-component structures. Refusing every
#' REML fit (the old behavior) refused that case too. `[glmmTMB#776]`
#'
#' @noRd
reml_comparable <- function(fits) {
  d1 <- reml_designs(fits[[1]])
  for (f in fits[-1]) {
    d2 <- reml_designs(f)
    if (!identical(names(d1), names(d2))) return(FALSE)
    for (k in names(d1)) {
      if (!same_column_space(d1[[k]], d2[[k]])) return(FALSE)
    }
  }
  TRUE
}

#' Every fixed-effect design of a fit, keyed by response and dpar.
#' `reml_designs()` keeps only the ones REML integrates out; nesting is
#' a question about all of them, a dpar's own predictor included.
#'
#' @noRd
fixef_designs <- function(fit) {
  parts <- list()
  for (lp in fit$frame[["linpreds"]]) {
    X <- if (is.null(lp[["X"]])) {
      matrix(numeric(0), fit$frame[["n_obs"]], 0L)
    } else {
      as.matrix(lp[["X"]])
    }
    parts[[linpred_key(lp[["resp"]], lp[["dpar"]])]] <- X
  }
  parts[order(names(parts))]
}

#' Every coefficient a fit estimates for its FIXED effects, by name.
#' The `beta` and `betad` components together, because a dpar's
#' predictor is a fixed effect too and dropping a term from it is
#' exactly the comparison this has to see.
#'
#' @noRd
fixef_coef_names <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  c(names(tpl[["beta"]]), names(tpl[["betad"]]))
}

#' Warn when two models being compared by likelihood ratio have
#' fixed-effect column sets neither of which contains the other.
#'
#' The LRT's null distribution needs the smaller model to be a
#' restriction of the larger. Deciding that in general is not possible
#' here - a nonlinear reparameterization is nesting that no comparison
#' of names or column spaces can see, and equality constraints across
#' dpars are not in the design at all - so this makes the cheap check
#' and says only what it checked. Names first, because renaming is what
#' a genuinely different term does; column space second, so that a
#' change of basis for the same span (`poly(x, 2)` against an
#' orthogonalized pair of columns) is not reported as a different
#' model.
#'
#' @noRd
warn_non_nested <- function(fits) {
  for (i in seq_len(length(fits) - 1L)) {
    a <- fixef_coef_names(fits[[i]])
    b <- fixef_coef_names(fits[[i + 1L]])
    if (all(a %in% b) || all(b %in% a)) next
    da <- fixef_designs(fits[[i]])
    db <- fixef_designs(fits[[i + 1L]])
    if (identical(names(da), names(db)) &&
        all(vapply(names(da),
                   function(k) column_space_within(da[[k]], db[[k]]) ||
                     column_space_within(db[[k]], da[[k]]),
                   TRUE))) {
      next
    }
    frm_warning("anova(): the fixed effects of ", model_label(fits[[i]]),
                " and ", model_label(fits[[i + 1L]]),
                " are not nested - neither model's coefficients are a ",
                "subset of the other's, and neither design sits inside the ",
                "other's column space. A likelihood-ratio test between ",
                "models that are not nested has no chi-square null ",
                "distribution; compare them by AIC instead",
                call. = FALSE)
  }
  invisible(NULL)
}

#' Whether every column of `A` lies in the span of `B`. One half of
#' `same_column_space()`, which is what nesting needs: containment, not
#' equality.
#'
#' @noRd
column_space_within <- function(A, B, tol = 1e-8) {
  if (nrow(A) != nrow(B)) return(FALSE)
  if (!ncol(A)) return(TRUE)
  if (!ncol(B)) return(FALSE)
  s <- max(1, maxabs(A), maxabs(B))
  maxabs(proj_resid(A, B)) <= tol * s
}

#' The cheapest correct REML -> ML conversion: reuse the assembled
#' design (no formula parsing, no frame assembly) and warm-start the
#' optimizer at the REML estimates. The parameter template is the same
#' list either way; REML only decides whether `beta` is integrated out,
#' so the REML estimates are a valid ML start.
#'
#' @noRd
anova_refit_ml <- function(fit) {
  ctl <- fit$control %||% frmtmb_control()
  # one anova() call can trigger several refits; a verbose original fit
  # must not make each of them narrate itself
  ctl$verbose <- FALSE
  fit_assembled(fit$spec, fit$frame, fit$bform, fit$call,
                REML = FALSE, start = NULL, control = ctl, se = FALSE,
                lower = fit$lower, upper = fit$upper,
                prior = fit$prior,
                quadrature = isTRUE(fit$quadrature),
                importance = fit$importance$draws %||% 0L,
                template = fit$estimates,
                data2 = fit$data2 %||% list())
}

#' Chi-square tail probability for a likelihood-ratio statistic, with a
#' zero degree-of-freedom difference reported as NA.
#'
#' `pchisq(0, df = 0, lower.tail = FALSE)` is 0, so two models of the
#' same dimension - a reparameterization, or the same model passed
#' twice - used to print "< 2.2e-16 ***" for a test that was never run.
#' A chi-square with no degrees of freedom is a point mass at zero and
#' has no p-value; NA says so, and drops the significance stars with it.
#'
#' @noRd
lrt_pvalue <- function(chisq, ddf) {
  p <- stats::pchisq(chisq, ddf, lower.tail = FALSE)
  p[!is.na(ddf) & ddf <= 0] <- NA_real_
  p
}

#' Likelihood-ratio tests between nested frmtmb fits
#'
#' ML fits compare freely. REML fits compare only with each other, and
#' only when their fixed-effect designs span the same column space: a
#' REML likelihood is a likelihood for the error contrasts of that
#' design, so two of them are on a common scale exactly when the design
#' is the same. That covers the usual REML use - testing
#' variance-component structures with the fixed effects held fixed - and
#' refuses the rest with the reason (glmmTMB#776).
#'
#' `refit = TRUE` is the lme4 convenience for the refused case: every
#' REML fit in the comparison is refit with `REML = FALSE` and the ML
#' fits are compared instead. lme4 does this silently by default; here
#' it is opt-in and the message names the models that were refit.
#'
#' When the smaller model removes a variance component, the null value
#' sits on the boundary of the parameter space and the usual chi-square
#' reference is wrong: the asymptotic null is a mixture (for one
#' component, half a point mass at zero and half a chi-square with one
#' df), so the reported p-value is conservative - up to a factor of two
#' for a single component. lme4 and glmmTMB report the same naive
#' p-value; halve it for the one-component case, or use
#' [frm_bootstrap()] for a simulation-based reference.
#'
#' @section Nesting is assumed, not verified:
#' A likelihood-ratio statistic has a chi-square null distribution only
#' when the smaller model is a restriction of the larger. `anova()`
#' cannot verify that in general: nesting through a nonlinear
#' reparameterization, or through a constraint that ties parameters
#' across distributional parameters, is invisible to anything the
#' fitted objects carry. What it does check is cheap and stated: if
#' neither model's fixed-effect coefficient names are a subset of the
#' other's, and neither fixed-effect design sits inside the other's
#' column space, it warns. A comparison that passes that check is not
#' thereby verified to be nested. Two models that are genuinely not
#' nested are compared by AIC, not by this table.
#'
#' @param object A `frmtmb_fit`.
#' @param ... Further `frmtmb_fit` objects, nested with `object`.
#' @param refit If `TRUE`, refit every REML fit in the comparison with
#'   ML and compare those, with a message naming what was refit. The
#'   refits reuse the assembled design and warm-start at the REML
#'   estimates. `FALSE` (the default) keeps the REML fits and refuses
#'   the comparisons a restricted likelihood cannot support.
#' @return An `anova` table.
#'
#' @srrstats {RE4.11} Goodness-of-fit statistics are available for the
#'   fitted model. `logLik()` reports the log-likelihood with its degrees
#'   of freedom and `nobs`, so `AIC()` and `BIC()` work through the
#'   `stats` defaults; `extractAIC()` and `deviance()` are implemented;
#'   `anova()` gives likelihood-ratio tests between nested fits and
#'   `drop1()` single-term deletions. Effect sizes with the coefficients
#'   come from `summary()` (estimate, standard error, z, p) and
#'   `confint()`. The boundary problem for variance-component tests is
#'   documented above rather than left implicit.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' dd <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:20, 10)))
#' u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.5))
#' dd$y <- rnorm(n, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
#'
#' m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' m1 <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
#' anova(m0, m1)
#'
#' # dropping a variance component puts the null on the boundary, so
#' # this p-value is conservative by up to a factor of two
#' m2 <- frm(bf(y ~ x) + gaussian(), data = dd)
#' anova(m2, m0)
#'
#' # REML fits compare only when the fixed-effect designs agree, which
#' # is the case for a variance-component test
#' r0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = TRUE)
#' r1 <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd, REML = TRUE)
#' anova(r0, r1)
#' # differing designs are refused; refit = TRUE compares ML fits instead
#' rz <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, REML = TRUE)
#' try(anova(r0, rz))
#' anova(r0, rz, refit = TRUE)
#' @export
anova.frmtmb_fit <- function(object, ..., refit = FALSE) {
  fits <- c(list(object), Filter(function(x) inherits(x, "frmtmb_fit"),
                                 list(...)))
  if (length(fits) < 2) {
    frm_stop("anova() needs at least two frmtmb fits to compare", call. = FALSE)
  }
  reml <- vapply(fits, `[[`, TRUE, "REML")
  if (any(reml) && isTRUE(refit)) {
    labs <- vapply(fits[reml], model_label, "")
    frm_message("anova(): refitting ", length(labs), " REML model",
                if (length(labs) != 1L) "s" else "", " with ML: ",
                paste(labs, collapse = "; "))
    fits[reml] <- lapply(fits[reml], anova_refit_ml)
    reml[] <- FALSE
  }
  if (any(reml)) {
    if (!all(reml)) {
      frm_stop("anova() cannot mix REML and ML fits: their likelihoods are ",
               "for different quantities. Refit them all with the same ",
               "REML setting, or pass refit = TRUE to compare them as ML ",
               "fits", call. = FALSE)
    }
    if (!reml_comparable(fits)) {
      frm_stop("REML likelihoods are comparable only between fits whose ",
               "fixed-effect designs span the same column space; these do ",
               "not. Pass refit = TRUE (or refit with REML = FALSE) to ",
               "compare fixed effects, or hold the fixed effects fixed to ",
               "compare random-effect structures", call. = FALSE)
    }
  }
  # Likelihoods computed on different data are not on a common scale, so
  # the LRT would be meaningless (and can come out negative). lme4 keys
  # its equivalent check off the `data` argument in the call, which both
  # false-positives on identical frames and misses NA-dropped rows;
  # comparing the response actually used catches the real cases.
  # [lme4#622]
  nobs_all <- vapply(fits, function(f) as.integer(f$frame[["n_obs"]]), 0L)
  if (length(unique(nobs_all)) > 1L) {
    frm_stop("anova() needs fits with the same number of observations (got ",
             paste(unique(nobs_all), collapse = ", "),
             "); models fit to different data or with different NA rows ",
             "dropped are not comparable", call. = FALSE)
  }
  ll <- vapply(fits, function(f) as.numeric(logLik(f)), 0)
  df <- vapply(fits, function(f) attr(logLik(f), "df"), 0L)
  ord <- order(df)
  fits <- fits[ord]; ll <- ll[ord]; df <- df[ord]
  warn_non_nested(fits)
  chisq <- c(NA, 2 * diff(ll))
  ddf <- c(NA, diff(df))
  p <- lrt_pvalue(chisq, ddf)
  tab <- data.frame(
    Df = df, logLik = ll, AIC = -2 * ll + 2 * df,
    Chisq = chisq, `Chi Df` = ddf, `Pr(>Chisq)` = p,
    check.names = FALSE
  )
  rownames(tab) <- make.unique(vapply(fits, model_label, ""))
  structure(tab, class = c("anova", "data.frame"),
            heading = paste0("Likelihood-ratio tests\n",
                             "Each test assumes the smaller model is ",
                             "nested in the larger; see ",
                             "?anova.frmtmb_fit\n"))
}

#' Single-term deletions
#'
#' Drops each fixed-effect term of the primary (`mu`) formula in turn,
#' refits, and tabulates AIC (and likelihood-ratio tests with
#' `test = "Chisq"`), following [stats::drop1()] and lme4's
#' `drop1.merMod`. Random-effect, smooth, and `mo()`/`mi()` terms are
#' not part of the deletion scope.
#'
#' @param object A `frmtmb_fit` from an ML fit (`REML = FALSE`) of a
#'   univariate model.
#' @param scope Terms to drop: a character vector or a right-hand-side
#'   formula. Defaults to all fixed-effect terms that marginality
#'   allows ([stats::drop.scope()]).
#' @param test `"Chisq"` adds likelihood-ratio tests.
#' @param k AIC penalty per parameter.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return An `anova` table with one row per dropped term.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), z = rnorm(100),
#'                  g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x * z + (1 | g)) + gaussian(), data = dd)
#'
#' # marginality keeps the main effects out of scope while x:z is in it
#' drop1(fit)
#' drop1(fit, test = "Chisq")
#'
#' # name the terms to override the default scope
#' drop1(fit, scope = ~ x + z, test = "Chisq")
#' @export
drop1.frmtmb_fit <- function(object, scope, test = c("none", "Chisq"),
                             k = 2, ...) {
  # `scale` and `trace` are not named here: they are in
  # `s3_contract_args`, with every other name R's own machinery passes
  # through a generic. step() reaches this one after nobs().
  frm_check_dots(...)
  test <- frm_match_arg(test)
  if (object$REML) {
    frm_stop("drop1() compares fixed effects; refit with REML = FALSE",
             call. = FALSE)
  }
  if (length(object$spec$responses) > 1) {
    frm_stop("drop1() is not supported for multivariate fits", call. = FALSE)
  }
  tt <- terms(object)
  labs <- attr(tt, "term.labels")
  if (missing(scope)) {
    scope <- stats::drop.scope(tt)
  } else if (!is.character(scope)) {
    scope <- attr(stats::terms(stats::update.formula(
      stats::formula(tt), scope)), "term.labels")
  }
  bad <- setdiff(scope, labs)
  if (length(bad)) {
    frm_stop("scope is not a subset of the term labels: ",
             paste(bad, collapse = ", "), call. = FALSE)
  }

  ll0 <- logLik(object)
  df0 <- attr(ll0, "df")
  aic0 <- -2 * as.numeric(ll0) + k * df0
  n_sc <- length(scope)
  ddf <- aic <- lrt <- rep(NA_real_, n_sc)
  for (i in seq_len(n_sc)) {
    # rebuild the bf() object with the term removed; the stored model
    # frame carries every variable, so the refit does not depend on
    # the original data still being visible
    nb <- object$bform
    nf <- stats::update.formula(nb$formula,
                                paste(". ~ . -", scope[i]))
    environment(nf) <- environment(nb$formula)
    nb$formula <- nf
    cl <- object$call
    cl$formula <- nb
    cl$data <- object$frame[["data_frame"]]
    # same reason for data2: the stored structural objects go in by
    # value, so the refit does not need the names the user passed to
    # still resolve where the call is evaluated
    if (length(object$data2)) cl$data2 <- object$data2
    fit_i <- eval(cl, environment(nb$formula) %||% parent.frame())
    ll_i <- logLik(fit_i)
    ddf[i] <- df0 - attr(ll_i, "df")
    aic[i] <- -2 * as.numeric(ll_i) + k * attr(ll_i, "df")
    lrt[i] <- 2 * (as.numeric(ll0) - as.numeric(ll_i))
  }
  tab <- data.frame(Df = c(NA, ddf), AIC = c(aic0, aic),
                    row.names = c("<none>", scope), check.names = FALSE)
  if (test == "Chisq") {
    tab$LRT <- c(NA, lrt)
    tab$`Pr(>Chi)` <- c(NA, lrt_pvalue(lrt, ddf))
  }
  structure(tab, class = c("anova", "data.frame"),
            heading = c("Single term deletions\n",
                        paste("Model:", model_label(object)), ""))
}

#' Does a formula carry a `.` term, which is what makes it a delta
#' against an existing formula rather than a complete one? Both
#' [stats::update()] spellings have one: the one-sided `~ . + z` and
#' the dotted-LHS `. ~ . + z`.
#'
#' @noRd
formula_has_dot <- function(f) {
  has <- function(e) {
    if (is.name(e)) return(identical(as.character(e), "."))
    if (!is.call(e)) return(FALSE)
    parts <- as.list(e)[-1L]
    for (p in parts) {
      if (identical(p, quote(expr = ))) next
      if (has(p)) return(TRUE)
    }
    FALSE
  }
  has(f)
}

#' Apply an update formula (`~ . + z` or `. ~ . + z`, the
#' [stats::update.formula] and brms spellings) to the stored model
#' formula. The whole `bf()` goes back into the call, so dpar formulas,
#' fixed dpar values and the family survive the update.
#'
#' @noRd
update_delta_formula <- function(object, f) {
  bform <- object$bform
  if (inherits(bform, "frmtmb_mvformula")) {
    frm_stop("An update formula written as a delta does not say which ",
             "response it changes. Pass the complete mvbf() as `formula`",
             call. = FALSE)
  }
  # nlf() on mu makes the response formula's right-hand side a body too,
  # whether or not nl = TRUE was written
  nl_mu <- isTRUE(bform$nl) ||
    length(intersect(all.vars(reformulas::RHSForm(bform$formula)),
                     names(bform$nlforms %||% list()))) > 0L
  if (nl_mu) {
    frm_stop("An update formula written as a delta cannot be applied to a ",
             "nonlinear formula, whose right-hand side is an expression ",
             "and not a sum of terms. Pass the complete ",
             "bf(..., nl = TRUE) as `formula`", call. = FALSE)
  }
  new <- stats::update.formula(bform$formula, f)
  environment(new) <- environment(bform$formula)
  bform$formula <- new
  bform
}

#' Update and refit a model
#'
#' Re-evaluates the stored [frm()] call with the given arguments
#' replaced. Any `frm()` argument can be updated by name.
#'
#' The formula argument is `formula.`, as in [stats::update()] and in
#' brms; `formula = ` reaches it by partial matching. A formula
#' carrying a `.` is a delta applied to the stored `mu` formula with
#' [stats::update.formula] semantics - one-sided `~ . + z`, dotted
#' `. ~ . + z`, or a changed response `z ~ . + x` - and keeps the dpar
#' formulas, the fixed dpar values and the family. A formula with no
#' `.` replaces the stored one. brms's `newdata` is accepted as a
#' synonym for `data`.
#'
#' @param object A `frmtmb_fit`.
#' @param formula. A complete formula or [bf()], or a delta such as
#'   `~ . + z` or `. ~ . + z`.
#' @param ... Arguments of [frm()] to replace, e.g. `data`, `family`,
#'   `REML`. `newdata` is accepted for `data`.
#' @param evaluate If `FALSE`, return the updated call instead of the
#'   refitted model.
#' @return A `frmtmb_fit`, or the updated call when
#'   `evaluate = FALSE`.
#' @examples
#' set.seed(3)
#' dd <- data.frame(x = rnorm(60), z = rnorm(60))
#' dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
#' fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
#'
#' # a delta on the stored formula, in either spelling
#' fit2 <- update(fit, ~ . + z)
#' formula(fit2)
#' formula(update(fit, . ~ . + z))
#' fit3 <- update(fit, formula. = ~ . - x, newdata = dd[1:40, ])
#' nobs(fit3)
#' @export
update.frmtmb_fit <- function(object, formula., ..., evaluate = TRUE) {
  cl <- object$call
  extras <- match.call(expand.dots = FALSE)$...
  # brms spells the data argument of update() `newdata`
  if ("newdata" %in% names(extras)) {
    if ("data" %in% names(extras)) {
      frm_stop("Give the updated data once: as `data`, or as brms's ",
               "`newdata`, but not both", call. = FALSE)
    }
    names(extras)[names(extras) == "newdata"] <- "data"
  }
  for (nm in names(extras)) cl[[nm]] <- extras[[nm]]
  if (!missing(formula.)) {
    # a complete formula goes in unevaluated, so the stored call stays
    # readable; a delta has to be resolved against the stored one
    cl$formula <- if (inherits(formula., "formula") &&
                      (length(formula.) == 2L ||
                         formula_has_dot(formula.))) {
      update_delta_formula(object, formula.)
    } else {
      substitute(formula.)
    }
  }
  if (!evaluate) return(cl)
  # the stored structural objects go into the call by value, so an
  # update in a session where the original data2 names are gone still
  # assembles; an explicit data2 = in the update wins
  if (length(object$data2) && !("data2" %in% names(extras))) {
    cl$data2 <- object$data2
  }
  eval(cl, parent.frame())
}

#' Likelihood profiles
#'
#' Wraps [TMB::tmbprofile()] per parameter. The returned objects have
#' `plot()` and `confint()` methods (from TMB).
#'
#' @param fitted A `frmtmb_fit`.
#' @param parm Parameter names or indices. Required; profiling is not
#'   free, so there is no all-parameters default. The names are the
#'   ones [confint()] takes, in any of its three spellings: internal
#'   (`theta_1`, `tarsus_(Intercept)`), parenthesis-free
#'   (`tarsus_Intercept`), or a one-to-one natural-scale alias
#'   (`sd_dam__Intercept`). The profile is of the internal parameter
#'   either way - a log standard deviation, not a standard deviation -
#'   and the returned element keeps the internal name. For a profile of
#'   a natural-scale quantity itself, including one that mixes several
#'   parameters, use `hypothesis(method = "profile")`.
#' @param ... Passed to [TMB::tmbprofile()].
#' @return A `tmbprofile` data frame, or a named list of them when
#'   `parm` has length above one.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # parameter names are the confint() row names, with or without the
#' # parentheses, or a one-to-one natural-scale alias
#' rownames(confint(fit))
#' pr <- profile(fit, "theta_1")
#' identical(profile(fit, "(Intercept)"), profile(fit, "Intercept"))
#' plot(pr)
#' # TMB's confint() reads the interval off the profile
#' confint(pr)
#'
#' # several parameters at once return a named list
#' prs <- profile(fit, c("x", "theta_1"))
#' names(prs)
#' @export
profile.frmtmb_fit <- function(fitted, parm, ...) {
  nm <- outer_par_names(fitted)
  idx <- resolve_par_index(fitted, parm, "profile")
  out <- lapply(idx, function(i) {
    TMB::tmbprofile(fitted$obj, name = i, trace = FALSE, ...)
  })
  names(out) <- nm[idx]
  if (length(out) == 1L) out[[1L]] else out
}

## hypothesis(): the expression environment and its parameter mapping.

#' Parameter values without any covariance machinery (usable on refits
#' inside a bootstrap without triggering sdreport).
#'
#' @noRd
hyp_vals_only <- function(fit) {
  est <- fit$estimates
  bd <- est[["betad"]]
  if (length(fx <- fit$frame[["betad_fixed_idx"]])) bd <- bd[-fx]
  vals <- c(est[["beta"]], bd,
            est[["theta"]], est[["thetaac"]], est[["thetar"]])
  comp <- c(rep("beta", length(est[["beta"]])), rep("betad", length(bd)),
            rep("theta", length(est[["theta"]])),
            rep("thetaac", length(est[["thetaac"]])),
            rep("thetar", length(est[["thetar"]])))
  # the ordinal thresholds and the category-specific coefficients: brms
  # reports both as parameters, and without them variables() listed 3 of
  # 9 names on an ordinal fit and hypothesis() could not reach a
  # threshold at all
  for (cp in ord_extra_comps(fit)) {
    v <- est[[cp]]
    vals <- c(vals, v)
    comp <- c(comp, rep(cp, length(v)))
  }
  list(vals = vals, comp = comp)
}

#' The parameter-template components that hold an ordinal fit's
#' thresholds and its `cs()` coefficients.
#'
#' They are `extra_names` components rather than coefficients, so the
#' coefficient machinery never saw them. They are parameters all the
#' same, and brms names them `b_Intercept[k]` and `bcs_<term>[k]`.
#'
#' @noRd
ord_extra_comps <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  out <- character(0)
  if (length(tpl[["tau_raw"]])) out <- "tau_raw"
  for (lp in fit$frame[["linpreds"]]) {
    for (ct in lp[["cs"]] %||% list()) out <- c(out, ct[["par"]])
  }
  intersect(unique(out), names(tpl))
}

#' The thresholds themselves, from the internal vector, through the map
#' the family declares.
#'
#' @noRd
ord_threshold_values <- function(fam, raw) {
  f <- fam[["post"]][["ord_thresholds"]]
  if (is.null(f)) raw else f(raw)
}

#' Values plus joint covariance of (beta, estimated betad, theta,
#' thetaac, thetar). ML: straight from cov.fixed in opt$par order
#' (outer_pos maps back into the full outer vector for tmbroot
#' lincombs). REML: beta is integrated out, so the blocks come from the
#' joint precision.
#'
#' @noRd
hyp_par_cov <- function(fit) {
  comps <- c("beta", "betad", "theta", "thetaac", "thetar",
             ord_extra_comps(fit))
  if (!fit$REML && !isTRUE(fit$control$profile)) {
    sdr <- sdr_of(fit)
    V <- sdr$cov.fixed
    rn <- rownames(V)
    keep <- which(rn %in% comps)
    # par.fixed equals opt$par at the optimum (this branch is never
    # taken under control profile = TRUE)
    list(vals = unname(sdr$par.fixed[keep]), comp = rn[keep],
         V = V[keep, keep, drop = FALSE], outer_pos = keep,
         n_outer = length(fit$opt$par))
  } else {
    Q <- sdr_of(fit)$jointPrecision
    Vall <- solve_joint_precision(Q, fit$cache, fit)
    rn <- rownames(Q)
    keep <- which(rn %in% comps)
    vo <- hyp_vals_only(fit)
    vals <- numeric(length(keep))
    cnt <- stats::setNames(integer(length(comps)), comps)
    for (i in seq_along(keep)) {
      k <- rn[keep[i]]
      cnt[k] <- cnt[k] + 1L
      vals[i] <- vo$vals[vo$comp == k][cnt[k]]
    }
    list(vals = vals, comp = rn[keep],
         V = as.matrix(Vall[keep, keep, drop = FALSE]), outer_pos = NULL,
         n_outer = length(fit$opt$par))
  }
}

#' Named list of every variable a hypothesis can name, under brms's
#' names, at one parameter vector:
#'
#' - the coefficients, `b_Intercept`, `b_sigma_Intercept`, `bs_sx_1`;
#' - a distributional parameter nobody wrote a formula for, on its
#'   natural scale, `sigma`, `shape`, `sigma_ya` (`brms_coef_table()`);
#' - the group-level standard deviations and correlations,
#'   `sd_<group>__<coef>` and `cor_<group>__<c1>__<c2>`;
#' - the residual autocorrelation parameters, `ar[1]`, `cosy`;
#' - the residual correlations of a multivariate model,
#'   `rescor__<resp1>__<resp2>`.
#'
#' Every name is brms's. brms refuses a model whose renaming gives two
#' coefficients one name, suffixes a clash across predictors with
#' `__1`, and refuses a duplicated group-level effect; the model frame
#' applies the same rules, so a name reaching here twice is a defect and
#' stops rather than keeping one of the two values.
#'
#' @noRd
hyp_env_vals <- function(fit, vals, comp) {
  env <- list()
  put <- function(nm, val) {
    if (!is.null(env[[nm]])) {
      frm_stop("Internal error: two parameters of this model share the ",
               "name '", nm, "'. Please report it with the model formula",
               call. = FALSE)
    }
    env[[nm]] <<- unname(val)
    invisible(NULL)
  }

  cf <- c(vals[comp == "beta"], vals[comp == "betad"])
  tab <- brms_coef_table(fit)
  inv <- attr(tab, "linkinv")
  for (i in which(!tab$natural)) put(tab$brms[i], cf[i])

  th <- vals[comp == "theta"]
  sds_nm <- brms_sds_names(fit)
  for (bi in seq_along(fit$frame[["re_blocks"]])) {
    bk <- fit$frame[["re_blocks"]][[bi]]
    # a smooth's one parameter is brms's smoothing sd, sds_sx_1
    if (!is.na(sds_nm[bi])) {
      put(sds_nm[bi], exp(th[bk[["theta_idx"]][1L]]))
      next
    }
    # Excluded: the structures whose theta segment is not a set of
    # standard deviations and correlations at all. `smooth` carries one
    # inverse smoothing parameter, `gp`/`hsgp` a marginal sd plus
    # lengthscales, `car` an sd plus a mixing proportion, `spde` a
    # precision and an inverse range. Their summaries live in
    # confint_varcorr() under their own names.
    #
    # Included: gr_cov, gr_prec and equalto, whose registry vcov() is the
    # WITHIN-level covariance brms names sd_<group>__<term>.
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "car", "spde")) next
    V <- covstruct_registry[[bk[["covstruct"]]]]$vcov(th[bk[["theta_idx"]]], bk)
    tn <- brms_re_rnames(fit, bk)
    g <- brms_group_name(bk)
    sds <- sqrt(diag(V))
    for (j in seq_along(sds)) put(paste0("sd_", g, "__", tn[j]), sds[j])
    if (nrow(V) > 1L) {
      C <- stats::cov2cor(V)
      for (j in seq_len(nrow(V) - 1L)) {
        for (k in seq(j + 1L, nrow(V))) {
          put(paste0("cor_", g, "__", tn[j], "__", tn[k]), C[j, k])
        }
      }
    }
  }

  # distributional parameters on their natural scale, after the group
  # summaries, where brms lists them
  smp <- attr(tab, "simplex")
  in_smp <- unlist(lapply(smp, `[[`, "pos"))
  for (i in setdiff(which(tab$natural), in_smp)) {
    put(tab$brms[i], inv[[i]](cf[i]))
  }
  # a mixture's weights, brms's theta1 ... thetaK, from all K - 1 log
  # ratios at once
  for (s in smp) {
    p <- s$to_simplex(cf[s$pos])
    for (k in seq_along(s$names)) put(s$names[k], p[1L, k])
  }

  thac <- vals[comp == "thetaac"]
  for (ac in fit$frame[["autocor"]] %||% list()) {
    nat <- autocor_natural(thac[ac[["theta_idx"]]], ac)
    for (j in seq_along(nat)) put(names(nat)[j], nat[j])
  }

  hyp_put_ordinal(fit, vals, comp, put)

  thr <- vals[comp == "thetar"]
  if (isTRUE(fit$spec$rescor) && length(thr)) {
    rs <- brms_stan_name(names(fit$spec$responses))
    C <- us_chol_cor(thr, length(rs))
    for (i in seq_along(rs)[-1L]) {
      for (j in seq_len(i - 1L)) {
        put(paste0("rescor__", rs[j], "__", rs[i]), C[i, j])
      }
    }
  }
  env
}

#' brms's names for an ordinal fit's thresholds and `cs()`
#' coefficients: `b_Intercept[k]` for the threshold and
#' `bcs_<term>[k]` for the category-specific coefficient, each prefixed
#' by the predictor the way every other coefficient name is.
#'
#' The thresholds of a model with more than one ordinal response are
#' NOT named: one template component holds them all and nothing here
#' says which belongs to which response, so naming them would be a
#' guess. The `cs()` coefficients are per predictor and are named
#' whatever the model looks like.
#'
#' @noRd
hyp_put_ordinal <- function(fit, vals, comp, put) {
  raw <- vals[comp == "tau_raw"]
  if (length(raw)) {
    ord_lps <- Filter(function(lp) {
      identical(brms_lp_family(fit, lp)[["type"]], "ordinal") &&
        identical(lp[["dpar"]], "mu")
    }, fit$frame[["linpreds"]])
    if (length(ord_lps) == 1L) {
      lp <- ord_lps[[1L]]
      th <- ord_threshold_values(brms_lp_family(fit, lp), raw)
      pre <- brms_lp_prefix(fit, lp)
      for (k in seq_along(th)) {
        put(paste0("b_", brms_usc(pre, "Intercept"), "[", k, "]"), th[k])
      }
    }
  }
  for (lp in fit$frame[["linpreds"]]) {
    for (ct in lp[["cs"]] %||% list()) {
      v <- vals[comp == ct[["par"]]]
      if (!length(v)) next
      lab <- brms_rename(sub("^cs", "", ct[["label"]]))
      pre <- brms_lp_prefix(fit, lp)
      for (k in seq_along(v)) {
        put(paste0("bcs_", brms_usc(pre, lab), "[", k, "]"), v[k])
      }
    }
  }
  invisible(NULL)
}

#' Split one hypothesis string the way `brms:::eval_hypothesis()` does:
#' whitespace removed, the two sides around the one sign, and the text
#' `(lhs)-(rhs)` with the right side dropped when it is `0`. Returns
#' that text and the alternative.
#'
#' A string with NO sign used to be read as the quantity tested against
#' zero, so `hypothesis(fit, "Age")` answered the row of `"Age = 0"`
#' with nothing said, where brms refuses: "Every hypothesis must be of
#' the form 'left (= OR < OR >) right'". `brms:::eval_hypothesis()`
#' requires exactly one sign and exactly two sides, and so does this now
#' (dev/adefects-findings.md, D3).
#'
#' @noRd
hyp_parse <- function(h) {
  h <- gsub("[ \t\r\n]", "", h)
  ops <- unlist(gregexpr("[<>]", h))
  ops <- ops[ops > 0L]
  if (length(ops)) {
    if (length(ops) > 1L) {
      frm_stop("A hypothesis has at most one '<' or '>': '", h, "'",
               call. = FALSE)
    }
    op <- substr(h, ops, ops)
    lhs <- substr(h, 1L, ops - 1L)
    # ">=" and "<=" read as ">" and "<": the boundary has probability
    # zero under every sampling distribution used here
    rhs <- sub("^=", "", substring(h, ops + 1L))
    if (grepl("=", lhs, fixed = TRUE) || grepl("=", rhs, fixed = TRUE)) {
      frm_stop("A hypothesis is directional ('<', '>') or an equality ",
               "('='), not both: '", h, "'", call. = FALSE)
    }
    return(list(text = hyp_two_sides(lhs, rhs, h),
                dir = if (op == ">") "greater" else "less"))
  }
  eq <- strsplit(h, "=", fixed = TRUE)[[1L]]
  if (length(eq) > 2L) {
    frm_stop("A hypothesis has at most one '=': '", h, "'", call. = FALSE)
  }
  if (length(eq) < 2L) {
    frm_stop("Every hypothesis must be of the form 'left (= OR < OR >) ",
             "right': '", h, "' states no relation. Write '", h, " = 0' ",
             "for the test against zero, which is what brms writes",
             call. = FALSE)
  }
  list(text = hyp_two_sides(eq[1L], eq[2L], h), dir = "two.sided")
}

#' brms's `(lhs)-(rhs)` text of a two-sided hypothesis.
#'
#' @noRd
hyp_two_sides <- function(lhs, rhs, h) {
  if (!nzchar(lhs) || !nzchar(rhs)) {
    frm_stop("Every hypothesis must be of the form 'left (= OR < OR >) ",
             "right': '", h, "'", call. = FALSE)
  }
  paste0("(", lhs, ")", if (rhs != "0") paste0("-(", rhs, ")"))
}

#' brms's `find_vars()`: the variable names of a hypothesis text. A name
#' may carry `:`, `.`, `_` and one bracketed index, `b_x:fe` and
#' `ar[1]`; a function name before `(` and the digits of a decimal
#' number are not variables.
#'
#' @noRd
hyp_find_vars <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  lead <- "([^([:digit:]|[:punct:])]|\\.)"
  pos_all <- gregexpr(paste0(lead, "[[:alnum:]_\\:\\.]*",
                             "(\\[[^],]+(,[^],]+)*\\])?"), x)[[1L]]
  pos_fun <- gregexpr(paste0(lead, "[[:alnum:]_\\.]*\\("), x)[[1L]]
  pos_dec <- gregexpr("\\.[[:digit:]]+", x)[[1L]]
  keep <- !pos_all %in% c(pos_fun, pos_dec)
  pos <- pos_all[keep]
  attr(pos, "match.length") <- attr(pos_all, "match.length")[keep]
  if (!length(pos)) return(character(0))
  unique(unlist(regmatches(x, list(pos))))
}

#' brms's hypothesis renaming, which makes a parameter name a legal R
#' name inside the expression: `:` to `___`, `[` and `]` to `.`, `,` to
#' `..`. Without it `b_x:fe` parses as R's sequence operator.
#'
#' @noRd
hyp_rename <- function(x) {
  brms_rename(x, c(":", "[", "]", ","), c("___", ".", ".", ".."))
}

#' Parse every hypothesis string of one call: the shared front end of
#' the `frmtmb_fit`, `frmtmb_draws` and `frmtmb_multiple` methods, and
#' of `hypothesis(scope = )` on draws.
#'
#' As `brms:::eval_hypothesis()`: the variables of the text are found,
#' each is prefixed with `class` and `group`, a prefixed name that is
#' not in `known` is refused with brms's message, and the text is
#' renamed before it is parsed. Each expression carries a `"vars"`
#' attribute mapping the renamed name it contains to the full parameter
#' name, which is what `hyp_eval()` binds.
#'
#' @noRd
hyp_parse_all <- function(hypothesis, known, class = "b", group = "") {
  if (!is.character(hypothesis) || !length(hypothesis) ||
        anyNA(hypothesis)) {
    frm_stop("Argument 'hypothesis' must be a character vector.", call. = FALSE)
  }
  prefix <- hyp_class_prefix(class, group)
  ps <- lapply(hypothesis, hyp_parse)
  exprs <- lapply(ps, function(p) {
    vars <- hyp_find_vars(p$text)
    full <- paste0(prefix, vars)
    miss <- setdiff(full, known)
    if (length(miss)) {
      frm_stop("Some parameters cannot be found in the model: \n",
               paste0("'", miss, "'", collapse = ", "),
               "\nvariables() lists every name; brms's default class = \"b\" ",
               "puts b_ before each name, so a name such as sigma or ",
               "sd_<group>__<coef> needs class = NULL", call. = FALSE)
    }
    ex <- str2lang(hyp_rename(p$text))
    # a backticked name such as `(Intercept)` is one symbol to R but not
    # a variable to hyp_find_vars(), so it would reach eval() unbound
    # and fail there with base R's "object not found"
    extra <- setdiff(all.vars(ex), hyp_rename(vars))
    extra <- extra[!vapply(extra, exists, NA, envir = environment())]
    if (length(extra)) {
      frm_stop("Some parameters cannot be found in the model: \n",
               paste0("'", extra, "'", collapse = ", "),
               "\nvariables() lists every name; write a name without ",
               "backticks, as brms does", call. = FALSE)
    }
    attr(ex, "vars") <- stats::setNames(full, hyp_rename(vars))
    ex
  })
  list(exprs = exprs, dir = vapply(ps, function(p) p$dir, ""))
}

#' The full parameter names one parsed hypothesis reads.
#'
#' @noRd
hyp_expr_vars <- function(ex) unname(attr(ex, "vars") %||% character(0))

#' Evaluate one parsed hypothesis over a named list of values keyed by
#' full parameter name: a number at a parameter vector, or a vector of
#' draws.
#'
#' @noRd
hyp_eval_in <- function(ex, values) {
  map <- attr(ex, "vars") %||% character(0)
  env <- lapply(unname(map), function(f) values[[f]])
  names(env) <- names(map)
  miss <- map[vapply(env, is.null, TRUE)]
  if (length(miss)) {
    frm_stop("Some parameters cannot be found in the model: \n",
             paste0("'", miss, "'", collapse = ", "), call. = FALSE)
  }
  attr(ex, "vars") <- NULL
  eval(ex, env, parent.frame())
}

#' Evaluate a parsed hypothesis at one parameter vector.
#'
#' @noRd
hyp_eval <- function(fit, ex, vals, comp) {
  hyp_eval_in(ex, hyp_env_vals(fit, vals, comp))
}

#' The name prefix implied by brms's `class` and `group` shorthand,
#' computed as `brms:::hypothesis.brmsfit()` computes it: `NULL` and
#' `""` are no prefix, a group makes `<class>_<group>__`, and a class
#' alone makes `<class>_`. So brms's default `class = "b"` reads a bare
#' `x` as the coefficient `b_x`, and a natural-scale name such as
#' `sd_g__Intercept` needs `class = NULL`, exactly as in brms.
#'
#' @noRd
hyp_class_prefix <- function(class = "b", group = "") {
  if (!length(class)) class <- ""
  if (!is.character(class) || length(class) != 1L || is.na(class)) {
    frm_stop("`class` must be a single string, or NULL for no prefix",
             call. = FALSE)
  }
  if (!length(group)) group <- ""
  if (!is.character(group) || length(group) != 1L || is.na(group)) {
    frm_stop("`group` must be a single string", call. = FALSE)
  }
  if (nzchar(group)) return(paste0(class, "_", group, "__"))
  if (nzchar(class)) return(paste0(class, "_"))
  ""
}

#' The `class` element brms stores on a hypothesis result: the prefix
#' with its trailing underscores dropped (`"b"`, `"sd_g"`, `""`).
#'
#' @noRd
hyp_class_label <- function(prefix) sub("_+$", "", prefix)

#' One-sided quantile bookkeeping: the interval bound and the p-value a
#' given alternative asks for, from an estimate, a standard error and a
#' normal (or t) reference.
#'
#' @noRd
hyp_wald_row <- function(est, se, dir, alpha, qfun, pfun) {
  z <- est / se
  q <- qfun(1 - if (dir == "two.sided") alpha / 2 else alpha)
  list(
    lwr = if (dir == "less") -Inf else est - q * se,
    upr = if (dir == "greater") Inf else est + q * se,
    stat = z,
    p = switch(dir,
               two.sided = 2 * pfun(-abs(z)),
               greater = pfun(-z),
               less = pfun(z))
  )
}

#' Tail proportion of a draws vector against zero, in the direction the
#' hypothesis asks for, with the (1 + k) / (1 + n) correction that
#' keeps a p-value away from exactly zero.
#'
#' @noRd
hyp_tail_p <- function(t, dir) {
  n <- length(t)
  switch(dir,
         two.sided = min(1, 2 * min((1 + sum(t <= 0)) / (1 + n),
                                    (1 + sum(t >= 0)) / (1 + n))),
         greater = (1 + sum(t <= 0)) / (1 + n),
         less = (1 + sum(t >= 0)) / (1 + n))
}

#' Central-difference gradient of a scalar function of the parameter
#' vector. The delta method needs a gradient, and a hypothesis is an
#' arbitrary R expression with no derivative available.
#'
#' @noRd
hyp_fd_grad <- function(f, v) {
  vapply(seq_along(v), function(i) {
    step <- max(1e-5, 1e-5 * abs(v[i]))
    vp <- v; vp[i] <- vp[i] + step
    vm <- v; vm[i] <- vm[i] - step
    (f(vp) - f(vm)) / (2 * step)
  }, numeric(1))
}

#' Hypothesis tests on parameter expressions
#'
#' The frequentist analog of brms's `hypothesis()`: evaluates
#' expressions of the model parameters at the estimates and tests them
#' against zero. A hypothesis is `"lhs = rhs"`, e.g. `"x1 - x2 = 0"` or
#' `"exp(Intercept) = 1"`, or brms's directional `"lhs > rhs"` /
#' `"lhs < rhs"`. Every hypothesis states a relation: a string with no
#' `=`, `<` or `>` is refused, as it is in brms, so write `"x1 = 0"`
#' rather than `"x1"`.
#'
#' @section The returned object:
#' brms's SHAPE under frmtmb's own class: a list of class
#' `"frmtmb_hypothesis"` with the elements brms has, in brms's order.
#' It does not carry brms's `brmshypothesis` class. frmtmb owns
#' `print()` and `plot()` for its own class and a frmtmb fit is not a
#' brms fit, so `is(x, "brmshypothesis")` in a ported script is a
#' rule-2 divergence like the other seventeen the port ledger records.
#'
#' - `hypothesis`: a data frame with brms's eight columns, one row per
#'   hypothesis. On a maximum-likelihood fit they mean:
#'   - `Hypothesis`: brms's label, `(lhs)-(rhs) > 0`, or the name given
#'     on the hypothesis vector.
#'   - `Estimate`: the expression at the estimates.
#'   - `Est.Error`: its delta-method standard error (`"wald"`,
#'     `"profile"`), the bootstrap standard deviation (`"boot"`), or the
#'     Rubin pooled standard error (a [frm_multiple()] result).
#'   - `CI.Lower`, `CI.Upper`: brms's interval, which is central at
#'     `1 - alpha` for `"="` and central at `1 - 2 * alpha` for a
#'     directional row, so that its relevant end is the one-sided bound.
#'     Wald, profile-likelihood or bootstrap percentile, per `method`.
#'   - `Evid.Ratio`, `Post.Prob`: `NA`. They are posterior quantities,
#'     and a fit has no posterior.
#'   - `Star`: `"*"` when a two-sided row's interval excludes 0, or when
#'     a directional row's one-sided test rejects at level `alpha`
#'     (brms stars a posterior probability above `1 - alpha` there).
#' - `samples`: brms's frame of draws, columns `H1`, `H2`, ...: the
#'   bootstrap replicates for `"boot"`, and no rows otherwise.
#' - `prior_samples`: the same columns, all `NA`.
#' - `class`: the name prefix `class` and `group` produced, without its
#'   trailing underscores, as brms stores it.
#' - `alpha`.
#'
#' What brms's frame has no column for rides on attributes, so the
#' frame keeps brms's shape: `attr(h, "test")` is the test statistic
#' and its p-value per row (`z`, or `t` with `df` for a pooled or
#' degrees-of-freedom-carrying covariance), `attr(h, "method")`, and the
#' method payload, `attr(h, "draws")` for the bootstrap matrix and
#' `attr(h, "profiles")` for the profile curves.
#'
#' @section Directional hypotheses:
#' `"lhs > rhs"` and `"lhs < rhs"` test the same difference
#' `(lhs) - (rhs)` against zero with a one-sided alternative, so the
#' reported `p` is the one-sided tail probability. `p` is `pnorm()` of
#' the signed z statistic for `"wald"` and, as in the two-sided case
#' where the standard error and the statistic stay Wald-based, for
#' `"profile"` too: the profile changes the INTERVAL and nothing else.
#' For `"boot"` `p` is the tail proportion of the replicates with the
#' `(1 + k) / (1 + n)` correction. Where brms reports the posterior
#' probability of the direction, this reports its frequentist
#' complement: small `p` is evidence for the stated direction.
#' `">="` and `"<="` read as `">"` and `"<"`.
#'
#' @section brms class and group shorthand:
#' The hypothesis is read as `brms:::eval_hypothesis()` reads it. Every
#' variable in the text gets the prefix `class` and `group` make, and a
#' prefixed name the model does not have is refused with brms's message,
#' "Some parameters cannot be found in the model". brms's default
#' `class = "b"` reads `"x1 - x2 = 0"` as `b_x1 - b_x2`, so it refuses
#' `"b_x1 = 0"`, which it reads as `b_b_x1`. `class = "sd", group =
#' "patient"` reads `"Intercept - age > 0"` as `sd_patient__Intercept -
#' sd_patient__age`. `class = NULL` (or `""`) takes every name as
#' written, which is what `sigma` or `sd_g__Intercept` needs, as in
#' brms.
#'
#' A name may carry what brms's names carry: `x:fe` is the interaction
#' coefficient `b_x:fe`, not R's `:` operator, and `ar[1]` is the
#' autocorrelation. brms's renaming (`:` to `___`, `[` and `]` to `.`,
#' `,` to `..`) is applied to the text before it is parsed, as in brms.
#'
#' Available names, which [variables()] lists, are brms's: the
#' coefficients (`b_Intercept`, `b_x`, `b_IxE2` for `I(x^2)`,
#' `b_sigma_Intercept` for a `sigma` formula, `bs_sx_1` for the
#' unpenalized part of `s(x)`); a distributional parameter nobody wrote
#' a formula for, on its natural scale (`sigma`, `shape`, `nu`, and
#' `sigma_ya` for response `y_a` of a multivariate model); the
#' group-level summaries `sd_<group>__<coef>` and
#' `cor_<group>__<coef1>__<coef2>`, where the coefficient of a
#' distributional or nonlinear parameter carries that parameter's name
#' (`sd_g__sigma_Intercept`); the autocorrelation parameters (`ar[1]`,
#' `cosy`); and the residual correlations `rescor__<resp1>__<resp2>`.
#' So an ICC is
#' `hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2) = 0",
#' class = NULL)`. The traffic runs the other way too: `confint(parm = )`
#' and `profile(parm = )` take these names, whenever one of them stands
#' for a single internal parameter.
#'
#' @section Which random-effect blocks contribute names:
#' Every block whose covariance parameters ARE standard deviations and
#' correlations: the plain structures (`us`, `diag`, `homdiag`, `cs`,
#' `ar1`, `toep`, the spatial and reduced-rank ones) and the
#' known-structure blocks `gr(cov = )`, `gr(prec = )` and `equalto()`,
#' whose `sd_`/`cor_` names describe the WITHIN-level covariance that
#' multiplies the fixed relationship matrix. That is what makes
#' heritability-as-ICC writable directly:
#' `"sd_id__Intercept^2 / (sd_id__Intercept^2 + sigma^2) = 0"` with
#' `class = NULL` on an animal model fitted with
#' `(1 | gr(id, cov = A))`. An `equalto()` block estimates nothing, so
#' its names are constants with zero variance.
#'
#' An `|ID|`-merged block is ONE block, so it contributes one name per
#' merged coefficient, and the names carry the predictor they came from
#' in brms's spelling. A two-trait animal model written
#' `(1 | q | gr(id, cov = A))` in both formulas of an [mvbf()] gives
#' `sd_id__y1_Intercept`, `sd_id__y2_Intercept` and
#' `cor_id__y1_Intercept__y2_Intercept`, the last being the genetic
#' correlation between the traits. [variables()] prints them.
#'
#' Two terms that give one grouping factor the same coefficient, an
#' animal model's `(1 | gr(id, cov = A)) + (1 | id)`, are refused when
#' the model is built, with brms's message "Duplicated group-level
#' effects are not allowed". Give the second term a copy of the factor
#' under another name, `(1 | gr(id, cov = A)) + (1 | id_pe)`, and the
#' two are `sd_id__Intercept` and `sd_id_pe__Intercept`.
#'
#' Excluded: `s()`/`t2()` smooths, `gp()`/`hsgp()`, `car()` and `spde()`.
#' Their theta segments are not standard deviations: an inverse
#' smoothing parameter, lengthscales, a mixing proportion, a precision
#' and an inverse range. There is no `sd_<group>__<coef>` to name. Read
#' those off [confint_varcorr()], which reports each under its own
#' label (`sd(gp)`, `range(gp)`, `sd(car)`, ...).
#'
#' @section Names that would collide:
#' brms's renaming can give two parameters one name, and frmtmb does
#' what brms does with each case:
#'
#' - Within one predictor, a covariate whose renamed column repeats
#'   another's, `y ~ Intercept + x` (`(Intercept)` and `Intercept` are
#'   both `Intercept`), is refused with brms's "Internal renaming led to
#'   duplicated names".
#' - Across predictors, the later name takes brms's `__1` suffix:
#'   in `bf(y ~ sigma_z, sigma ~ z)` the covariate `sigma_z` of `mu` and
#'   the coefficient `z` of `sigma` are `b_sigma_z` and `b_sigma_z__1`.
#' - A group-level coefficient given twice on one group is refused, see
#'   above.
#' - Two responses brms spells alike, `y_a` and `ya`, are refused with
#'   brms's "Cannot use the same response variable twice".
#' - A group-level label given twice on draws, from levels `lvl 1` and
#'   `lvl.1`, takes brms's `__1` on the later level.
#' - An interaction group two of whose levels brms joins to one string
#'   (`1_2:3` and `1:2_3`) is refused: brms pools the two levels.
#'
#' A coefficient and a natural-scale name cannot meet: every coefficient
#' starts `b_`, `bs_` or `bsp_`, and `sigma` does not.
#'
#' A mixture's weights with no theta formula are brms's simplex,
#' `theta1 ... thetaK`, the mixing probabilities, computed together from
#' the `K - 1` estimated log ratios against the last component.
#'
#' @seealso [vcov.frmtmb_fit()] with `full = TRUE` for the same joint
#'   covariance (fixed effects plus covariance parameters, on their
#'   internal scale) as a matrix, which is what the `"wald"` method uses
#'   here.
#'
#' Methods:
#' - `"wald"` (default): delta-method z-test, finite-difference
#'   gradient against the joint parameter covariance (under REML, from
#'   the joint precision).
#' - `"profile"`: profile-likelihood interval via [TMB::tmbroot()] with
#'   a `lincomb` direction. Only for hypotheses that are linear in the
#'   parameters, and only for ML fits; the standard error and the test
#'   stay Wald-based, because the method changes the interval.
#' - `"boot"`: parametric bootstrap through [frm_bootstrap()]
#'   (percentile interval; `p` is the two-sided percentile p-value,
#'   whose resolution is limited by `nsim`; `Est.Error` is the bootstrap
#'   SD). Handles any expression, including the variance-component
#'   names, whose sampling distributions Wald approximates poorly.
#'
#' For a [frm_multiple()] result the Wald estimate and delta-method
#' variance are computed per imputation and pooled by Rubin's rules
#' with Barnard-Rubin degrees of freedom; the test attribute carries
#' `t` and `df` in place of `z`, and only Wald inference is available.
#'
#' @param x A `frmtmb_fit`, or a `frmtmb_multiple` for pooled tests.
#' @param hypothesis Character vector of hypotheses. Names on it become
#'   the `Hypothesis` labels, as in brms.
#' @param class,group brms's name prefix; see *brms class and group
#'   shorthand*. The default `class = "b"` is brms's.
#' @param scope brms's `"standard"` is the only scope a fit supports:
#'   `"ranef"` and `"coef"` evaluate the hypothesis on each group
#'   level's draws, which a maximum-likelihood fit does not have, and
#'   are refused by name.
#' @param alpha Test level; the reported interval covers `1 - alpha`
#'   for a two-sided row and `1 - 2 * alpha` for a directional one, as
#'   in brms.
#' @param robust brms's median-and-MAD switch. It summarizes draws and
#'   a fit has none, so `TRUE` is refused by name.
#' @param seed Optional seed for `method = "boot"`.
#' @param method `"wald"`, `"profile"`, or `"boot"`.
#' @param nsim Bootstrap draws for `method = "boot"`; all hypotheses
#'   share one bootstrap run.
#' @param vcov `method = "wald"` only: a covariance matrix over the
#'   whole outer parameter vector to use in place of the model-based
#'   one - [vcov_cluster()] with `full = TRUE`, or a function of the
#'   fit returning such a matrix. The delta-method standard error is
#'   then the cluster-robust one, and a matrix carrying reference
#'   degrees of freedom switches the test to a `t` reference.
#' @param ... Backend controls: passed to [TMB::tmbprofile()] for
#'   `method = "profile"` (e.g. `ytol`, `ystep`, `maxit`,
#'   `parm.range`) and to [frm_bootstrap()] for `method = "boot"`
#'   (e.g. `re_formula = NULL` for a conditional bootstrap). Refused for
#'   `"wald"`.
#' @return A `frmtmb_hypothesis` list in brms's shape; see
#'   *The returned object*.
#'   `plot()` shows the bootstrap distribution, the profile curve, or
#'   the implied Wald normal density, one panel per hypothesis.
#' @examples
#' set.seed(4)
#' dd <- data.frame(x1 = rnorm(120), x2 = rnorm(120),
#'                  g = factor(rep(1:10, 12)))
#' dd$y <- rnorm(120, 1 + 0.6 * dd$x1 + 0.4 * dd$x2 +
#'                 rnorm(10, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x1 + x2 + (1 | g)) + gaussian(), data = dd)
#' h <- hypothesis(fit, c("x1 - x2 = 0", "exp(Intercept) = 1"))
#' h
#' h$hypothesis$Est.Error
#' attr(h, "test")$p
#' # brms's directional form: one-sided p, and a 90% interval
#' hypothesis(fit, "x1 > x2")
#' # class/group name the natural-scale random-effect summaries
#' hypothesis(fit, "Intercept > 0", class = "sd", group = "g")
#' # variance-component expressions need class = NULL, as in brms: an
#' # ICC with bootstrap intervals
#' hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2) = 0",
#'            class = NULL, method = "boot", nsim = 20, seed = 1)
#' @export
hypothesis <- function(x, ...) UseMethod("hypothesis")

#' Usable parameter names
#'
#' brms's names for the parameters of a fit, which are the names
#' [hypothesis()] expressions accept: fixed-effect coefficients with
#' brms's `b_` prefix (`b_Intercept`, `b_x`, and `b_sigma_Intercept`
#' for a coefficient of a `sigma` formula), natural-scale random-effect
#' summaries (`sd_<group>__<coef>`, `cor_<group>__<coef1>__<coef2>`,
#' with a distributional or nonlinear parameter's name in the
#' coefficient part, `sd_g__sigma_Intercept`), and a distributional
#' parameter nobody wrote a formula for, on its natural scale (`sigma`,
#' `shape`, `sigma_ya`). Every name is spelled through brms's renaming:
#' `b_IxE2` for `I(x^2)`, `sd_g:h__Intercept` for `(1 | g:h)`. For sampled fits,
#' `variables()` on the `frmtmb.sample::frm_sample()` result lists the
#' draw columns, which follow the same convention.
#'
#' brms's `variables()` also lists what a fit has no counterpart of:
#' group-level coefficients `r_<group>[<level>,<coef>]`, the centered
#' `Intercept`, `lprior` and `lp__`. A maximum-likelihood fit has no
#' draws of those, and [ranef()] reports the conditional modes.
#'
#' A residual correlation term ([frmtmb-autocor]) contributes its
#' natural-scale parameters under brms's names: `ar[1]`, `ar[2]`,
#' `ma[1]`, `cosy`, `cortime__<t1>__<t2>`.
#'
#' `gr(cov = )`, `gr(prec = )` and `equalto()` blocks contribute
#' `sd_`/`cor_` names for their within-level covariance. Smooths,
#' `gp()`/`hsgp()`, `car()` and `spde()` blocks contribute none: their
#' parameters are not standard deviations. See the "Which random-effect
#' blocks contribute names" section of [hypothesis()].
#'
#' @param x A `frmtmb_fit` or `frmtmb_draws`.
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
#' @return A character vector.
#' @examples
#' dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#' dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' variables(fit)
#' @export
variables <- function(x, ...) UseMethod("variables")

#' @rdname variables
#' @exportS3Method posterior::variables
#' @export
variables.frmtmb_fit <- function(x, ...) {
  frm_check_dots(...)
  vo <- hyp_vals_only(x)
  names(hyp_env_vals(x, vo$vals, vo$comp))
}

#' The `Hypothesis` label brms writes for one hypothesis string:
#' whitespace removed, the two sides written `(lhs)-(rhs)`, the right
#' side dropped when it is `0`, then the sign and `0`
#' (`brms:::eval_hypothesis()`). A name given on the hypothesis vector
#' replaces it, as in brms. The bare-quantity spelling brms does not
#' have is labeled `(expr)`.
#'
#' @noRd
hyp_labels <- function(hypothesis) {
  nms <- names(hypothesis)
  vapply(seq_along(hypothesis), function(i) {
    if (length(nms) && !is.na(nms[i]) && nzchar(nms[i])) return(nms[i])
    h <- gsub("[ \t\r\n]", "", hypothesis[i])
    sign <- regmatches(h, regexpr("[<>]=?|=", h))
    if (!length(sign)) return(paste0("(", h, ")"))
    lr <- strsplit(h, "[<>]=?|=")[[1L]]
    lhs <- paste0("(", lr[1L], ")")
    if (length(lr) > 1L && lr[2L] != "0") {
      lhs <- paste0(lhs, "-(", lr[2L], ")")
    }
    paste(lhs, substr(sign, 1L, 1L), "0")
  }, "")
}

#' brms's hypothesis result, in brms's shape: a list with the
#' eight-column `hypothesis` frame, `samples`, `prior_samples`, `class`
#' and `alpha`, in brms's order (`brms:::combine_hlist()`). Anything a
#' method knows beyond brms's columns rides on attributes, so the list
#' and the frame stay exactly the shape a ported script indexes.
#'
#' `star` follows brms: for a two-sided row the interval excludes 0,
#' for a one-sided row `onesided_ok` says whether the claim passes at
#' `1 - alpha`.
#'
#' @noRd
hyp_brms_result <- function(labels, estimate, error, lower, upper,
                            evid_ratio, post_prob, dir, onesided_ok,
                            samples, prior_samples, prefix, alpha,
                            attrs = list()) {
  star <- ifelse(dir == "two.sided", !(lower <= 0 & 0 <= upper),
                 onesided_ok)
  star <- ifelse(!is.na(star) & star, "*", "")
  hs <- data.frame(Hypothesis = labels, Estimate = unname(estimate),
                   Est.Error = unname(error), CI.Lower = unname(lower),
                   CI.Upper = unname(upper),
                   Evid.Ratio = unname(evid_ratio),
                   Post.Prob = unname(post_prob), Star = star,
                   stringsAsFactors = FALSE)
  out <- list(hypothesis = hs, samples = samples,
              prior_samples = prior_samples,
              class = hyp_class_label(prefix), alpha = alpha)
  for (nm in names(attrs)) attr(out, nm) <- attrs[[nm]]
  attr(out, "direction") <- dir
  # NOT brms's class as well. The stated reason for carrying it was
  # that print() and plot() would dispatch differently without it, and
  # that was measured false: frmtmb exports both generics for
  # frmtmb_hypothesis, which comes first, so stripping brmshypothesis
  # leaves the printed output identical and plot() working
  # (dev/reviews/20260918-shapes.md). What a frmtmb object must not do
  # is answer is(x, "<brms class>") TRUE, which is the rule the port
  # ledger cites eighteen times.
  class(out) <- "frmtmb_hypothesis"
  out
}

#' A draws frame in brms's `samples` layout: one column per hypothesis,
#' named `H1`, `H2`, ...
#'
#' @noRd
hyp_samples_frame <- function(m, k) {
  if (is.null(m)) m <- matrix(numeric(0), 0L, k)
  m <- as.data.frame(unname(as.matrix(m)))
  names(m) <- paste0("H", seq_len(k))
  m
}

#' The brms arguments a method without posterior draws cannot honor,
#' refused by name rather than ignored.
#'
#' @noRd
hyp_refuse_draws_args <- function(scope, robust, what) {
  scope <- frm_match_arg(scope, c("standard", "ranef", "coef"))
  if (!identical(scope, "standard")) {
    frm_stop(what, " cannot honor scope = \"", scope, "\": brms evaluates ",
             "the hypothesis on each group level's draws from ", scope,
             "(summary = FALSE), and a maximum-likelihood fit has no draws ",
             "of a group level. Sample with frmtmb.sample::frm_sample() for ",
             "that, or write the level's own parameter out", call. = FALSE)
  }
  check_flag(robust, "robust")
  if (robust) {
    frm_stop(what, " cannot honor robust = TRUE: brms's robust summary is ",
             "the median and MAD of the draws, and a maximum-likelihood ",
             "fit reports one estimate and its standard error. Sample with ",
             "frmtmb.sample::frm_sample() for that", call. = FALSE)
  }
  invisible(scope)
}

#' @rdname hypothesis
#' @exportS3Method brms::hypothesis
#' @export
hypothesis.frmtmb_fit <- function(x, hypothesis, class = "b", group = "",
                                  scope = c("standard", "ranef", "coef"),
                                  alpha = 0.05, robust = FALSE,
                                  seed = NULL,
                                  method = c("wald", "profile", "boot"),
                                  nsim = 500, vcov = NULL, ...) {
  method <- frm_match_arg(method)
  hyp_refuse_draws_args(scope, robust, "hypothesis()")
  check_probability(alpha, "alpha")
  if (!is.null(vcov) && method != "wald") {
    frm_stop("hypothesis(vcov = ) applies to method = 'wald' only: ",
             "method = '", method, "' does not go through a covariance ",
             "matrix", call. = FALSE)
  }
  # Same rule as confint(): checked against what THIS method forwards
  # them to, and refused rather than warned about, because a warning
  # returns a table built as if the argument had not been given.
  frm_check_dots(..., .allow = switch(method,
    boot = names(formals(frm_bootstrap)),
    profile = names(formals(TMB::tmbprofile)),
    NULL))
  vo <- hyp_vals_only(x)
  known <- names(hyp_env_vals(x, vo$vals, vo$comp))
  prefix <- hyp_class_prefix(class, group)
  hp <- hyp_parse_all(hypothesis, known, class, group)
  exs <- hp$exprs
  labels <- hyp_labels(hypothesis)
  k_n <- length(exs)
  vals0 <- vapply(seq_along(exs), function(i) {
    val <- hyp_eval(x, exs[[i]], vo$vals, vo$comp)
    if (!is.numeric(val) || length(val) != 1L) {
      frm_stop("Hypothesis '", hypothesis[i], "' must evaluate to a single ",
               "number at the fitted estimates", call. = FALSE)
    }
    val
  }, numeric(1))
  # brms's interval: central 1 - alpha for "=", central 1 - 2 alpha for a
  # directional row, whose relevant end is then the one-sided bound
  lo_p <- ifelse(hp$dir == "two.sided", alpha / 2, alpha)
  na_col <- rep(NA_real_, k_n)
  no_prior <- hyp_samples_frame(matrix(NA_real_, 0L, k_n), k_n)

  if (method == "boot") {
    FUN <- function(ft) {
      w <- hyp_vals_only(ft)
      vapply(exs, function(ex) hyp_eval(ft, ex, w$vals, w$comp),
             numeric(1))
    }
    bs <- frm_bootstrap(x, FUN, nsim = nsim, seed = seed, ...)
    colnames(bs$t) <- hypothesis
    fin <- lapply(seq_len(k_n), function(i) {
      t_i <- bs$t[, i]
      t_i[is.finite(t_i)]
    })
    se <- vapply(fin, stats::sd, 1)
    lwr <- vapply(seq_len(k_n), function(i) {
      unname(stats::quantile(fin[[i]], lo_p[i]))
    }, 1)
    upr <- vapply(seq_len(k_n), function(i) {
      unname(stats::quantile(fin[[i]], 1 - lo_p[i]))
    }, 1)
    p <- vapply(seq_len(k_n), function(i) hyp_tail_p(fin[[i]], hp$dir[i]),
                1)
    samp <- hyp_samples_frame(bs$t, k_n)
    return(hyp_brms_result(
      labels, vals0, se, lwr, upr, na_col, na_col, hp$dir, p < alpha,
      samp, hyp_samples_frame(matrix(NA_real_, nrow(samp), k_n), k_n),
      prefix, alpha,
      list(method = method,
           test = data.frame(Hypothesis = labels, z = vals0 / se, p = p),
           draws = bs$t, nsim = nsim, converged = bs$converged)))
  }

  pc <- hyp_par_cov(x)
  qfun <- stats::qnorm
  pfun <- stats::pnorm
  stat_name <- "z"
  if (!is.null(vcov)) {
    rv <- resolve_vcov_arg(x, vcov, "hypothesis")
    if (is.null(pc$outer_pos)) {
      frm_stop("hypothesis(vcov = ) needs a plain maximum-likelihood fit ",
               "(the REML / profile branch reads the joint precision, ",
               "which a supplied covariance does not replace)",
               call. = FALSE)
    }
    pc$V <- rv$V[pc$outer_pos, pc$outer_pos, drop = FALSE]
    if (!is.null(rv$df)) {
      qfun <- function(p) stats::qt(p, rv$df)
      pfun <- function(q) stats::pt(q, rv$df)
      stat_name <- "t"
    }
  }
  profiles <- vector("list", k_n)
  se <- lwr <- upr <- stat <- p <- numeric(k_n)
  for (i in seq_len(k_n)) {
    ex <- exs[[i]]
    fn <- function(v) hyp_eval(x, ex, v, pc$comp)
    g <- hyp_fd_grad(fn, pc$vals)
    se[i] <- sqrt(max(0, drop(t(g) %*% pc$V %*% g)))
    dir <- hp$dir[i]
    wr <- hyp_wald_row(vals0[i], se[i], dir, alpha, qfun, pfun)
    stat[i] <- wr$stat
    p[i] <- wr$p
    q <- qfun(1 - lo_p[i])
    lwr[i] <- vals0[i] - q * se[i]
    upr[i] <- vals0[i] + q * se[i]
    if (method == "profile") {
      if (x$REML) {
        frm_stop("method = 'profile' requires an ML fit (REML integrates ",
                 "the fixed effects out of the outer problem)",
                 call. = FALSE)
      }
      if (isTRUE(x$control$profile)) {
        frm_stop("hypothesis(method = 'profile') needs a fit without ",
                 "frmtmb_control(profile = TRUE)", call. = FALSE)
      }
      g2 <- hyp_fd_grad(fn, pc$vals + 0.1 * (1 + abs(pc$vals)))
      if (max(abs(g - g2)) > 1e-4 * max(1, max(abs(g)))) {
        frm_stop("Hypothesis '", hypothesis[i], "' is not linear in the ",
                 "parameters; use method = 'boot'", call. = FALSE)
      }
      v <- numeric(pc$n_outer)
      v[pc$outer_pos] <- g
      const <- vals0[i] - sum(g * pc$vals)
      lev <- 1 - 2 * lo_p[i]
      if (lev <= 0) {
        frm_stop("A one-sided profile bound needs alpha below 0.5",
                 call. = FALSE)
      }
      pargs <- utils::modifyList(
        list(obj = x$obj, lincomb = v, trace = FALSE,
             ytol = 0.5 * stats::qchisq(lev, 1) + 1),
        list(...)
      )
      pr <- do.call(TMB::tmbprofile, pargs)
      ci <- stats::confint(pr, level = lev)
      pr[[1L]] <- pr[[1L]] + const
      profiles[[i]] <- pr
      lwr[i] <- unname(ci[1]) + const
      upr[i] <- unname(ci[2]) + const
    }
  }
  test <- data.frame(Hypothesis = labels, stat = stat, p = p)
  names(test)[2L] <- stat_name
  # a directional claim passes when its one-sided bound excludes 0,
  # which for a Wald interval is the same event as p < alpha
  ok <- ifelse(hp$dir == "greater", lwr > 0, upr < 0)
  extra <- list(method = method, test = test)
  if (method == "profile") {
    extra$profiles <- stats::setNames(profiles, labels)
  }
  hyp_brms_result(labels, vals0, se, lwr, upr, na_col, na_col, hp$dir, ok,
                  hyp_samples_frame(NULL, k_n), no_prior, prefix, alpha,
                  extra)
}

#' @export
print.frmtmb_hypothesis <- function(x, digits = 2, ...) {
  frm_check_dots(...)
  method <- attr(x, "method") %||% "posterior"
  cat("Hypothesis Tests for class ", x$class, ":\n", sep = "")
  hs <- x$hypothesis
  num <- vapply(hs, is.numeric, TRUE)
  hs[num] <- lapply(hs[num], round, digits)
  print(hs, quote = FALSE)
  pone <- (1 - x$alpha * 2) * 100
  ptwo <- (1 - x$alpha) * 100
  cat("---\n'CI': ", pone, "%-CI for one-sided and ", ptwo,
      "%-CI for two-sided hypotheses.\n", sep = "")
  if (identical(method, "posterior")) {
    cat("'*': For one-sided hypotheses, the posterior probability ",
        "exceeds ", ptwo, "%;\nfor two-sided hypotheses, the value ",
        "tested against lies outside the ", ptwo, "%-CI.\n", sep = "")
    cat("Posterior probabilities of point hypotheses assume equal ",
        "prior probabilities.\n", sep = "")
    return(invisible(x))
  }
  cat("Method: ", method, if (identical(method, "boot"))
        paste0(" (", attr(x, "nsim"), " bootstrap draws, ",
               sum(!attr(x, "converged")), " failed or not converged)"),
      ". Est.Error is the ", if (identical(method, "boot"))
        "bootstrap SD" else "standard error",
      "; Evid.Ratio and\nPost.Prob are NA, because a maximum-likelihood ",
      "fit has no posterior.\n", sep = "")
  cat("'*': For one-sided hypotheses, the one-sided test rejects at ",
      "level ", x$alpha, ";\nfor two-sided hypotheses, the value tested ",
      "against lies outside the ", ptwo, "%-CI.\n", sep = "")
  tst <- attr(x, "test")
  if (!is.null(tst)) {
    tst[-1L] <- lapply(tst[-1L], signif, 4)
    print(tst, row.names = FALSE)
  }
  invisible(x)
}

#' @export
plot.frmtmb_hypothesis <- function(x, ask = NULL, ...) {
  frm_check_dots(...)
  method <- attr(x, "method") %||% "posterior"
  alpha <- x$alpha %||% 0.05
  hs <- x$hypothesis
  n <- nrow(hs)
  ask <- ask %||% (n > 1L && grDevices::dev.interactive())
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask), add = TRUE)
  }
  mark <- function(i) {
    graphics::abline(v = hs$Estimate[i], lwd = 2)
    graphics::abline(v = c(hs$CI.Lower[i], hs$CI.Upper[i]), lty = 2)
    graphics::abline(v = 0, col = 2)
  }
  for (i in seq_len(n)) {
    h <- hs$Hypothesis[i]
    if (method %in% c("boot", "posterior")) {
      d <- x$samples[[i]]
      d <- d[is.finite(d)]
      graphics::hist(d, freq = FALSE, breaks = "FD", main = h,
                     xlab = if (method == "boot") "bootstrap value" else
                       "posterior value",
                     col = "gray90", border = "gray60")
      if (length(unique(d)) > 1L) {
        graphics::lines(stats::density(d), lwd = 2)
      }
      mark(i)
    } else if (method == "profile") {
      pr <- attr(x, "profiles")[[i]]
      dnll <- pr$value - min(pr$value, na.rm = TRUE)
      graphics::plot(pr[[1L]], dnll, type = "l", lwd = 2, main = h,
                     xlab = "value",
                     ylab = "profile neg. log-likelihood change")
      graphics::abline(h = 0.5 * stats::qchisq(1 - 2 * alpha, 1), lty = 3)
      mark(i)
    } else {
      est <- hs$Estimate[i]
      se <- hs$Est.Error[i]
      xs <- seq(est - 4 * se, est + 4 * se, length.out = 200)
      graphics::plot(xs, stats::dnorm(xs, est, se),
                     type = "l", lwd = 2, main = h, xlab = "value",
                     ylab = "Wald (normal) density")
      mark(i)
    }
  }
  invisible(x)
}
