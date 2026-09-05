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
      stop("Parameter name '", x[k], "' is ambiguous once parentheses ",
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
    g <- hyp_san(bk[["group_name"]])
    tn <- hyp_san(bk[["cnms"]])
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
    add(hyp_san(names(nat)[1L]),
        if (p >= 1L && p <= length(ac_pos)) ac_pos[p] else NA_integer_)
  }
  stats::setNames(pos, nms)
}

#' Bare nonlinear-parameter names, as positions in `outer_par_map()`.
#' An `nl` parameter declared `la ~ 1` has exactly one coefficient,
#' `la_(Intercept)`, so the bare `la` names it without ambiguity. That
#' is the one-to-one rule the `sd_`/`ar1` aliases already follow. A
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
    stop("Nonlinear parameter '", bad[1L], "' has more than one ",
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
    alias <- par_alias_index(fit)
    hit <- match(parm, names(alias))
    took <- which(is.na(idx) & !is.na(hit))
    idx[took] <- alias[hit[took]]
    if (length(took)) {
      message(what, "(): ",
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
      stop("'", b, "' is a natural-scale summary rather than a fitted ",
           "parameter, and it does not stand for a single internal one ",
           "here (a correlation of a wider us() block mixes several, ",
           "and a response-scale summary such as sigma is a transform ",
           "of one). Use hypothesis(fit, \"", b, "\", method = ",
           "'profile'), which profiles the combination itself, read ",
           "confint_varcorr() for natural-scale variance components, ",
           "or name the internal parameter: confint(fit) lists them",
           call. = FALSE)
    }
    stop("Unknown parameter(s) in ", what, "(parm =): ",
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
#'   [frm_bootstrap()] for `method = "boot"` (e.g. `re.form`).
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
  method <- match.arg(method)
  if (method == "Wald") method <- "wald"
  # A length-2 level makes a length-2 quantile, which then RECYCLES
  # against the parameter vector: rows 1 and 3 of one table came back at
  # 90% and rows 2 and 4 at 95%, with nothing in the output recording
  # it. A level outside (0, 1) produced a table of NaN or of Inf.
  check_probability(level, "level")
  check_count(nsim, "nsim", min = 1L)
  if (!is.null(parm) && !is.character(parm)) {
    stop("`parm` must be a character vector of parameter names, or NULL ",
         "for all of them, not ", arg_desc(parm), call. = FALSE)
  }
  if (!is.null(vcov) && method != "wald") {
    stop("confint(vcov = ) applies to method = 'wald' only: ",
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
    stop("confint(method = '", method, "') needs a fit without ",
         "frmtmb_control(profile = TRUE)", call. = FALSE)
  }
  if (is.null(parm)) {
    stop("`parm` is required for method = '", method, "'", call. = FALSE)
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
    warning("No interval for ", length(bad), " component",
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
    warning("Uninformative interval for ", length(wide), " component",
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
  vc <- tryCatch(as.data.frame(VarCorr(fit)), error = function(e) NULL)
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
#' Reports the optimizer's own verdict plus four checks that a converged
#' fit can still fail: non-finite standard errors, complete separation
#' in a binomial-type fit, predictor columns scaled far from one, and
#' variance components on the boundary of their parameter space
#' (lme4's `isSingular()`, read off the estimates rather than the
#' Hessian).
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
#'   of the covariance, boundary (singular) variance components,
#'   separation, and predictor scaling. `frm_allfit()` refits across
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
  stopifnot(inherits(fit, "frmtmb_fit"))
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
    max_grad = if (length(gr)) max(abs(gr)) else NA_real_,
    worst_grad = if (length(gr)) nm[which.max(abs(gr))] else NA_character_,
    pdHess = isTRUE(sdr_of(fit)$pdHess),
    bad_se = nm[!is.finite(se)],
    min_cov_eigenvalue = if (!is.null(ev) && length(ev)) min(ev),
    extreme_theta = sd_i[abs(th[sd_i]) > 8],
    separation = if (!is.null(ps)) diagnose_separation(fit, ps),
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
    if (length(out$bad_se)) {
      cat("Non-finite standard errors:",
          paste(out$bad_se, collapse = ", "), "\n")
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
    if (!is.null(out$predictor_scale)) {
      cat("Badly scaled predictors: ",
          paste(paste0(out$predictor_scale$column, " (sd ",
                       format(out$predictor_scale$sd, digits = 3), ")"),
                collapse = "; "),
          "\n  Rescale the column, or refit with ",
          "frmtmb_control(autoscale = TRUE).\n", sep = "")
    }
    clean <- out$convergence == 0 && out$pdHess && !length(out$bad_se) &&
      is.null(out$singular) && is.null(out$separation) &&
      is.null(out$predictor_scale) &&
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
    warning("anova(): the fixed effects of ", model_label(fits[[i]]),
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
    stop("anova() needs at least two frmtmb fits to compare", call. = FALSE)
  }
  reml <- vapply(fits, `[[`, TRUE, "REML")
  if (any(reml) && isTRUE(refit)) {
    labs <- vapply(fits[reml], model_label, "")
    message("anova(): refitting ", length(labs), " REML model",
            if (length(labs) != 1L) "s" else "", " with ML: ",
            paste(labs, collapse = "; "))
    fits[reml] <- lapply(fits[reml], anova_refit_ml)
    reml[] <- FALSE
  }
  if (any(reml)) {
    if (!all(reml)) {
      stop("anova() cannot mix REML and ML fits: their likelihoods are ",
           "for different quantities. Refit them all with the same ",
           "REML setting, or pass refit = TRUE to compare them as ML ",
           "fits", call. = FALSE)
    }
    if (!reml_comparable(fits)) {
      stop("REML likelihoods are comparable only between fits whose ",
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
    stop("anova() needs fits with the same number of observations (got ",
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
#' @param ... Unused.
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
  test <- match.arg(test)
  if (object$REML) {
    stop("drop1() compares fixed effects; refit with REML = FALSE",
         call. = FALSE)
  }
  if (length(object$spec$responses) > 1) {
    stop("drop1() is not supported for multivariate fits", call. = FALSE)
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
    stop("scope is not a subset of the term labels: ",
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
    stop("An update formula written as a delta does not say which ",
         "response it changes. Pass the complete mvbf() as `formula`",
         call. = FALSE)
  }
  # nlf() on mu makes the response formula's right-hand side a body too,
  # whether or not nl = TRUE was written
  nl_mu <- isTRUE(bform$nl) ||
    length(intersect(all.vars(reformulas::RHSForm(bform$formula)),
                     names(bform$nlforms %||% list()))) > 0L
  if (nl_mu) {
    stop("An update formula written as a delta cannot be applied to a ",
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
      stop("Give the updated data once: as `data`, or as brms's ",
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

#' Make a term or group name usable as a variable name in a hypothesis
#' expression: drop parentheses and every other character an R name
#' cannot carry.
#'
#' @noRd
hyp_san <- function(s) gsub("[^[:alnum:]_.]", "", gsub("[()]", "", s))

#' Parameter values without any covariance machinery (usable on refits
#' inside a bootstrap without triggering sdreport).
#'
#' @noRd
hyp_vals_only <- function(fit) {
  est <- fit$estimates
  bd <- est[["betad"]]
  if (length(fx <- fit$frame[["betad_fixed_idx"]])) bd <- bd[-fx]
  list(
    vals = c(est[["beta"]], bd,
       est[["theta"]], est[["thetaac"]], est[["thetar"]]),
    comp = c(rep("beta", length(est[["beta"]])), rep("betad", length(bd)),
             rep("theta", length(est[["theta"]])),
             rep("thetaac", length(est[["thetaac"]])),
             rep("thetar", length(est[["thetar"]])))
  )
}

#' Values plus joint covariance of (beta, estimated betad, theta,
#' thetaac, thetar). ML: straight from cov.fixed in opt$par order
#' (outer_pos maps back into the full outer vector for tmbroot
#' lincombs). REML: beta is integrated out, so the blocks come from the
#' joint precision.
#'
#' @noRd
hyp_par_cov <- function(fit) {
  comps <- c("beta", "betad", "theta", "thetaac", "thetar")
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
    Vall <- solve_joint_precision(Q, fit$cache)
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

# Shadowing note bookkeeping. hyp_env_vals() is rebuilt for every
# hypothesis, and once per finite-difference step inside the delta
# method, so the note has to be armed and deduplicated by the call that
# a user actually typed rather than emitted where it is detected.
hyp_shadow_state <- new.env(parent = emptyenv())

#' Arm the shadowing note for one user-level call and return the state
#' to restore afterwards (nested calls therefore stay one-shot too).
#'
#' @noRd
hyp_shadow_arm <- function() {
  old <- list(armed = hyp_shadow_state$armed, seen = hyp_shadow_state$seen)
  hyp_shadow_state$armed <- TRUE
  hyp_shadow_state$seen <- character(0)
  old
}

#' @noRd
hyp_shadow_disarm <- function(old) {
  hyp_shadow_state$armed <- old$armed
  hyp_shadow_state$seen <- old$seen
  invisible(NULL)
}

#' Report the natural-scale names a fixed-effect coefficient has taken
#' over, once per armed call and once per name.
#'
#' @noRd
hyp_shadow_note <- function(shadow) {
  if (!length(shadow) || !isTRUE(hyp_shadow_state$armed)) {
    return(invisible(NULL))
  }
  new <- setdiff(names(shadow), hyp_shadow_state$seen)
  if (!length(new)) return(invisible(NULL))
  hyp_shadow_state$seen <- c(hyp_shadow_state$seen, new)
  parts <- vapply(new, function(nm) {
    sh <- shadow[[nm]]
    tail <- if (isTRUE(sh$dot)) {
      paste0("; that quantity is available as '.", nm, "'")
    } else {
      paste0("; '.", nm, "' is a coefficient too, so read the shadowed ",
             "quantity from summary() or VarCorr() instead")
    }
    paste0("'", nm, "' as the coefficient of the model term of that ",
           "name, not ", sh$meaning, tail)
  }, character(1))
  message("hypothesis() reads ", paste0(parts, collapse = ", and "), ".")
}

#' Named list the hypothesis expressions are evaluated in: fixed
#' coefficients under their vcov() names (parentheses stripped),
#' natural-scale random-effect summaries (`sd_<group>__<term>`,
#' `cor_<group>__<t1>__<t2>`), and `sigma` when it is a scalar.
#'
#' A coefficient name wins any collision with a natural-scale name (the
#' v0.21 guard: a covariate literally named `sigma` must stay
#' addressable). The shadowed quantity is then registered under a
#' leading dot (`.sigma`, `.sd_g__Intercept`) and named in a one-time
#' message, so the other meaning is reachable rather than merely
#' documented.
#'
#' @noRd
hyp_env_vals <- function(fit, vals, comp) {
  env <- list()
  cf <- c(vals[comp == "beta"], vals[comp == "betad"])
  raw <- estimated_coef_names(fit)
  cn <- gsub("[()]", "", raw)
  for (i in seq_along(cn)) env[[cn[i]]] <- cf[i]
  # the internal spelling as an alias, so a name copied from confint()
  # or vcov() also resolves when it is backquoted in the expression.
  # variables() keeps listing the parenthesis-free names, which are the
  # ones an expression can carry unquoted
  for (i in seq_along(raw)) {
    if (raw[i] != cn[i] && is.null(env[[raw[i]]])) env[[raw[i]]] <- cf[i]
  }

  # every name in the environment at this point came from a coefficient,
  # so a later collision is a shadow rather than two natural-scale
  # summaries competing (the latter keeps the first writer, as before)
  coef_names <- names(env)
  shadow <- list()
  put <- function(nm, val, meaning) {
    if (is.null(env[[nm]])) {
      env[[nm]] <<- val
      return(invisible(NULL))
    }
    if (!nm %in% coef_names) return(invisible(NULL))
    dn <- paste0(".", nm)
    # the dot slot itself can be taken by a coefficient literally named
    # `.sigma`; the note must then not claim the quantity is reachable
    filed <- !dn %in% coef_names
    if (filed && is.null(env[[dn]])) env[[dn]] <<- val
    if (is.null(shadow[[nm]])) {
      shadow[[nm]] <<- list(meaning = meaning, dot = filed)
    }
    invisible(NULL)
  }

  th <- vals[comp == "theta"]
  for (bk in fit$frame[["re_blocks"]]) {
    # Excluded: the structures whose theta segment is not a set of
    # standard deviations and correlations at all. `smooth` carries one
    # inverse smoothing parameter, `gp`/`hsgp` a marginal sd plus
    # lengthscales, `car` an sd plus a mixing proportion, `spde` a
    # precision and an inverse range. Their summaries live in
    # confint_varcorr() under their own names.
    #
    # Included (since v0.29): gr_cov, gr_prec and equalto. Their
    # registry vcov() is the WITHIN-level covariance - a plain sd for a
    # scalar block, sds plus correlations for the correlated-slopes
    # Kronecker path - which is exactly the quantity brms names
    # sd_<group>__<term>. equalto contributes fixed constants (it
    # estimates nothing), so its names read as knowns.
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "car", "spde")) next
    V <- covstruct_registry[[bk[["covstruct"]]]]$vcov(th[bk[["theta_idx"]]], bk)
    tn <- hyp_san(bk[["cnms"]])
    g <- hyp_san(bk[["group_name"]])
    sds <- sqrt(diag(V))
    for (j in seq_along(sds)) {
      nm <- paste0("sd_", g, "__", tn[j])
      put(nm, sds[j], "the random-effect standard deviation")
    }
    if (nrow(V) > 1L) {
      C <- stats::cov2cor(V)
      for (j in seq_len(nrow(V) - 1L)) {
        for (k in seq(j + 1L, nrow(V))) {
          nm <- paste0("cor_", g, "__", tn[j], "__", tn[k])
          put(nm, C[j, k], "the random-effect correlation")
        }
      }
    }
  }

  # R-side residual correlation: brms's own names, sanitized the way
  # every other name here is (ar[1] -> ar1, cortime__1__2 unchanged)
  thac <- vals[comp == "thetaac"]
  for (ac in fit$frame[["autocor"]] %||% list()) {
    nat <- autocor_natural(thac[ac[["theta_idx"]]], ac)
    for (j in seq_along(nat)) {
      nm <- hyp_san(names(nat)[j])
      put(nm, unname(nat[j]), "the residual autocorrelation")
    }
  }

  if (length(fit$spec$responses) == 1L) {
    # put() keeps a covariate literally named `sigma` visible under that
    # name (the v0.21 guard) and files the residual SD under `.sigma`
    for (lp in fit$frame[["linpreds"]]) {
      if (lp[["dpar"]] != "sigma") next
      if (!is.null(lp[["constant"]])) {
        put("sigma", lp[["constant"]], "the residual standard deviation")
      } else if (ncol(lp[["X"]]) == 1L &&
                 identical(colnames(lp[["X"]]), "(Intercept)") &&
                 is.null(lp[["Z"]]) && lp[["par"]] == "betad") {
        tpl_len <- length(fit$frame[["par_template"]][["betad"]])
        rk <- match(lp[["idx"]], setdiff(seq_len(tpl_len),
                                    fit$frame[["betad_fixed_idx"]]))
        bd <- vals[comp == "betad"]
        if (!is.na(rk)) {
          put("sigma", lp[["link"]]$linkinv(bd[rk]),
              "the residual standard deviation")
        }
      }
    }
  }
  hyp_shadow_note(shadow)
  env
}

#' Turn one hypothesis string into the language object to evaluate plus
#' the alternative it asks for. An `"lhs = rhs"` hypothesis becomes the
#' difference of the two sides, so every hypothesis is then tested
#' against zero; brms's directional `"lhs > rhs"` and `"lhs < rhs"`
#' become the same difference with a one-sided alternative.
#'
#' @noRd
hyp_parse <- function(h) {
  ops <- unlist(gregexpr("[<>]", h))
  ops <- ops[ops > 0L]
  if (length(ops)) {
    if (length(ops) > 1L) {
      stop("A hypothesis has at most one '<' or '>': '", h, "'",
           call. = FALSE)
    }
    op <- substr(h, ops, ops)
    lhs <- substr(h, 1L, ops - 1L)
    # ">=" and "<=" read as ">" and "<": the boundary has probability
    # zero under every sampling distribution used here
    rhs <- sub("^[[:space:]]*=", "", substring(h, ops + 1L))
    if (grepl("=", lhs, fixed = TRUE) || grepl("=", rhs, fixed = TRUE)) {
      stop("A hypothesis is directional ('<', '>') or an equality ",
           "('='), not both: '", h, "'", call. = FALSE)
    }
    return(list(expr = str2lang(paste0("(", lhs, ") - (", rhs, ")")),
                dir = if (op == ">") "greater" else "less"))
  }
  eq <- strsplit(h, "=", fixed = TRUE)[[1L]]
  txt <- if (length(eq) == 2L) {
    paste0("(", eq[1L], ") - (", eq[2L], ")")
  } else if (length(eq) == 1L) {
    h
  } else {
    stop("A hypothesis has at most one '=': '", h, "'", call. = FALSE)
  }
  list(expr = str2lang(txt), dir = "two.sided")
}

#' The name prefix implied by brms's `class` and `group` shorthand.
#' Class `"b"` (brms's default) and no class both mean the plain
#' coefficient names; a class with a group is the `sd_<group>__` /
#' `cor_<group>__` naming of the natural-scale random-effect
#' summaries; a class without a group is a dpar prefix such as
#' `sigma_`.
#'
#' @noRd
hyp_class_prefix <- function(class = NULL, group = NULL) {
  if (is.null(class) || !nzchar(class) || identical(class, "b")) return("")
  if (!is.null(group) && nzchar(group)) {
    return(paste0(class, "_", group, "__"))
  }
  paste0(class, "_")
}

#' Rewrite the bare names of a hypothesis expression under a `class` /
#' `group` prefix. A name already written in full keeps its spelling,
#' and a name that is neither is left for the evaluator to report; but
#' a name that exists WITHOUT the prefix and not with it is refused
#' rather than silently tested under the wrong class, which is the one
#' way this shorthand can quietly answer a different question.
#'
#' @noRd
hyp_prefix_names <- function(ex, prefix, known) {
  if (!nzchar(prefix)) return(ex)
  rec <- function(e) {
    if (is.name(e)) {
      s <- as.character(e)
      cand <- paste0(prefix, s)
      if (cand %in% known) return(as.name(cand))
      if (s %in% known && !startsWith(s, prefix)) {
        stop("'", cand, "' is not a parameter of this model, while '",
             s, "' is: the `class`/`group` shorthand would be ignored ",
             "for it. Drop them, or correct them; variables() lists ",
             "every usable name", call. = FALSE)
      }
      return(e)
    }
    if (is.call(e)) {
      for (i in seq_along(e)[-1L]) {
        if (!identical(e[[i]], quote(expr = ))) e[[i]] <- rec(e[[i]])
      }
    }
    e
  }
  rec(ex)
}

#' Parse every hypothesis string of one call: the shared front end of
#' the `frmtmb_fit`, `frmtmb_draws` and `frmtmb_multiple` methods.
#' Returns the expressions and their alternatives.
#'
#' @noRd
hyp_parse_all <- function(hypothesis, known, class = NULL, group = NULL) {
  prefix <- hyp_class_prefix(class, group)
  ps <- lapply(hypothesis, hyp_parse)
  list(exprs = lapply(ps, function(p) {
         hyp_prefix_names(p$expr, prefix, known)
       }),
       dir = vapply(ps, function(p) p$dir, ""))
}

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

#' The names of a hypothesis environment that are worth listing: the
#' parenthesis-free spelling of every parameter. The internal spellings
#' are in the environment as aliases (they need backquotes in an
#' expression), and listing both would double the vocabulary a reader
#' has to scan.
#'
#' @noRd
hyp_public_names <- function(ev) grep("[()]", names(ev), invert = TRUE,
                                      value = TRUE)

#' Evaluate a parsed hypothesis at one parameter vector. A failure lists
#' the available names, because an unknown name is the usual cause.
#'
#' @noRd
hyp_eval <- function(fit, ex, vals, comp) {
  ev <- hyp_env_vals(fit, vals, comp)
  tryCatch(eval(ex, ev), error = function(e) {
    stop(conditionMessage(e), "\nAvailable names: ",
         paste(hyp_public_names(ev), collapse = ", "),
         "\n(the internal spellings of confint() work too, backquoted; ",
         "a natural-scale name a coefficient has taken over carries a ",
         "leading dot)", call. = FALSE)
  })
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
#' against zero. A hypothesis is `"expr"` (tested against 0),
#' `"expr = rhs"`, e.g. `"x1 - x2 = 0"` or `"exp(Intercept) = 1"`, or
#' brms's directional `"lhs > rhs"` / `"lhs < rhs"`.
#'
#' @section Directional hypotheses:
#' `"lhs > rhs"` and `"lhs < rhs"` test the same difference
#' `(lhs) - (rhs)` against zero with a one-sided alternative, so the
#' reported `p` is the one-sided tail probability and the interval is
#' one-sided at level `1 - alpha`: the unbounded end prints as `Inf`
#' or `-Inf`. `p` is `pnorm()` of the signed z statistic for `"wald"`
#' and, as in the two-sided case where `se`, `z` and `p` stay
#' Wald-based, for `"profile"` too - the profile changes the BOUND,
#' which is the matching endpoint of the two-sided `1 - 2 * alpha`
#' profile interval, and nothing else. For `"boot"` and the draws
#' method `p` is the tail proportion of the draws with the
#' `(1 + k) / (1 + n)` correction. Where brms reports the posterior
#' probability of the direction, this reports its frequentist
#' complement: small `p` is evidence for the stated direction.
#' `">="` and `"<="` read as `">"` and `"<"`.
#'
#' @section brms class and group shorthand:
#' `class` and `group` prefix the bare names in the hypothesis, so
#' `hypothesis(fit, "Intercept - age > 0", class = "sd",
#' group = "patient")` tests
#' `sd_patient__Intercept - sd_patient__age`. `class = "b"` (brms's
#' default) and `class = NULL` leave the names alone; a `class`
#' without a `group` prefixes `<class>_`, which is how a
#' distributional coefficient such as `sigma_Intercept` is named. A
#' name already written in full keeps its spelling, so the two can be
#' mixed; but a name that exists only WITHOUT the prefix is an error
#' rather than a test of the unprefixed parameter, so a wrong `class`
#' or `group` cannot quietly answer a different question.
#'
#' Available names: the fixed-effect coefficients under their `vcov()`
#' row names with parentheses stripped (`Intercept`, `x`,
#' `sigma_Intercept`, ...), natural-scale random-effect summaries
#' `sd_<group>__<term>` and `cor_<group>__<t1>__<t2>` (brms naming),
#' and `sigma` when the residual SD is a scalar. So an ICC is
#' `"sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)"`.
#' [variables()] lists every usable name for a fit. The internal
#' spelling that [confint()] and [vcov()] print is accepted as well,
#' backquoted because it carries parentheses:
#' `` "`(Intercept)` - x" ``. The traffic runs the other way too:
#' `confint(parm = )` and `profile(parm = )` take these names,
#' whenever one of them stands for a single internal parameter.
#'
#' @section Which random-effect blocks contribute names:
#' Every block whose covariance parameters ARE standard deviations and
#' correlations: the plain structures (`us`, `diag`, `homdiag`, `cs`,
#' `ar1`, `toep`, the spatial and reduced-rank ones) and the
#' known-structure blocks `gr(cov = )`, `gr(prec = )` and `equalto()`,
#' whose `sd_`/`cor_` names describe the WITHIN-level covariance that
#' multiplies the fixed relationship matrix. That is what makes
#' heritability-as-ICC writable directly:
#' `"sd_id__Intercept^2 / (sd_id__Intercept^2 + sigma^2)"` on an animal
#' model fitted with `(1 | gr(id, cov = A))`. An `equalto()` block
#' estimates nothing, so its names are constants with zero variance.
#'
#' An `|ID|`-merged block is ONE block, so it contributes one name per
#' merged coefficient, and the names carry the linear predictor they
#' came from just as correlated slopes do. A two-trait animal model
#' written `(1 | q | gr(id, cov = A))` in both formulas of an [mvbf()]
#' gives `sd_id__y1.muIntercept`, `sd_id__y2.muIntercept` and
#' `cor_id__y1.muIntercept__y2.muIntercept` - the last being the
#' genetic correlation between the traits. [variables()] prints them.
#'
#' Two blocks on the same grouping factor with the same term name - an
#' animal model's `(1 | gr(id, cov = A)) + (1 | id)`, where the genetic
#' and permanent-environment terms both name the group `id` - collide
#' on one `sd_id__Intercept`, and the first block in formula order
#' claims it. Give the second term its own grouping column (a copy of
#' the factor under another name) when both are wanted by name.
#'
#' Excluded: `s()`/`t2()` smooths, `gp()`/`hsgp()`, `car()` and `spde()`.
#' Their theta segments are not standard deviations - an inverse
#' smoothing parameter, lengthscales, a mixing proportion, a precision
#' and an inverse range - so there is no `sd_<group>__<term>` to name.
#' Read those off [confint_varcorr()], which reports each under its own
#' label (`sd(gp)`, `range(gp)`, `sd(car)`, ...).
#'
#' @section When a coefficient shadows a natural-scale name:
#' Model terms and natural-scale summaries share one namespace here, and
#' the coefficient wins: a covariate literally named `sigma` makes
#' `"sigma = 0"` a test on ITS coefficient, not on the residual standard
#' deviation. The same holds for a coefficient that spells out
#' `sd_<group>__<term>`, `cor_...` or an autocorrelation name such as
#' `ar1`. The shadowed quantity keeps a name: prefix it with a dot,
#' `.sigma`, `.sd_g__Intercept`, `.ar1`. The dot spelling exists only
#' where a collision does, and `hypothesis()` says so once per call when
#' one is in play. `variables()` lists both names in that case.
#'
#' @seealso [vcov.frmtmb_fit()] with `full = TRUE` for the same joint covariance
#'   (fixed effects plus covariance parameters, on their internal
#'   scale) as a matrix, which is what the `"wald"` method uses here.
#'
#' Methods:
#' - `"wald"` (default): delta-method z-test, finite-difference
#'   gradient against the joint parameter covariance (under REML, from
#'   the joint precision).
#' - `"profile"`: profile-likelihood interval via [TMB::tmbroot()] with
#'   a `lincomb` direction. Only for hypotheses that are linear in the
#'   parameters, and only for ML fits; `se`, `z`, and `p` stay
#'   Wald-based - the method changes the interval.
#' - `"boot"`: parametric bootstrap through [frm_bootstrap()]
#'   (percentile interval; `p` is the two-sided percentile p-value,
#'   whose resolution is limited by `nsim`; `se` is the bootstrap SD).
#'   Handles any expression, including the variance-component names,
#'   whose sampling distributions Wald approximates poorly.
#'
#' For a [frm_multiple()] result the Wald estimate and delta-method
#' variance are computed per imputation and pooled by Rubin's rules
#' with Barnard-Rubin degrees of freedom; the returned table carries
#' `t` and `df` columns in place of `z` (reference t distribution, not
#' normal), and only Wald inference is available.
#'
#' @param x A `frmtmb_fit`, or a `frmtmb_multiple` for pooled tests.
#' @param hypothesis Character vector of hypotheses.
#' @param alpha Test level; the reported interval covers `1 - alpha`
#'   (brms spelling).
#' @param method `"wald"`, `"profile"`, or `"boot"`.
#' @param nsim Bootstrap draws for `method = "boot"`; all hypotheses
#'   share one bootstrap run.
#' @param seed Optional seed for `method = "boot"`.
#' @param class,group brms shorthand for the parameter names: the
#'   hypothesis is written with bare names and `class` (and `group`,
#'   for the `sd_`/`cor_` summaries) supplies the prefix. The default
#'   `NULL` (like brms's `class = "b"`) takes the names as written.
#' @param vcov `method = "wald"` only: a covariance matrix over the
#'   whole outer parameter vector to use in place of the model-based
#'   one - [vcov_cluster()] with `full = TRUE`, or a function of the
#'   fit returning such a matrix. The delta-method standard error is
#'   then the cluster-robust one, and a matrix carrying reference
#'   degrees of freedom switches the test to a `t` reference.
#' @param ... Backend controls: passed to [TMB::tmbprofile()] for
#'   `method = "profile"` (e.g. `ytol`, `ystep`, `maxit`,
#'   `parm.range`) and to [frm_bootstrap()] for `method = "boot"`
#'   (e.g. `re.form = NULL` for a conditional bootstrap). Unused for
#'   `"wald"` (a warning).
#' @return A `frmtmb_hypothesis` object: a data frame with one row per
#'   hypothesis (`estimate`, `se`, `lwr`, `upr`, `z`, `p`) carrying the
#'   method payload in attributes - the bootstrap draws matrix
#'   (`attr(., "draws")`) or the profile curves (`attr(., "profiles")`).
#'   `plot()` shows the bootstrap distribution, the profile curve, or
#'   the implied Wald normal density, one panel per hypothesis.
#'   Subsetting with `[` drops the attributes; keep the full object for
#'   plotting.
#' @examples
#' set.seed(4)
#' dd <- data.frame(x1 = rnorm(120), x2 = rnorm(120),
#'                  g = factor(rep(1:10, 12)))
#' dd$y <- rnorm(120, 1 + 0.6 * dd$x1 + 0.4 * dd$x2 +
#'                 rnorm(10, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x1 + x2 + (1 | g)) + gaussian(), data = dd)
#' hypothesis(fit, c("x1 - x2 = 0", "exp(Intercept)"))
#' # brms's directional form: one-sided p, one-sided interval
#' hypothesis(fit, "x1 > x2")
#' # class/group name the natural-scale random-effect summaries
#' hypothesis(fit, "Intercept > 0", class = "sd", group = "g")
#' # variance-component expressions: an ICC with bootstrap intervals
#' hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)",
#'            method = "boot", nsim = 20, seed = 1)
#' @export
hypothesis <- function(x, ...) {
  # the shadowing note belongs to the call the user typed, not to any of
  # the many environment rebuilds it triggers, so it is armed here and
  # restored when the method returns (on.exit survives UseMethod)
  old <- hyp_shadow_arm()
  on.exit(hyp_shadow_disarm(old), add = TRUE)
  UseMethod("hypothesis")
}

#' Usable parameter names
#'
#' The names that [hypothesis()] expressions (and `set_prior()`
#' targeting) accept: fixed-effect coefficients under their `vcov()`
#' names with parentheses stripped, natural-scale random-effect
#' summaries (`sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`), and
#' `sigma` when the residual SD is a scalar. The brms spelling; for
#' sampled fits, `variables()` on the `frmtmb.sample::frm_sample()`
#' result lists the draw columns instead.
#'
#' A residual correlation term ([frmtmb-autocor]) contributes its
#' natural-scale parameters under brms's names, sanitized the same way:
#' `ar1`, `ar2`, `ma1`, `cosy`, `cortime__<t1>__<t2>`.
#'
#' `gr(cov = )`, `gr(prec = )` and `equalto()` blocks contribute
#' `sd_`/`cor_` names for their within-level covariance. Smooths,
#' `gp()`/`hsgp()`, `car()` and `spde()` blocks contribute none: their
#' parameters are not standard deviations. See the "Which random-effect
#' blocks contribute names" section of [hypothesis()].
#'
#' @param x A `frmtmb_fit` or `frmtmb_draws`.
#' @param ... Unused.
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
  vo <- hyp_vals_only(x)
  hyp_public_names(hyp_env_vals(x, vo$vals, vo$comp))
}

#' @rdname hypothesis
#' @exportS3Method brms::hypothesis
#' @export
hypothesis.frmtmb_fit <- function(x, hypothesis, alpha = 0.05,
                                  method = c("wald", "profile", "boot"),
                                  nsim = 500, seed = NULL, class = NULL,
                                  group = NULL, vcov = NULL, ...) {
  method <- match.arg(method)
  if (!is.null(vcov) && method != "wald") {
    stop("hypothesis(vcov = ) applies to method = 'wald' only: ",
         "method = '", method, "' does not go through a covariance ",
         "matrix", call. = FALSE)
  }
  if (method == "wald" && ...length()) {
    warning("ignoring arguments unused by method = 'wald': ",
            paste(...names(), collapse = ", "), call. = FALSE)
  }
  vo <- hyp_vals_only(x)
  known <- names(hyp_env_vals(x, vo$vals, vo$comp))
  hp <- hyp_parse_all(hypothesis, known, class, group)
  exs <- hp$exprs
  vals0 <- vapply(seq_along(exs), function(i) {
    val <- hyp_eval(x, exs[[i]], vo$vals, vo$comp)
    if (!is.numeric(val) || length(val) != 1L) {
      stop("Hypothesis '", hypothesis[i], "' must evaluate to a single ",
           "number at the fitted estimates", call. = FALSE)
    }
    val
  }, numeric(1))

  hyp_result <- function(out, extra = list()) {
    rownames(out) <- NULL
    attr(out, "method") <- method
    attr(out, "alpha") <- alpha
    attr(out, "direction") <- hp$dir
    for (nm in names(extra)) attr(out, nm) <- extra[[nm]]
    class(out) <- c("frmtmb_hypothesis", "data.frame")
    out
  }

  if (method == "boot") {
    FUN <- function(ft) {
      w <- hyp_vals_only(ft)
      vapply(exs, function(ex) hyp_eval(ft, ex, w$vals, w$comp),
             numeric(1))
    }
    bs <- frm_bootstrap(x, FUN, nsim = nsim, seed = seed, ...)
    colnames(bs$t) <- hypothesis
    rows <- lapply(seq_along(exs), function(i) {
      t_i <- bs$t[, i]
      t_i <- t_i[is.finite(t_i)]
      se <- stats::sd(t_i)
      dir <- hp$dir[i]
      lo <- if (dir == "two.sided") alpha / 2 else alpha
      data.frame(hypothesis = hypothesis[i], estimate = vals0[i],
                 se = se,
                 lwr = if (dir == "less") -Inf else
                   unname(stats::quantile(t_i, lo)),
                 upr = if (dir == "greater") Inf else
                   unname(stats::quantile(t_i, 1 - lo)),
                 z = vals0[i] / se, p = hyp_tail_p(t_i, dir))
    })
    return(hyp_result(do.call(rbind, rows),
                      list(draws = bs$t, nsim = nsim,
                           converged = bs$converged)))
  }

  pc <- hyp_par_cov(x)
  qfun <- stats::qnorm
  pfun <- stats::pnorm
  if (!is.null(vcov)) {
    rv <- resolve_vcov_arg(x, vcov, "hypothesis")
    if (is.null(pc$outer_pos)) {
      stop("hypothesis(vcov = ) needs a plain maximum-likelihood fit ",
           "(the REML / profile branch reads the joint precision, ",
           "which a supplied covariance does not replace)",
           call. = FALSE)
    }
    pc$V <- rv$V[pc$outer_pos, pc$outer_pos, drop = FALSE]
    if (!is.null(rv$df)) {
      qfun <- function(p) stats::qt(p, rv$df)
      pfun <- function(q) stats::pt(q, rv$df)
    }
  }
  profiles <- vector("list", length(exs))
  rows <- vector("list", length(exs))
  for (i in seq_along(exs)) {
    ex <- exs[[i]]
    fn <- function(v) hyp_eval(x, ex, v, pc$comp)
    g <- hyp_fd_grad(fn, pc$vals)
    se <- sqrt(max(0, drop(t(g) %*% pc$V %*% g)))
    dir <- hp$dir[i]
    wr <- hyp_wald_row(vals0[i], se, dir, alpha, qfun, pfun)
    zv <- wr$stat
    lwr <- wr$lwr
    upr <- wr$upr
    if (method == "profile") {
      if (x$REML) {
        stop("method = 'profile' requires an ML fit (REML integrates ",
             "the fixed effects out of the outer problem)",
             call. = FALSE)
      }
      if (isTRUE(x$control$profile)) {
        stop("hypothesis(method = 'profile') needs a fit without ",
             "frmtmb_control(profile = TRUE)", call. = FALSE)
      }
      g2 <- hyp_fd_grad(fn, pc$vals + 0.1 * (1 + abs(pc$vals)))
      if (max(abs(g - g2)) > 1e-4 * max(1, max(abs(g)))) {
        stop("Hypothesis '", hypothesis[i], "' is not linear in the ",
             "parameters; use method = 'boot'", call. = FALSE)
      }
      v <- numeric(pc$n_outer)
      v[pc$outer_pos] <- g
      const <- vals0[i] - sum(g * pc$vals)
      # a one-sided bound at level 1 - alpha is the matching endpoint of
      # the two-sided 1 - 2 * alpha profile interval
      lev <- if (dir == "two.sided") 1 - alpha else 1 - 2 * alpha
      if (lev <= 0) {
        stop("A one-sided profile bound needs alpha below 0.5",
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
      lwr <- if (dir == "less") -Inf else unname(ci[1]) + const
      upr <- if (dir == "greater") Inf else unname(ci[2]) + const
    }
    rows[[i]] <- data.frame(hypothesis = hypothesis[i],
                            estimate = vals0[i], se = se,
                            lwr = lwr, upr = upr, z = zv,
                            p = wr$p)
  }
  hyp_result(do.call(rbind, rows),
             if (method == "profile") {
               list(profiles = stats::setNames(profiles, hypothesis))
             } else {
               list()
             })
}

#' @export
print.frmtmb_hypothesis <- function(x, digits = 4, ...) {
  method <- attr(x, "method")
  cat("Hypothesis tests (method = ", method, ")\n", sep = "")
  if (identical(method, "boot")) {
    cat("  bootstrap draws: ", attr(x, "nsim"), " (",
        sum(!attr(x, "converged")), " failed or not converged)\n",
        sep = "")
  }
  df <- x
  class(df) <- "data.frame"
  df[-1] <- lapply(df[-1], signif, digits)
  print(df, row.names = FALSE)
  dir <- attr(x, "direction") %||% rep("two.sided", nrow(x))
  if (any(dir != "two.sided")) {
    cat("  rows written with '<' or '>' are one-sided: p and the ",
        "finite interval\n  bound hold at level ",
        signif(1 - (attr(x, "alpha") %||% 0.05), 3), "\n", sep = "")
  }
  invisible(x)
}

#' @export
plot.frmtmb_hypothesis <- function(x, ask = NULL, ...) {
  method <- attr(x, "method") %||% "wald"
  alpha <- attr(x, "alpha") %||% 0.05
  n <- nrow(x)
  ask <- ask %||% (n > 1L && grDevices::dev.interactive())
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask), add = TRUE)
  }
  mark <- function(i) {
    graphics::abline(v = x$estimate[i], lwd = 2)
    graphics::abline(v = c(x$lwr[i], x$upr[i]), lty = 2)
    graphics::abline(v = 0, col = 2)
  }
  for (i in seq_len(n)) {
    h <- x$hypothesis[i]
    if (method %in% c("boot", "posterior")) {
      d <- attr(x, "draws")[, i]
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
      graphics::abline(h = 0.5 * stats::qchisq(1 - alpha, 1), lty = 3)
      mark(i)
    } else {
      xs <- seq(x$estimate[i] - 4 * x$se[i], x$estimate[i] + 4 * x$se[i],
                length.out = 200)
      graphics::plot(xs, stats::dnorm(xs, x$estimate[i], x$se[i]),
                     type = "l", lwd = 2, main = h, xlab = "value",
                     ylab = "Wald (normal) density")
      mark(i)
    }
  }
  invisible(x)
}
