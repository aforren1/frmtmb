# Multiply-imputed data: fit per imputation, pool by Rubin's rules.

#' Rubin's rules with Barnard-Rubin small-sample df. `Q` and `U` are
#' quantity-by-imputation matrices of estimates and squared standard
#' errors; `dfcom` is the complete-data residual df.
#'
#' @noRd
rubin_pool <- function(Q, U, dfcom) {
  m <- ncol(Q)
  qbar <- rowMeans(Q)
  ubar <- rowMeans(U)
  bvar <- apply(Q, 1, stats::var)
  tvar <- ubar + (1 + 1 / m) * bvar
  lam <- (1 + 1 / m) * bvar / tvar
  df_old <- (m - 1) / pmax(lam, 1e-12)^2
  df_obs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lam)
  data.frame(estimate = qbar, se = sqrt(tvar),
             df = df_old * df_obs / (df_old + df_obs), fmi = lam)
}

#' Row names of the pooled coefficient table: brms's, as `fixef()` and
#' `vcov()` of each fit name the same coefficients (`Intercept`,
#' `sigma_x`). A distributional parameter nobody wrote a formula for is
#' pooled here on its link scale, where brms reports the parameter
#' itself under its bare name, so it keeps its predictor's name,
#' `sigma_Intercept`, which says what the number is.
#'
#' @noRd
multiple_row_names <- function(fit) {
  tab <- brms_coef_table(fit)
  nm <- sub("^(b|bs|bsp|bcs)_", "", tab$brms)
  nat <- tab$natural
  # brms's name of the intercept that predictor would have, had a
  # formula been written for it (sigma_Intercept, sigma_bmi_Intercept)
  nm[nat] <- paste0(tab$brms[nat], "_Intercept")
  nm
}

#' Fit a model across multiply imputed datasets
#'
#' The frequentist counterpart of brms's `brm_multiple()`: fits the
#' model on every imputed dataset and pools by Rubin's rules with
#' Barnard-Rubin adjusted degrees of freedom. Accepts a plain list of
#' data frames or a `mice::mids` object.
#'
#' Fixed effects are pooled on the link scale, so distributional
#' coefficients like `sigma` appear in `$pooled` on their link (log)
#' scale, exactly as in `vcov()`. Random-effect SDs and correlations
#' are pooled on the scales where a Wald argument is defensible (log
#' for SDs and GP ranges, Fisher z for correlations, the
#' [confint_varcorr()] convention) and back-transformed, giving
#' `$pooled_varcorr`. For missing-data mechanisms beyond imputation,
#' see in-model `mi()`.
#'
#' [hypothesis()] pools arbitrary functions of the parameters by the
#' same rules, and [anova.frmtmb_multiple()] compares two nested
#' `frm_multiple()` fits with the D1, D2 and D3 rules of `mice`.
#' `frmtmb.sample::frm_sample()` samples the fit of every imputation
#' and combines the draws, as brms's `brm_multiple()` does.
#'
#' @param formula,... As in [frm()].
#' @param data A list of completed data frames, or a `mice::mids`.
#' @param level Confidence level for the `$pooled_varcorr` interval.
#' @return A `frmtmb_multiple` object: `pooled` (the Rubin table for
#'   the fixed effects, its rows named as brms names them, `Intercept`
#'   and `sigma_Intercept`), `pooled_varcorr` (grp/term/type/estimate/
#'   lwr/upr/df/fmi for the random-effect SDs and correlations; `NULL`
#'   without random effects), and `fits` (the per-imputation fits).
#'   `fixef()`, `summary()`, `vcov()`, `fitted()`, `conditional_effects()`
#'   and the other methods of [frm_multiple-methods] pool across the
#'   imputations.
#' @examples
#' set.seed(8)
#' n <- 80
#' x <- rnorm(n)
#' y <- rnorm(n, 1 + 0.5 * x, 1)
#' x[sample(n, 15)] <- NA
#' imps <- lapply(1:3, function(i) {
#'   xi <- x
#'   xi[is.na(xi)] <- sample(x[!is.na(x)], sum(is.na(xi)), TRUE)
#'   data.frame(y = y, x = xi)
#' })
#' frm_multiple(bf(y ~ x) + gaussian(), data = imps)
#' @srrstats {G2.1} `data` is asserted to be a list of at least two data
#'   frames, checked element by element, and errors with that requirement
#'   spelled out. A `mice::mids` object is converted to that form first,
#'   after checking that mice is installed.
#' @srrstats {G2.14c} This is the second route for replacing missing data
#'   with imputed values: imputations produced outside the model (by
#'   mice, or by any procedure that yields a list of completed data
#'   frames) are each fitted and then pooled by Rubin's rules, covering
#'   the coefficients, the variance components, and `hypothesis()` and
#'   `anova()` on the pooled fit. `bf(x | mi() ~ ...)` is the in-model
#'   alternative.
#'
#' @export
frm_multiple <- function(formula, data, level = 0.95, ...) {
  # the pooled interval is built from one normal quantile; a length-2
  # level recycles it down the coefficient table
  check_probability(level, "level")
  if (inherits(data, "mids")) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      frm_stop("A mids object needs the 'mice' package", call. = FALSE)
    }
    data <- mice::complete(data, action = "all")
  }
  if (!is.list(data) || length(data) < 2L ||
      !all(vapply(data, is.data.frame, TRUE))) {
    frm_stop("data must be a list of at least two data frames (or a ",
             "mice::mids object)", call. = FALSE)
  }
  fits <- lapply(data, function(d) frm(formula, data = d, ...))
  m <- length(fits)
  dfcom <- df.residual(fits[[1]])

  nm <- estimated_coef_names(fits[[1]])
  for (j in seq_len(m)[-1L]) {
    # vapply() below reported this as "values must be length 2"
    if (!identical(estimated_coef_names(fits[[j]]), nm)) {
      frm_stop("frm_multiple(): imputation ", j, " has other coefficients ",
               "than the first (", paste(setdiff(union(
                 estimated_coef_names(fits[[j]]), nm), intersect(
                   estimated_coef_names(fits[[j]]), nm)), collapse = ", "),
               "), so there is nothing to pool: a factor level, say, is ",
               "missing from one of the imputed data sets", call. = FALSE)
    }
  }
  cf <- vapply(fits, function(f) {
    bd <- f$estimates[["betad"]]
    if (length(fx <- f$frame[["betad_fixed_idx"]])) bd <- bd[-fx]
    c(f$estimates[["beta"]], bd)
  }, numeric(length(nm)))
  cf <- matrix(cf, nrow = length(nm), dimnames = list(nm, NULL))
  # vcov_estimated(), not vcov(): the pooling runs over EVERY estimated
  # coefficient, which is what `nm` and `cf` above are, and vcov() is
  # brms's population-level block, one row shorter on a fit with an
  # intercept-only dispersion parameter
  us <- vapply(fits, function(f) diag(vcov_estimated(f)),
               numeric(length(nm)))
  pl <- rubin_pool(cf, matrix(us, nrow = length(nm)), dfcom)
  tstat <- pl$estimate / pl$se
  # brms's names, which fixef() of each fit uses: `Intercept`, not the
  # design column's `(Intercept)`
  tab <- data.frame(
    estimate = pl$estimate, se = pl$se, df = pl$df, t = tstat,
    p = 2 * stats::pt(-abs(tstat), pl$df), fmi = pl$fmi,
    row.names = multiple_row_names(fits[[1]])
  )

  # variance components, pooled on the log/atanh scales and
  # back-transformed (interval via the Barnard-Rubin t quantile)
  vc <- lapply(fits, varcorr_trans_rows)
  pooled_varcorr <- NULL
  if (!is.null(vc[[1]])) {
    nq <- nrow(vc[[1]])
    Q <- matrix(vapply(vc, function(v) v$est_t, numeric(nq)), nrow = nq)
    U <- matrix(vapply(vc, function(v) v$se_t^2, numeric(nq)), nrow = nq)
    pv <- rubin_pool(Q, U, dfcom)
    tq <- stats::qt(1 - (1 - level) / 2, pv$df)
    pooled_varcorr <- data.frame(
      grp = vc[[1]]$block, term = vc[[1]]$term, type = vc[[1]]$type,
      estimate = varcorr_untrans(vc[[1]]$type, pv$estimate),
      lwr = varcorr_untrans(vc[[1]]$type, pv$estimate - tq * pv$se),
      upr = varcorr_untrans(vc[[1]]$type, pv$estimate + tq * pv$se),
      df = pv$df, fmi = pv$fmi
    )
  }
  structure(list(pooled = tab, pooled_varcorr = pooled_varcorr,
                 fits = fits, m = m, level = level),
            class = "frmtmb_multiple")
}

#' @export
print.frmtmb_multiple <- function(x, digits = 4, ...) {
  frm_check_dots(...)
  cat("Pooled over", x$m, "imputations (Rubin's rules):\n\n")
  print(signif(x$pooled, digits))
  if (!is.null(x$pooled_varcorr)) {
    cat("\nPooled variance components (log/Fisher-z scale, ",
        "back-transformed):\n\n", sep = "")
    vc <- x$pooled_varcorr
    vc[-(1:3)] <- lapply(vc[-(1:3)], signif, digits)
    print(vc, row.names = FALSE)
  }
  invisible(x)
}

#' Pooled model comparison across imputations (D1, D2, D3)
#'
#' The multiply-imputed counterpart of [anova.frmtmb_fit()]. A
#' likelihood-ratio test per imputation is not a test: the m statistics
#' have to be combined, and the three standard combining rules
#' (`mice::D1()`, `mice::D2()`, `mice::D3()`) each reference an F
#' distribution whose denominator degrees of freedom absorb the extra
#' between-imputation uncertainty.
#'
#' \describe{
#'   \item{`"D1"`}{Multivariate Wald (Li, Raghunathan and Rubin 1991).
#'     Let `Qbar` be the pooled estimate of the `k` coefficients that
#'     the larger model adds, `Ubar` their average within-imputation
#'     covariance and `B` the between-imputation covariance. With the
#'     average relative increase in variance r, which is
#'     `(1 + 1/m) tr(B Ubar^-1) / k`, the statistic is
#'     `Qbar' ((1 + r) Ubar)^-1 Qbar / k` on `k` and `v` degrees of
#'     freedom. `v` uses the Reiter (2007) small-sample form when
#'     `dfcom` is finite and the Li et al form otherwise.}
#'   \item{`"D2"`}{Pooled test statistics (Li, Meng, Raghunathan and
#'     Rubin 1991). From the `m` per-imputation statistics `d` - either
#'     likelihood-ratio chi-squares (`use = "likelihood"`) or Wald
#'     chi-squares (`use = "wald"`) - and r set to
#'     `(1 + 1/m) var(sqrt(d))`, the statistic is
#'     `(mean(d)/k - (m + 1)/(m - 1) r) / (1 + r)` on `k` and
#'     `k^(-3/m) (m - 1) (1 + 1/r)^2` degrees of freedom. D2 needs no
#'     covariance matrix, so it is the cheapest rule and the least
#'     efficient.}
#'   \item{`"D3"`}{Pooled likelihood ratio (Meng and Rubin 1992). The
#'     average deviance difference is computed twice: at the
#'     imputation-specific estimates (`dbar`) and at the pooled
#'     estimates (`dtilde`), the second by evaluating each imputation's
#'     own likelihood at one common parameter vector. With r set to
#'     `(m + 1) / (k (m - 1)) (dbar - dtilde)`, the statistic is
#'     `dtilde / (k (1 + r))`. D3 is the default: it needs no
#'     covariance matrix either, and unlike D2 it stays valid when the
#'     coefficients are far from normal.}
#' }
#'
#' Pooled parameter values are the plain across-imputation means of the
#' optimizer's own parameter vector, so variance components are pooled
#' on the internal (log / Cholesky) scale, matching `$pooled`. No model
#' is refit: each imputation's objective is re-evaluated at the pooled
#' vector, which costs one Laplace solve.
#'
#' `mice::D3()` instead fixes the pooled coefficients as an offset and
#' re-estimates every remaining parameter, so it can differ by a few
#' percent; `mitml::testModels(method = "D3")` uses the plug-in form
#' implemented here. Note also that `df.residual()` counts the
#' dispersion parameter of a gaussian fit, which `lm()` does not, so
#' `dfcom` is one smaller than mice's default; pass `dfcom` explicitly
#' to reproduce mice exactly.
#'
#' @param object A `frmtmb_multiple` from [frm_multiple()].
#' @param ... A second `frmtmb_multiple`, nested with `object` and fit
#'   to the same imputed datasets.
#' @param method Combining rule: `"D3"`, `"D1"` or `"D2"`.
#' @param use For `method = "D2"`, whether the per-imputation
#'   statistics are likelihood-ratio (`"likelihood"`) or Wald
#'   (`"wald"`) chi-squares.
#' @param constraint Instead of a second model, the names of the
#'   coefficients to test jointly against zero (`"D1"` and
#'   `method = "D2", use = "wald"` only).
#' @param dfcom Complete-data residual degrees of freedom for the
#'   `"D1"` reference distribution. Defaults to `df.residual()` of the
#'   larger model's first fit; `Inf` selects the unadjusted form.
#' @return A one-row `frmtmb_pooled_anova` table: `statistic`, `df1`,
#'   `df2`, `p`, `riv` (the relative increase in variance).
#' @references
#' Li, K. H., Raghunathan, T. E. and Rubin, D. B. (1991).
#' Large-sample significance levels from multiply imputed data using
#' moment-based statistics and an F reference distribution.
#' *JASA* 86, 1065-1073.
#'
#' Li, K. H., Meng, X.-L., Raghunathan, T. E. and Rubin, D. B. (1991).
#' Significance levels from repeated p-values with multiply-imputed
#' data. *Statistica Sinica* 1, 65-92.
#'
#' Meng, X.-L. and Rubin, D. B. (1992). Performing likelihood ratio
#' tests with multiply-imputed data sets. *Biometrika* 79, 103-111.
#' @seealso [frm_multiple()], [anova.frmtmb_fit()].
#' @examples
#' set.seed(4)
#' n <- 60
#' imps <- lapply(1:4, function(i) {
#'   x <- rnorm(n)
#'   data.frame(y = rnorm(n, 1 + 0.5 * x), x = x, z = rnorm(n))
#' })
#' m1 <- frm_multiple(bf(y ~ x + z) + gaussian(), data = imps)
#' m0 <- frm_multiple(bf(y ~ x) + gaussian(), data = imps)
#' anova(m1, m0)
#' anova(m1, m0, method = "D1")
#' @export
anova.frmtmb_multiple <- function(object, ...,
                                  method = c("D3", "D1", "D2"),
                                  use = c("likelihood", "wald"),
                                  constraint = NULL, dfcom = NULL) {
  method <- frm_match_arg(method)
  use <- frm_match_arg(use)
  others <- Filter(function(x) inherits(x, "frmtmb_multiple"), list(...))
  if (length(others) > 1L) {
    frm_stop("anova() compares two frmtmb_multiple fits at a time",
             call. = FALSE)
  }
  # Same refusal as anova.frmtmb_fit: the restricted likelihood depends
  # on the fixed-effect design, so neither its value nor the Wald
  # covariance built from it compares across models.
  for (mf in c(list(object), others)) {
    if (any(vapply(mf$fits, `[[`, TRUE, "REML"))) {
      frm_stop("Pooled model comparison requires ML fits (REML = FALSE)",
               call. = FALSE)
    }
  }
  wald <- method == "D1" || (method == "D2" && use == "wald")
  if (!length(others)) {
    if (is.null(constraint)) {
      frm_stop("anova() needs a second frmtmb_multiple fit, or a ",
               "`constraint` naming the coefficients to test", call. = FALSE)
    }
    if (!wald) {
      frm_stop("`constraint` needs a Wald rule: method = \"D1\", or ",
               "method = \"D2\" with use = \"wald\". Likelihood rules ",
               "need a fitted null model", call. = FALSE)
    }
    parts <- pooled_wald_parts(object, NULL, constraint)
  } else {
    ord <- pooled_order(object, others[[1]])
    if (wald) {
      parts <- pooled_wald_parts(ord$big, ord$small, NULL)
    } else {
      parts <- pooled_lrt_parts(ord$big, ord$small, method)
    }
  }
  m <- object$m

  if (method == "D1") {
    if (is.null(dfcom)) dfcom <- df.residual(parts$big_fit)
    st <- d1_stat(parts$Q, parts$U, m, dfcom, parts$K)
  } else if (method == "D2") {
    st <- d2_stat(parts$d, parts$k, m)
  } else {
    st <- d3_stat(parts$dbar, parts$dtilde, parts$k, m)
  }
  # a Wald rule through the factors whose inverse does not exist
  if (wald && !is.null(parts$K) && !is.finite(st$statistic)) {
    warn_cov_singular(parts$big_fit, "anova()")
  }
  p <- stats::pf(st$statistic, st$df1, st$df2, lower.tail = FALSE)
  structure(
    data.frame(statistic = st$statistic, df1 = st$df1, df2 = st$df2,
               p = p, riv = st$riv, row.names = NULL),
    class = c("frmtmb_pooled_anova", "data.frame"),
    method = method, use = if (method == "D2") use else NULL,
    m = m, models = parts$labels,
    dfcom = if (method == "D1") dfcom else NULL
  )
}

#' @export
print.frmtmb_pooled_anova <- function(x, digits = 4, ...) {
  frm_check_dots(...)
  cat("Pooled model comparison over ", attr(x, "m"), " imputations (",
      attr(x, "method"),
      if (!is.null(attr(x, "use"))) paste0(", ", attr(x, "use")),
      ")\n", sep = "")
  lab <- attr(x, "models")
  for (i in seq_along(lab)) {
    cat("  Model ", i, ": ", lab[i], "\n", sep = "")
  }
  cat("\n")
  print(as.data.frame(lapply(unclass(x), signif, digits)),
        row.names = FALSE)
  invisible(x)
}

#' Order two pooled fits by model size (logLik df of the first fit) and
#' check that they were fit to the same imputed datasets.
#'
#' @noRd
pooled_order <- function(a, b) {
  if (a$m != b$m) {
    frm_stop("anova() needs both fits pooled over the same imputations ",
             "(got m = ", a$m, " and ", b$m, ")", call. = FALSE)
  }
  # Same guard as anova.frmtmb_fit, applied imputation by imputation:
  # likelihoods computed on different data are not on a common scale.
  na <- vapply(a$fits, function(f) as.integer(f$frame[["n_obs"]]), 0L)
  nb <- vapply(b$fits, function(f) as.integer(f$frame[["n_obs"]]), 0L)
  if (!identical(na, nb)) {
    frm_stop("anova() needs both fits fit to the same imputed datasets; ",
             "imputation ", which(na != nb)[1], " has ", na[which(na != nb)[1]],
             " observations in one fit and ", nb[which(na != nb)[1]],
             " in the other", call. = FALSE)
  }
  dfa <- attr(logLik(a$fits[[1]]), "df")
  dfb <- attr(logLik(b$fits[[1]]), "df")
  if (dfa == dfb) {
    frm_stop("anova() needs nested fits of different size (both have ",
             dfa, " parameters)", call. = FALSE)
  }
  if (dfa > dfb) list(big = a, small = b) else list(big = b, small = a)
}

#' `Qhat (k x m)` and `Uhat (k x k x m)` for the coefficients under
#' test: those the larger model adds, or the named `constraint` set.
#'
#' @noRd
pooled_wald_parts <- function(big, small, constraint) {
  nm1 <- estimated_coef_names(big$fits[[1]])
  if (is.null(small)) {
    bad <- setdiff(constraint, nm1)
    if (length(bad)) {
      frm_stop("constraint names no coefficient of the fit: ",
               paste(bad, collapse = ", "), call. = FALSE)
    }
    tested <- nm1[nm1 %in% constraint]
    labels <- c(model_label(big$fits[[1]]),
                paste0("constraint: ",
                       paste(paste(tested, "= 0"), collapse = ", ")))
  } else {
    nm0 <- estimated_coef_names(small$fits[[1]])
    if (length(setdiff(nm0, nm1))) {
      frm_stop("the smaller fit is not nested in the larger one: ",
               paste(setdiff(nm0, nm1), collapse = ", "),
               " has no counterpart", call. = FALSE)
    }
    # A Wald test on the coefficients cannot express a difference in
    # the covariance parameters, so refuse rather than test the wrong
    # hypothesis. Checked before the empty-difference case, which the
    # same situation would otherwise trip first.
    nth <- function(mf) {
      pn <- names(mf$fits[[1]]$opt$par)
      length(pn[!pn %in% c("beta", "betad")])
    }
    if (nth(big) != nth(small)) {
      frm_stop("the fits differ in their covariance parameters, which a ",
               "Wald rule cannot test. Use method = \"D3\" or ",
               "method = \"D2\" with use = \"likelihood\"", call. = FALSE)
    }
    tested <- setdiff(nm1, nm0)
    if (!length(tested)) {
      frm_stop("the two fits have the same coefficients; a Wald rule has ",
               "nothing to test. Use method = \"D3\"", call. = FALSE)
    }
    labels <- c(model_label(big$fits[[1]]), model_label(small$fits[[1]]))
  }
  k <- length(tested)
  m <- big$m
  Q <- matrix(NA_real_, k, m)
  U <- array(NA_real_, c(k, k, m))
  for (j in seq_len(m)) {
    f <- big$fits[[j]]
    cf <- c(f$estimates[["beta"]],
            if (length(fx <- f$frame[["betad_fixed_idx"]])) {
              f$estimates[["betad"]][-fx]
            } else {
              f$estimates[["betad"]]
            })
    names(cf) <- estimated_coef_names(f)
    Q[, j] <- cf[tested]
    # `tested` are estimated_coef_names(), so the covariance has to be
    # the one over that vector rather than brms's block
    U[, , j] <- vcov_estimated(f)[tested, tested, drop = FALSE]
  }
  # per-imputation Wald chi-squares, for D2(use = "wald"), through the
  # factors where the fits hold them: solve() of U stops as
  # computationally singular when `tested` are the powers of a degree-7
  # raw polynomial (dev/illse2-findings.md, item 1)
  K <- pooled_wald_roots(big$fits, tested, U)
  d <- if (is.null(K)) {
    vapply(seq_len(m), function(j) {
      drop(t(Q[, j]) %*% solve(U[, , j], Q[, j]))
    }, 0)
  } else {
    vapply(seq_len(m), function(j) {
      root_inv_quad(K[[j]], Q[, j]) %||% NaN
    }, 0)
  }
  list(Q = Q, U = U, d = d, k = k, labels = labels,
       big_fit = big$fits[[1]], K = K)
}

#' Factors `K_j` of the per-imputation covariances `U[, , j]` of the
#' coefficients `tested`, `U_j = K_j K_j'`, when at least one fit holds
#' its covariance as a factor (se_reparam_cov()); NULL otherwise, and
#' then the Wald rules invert `U` as they did. A fit without a factor
#' gets the pivoted Cholesky factor of its `U_j`, whose condition the
#' fit's own SE check accepted; NULL when one has none.
#'
#' @noRd
pooled_wald_roots <- function(fits, tested, U) {
  if (!all(is.finite(U))) return(NULL)
  K <- lapply(fits, function(f) {
    Kf <- coef_vcov_root(f, extra = FALSE)
    if (is.null(Kf) || !all(tested %in% rownames(Kf))) return(NULL)
    Kf[tested, , drop = FALSE]
  })
  have <- !vapply(K, is.null, NA)
  if (!any(have)) return(NULL)
  for (j in which(!have)) {
    R <- tryCatch(suppressWarnings(chol(U[, , j], pivot = TRUE)),
                  error = function(e) NULL)
    if (is.null(R) || attr(R, "rank") < nrow(R)) return(NULL)
    Kj <- matrix(0, nrow(R), ncol(R))
    Kj[attr(R, "pivot"), ] <- t(R)
    K[[j]] <- Kj
  }
  K
}

#' Per-imputation deviance differences, at the imputation-specific
#' estimates and (for D3) at the pooled parameter vector.
#'
#' @noRd
pooled_lrt_parts <- function(big, small, method) {
  m <- big$m
  k <- attr(logLik(big$fits[[1]]), "df") -
    attr(logLik(small$fits[[1]]), "df")
  dev1 <- vapply(big$fits, function(f) 2 * f$opt$objective, 0)
  dev0 <- vapply(small$fits, function(f) 2 * f$opt$objective, 0)
  out <- list(k = k, d = pmax(dev0 - dev1, 0),
              labels = c(model_label(big$fits[[1]]),
                         model_label(small$fits[[1]])),
              big_fit = big$fits[[1]])
  if (method == "D3") {
    out$dbar <- mean(dev0 - dev1)
    out$dtilde <- mean(dev_at_pooled(small$fits) - dev_at_pooled(big$fits))
  }
  out
}

#' Deviance of every imputation's own likelihood at the across-imputation
#' mean of the optimizer's parameter vector. The stored objective is
#' reusable: `obj$fn()` reruns the inner Laplace solve at the new outer
#' values, so this costs one likelihood evaluation per imputation rather
#' than a refit. `fn()` leaves its argument in the tape's `last.par`,
#' which `sdreport()` and the conditional modes read, so restore the
#' optimum.
#'
#' @noRd
dev_at_pooled <- function(fits) {
  p1 <- fits[[1]]$opt$par
  P <- vapply(fits, function(f) {
    p <- f$opt$par
    if (!identical(names(p), names(p1))) {
      frm_stop("the imputations produced different parameter vectors; ",
               "method = \"D3\" needs one common parameterization",
               call. = FALSE)
    }
    p
  }, numeric(length(p1)))
  pbar <- rowMeans(matrix(P, nrow = length(p1)))
  vapply(fits, function(f) {
    val <- 2 * as.numeric(f$obj$fn(pbar))
    f$obj$fn(f$opt$par)
    val
  }, 0)
}

#' Denominator df of the F reference when the complete-data df is not
#' used (Li, Raghunathan and Rubin 1991 eq. 2.2). A relative increase
#' `r` of zero gives `Inf`, the complete-data chi-square limit.
#'
#' @noRd
d_denom_df <- function(r, k, m) {
  tdf <- k * (m - 1)
  if (tdf > 4) {
    4 + (tdf - 4) * (1 + (1 - 2 / tdf) / r)^2
  } else {
    tdf * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  }
}

#' D1 denominator df: the Reiter (2007) small-sample form when the
#' complete-data df is finite, which is what `mice::D1()` uses by
#' default. Every one of its terms carries a factor `a = r t / (t - 2)`,
#' so a zero ARIV collapses it to the complete-data df. Below `t = 4`
#' its `c0 = 1 / (t - 4)` turns negative and can carry the result past
#' zero, where `pf()` has no answer; fall back to the large-sample form
#' there and on any other degenerate value.
#'
#' @noRd
d1_denom_df <- function(r, k, m, dfcom) {
  tdf <- k * (m - 1)
  if (!is.finite(dfcom) || tdf < 4) return(d_denom_df(r, k, m))
  vstar <- ((dfcom + 1) / (dfcom + 3)) * dfcom
  if (r <= 0) return(vstar)
  a <- r * tdf / (tdf - 2)
  c0 <- 1 / (tdf - 4)
  c1 <- vstar - 2 * (1 + a)
  c2 <- vstar - 4 * (1 + a)
  z <- 1 / c2 +
    c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
    c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
    c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) +
    c0 * (8 * a^2 / c2^2)
  v <- 4 + 1 / z
  if (is.na(v) || v <= 0) d_denom_df(r, k, m) else v
}

#' D1, the multivariate Wald rule of Li, Raghunathan and Rubin (1991).
#' Pools the coefficients `Q` and their within-imputation covariances
#' `U`, then inflates the pooled covariance by the average relative
#' increase in variance. The result is referred to an F distribution on
#' `k` and the `d1_denom_df()` denominator degrees of freedom. `K`, the
#' factors of `U` (pooled_wald_roots()), replaces `U` where it is given.
#'
#' @noRd
d1_stat <- function(Q, U, m, dfcom, K = NULL) {
  k <- nrow(Q)
  qbar <- rowMeans(Q)
  if (!is.null(K)) {
    # the mean of the factored covariances is factored too, so neither
    # the trace nor the form inverts an ill-conditioned ubar
    Kb <- do.call(cbind, K) / sqrt(m)
    C <- (Q - qbar) / sqrt(m - 1)
    tr <- root_inv_quad(Kb, C)
    qf <- root_inv_quad(Kb, qbar)
    if (is.null(tr) || is.null(qf)) {
      return(list(statistic = NaN, df1 = k, df2 = NaN, riv = NaN))
    }
    r <- max(0, (1 + 1 / m) * sum(tr) / k)
    return(list(statistic = qf / ((1 + r) * k), df1 = k,
                df2 = d1_denom_df(r, k, m, dfcom), riv = r))
  }
  ubar <- apply(U, c(1, 2), mean)
  bvar <- stats::var(t(Q))
  # ARIV clamped at zero: the moment estimator can go negative on few
  # imputations, and a negative r makes the reference df negative too
  # (mitml's ariv = "positive").
  r <- max(0, (1 + 1 / m) * sum(diag(bvar %*% solve(ubar))) / k)
  stat <- drop(t(qbar) %*% solve((1 + r) * ubar, qbar)) / k
  df2 <- d1_denom_df(r, k, m, dfcom)
  list(statistic = stat, df1 = k, df2 = df2, riv = r)
}

#' D2, the pooled test statistic rule of Li, Meng, Raghunathan and Rubin
#' (1991). Combines the `m` per-imputation chi-squares `d`, so it needs
#' no covariance matrix. The result is referred to an F distribution on
#' `k` and `k^(-3/m) (m - 1) (1 + 1/r)^2` degrees of freedom.
#'
#' @noRd
d2_stat <- function(d, k, m) {
  r <- max(0, (1 + 1 / m) * stats::var(sqrt(d)))
  stat <- (mean(d) / k - (m + 1) / (m - 1) * r) / (1 + r)
  list(statistic = stat, df1 = k,
       df2 = k^(-3 / m) * (m - 1) * (1 + 1 / r)^2, riv = r)
}

#' D3, the pooled likelihood ratio rule of Meng and Rubin (1992). Takes
#' the average deviance difference at the imputation-specific estimates
#' (`dbar`) and at the pooled estimates (`dtilde`). The result is
#' referred to an F distribution on `k` and the `d_denom_df()`
#' denominator degrees of freedom.
#'
#' @noRd
d3_stat <- function(dbar, dtilde, k, m) {
  r <- max(0, (m + 1) / (k * (m - 1)) * (dbar - dtilde))
  list(statistic = dtilde / (k * (1 + r)), df1 = k,
       df2 = d_denom_df(r, k, m), riv = r)
}

#' @rdname hypothesis
#' @exportS3Method brms::hypothesis
#' @export
hypothesis.frmtmb_multiple <- function(x, hypothesis, class = "b",
                                       group = "",
                                       scope = c("standard", "ranef",
                                                 "coef"),
                                       alpha = 0.05, robust = FALSE,
                                       seed = NULL, ...) {
  # refused, not warned about: a pooled table built as if the argument
  # had not been given is a wrong answer with a note beside it
  frm_check_dots(...)
  hyp_refuse_draws_args(scope, robust, "hypothesis()")
  check_probability(alpha, "alpha")
  vo <- hyp_vals_only(x$fits[[1]])
  prefix <- hyp_class_prefix(class, group)
  hp <- hyp_parse_all(hypothesis,
                      names(hyp_env_vals(x$fits[[1]], vo$vals, vo$comp)),
                      class, group)
  exs <- hp$exprs
  labels <- hyp_labels(hypothesis)
  m <- x$m
  Q <- matrix(NA_real_, length(exs), m)
  U <- matrix(NA_real_, length(exs), m)
  for (j in seq_len(m)) {
    fit <- x$fits[[j]]
    pc <- hyp_par_cov(fit)
    for (i in seq_along(exs)) {
      ex <- exs[[i]]
      val <- hyp_eval(fit, ex, pc$vals, pc$comp)
      if (!is.numeric(val) || length(val) != 1L) {
        frm_stop("Hypothesis '", hypothesis[i], "' must evaluate to a ",
                 "single number in imputation ", j, " of ", m, call. = FALSE)
      }
      Q[i, j] <- val
      # NaN along a direction this imputation's Hessian lost (RB2)
      U[i, j] <- as.numeric(hyp_delta_var(
        pc, function(v) hyp_eval(fit, ex, v, pc$comp), pc$vals))
    }
  }
  se_pred_warn(rowSums(is.nan(U)) > 0, "hypotheses")
  pl <- rubin_pool(Q, U, df.residual(x$fits[[1]]))
  k_n <- length(exs)
  lwr <- upr <- stat <- p <- numeric(k_n)
  for (i in seq_len(k_n)) {
    # the reference is the Barnard-Rubin t of this row, so the
    # quantile and tail functions go in per row
    df_i <- pl$df[i]
    wr <- hyp_wald_row(pl$estimate[i], pl$se[i], hp$dir[i], alpha,
                       function(p) stats::qt(p, df_i),
                       function(q) stats::pt(q, df_i))
    stat[i] <- wr$stat
    p[i] <- wr$p
    # brms's interval: central 1 - 2 alpha for a directional row
    lo <- if (hp$dir[i] == "two.sided") alpha / 2 else alpha
    q <- stats::qt(1 - lo, df_i)
    lwr[i] <- pl$estimate[i] - q * pl$se[i]
    upr[i] <- pl$estimate[i] + q * pl$se[i]
  }
  na_col <- rep(NA_real_, k_n)
  hyp_brms_result(
    labels, pl$estimate, pl$se, lwr, upr, na_col, na_col, hp$dir,
    ifelse(hp$dir == "greater", lwr > 0, upr < 0),
    hyp_samples_frame(NULL, k_n),
    hyp_samples_frame(matrix(NA_real_, 0L, k_n), k_n), prefix, alpha,
    list(method = "wald",
         test = data.frame(Hypothesis = labels, t = stat, df = pl$df,
                           p = p)))
}

## Post-processing entry points that a pooled fit cannot serve.
##
## A frm_multiple() result is m maximum-likelihood fits plus the Rubin
## pooling of their coefficients. It has no draws, no chains, and no
## pooled predictive distribution, so the brms/posterior entry points
## that assume one of those are refused by name rather than left to
## fail inside base graphics or a posterior-package default method.

#' The one refusal for the posterior-package entry points: name the
#' function the user called, say why there is nothing to hand over, and
#' name what to call instead.
#'
#' @noRd
multiple_no_draws <- function(fn) {
  frm_stop(fn, "() needs draws, and a frm_multiple() result has none: it ",
           "is m maximum-likelihood fits pooled by Rubin's rules, with no ",
           "chains. Read the pooled tables from `x$pooled` and ",
           "`x$pooled_varcorr` or test with hypothesis(), and use ",
           "frmtmb.sample::frm_sample(x) for the draws of every imputation ",
           "combined, as brm_multiple() combines them",
           call. = FALSE)
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws
#' @export
as_draws.frmtmb_multiple <- function(x, ...) {
  multiple_no_draws("as_draws")
}

#' @exportS3Method posterior::as_draws_array
as_draws_array.frmtmb_multiple <- function(x, ...) {
  multiple_no_draws("as_draws_array")
}

#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.frmtmb_multiple <- function(x, ...) {
  multiple_no_draws("as_draws_matrix")
}

#' @exportS3Method posterior::as_draws_df
as_draws_df.frmtmb_multiple <- function(x, ...) {
  multiple_no_draws("as_draws_df")
}

# posterior's nchains() and ndraws() generics take x alone, so these
# methods do too
#' @exportS3Method posterior::nchains
nchains.frmtmb_multiple <- function(x) multiple_no_draws("nchains")

#' @exportS3Method posterior::ndraws
ndraws.frmtmb_multiple <- function(x) multiple_no_draws("ndraws")

#' @export
plot.frmtmb_multiple <- function(x, ...) {
  # brms's plot() of a brm_multiple() fit draws the trace and the
  # histogram of the combined draws, chosen with variable =, regex =
  # and N =; there are no draws here to choose from
  frm_stop("plot() has no pooled display for a frm_multiple() result: ",
           "brms's display is the posterior draws of each parameter, and ",
           "this is m maximum-likelihood fits pooled by Rubin's rules. ",
           "summary() and fixef() give the pooled coefficients, ",
           "conditional_effects() the pooled effect curves, and ",
           "`plot(x$fits[[1]])` one imputation's residual diagnostics",
           call. = FALSE)
}

## Pooled post-processing: the brms_missings vignette's calls on a
## brm_multiple() fit, with Rubin's rules where brms combines the draws
## of the imputations.

#' Pooled post-processing of a frm_multiple() fit
#'
#' brms's `brm_multiple()` combines the posterior draws of the m
#' imputations, so `fixef()`, `summary()`, `conditional_effects()` and
#' the other methods of its result summarize the combined draws. A
#' [frm_multiple()] result is m maximum-likelihood fits, and these
#' methods pool them by Rubin's rules with the Barnard-Rubin degrees of
#' freedom instead:
#'
#' * `fixef()` pools `fixef()` of each imputation, so it has brms's rows,
#'   names and order, ordinal thresholds included. `Estimate` is the
#'   mean over the imputations, `Est.Error` the square root of the total
#'   variance (within plus `1 + 1/m` times between), and the quantile
#'   columns a t interval on the pooled degrees of freedom.
#' * `summary()` prints that table, the distributional parameters
#'   without a formula pooled on their link scale and transformed back,
#'   and the pooled variance components, with the degrees of freedom
#'   and the fraction of missing information (`fmi`) of each row where
#'   brms prints R-hat and effective sample sizes.
#' * `conditional_effects()` computes the curves of every imputation on
#'   the grid of the first (brms's grid, which `brm_multiple()` reads
#'   from the first data set), pools each grid point on the scale its
#'   Wald band is symmetric on (the link, or the logit of an ordinal
#'   category probability), and transforms back. A category probability
#'   is the exception for its estimate: that is the mean of the
#'   imputations' probabilities, so the categories sum to 1, and the
#'   band is moved out to it where needed. Only the Wald band pools this
#'   way: `method = "predict"` and the bootstrap and profile bands are
#'   refused by name.
#' * `vcov()` is the total covariance of the `fixef()` rows: the mean of
#'   each imputation's `vcov()` plus `1 + 1/m` times the covariance of
#'   the estimates between imputations, which is `mice::pool()`'s.
#' * `fitted()` pools `fitted()` of each imputation row by row. It pools
#'   on the link scale of the predictor, as `conditional_effects()`
#'   does: each estimate goes through the link and its standard error
#'   through the delta method. The estimate and a t interval on the
#'   pooled degrees of freedom come back through the inverse link, so
#'   the interval stays inside the response's range, and `Est.Error` is
#'   the pooled standard error carried back by the delta method. A
#'   category probability takes its interval from the logit scale, but
#'   its estimate and `Est.Error` from Rubin's rules on the probability
#'   scale, so the categories sum to 1; an end on the wrong side of the
#'   estimate is moved to it. These are the values of the pooled
#'   `conditional_effects()`. `scale = "linear"` pools
#'   the linear predictor itself. In sample, each imputation predicts
#'   its own rows, so an imputed covariate enters at the value of each
#'   imputation; brms answers from the first imputed data set alone and
#'   warns that it does. With `newdata`, every imputation predicts the
#'   same rows.
#' * `coef()` and `ranef()` are the mean over the imputations of the
#'   coefficients and conditional modes of each level, which is Rubin's
#'   point estimate. `ranef(condVar = TRUE)` pools the conditional
#'   standard deviations as Rubin's total variance: the mean
#'   conditional variance plus `1 + 1/m` times the variance of the modes.
#' * `VarCorr()` pools each standard deviation on the log scale and each
#'   correlation on the Fisher z scale, the scales `frm_multiple()`
#'   pools its variance components on, and transforms back, with
#'   `Est.Error` by the delta method. A covariance pools on its own
#'   scale. A component at its boundary (a standard deviation of zero in
#'   one imputation) has no log, and its row is `NA`.
#' * `emmeans::emmeans()` pools the estimated marginal means of each
#'   imputation as emmeans pools a `mice::mira`: at the covariate means
#'   over the imputations, with the Barnard-Rubin degrees of freedom on
#'   the asymptotic reference that emmeans uses for one fit.
#' * `nobs()` is the number of observations, which every imputation
#'   shares, and `formula()` and `family()` are those of the model every
#'   imputation is fitted with.
#'
#' `predict()`, `residuals()`, `pp_check()`, `logLik()` (and so `AIC()`
#' and `BIC()`), `loo()`, `waic()`, `update()` and `model.frame()` have
#' no Rubin's-rules meaning. Each is refused by name, with what to call instead:
#' `anova()` pools a model comparison, and `x$fits[[j]]` is the fit of
#' imputation `j`.
#'
#' brms's draws arguments have no meaning on pooled fits and are
#' refused by name, as on one fit.
#'
#' @param object,x A `frmtmb_multiple` from [frm_multiple()].
#' @param summary,robust brms's arguments, refused when not the
#'   default: there are no draws.
#' @param probs Probabilities of the quantile columns of `fixef()`,
#'   `fitted()` and `VarCorr()`.
#' @param pars Row names to keep, in the order given, as in brms. For
#'   `ranef()`, the coefficients to keep, as in [ranef()].
#' @param priors If `TRUE`, `summary()` also prints the prior of the
#'   first imputation's fit.
#' @param prob The coverage of the summary's and the curves' intervals.
#' @param mc_se brms's argument, refused when `TRUE`.
#' @param waic Old brms printed WAIC in `summary()` when it was `TRUE`.
#'   Maximum-likelihood fits have no WAIC, so `TRUE` is refused by name
#'   (brms 2.23.0 ignores the argument).
#' @param effects,... For `conditional_effects()`, passed to
#'   [conditional_effects()] on each imputation's fit. Otherwise refused
#'   by name.
#' @return `fixef()`: brms's coefficient matrix. `summary()`: a
#'   `summary.frmtmb_multiple` object, printed in brms's layout.
#'   `conditional_effects()`: a conditional-effects object, which
#'   `plot()` draws. `vcov()`, `coef()`, `ranef()`, `VarCorr()` and
#'   `fitted()`: the shape the same method gives on one fit. `nobs()`: a
#'   number.
#' @seealso [frm_multiple()], [hypothesis()] for pooled functions of
#'   the parameters, [anova.frmtmb_multiple()].
#' @examples
#' set.seed(8)
#' n <- 80
#' x <- rnorm(n)
#' y <- rnorm(n, 1 + 0.5 * x, 1)
#' x[sample(n, 15)] <- NA
#' imps <- lapply(1:3, function(i) {
#'   xi <- x
#'   xi[is.na(xi)] <- sample(x[!is.na(x)], sum(is.na(xi)), TRUE)
#'   data.frame(y = y, x = xi)
#' })
#' fm <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps)
#' fixef(fm)
#' summary(fm)
#' # the total covariance, whose diagonal is the squared Est.Error
#' vcov(fm)
#' fitted(fm, newdata = data.frame(x = c(-1, 0, 1)))
#' ce <- conditional_effects(fm, "x")
#' head(ce$x[, c("x", "estimate__", "lower__", "upper__")])
#' @name frm_multiple-methods
NULL

#' @rdname frm_multiple-methods
#' @export
fixef.frmtmb_multiple <- function(object, summary = TRUE, robust = FALSE,
                                  probs = c(0.025, 0.975), pars = NULL,
                                  ...) {
  frm_check_dots(...)
  fit_refuse_draws_args("fixef()", summary = summary, robust = robust)
  pl <- multiple_pool_fixef(object)
  qs <- lapply(probs %||% numeric(0), function(p) {
    pl$estimate + stats::qt(p, pl$df) * pl$se
  })
  out <- cbind(Estimate = pl$estimate, `Est.Error` = pl$se)
  if (length(qs)) out <- cbind(out, do.call(cbind, qs))
  colnames(out) <- c("Estimate", "Est.Error", brms_prob_cols(probs))
  rownames(out) <- pl$names
  brms_pars_filter(out, pars, "fixef()")
}

#' Rubin's rules over `fixef()` of each imputation: brms's rows, names
#' and order, with the pooled degrees of freedom and fmi.
#'
#' @noRd
multiple_pool_fixef <- function(object) {
  fx <- lapply(object$fits, function(f) fixef(f))
  rn <- rownames(fx[[1L]])
  for (j in seq_along(fx)[-1L]) {
    if (!identical(rownames(fx[[j]]), rn)) {
      frm_stop("the imputations' fits have different coefficients (",
               "imputation ", j, " differs from the first), so there is ",
               "no pooled table: a factor level, say, is missing from ",
               "one of the imputed data sets", call. = FALSE)
    }
  }
  Q <- matrix(vapply(fx, function(m) m[, "Estimate"], numeric(length(rn))),
              nrow = length(rn))
  U <- matrix(vapply(fx, function(m) m[, "Est.Error"]^2,
                     numeric(length(rn))), nrow = length(rn))
  pl <- rubin_pool(Q, U, df.residual(object$fits[[1L]]))
  pl$names <- rn
  pl
}

#' A pooled summary block: brms's four columns with a t interval on the
#' pooled degrees of freedom, and `df` and `fmi` where brms prints its
#' convergence columns.
#'
#' @noRd
multiple_block <- function(est, se, df, fmi, prob, names,
                           linkinv = identity, mu_eta = NULL) {
  q <- stats::qt(1 - (1 - prob) / 2, df)
  lo <- linkinv(est - q * se)
  hi <- linkinv(est + q * se)
  e <- if (is.null(mu_eta)) se else abs(mu_eta(est)) * se
  pct <- format(prob * 100, drop0trailing = TRUE)
  out <- data.frame(linkinv(est), e, pmin(lo, hi), pmax(lo, hi), df, fmi,
                    row.names = names)
  names(out) <- c("Estimate", "Est.Error", paste0("l-", pct, "% CI"),
                  paste0("u-", pct, "% CI"), "df", "fmi")
  out
}

#' @rdname frm_multiple-methods
#' @export
summary.frmtmb_multiple <- function(object, priors = FALSE, prob = 0.95,
                                    robust = FALSE, mc_se = FALSE, ...,
                                    waic = FALSE) {
  frm_check_dots(..., .unsupported = brms_summary_args)
  summary_refuse_waic(waic)
  check_flag(priors, "priors")
  check_probability(prob, "prob")
  fit_refuse_draws_args("summary()", robust = robust)
  check_flag(mc_se, "mc_se")
  if (mc_se) {
    frm_stop("summary() of a frm_multiple() result cannot honor ",
             "mc_se = TRUE: ", brms_summary_args[["mc_se"]], call. = FALSE)
  }
  f1 <- object$fits[[1L]]
  s1 <- summary(f1, priors = priors)
  pl <- multiple_pool_fixef(object)
  fixed <- multiple_block(pl$estimate, pl$se, pl$df, pl$fmi, prob,
                          pl$names)
  # the distributional parameters nobody wrote a formula for: pooled on
  # their link scale in `$pooled`, reported on their own, as brms does
  tab <- brms_coef_table(f1)
  nat <- which(tab$natural)
  spec <- NULL
  if (length(nat)) {
    pooled <- object$pooled
    spec <- do.call(rbind, lapply(nat, function(i) {
      lp <- f1$frame[["linpreds"]][[linpred_key(tab$resp[i], tab$dpar[i])]]
      multiple_block(pooled$estimate[i], pooled$se[i], pooled$df[i],
                     pooled$fmi[i], prob, tab$brms[i],
                     linkinv = lp[["link"]]$linkinv,
                     mu_eta = lp[["link"]]$mu_eta)
    }))
  }
  structure(list(m = object$m, family = s1$family, links = s1$links,
                 formula = s1$formula, formulas = s1$formulas,
                 data_name = "the imputed data sets", nobs = s1$nobs,
                 ngrps = s1$ngrps, fixed = fixed, spec_pars = spec,
                 varcorr = object$pooled_varcorr, level = object$level,
                 prob = prob, priors = priors, prior = s1$prior),
            class = "summary.frmtmb_multiple")
}

#' @export
print.summary.frmtmb_multiple <- function(x, digits = 2, ...) {
  frm_check_dots(...)
  mv_forms <- x$formulas %||% list()
  if (length(mv_forms) > 1L) {
    cat(" Family: MV(", paste(vapply(x$family, `[[`, "", "family"),
                              collapse = ", "), ") \n", sep = "")
  } else {
    cat(" Family:", x$family[["family"]], "\n")
  }
  cat_family_links(x$links %||% family_link_str(x$family))
  if (length(mv_forms) > 1L) {
    cat("Formula: ", paste(vapply(mv_forms, deparse1, ""),
                           collapse = " \n         "), " \n", sep = "")
  } else {
    cat("Formula:", deparse1(x$formula), "\n")
  }
  cat("   Data:", x$data_name,
      paste0("(Number of observations: ", x$nobs, ")"), "\n")
  cat(" Method: ML on each of ", x$m, " imputations, pooled by Rubin's ",
      "rules (Barnard-Rubin df)\n", sep = "")
  if (NROW(x$varcorr)) {
    cat("\nMultilevel Hyperparameters (pooled on the log and Fisher-z ",
        "scales, ", format(x$level * 100), "% interval):\n", sep = "")
    vc <- x$varcorr
    vc[-(1:3)] <- lapply(vc[-(1:3)], signif, digits + 2L)
    print(vc, row.names = FALSE)
  }
  cat("\nRegression Coefficients:\n")
  print(round(x$fixed, digits))
  if (NROW(x$spec_pars)) {
    cat("\nFurther Distributional Parameters:\n")
    print(round(x$spec_pars, digits))
  }
  if (isTRUE(x$priors)) {
    cat("\nPriors (of the first imputation's fit):\n")
    if (is.null(x$prior)) {
      cat("No priors were set (plain maximum likelihood).\n")
    } else {
      print(x$prior)
    }
  }
  invisible(x)
}

#' @rdname frm_multiple-methods
#' @exportS3Method brms::conditional_effects
#' @export
conditional_effects.frmtmb_multiple <- function(x, effects = NULL, ...,
                                                prob = 0.95) {
  dots <- list(...)
  check_probability(prob, "prob")
  if (!is.null(dots$method) &&
        !identical(ce_method(dots$method), "epred")) {
    frm_stop("conditional_effects() of a frm_multiple() result cannot ",
             "honor method = \"predict\": a predictive interval is a ",
             "quantile of simulated responses, and the simulations of the ",
             "imputations would have to be pooled as draws, which is not ",
             "implemented. Use method = \"epred\" (the default), or ",
             "conditional_effects(x$fits[[1]], method = \"predict\") for ",
             "one imputation", call. = FALSE)
  }
  if (!is.null(dots$band) && !identical(dots$band, "wald")) {
    frm_stop("conditional_effects() of a frm_multiple() result cannot ",
             "honor band = \"", dots$band, "\": Rubin's rules pool an ",
             "estimate and its standard error, which the Wald band has ",
             "and a ", dots$band, " band does not. Use band = \"wald\" ",
             "(the default)", call. = FALSE)
  }
  for (a in c("spaghetti", "surface")) {
    if (isTRUE(dots[[a]])) {
      frm_stop("conditional_effects() of a frm_multiple() result cannot ",
               "honor ", a, " = TRUE: it pools curves point by point, ",
               "and ", if (a == "spaghetti") "a spaghetti display draws " else
                 "a surface is a grid of ", "one ",
               if (a == "spaghetti") "curve per draw" else
                 "fit's points, which pool the same way only on line plots",
               ". Call it on one imputation's fit for that",
               call. = FALSE)
    }
  }
  ce1 <- conditional_effects(x$fits[[1L]], effects = effects, prob = prob,
                             ...)
  ps <- attr(ce1, "pool_scale")
  if (is.null(ps)) {
    frm_stop("conditional_effects() of a frm_multiple() result pools a ",
             "Wald band on the scale it is symmetric on, and this display ",
             "has none: its band scale is chosen grid point by grid point ",
             "(a reported dpar such as a mixing weight). Call it on one ",
             "imputation's fit, conditional_effects(x$fits[[1]])",
             call. = FALSE)
  }
  link <- ps$link
  q <- function(df) stats::qt(1 - (1 - prob) / 2, df)
  internal <- c("cond__", "effect1__", "effect2__", "estimate__", "se__",
                "lower__", "upper__", "cats__")
  resp_vars <- names(x$fits[[1L]]$spec$responses)
  out <- ce1
  for (e in names(ce1)) {
    d1 <- ce1[[e]]
    # the per-category display keys its effect "x:cats__", which names
    # the categories and not a variable of the call
    ev <- setdiff(attr(d1, "effects"), "cats__")
    e_call <- paste(ev, collapse = ":")
    covs <- setdiff(names(d1), c(internal, ev, resp_vars))
    # the first imputation's grid, given to every other fit exactly:
    # the effect values as int_conditions and the other covariates as
    # condition sets
    ic <- lapply(ev, function(v) {
      u <- unique(d1[[v]])
      if (is.factor(u)) as.character(u) else u
    })
    names(ic) <- ev
    conds <- unique(d1[covs])
    key <- function(d) {
      do.call(paste, c(lapply(d[c(ev, covs, intersect("cats__", names(d)))],
                              as.character), sep = "\r"))
    }
    k1 <- key(d1)
    tr <- function(d) {
      t <- link$linkfun(d$estimate__)
      s <- if (identical(ps$se_on, "link")) d$se__ else
        d$se__ / abs(link$mu_eta(t))
      # the response-scale estimate and standard error too, which a
      # category's estimate pools on
      r <- if (identical(ps$se_on, "link")) {
        d$se__ * abs(link$mu_eta(t))
      } else {
        d$se__
      }
      list(t = t, s = s, e = d$estimate__, r = r)
    }
    T1 <- tr(d1)
    Q <- matrix(NA_real_, nrow(d1), x$m)
    U <- Q
    R <- Q
    V <- Q
    Q[, 1L] <- T1$t
    U[, 1L] <- T1$s^2
    R[, 1L] <- T1$e
    V[, 1L] <- T1$r^2
    for (j in seq_len(x$m)[-1L]) {
      cej <- conditional_effects(x$fits[[j]], effects = e_call,
                                 prob = prob,
                                 int_conditions = ic,
                                 conditions = if (length(covs)) {
                                   conds
                                 } else {
                                   list()
                                 },
                                 ...)
      dj <- cej[[1L]]
      pos <- match(k1, key(dj))
      if (anyNA(pos)) {
        frm_stop("conditional_effects() could not evaluate imputation ",
                 j, " on the first imputation's grid for '", e, "', so ",
                 "the curves do not pool point by point. Give the grid ",
                 "yourself with int_conditions = and conditions =",
                 call. = FALSE)
      }
      Tj <- tr(dj[pos, , drop = FALSE])
      Q[, j] <- Tj$t
      U[, j] <- Tj$s^2
      R[, j] <- Tj$e
      V[, j] <- Tj$r^2
    }
    pl <- rubin_pool(Q, U, df.residual(x$fits[[1L]]))
    lo <- link$linkinv(pl$estimate - q(pl$df) * pl$se)
    hi <- link$linkinv(pl$estimate + q(pl$df) * pl$se)
    d1$estimate__ <- link$linkinv(pl$estimate)
    d1$se__ <- if (identical(ps$se_on, "link")) pl$se else
      pl$se * abs(link$mu_eta(pl$estimate))
    d1$lower__ <- pmin(lo, hi)
    d1$upper__ <- pmax(lo, hi)
    if ("cats__" %in% names(d1)) {
      # category probabilities pool their estimate on the response
      # scale, so the categories sum to 1 as each imputation's do, and
      # only the band's ends on the link scale, as fitted() of a pooled
      # fit does (multiple_pool_summaries())
      plr <- rubin_pool(R, V, df.residual(x$fits[[1L]]))
      d1$estimate__ <- plr$estimate
      d1$se__ <- if (identical(ps$se_on, "link")) {
        plr$se / abs(link$mu_eta(link$linkfun(plr$estimate)))
      } else {
        plr$se
      }
      d1$lower__ <- pmin(d1$lower__, d1$estimate__)
      d1$upper__ <- pmax(d1$upper__, d1$estimate__)
    }
    out[[e]] <- d1
  }
  attr(out, "pool_scale") <- NULL
  out
}

## Pooled point summaries beyond the coefficients: vcov(), coef(),
## ranef(), VarCorr() and fitted() pool each imputation's own answer by
## Rubin's rules, and the methods without a Rubin's-rules meaning are
## refused by name.

#' The refusal of an answer that differs in layout between imputations,
#' which pooling entry by entry cannot read.
#'
#' @noRd
multiple_shape_stop <- function(what, j) {
  frm_stop(what, " of a frm_multiple() result pools the imputations ",
           "entry by entry, and imputation ", j, " gives a differently ",
           "shaped answer from the first, so there is nothing to pool: a ",
           "factor or group level, say, is missing from one of the imputed ",
           "data sets. Call ", what, " on each fit, x$fits[[j]]",
           call. = FALSE)
}

#' Whether two answers have the same layout throughout: the same nesting,
#' names, dimensions and dimnames.
#'
#' @noRd
multiple_same_layout <- function(a, b) {
  if (is.list(a) && !is.data.frame(a)) {
    a <- unclass(a)
    b <- unclass(b)
    return(is.list(b) && !is.data.frame(b) && length(a) == length(b) &&
             identical(names(a), names(b)) &&
             all(vapply(seq_along(a), function(k) {
               multiple_same_layout(a[[k]], b[[k]])
             }, NA)))
  }
  identical(class(a), class(b)) && identical(dim(a), dim(b)) &&
    identical(dimnames(a), dimnames(b)) && identical(names(a), names(b)) &&
    length(a) == length(b)
}

#' Each imputation's answer to `f`, refused when the layouts differ.
#'
#' @noRd
multiple_answers <- function(object, what, f) {
  xs <- lapply(object$fits, f)
  for (j in seq_along(xs)[-1L]) {
    if (!multiple_same_layout(xs[[1L]], xs[[j]])) {
      multiple_shape_stop(what, j)
    }
  }
  xs
}

#' Rubin's rules over brms-shaped summaries. Each element of `xs` is one
#' imputation's `n x C` matrix or `n x C x K` array whose first two
#' columns are `Estimate` and `Est.Error`. The pooled object has the
#' same shape.
#'
#' `links`, one link object (or `NULL` for none) per layer, is the scale
#' each layer pools on, as `conditional_effects()` of a pooled fit pools
#' its curves: each imputation's estimate goes through the link, its
#' standard error through the derivative (the delta method), Rubin's
#' rules pool the two there, and the estimate and the t interval on the
#' Barnard-Rubin degrees of freedom come back through the inverse link,
#' which keeps them inside the range the link maps. `Est.Error` is the
#' pooled standard error carried back by the delta method. A cell whose
#' estimate has no value on the link scale (a probability of exactly 0)
#' pools on the response scale instead.
#'
#' A category probability (a layer of `links` with the attribute
#' `simplex = TRUE`) pools its estimate and `Est.Error` on the response
#' scale, where the categories of each imputation sum to 1 and so do
#' their means; only the interval ends come from the logit scale. The
#' two scales put the centre in slightly different places (the mean of
#' the probabilities against the inverse logit of the mean logit), so
#' an end that would fall on the wrong side of the estimate is moved to
#' it. Pooled on the logit alone, the categories of a cumulative fit
#' over 3 imputations summed to 0.99669 to 0.99983
#' (dev/postfix-rev2-rubin.R).
#'
#' @noRd
multiple_pool_summaries <- function(xs, probs, dfcom, links = NULL) {
  two <- length(dim(xs[[1L]])) == 2L
  A <- lapply(xs, function(x) if (two) array(x, c(dim(x), 1L)) else x)
  n <- dim(A[[1L]])[1L]
  K <- dim(A[[1L]])[3L]
  m <- length(A)
  Q <- matrix(vapply(A, function(a) as.vector(a[, 1L, ]), numeric(n * K)),
              ncol = m)
  S <- matrix(vapply(A, function(a) as.vector(a[, 2L, ]),
                     numeric(n * K)), ncol = m)
  # the layer of each cell, in the column-major order of as.vector()
  lay <- rep(seq_len(K), each = n)
  lk <- rep(list(NULL), K)
  if (!is.null(links)) lk[seq_along(links)] <- links
  Qt <- Q
  St <- S
  on_link <- logical(n * K)
  for (k in seq_len(K)) {
    l <- lk[[k]]
    if (is.null(l)) next
    rows <- which(lay == k)
    # a link's functions need not keep the matrix shape (identity's
    # mu_eta is a plain vector of ones), so it is put back
    t <- matrix(l$linkfun(Q[rows, , drop = FALSE]), length(rows), m)
    d <- matrix(abs(l$mu_eta(t)), length(rows), m)
    ok <- rowSums(!is.finite(t) | !is.finite(d) | d <= 0) == 0L
    rows <- rows[ok]
    Qt[rows, ] <- t[ok, , drop = FALSE]
    St[rows, ] <- S[rows, , drop = FALSE] / d[ok, , drop = FALSE]
    on_link[rows] <- TRUE
  }
  pl <- rubin_pool(Qt, St^2, dfcom)
  simplex <- vapply(lk, function(l) isTRUE(attr(l, "simplex")), NA)
  plr <- if (any(simplex)) rubin_pool(Q, S^2, dfcom)
  cn <- c("Estimate", "Est.Error", brms_prob_cols(probs))
  est <- pl$estimate
  se <- pl$se
  qs <- lapply(probs, function(p) pl$estimate + stats::qt(p, pl$df) * pl$se)
  for (k in seq_len(K)) {
    l <- lk[[k]]
    rows <- which(lay == k & on_link)
    if (is.null(l) || !length(rows)) next
    est[rows] <- l$linkinv(pl$estimate[rows])
    se[rows] <- pl$se[rows] * abs(l$mu_eta(pl$estimate[rows]))
    ends <- lapply(qs, function(q) l$linkinv(q[rows]))
    # a decreasing link (1/mu) reverses the two ends
    if (length(ends) == 2L) {
      lo <- pmin(ends[[1L]], ends[[2L]])
      hi <- pmax(ends[[1L]], ends[[2L]])
      ends <- if (probs[1L] <= probs[2L]) list(lo, hi) else list(hi, lo)
    }
    if (simplex[k]) {
      # the estimate a category's interval is built around is the
      # response-scale mean; an end on its wrong side is moved to it
      est[rows] <- plr$estimate[rows]
      se[rows] <- plr$se[rows]
      for (i in seq_along(ends)) {
        ends[[i]] <- if (probs[i] < 0.5) pmin(ends[[i]], est[rows]) else
          pmax(ends[[i]], est[rows])
      }
    }
    for (i in seq_along(qs)) qs[[i]][rows] <- ends[[i]]
  }
  out <- array(NA_real_, c(n, length(cn), K))
  out[, 1L, ] <- est
  out[, 2L, ] <- se
  for (i in seq_along(probs)) out[, 2L + i, ] <- qs[[i]]
  dn <- dimnames(xs[[1L]])
  if (two) {
    return(matrix(out, n, length(cn), dimnames = list(dn[[1L]], cn)))
  }
  dimnames(out) <- list(dn[[1L]], cn, dn[[3L]])
  out
}

#' The mean over the imputations of every number in a nested list of
#' numeric vectors, matrices or data frames, which is Rubin's point
#' estimate of each. The layout is the first answer's, attributes
#' included.
#'
#' @noRd
multiple_mean_struct <- function(xs) {
  x1 <- xs[[1L]]
  if (is.list(x1)) {
    out <- x1
    for (k in seq_along(x1)) {
      out[[k]] <- multiple_mean_struct(lapply(xs, function(x) {
        unclass(x)[[k]]
      }))
    }
    return(out)
  }
  out <- Reduce(`+`, lapply(xs, as.numeric)) / length(xs)
  attributes(out) <- attributes(x1)
  out
}

#' @rdname frm_multiple-methods
#' @param correlation If `TRUE`, `vcov()` returns the pooled correlation
#'   matrix.
#' @export
vcov.frmtmb_multiple <- function(object, correlation = FALSE, pars = NULL,
                                 ...) {
  frm_check_dots(...)
  check_flag(correlation, "correlation")
  vs <- multiple_answers(object, "vcov()", function(f) vcov(f))
  rn <- rownames(vs[[1L]])
  Q <- matrix(vapply(object$fits, function(f) {
    fixef(f)[rn, "Estimate"]
  }, numeric(length(rn))), nrow = length(rn))
  m <- object$m
  # Rubin's total covariance: the mean within-imputation covariance plus
  # (1 + 1/m) times the covariance of the estimates between imputations
  V <- Reduce(`+`, vs) / m + (1 + 1 / m) * stats::cov(t(Q))
  dimnames(V) <- list(rn, rn)
  if (correlation) V <- stats::cov2cor(V)
  brms_pars_filter(t(brms_pars_filter(t(V), pars, "vcov()")), pars,
                   "vcov()")
}

#' @rdname frm_multiple-methods
#' @export
coef.frmtmb_multiple <- function(object, summary = TRUE, robust = FALSE,
                                 probs = c(0.025, 0.975), ...) {
  frm_check_dots(..., .unsupported = brms_draws_summary_args["pars"])
  fit_refuse_draws_args("coef()", summary = summary, robust = robust,
                        probs = probs)
  multiple_mean_struct(multiple_answers(object, "coef()",
                                        function(f) coef(f)))
}

#' @rdname frm_multiple-methods
#' @param groups For `ranef()`, the grouping factors to keep, as in
#'   [ranef()].
#' @param condVar For `ranef()`: if `TRUE`, each matrix carries the
#'   pooled conditional standard deviations as its `"condSD"` attribute.
#' @export
ranef.frmtmb_multiple <- function(object, summary = TRUE, robust = FALSE,
                                  probs = c(0.025, 0.975), pars = NULL,
                                  groups = NULL, ..., condVar = FALSE) {
  frm_check_dots(...)
  fit_refuse_draws_args("ranef()", summary = summary, robust = robust,
                        probs = probs)
  check_flag(condVar, "condVar")
  xs <- multiple_answers(object, "ranef()", function(f) {
    ranef(f, pars = pars, groups = groups, condVar = condVar)
  })
  out <- unclass(xs[[1L]])
  m <- object$m
  for (k in seq_along(out)) {
    M <- lapply(xs, function(x) unclass(x)[[k]])
    E <- matrix(vapply(M, as.numeric, numeric(length(M[[1L]]))), ncol = m)
    v <- out[[k]]
    v[] <- rowMeans(E)
    if (condVar) {
      # Rubin's total variance of a predicted effect: the mean
      # conditional variance plus (1 + 1/m) times the variance of the
      # modes between imputations
      W <- matrix(vapply(M, function(x) as.numeric(attr(x, "condSD"))^2,
                         numeric(length(M[[1L]]))), ncol = m)
      S <- attr(out[[k]], "condSD")
      S[] <- sqrt(rowMeans(W) + (1 + 1 / m) * apply(E, 1, stats::var))
      attr(v, "condSD") <- S
    }
    out[[k]] <- v
  }
  structure(out, class = class(xs[[1L]]))
}

#' One block of `VarCorr()` entries, pooled by Rubin's rules on the
#' scale `frm_multiple()` pools variance components on: the log of a
#' standard deviation (`type = "sd"`) and the Fisher z of a correlation
#' (`"cor"`); a covariance (`"raw"`) on its own scale, since it has no
#' sign. `est` and `se` are quantity-by-imputation matrices of the
#' natural-scale estimates and standard errors, which the derivative of
#' the transform carries to the pooling scale. A standard deviation at
#' zero, or a correlation at one, has no value there and pools to `NA`.
#'
#' @noRd
multiple_pool_vc <- function(est, se, type, probs, dfcom) {
  tr <- switch(type, sd = log, cor = atanh, raw = identity)
  inv <- switch(type, sd = exp, cor = tanh, raw = identity)
  dtr <- switch(type,
                sd = function(v) 1 / v,
                cor = function(v) 1 / (1 - v^2),
                raw = function(v) 1 + 0 * v)
  Q <- tr(est)
  U <- (se * dtr(est))^2
  bad <- rowSums(!is.finite(Q) | !is.finite(U)) > 0L
  Q[bad, ] <- 0
  U[bad, ] <- 1
  pl <- rubin_pool(Q, U, dfcom)
  e <- inv(pl$estimate)
  out <- cbind(e, pl$se / abs(dtr(e)))
  for (p in probs) {
    out <- cbind(out, inv(pl$estimate + stats::qt(p, pl$df) * pl$se))
  }
  out[bad, ] <- NA_real_
  out
}

#' @rdname frm_multiple-methods
#' @param sigma Ignored, as in brms; nlme's generic carries it.
#' @export
VarCorr.frmtmb_multiple <- function(x, sigma = 1, summary = TRUE,
                                    robust = FALSE,
                                    probs = c(0.025, 0.975), ...) {
  frm_check_dots(...)
  fit_refuse_draws_args("VarCorr()", summary = summary, robust = robust)
  xs <- multiple_answers(x, "VarCorr()", function(f) {
    VarCorr(f, probs = probs)
  })
  dfcom <- df.residual(x$fits[[1L]])
  m <- x$m
  out <- xs[[1L]]
  for (g in seq_along(out)) {
    for (el in names(out[[g]])) {
      A <- lapply(xs, function(v) v[[g]][[el]])
      if (el == "sd") {
        col <- function(s) {
          matrix(vapply(A, function(a) a[, s], numeric(nrow(A[[1L]]))),
                 ncol = m)
        }
        out[[g]][[el]][] <- multiple_pool_vc(col("Estimate"),
                                             col("Est.Error"), "sd", probs,
                                             dfcom)
        next
      }
      # cor and cov are coefficient x statistic x coefficient arrays;
      # the diagonal of cor is 1 with error 0 and stays so
      K <- dim(A[[1L]])[1L]
      arr <- A[[1L]]
      for (i in seq_len(K)) {
        for (j in seq_len(K)) {
          if (el == "cor" && i == j) next
          cell <- function(s) {
            matrix(vapply(A, function(a) a[i, s, j], 0), nrow = 1L)
          }
          arr[i, , j] <- multiple_pool_vc(cell("Estimate"),
                                          cell("Est.Error"),
                                          if (el == "cor") "cor" else "raw",
                                          probs, dfcom)
        }
      }
      out[[g]][[el]] <- arr
    }
  }
  out
}

#' @rdname frm_multiple-methods
#' @param newdata,re_formula,scale,resp,dpar,nlpar For `fitted()`: as in
#'   [fitted.frmtmb_fit()], given to every imputation's fit. `resp` also
#'   names the response whose family `family()` returns.
#' @param allow_new_levels,sample_new_levels For `fitted()`: as in
#'   [fitted.frmtmb_fit()], given to every imputation's fit.
#' @param ndraws,draw_ids,sort For `fitted()`: brms's draws arguments, as
#'   on one fit.
#' @export
fitted.frmtmb_multiple <- function(object, newdata = NULL, re_formula = NULL,
                                   scale = c("response", "linear"),
                                   resp = NULL, dpar = NULL, nlpar = NULL,
                                   ndraws = NULL, draw_ids = NULL,
                                   sort = FALSE, summary = TRUE,
                                   robust = FALSE, probs = c(0.025, 0.975),
                                   ..., allow_new_levels = FALSE,
                                   sample_new_levels = NULL) {
  frm_check_dots(..., .unsupported = fitted_no_draws)
  scale <- frm_match_arg(scale)
  brms_prob_cols(probs)
  # once for the pooled call, so `ndraws` is said once and not per fit
  fitted_refuse_draws_args("fitted()", ndraws, draw_ids, sort, summary,
                           robust)
  xs <- multiple_answers(object, "fitted()", function(f) {
    fitted(f, newdata = newdata, re_formula = re_formula, scale = scale,
           resp = resp, dpar = dpar, nlpar = nlpar,
           allow_new_levels = allow_new_levels,
           sample_new_levels = sample_new_levels)
  })
  links <- if (identical(scale, "response")) {
    multiple_fitted_links(object$fits[[1L]], xs[[1L]], resp,
                          dpar %||% nlpar)
  }
  multiple_pool_summaries(xs, probs, df.residual(object$fits[[1L]]), links)
}

#' The scale each layer of a pooled `fitted()` pools on: the link of the
#' predictor the layer reports (the response's `mu`, or the `dpar`
#' asked for), and the logit for a category probability, the scales
#' `conditional_effects()` of a pooled fit uses. A dpar the family
#' reports on a scale of its own (a mixture's weight) pools as it is
#' reported.
#'
#' @noRd
multiple_fitted_links <- function(fit, x1, resp, dpar) {
  rn <- names(fit$spec$responses)
  one <- function(r) {
    rspec <- fit$spec$responses[[r]]
    dp <- dpar %||% (if ("mu" %in% names(rspec$dpars)) "mu" else
      rspec$primary_dpars[1L])
    if (!is.null(dpar_report_hook(rspec$family, dp, rspec))) return(NULL)
    lp <- fit$frame[["linpreds"]][[linpred_key(r, dp)]]
    lp[["link"]]
  }
  if (length(dim(x1)) == 2L) return(list(one(resp %||% rn[1L])))
  lapply(dimnames(x1)[[3L]], function(k) {
    if (grepl("^P[(]Y = ", k)) {
      return(structure(get_link("logit"), simplex = TRUE))
    }
    if (k %in% rn) return(one(k))
    NULL
  })
}

# The reference grid of a pooled fit, as emmeans builds it for a
# `mice::mira`: each imputation's grid data, with every numeric column
# averaged over the imputations. A refusal comes back as a string,
# which is emmeans's convention for one.
#' @exportS3Method emmeans::recover_data
recover_data.frmtmb_multiple <- function(object, ...) {
  rds <- lapply(object$fits, function(f) emmeans::recover_data(f, ...))
  rd <- rds[[1L]]
  if (is.character(rd)) return(rd)
  for (j in which(vapply(rd, is.numeric, NA))) {
    rd[[j]] <- Reduce(`+`, lapply(rds, `[[`, j)) / length(rds)
  }
  rd
}

# emmeans's basis of a pooled fit, by Rubin's rules over each
# imputation's basis on the one grid, as emmeans pools a `mice::mira`:
# the mean coefficient vector, the total covariance, and the
# Barnard-Rubin degrees of freedom of each contrast. The complete-data
# degrees of freedom are those one fit's basis reports, which is `Inf`
# (emmeans's asymptotic reference for a frmtmb fit), so the result is
# the one emmeans gives the same imputations fitted by `glm()`. The
# pooled table of `fixef()` takes `df.residual()` instead, as
# `mice::pool()` does.
#' @exportS3Method emmeans::emm_basis
emm_basis.frmtmb_multiple <- function(object, trms, xlev, grid, ...) {
  bas <- lapply(object$fits, function(f) {
    emmeans::emm_basis(f, trms, xlev, grid, ...)
  })
  m <- length(bas)
  b <- matrix(vapply(bas, function(z) as.numeric(z$bhat),
                     numeric(length(bas[[1L]]$bhat))), ncol = m)
  ok <- which(stats::complete.cases(b))
  B <- matrix(0, nrow(b), nrow(b))
  B[ok, ok] <- stats::cov(t(b[ok, , drop = FALSE]))
  W <- Reduce(`+`, lapply(bas, `[[`, "V")) / m
  Tv <- W + (1 + 1 / m) * B
  out <- bas[[1L]]
  out$bhat <- rowMeans(b)
  out$V <- Tv
  out$dfargs <- list(B = B, Tv = Tv, m = m, df1 = bas[[1L]]$dffun,
                     dfargs1 = bas[[1L]]$dfargs)
  out$dffun <- function(k, dfargs) {
    dfcom <- dfargs$df1(k, dfargs$dfargs1)
    bk <- sum(k * (dfargs$B %*% k))
    tk <- sum(k * (dfargs$Tv %*% k))
    lam <- max((1 + 1 / dfargs$m) * bk / tk, 1e-12)
    df_old <- (dfargs$m - 1) / lam^2
    if (!is.finite(dfcom)) return(df_old)
    df_obs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lam)
    df_old * df_obs / (df_old + df_obs)
  }
  out
}

#' @rdname frm_multiple-methods
#' @export
nobs.frmtmb_multiple <- function(object, ...) {
  frm_check_dots(...)
  n <- vapply(object$fits, function(f) as.numeric(nobs(f)), 0)
  if (length(unique(n)) != 1L) {
    frm_stop("nobs(): the imputations' fits have different numbers of ",
             "observations (", paste(n, collapse = ", "), "), so there is ",
             "no one count. Read nobs(x$fits[[j]])", call. = FALSE)
  }
  unname(n[1L])
}

#' @rdname frm_multiple-methods
#' @export
formula.frmtmb_multiple <- function(x, ...) {
  frm_check_dots(...)
  # every imputation is fitted with the one formula
  formula(x$fits[[1L]])
}

#' @rdname frm_multiple-methods
#' @export
family.frmtmb_multiple <- function(object, resp = NULL, ...) {
  frm_check_dots(...)
  family(object$fits[[1L]], resp = resp)
}

#' @export
model.frame.frmtmb_multiple <- function(formula, ...) {
  # stats' default would rebuild a frame from formula() and whatever
  # its environment holds, which is none of the imputed data sets
  multiple_refuse(
    "model.frame()",
    "each imputed data set has its own model frame",
    "Call model.frame(x$fits[[j]]) for imputation j")
}

#' The refusal of a method that has no Rubin's-rules meaning on a
#' pooled fit: what was asked, why the imputations do not pool into one
#' answer to it, and what to call instead.
#'
#' @noRd
multiple_refuse <- function(what, why, instead) {
  frm_stop(what, " has no pooled answer for a frm_multiple() result: ",
           why, ". ", instead, call. = FALSE)
}

# brms answers the in-sample calls below from the first imputation's
# data alone and warns "Using only the first imputed data set"; one
# imputation is what each refusal names instead

#' @export
predict.frmtmb_multiple <- function(object, ...) {
  multiple_refuse(
    "predict()",
    paste0("a predictive interval is a quantile of simulated responses, ",
           "and Rubin's rules pool an estimate and its standard error, ",
           "not a distribution of draws"),
    paste0("fitted() gives the pooled expected response, and ",
           "predict(x$fits[[j]]) one imputation's predictive interval"))
}

#' @export
residuals.frmtmb_multiple <- function(object, ...) {
  multiple_refuse(
    "residuals()",
    paste0("a residual is an observed response less a fitted value, and ",
           "the imputed data sets observe different values wherever one ",
           "was imputed"),
    "Call residuals(x$fits[[j]]) for imputation j")
}

#' @export
logLik.frmtmb_multiple <- function(object, ...) {
  multiple_refuse(
    "logLik() (and so AIC() and BIC())",
    paste0("the m log-likelihoods are of m different data sets, and ",
           "their average is the likelihood of no model"),
    paste0("Compare nested pooled fits with anova(m1, m0), which pools ",
           "the likelihood-ratio test (D3, or D1 and D2), or read ",
           "logLik(x$fits[[j]])"))
}

#' @exportS3Method bayesplot::pp_check
#' @export
pp_check.frmtmb_multiple <- function(object, ...) {
  multiple_refuse(
    "pp_check()",
    paste0("the check compares simulated responses with the observed ",
           "ones, and the imputed data sets observe different values ",
           "wherever one was imputed"),
    "Call pp_check(x$fits[[j]]) for imputation j")
}

#' @exportS3Method loo::loo
#' @export
loo.frmtmb_multiple <- function(x, ...) {
  multiple_refuse(
    "loo()",
    paste0("an elpd averages the likelihood over posterior draws, and ",
           "these are m maximum-likelihood fits"),
    paste0("Compare nested pooled fits with anova(m1, m0); for one ",
           "imputation, loo(frmtmb.sample::frm_sample(x$fits[[j]]))"))
}

#' @exportS3Method loo::waic
#' @export
waic.frmtmb_multiple <- function(x, ...) {
  multiple_refuse(
    "waic()",
    paste0("it averages the likelihood over posterior draws, and these ",
           "are m maximum-likelihood fits"),
    paste0("Compare nested pooled fits with anova(m1, m0); for one ",
           "imputation, waic(frmtmb.sample::frm_sample(x$fits[[j]]))"))
}

#' @export
update.frmtmb_multiple <- function(object, ...) {
  multiple_refuse(
    "update()",
    paste0("each imputation's fit stores a call that reads its data ",
           "under a name local to frm_multiple(), so no stored call ",
           "refits it"),
    paste0("Call frm_multiple() again with the new formula and the same ",
           "imputed data sets"))
}
