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
#'
#' @param formula,... As in [frm()].
#' @param data A list of completed data frames, or a `mice::mids`.
#' @param level Confidence level for the `$pooled_varcorr` interval.
#' @return A `frmtmb_multiple` object: `pooled` (the Rubin table for
#'   the fixed effects), `pooled_varcorr` (grp/term/type/estimate/
#'   lwr/upr/df/fmi for the random-effect SDs and correlations; `NULL`
#'   without random effects), and `fits` (the per-imputation fits).
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
  if (inherits(data, "mids")) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("A mids object needs the 'mice' package", call. = FALSE)
    }
    data <- mice::complete(data, action = "all")
  }
  if (!is.list(data) || length(data) < 2L ||
      !all(vapply(data, is.data.frame, TRUE))) {
    stop("data must be a list of at least two data frames (or a ",
         "mice::mids object)", call. = FALSE)
  }
  fits <- lapply(data, function(d) frm(formula, data = d, ...))
  m <- length(fits)
  dfcom <- df.residual(fits[[1]])

  nm <- estimated_coef_names(fits[[1]])
  cf <- vapply(fits, function(f) {
    bd <- f$estimates$betad
    if (length(fx <- f$frame$betad_fixed_idx)) bd <- bd[-fx]
    c(f$estimates$beta, bd)
  }, numeric(length(nm)))
  cf <- matrix(cf, nrow = length(nm), dimnames = list(nm, NULL))
  us <- vapply(fits, function(f) diag(vcov(f)), numeric(length(nm)))
  pl <- rubin_pool(cf, matrix(us, nrow = length(nm)), dfcom)
  tstat <- pl$estimate / pl$se
  tab <- data.frame(
    estimate = pl$estimate, se = pl$se, df = pl$df, t = tstat,
    p = 2 * stats::pt(-abs(tstat), pl$df), fmi = pl$fmi,
    row.names = nm
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
  method <- match.arg(method)
  use <- match.arg(use)
  others <- Filter(function(x) inherits(x, "frmtmb_multiple"), list(...))
  if (length(others) > 1L) {
    stop("anova() compares two frmtmb_multiple fits at a time",
         call. = FALSE)
  }
  # Same refusal as anova.frmtmb_fit: the restricted likelihood depends
  # on the fixed-effect design, so neither its value nor the Wald
  # covariance built from it compares across models.
  for (mf in c(list(object), others)) {
    if (any(vapply(mf$fits, `[[`, TRUE, "REML"))) {
      stop("Pooled model comparison requires ML fits (REML = FALSE)",
           call. = FALSE)
    }
  }
  wald <- method == "D1" || (method == "D2" && use == "wald")
  if (!length(others)) {
    if (is.null(constraint)) {
      stop("anova() needs a second frmtmb_multiple fit, or a ",
           "`constraint` naming the coefficients to test", call. = FALSE)
    }
    if (!wald) {
      stop("`constraint` needs a Wald rule: method = \"D1\", or ",
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
    st <- d1_stat(parts$Q, parts$U, m, dfcom)
  } else if (method == "D2") {
    st <- d2_stat(parts$d, parts$k, m)
  } else {
    st <- d3_stat(parts$dbar, parts$dtilde, parts$k, m)
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
    stop("anova() needs both fits pooled over the same imputations ",
         "(got m = ", a$m, " and ", b$m, ")", call. = FALSE)
  }
  # Same guard as anova.frmtmb_fit, applied imputation by imputation:
  # likelihoods computed on different data are not on a common scale.
  na <- vapply(a$fits, function(f) as.integer(f$frame$n_obs), 0L)
  nb <- vapply(b$fits, function(f) as.integer(f$frame$n_obs), 0L)
  if (!identical(na, nb)) {
    stop("anova() needs both fits fit to the same imputed datasets; ",
         "imputation ", which(na != nb)[1], " has ", na[which(na != nb)[1]],
         " observations in one fit and ", nb[which(na != nb)[1]],
         " in the other", call. = FALSE)
  }
  dfa <- attr(logLik(a$fits[[1]]), "df")
  dfb <- attr(logLik(b$fits[[1]]), "df")
  if (dfa == dfb) {
    stop("anova() needs nested fits of different size (both have ",
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
      stop("constraint names no coefficient of the fit: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
    tested <- nm1[nm1 %in% constraint]
    labels <- c(model_label(big$fits[[1]]),
                paste0("constraint: ",
                       paste(paste(tested, "= 0"), collapse = ", ")))
  } else {
    nm0 <- estimated_coef_names(small$fits[[1]])
    if (length(setdiff(nm0, nm1))) {
      stop("the smaller fit is not nested in the larger one: ",
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
      stop("the fits differ in their covariance parameters, which a ",
           "Wald rule cannot test. Use method = \"D3\" or ",
           "method = \"D2\" with use = \"likelihood\"", call. = FALSE)
    }
    tested <- setdiff(nm1, nm0)
    if (!length(tested)) {
      stop("the two fits have the same coefficients; a Wald rule has ",
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
    cf <- c(f$estimates$beta,
            if (length(fx <- f$frame$betad_fixed_idx)) {
              f$estimates$betad[-fx]
            } else {
              f$estimates$betad
            })
    names(cf) <- estimated_coef_names(f)
    Q[, j] <- cf[tested]
    U[, , j] <- vcov(f)[tested, tested, drop = FALSE]
  }
  # per-imputation Wald chi-squares, for D2(use = "wald")
  d <- vapply(seq_len(m), function(j) {
    drop(t(Q[, j]) %*% solve(U[, , j], Q[, j]))
  }, 0)
  list(Q = Q, U = U, d = d, k = k, labels = labels,
       big_fit = big$fits[[1]])
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
      stop("the imputations produced different parameter vectors; ",
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
#' `k` and the `d1_denom_df()` denominator degrees of freedom.
#'
#' @noRd
d1_stat <- function(Q, U, m, dfcom) {
  k <- nrow(Q)
  qbar <- rowMeans(Q)
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
hypothesis.frmtmb_multiple <- function(x, hypothesis, alpha = 0.05,
                                       class = NULL, group = NULL, ...) {
  if (...length()) {
    warning("ignoring arguments unused by pooled hypothesis tests: ",
            paste(...names(), collapse = ", "), call. = FALSE)
  }
  vo <- hyp_vals_only(x$fits[[1]])
  hp <- hyp_parse_all(hypothesis,
                      names(hyp_env_vals(x$fits[[1]], vo$vals, vo$comp)),
                      class, group)
  exs <- hp$exprs
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
        stop("Hypothesis '", hypothesis[i], "' must evaluate to a ",
             "single number in imputation ", j, " of ", m, call. = FALSE)
      }
      g <- hyp_fd_grad(function(v) hyp_eval(fit, ex, v, pc$comp),
                       pc$vals)
      Q[i, j] <- val
      U[i, j] <- max(0, drop(t(g) %*% pc$V %*% g))
    }
  }
  pl <- rubin_pool(Q, U, df.residual(x$fits[[1]]))
  rows <- lapply(seq_along(exs), function(i) {
    # the reference is the Barnard-Rubin t of this row, so the
    # quantile and tail functions go in per row
    df_i <- pl$df[i]
    wr <- hyp_wald_row(pl$estimate[i], pl$se[i], hp$dir[i], alpha,
                       function(p) stats::qt(p, df_i),
                       function(q) stats::pt(q, df_i))
    data.frame(hypothesis = hypothesis[i], estimate = pl$estimate[i],
               se = pl$se[i], lwr = wr$lwr, upr = wr$upr,
               t = wr$stat, df = df_i, p = wr$p)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "method") <- "wald"
  attr(out, "alpha") <- alpha
  attr(out, "direction") <- hp$dir
  class(out) <- c("frmtmb_hypothesis", "data.frame")
  out
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
  stop(fn, "() needs draws, and a frm_multiple() result has none: it ",
       "is m maximum-likelihood fits pooled by Rubin's rules, with no ",
       "chains. Read the pooled tables from `x$pooled` and ",
       "`x$pooled_varcorr` or test with hypothesis(), and use ",
       "frm_sample() on one imputation's fit (`x$fits[[1]]`) for draws",
       call. = FALSE)
}

#' @rdname as_draws
#' @exportS3Method posterior::as_draws
#' @export
as_draws.frmtmb_multiple <- function(x, ...) multiple_no_draws("as_draws")

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
  stop("plot() has no pooled display for a frm_multiple() result. ",
       "Plot one imputation's fit, `plot(x$fits[[1]])`, or print the ",
       "object for the pooled coefficients and variance components",
       call. = FALSE)
}

#' @exportS3Method brms::conditional_effects
#' @export
conditional_effects.frmtmb_multiple <- function(x, ...) {
  stop("conditional_effects() has no pooled version for a ",
       "frm_multiple() result: an effect curve would have to be ",
       "pooled across imputations, which is not implemented. Compute ",
       "it on one imputation's fit, `conditional_effects(x$fits[[1]])`",
       call. = FALSE)
}
