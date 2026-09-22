# brms's return shapes for a maximum-likelihood fit (item 2.6f).
#
# brms summarizes posterior draws and reports `Estimate`, `Est.Error`
# and two quantile columns. A maximum-likelihood fit has no draws, so
# the columns are kept and their CONTENT is the frequentist analogue:
# `Est.Error` is a standard error and the `Q` columns are the ends of a
# Wald interval at those probabilities. The alternative, inventing new
# column names, would break every ported brms script for a difference
# that the documentation can state once.

#' Validate brms's `probs` and give the column names it implies.
#'
#' brms writes `Q2.5`, `Q97.5`, `Q10`, `Q65`: the letter Q and the
#' probability in percent, with no padding and no trailing zeros. Two
#' probabilities is brms's default and not a rule; `probs = 0.65` is one
#' column and `probs = NULL` is none.
#'
#' @noRd
brms_prob_cols <- function(probs) {
  if (is.null(probs) || !length(probs)) return(character(0))
  if (!is.numeric(probs) || anyNA(probs) || any(probs < 0) ||
      any(probs > 1)) {
    frm_stop("`probs` must be probabilities in [0, 1], not ",
             arg_desc(probs), call. = FALSE)
  }
  paste0("Q", format(probs * 100, trim = TRUE, scientific = FALSE,
                     drop0trailing = TRUE))
}

#' brms's summary columns from a point estimate and a standard error.
#'
#' The Q columns are a Wald interval, `est + qnorm(p) * se`, which is
#' the interval [confint.frmtmb_fit()] reports for a coefficient. A row
#' whose standard error could not be recovered gets `NA` there rather
#' than a bound built from `NaN`.
#'
#' @noRd
brms_summary_matrix <- function(est, se, probs = c(0.025, 0.975),
                                rownames = NULL) {
  est <- as.numeric(est)
  se <- if (is.null(se)) rep(NA_real_, length(est)) else as.numeric(se)
  qs <- lapply(probs %||% numeric(0),
               function(p) est + stats::qnorm(p) * se)
  out <- cbind(Estimate = est, `Est.Error` = se)
  if (length(qs)) out <- cbind(out, do.call(cbind, qs))
  colnames(out) <- c("Estimate", "Est.Error", brms_prob_cols(probs))
  rownames(out) <- rownames
  out
}

#' The same summary for a matrix-valued prediction: one `Estimate`,
#' `Est.Error` and quantile layer per column of `m`, stacked into brms's
#' `n x ncol x K` array.
#'
#' @noRd
brms_summary_array <- function(m, se, probs = c(0.025, 0.975),
                               third = NULL) {
  m <- as.matrix(m)
  if (is.null(se)) se <- array(NA_real_, dim(m))
  cn <- c("Estimate", "Est.Error", brms_prob_cols(probs))
  out <- array(NA_real_, c(nrow(m), length(cn), ncol(m)),
               dimnames = list(rownames(m), cn, third %||% colnames(m)))
  for (k in seq_len(ncol(m))) {
    out[, , k] <- brms_summary_matrix(m[, k], se[, k], probs)
  }
  out
}

#' Summarize simulated draws the way brms summarizes posterior draws:
#' the mean (or median under `robust`), the spread, and the empirical
#' quantiles of the draws themselves.
#'
#' The quantiles here are NOT a Wald interval: `predict()` has real
#' draws, because it simulates, so its interval is the simulation's own
#' and carries the family's skew and its discreteness.
#'
#' @noRd
brms_summarize_draws <- function(d, probs = c(0.025, 0.975),
                                 robust = FALSE) {
  d <- as.matrix(d)
  # na.rm throughout: a simulated replicate that drew a non-finite
  # distributional parameter is NA in the draws matrix, and the
  # summary is taken over the replicates that are there. With no NAs
  # these are the same numbers.
  ctr <- if (robust) {
    function(x) stats::median(x, na.rm = TRUE)
  } else {
    function(x) mean(x, na.rm = TRUE)
  }
  spr <- if (robust) {
    function(x) stats::mad(x, na.rm = TRUE)
  } else {
    function(x) stats::sd(x, na.rm = TRUE)
  }
  est <- apply(d, 2L, ctr)
  err <- apply(d, 2L, spr)
  out <- cbind(Estimate = est, `Est.Error` = err)
  for (p in probs %||% numeric(0)) {
    out <- cbind(out, apply(d, 2L, stats::quantile, probs = p,
                            names = FALSE, na.rm = TRUE))
  }
  colnames(out) <- c("Estimate", "Est.Error", brms_prob_cols(probs))
  # brms leaves the row dimnames NULL on fitted(), residuals() and
  # predict() (measured, dev/shapes-rev-brmsref.rds)
  rownames(out) <- NULL
  out
}

#' The population-level parameters that are NOT rows of
#' `brms_coef_table()`: an ordinal fit's thresholds and its `cs()`
#' coefficients.
#'
#' Both live in `extra_names` components of the parameter template
#' rather than in `beta` or `betad`, so the coefficient machinery never
#' saw them and `fixef()` on `bf(ord ~ x) + cumulative()` reported one
#' row where brms reports three. brms treats them as population-level
#' parameters: `fixef()`, `vcov()` and `summary()$fixed` all carry
#' `Intercept[1]`, `Intercept[2]`, `x`, and `z[1]`, `z[2]` for `cs(z)`
#' (measured, `dev/shapes-p1-brmsord.R`).
#'
#' Each block names the template component it reads, the reported
#' values, and the map from the internal vector to them, which is the
#' identity for a `cs()` coefficient and the family's threshold map
#' otherwise. The covariance needs that map's derivative, so it is
#' carried rather than the values alone.
#'
#' The thresholds of a model with more than one ordinal response are
#' left out, as in `hyp_put_ordinal()`: one component holds them all
#' and nothing says which belongs to which response.
#'
#' @noRd
brms_extra_fixef <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  est <- fit$estimates
  out <- list()
  lp_dp <- function(lp) {
    paste(lp[["resp"]] %||% NA_character_, lp[["dpar"]], sep = ":")
  }
  ord_lps <- Filter(function(lp) {
    identical(brms_lp_family(fit, lp)[["type"]], "ordinal") &&
      identical(lp[["dpar"]], "mu")
  }, fit$frame[["linpreds"]])
  raw <- est[["tau_raw"]]
  if (length(tpl[["tau_raw"]]) && length(raw) && length(ord_lps) == 1L) {
    lp <- ord_lps[[1L]]
    fam <- brms_lp_family(fit, lp)
    map <- function(r) ord_threshold_values(fam, r)
    v <- map(as.numeric(raw))
    pre <- brms_lp_prefix(fit, lp)
    out[[length(out) + 1L]] <- list(
      comp = "tau_raw", cls = "b", is_int = TRUE, dp = lp_dp(lp),
      key = coef_block_key(fit, lp),
      names = paste0(brms_usc(pre, "Intercept"), "[", seq_along(v), "]"),
      values = v, raw = as.numeric(raw), map = map)
  }
  for (lp in fit$frame[["linpreds"]]) {
    for (ct in lp[["cs"]] %||% list()) {
      v <- as.numeric(est[[ct[["par"]]]])
      if (!length(v)) next
      lab <- brms_rename(sub("^cs", "", ct[["label"]]))
      pre <- brms_lp_prefix(fit, lp)
      out[[length(out) + 1L]] <- list(
        comp = ct[["par"]], cls = "bcs", is_int = FALSE, dp = lp_dp(lp),
        key = coef_block_key(fit, lp),
        names = paste0(brms_usc(pre, lab), "[", seq_along(v), "]"),
        values = v, raw = v, map = identity)
    }
  }
  out
}

#' brms's population-level rows: which rows of `brms_coef_table()` they
#' are, which come from `brms_extra_fixef()` instead, their brms names
#' with the class prefix dropped, and brms's own ORDER.
#'
#' The order is brms's Stan parameter order and is not frmtmb's
#' linear-predictor order: brms holds each predictor's intercept in a
#' scalar of its own ahead of the `b` vector, so every intercept comes
#' first, then the ordinary coefficients predictor by predictor, then
#' the spline (`bs_`), monotonic (`bsp_`) and category-specific
#' (`bcs_`) blocks. A ported script indexes `fixef(fit)[3, ]`, so the
#' order is part of the shape.
#'
#' `idx` indexes the coefficient table and is `NA` on an extra row;
#' `blk` and `pos` index `extra` and are `NA` on a coefficient row.
#'
#' @noRd
brms_fixef_rows <- function(fit) {
  tab <- brms_coef_table(fit)
  keep <- which(!tab$natural)
  ex <- brms_extra_fixef(fit)
  nm <- sub("^(b|bs|bsp|bcs)_", "", tab$brms[keep])
  cls <- sub("^(b|bs|bsp|bcs)_.*$", "\\1", tab$brms[keep])
  cls[!cls %in% c("b", "bs", "bsp", "bcs")] <- "b"
  # the internal name is the design column, which is where "(Intercept)"
  # is written; the brms name would also match a covariate called
  # Intercept
  is_int <- grepl("[(]Intercept[)]$", tab$internal[keep])
  dp <- paste(tab$resp[keep], tab$dpar[keep], sep = ":")
  idx <- keep
  blk <- rep(NA_integer_, length(keep))
  pos <- rep(NA_integer_, length(keep))
  # a threshold is an intercept of the predictor it belongs to and
  # sorts ahead of that predictor's other intercepts, which is where
  # brms puts it; the tie-break below is what places it
  tie <- rep(0L, length(keep))
  for (b in seq_along(ex)) {
    e <- ex[[b]]
    k <- length(e$names)
    nm <- c(nm, e$names)
    cls <- c(cls, rep(e$cls, k))
    is_int <- c(is_int, rep(e$is_int, k))
    dp <- c(dp, rep(e$dp, k))
    idx <- c(idx, rep(NA_integer_, k))
    blk <- c(blk, rep(b, k))
    pos <- c(pos, seq_len(k))
    tie <- c(tie, rep(-1L, k))
  }
  if (!length(nm)) {
    return(list(idx = integer(0), names = character(0),
                blk = integer(0), pos = integer(0), extra = ex))
  }
  rank <- match(cls, c("b", "bs", "bsp", "bcs"))
  dpr <- match(dp, unique(dp))
  ord <- order(rank, !is_int, dpr, tie, seq_along(nm))
  list(idx = idx[ord], names = nm[ord], blk = blk[ord], pos = pos[ord],
       extra = ex)
}

#' The estimates of every `brms_fixef_rows()` row, in that order.
#'
#' @noRd
brms_fixef_values <- function(fit, rows = brms_fixef_rows(fit)) {
  cf <- fixef_estimated(fit)
  out <- numeric(length(rows$names))
  co <- !is.na(rows$idx)
  out[co] <- cf[rows$idx[co]]
  for (i in which(!co)) {
    out[i] <- rows$extra[[rows$blk[i]]]$values[rows$pos[i]]
  }
  stats::setNames(out, rows$names)
}

#' The covariance of the `brms_fixef_rows()` rows, when some of them
#' are `brms_extra_fixef()` blocks.
#'
#' `vcov_estimated()` covers `beta` and `betad` only, so a joint matrix
#' has to come from somewhere that also holds `tau_raw` and the `cs()`
#' components. `hyp_par_cov()` is that place: it already assembles them
#' under both the maximum-likelihood and the REML parameterization, and
#' `hypothesis()` reads the same matrix.
#'
#' A threshold is a nonlinear function of the internal vector when the
#' family estimates increments, so the block is a delta method rather
#' than a subset. The derivative is taken numerically from the map the
#' FAMILY declares, so an ordinal family shipped by another package is
#' covered; a hand-written derivative for the four built-in ones would
#' not be.
#'
#' @noRd
brms_fixef_extra_vcov <- function(fit, rows) {
  pc <- hyp_par_cov(fit)
  tab <- brms_coef_table(fit)
  n_beta <- attr(tab, "n_beta")
  cpos <- function(cp) which(pc$comp == cp)
  A <- matrix(0, length(rows$names), length(pc$comp))
  for (i in seq_along(rows$names)) {
    j <- rows$idx[i]
    if (!is.na(j)) {
      p <- if (j <= n_beta) cpos("beta")[j] else cpos("betad")[j - n_beta]
      if (is.na(p)) return(NULL)
      A[i, p] <- 1
      next
    }
    e <- rows$extra[[rows$blk[i]]]
    p <- cpos(e$comp)
    if (length(p) != length(e$raw)) return(NULL)
    A[i, p] <- fd_gradient_row(e$map, e$raw, rows$pos[i])
  }
  V <- A %*% pc$V %*% t(A)
  dimnames(V) <- list(rows$names, rows$names)
  V
}

#' Row `k` of the Jacobian of `f` at `x`, by central differences.
#'
#' @noRd
fd_gradient_row <- function(f, x, k, eps = 1e-6) {
  g <- numeric(length(x))
  for (j in seq_along(x)) {
    h <- eps * max(1, abs(x[j]))
    xp <- xm <- x
    xp[j] <- x[j] + h
    xm[j] <- x[j] - h
    g[j] <- (f(xp)[k] - f(xm)[k]) / (2 * h)
  }
  g
}

#' brms's `pars` filter: keep the named rows, in the order given, and
#' error on a name the object does not have.
#'
#' brms matches exactly, so `pars = "Trt1"` does not also take
#' `Trt1:Age`.
#'
#' @noRd
brms_pars_filter <- function(x, pars, what) {
  if (is.null(pars)) return(x)
  if (!is.character(pars) || anyNA(pars)) {
    frm_stop("`pars` must be a character vector of row names, or NULL ",
             "for all of them, not ", arg_desc(pars), call. = FALSE)
  }
  hit <- match(pars, rownames(x))
  if (anyNA(hit)) {
    frm_stop(what, ": no such parameter(s): ",
             paste0("'", pars[is.na(hit)], "'", collapse = ", "),
             ". Available: ", paste(rownames(x), collapse = ", "),
             call. = FALSE)
  }
  x[hit, , drop = FALSE]
}

#' The outer parameter vector of a fit, and the inverse.
#'
#' `outer_par_map()` already names the outer vector in the order
#' `vcov(full = TRUE)` and `confint()` use, so a value vector in that
#' order can be written back into the estimate list component by
#' component. This is the seam a delta method or a joint draw needs: it
#' perturbs the parameters and re-reads whatever function of them is
#' wanted, instead of a second, hand-written derivative for each one.
#'
#' The cache is dropped with the values, because a stored `sdreport`
#' belongs to the estimates it was computed at.
#'
#' @noRd
fit_outer_vector <- function(fit) {
  map <- outer_par_map(fit)
  est <- fit$estimates
  out <- numeric(length(map$names))
  for (cp in unique(map$comp)) {
    k <- which(map$comp == cp)
    v <- est[[cp]]
    if (cp == "betad" && length(fx <- fit$frame[["betad_fixed_idx"]])) {
      v <- v[-fx]
    }
    out[k] <- as.numeric(v)
  }
  stats::setNames(out, map$names)
}

#' @noRd
fit_set_outer <- function(fit, v, map = outer_par_map(fit)) {
  est <- fit$estimates
  for (cp in unique(map$comp)) {
    k <- which(map$comp == cp)
    if (cp == "betad" && length(fx <- fit$frame[["betad_fixed_idx"]])) {
      pos <- setdiff(seq_along(est[[cp]]), fx)
      est[[cp]][pos] <- v[k]
    } else {
      est[[cp]][] <- v[k]
    }
  }
  fit$estimates <- est
  fit$cache <- new.env(parent = emptyenv())
  fit
}

#' Delta-method standard errors of a matrix-valued function of the outer
#' parameters, by central differences.
#'
#' Used where no analytic gradient exists: the ordinal and categorical
#' category probabilities, which depend on the thresholds and the
#' category-specific coefficients as well as on the linear predictor.
#' The step is relative to each parameter's own scale so that a
#' threshold and a log standard deviation are perturbed comparably, and
#' it falls back to the absolute step at a parameter estimated near
#' zero.
#'
#' Returns `NULL` when the covariance is not usable, so the caller
#' reports `NA` standard errors rather than bounds built from `NaN`.
#'
#' @noRd
fit_fd_se <- function(fit, f, eps = 1e-5) {
  V <- tryCatch(suppressWarnings(vcov(fit, full = TRUE)),
                error = function(e) NULL)
  if (is.null(V) || !all(is.finite(V))) return(NULL)
  map <- outer_par_map(fit)
  v0 <- fit_outer_vector(fit)
  m0 <- as.matrix(f(fit))
  p <- length(v0)
  # J is (n * K) x p: one column per outer parameter, flattened the way
  # as.vector() flattens the prediction, so the quadratic form below is
  # elementwise over the same flattening
  J <- matrix(NA_real_, length(m0), p)
  for (j in seq_len(p)) {
    h <- eps * max(1, abs(v0[j]))
    vp <- v0; vp[j] <- vp[j] + h
    vm <- v0; vm[j] <- vm[j] - h
    up <- as.vector(as.matrix(f(fit_set_outer(fit, vp, map))))
    dn <- as.vector(as.matrix(f(fit_set_outer(fit, vm, map))))
    J[, j] <- (up - dn) / (2 * h)
  }
  se <- sqrt(pmax(0, rowSums((J %*% V) * J)))
  matrix(se, nrow(m0), ncol(m0), dimnames = dimnames(m0))
}

#' brms's `P(Y = k)` labels: the response's own categories, or their
#' positions when the response carries no labels.
#'
#' @noRd
brms_category_labels <- function(labels, k) {
  paste0("P(Y = ", labels %||% as.character(seq_len(k)), ")")
}
