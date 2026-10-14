# --- simplex responses: dirichlet(), logistic_normal(), dirichlet2 -------
#
# brms's three families whose response is a row of proportions on the
# simplex. They share the categorical machinery's shape: one linear
# predictor per non-reference category, all receiving the main formula,
# category names read off the response. They differ from it in the
# response, an n x K matrix of proportions rather than a category index,
# which is why they are continuous families with a matrix response and
# carry their categories in `fam[["simplex"]]`.

#' brms's `make_stan_names()`, which turns `paste0(dpar, category)` into
#' the dpar names: dots and underscores removed after `make.names()`.
#'
#' @noRd
simplex_dpar_names <- function(prefix, cats) {
  out <- gsub("[._]", "", make.names(paste0(prefix, cats), unique = TRUE))
  if (anyDuplicated(out)) {
    frm_stop("Invalid response category names. Please avoid using any ",
             "special characters in the names. The categories ",
             paste(cats, collapse = ", "), " give the dpar names ",
             paste(out, collapse = ", "), call. = FALSE)
  }
  out
}

#' The response categories of a simplex family: the column names of the
#' response matrix, or `1..K` when it has none, as brms's
#' `extract_cat_names()` reads them.
#'
#' @noRd
simplex_categories <- function(formula, data, fname) {
  ri <- parse_response(formula)
  y <- tryCatch(eval(ri$resp, data, environment(formula) %||% globalenv()),
                error = function(e) NULL)
  if (is.null(y)) {
    frm_stop(fname, "(): the response '", deparse1(ri$resp), "' could not ",
             "be evaluated on the data to find its categories", call. = FALSE)
  }
  y <- as.matrix(y)
  cats <- colnames(y) %||% as.character(seq_len(ncol(y)))
  if (length(cats) < 2L) {
    frm_stop("At least 2 response categories are required. ", fname,
             "() takes an n x K matrix of proportions, each row summing ",
             "to 1, as cbind(y1, y2, y3) or a matrix column of the data; '",
             deparse1(ri$resp), "' has ", length(cats), " column",
             call. = FALSE)
  }
  cats
}

#' The position of the reference category among `cats`, or `NA` for none.
#' brms's rule: `NULL` is the first category, a name must be one of
#' them, and `NA` (dirichlet() only) predicts every category.
#'
#' @noRd
simplex_refcat <- function(refcat, cats) {
  if (is.null(refcat)) return(1L)
  if (is.na(refcat)) return(NA_integer_)
  if (!refcat %in% cats) {
    frm_stop("The reference response category must be one of ",
             paste(cats, collapse = ", "), ". Got refcat = '", refcat, "'",
             call. = FALSE)
  }
  match(refcat, cats)
}

#' Check a constructor's `refcat`: one string, or `NA` where the family
#' allows predicting every category. A number or a factor is taken as
#' its string, as brms's `as_one_character()` takes it, so `refcat = 2`
#' names the category `"2"` of an unnamed response.
#'
#' @noRd
check_simplex_refcat <- function(refcat, fname, allow_na) {
  if (is.null(refcat)) return(NULL)
  if (length(refcat) == 1L && !is.na(refcat) &&
        (is.numeric(refcat) || is.factor(refcat))) {
    refcat <- as.character(refcat)
  }
  ok <- length(refcat) == 1L &&
    (is.character(refcat) || (allow_na && is.logical(refcat))) &&
    (allow_na || !is.na(refcat))
  if (!ok) {
    frm_stop(fname, "(refcat =) names the reference category with one ",
             "string", if (allow_na) {
               ", or is NA to predict every category"
             } else {
               paste0("; ", fname, "() needs a reference category, because ",
                      "its log ratios are taken against one")
             }, ". Got ", arg_desc(refcat), call. = FALSE)
  }
  if (is.na(refcat)) NA else as.character(refcat)
}

#' The response check every simplex family makes, in brms's order and
#' with brms's sentences: a matrix of K columns, every row on the
#' simplex, every value strictly inside (0, 1).
#'
#' brms tests the row sums with `all.equal(rowSums(Y), 1)`, the MEAN
#' absolute deviation across rows against 1.5e-8, so one row far off the
#' simplex passes in a large enough data set. The test here is per row
#' at the same tolerance: it accepts no data set brms refuses, and it
#' refuses a row that is not a composition however many rows are.
#'
#' @noRd
simplex_valid_y <- function(y, K, fname) {
  if (!is.matrix(y)) {
    frm_stop("This model requires a response matrix. ", fname, "() takes ",
             "an n x ", K, " matrix of proportions, as cbind(y1, y2, ...) ",
             "or a matrix column of the data", call. = FALSE)
  }
  if (ncol(y) != K) {
    frm_stop(fname, "(): the response has ", ncol(y), " columns where the ",
             "family has ", K, " categories", call. = FALSE)
  }
  ok <- stats::complete.cases(y)
  dev <- abs(rowSums(y[ok, , drop = FALSE]) - 1)
  tol <- sqrt(.Machine$double.eps)
  if (any(dev > tol)) {
    bad <- which(ok)[dev > tol]
    frm_stop("Response values in simplex models must sum to 1. ",
             length(bad), " row", if (length(bad) > 1L) "s do" else " does",
             " not (row ", paste(utils::head(bad, 5L), collapse = ", "),
             if (length(bad) > 5L) ", ...", "; the largest departure is ",
             signif(max(dev), 3), "). Each row is held to 1.5e-8 of 1",
             call. = FALSE)
  }
  if (any(y[ok, ] <= 0)) {
    frm_stop("Family '", fname, "' requires response greater than 0. A ",
             "proportion of exactly 0 has density 0 under it (its log is ",
             "-Inf); replace the zeros, or model the counts behind them",
             call. = FALSE)
  }
  if (any(y[ok, ] >= 1)) {
    frm_stop("Family '", fname, "' requires response smaller than 1.",
             call. = FALSE)
  }
  invisible(NULL)
}

#' The link-scale predictor of every category as a list of K entries,
#' `NULL` at the reference category, whose predictor is the pinned 0.
#'
#' @noRd
simplex_eta_list <- function(dpars, sx) {
  out <- vector("list", length(sx$cats))
  for (j in seq_along(sx$pred)) out[[sx$pred[j]]] <- dpars[[sx$mu[j]]]
  out
}

#' Log of the softmax denominator over the categories, accumulated in
#' log space from the reference category's implicit zero, as
#' `cat_logit_lpdf()` does, so a wide predictor cannot overflow it.
#'
#' @noRd
simplex_log_denominator <- function(E) {
  lden <- NULL
  for (e in E) {
    if (is.null(e)) next
    lden <- if (is.null(lden)) e else RTMB::logspace_add(lden, e)
  }
  if (length(lden) && any(vapply(E, is.null, NA))) {
    lden <- RTMB::logspace_add(lden, 0 * lden)
  }
  lden
}

#' The `n x K` category probabilities from a NUMERIC dpar list, columns
#' named by category. The row maximum comes out before `exp()`.
#'
#' @noRd
simplex_probs <- function(dpars, sx, n = NULL) {
  E <- simplex_eta_list(dpars, sx)
  n <- n %||% max(vapply(E, length, 1L))
  M <- vapply(E, function(e) if (is.null(e)) numeric(n) else {
    rep(as.numeric(e), length.out = n)
  }, numeric(n))
  M <- matrix(M, n, length(sx$cats))
  P <- exp(M - apply(M, 1L, max))
  P <- P / rowSums(P)
  colnames(P) <- sx$cats
  P
}

#' Dirichlet draws by gamma normalization, one row per row of `alpha`.
#'
#' @noRd
simplex_rdirichlet <- function(alpha, cats) {
  G <- matrix(stats::rgamma(length(alpha), shape = alpha), nrow(alpha))
  out <- G / rowSums(G)
  colnames(out) <- cats
  out
}

#' Starting values on the response scale, read off the column means of
#' the response: additive log ratios against the reference category, a
#' moment estimate of the Dirichlet precision.
#'
#' @noRd
simplex_mean_alr <- function(y, k, ref) {
  m <- colMeans(y, na.rm = TRUE)
  if (is.na(ref)) log(m[k]) - mean(log(m)) else log(m[k]) - log(m[ref])
}

#' @noRd
simplex_phi_start <- function(y) {
  m <- colMeans(y, na.rm = TRUE)
  v <- apply(y, 2L, stats::var, na.rm = TRUE)
  ph <- stats::median(m * (1 - m) / v - 1)
  if (!is.finite(ph)) ph <- 1
  max(ph, 0.5)
}

#' The common part of a simplex family: its category bookkeeping and the
#' fields every route reads.
#'
#' @noRd
simplex_info <- function(kind, cats, ref) {
  K <- length(cats)
  pred <- if (is.na(ref)) seq_len(K) else seq_len(K)[-ref]
  list(kind = kind, cats = cats, ref = ref, pred = pred,
       mu = simplex_dpar_names("mu", cats[pred]))
}

#' Attach the simplex bookkeeping to a built family.
#'
#' @noRd
simplex_finish <- function(fam, sx, joint_link) {
  fam[["simplex"]] <- sx
  fam[["cat_levels"]] <- sx$cats
  fam[["cat_K"]] <- length(sx$cats)
  if (!is.null(joint_link)) fam[["joint_link"]] <- joint_link
  fam
}

#' Dirichlet family over the categories `cats`, the reference at position
#' `ref` (`NA` for none). brms's `dirichlet_logit_lpdf()`:
#' `y ~ Dirichlet(phi * softmax(mu))`.
#'
#' @noRd
fam_dirichlet_impl <- function(cats, ref, link_phi = "log") {
  sx <- simplex_info("dirichlet", cats, ref)
  K <- length(cats)
  mun <- sx$mu
  lk_phi <- dpar_link(link_phi, "phi", "dirichlet", dpar_links_positive)
  init <- stats::setNames(lapply(sx$pred, function(k) {
    force(k)
    function(y, aterms) simplex_mean_alr(y, k, ref)
  }), mun)
  init[["phi"]] <- function(y, aterms) simplex_phi_start(y)
  fam <- frmtmb_family(
    "dirichlet",
    accepts_aterms = "weights",
    dpars = c(mun, "phi"),
    links = c(stats::setNames(rep(list("identity"), length(mun)), mun),
              list(phi = lk_phi)),
    lpdf = function(y, dpars, aterms) {
      ly <- log(y)
      E <- simplex_eta_list(dpars, sx)
      lden <- simplex_log_denominator(E)
      phi <- dpars[["phi"]]
      ll <- lgamma(phi)
      for (k in seq_len(K)) {
        lp <- if (is.null(E[[k]])) -lden else E[[k]] - lden
        a <- phi * exp(lp)
        ll <- ll - lgamma(a) + (a - 1) * ly[, k]
      }
      ll
    },
    valid_y = function(y, aterms) simplex_valid_y(y, K, "dirichlet"),
    init_dpars = init,
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) simplex_probs(dpars, sx)
    ),
    sim = function(dpars, aterms, n) {
      P <- simplex_probs(dpars, sx, n)
      simplex_rdirichlet(P * rep(dpars[["phi"]], length.out = n), cats)
    },
    primary_dpars = mun
  )
  simplex_finish(fam, sx, "logit")
}

#' brms's internal `dirichlet2`, reached as `brmsfamily("dirichlet2")`:
#' `y ~ Dirichlet(alpha)`, one predictor per category on the scale of its
#' concentration, no reference category and no precision.
#'
#' @noRd
fam_dirichlet2_impl <- function(cats, link = "log") {
  sx <- simplex_info("dirichlet2", cats, NA_integer_)
  K <- length(cats)
  mun <- sx$mu
  lk <- mu_link(link, "dirichlet2")
  init <- stats::setNames(lapply(seq_len(K), function(k) {
    force(k)
    function(y, aterms) {
      colMeans(y, na.rm = TRUE)[k] * simplex_phi_start(y)
    }
  }), mun)
  fam <- frmtmb_family(
    "dirichlet2",
    accepts_aterms = "weights",
    dpars = mun,
    links = stats::setNames(rep(list(lk), K), mun),
    lpdf = function(y, dpars, aterms) {
      ly <- log(y)
      a0 <- 0
      ll <- 0
      for (k in seq_len(K)) {
        a <- dpars[[mun[k]]]
        a0 <- a0 + a
        ll <- ll - lgamma(a) + (a - 1) * ly[, k]
      }
      ll + lgamma(a0)
    },
    valid_y = function(y, aterms) simplex_valid_y(y, K, "dirichlet2"),
    init_dpars = init,
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        A <- dirichlet2_alpha(dpars, mun)
        P <- A / rowSums(A)
        colnames(P) <- cats
        P
      }
    ),
    sim = function(dpars, aterms, n) {
      simplex_rdirichlet(dirichlet2_alpha(dpars, mun, n), cats)
    },
    primary_dpars = mun
  )
  simplex_finish(fam, sx, NULL)
}

#' The `n x K` concentrations of a numeric `dirichlet2` dpar list.
#'
#' @noRd
dirichlet2_alpha <- function(dpars, mun, n = NULL) {
  n <- n %||% max(vapply(mun, function(m) length(dpars[[m]]), 1L))
  matrix(vapply(mun, function(m) rep(as.numeric(dpars[[m]]), length.out = n),
                numeric(n)), n, length(mun))
}

#' Logistic-normal family: the additive log ratios against the reference
#' category, `log(y_k / y_ref)`, are multivariate normal with means
#' `mu<k>`, scales `sigma<k>` and correlation matrix `lncor`. brms's
#' `logistic_normal_cholesky_cor_lpdf()`, which adds the Jacobian
#' `-sum(log(y))` of the log-ratio map.
#'
#' The correlations are a family-level parameter, `lncor_raw`, held as the
#' row-normalized Cholesky chart the residual correlations of
#' `set_rescor(TRUE)` use, so the LKJ prior and its Jacobian carry over.
#'
#' @noRd
fam_logistic_normal_impl <- function(cats, ref, link_sigma = "log") {
  sx <- simplex_info("logistic_normal", cats, ref)
  K <- length(cats)
  D <- K - 1L
  mun <- sx$mu
  sgn <- simplex_dpar_names("sigma", cats[sx$pred])
  sx$sigma <- sgn
  lk_sigma <- dpar_link(link_sigma, "sigma", "logistic_normal",
                        dpar_links_positive)
  alr <- function(y) {
    ly <- log(y)
    ly[, sx$pred, drop = FALSE] - ly[, ref]
  }
  init <- stats::setNames(lapply(seq_len(D), function(j) {
    force(j)
    function(y, aterms) mean(alr(y)[, j], na.rm = TRUE)
  }), mun)
  init <- c(init, stats::setNames(lapply(seq_len(D), function(j) {
    force(j)
    function(y, aterms) max(stats::sd(alr(y)[, j], na.rm = TRUE), 1e-3)
  }), sgn))
  fam <- frmtmb_family(
    "logistic_normal",
    accepts_aterms = "weights",
    dpars = c(mun, sgn),
    links = c(stats::setNames(rep(list("identity"), D), mun),
              stats::setNames(rep(list(lk_sigma), D), sgn)),
    lpdf = function(y, dpars, aterms, extra) {
      "c" <- RTMB::ADoverload("c")
      n <- nrow(y)
      ly <- log(y)
      zs <- NULL
      lsig <- 0
      for (j in seq_len(D)) {
        s <- dpars[[sgn[j]]]
        z <- (ly[, sx$pred[j]] - ly[, ref] - dpars[[mun[j]]]) / s
        if (length(z) == 1L) z <- z + numeric(n)
        zs <- if (is.null(zs)) z else c(zs, z)
        lsig <- lsig + log(s)
      }
      ll <- if (D == 1L) {
        RTMB::dnorm(zs, 0, 1, log = TRUE)
      } else {
        C <- us_chol_cor(extra[["lncor_raw"]], D)
        RTMB::dmvnorm(RTMB::matrix(zs, n, D), 0, C, log = TRUE)
      }
      ll - lsig - rowSums(ly)
    },
    valid_y = function(y, aterms) {
      simplex_valid_y(y, K, "logistic_normal")
    },
    init_dpars = init,
    type = "continuous",
    post = list(
      # brms refuses it in these words: the mean of a logistic-normal
      # composition has no closed form
      mean_fn = function(dpars, aterms) {
        frm_stop("Cannot compute expected values of the posterior ",
                 "predictive distribution for family 'logistic_normal'. ",
                 "Its mean has no closed form, as in brms; predict() and ",
                 "simulate() give the predictive distribution, and ",
                 "frm_linpred(dpar = ) a category's log-ratio mean",
                 call. = FALSE)
      }
    ),
    extra_pars = if (D > 1L) {
      function(y, aterms) {
        X <- alr(y)
        X <- X[stats::complete.cases(X), , drop = FALSE]
        sds <- if (nrow(X) > D) apply(X, 2L, stats::sd) else 0
        # a log ratio that does not vary (a template's placeholder
        # response) has no correlation to start from
        C <- if (all(is.finite(sds) & sds > 0)) stats::cor(X) else diag(D)
        th <- tryCatch(us_theta_cor(C), error = function(e) NULL)
        if (is.null(th) || any(!is.finite(th))) th <- numeric(D * (D - 1) / 2)
        list(lncor_raw = th)
      }
    },
    sim = function(dpars, aterms, n, extra) {
      C <- if (D == 1L) matrix(1, 1, 1) else {
        lncor_matrix(extra[["lncor_raw"]], D)
      }
      W <- matrix(stats::rnorm(n * D), n, D) %*% chol(C)
      E <- matrix(0, n, K)
      for (j in seq_len(D)) {
        E[, sx$pred[j]] <- rep(as.numeric(dpars[[mun[j]]]), length.out = n) +
          rep(as.numeric(dpars[[sgn[j]]]), length.out = n) * W[, j]
      }
      P <- exp(E - apply(E, 1L, max))
      P <- P / rowSums(P)
      colnames(P) <- cats
      P
    },
    primary_dpars = mun
  )
  simplex_finish(fam, sx, "identity")
}

#' The logistic-normal correlation matrix from its chart, numerically.
#'
#' @noRd
lncor_matrix <- function(theta, D) {
  if (D < 2L) return(matrix(1, 1, 1))
  C <- us_chol_cor(as.numeric(theta), D)
  C <- (C + t(C)) / 2
  diag(C) <- 1
  C
}

#' The correlations of a logistic-normal fit in brms's order and under
#' brms's names, `lncor__b__c`: the upper triangle read column by column,
#' which is the loop brms's generated quantities write it with.
#'
#' @noRd
lncor_values <- function(theta, sx) {
  D <- length(sx$pred)
  if (D < 2L) return(stats::setNames(numeric(0), character(0)))
  C <- lncor_matrix(theta, D)
  lab <- sub("^mu", "", sx$mu)
  idx <- which(upper.tri(C), arr.ind = TRUE)
  idx <- idx[order(idx[, "col"], idx[, "row"]), , drop = FALSE]
  stats::setNames(C[idx], paste0("lncor__", lab[idx[, "row"]], "__",
                                 lab[idx[, "col"]]))
}

#' Placeholder a bare simplex constructor returns: the categories are a
#' property of the data, so `frm()` swaps this for the real family
#' through `resolve_deferred_families()` once it can see the response,
#' as it does for `categorical()`.
#'
#' `links` names one stand-in per kind of parameter, `mu` and the
#' family's `phi` or `sigma`, so that the placeholder answers brms's
#' `$link` and `$link_phi` before the categories are known.
#'
#' @noRd
fam_simplex_deferred <- function(fname, build, links, joint_link) {
  fam <- frmtmb_family(
    fname,
    accepts_aterms = "weights",
    dpars = names(links),
    links = links,
    lpdf = function(y, dpars, aterms) {
      frm_stop(fname, "(): the response categories were never resolved. ",
               "They are read from the response matrix's column names by ",
               "frm(); pass data on the entry point you are using",
               call. = FALSE)
    },
    type = "continuous"
  )
  fam[["defer"]] <- function(formula, data) {
    build(simplex_categories(formula, data, fname))
  }
  fam[["simplex"]] <- list(kind = fname, cats = NULL)
  if (!is.null(joint_link)) fam[["joint_link"]] <- joint_link
  fam
}

#' @noRd
fam_dirichlet_deferred <- function(link = "logit", link_phi = "log",
                                   refcat = NULL) {
  mu_link(link, "dirichlet")
  refcat <- check_simplex_refcat(refcat, "dirichlet", allow_na = TRUE)
  # checked now and not only once the categories are known, so a bad link
  # is not masked by a data problem
  dpar_link(link_phi, "phi", "dirichlet", dpar_links_positive)
  fam_simplex_deferred("dirichlet", function(cats) {
    fam_dirichlet_impl(cats, simplex_refcat(refcat, cats), link_phi)
  }, list(mu = "identity", phi = link_phi), "logit")
}

#' @noRd
fam_dirichlet2_deferred <- function(link = "log") {
  mu_link(link, "dirichlet2")
  fam_simplex_deferred("dirichlet2", function(cats) {
    fam_dirichlet2_impl(cats, link)
  }, list(mu = link), NULL)
}

#' @noRd
fam_logistic_normal_deferred <- function(link = "identity",
                                         link_sigma = "log", refcat = NULL) {
  mu_link(link, "logistic_normal")
  refcat <- check_simplex_refcat(refcat, "logistic_normal", allow_na = FALSE)
  dpar_link(link_sigma, "sigma", "logistic_normal", dpar_links_positive)
  fam_simplex_deferred("logistic_normal", function(cats) {
    fam_logistic_normal_impl(cats, simplex_refcat(refcat, cats), link_sigma)
  }, list(mu = "identity", sigma = link_sigma), "identity")
}

#' Whether a family's response is a composition on the simplex.
#'
#' @noRd
fam_is_simplex <- function(fam) !is.null(fam[["simplex"]])

#' @section Compositional responses:
#' `dirichlet()` and `logistic_normal()` model a row of proportions that
#' sums to one, brms's simplex families. Write the response as an
#' `n x K` matrix: `cbind(y1, y2, y3)`, or a matrix column of the data
#' (`dd$Y <- Y`). Every row must sum to one, to within 1.5e-8, and every
#' proportion must lie strictly between 0 and 1. The column names are
#' the categories, and the FIRST is the reference category unless
#' `refcat` names another, as in brms.
#'
#' Each non-reference category `k` gets its own linear predictor
#' `mu<k>`, and the main formula applies to every one of them;
#' `bf(Y ~ x, muc ~ z)` gives category `c` a predictor of its own. A
#' group term `(1 | g)` gives each category its own effects, and
#' `(1 | ID | g)` correlates them across categories, as for
#' `categorical()`.
#'
#' `dirichlet()` is `y ~ Dirichlet(phi * softmax(mu))`: the category
#' means are the softmax of the predictors, with the reference
#' category's held at zero, and `phi` is the precision, which
#' `bf(phi ~ z)` models. `refcat = NA` predicts every category, which
#' brms allows; without a prior the predictors are then identified
#' only up to a common shift, and the fit says so.
#'
#' `logistic_normal()` makes the log ratios against the reference
#' category, `log(y_k / y_ref)`, multivariate normal with means
#' `mu<k>`, standard deviations `sigma<k>` and correlations `lncor`.
#' Each `sigma<k>` takes its own formula; `lncor` is one correlation
#' matrix for the model, reported as `lncor(b,c)` in the summary, and
#' takes an `lkj()` prior under `class = "lncor"`.
#'
#' `fitted()` on a `dirichlet()` fit returns brms's `n x 4 x K` array of
#' category means, layers named `P(Y = k)`. A logistic-normal mean has
#' no closed form, so `fitted()` and `conditional_effects()` are refused
#' on it, as in brms. On both, `predict()` and `simulate()` draw
#' compositions (rows on the simplex), `frm_linpred(dpar = "mub")` gives
#' one category's predictor, and `residuals()` is refused, as brms
#' refuses it for every family whose response is a set of categories.
#'
#' brms also has `dirichlet2`, not exported there, reached as
#' `brmsfamily("dirichlet2")`: one predictor per category on the scale
#' of its Dirichlet concentration (log link), and no precision. The
#' same spelling reaches it here.
#' @name frmtmb-families
NULL

#' @rdname frmtmb-families
#' @param refcat For `dirichlet()` and `logistic_normal()`: the name of
#'   the reference category, whose predictor is held at zero. `NULL`
#'   (the default) takes the first column of the response. `NA` predicts
#'   every category, which `dirichlet()` alone allows.
#' @export
dirichlet <- function(link = "logit", link_phi = "log", refcat = NULL) {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["dirichlet"]], "logit")
  fam_dirichlet_deferred(link, link_phi, refcat)
}

#' @rdname frmtmb-families
#' @export
logistic_normal <- function(link = "identity", link_sigma = "log",
                            refcat = NULL) {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["logistic_normal"]], "identity")
  fam_logistic_normal_deferred(link, link_sigma, refcat)
}

#' The logistic-normal correlation blocks of a fit: the template
#' component holding each one's chart and the family's categories.
#' `extra_pars` refuses a family with extras in a multivariate model, so
#' there is at most one, under its own name.
#'
#' @noRd
lncor_blocks <- function(fit) {
  tpl <- fit$frame[["par_template"]]
  out <- list()
  for (r in names(fit$spec$responses)) {
    sx <- fit$spec$responses[[r]]$family[["simplex"]]
    if (!identical(sx$kind, "logistic_normal")) next
    nm <- extra_tpl_name(fit$frame, r, "lncor_raw")
    if (!length(tpl[[nm]])) next
    out[[length(out) + 1L]] <- list(resp = r, comp = nm, sx = sx)
  }
  out
}

#' The template components a hypothesis and the summary read beside the
#' coefficients: the logistic-normal correlation charts.
#'
#' @noRd
simplex_extra_comps <- function(fit) {
  vapply(lncor_blocks(fit), `[[`, "", "comp")
}

#' Put the logistic-normal correlations under brms's names,
#' `lncor__b__c`, for `hyp_env_vals()`.
#'
#' @noRd
lncor_put <- function(fit, vals, comp, put) {
  for (b in lncor_blocks(fit)) {
    v <- lncor_values(vals[comp == b$comp], b$sx)
    for (j in seq_along(v)) put(names(v)[j], v[[j]])
  }
  invisible(NULL)
}

#' The summary rows of the logistic-normal correlations, brms's
#' `lncor(b,c)` under Further Distributional Parameters. The standard
#' error is the delta method through the chart; the interval is the
#' Wald interval of `atanh(r)` mapped back, so it stays inside (-1, 1).
#'
#' @noRd
lncor_summary_rows <- function(object, q) {
  bl <- lncor_blocks(object)
  if (!length(bl)) return(NULL)
  pc <- tryCatch(suppressWarnings(hyp_par_cov(object)),
                 error = function(e) NULL)
  out <- list(est = numeric(0), err = numeric(0), lo = numeric(0),
              hi = numeric(0), nm = character(0))
  for (b in bl) {
    th <- object$estimates[[b$comp]]
    r <- lncor_values(th, b$sx)
    G <- matrix(0, length(r), length(th))
    for (k in seq_along(th)) {
      h <- 1e-6 * max(1, abs(th[k]))
      tp <- th
      tm <- th
      tp[k] <- tp[k] + h
      tm[k] <- tm[k] - h
      G[, k] <- (lncor_values(tp, b$sx) - lncor_values(tm, b$sx)) / (2 * h)
    }
    se <- rep(NA_real_, length(r))
    if (!is.null(pc)) {
      pos <- which(pc$comp == b$comp)
      if (length(pos) == length(th)) {
        Gf <- matrix(0, length(r), length(pc$comp))
        Gf[, pos] <- G
        se <- sqrt(as.numeric(hyp_prop_var(pc, Gf)))
      }
    }
    z <- atanh(r)
    sz <- se / (1 - r^2)
    lab <- sub("^lncor__", "", names(r))
    out$est <- c(out$est, unname(r))
    out$err <- c(out$err, se)
    out$lo <- c(out$lo, tanh(z - q * sz))
    out$hi <- c(out$hi, tanh(z + q * sz))
    out$nm <- c(out$nm, paste0("lncor(", sub("__", ",", lab, fixed = TRUE),
                               ")"))
  }
  out
}

#' The one refusal of `residuals()` on a simplex family, for every type
#' and for new data: brms defines no residual of a composition.
#'
#' @noRd
simplex_residuals_refusal <- function(fam) {
  frm_stop("residuals() is not defined for the '", fam[["family"]],
           "' family, as in brms: the response is a row of proportions ",
           "over categories, which has no scalar residual. Compare ",
           "predict(fit), a summary of simulated compositions, with the ",
           "observed proportions instead", call. = FALSE)
}
