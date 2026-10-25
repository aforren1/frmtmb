# sar(), brms's spatial simultaneous autoregressive term, registered
# with core through frmtmb_register_autocor() (R/zzz.R). It is a
# residual correlation term, not a covariance structure: it changes the
# density of the response and, in its lag form, the mean, so it cannot
# be a latent field added to the linear predictor
# (frmtmb_register_covstruct()).

#' Spatial simultaneous autoregressive terms: `sar()`
#'
#' `sar(M, type = )` in a formula makes the residuals of a `gaussian()`
#' or `student()` response one spatially autoregressive draw over a
#' spatial weights matrix `M`, with brms's spelling and its two types.
#' It is available once frmtmb.spatial is loaded:
#' `library(frmtmb.spatial)`. Without the package, [frmtmb::frm()]
#' refuses the term and names this package.
#'
#' With \eqn{A = I - \rho W}, the two types are, as in brms 2.23.0:
#' \describe{
#'   \item{`type = "lag"` (the default)}{SAR of the response,
#'     \eqn{y = \rho W y + \eta + e}, so \eqn{A y - \eta = e}. The mean of
#'     the response is \eqn{A^{-1} \eta}, which [fitted()] reports, as
#'     brms's `posterior_epred()` does. The parameter is `lagsar`.}
#'   \item{`type = "error"`}{SAR of the residuals,
#'     \eqn{y = \eta + u}, \eqn{u = \rho W u + e}, so
#'     \eqn{A (y - \eta) = e}. The mean is \eqn{\eta}. The parameter is
#'     `errorsar`.}
#' }
#' `e` holds independent normal residuals with standard deviation
#' `sigma`, or for `student()` one multivariate-t draw with scale
#' `sigma` and `nu`. The log-likelihood has the log-determinant
#' \eqn{\log|I - \rho W| = \sum_i \log(1 - \rho \lambda_i)}, with
#' \eqn{\lambda_i} the eigenvalues of `W`, computed once.
#'
#' @section Arguments of the term:
#' \describe{
#'   \item{`M`}{The spatial weights matrix, one row and column per row
#'     of the data, in `data2 = list(M = M)`: a matrix, a sparse
#'     `Matrix`, or a `listw` or `nb` object of the spdep package, which
#'     is read as brms reads it (`spdep::listw2mat()`, and
#'     `spdep::nb2mat()` with its row-standardized default) without
#'     needing spdep. The rows the model drops (missing values) are
#'     dropped from `M`, as brms drops them. `M` is not checked for row
#'     standardization or symmetry, as in brms. Its eigenvalues must be
#'     real, with at least one negative and one positive.}
#'   \item{`type`}{`"lag"` or `"error"`. A unique abbreviation is read
#'     as brms reads it (`match.arg()`).}
#' }
#'
#' @section The range of rho:
#' \eqn{I - \rho W} is invertible and its determinant positive on the
#' interval \eqn{(1 / \lambda_{min}, 1 / \lambda_{max})}, the default
#' interval of `spatialreg::lagsarlm()`. frmtmb estimates \eqn{\rho} on
#' that interval through a logistic map, starting at 0. brms 2.23.0
#' bounds \eqn{\rho} by \eqn{(\lambda_{min}, \lambda_{max})} instead.
#' For a row-standardized `W` (\eqn{\lambda_{max} = 1}) that is a
#' narrower interval with the same upper end, which cuts off the
#' negative values below \eqn{\lambda_{min}}. For a `W` that is not
#' row-standardized, brms's interval reaches past
#' \eqn{1 / \lambda_{max}}, where the determinant changes sign; there
#' brms's log-determinant is the log of a negative number and Stan
#' rejects the proposal, so its support is frmtmb's interval there.
#'
#' The consequence on a row-standardized `W`: when the data favor rho
#' below the smallest eigenvalue, frmtmb and brms fit different models.
#' On a row-standardized 8 x 8 queen lattice (eigenvalues -0.5006 to 1)
#' with data drawn at rho = -0.8, frmtmb's estimate is -1.103, with a
#' log-likelihood 3.74 above the best inside brms's bounds, and brms's
#' posterior piles at its bound, -0.5006 (median -0.443; a fifth of the
#' draws within 0.02 of it).
#'
#' A `W` whose extreme eigenvalue is repeated in all but one direction,
#' as on a complete graph, gives a likelihood with no maximum at that
#' end when the model has an intercept; frmtmb warns.
#'
#' @section What else the term does:
#' Only `gaussian()` and `student()` take it, with a constant `sigma`
#' and `nu`, as in brms ("SAR models are not implemented when
#' predicting 'sigma' or 'nu'"). It is written in the formula of `mu`.
#' The likelihood is one joint density over all the rows, so
#' `weights()`, `cens()`, `trunc()`, `se()`, `mi()`, `rescor = TRUE`,
#' mixtures and `quadrature` are refused, and a response carries one
#' residual correlation term. Random effects work.
#'
#' `summary()` and [frmtmb::confint_varcorr()] report `lagsar`
#' or `errorsar` with a delta-method interval. [fitted()] and
#' `residuals()` use the mean above; pearson residuals divide by the
#' marginal standard deviation, \eqn{\sigma \sqrt{[(A'A)^{-1}]_{ii}}}.
#' [simulate()] and [predict()] draw the whole response through the
#' model. [frmtmb::autocor_matrix()] returns \eqn{(A'A)^{-1}}.
#' [fitted()], [predict()] and [simulate()] take `newdata` with one row
#' per row of `M`, read as the locations in their order, as brms takes
#' it: a counterfactual such as `x + 1` everywhere gives the lag model's
#' spillover, \eqn{A^{-1} \eta_{new}}. When the fit dropped rows for
#' missing values, `newdata` covers them too, and A is then over every
#' row of `M`, at the fitted rho; a missing predictor gives `NA`, which
#' the lag form spreads to every row, as in brms. Other `newdata` is
#' refused; brms refuses it too unless a matrix over the new rows comes
#' in its `newdata2`.
#' [frmtmb::frm_linpred()] gives the linear predictor \eqn{\eta} on any
#' new data.
#'
#' brms's class priors `lagsar` and `errorsar` work:
#' `set_prior("normal(0, 0.5)", class = "lagsar")` is a density on rho,
#' carried to its internal parameter with the log-Jacobian of the map.
#' frmtmb.sample's default is brms's, flat on rho over its interval.
#' The pointwise `log_lik()` is refused, as for frmtmb's other residual
#' covariance terms; brms's is a leave-one-out conditional density.
#'
#' @return `sar()` returns the term as an object of class
#'   `frmtmb_ac_term`. Its use is inside a formula, where [frmtmb::frm()]
#'   reads it.
#' @param M The spatial weights matrix; see "Arguments of the term".
#' @param type `"lag"` or `"error"`.
#' @seealso [frmtmb::frmtmb_register_autocor()] for the seam the term is
#'   registered through, and [frmtmb::frmtmb-autocor] for frmtmb's other
#'   residual correlation terms, `fcor()` among them.
#' @examples
#' set.seed(2)
#' grid <- expand.grid(r = 1:6, c = 1:6)
#' n <- nrow(grid)
#' A <- outer(seq_len(n), seq_len(n), function(i, j) {
#'   as.numeric(abs(grid$r[i] - grid$r[j]) + abs(grid$c[i] - grid$c[j]) == 1)
#' })
#' W <- A / rowSums(A)
#' d <- data.frame(x = rnorm(n))
#' d$y <- as.vector(solve(diag(n) - 0.5 * W, 1 + 0.5 * d$x + rnorm(n)))
#' fit <- frm(bf(y ~ x + sar(W, type = "lag")), family = gaussian(),
#'            data = d, data2 = list(W = W))
#' summary(fit)
#' @export
sar <- function(M, type = "lag") {
  cl <- sys.call()
  cl[[1L]] <- as.name("sar")
  structure(list(term = cl, label = deparse1(cl)),
            class = "frmtmb_ac_term")
}

#' The spec of a `sar()` call: `M` as written, and the type, matched as
#' brms matches it.
#'
#' @noRd
sar_parse <- function(args, env) {
  if (is.null(args[["M"]])) {
    frm_stop("Argument 'M' is missing in sar().", call. = FALSE)
  }
  type <- args[["type"]]
  if (is.null(type)) {
    type <- "lag"
  } else if (!is.character(type)) {
    type <- tryCatch(eval(type, env), error = function(e) {
      frm_stop("sar(): type = ", deparse1(args[["type"]]), " could not be ",
               "evaluated: ", conditionMessage(e), call. = FALSE)
    })
  }
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    frm_stop("sar(): type must be \"lag\" or \"error\"", call. = FALSE)
  }
  type <- tryCatch(match.arg(type, c("lag", "error")), error = function(e) {
    frm_stop("sar(): type must be \"lag\" or \"error\", not \"", type, "\"",
             call. = FALSE)
  })
  list(M = args[["M"]], type = type)
}

#' The weights matrix of a `sar()` term over the kept rows, with its
#' eigenvalues and the interval of rho they give.
#'
#' @noRd
sar_build <- function(spec, rows, n_all, lookup) {
  M <- lookup(spec[["M"]], "sar()")
  if (inherits(M, "listw")) {
    M <- sar_listw_matrix(M)
  } else if (inherits(M, "nb")) {
    M <- sar_nb_matrix(M)
  }
  if (length(dim(M)) != 2L || is.data.frame(M)) {
    frm_stop("'M' for sar() must be of class 'matrix', 'listw', or 'nb', ",
             "as in brms", call. = FALSE)
  }
  M <- methods::as(Matrix::Matrix(M, sparse = TRUE), "generalMatrix")
  if (nrow(M) != ncol(M) || nrow(M) != n_all) {
    frm_stop("Dimensions of 'M' for sar() must be equal to the number of ",
             "observations: M is ", nrow(M), " x ", ncol(M), " and the ",
             "data have ", n_all, " rows. As in brms, M covers every row ",
             "of the data, and the rows the model drops are dropped from ",
             "M", call. = FALSE)
  }
  M <- M[rows, rows, drop = FALSE]
  if (anyNA(M@x) || any(!is.finite(M@x))) {
    frm_stop("'M' for sar() has missing or infinite values", call. = FALSE)
  }
  dimnames(M) <- list(NULL, NULL)
  ev <- eigen(as.matrix(M), only.values = TRUE)$values
  if (is.complex(ev)) {
    # rounding leaves a few ulps of imaginary part on a W that is
    # similar to a symmetric matrix, such as a row-standardized
    # adjacency matrix; a real imaginary part makes log|I - rho W| a
    # sum over complex pairs, which brms does not compute either
    if (max(abs(Im(ev))) > 1e-8 * max(Mod(ev))) {
      frm_stop("sar(): the eigenvalues of 'M' are complex (largest ",
               "imaginary part ", format(max(abs(Im(ev))), digits = 3),
               "), so log|I - rho M| has no real form over them. A ",
               "spatial weights matrix from an undirected neighbor graph, ",
               "symmetric or row-standardized, has real eigenvalues",
               call. = FALSE)
    }
    ev <- Re(ev)
  }
  lo <- min(ev)
  hi <- max(ev)
  tol <- 1e-8 * max(abs(ev))
  if (!(lo < -tol && hi > tol)) {
    frm_stop("sar(): 'M' needs a negative and a positive eigenvalue, ",
             "which a weights matrix with a zero diagonal has; its range ",
             "is [", format(lo, digits = 4), ", ", format(hi, digits = 4),
             "]", call. = FALSE)
  }
  sar_warn_repeated(ev, tol)
  list(W = M, ev = ev, lo = 1 / lo, hi = 1 / hi, type = spec[["type"]])
}

#' The weights matrix of an spdep `nb` object, as brms takes it:
#' `spdep::nb2mat()` with its defaults, the row-standardized binary
#' matrix (style "W"), which refuses a location without neighbors. An
#' `nb` object is a list holding each location's neighbor indices, with
#' `0L` for none, so no spdep call is needed to read it.
#'
#' @noRd
sar_nb_matrix <- function(nb) {
  n <- length(nb)
  empty <- vapply(nb, function(v) length(v) == 1L && v[1L] == 0L, TRUE)
  if (any(empty)) {
    frm_stop("sar(): the 'nb' object has ", sum(empty), " location(s) ",
             "without neighbors, which spdep::nb2mat(), the conversion ",
             "brms uses, refuses (\"Empty neighbour sets found\")",
             call. = FALSE)
  }
  M <- matrix(0, n, n)
  for (i in seq_len(n)) M[i, nb[[i]]] <- 1 / length(nb[[i]])
  M
}

#' The weights matrix of an spdep `listw` object, as
#' `spdep::listw2mat()` gives it: each location's weights at its
#' neighbors' columns.
#'
#' @noRd
sar_listw_matrix <- function(lw) {
  nb <- lw[["neighbours"]]
  wt <- lw[["weights"]]
  n <- length(nb)
  if (!is.list(nb) || !is.list(wt) || length(wt) != n) {
    frm_stop("sar(): the 'listw' object needs neighbours and weights, ",
             "one entry per location", call. = FALSE)
  }
  M <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (length(nb[[i]]) == 1L && nb[[i]][1L] == 0L) next
    M[i, nb[[i]]] <- wt[[i]]
  }
  M
}

#' Warn when an extreme eigenvalue of `W` is repeated in all but one
#' direction, as on a complete graph. At that end of rho's interval
#' `I - rho W` keeps one direction, the constant one for a
#' row-standardized `W`, which an intercept spans: the residuals vanish
#' in every other direction and the likelihood of both types grows
#' without bound, so rho runs to the end. A warning and not a refusal:
#' a model without an intercept keeps its maximum, and `build` cannot
#' see the formula.
#'
#' @noRd
sar_warn_repeated <- function(ev, tol) {
  n <- length(ev)
  if (n < 3L) return(invisible(NULL))
  for (end in c(min(ev), max(ev))) {
    k <- sum(abs(ev - end) <= tol)
    if (k >= n - 1L) {
      frm_warning("sar(): the eigenvalue ", format(end, digits = 4),
                  " of 'M' is repeated ", k, " times in ", n, " rows ",
                  "(a complete graph, for one). At rho = ",
                  format(1 / end, digits = 4), ", the end of its ",
                  "interval, I - rho M keeps one direction, and with an ",
                  "intercept the likelihood has no maximum there: rho ",
                  "runs to that end. brms bounds rho by the eigenvalues ",
                  "themselves, which hides it. A weights matrix whose ",
                  "rows differ avoids it", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' rho from its unbounded parameter: a logistic map onto the interval
#' where `I - rho W` has a positive determinant.
#'
#' @noRd
sar_rho <- function(theta, aux) {
  aux[["lo"]] + (aux[["hi"]] - aux[["lo"]]) / (1 + exp(-theta[1]))
}

#' The innovations `(A y - mu) / sigma` (lag) or `A (y - mu) / sigma`
#' (error) and the log-Jacobian `log|A| - n log(sigma)`. AD-safe: `W` is
#' data, and a sparse `Matrix` times an advector stays on the tape.
#'
#' @noRd
sar_whiten <- function(y, mu, sigma, theta, aux) {
  rho <- sar_rho(theta, aux)
  W <- aux[["W"]]
  r <- if (identical(aux[["type"]], "lag")) {
    y - rho * as.vector(W %*% y) - mu
  } else {
    e <- y - mu
    e - rho * as.vector(W %*% e)
  }
  # one sigma: sigma_varies = FALSE
  list(u = r / sigma,
       logjac = sum(log1p(-rho * aux[["ev"]])) - length(y) * log(sigma))
}

#' The response at innovations `u`: `A^-1 (mu + sigma u)` (lag) or
#' `mu + A^-1 sigma u` (error). On numbers. `u` may be a matrix, one
#' column per set of innovations: one sparse factorization of `A` then
#' serves every column.
#'
#' @noRd
sar_unwhiten <- function(u, mu, sigma, theta, aux) {
  rho <- sar_rho(theta, aux)
  A <- Matrix::Diagonal(nrow(aux[["W"]])) - rho * aux[["W"]]
  out <- if (identical(aux[["type"]], "lag")) {
    as.matrix(Matrix::solve(A, mu + sigma * u))
  } else {
    mu + as.matrix(Matrix::solve(A, sigma * u))
  }
  if (is.matrix(u)) out else as.vector(out)
}

#' The interval of rho, for brms's class prior on it and frmtmb.sample's
#' flat default (`frmtmb::frmtmb_autocor_pars()`).
#'
#' @noRd
sar_bounds <- function(aux) {
  cbind(lower = aux[["lo"]], upper = aux[["hi"]])
}

#' rho under brms's name for the type.
#'
#' @noRd
sar_natural <- function(theta, aux) {
  stats::setNames(sar_rho(theta, aux), paste0(aux[["type"]], "sar"))
}

#' Start at rho = 0, the independent model.
#'
#' @noRd
sar_start <- function(aux) {
  stats::qlogis((0 - aux[["lo"]]) / (aux[["hi"]] - aux[["lo"]]))
}

# What frm_compat() says for sar() under gaussian() and student().
sar_compat_note <- "sar(M, type = \"lag\") is brms's SAR of the response, A y - mu ~ N(0, sigma^2 I) with A = I - rho W, whose mean is A^-1 mu; type = \"error\" is SAR of the residuals, A (y - mu) ~ N(0, sigma^2 I). student() gets the multivariate-t analog. The log-likelihood carries log|A| from the eigenvalues of W. rho is estimated on (1 / min eigenvalue, 1 / max eigenvalue), where |A| > 0; brms bounds it by (min eigenvalue, max eigenvalue). M is the weights matrix over every row of the data (a matrix, a sparse Matrix, or spdep's listw or nb), and the rows the model drops are dropped from it, as brms drops them. sigma and nu must be constant, as in brms."

#' Register `sar()` with core.
#'
#' @noRd
register_sar <- function() {
  frmtmb_register_autocor(
    "sar", args = c("M", "type"), parse = sar_parse, build = sar_build,
    whiten = sar_whiten, unwhiten = sar_unwhiten, npar = 1L,
    start = sar_start, natural = sar_natural, sigma_varies = FALSE,
    bounds = sar_bounds, usage = "sar(M, type = \"lag\")",
    compat_note = sar_compat_note)
}
