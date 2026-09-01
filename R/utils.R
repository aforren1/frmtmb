`%||%` <- function(x, y) if (is.null(x)) y else x

# Split `a + b + c` into list(a, b, c), left-associatively.
split_plus <- function(expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    c(split_plus(expr[[2]]), list(expr[[3]]))
  } else {
    list(expr)
  }
}

linpred_key <- function(resp, dpar) paste(resp, dpar, sep = ".")

#' Factor with numeric-coded levels for coordinate covariance structures
#'
#' `ou()` and the spatial structures (`exp()`, `gau()`, `mat()`) need
#' the positions of the term levels. `num_factor(x)` (one dimension) or
#' `num_factor(x, y)` (planar coordinates) encodes them in the level
#' labels the same way `glmmTMB::numFactor()` does, so factors created
#' by either function work.
#'
#' @param x Numeric positions (times, coordinates).
#' @param y Optional second coordinate.
#' @return A factor whose levels encode the unique positions.
#' @export
num_factor <- function(x, y = NULL) {
  if (is.null(y)) {
    ux <- sort(unique(as.numeric(x)))
    return(factor(sprintf("(%g)", as.numeric(x)),
                  levels = sprintf("(%g)", ux)))
  }
  lab <- sprintf("(%g,%g)", as.numeric(x), as.numeric(y))
  ord <- order(as.numeric(x), as.numeric(y))
  factor(lab, levels = unique(lab[ord]))
}

# Recover coordinates from num_factor / glmmTMB::numFactor levels
# (vector for 1-D, matrix for 2-D), or from plainly numeric labels.
parse_num_levels <- function(lv) {
  s <- gsub("[()]", "", lv)
  if (any(grepl(",", s, fixed = TRUE))) {
    parts <- strsplit(s, ",", fixed = TRUE)
    out <- t(vapply(parts, function(p) suppressWarnings(as.numeric(p)),
                    numeric(length(parts[[1]]))))
    if (anyNA(out)) {
      stop("Levels must encode numeric positions; build the factor ",
           "with num_factor()", call. = FALSE)
    }
    return(out)
  }
  out <- suppressWarnings(as.numeric(s))
  if (anyNA(out)) {
    stop("Levels must encode numeric positions; build the factor with ",
         "num_factor()", call. = FALSE)
  }
  out
}

# The verdict both covariance paths give when the standard errors come
# back non-finite. It is a separate verdict from "the Hessian is not
# positive definite", because the two do not coincide: sdreport reports
# pdHess from a Cholesky, which succeeds on a matrix LAPACK's solver
# then refuses as computationally singular (reciprocal condition number
# below its tolerance), and cov.fixed comes back filled with NaN on a
# fit that converged with a small gradient and a positive definite
# Hessian. diagnose() names the parameters; nothing else did.
#
# The verdict is a property of the FIT, not of the call, so it is worth
# saying once per fit: vcov(), summary(), confint() and predict() all
# reach it on the same degenerate object, and frm_multiple() and
# influence() loop over those. The flag lives on the fit's cache
# environment, which is replaced whenever the fit is re-estimated, so a
# refit that is still degenerate warns again.
warn_nonfinite_cov <- function(cache = NULL) {
  if (is.environment(cache)) {
    if (isTRUE(cache$warned_nonfinite_cov)) return(invisible(NULL))
    cache$warned_nonfinite_cov <- TRUE
  }
  warning("Some standard errors are not finite, so vcov() and ",
          "summary() report NaN: the covariance could not be recovered ",
          "from the Hessian and the model is probably ",
          "overparameterized. diagnose() names the offending ",
          "parameters; see the 'Convergence problems' section of ",
          "vignette('diagnostics')", call. = FALSE)
}

# Invert the joint precision of a REML / profile fit.
#
# The ML branch of vcov() reads an already-inverted cov.fixed, which
# sdreport fills with NaN when the Hessian is singular; inverting by
# hand here would instead throw a raw LAPACK message that names neither
# the model nor the remedy. Degrade the same way ML does: NaN entries
# plus one warning pointing at diagnose().
solve_joint_precision <- function(Q, cache = NULL) {
  # Matrix::solve, not base solve: the joint precision of a GLMM is
  # sparse and often badly conditioned, and base's dense LAPACK path
  # refuses it on a reciprocal-condition-number test that the sparse
  # Cholesky has no reason to apply. Rejecting an invertible precision
  # would turn every standard error on such a fit into NaN.
  # CHOLMOD narrates a failed factorization from C before the error
  # reaches R; the message below is the one that names the remedy
  V <- try(suppressWarnings(Matrix::solve(Q)), silent = TRUE)
  if (!inherits(V, "try-error")) {
    # A precision matrix can factor and still invert into non-finite
    # entries. That is the same verdict as a failed solve() from the
    # user's side - NaN standard errors - so it earns the same warning,
    # or a profile fit reports NaN as quietly as the ML branch used to.
    # The x slot is the stored values of a Matrix; a base matrix has no
    # such slot and is read whole.
    xs <- tryCatch(V@x, error = function(e) as.numeric(as.matrix(V)))
    if (any(!is.finite(xs))) warn_nonfinite_cov(cache)
    return(V)
  }
  if (is.environment(cache)) {
    # same once-per-fit rule as warn_nonfinite_cov(): a singular joint
    # precision stays singular for every caller that asks
    if (isTRUE(cache$warned_singular_precision)) {
      return(matrix(NaN, nrow(Q), ncol(Q), dimnames = dimnames(Q)))
    }
    cache$warned_singular_precision <- TRUE
  }
  warning("The joint precision matrix is singular, so standard errors ",
          "are NaN; the model is probably overparameterized. ",
          "diagnose() names the offending parameter; see the ",
          "'Convergence problems' section of vignette('diagnostics')",
          call. = FALSE)
  matrix(NaN, nrow(Q), ncol(Q), dimnames = dimnames(Q))
}

# One string key per coordinate row, used to match gp() prediction
# positions against fitted positions. Defined once so frame assembly
# and kriging can never disagree on the separator.
pos_rowkey <- function(M) {
  do.call(paste, c(as.data.frame(M), sep = "\r"))
}
