#' Null-coalescing operator: `y` when `x` is `NULL`, else `x`.
#'
#' Defined here rather than imported so the package keeps working on the
#' oldest supported R, where base does not export it.
#'
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

# --- argument checks shared by the whole export surface ---------------
#
# These exist because the idioms they replace fail SILENTLY. `isTRUE(x)`
# reads a length-2 logical, a string, an integer and `NA` all as FALSE,
# so a flag set by mistake used to fit the OTHER model and report it as
# the user's. A length-2 `level` recycles against the parameter vector
# and puts DIFFERENT coverages on different rows of one table, with no
# column heading to record it. Each check turns one of those into a
# refusal that names the argument and says what it had to be.
#
# Every check raises exactly one message, so the whole family adds a
# handful of condition messages rather than one per call site, and the
# argument name that varies at run time comes from the caller.

#' What a rejected argument actually was, as a phrase a message can end
#' with. Deliberately short: the caller's message has already said what
#' the argument had to be, so this only has to identify what arrived.
#'
#' @noRd
arg_desc <- function(x) {
  if (is.null(x)) return("NULL")
  cls <- class(x)[1L]
  if (is.atomic(x) && length(x) == 1L) {
    if (is.na(x)) return(paste0(cls, " NA"))
    return(paste0(cls, " ", encodeString(format(x), quote = "\"")))
  }
  art <- if (substr(cls, 1L, 1L) %in% c("a", "e", "i", "o", "u")) "an" else "a"
  paste0(art, " ", cls, " of length ", length(x))
}

#' The one refusal every TRUE/FALSE argument in the package shares.
#'
#' @noRd
check_flag <- function(x, arg, what = NULL) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE, not ", arg_desc(x),
         if (is.null(what)) "" else paste0(". ", what), call. = FALSE)
  }
  invisible(x)
}

#' A single whole number at or above `min`: counts of draws, restarts
#' and grid points. A length-2 or fractional count otherwise reaches
#' `seq_len()` or `numeric()` as a size, where the message names neither
#' the argument nor the function the user called.
#'
#' @noRd
check_count <- function(x, arg, min = 0L) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x == round(x) && x >= min
  if (!ok) {
    stop("`", arg, "` must be a single whole number of at least ", min,
         ", not ", arg_desc(x), call. = FALSE)
  }
  invisible(as.integer(x))
}

#' A single number strictly inside (0, 1): every coverage argument
#' (`level`, `prob`).
#'
#' @noRd
check_probability <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x <= 0 || x >= 1) {
    stop("`", arg, "` must be a single number strictly between 0 and 1, ",
         "not ", arg_desc(x), call. = FALSE)
  }
  invisible(as.numeric(x))
}

#' A single finite positive number: scale parameters and tolerances.
#'
#' @noRd
check_positive <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be a single finite positive number, not ",
         arg_desc(x), call. = FALSE)
  }
  invisible(as.numeric(x))
}

#' A single finite number with no sign constraint: the location
#' arguments of the prior constructors.
#'
#' @noRd
check_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("`", arg, "` must be a single finite number, not ", arg_desc(x),
         call. = FALSE)
  }
  invisible(as.numeric(x))
}

#' A named list, for the argument slots keyed by parameter name. An
#' unnamed value used to be iterated over an empty `names()` and so
#' ignored in silence, and the fit then ran from the default and
#' reported it as the user's.
#'
#' @noRd
check_named_list <- function(x, arg, example) {
  # every element must carry a name: a PARTIALLY named list would have
  # its unnamed elements silently dropped by the by-name consumers
  if (!is.list(x) ||
        (length(x) && (is.null(names(x)) || !all(nzchar(names(x)))))) {
    stop("`", arg, "` must be a named list, e.g. ", example, ", not ",
         arg_desc(x), call. = FALSE)
  }
  invisible(x)
}

#' A single string drawn from a fixed set. `match.arg()` covers the
#' arguments whose default IS the set; this covers the ones whose
#' default is one value of it, where `match.arg()` has nothing to match
#' against and any string is accepted.
#'
#' @noRd
check_string_choice <- function(x, arg, choices) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
        !x %in% choices) {
    stop("`", arg, "` must be one of ",
         paste0("\"", choices, "\"", collapse = ", "), ", not ",
         arg_desc(x), call. = FALSE)
  }
  invisible(x)
}

#' Split `a + b + c` into list(a, b, c), left-associatively.
#'
#' @noRd
split_plus <- function(expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    c(split_plus(expr[[2]]), list(expr[[3]]))
  } else {
    list(expr)
  }
}

#' The `"response.dpar"` key that names one linear predictor throughout
#' the package. Written once so every table that is indexed by linear
#' predictor agrees on the separator.
#'
#' @noRd
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
#' @examples
#' # unequally spaced observation times, kept as distances
#' num_factor(c(0, 1.5, 4))
#'
#' # planar coordinates for a spatial covariance
#' levels(num_factor(rep(1:3, 3), rep(1:3, each = 3)))
#'
#' # ou() reads the distances out of the level labels; a plain factor
#' # would only give it an ordering
#' set.seed(1)
#' tim <- c(0, 1, 1.5, 3)
#' n_g <- 40
#' S <- 0.9^2 * exp(-1.2 * abs(outer(tim, tim, "-")))
#' u <- matrix(rnorm(n_g * length(tim)), n_g) %*% chol(S)
#' dd <- data.frame(
#'   y = 1 + as.vector(t(u)) + rnorm(n_g * length(tim), 0, 0.4),
#'   g = factor(rep(seq_len(n_g), each = length(tim))),
#'   tim = num_factor(rep(tim, n_g))
#' )
#' fit <- frm(bf(y ~ 1 + ou(tim + 0 | g)) + gaussian(), data = dd)
#' round(VarCorr(fit)[[1]], 3)
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

#' Recover coordinates from num_factor / glmmTMB::numFactor levels
#' (vector for 1-D, matrix for 2-D), or from plainly numeric labels.
#' Errors when a level does not encode a position.
#'
#' @noRd
parse_num_levels <- function(lv) {
  s <- gsub("[()]", "", lv)
  if (any(grepl(",", s, fixed = TRUE))) {
    parts <- strsplit(s, ",", fixed = TRUE)
    out <- t(vapply(parts, function(p) suppressWarnings(as.numeric(p)),
                    numeric(length(parts[[1]]))))
    if (anyNA(out)) {
      stop("Levels must encode numeric coordinates like '(1,2)'; build ",
           "the factor with num_factor(x, y)", call. = FALSE)
    }
    return(out)
  }
  out <- suppressWarnings(as.numeric(s))
  if (anyNA(out)) {
    stop("Levels must encode numeric positions on one axis; build the ",
         "factor with num_factor(x)", call. = FALSE)
  }
  out
}

#' The verdict both covariance paths give when the standard errors come
#' back non-finite. It is a separate verdict from "the Hessian is not
#' positive definite", because the two do not coincide: sdreport reports
#' pdHess from a Cholesky, which succeeds on a matrix LAPACK's solver
#' then refuses as computationally singular (reciprocal condition number
#' below its tolerance), and cov.fixed comes back filled with NaN on a
#' fit that converged with a small gradient and a positive definite
#' Hessian. diagnose() names the parameters; nothing else did.
#'
#' The verdict is a property of the FIT, not of the call, so it is worth
#' saying once per fit: vcov(), summary(), confint() and predict() all
#' reach it on the same degenerate object, and frm_multiple() and
#' influence() loop over those. The flag lives on the fit's cache
#' environment, which is replaced whenever the fit is re-estimated, so a
#' refit that is still degenerate warns again.
#'
#' @noRd
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

#' Invert the joint precision of a REML / profile fit.
#'
#' The ML branch of vcov() reads an already-inverted cov.fixed, which
#' sdreport fills with NaN when the Hessian is singular; inverting by
#' hand here would instead throw a raw LAPACK message that names neither
#' the model nor the remedy. Degrade the same way ML does: NaN entries
#' plus one warning pointing at diagnose().
#'
#' @noRd
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

#' One string key per coordinate row, used to match gp() prediction
#' positions against fitted positions. Defined once so frame assembly
#' and kriging can never disagree on the separator.
#'
#' @noRd
pos_rowkey <- function(M) {
  do.call(paste, c(as.data.frame(M), sep = "\r"))
}
