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
#' `ou()` (and future spatial structures) need the positions of the term
#' levels. `num_factor(x)` encodes them in the level labels the same way
#' `glmmTMB::numFactor()` does, so factors created by either function
#' work.
#'
#' @param x Numeric positions (times, coordinates).
#' @return A factor whose levels encode the sorted unique positions.
#' @export
num_factor <- function(x) {
  ux <- sort(unique(as.numeric(x)))
  factor(sprintf("(%g)", as.numeric(x)),
         levels = sprintf("(%g)", ux))
}

# Recover coordinates from num_factor / glmmTMB::numFactor levels, or
# from plainly numeric level labels.
parse_num_levels <- function(lv) {
  out <- suppressWarnings(as.numeric(gsub("[()]", "", lv)))
  if (anyNA(out)) {
    stop("Levels must encode numeric positions; build the factor with ",
         "num_factor()", call. = FALSE)
  }
  out
}
