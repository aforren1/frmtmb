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

# One string key per coordinate row, used to match gp() prediction
# positions against fitted positions. Defined once so frame assembly
# and kriging can never disagree on the separator.
pos_rowkey <- function(M) {
  do.call(paste, c(as.data.frame(M), sep = "\r"))
}
