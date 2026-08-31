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
