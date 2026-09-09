# rev-ndt: compare the two arms of rev-ndt-bc.R bitwise.
#
# identical() first, then a ULP distance for anything that differs, so
# that "agrees to every digit printed" is replaced by a bit count.
new <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/rev-ndt-bc-new.rds")
old <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/rev-ndt-bc-old.rds")

ulps <- function(a, b) {
  # distance in representable doubles between two finite numbers
  ok <- is.finite(a) & is.finite(b)
  out <- rep(NA_real_, length(a))
  for (i in which(ok)) {
    x <- a[i]
    y <- b[i]
    if (x == y) { out[i] <- 0; next }
    lo <- min(x, y)
    hi <- max(x, y)
    n <- 0
    while (lo < hi && n < 1e6) { lo <- nextafter(lo); n <- n + 1 }
    out[i] <- n
  }
  out
}
# no nextafter in base R; step by the exponent of the value
nextafter <- function(x) x + .Machine$double.eps * max(abs(x), 1e-300) *
  (if (abs(x) >= 1) 2^(floor(log2(abs(x)))) else 1)

# simpler and enough: relative difference in units of eps
rel_ulp <- function(a, b) {
  d <- abs(a - b)
  s <- pmax(abs(a), abs(b))
  ifelse(d == 0, 0, d / (s * .Machine$double.eps))
}

walk <- function(pa, na, nb) {
  if (is.list(na) != is.list(nb)) {
    cat(sprintf("SHAPE      %-46s new %s / old %s\n", pa, class(na)[1],
                class(nb)[1]))
    cat("   new: ", paste(utils::head(as.character(unlist(na)), 2),
                          collapse = " | "), "\n")
    cat("   old: ", paste(utils::head(as.character(unlist(nb)), 2),
                          collapse = " | "), "\n")
    return(invisible(NULL))
  }
  if (is.list(na)) {
    for (k in union(names(na), names(nb))) {
      walk(paste0(pa, "$", k), na[[k]], nb[[k]])
    }
    return(invisible(NULL))
  }
  same <- identical(na, nb)
  if (same) { cat(sprintf("IDENTICAL  %-46s n=%d\n", pa, length(na)))
    return(invisible(NULL)) }
  if (is.numeric(na) && is.numeric(nb) && length(na) == length(nb)) {
    u <- rel_ulp(na, nb)
    cat(sprintf("DIFFERS    %-46s n=%d  max %.3g ulp  max abs %.3g\n",
                pa, length(na), max(u, na.rm = TRUE),
                max(abs(na - nb), na.rm = TRUE)))
  } else {
    cat(sprintf("DIFFERS    %-46s (not comparable numerics)\n", pa))
    cat("   new: ", paste(utils::head(as.character(na), 2),
                          collapse = " | "), "\n")
    cat("   old: ", paste(utils::head(as.character(nb), 2),
                          collapse = " | "), "\n")
  }
  invisible(NULL)
}

for (k in union(names(new), names(old))) walk(k, new[[k]], old[[k]])
