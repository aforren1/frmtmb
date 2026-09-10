# identical() over the two smoke records, quantity by quantity, plus a
# ulp count where the two differ. A printed zero is not a measured zero,
# so nothing here is compared with a tolerance.
#
#   Rscript dev/rlddm-scripts/rlddm-compare.R <ref.rds> <new.rds>

args <- commandArgs(trailingOnly = TRUE)
a <- readRDS(args[[1L]])
b <- readRDS(args[[2L]])

ulps <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  d <- abs(x - y)
  e <- pmax(abs(x), abs(y))
  out <- ifelse(d == 0, 0, d / (e * .Machine$double.eps))
  max(out[is.finite(out)], 0)
}

nm <- union(names(a), names(b))
cat(sprintf("%-20s %-10s %s\n", "quantity", "identical", "max ulp"))
for (k in nm) {
  x <- a[[k]]
  y <- b[[k]]
  id <- identical(x, y)
  u <- if (is.numeric(x) && is.numeric(y) && length(x) == length(y)) {
    formatC(ulps(x, y), digits = 4, format = "g")
  } else {
    ""
  }
  cat(sprintf("%-20s %-10s %s\n", k, if (id) "YES" else "no", u))
  if (!id && is.character(x) != is.character(y)) {
    cat("   ref: ", substr(paste(format(x), collapse = " "), 1, 100),
        "\n   new: ", substr(paste(format(y), collapse = " "), 1, 100),
        "\n", sep = "")
  }
}
