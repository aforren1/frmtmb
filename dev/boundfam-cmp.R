# Lane boundfam: compare the standard errors of every analysis between
# two traced suite runs (dev/boundfam-par.sh base / lane).
#   Rscript dev/boundfam-cmp.R <dir a> <dir b>
# One line per test file that both arms traced, then the totals: COV
# records equal bitwise, equal to 1e-8 relative (cross-process ulps,
# dev/rtmb-pitfalls.md item 21), moved, and NaN pattern changed; and
# the same for VCT (covariance-interval) records.
a <- commandArgs(trailingOnly = TRUE)
da <- a[1]
db <- a[2]
rd <- function(f, tag) {
  if (!file.exists(f)) return(NULL)
  x <- readLines(f, warn = FALSE)
  x <- x[startsWith(x, paste0(tag, " "))]
  key <- if (tag == "COV") "se=" else "se_t="
  lapply(x, function(l) {
    v <- sub(paste0(".*", key), "", l)
    v <- suppressWarnings(as.numeric(strsplit(v, ",", fixed = TRUE)[[1]]))
    list(v = v, lost = if (tag == "COV") sub(" se=.*", "", sub(".*lost=", "", l)))
  })
}
cmp1 <- function(x, y) {
  if (length(x$v) != length(y$v)) return("length")
  na_x <- is.na(x$v)
  na_y <- is.na(y$v)
  if (!identical(na_x, na_y)) return("nan")
  if (identical(x$v, y$v)) return("bitwise")
  ok <- !na_x
  r <- max(abs(x$v[ok] - y$v[ok]) / pmax(abs(x$v[ok]), 1e-300))
  if (r <= 1e-8) "ulps" else "moved"
}
files <- intersect(basename(Sys.glob(file.path(da, "*.trace"))),
                   basename(Sys.glob(file.path(db, "*.trace"))))
tot <- list(COV = integer(0), VCT = integer(0))
for (f in files) {
  for (tag in c("COV", "VCT")) {
    xa <- rd(file.path(da, f), tag)
    xb <- rd(file.path(db, f), tag)
    n <- min(length(xa), length(xb))
    res <- character(0)
    if (n) res <- vapply(seq_len(n), function(i) cmp1(xa[[i]], xb[[i]]), "")
    if (length(xa) != length(xb)) res <- c(res, "count")
    tot[[tag]] <- c(tot[[tag]], res)
    bad <- which(!res %in% c("bitwise", "ulps"))
    if (length(bad)) {
      cat(sprintf("%s %s: %d of %d records differ (%s)\n", f, tag,
                  length(bad), length(res),
                  paste(unique(res[bad]), collapse = ",")))
      for (i in bad[bad <= n][seq_len(min(4, sum(bad <= n)))]) {
        cat("   #", i, " a:", signif(xa[[i]]$v, 4), " lost", xa[[i]]$lost, "\n")
        cat("   #", i, " b:", signif(xb[[i]]$v, 4), " lost", xb[[i]]$lost, "\n")
      }
    }
  }
}
cat("files compared", length(files), "\n")
for (tag in names(tot)) {
  t <- table(factor(tot[[tag]],
                    levels = c("bitwise", "ulps", "moved", "nan", "length",
                               "count")))
  cat(tag, paste(names(t), t, sep = "=", collapse = " "), "\n")
}
