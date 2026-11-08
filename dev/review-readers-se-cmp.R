# Review of lane readers: compare se-before.rds and se-after.rds from
# dev/review-readers-se.R. The two arms ran in separate processes, so a
# difference of a few ulps is the process (dev/rtmb-pitfalls.md item 21).
a <- readRDS("dev/review-readers-out/se-before.rds")
b <- readRDS("dev/review-readers-out/se-after.rds")
for (nm in names(b)) {
  x <- a[[nm]]$ce; y <- b[[nm]]$ce
  if (is.character(x) || is.character(y)) {
    cat(sprintf("%-16s before: %s | after: %s\n", nm,
                if (is.character(x)) x else "ok",
                if (is.character(y)) y else "ok"))
    next
  }
  rd <- function(col) {
    if (!col %in% names(x)) return(NA)
    max(abs(y[[col]] / x[[col]] - 1), na.rm = TRUE)
  }
  ratio <- y$se__ / x$se__
  fe <- b[[nm]]$fitted
  vs_fit <- if (!is.null(fe) && length(dim(fe)) == 2L &&
                  nrow(fe) == nrow(y)) {
    format(max(abs(y$se__ / fe[, "Est.Error"] - 1)), digits = 3)
  } else "-"
  cat(sprintf(paste("%-16s est %.1e lower %.1e upper %.1e | se__ after /",
                    "before %s | se__ / fitted Est.Error - 1: %s\n"),
              nm, rd("estimate__"), rd("lower__"), rd("upper__"),
              paste(format(range(ratio), digits = 4), collapse = " to "),
              vs_fit))
}
