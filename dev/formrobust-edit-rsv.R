# One-off edit of R/parse.R: the lower-case intercept. Kept as the record.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/parse.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("  list(fixed = new, pos = k - 1L)\n}", "  list(fixed = new, pos = k - 1L, lower = lower)\n}")
rep1("#' than guessed. So is brms's deprecated lower-case `intercept`.",
"#' than guessed. brms's deprecated lower-case `intercept` is read as\n#' brms reads it, a data column of ones with a warning; `lower` says so,\n#' and is all the result holds when `Intercept` itself is absent.")
rep1("  rsv <- rsv_intercept_fixed(fixed)
  if (!is.null(rsv)) {",
"  rsv <- rsv_intercept_fixed(fixed)
  rsv_lower <- isTRUE(rsv[[\"lower\"]])
  if (!is.null(rsv[[\"fixed\"]])) {")
rep1("  if (!is.null(rsv)) {
    out$center <- FALSE
    out$rsv_intercept <- rsv$pos
  }
  out
}",
"  if (!is.null(rsv[[\"fixed\"]])) {
    out$center <- FALSE
    out$rsv_intercept <- rsv$pos
  }
  if (rsv_lower) out$rsv_lower <- TRUE
  out
}")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
