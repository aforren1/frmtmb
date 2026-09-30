# Lane sampfix last round (R1, R3): the watch compares NA patterns and
# builds no names.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/R/methods-draws.R"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("found ", n, ":\n", old)
  s <<- sub(old, new, s, fixed = TRUE)
}
rep1("#' draws 2 to 6 were NaN. So on laplace draws each draw's non-finite
#' cells must be the first draw's; a draw with others read an integrated
#' value, and the call is refused rather than returning NaN. The first
#' draw's pattern is the reference, not \"none\", so a row that is `NA` at
#' every draw (a missing covariate) is not taken for one.",
"#' draws 2 to 6 were NaN. So on laplace draws each draw's `NA` cells
#' must be the first draw's; a draw with others read an integrated
#' value, and the call is refused rather than returning NaN. The first
#' draw's pattern is the reference, not \"none\", so a row that is `NA` at
#' every draw (a missing covariate) is not taken for one.
#'
#' `is.na()`, not `!is.finite()`: a read of the `NA` fill came out `NA`
#' or NaN on every watched path and never `Inf` (38 path and model
#' combinations, dev/sampfix-12-fillpaths.R), while a later draw that
#' overflows to `Inf` on its own was refused as a read. The names are
#' dropped because `unlist()` building one per cell was nearly all of
#' the watch's cost, about as much as the prediction itself.")
rep1("    nf <- !is.finite(suppressWarnings(as.numeric(unlist(v))))",
     "    nf <- is.na(unlist(v, use.names = FALSE))")
writeLines(s, f)
cat("ok\n")
