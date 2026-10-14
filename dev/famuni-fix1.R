p <- "extensions/frmtmb.eam/R/wiener-family.R"
x <- paste(readLines(p), collapse = "\n")
old <- "#' not a density of anything. `0 / (x + |x|)` is exactly 0 for a
#' positive `x` and `0 / 0`, NaN, otherwise, which the optimizer reads
#' as a rejected step, as frmtmb's own families do; unlike `0 * log(x)`
#' it raises no warning when a post-fit method evaluates it on doubles."
new <- "#' not a density of anything. With `u = x + |x|`, `u / u - 1` is exactly
#' 0 for a positive `x`, with derivatives 0, and `0 / 0`, NaN,
#' otherwise, which the optimizer reads as a rejected step, as frmtmb's
#' own families do. It raises no warning on doubles, where `log(x)` would.
#'
#' The guard must not have a constant 0 as an operand. RTMB folds `0 * x`
#' and `0 / x` to the constant 0 when it records the tape, so the first
#' spelling, `0 / (x + |x|)`, was NaN on doubles and 0 in every fit: at a
#' bias of -0.5 the objective was finite, with a log-likelihood above the
#' maximum (dev/reviews/2026-10-13-famuni.md, fix 1;
#' dev/famuni-out/fold-probe.txt)."
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
old2 <- "  pos <- function(x) 0 / (x + abs(x))"
new2 <- "  pos <- function(x) {\n    u <- x + abs(x)\n    u / u - 1\n  }"
stopifnot(grepl(old2, x, fixed = TRUE))
x <- sub(old2, new2, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
