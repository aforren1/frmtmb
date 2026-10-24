# Fix round 2: the exported docs of frmtmb_register_covstruct() for F1
# and F2, on the tree named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
p <- file.path(root, "R/covstruct-seam.R")
s <- paste(readLines(p, warn = FALSE), collapse = "\n")
rep1 <- function(s, old, new) {
  if (!grepl(old, s, fixed = TRUE)) stop("no match: ", substr(old, 1, 70))
  sub(old, new, s, fixed = TRUE)
}
s <- rep1(s, "#'   as another extension's (only that package registers it), and not a
#'   function of base, stats, splines, utils, methods or mgcv.",
"#'   as another extension's (only that package registers it), and not a
#'   function that formulas call as a term, such as `poly`, `ns`, `log`
#'   or `factor`.")
s <- rep1(s, "#' @param lpdf `function(b, theta, blk)`: the log density of the block's
#'   `b`, summed over levels. AD-safe.",
"#' @param lpdf `function(b, theta, blk)`: the log density of the block's
#'   `b`, summed over levels. AD-safe. The frame checks only that it is
#'   finite at `b = 0` and the start: a correct density can be flat along
#'   a direction, centered away from zero or bimodal, so its shape tells
#'   nothing, and a density with the wrong sign is for the structure's
#'   own test to catch (\"Testing a structure\").")
writeLines(s, p)
cat("fix-r2b applied to", root, "\n")
