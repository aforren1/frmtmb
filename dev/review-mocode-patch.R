## Reviewer prototype of the fix proposed for findings on fractional
## mo() new data and conditional_effects()' reference values. Patches a
## COPY of the lane's R/ (path in the first argument), never the
## worktree.
##   Rscript dev/review-mocode-patch.R <copy root>
root <- commandArgs(trailingOnly = TRUE)[1L]
sub1 <- function(file, old, new) {
  f <- file.path(root, "R", file)
  x <- paste(readLines(f), collapse = "\n")
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <- sub(old, new, x, fixed = TRUE)
  writeLines(x, f)
}
sub1("predict.R",
"  lo <- mo_min(mi)
  codes <- as.integer(round(v - lo))   # grids may land between categories
  if (anyNA(codes) || any(codes < 0L | codes > mi$D)) {",
"  lo <- mo_min(mi)
  if (!is.numeric(v) || anyNA(v) || any(v != round(v)) ||
      any(v < lo | v > lo + mi$D)) {")
sub1("predict.R",
"             \"range \", lo, \"..\", lo + mi$D, call. = FALSE)
  }
  codes
}",
"             \"range \", lo, \"..\", lo + mi$D, call. = FALSE)
  }
  as.integer(v - lo)
}")
sub1("conditional-effects.R",
"ce_ref_value <- function(col) {",
"ce_ref_value <- function(col, step = FALSE) {
  if (step && is.numeric(col) && !is.matrix(col)) {
    return(round(stats::median(col, na.rm = TRUE)))
  }")
sub1("conditional-effects.R",
"ce_second_values <- function(col, int_cond = NULL) {
  if (!is.null(int_cond)) return(ce_int_cond(int_cond, col))",
"ce_second_values <- function(col, int_cond = NULL, step = FALSE) {
  if (!is.null(int_cond)) return(ce_int_cond(int_cond, col))
  if (step && is.numeric(col) && !is.matrix(col)) {
    return(round(stats::median(col, na.rm = TRUE) +
                   c(-1, 0, 1) * stats::mad(col, na.rm = TRUE)))
  }")
sub1("conditional-effects.R",
"                        extend_vars = character(0)) {
  nd <- data.frame(.ce_row = seq_len(n))",
"                        extend_vars = character(0),
                        step_vars = character(0)) {
  nd <- data.frame(.ce_row = seq_len(n))")
sub1("conditional-effects.R",
"      ce_ref_value(base[[nm]])",
"      ce_ref_value(base[[nm]], step = nm %in% step_vars)")
sub1("conditional-effects.R",
"      ce_second_values(base[[ev[2L]]], int_conditions[[ev[2L]]])",
"      ce_second_values(base[[ev[2L]]], int_conditions[[ev[2L]]],
                       step = ev[2L] %in% step_vars)")
sub1("conditional-effects.R",
"                        na_vars, setdiff(na_vars, ce_locked_vars(x)))",
"                        na_vars, setdiff(na_vars, ce_locked_vars(x)),
                        step_vars)")
cat("patched\n")
