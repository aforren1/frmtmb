# One-off edit of R/bf.R: bf(autocor = ). Kept as the record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/bf.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(n, " matches for ", old)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("bf <- function(formula, ..., family = NULL, nl = NULL, center = NULL,
               cmc = NULL) {",
"bf <- function(formula, ..., family = NULL, nl = NULL, center = NULL,
               cmc = NULL, autocor = NULL) {")
rep1("    return(bf_update(formula, ..., family = family, nl = nl,
                     center = center, cmc = cmc))",
"    return(bf_update(formula, ..., family = family, nl = nl,
                     center = center, cmc = cmc, autocor = autocor))")
rep1("      bf(f1, ..., family = family, nl = nl, center = center, cmc = cmc)",
"      bf(f1, ..., family = family, nl = nl, center = center, cmc = cmc,
         autocor = autocor)")
rep1("  parsed <- bf_dots(list(...))
  pforms <- parsed[[\"pforms\"]]",
"  # brms's autocor argument: the terms join the formula of mu
  if (!is.null(autocor)) formula <- add_ac_terms(formula, autocor)
  parsed <- bf_dots(list(...))
  pforms <- parsed[[\"pforms\"]]")
rep1("bf_update <- function(formula, ..., family = NULL, nl = NULL,
                      center = NULL, cmc = NULL) {
  dots <- list(...)
  if (!length(dots) && is.null(family) && is.null(nl) && is.null(center) &&
        is.null(cmc)) {
    return(formula)
  }",
"bf_update <- function(formula, ..., family = NULL, nl = NULL,
                      center = NULL, cmc = NULL, autocor = NULL) {
  dots <- list(...)
  if (!length(dots) && is.null(family) && is.null(nl) && is.null(center) &&
        is.null(cmc) && is.null(autocor)) {
    return(formula)
  }
  if (!is.null(autocor)) {
    formula[[\"formula\"]] <- add_ac_terms(formula[[\"formula\"]], autocor)
  }")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
