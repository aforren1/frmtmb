# Fix round 2 (review "Short check of fix round 1"), on the tree named
# as the argument.
# F1: cs_check_block() keeps the finiteness of lpdf and drops the level
#     comparison, which refused correct densities (an improper RW1, a
#     shifted mean) and would refuse a bimodal one under any curvature
#     test.
# F2: the names refused as formula functions are a fixed list of the
#     functions formulas call as terms, not every function of six
#     namespaces, and registration loads no package to decide.
root <- commandArgs(trailingOnly = TRUE)[1]
p <- file.path(root, "R/covstruct-seam.R")
s <- paste(readLines(p, warn = FALSE), collapse = "\n")
rep1 <- function(s, old, new) {
  if (!grepl(old, s, fixed = TRUE)) stop("no match: ", substr(old, 1, 70))
  sub(old, new, s, fixed = TRUE)
}

## F1 ----------------------------------------------------------------
s <- rep1(s, "#' registers, in addition: `lpdf` at the start is one finite number,
#' and larger at `b = 0` than far from it, which a log density centered
#' at zero is and a negative log density is not. The built-in
#' structures skip the `lpdf` evaluation, which for a dense field costs
#' a factorization, because their densities are tested directly.",
"#' registers, in addition: `lpdf` at `b = 0` and the start is one finite
#' number. Nothing is asked of its shape. A correct density can be flat
#' along a direction (an intrinsic RW1 or ICAR along the constant
#' vector), centered away from zero, or bimodal with a minimum at zero,
#' so no test of its level or curvature at zero separates it from a
#' density with the wrong sign; the author's own test of the structure
#' does. The built-in structures skip the `lpdf` evaluation, which for
#' a dense field costs a factorization, because their densities are
#' tested directly.")
s <- rep1(s, "  l1 <- e$nll(rep(1e3, nb), st, bk)
  if (is.numeric(l1) && length(l1) == 1L && !is.na(l1) && l1 >= l0) {
    bad(\"lpdf\", paste0(\"is no smaller at b = 1000 than at b = 0 (\",
                       format(l1, digits = 4), \" against \",
                       format(l0, digits = 4), \"); it has to be the LOG \",
                       \"DENSITY of b, not its negative\"))
  }
", "")

## F2 ----------------------------------------------------------------
s <- rep1(s, "#' registered in this session. `covstruct_formula_function()` adds the
#' functions of the packages a formula calls into.",
"#' registered in this session. `covstruct_formula_terms` adds the
#' functions that formulas call as terms.")
s <- rep1(s, "  unique(c(brms_whole_term_specials, \"mo\", \"mi\", \"me\", \"gr\", \"mm\",
           \"mmc\", \"offset\", \"I\", \"propto\", \"rr\", \"us_t\", \"diag_t\",
           \"smooth\", \"s\", \"t2\", \"te\", \"ti\", core_aterms,",
"  unique(c(brms_whole_term_specials, \"mo\", \"mi\", \"me\", \"gr\", \"mm\",
           \"mmc\", \"propto\", \"rr\", \"us_t\", \"diag_t\", \"smooth\",
           core_aterms,")
i0 <- regexpr("#' The package whose function `name` is", s, fixed = TRUE)
i1 <- regexpr("#' The package that owns a reserved name", s, fixed = TRUE)
stopifnot(i0 > 0, i1 > i0)
s <- paste0(substr(s, 1, i0 - 1), "#' The functions that formulas call as terms or around a term's
#' variable, which a structure of the same name would take from every
#' formula: `poly(x, 2)` stops working once `poly` is a structure. It
#' is a fixed list, and not every function of base, stats, splines or
#' mgcv, because a structure conflicts only with a call a formula
#' makes: `toeplitz`, `kernel` or `dist` is free. Being fixed, it also
#' loads no package to decide.
#'
#' @noRd
covstruct_formula_terms <- c(
  # splines, polynomials and mgcv's smooths
  \"poly\", \"polym\", \"ns\", \"bs\", \"s\", \"t2\", \"te\", \"ti\",
  # terms with a meaning of their own in a model formula
  \"offset\", \"I\", \"C\", \"cbind\", \"Surv\", \"strata\", \"cluster\",
  # transformations of a variable
  \"log\", \"log1p\", \"log2\", \"log10\", \"exp\", \"expm1\", \"sqrt\", \"abs\",
  \"sin\", \"cos\", \"tan\", \"plogis\", \"qlogis\", \"pnorm\", \"qnorm\",
  \"round\", \"floor\", \"ceiling\", \"pmin\", \"pmax\", \"ifelse\", \"scale\",
  # recodings of a variable
  \"factor\", \"ordered\", \"as.factor\", \"as.ordered\", \"as.numeric\",
  \"as.integer\", \"as.logical\", \"interaction\", \"cut\", \"relevel\"
)

", substr(s, i1, nchar(s)))
s <- rep1(s, "  fpkg <- covstruct_formula_function(name)
  if (!is.null(fpkg)) {
    frm_stop(\"frmtmb_register_covstruct() cannot register '\", name,
             \"', because \", name, \"() is a function of \", fpkg,
             \" that formulas call, and a structure of that name would \",
             \"take it from every formula. Choose another name.\",
             call. = FALSE)
  }",
"  if (name %in% covstruct_formula_terms) {
    frm_stop(\"frmtmb_register_covstruct() cannot register '\", name,
             \"', because formulas call \", name, \"() as a term (as in \",
             \"poly(x, 2) or log(x)), and a structure of that name would \",
             \"take it from every formula. Choose another name.\",
             call. = FALSE)
  }")
writeLines(s, p)

## NEWS --------------------------------------------------------------
p <- file.path(root, "NEWS.md")
s <- paste(readLines(p, warn = FALSE), collapse = "\n")
s <- rep1(s, "  registered structure that `lpdf` is a finite log density, larger at
  zero than far from it. A name the rest of frmtmb, brms or another
  extension owns is refused (`poly`, `ns`, `gaussian`, `dec`, ...), and",
"  registered structure that `lpdf` is finite at zero. A name the rest of
  frmtmb, brms or another extension owns, or a function formulas call
  as a term, is refused (`poly`, `ns`, `log`, `gaussian`, `dec`, ...), and")
writeLines(s, p)
cat("fix-r2 applied to", root, "\n")
