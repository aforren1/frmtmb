# Phase 2.5a: drop frmtmb's rival generics for the four names that are
# taken by importFrom() instead, and move each documentation block onto
# the method it now documents.  Done as a script so the exact before
# and after text is on the record.
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"

sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (n != 1L) stop(n, " matches in ", path)
  txt <- sub(old, new, txt, fixed = TRUE)
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

sub1("R/methods-fit.R",
"#' @export\nfixef <- function(object, ...) UseMethod(\"fixef\")\n\n#' @rdname fixef\n#' @exportS3Method nlme::fixef\n#' @export\n",
"#' @rdname fixef\n#' @aliases fixef\n#' @export\n")

sub1("R/methods-fit.R",
"#' @export\nranef <- function(object, ...) UseMethod(\"ranef\")\n\n#' @rdname ranef\n#' @exportS3Method nlme::ranef\n#' @export\n",
"#' @rdname ranef\n#' @aliases ranef\n#' @export\n")

sub1("R/methods-fit.R",
"#' @export\nVarCorr <- function(x, ...) UseMethod(\"VarCorr\")\n\n#' @rdname VarCorr\n#' @exportS3Method nlme::VarCorr\n#' @export\nVarCorr.frmtmb_fit <- function(x, ...) {",
"#' @rdname VarCorr\n#' @aliases VarCorr\n#' @export\nVarCorr.frmtmb_fit <- function(x, sigma = 1, ...) {")

sub1("R/sugar.R",
"#' @export\nrefit <- function(object, newresp, ...) UseMethod(\"refit\")\n",
"#' @name refit\n#' @aliases refit\nNULL\n")

# nvariables: posterior's generic takes (x, ...) and frmtmb's took (x).
# The comment above it asserted the opposite; dev/generics-audit2.R
# measured it.
sub1("R/draws-generics.R",
"# posterior's nchains()/ndraws()/niterations()/nvariables() generics take\n# x alone, so these do too",
"# posterior's nchains(), ndraws() and niterations() generics take x\n# alone and nvariables() takes (x, ...); dev/generics-audit2.R measured\n# it.  The signatures have to agree because frmtmb hands these generics\n# back to posterior when posterior is loaded.")

sub1("R/draws-generics.R",
"nvariables <- function(x) UseMethod(\"nvariables\")",
"nvariables <- function(x, ...) UseMethod(\"nvariables\")")

# posterior_summary: brms's generic takes (x, ...), frmtmb's took
# (object, ...).  Same reason.
sub1("R/draws-generics.R",
"#' @param object A matrix of draws, variables in columns.",
"#' @param x A matrix of draws, variables in columns.")

sub1("R/draws-generics.R",
"posterior_summary <- function(object, ...) UseMethod(\"posterior_summary\")",
"posterior_summary <- function(x, ...) UseMethod(\"posterior_summary\")")

sub1("R/draws-generics.R",
"posterior_summary.default <- function(object, probs = c(0.025, 0.975),\n                                      robust = FALSE, ...) {\n  m <- as.matrix(object)",
"posterior_summary.default <- function(x, probs = c(0.025, 0.975),\n                                      robust = FALSE, ...) {\n  m <- as.matrix(x)")

cat("DONE\n")
