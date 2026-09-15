root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (n != 1L) stop(n, " matches in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

sub1("R/draws-generics.R",
paste0("#' @param x An object holding draws.\n",
       "#' @return A single integer.\n",
       "#' @examples\n#' dd <- data.frame(y = rnorm(40), x = rnorm(40))"),
paste0("#' @param x An object holding draws.\n",
       "#' @param ... Unused. Carried because posterior's `nvariables()`\n",
       "#'   generic has it and frmtmb hands these generics back to\n",
       "#'   posterior when posterior is loaded.\n",
       "#' @return A single integer.\n",
       "#' @examples\n#' dd <- data.frame(y = rnorm(40), x = rnorm(40))"))

# ---- sibling package, kept compilable against the new signatures ----
# These three are NOT a widening of scope: core's generics now match
# the owner's formals exactly (posterior's nvariables() takes (x, ...),
# brms's posterior_summary() takes (x, ...), nlme's VarCorr() takes
# (x, sigma, ...)), and an S3 method must carry every argument of its
# generic. Without these edits R CMD check on frmtmb.sample reports a
# generic/method mismatch that this lane created.
S <- "extensions/frmtmb.sample"
sub1(file.path(S, "R/methods-draws.R"),
  "nvariables.frmtmb_draws <- function(x) ncol(x$draws)",
  "nvariables.frmtmb_draws <- function(x, ...) ncol(x$draws)")
sub1(file.path(S, "R/methods-draws.R"),
paste0("posterior_summary.frmtmb_draws <- function(object, ",
       "probs = c(0.025, 0.975),"),
paste0("posterior_summary.frmtmb_draws <- function(x, ",
       "probs = c(0.025, 0.975),"))
sub1(file.path(S, "R/methods-draws.R"),
  "  posterior_summary(draws_columns(object, variable),",
  "  posterior_summary(draws_columns(x, variable),")
sub1(file.path(S, "R/methods-draws.R"),
  "VarCorr.frmtmb_draws <- function(x, ...) {",
  "VarCorr.frmtmb_draws <- function(x, sigma = 1, ...) {")
cat("DONE\n")
