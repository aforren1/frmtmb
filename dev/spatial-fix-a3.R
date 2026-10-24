# Fix round 1: cs_entry()'s refusal goes through require_extension()
# (A3), the method heads check the owners (A3), the stale comment (A8),
# on the tree named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
rep_file <- function(p, old, new) {
  p <- file.path(root, p)
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, s, fixed = TRUE)[[1]])
  if (regexpr(old, s, fixed = TRUE) < 0 || n != 1) {
    stop(p, ": count ", n, ": ", substr(old, 1, 70))
  }
  writeLines(sub(old, new, s, fixed = TRUE), p)
}
rep_file("R/covstruct.R", '  home <- extension_special_home(name)
  if (!is.null(home)) {
    frm_stop("The ", name, "() structure is in the extension package ",
             home, ", which is not loaded. Install it and run library(",
             home, ")", call. = FALSE)
  }
  frm_stop("internal: no covariance structure \'", name, "\' is registered",
           call. = FALSE)', '  # loading the extension that owns the name registers it
  home <- extension_special_home(name)
  if (!is.null(home)) {
    require_extension(home, paste0("The ", name, "() structure"))
    e <- frmtmb_covstruct_registry$reg[[name]]
    if (!is.null(e)) return(e)
  }
  frm_stop("The ", name, "() structure is not registered in this session. ",
           "A fit that uses it needs the package that registered it ",
           "loaded, or the registration repeated",
           class = "frmtmb_extension_missing", call. = FALSE)')
rep_file("R/covstruct.R", "# The three readers below exist so that `covstruct_registry` and
# `lkj_refusals` can stay internal", "# The three readers below exist so that the registry and
# `lkj_refusals` can stay internal")
# the methods that read a block's entry without going through sdr_of()
# or the objective check the owners first
heads <- list(
  c("R/methods-fit.R", "print.frmtmb_fit <- function(x, ...) {\n", "x"),
  c("R/methods-fit.R", "ranef.frmtmb_fit <- function(object, summary = TRUE, robust = FALSE,", "object"),
  c("R/methods-fit.R", "VarCorr.frmtmb_fit <- function(x, sigma = 1, summary = TRUE,", "x"),
  c("R/predict-brms.R", "predict.frmtmb_fit <- function(object, newdata = NULL, re_formula = NULL,", "object"),
  c("R/predict.R", "fitted.frmtmb_fit <- function(object, newdata = NULL, re_formula = NULL,", "object"),
  c("R/predict.R", "frm_linpred <- function(object, newdata = NULL,", "object"),
  c("R/predict.R", "simulate.frmtmb_fit <- function(object, nsim = 1, seed = NULL,", "object"),
  c("R/confint.R", "update.frmtmb_fit <- function(object, formula., ..., evaluate = TRUE) {\n", "object"),
  c("R/confint.R", "confint_varcorr <- function(fit, level = 0.95) {\n", "fit"))
for (h in heads) {
  p <- file.path(root, h[1])
  x <- readLines(p, warn = FALSE)
  i <- grep(sub("\n$", "", h[2]), x, fixed = TRUE)
  stopifnot(length(i) == 1L)
  # the line that closes the signature
  j <- i
  while (!grepl("\\{\\s*$", x[j])) j <- j + 1L
  if (grepl("require_structure_packages", x[j + 1L], fixed = TRUE)) next
  x <- append(x, paste0("  require_structure_packages(", h[3], ")"), j)
  writeLines(x, p)
}
cat("fix-a3 applied to", root, "\n")
