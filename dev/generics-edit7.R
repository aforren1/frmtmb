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

# put the measured agreement on the page rather than an adjective
sub1("R/scales.R",
paste0("#' `fitted()`, `residuals()`, `simulate()`, `ngrps()`, `ranef()` and the\n",
       "#' coefficients of a linear predictor whose link is the identity. For a\n",
       "#' lognormal fit, `fitted()` is `exp(mu + sigma^2 / 2)`, the mean of the\n",
       "#' distribution rather than its median `exp(mu)`; the ratio between the\n",
       "#' two is `exp(sigma^2 / 2)`, which is the whole of the discrepancy a\n",
       "#' ported script sees.\n"),
paste0("#' `fitted()`, `residuals()`, `simulate()`, `ngrps()`, `ranef()` and the\n",
       "#' coefficients of a linear predictor whose link is the identity. For a\n",
       "#' lognormal fit, `fitted()` is `exp(mu + sigma^2 / 2)`, the mean of the\n",
       "#' distribution rather than its median `exp(mu)`; the ratio between the\n",
       "#' two is `exp(sigma^2 / 2)`, which is the whole of the discrepancy a\n",
       "#' ported script sees.\n",
       "#'\n",
       "#' Measured on one 400-row lognormal fit against a two-chain brms fit\n",
       "#' of the same model (`dev/generics-scale-brms.R`, seed 2026): the\n",
       "#' largest relative disagreement in `fitted()` over the 400 rows is\n",
       "#' 0.0106, which is 0.113 of brms's own posterior standard deviation\n",
       "#' for that row and 0.058 of it at the median. `predict()` disagrees\n",
       "#' by a factor of 443, because it is a different scale.\n"))

# make the contract reachable from the pages it constrains
sub1("R/predict.R",
"#' @seealso [predict.frmtmb_fit()], [residuals.frmtmb_fit()]",
"#' @seealso [predict.frmtmb_fit()], [residuals.frmtmb_fit()],\n#'   [frmtmb-scales] for which scale each method reports")
cat("DONE\n")
