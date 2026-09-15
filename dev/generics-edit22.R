root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (length(gregexpr(old, txt, fixed = TRUE)[[1]]) != 1L)
    stop("many in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

# ---- the table: the missing rows, measured in dev/generics-scale2.R --
sub1("R/scales.R",
paste0("#' | `predict(dpar = \"sigma\")` | that predictor | link, or response |\n",
       "#' | `fitted()` | the conditional mean | response |\n"),
paste0("#' | `predict(dpar = \"sigma\")` | that predictor | link; `type = \"response\"` gives response |\n",
       "#' | `fitted()` | the conditional mean | response |\n",
       "#' | `conditional_effects()` | `estimate__` and its band | response |\n"))

sub1("R/scales.R",
paste0("#' | `sigma()` | the residual standard deviation | response |\n"),
paste0("#' | `sigma()` | the residual SD, `NA` if it varies by row | response |\n",
       "#' | `hypothesis()` | `estimate` and its interval | link, per dpar |\n",
       "#' | `posterior_summary()` | a summary of a draws MATRIX | that matrix's own |\n",
       "#' | `pp_check()` | the outcome against replicates | response |\n",
       "#' | `bayes_R2()` | refuses on a `frmtmb_fit` | neither |\n"))

# ---- the identity's two conditions ----------------------------------
sub1("R/scales.R",
paste0("#' `fitted()`, `residuals()`, `simulate()`, `ngrps()`, `ranef()` and the\n",
       "#' coefficients of a linear predictor whose link is the identity. For a\n",
       "#' lognormal fit, `fitted()` is `exp(mu + sigma^2 / 2)`, the mean of the\n",
       "#' distribution rather than its median `exp(mu)`; the ratio between the\n",
       "#' two is `exp(sigma^2 / 2)`, which is the whole of the discrepancy a\n",
       "#' ported script sees.\n"),
paste0("#' `fitted()`, `residuals()`, `simulate()`, `ngrps()`, `ranef()`,\n",
       "#' `conditional_effects()` and the coefficients of a linear predictor\n",
       "#' whose link is the identity.\n",
       "#'\n",
       "#' For a lognormal fit with a CONSTANT sigma and no truncation,\n",
       "#' `fitted()` is `exp(mu + sigma^2 / 2)`, the mean of the distribution\n",
       "#' rather than its median `exp(mu)`; the ratio between the two is\n",
       "#' `exp(sigma^2 / 2)`, which is the whole of the discrepancy a ported\n",
       "#' script sees. **Both conditions are load-bearing**, and neither is\n",
       "#' obvious from the formula:\n",
       "#'\n",
       "#' * With a distributional sigma (`bf(y ~ x, sigma ~ x)`), `sigma()`\n",
       "#'   returns `NA` with a warning, because there is no one number to\n",
       "#'   return, and the formula gives `NA` with it. Use\n",
       "#'   `predict(dpar = \"sigma\", type = \"response\")`, which reproduces\n",
       "#'   `fitted()` exactly (`identical()` is `TRUE`).\n",
       "#' * Under truncation, `fitted()` is the TRUNCATED mean and the\n",
       "#'   formula is simply wrong. On a `y | trunc(lb = 2000)` lognormal\n",
       "#'   fit of the same design the naive formula is out by a relative\n",
       "#'   0.1150 at the worst row (`dev/generics-scale2.R`).\n"))

# ---- N10: 443 is a median over rows ---------------------------------
sub1("R/scales.R",
paste0("#' largest relative disagreement in `fitted()` over the 400 rows is\n",
       "#' 0.0106, which is 0.113 of brms's own posterior standard deviation\n",
       "#' for that row and 0.058 of it at the median. `predict()` disagrees\n",
       "#' by a factor of 443, because it is a different scale.\n"),
paste0("#' largest relative disagreement in `fitted()` over the 400 rows is\n",
       "#' 0.0106, which is 0.113 of brms's own posterior standard deviation\n",
       "#' for that row and 0.058 of it at the median. `predict()` disagrees\n",
       "#' by a factor of 443 at the MEDIAN row, and 768.7 at row 1, because\n",
       "#' it is a different scale rather than a different answer.\n"))

# ---- N1: the page has to be findable --------------------------------
sub1("R/scales.R",
"#' @name frmtmb-scales\n#' @keywords internal\nNULL\n",
"#' @name frmtmb-scales\nNULL\n")
sub1("_pkgdown.yml",
"  - sigma.frmtmb_fit\n",
"  - sigma.frmtmb_fit\n  - frmtmb-scales\n")
cat("DONE\n")
