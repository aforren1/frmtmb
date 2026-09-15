root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl(old, txt, fixed = TRUE)) stop("no match: ", old)
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (n != 1L) stop(n, " matches: ", old)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
}
# markdown table rows cannot be wrapped, so the cells are shortened
# instead to keep R/ inside 80 columns
sub1("R/scales.R",
"#' | `predict(dpar = \"sigma\")` | that parameter's predictor | link, or response with `type = \"response\"` |",
"#' | `predict(dpar = \"sigma\")` | that predictor | link, or response |")
sub1("R/scales.R",
"#' | `coef()`, `fixef()`, `ranef()` | coefficients and modes | link, per parameter |",
"#' | `coef()`, `fixef()`, `ranef()` | coefficients, modes | link, per dpar |")
sub1("R/scales.R",
"#' | `confint()`, `vcov()`, `summary()` | the same coefficients | link, per parameter |",
"#' | `confint()`, `vcov()`, `summary()` | the same | link, per dpar |")
sub1("R/scales.R",
"#' | `VarCorr()` | covariance of a linear predictor | link, per parameter |",
"#' | `VarCorr()` | covariance of a predictor | link, per dpar |")
sub1("R/scales.R",
"#' | `logLik()`, `AIC()`, `BIC()`, `deviance()` | the fitted likelihood | neither; the data's own |",
"#' | `logLik()`, `AIC()`, `BIC()` | the fitted likelihood | the data's own |")
cat("DONE\n")
