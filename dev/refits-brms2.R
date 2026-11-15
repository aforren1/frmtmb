# brms 2.23.0: data that lack a column the formula removes
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
go <- function(lab, expr) {
  cat("\n-- ", lab, " --\n", sep = "")
  r <- tryCatch(expr, error = function(e)
    paste("ERROR:", gsub("[\r\n]+", " ", conditionMessage(e))))
  print(r)
  invisible(r)
}
set.seed(3)
d4 <- data.frame(x = rnorm(80), z = rnorm(80), w = rnorm(80))
d4$y <- d4$x + rnorm(80)
d5 <- d4; d5$w <- NULL
go("y ~ (x + z + w) - w without w", colnames(standata(
  y ~ (x + z + w) - w, data = d5)$X))
go("y ~ x + w - w without w", colnames(standata(
  y ~ x + w - w, data = d5)$X))
go("bf(y ~ x, sigma ~ (x + z + w) - w) without w", names(standata(
  bf(y ~ x, sigma ~ (x + z + w) - w), data = d5)))
go("y ~ . - w with w: X columns", colnames(standata(
  y ~ . - w, data = d4)$X))
go("y ~ . - w + w:x with w", colnames(standata(
  y ~ . - w + w:x, data = d4)$X))
