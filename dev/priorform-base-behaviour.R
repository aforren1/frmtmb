# What the base commit DOES with the three formulas this lane refuses:
# not whether it errors, but which model it fits instead.
#   Rscript dev/priorform-base-behaviour.R ref|lane
args <- commandArgs(trailingOnly = TRUE)
lib <- switch(if (length(args)) args[1] else "ref",
  ref = "C:/Users/adf44/source/r/rellib-r3",
  lane = "C:/Users/adf44/source/r/priorform-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
`%||%` <- function(a, b) if (is.null(a)) b else a
cat("frmtmb from", find.package("frmtmb"), "\n\n")

set.seed(20260916)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:12, 20)))
d$y <- 1 + 0.5 * d$x + rnorm(12, 0, 0.8)[d$g] +
  rnorm(12, 0, 0.5)[d$g] * d$x + rnorm(n)

try_show <- function(label, expr) {
  cat("---", label, "---\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  print(r)
  cat("\n")
}

try_show("y ~ ~x: fixed effects of the fit", {
  f <- suppressWarnings(frm(y ~ ~x, d))
  list(fixef = fixef(f), logLik = as.numeric(logLik(f)))
})
try_show("y ~ x, for comparison", {
  f <- frm(y ~ x, d)
  list(fixef = fixef(f), logLik = as.numeric(logLik(f)))
})

try_show("(1 | g) + (x | g): VarCorr and logLik", {
  f <- suppressWarnings(frm(y ~ x + (1 | g) + (x | g), d))
  list(VarCorr = VarCorr(f), logLik = as.numeric(logLik(f)))
})
try_show("(x | g), for comparison", {
  f <- frm(y ~ x + (x | g), d)
  list(VarCorr = VarCorr(f), logLik = as.numeric(logLik(f)))
})

d$yo <- cut(d$y, quantile(d$y, 0:3 / 3), include.lowest = TRUE,
            ordered_result = TRUE)
d$c3 <- factor(rep(c("a", "b", "c"), length.out = n))
try_show("yo ~ x * cs(c3), sratio: design columns", {
  fr <- frm(yo ~ x * cs(c3), d, family = sratio(), dry_run = "frame")
  lapply(fr[["linpreds"]], function(lp) {
    list(X = colnames(lp[["X"]]),
         cs = vapply(lp[["cs"]] %||% list(), `[[`, "", "label"))
  })
})
try_show("yo ~ x + cs(c3), for comparison", {
  fr <- frm(yo ~ x + cs(c3), d, family = sratio(), dry_run = "frame")
  lapply(fr[["linpreds"]], function(lp) {
    list(X = colnames(lp[["X"]]),
         cs = vapply(lp[["cs"]] %||% list(), `[[`, "", "label"))
  })
})
