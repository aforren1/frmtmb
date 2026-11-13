# Review of lane quadseam: is the 3466 runaway of quadseam's test the case
# lane wcheck refuses? The test's data, fitted on wcheck's library.
#   Rscript quadseam-review-wcheck.R <libpath>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(20261111)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
d$y1 <- 1 + 0.5 * d$x + rnorm(n)
d$y2 <- -1 + 0.3 * d$x + rnorm(n)
d$w <- 1.5
fo <- mvbf(bf(y1 | weights(w) ~ x + (1 | obs)),
           bf(y2 | weights(w) ~ x + (1 | obs)))
r <- tryCatch(suppressWarnings(suppressMessages(
  frm(fo, family = gaussian(), data = d, quadrature = TRUE))),
  error = function(e) e)
if (inherits(r, "error")) {
  cat("REFUSED:", class(r)[1], substr(conditionMessage(r), 1, 300), "\n")
} else {
  cat("FIT logLik", as.numeric(logLik(r)), "betad", r$estimates$betad, "\n")
}
