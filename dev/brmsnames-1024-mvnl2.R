## Two responses, each nonlinear, with exp() in the second response's
## nonlinear formula: the review's case 15 (dev/review-brmsnames-cases.R,
## rv_data(), seed 2410). On 0.74.0 it does not fit.
##   Rscript dev/brmsnames-1024-mvnl2.R base|lane > dev/brmsnames-1024-log/mvnl2-<arm>.txt
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
src <- readLines("dev/review-brmsnames-cases.R")
i <- grep("^rv_data <- function", src)
j <- i + which(src[(i + 1):length(src)] == "}")[1L]
eval(parse(text = src[i:j]))
d <- rv_data()
form <- mvbf(bf(yn ~ a * exp(b * x), a ~ 1 + z + (1 | g), b ~ 1,
                sigma ~ z, nl = TRUE),
             bf(y ~ c + d * x, c ~ exp(z), d ~ 1 + (1 | g),
                sigma ~ x, nl = TRUE), rescor = FALSE)
fit <- tryCatch(q(frm(form, data = d)), error = function(e) e)
if (inherits(fit, "error")) {
  cat("frm:", conditionMessage(fit), "\n")
} else {
  cat("fits; logLik", format(as.numeric(logLik(fit)), digits = 10), "\n")
  cat("fixef:", rownames(fixef(fit)), "\n")
}
