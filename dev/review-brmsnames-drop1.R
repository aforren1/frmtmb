## drop1() and step() on a plain exp() term. Seed 3.
##   Rscript dev/review-brmsnames-drop1.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
set.seed(3)
n <- 200
d <- data.frame(x = runif(n), z = rnorm(n))
d$y <- 1 + 0.8 * exp(d$x) + 0.05 * d$z + rnorm(n)
fit <- q(frm(y ~ exp(x) + z, data = d))
print(q(drop1(fit)))
red <- q(frm(y ~ z, data = d))
cat("AIC of y ~ z refitted:", AIC(red), "\n")
s <- tryCatch(q(step(fit, trace = 0)), error = function(e) e)
cat("step() keeps:", if (inherits(s, "error")) conditionMessage(s) else
  deparse1(formula(s)), "\n")
