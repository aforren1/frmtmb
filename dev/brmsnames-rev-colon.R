## Reviewer, claim 3: an interaction name in a hypothesis string, which
## brms rewrites (`:` to `___`) and frmtmb parses as R's `:` operator.
##   Rscript dev/brmsnames-rev-colon.R base|lane     data seed 8
arm <- commandArgs(trailingOnly = TRUE)[1L]
lib <- if (arm == "lane") c("C:/Users/adf44/source/r/brmsnames-lib",
                            "C:/Users/adf44/source/r/rellib-r3") else
  "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
set.seed(8)
n <- 400
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "e"), n, TRUE)))
d$y <- 1 + 0.3 * d$x + 0.8 * (d$f == "e") + 1.5 * d$x * (d$f == "e") + rnorm(n)
fit <- q(frm(bf(y ~ x * f), family = gaussian(), data = d))
print(round(fixef(fit)$mu, 4))
est <- function(h) {
  r <- tryCatch(q(hypothesis(fit, h)), error = function(e) e)
  if (inherits(r, "error")) return(paste("ERROR:", substr(conditionMessage(r), 1, 80)))
  if (inherits(r, "brmshypothesis")) r$hypothesis$Estimate else r$estimate
}
for (h in c("x:fe > 0", "`x:fe` > 0", "b_x:fe > 0")) cat(arm, h, "->", est(h), "\n")
