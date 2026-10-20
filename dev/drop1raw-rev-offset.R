# Reviewer: a variable read only inside offset() is left out of
# raw_vars; does an update() that reads it outside offset() work?
# Arg: lane. Log: dev/drop1raw-rev-out/offset-lane.txt (seed 5, n = 200)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, expo = runif(n, 0.5, 3))
d$cnt <- rpois(n, d$expo * exp(0.3 + 0.4 * d$x + 0.3 * log(d$x2)))
fo <- q(frm(bf(cnt ~ x + log(x2) + offset(log(expo))), family = poisson(),
            data = d))
g <- glm(cnt ~ x + log(x2) + offset(log(expo)) + log(expo), poisson, d)
r <- tryCatch(q(update(fo, ~ . + log(expo))), error = function(e) e)
cat("data in place: update(fo, ~ . + log(expo)):",
    if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
      sprintf("logLik %.8f | glm %.8f", as.numeric(logLik(r)),
              as.numeric(logLik(g))), "\n")
fo_gone <- local({
  dg <- d
  q(frm(bf(cnt ~ x + log(x2) + offset(log(expo))), family = poisson(),
        data = dg))
})
r <- tryCatch(q(update(fo_gone, ~ . + log(expo))), error = function(e) e)
cat("data gone: update(fo, ~ . + log(expo)):",
    if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
      sprintf("logLik %.8f", as.numeric(logLik(r))), "\n")
expo <- rep(2, n)
r <- tryCatch(q(update(fo_gone, ~ . + log(expo))), error = function(e) e)
cat("data gone, global expo = 2: update(fo, ~ . + log(expo)):",
    if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
      sprintf("logLik %.8f", as.numeric(logLik(r))), "\n")
