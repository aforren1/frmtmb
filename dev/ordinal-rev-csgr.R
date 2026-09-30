# Reviewer: cs() beside thres(gr = ), where frame.R now reads the family's
# threshold count (a vector under grouping). Usage: ... <base|lane>
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r4", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261011); n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), h = factor(sample(c("p", "q"), n, TRUE)))
d$y <- 1L + (stats::rlogis(n) + d$x > -1) + (stats::rlogis(n) + d$x > 0) +
  (stats::rlogis(n) + d$x > 1)
d$y[d$h == "p"] <- pmin(d$y[d$h == "p"], 3L)
r <- tryCatch({
  f <- frm(y | thres(gr = h) ~ x + cs(z), family = sratio(), data = d)
  cat(arm, "fits; logLik", as.numeric(logLik(f)), "; cs coefs:",
      grep("^z", rownames(fixef(f)), value = TRUE), "\n")
}, error = function(e) cat(arm, "ERROR:", conditionMessage(e), "\n"))
