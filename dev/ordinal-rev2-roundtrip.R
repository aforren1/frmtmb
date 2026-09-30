# Reviewer re-check: refit from par_template() at the optimum, base vs
# lane, flexible cumulative; and the grouped sum-to-zero disc model with
# its own formula. Usage: ... <base|lane>. Data seed 20261015.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r4", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261015)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), h = factor(sample(c("p", "q"), n, TRUE)))
d$yb <- 1L + (stats::rlogis(n) / exp(0.3 * d$z) + d$x > 0)
d$y <- 1L + findInterval(stats::rlogis(n) + d$x, c(-1, 0, 1))
f <- frm(y ~ x, family = cumulative(), data = d)
g <- frm(y ~ x, family = cumulative(), data = d, start = par_template(f))
cat(arm, "flexible: refit from template, logLik rel diff",
    abs(as.numeric(logLik(g)) - as.numeric(logLik(f))) / abs(as.numeric(logLik(f))),
    " par identical:", identical(g$opt$par, f$opt$par), "\n")
if (arm == "lane") {
  fo <- bf(yb | thres(gr = h) ~ x, disc ~ 0 + z)
  f2 <- suppressMessages(frm(fo, data = d, family = sratio(threshold = "sum_to_zero")))
  g2 <- suppressMessages(frm(fo, data = d, family = sratio(threshold = "sum_to_zero"),
                             start = par_template(f2)))
  cat("gr_stz_disc: refit from template, logLik rel diff",
      abs(as.numeric(logLik(g2)) - as.numeric(logLik(f2))) / abs(as.numeric(logLik(f2))), "\n")
}
