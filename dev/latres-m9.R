# m9: simulate(fit, newdata) on one series, against the review's 2.25 s at 20000
# rows and 49 s at 100000. Seed 1.
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
d <- expand.grid(t = 1:6, g = factor(1:10)); d$x <- rnorm(60)
d$y <- rpois(60, exp(.5 + .3 * d$x))
f <- suppressMessages(suppressWarnings(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = poisson(), data = d)))
for (n in c(20000L, 100000L)) {
  nd <- data.frame(t = seq_len(n), g = factor(rep("a", n)), x = 0)
  tt <- system.time(s <- simulate(f, newdata = nd, allow_new_levels = TRUE))
  cat(n, "rows:", tt[["elapsed"]], "s\n")
}
