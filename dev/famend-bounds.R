# Lane famend: where a fit keeps the bounds it ran under (for the
# zero-density-end assertion of test-prior-table-brms.R)
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + 0.3 * d$z + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
for (pr in list(set_prior("beta(2, 2)", class = "b", coef = "z"),
                set_prior("beta(2, 2)", class = "b", coef = "z", lb = 0, ub = 1))) {
  f <- suppressWarnings(frm(yneg ~ z, data = d, prior = pr))
  str(f$cache$bounds)
  print(ls(f$cache))
}
