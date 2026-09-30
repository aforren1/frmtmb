# Reviewer: brms's data declarations of the threshold counts, grouped
# and ungrouped, sum_to_zero; and the tau_raw escape hatch spelled right.
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261005)
n <- 300
d <- data.frame(x = rnorm(n), h = factor(sample(c("p", "q"), n, TRUE)))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yb <- 1L + (u > 0)
for (f in list(brms::bf(yb | thres(gr = h) ~ x), brms::bf(y ~ x))) {
  sc <- brms::stancode(f, data = d,
                       family = brms::cumulative(threshold = "sum_to_zero"))
  L <- strsplit(sc, "\n")[[1]]
  cat(grep("int<lower=[0-9]+> nthres|nthres;", L, value = TRUE), sep = "\n")
  cat("--\n")
}
f <- frm(y ~ x, family = sratio(threshold = "sum_to_zero"), data = d,
         prior = list(tau_raw = prior_normal(0, 5)))
cat("escape hatch fits; thresholds:",
    frmtmb:::ord_threshold_values(family(f), f$estimates$tau_raw), "\n")
