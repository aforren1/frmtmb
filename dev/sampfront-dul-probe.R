# Lane sampfront: what drop_unused_levels = FALSE leaves in a fit and draws.
.libPaths(c("C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1655)
dd <- data.frame(g = factor(rep(1:6, each = 6L), levels = 1:7),
                 x = stats::rnorm(36))
dd$y <- 1 + 0.5 * dd$x + stats::rnorm(6L, 0, 0.5)[dd$g] + stats::rnorm(36)
f <- frm(bf(y ~ x + (1 | g)), data = dd, drop_unused_levels = FALSE)
print(length(f$estimates$b))
print(nlevels(f$data$g))
print(grep("^r_g", brms_par_labels(f), value = TRUE))
f0 <- frm(bf(y ~ x + (1 | g)), data = dd)
print(length(f0$estimates$b))
