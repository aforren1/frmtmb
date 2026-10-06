# Punch round 1, m7: summary()$gp row labels and get_prior() rows for a
# gp() fit on each arm, for the NEWS entries that mark them breaking.
# Usage: Rscript dev/gpby-p1-sumgp.R <lane|base>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
d <- data.frame(x = stats::runif(50, 0, 5), z = stats::runif(50, 0, 3))
d$y <- sin(d$x) + stats::rnorm(50, 0, 0.3)
for (f in list(bf(y ~ gp(x)), bf(y ~ gp(x, z, iso = FALSE)))) {
  fit <- suppressWarnings(frm(f, data = d))
  cat("summary()$gp rows:", rownames(summary(fit)$gp), "\n")
  g <- as.data.frame(get_prior(f, data = d))
  cat("get_prior() classes:", unique(g$class), "\n")
}
