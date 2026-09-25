# brms 2.23.0 default priors on reviewer data, for the sample-route check.
.libPaths(c("C:/Users/adf44/source/r/pinlib", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-cases.R")
d <- r2_data()
sh <- function(nm, f, fam) {
  cat("==", nm, "\n")
  p <- as.data.frame(default_prior(f, data = d, family = fam))
  p <- p[p$prior != "", ]
  print(p[, c("prior", "class", "coef", "group", "resp", "dpar")], row.names = FALSE)
}
sh("cat4re", bf(cat4 ~ x + (1 | g)), categorical())
sh("mix3", bf(ym3 ~ x), mixture(gaussian, gaussian, gaussian))
sh("mixdiff_re", bf(ym ~ x + (1 | g)), mixture(gaussian, student))
cat("mad(ym3) =", mad(d$ym3), " median =", median(d$ym3), "\n")
cat("mad(ym) =", mad(d$ym), " median =", median(d$ym), "\n")
