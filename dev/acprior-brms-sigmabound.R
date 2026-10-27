# brms 2.23.0: does a user lb/ub on sigma under ar(cov = TRUE) bound its
# innovation sd? Log: dev/acprior-out/brms-sigmabound.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- data.frame(y = rnorm(40), x = rnorm(40), t = rep(1:10, 4),
                g = rep(1:4, each = 10))
for (pr in list(set_prior("", class = "sigma", lb = 0.2),
                set_prior("normal(1, 1)", class = "sigma", ub = 3))) {
  sc <- stancode(y ~ x + ar(time = t, gr = g, cov = TRUE), data = d,
                 prior = pr)
  cat(grep("sigma;|lprior.*sigma|cholesky_cor_ar1[(]ar|time_hom_lpdf",
           strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
  cat("--\n")
}
