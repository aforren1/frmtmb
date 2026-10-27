# Reviewer, fix round 2: brms 2.23.0's sigma ~ x under ar(cov = TRUE).
# Is its Intercept_sigma the log of the innovation sd?
# Log: dev/review-acprior-log/r2-brmssigx.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- data.frame(y = rnorm(40), x = rnorm(40), t = rep(1:10, 4),
                g = rep(1:4, each = 10))
sc <- stancode(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), data = d,
               prior = set_prior("normal(0, 1)", class = "Intercept",
                                 dpar = "sigma"))
ln <- strsplit(sc, "\n")[[1]]
cat(grep("Intercept_sigma|sigma =|exp[(]sigma|_het_lpdf[(]Y|cholesky_cor_ar1[(]ar|mat [.][/]", ln,
         value = TRUE), sep = "\n")
