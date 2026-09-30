# Reviewer: brms's Stan code for an equation across two sigma links, and
# brms's reading of an equation plus a formula on the same dpar.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(11)
d <- data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50))
sc <- stancode(bf(y ~ x, sigma1 = "sigma2"), data = d,
               family = mixture(brmsfamily("gaussian", link_sigma = "softplus"),
                                gaussian()))
cat(grep("sigma", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
for (e in list(quote(bf(y ~ x, sigma1 = "sigma2", sigma1 ~ z)),
               quote(bf(y ~ x, sigma1 ~ z, sigma1 = "sigma2")),
               quote(bf(y ~ x, sigma1 = "sigma2", sigma2 = "sigma3")),
               quote(bf(y ~ x, sigma1 = "sigma2") + lf(sigma2 = "sigma3")),
               quote(bf(y ~ x, sigma1 = "sigma2") + lf(sigma1 ~ z)),
               quote(bf(cat ~ x, mub = "muc")))) {
  r <- tryCatch({eval(e); "ok (no error at bf())"},
                error = function(err) conditionMessage(err))
  cat(deparse1(e), "->", r, "\n")
}
