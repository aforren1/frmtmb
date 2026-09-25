# Lane wt-mvprior: brms 2.23.0's default_prior() rows for a nonlinear
# sigma, the shape of test-nlf.R:176, beside frmtmb's on both arms.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
set.seed(3)
d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60))
if (identical(Sys.getenv("MVPRIOR_WHO"), "brms")) {
  suppressMessages(library(brms))
  tab <- default_prior(bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1),
                       data = d)
} else {
  suppressMessages(library(frmtmb))
  cat("frmtmb from", find.package("frmtmb"), "\n")
  tab <- get_prior(bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1) +
                     gaussian(), data = d)
}
tab <- as.data.frame(tab)
print(tab[tab$class != "theta", c("prior", "class", "coef", "dpar", "nlpar")],
      row.names = FALSE)
