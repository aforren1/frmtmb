# Review of lane names, 1543: brms 2.23.0's printed formula lines for
# the shapes of dev/reviews/names-review-dataname.R (no compile).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
sh <- function(nm, f) {
  cat("--", nm, "\n")
  print(f)
}
sh("dist", bf(y ~ x + (1 | g), sigma ~ z + (1 | g)))
sh("nl", bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE))
sh("pfix nu", bf(y ~ x, nu = 4))
sh("mix", bf(y ~ x, sigma1 = "sigma2"))
sh("mixth", bf(y ~ x, theta1 = 0.3))
sh("mv_pf", bf(y ~ x, sigma ~ x) + bf(x ~ 1) + set_rescor(FALSE))
sh("nlf", bf(y ~ a + b * x, nlf(a ~ 1 + 0 * x), b ~ 1, nl = TRUE))
sh("lf", bf(y ~ x) + lf(sigma ~ x))
sh("ar", bf(y ~ x + ar(p = 1)))
sh("pfix 1/3", bf(y ~ x, nu = 1 / 3))
# the summary of a fit: the Formula block is print(x$formula) inside
# summary's print, so the lines above are what summary() shows
