# Lane meme, fix round 2: frm() with set_prior(class = "Intercept",
# nlpar =, lb =) is unchanged. Rscript dev/meme-fix2-frm.R <base|lane>
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r15", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(2)
d <- data.frame(x = runif(60, 0, 3), w = rnorm(60))
d$y <- 2 * exp(-0.7 * d$x) + rnorm(60, 0, 0.1)
for (f in list(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE),
               bf(y ~ a * exp(-b * x), a ~ 1 + w, b ~ 1, nl = TRUE))) {
  fit <- suppressWarnings(frm(f, data = d,
    prior = set_prior("normal(1, 1)", class = "Intercept", nlpar = "a", lb = 0) +
      set_prior("", class = "Intercept", nlpar = "b", lb = 0.75)))
  cat(sprintf("%.15g\n", c(as.numeric(logLik(fit)), unlist(fit$estimates))))
}
