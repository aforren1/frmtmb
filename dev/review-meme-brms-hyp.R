# Review of lane meme: brms 2.23.0's hypothesis() Evid.Ratio on me()
# hyperparameters, and its update(prior =) merge. One small fit.
#   Rscript dev/review-meme-brms-hyp.R > dev/review-meme-out/brms-hyp.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
fit <- suppressWarnings(suppressMessages(brm(
  bf(y ~ me(x, sx) + me(z, sz) + w), data = d,
  prior = c(set_prior("normal(0, 1)", class = "meanme"),
            set_prior("exponential(1)", class = "sdme")),
  sample_prior = "yes", chains = 2, iter = 1000, seed = 3, refresh = 0)))
cat("prior draws:", names(prior_draws(fit)), "\n")
cat("prior_draws(variable = meanme_mex):",
    names(prior_draws(fit, variable = "meanme_mex")), "\n")
for (h in c("meanme_mex = 0", "meanme_mex > 0", "sdme_mex = 1",
            "corme__mex__mez = 0", "sigma = 1")) {
  r <- tryCatch(hypothesis(fit, h, class = NULL)$hypothesis[
    , c("Estimate", "Evid.Ratio", "Post.Prob")],
    error = function(e) conditionMessage(e))
  cat("--", h, "\n")
  print(r)
}
