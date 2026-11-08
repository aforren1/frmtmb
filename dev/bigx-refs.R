# Lane bigx: what the reference packages do with y ~ (z + m), poisson,
# data seed 1: brms 2.23.0's design (standata(), no compile) and
# whether its Stan code centers it; lme4 2.1-0's glmer() and lm()'s
# column drop.
#   BIGX_ARM=before Rscript dev/bigx-refs.R
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
suppressPackageStartupMessages({
  library(brms)
  library(lme4)
})
cat("brms", format(packageVersion("brms")), "lme4",
    format(packageVersion("lme4")), "\n")
for (m in c(1e5, 1e6, 1e7, 1e8)) {
  d <- bigx_data("poisson", m, 1)
  sd <- standata(y ~ w + (1 | g), data = d, family = poisson())
  sc <- stancode(y ~ w + (1 | g), data = d, family = poisson())
  centers <- grepl("Xc", sc) && grepl("means_X", sc)
  r <- tryCatch({
    fit <- glmer(y ~ w + (1 | g), data = d, family = poisson())
    sprintf("glmer slope %.5f", fixef(fit)[["w"]])
  }, error = function(e) paste("glmer error:", conditionMessage(e)),
  warning = function(w) paste("glmer warning:", conditionMessage(w)))
  lmc <- coef(lm(y ~ w, data = d))
  cat(sprintf("m %g: brms K %d (columns of X with the intercept), Stan code centers X: %s\n  %s\n  lm() slope %s\n",
              m, sd$K, centers, substr(r, 1, 150), format(lmc[["w"]])))
}
