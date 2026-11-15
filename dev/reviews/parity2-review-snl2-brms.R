# Review of lane parity2, item 1497: brms 2.23.0 on snl2-data.rds of
# dev/reviews/parity2-review-snl2.R.
#   Rscript dev/reviews/parity2-review-snl2-brms.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
cat("brms", as.character(packageVersion("brms")), "\n")
d <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/reviews/parity2-review-out/snl2-data.rds")
fit <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g), data = d, chains = 2, iter = 2000, seed = 1,
      refresh = 0, silent = 2)))
ic <- list(x = c(-1, 0, 1))
wid <- function(ce) ce$x$upper__ - ce$x$lower__
for (snl in c("uncertainty", "gaussian", "old_levels")) {
  set.seed(2)
  ce <- conditional_effects(fit, "x", int_conditions = ic, re_formula = NULL,
                            sample_new_levels = snl)
  cat("brms", snl, "width:", format(wid(ce), digits = 4), "\n")
}
r <- as.matrix(fit, variable = "^r_g\\[", regex = TRUE)
cat("sd over draws and levels of r_g draws:", sd(as.vector(r)), "\n")
cat("sd of posterior means of r_g:", sd(colMeans(r)), "\n")
cat("sd_g posterior median:",
    median(as.matrix(fit, variable = "sd_g__Intercept")), "\n")
