# Review of lane parity2, item 1497: brms 2.23.0 on the strong-shrinkage
# data of dev/reviews/parity2-review-snl.R part B (saved there as
# snl-data.rds): band widths at x = -1, 0, 1 for a new level.
#   Rscript dev/reviews/parity2-review-snl-brms.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
cat("brms", as.character(packageVersion("brms")), "\n")
set.seed(31)
G <- 30
ds_ <- data.frame(x = rnorm(2 * G), g = factor(rep(1:G, each = 2)))
ds_$y <- 1 + 0.5 * ds_$x + rnorm(G)[ds_$g] + rnorm(2 * G, sd = 2.5)
fit <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g), data = ds_, chains = 2, iter = 2000, seed = 1,
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
cat("sd_g posterior median:", median(as.matrix(fit, variable = "sd_g__Intercept")), "\n")
