# Reviewer, lane dblbar: brms 2.23.0's re_formula reading on a
# (1 + f || g) fit, and its new-level predictions. Data seed 20261013;
# 2 chains of 1000, sampler seed 1. Output: dev/dblbar-rev-log/brms-refo.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(brms))
dd <- dblbar_data()
fit <- brm(bf(y ~ f + (1 + f || g)), data = dd, family = gaussian(),
           chains = 2, iter = 1000, seed = 1, refresh = 0, cores = 2,
           backend = "rstan")
sink(file.path(wt, "dev/dblbar-rev-log/brms-refo.txt"), split = TRUE)
nd_old <- dd[c(1, 5, 9, 40), ]
tc <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
full <- colMeans(posterior_epred(fit, newdata = nd_old))
for (rf in c("~ (0 + f || g)", "~ (0 + f | g)", "~ (1 + f | g)",
             "~ (f | g)", "~ (1 | g)", "~ (1 + x || g)",
             "~ (1 + f || gr(g, cor = TRUE))")) {
  r <- tc(colMeans(posterior_epred(fit, newdata = nd_old,
                                   re_formula = as.formula(rf))))
  cat(sprintf("re_formula %-32s ", rf))
  if (is.character(r)) cat(substr(r, 1, 160), "\n") else
    cat("max |diff from full|", sprintf("%.3e", max(abs(r - full))),
        " means", sprintf("%.3f", r), "\n")
}
nd_new <- nd_old
nd_new$g <- c("NEW1", "NEW2", "NEW1", "NEW3")
for (snl in c("gaussian", "old_levels", "uncertainty")) {
  pe <- posterior_epred(fit, newdata = nd_new, allow_new_levels = TRUE,
                        sample_new_levels = snl, seed = 3)
  cat(sprintf("new levels %-12s mean %s sd %s\n", snl,
              paste(sprintf("%.3f", colMeans(pe)), collapse = " "),
              paste(sprintf("%.3f", apply(pe, 2, sd)), collapse = " ")))
}
sink()
