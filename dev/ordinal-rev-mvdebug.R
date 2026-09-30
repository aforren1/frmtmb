# Reviewer: locate the multivariate flexible ordinal draws failure.
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib", "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261009)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 1L + (stats::rlogis(n) + d$x > -1) + (stats::rlogis(n) + d$x > 0) +
  (stats::rlogis(n) + d$x > 1)
d$y2 <- 1L + (stats::rlogis(n) + d$x > 0) + (stats::rlogis(n) + d$x > 1)
fit <- frm(bf(y ~ x) + bf(y2 ~ x), data = d, family = list(cumulative(), sratio()))
ds <- suppressWarnings(suppressMessages(frm_sample(fit, chains = 1, iter = 100,
                                                   refresh = 0, seed = 5)))
nc <- frmtmb.sample:::draws_natural_cols(ds$fit)
for (o in nc$ordinal) {
  cat("internal:", o$internal, "\n names:", o$names, "\n")
  r <- seq_along(o$internal) / 10
  cat(" map(r):", o$map(r), "\n")
}
cat("extra:", nc$extra, "\n")
cat("colnames(draws):", colnames(ds$draws), "\n")
fe <- frmtmb::brms_fixef_rows(ds$fit)$extra
for (e in fe) cat("extra row comp", e$comp, " cls", e$cls, " names", e$names, " raw len", length(e$raw), "\n")
cat("ord_delta_info:\n"); str(frmtmb::ord_delta_info(ds$fit))
