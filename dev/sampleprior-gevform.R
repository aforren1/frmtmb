# Lane sampleprior, fix round 1: frm_sample() from a formula on a
# gen_extreme_value() model with a scalar xi, which the apply sweep of
# dev/sampleprior-table.R could not tape (data seed 20261022).
#   Rscript dev/sampleprior-gevform.R <lane|r9>
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "r9")) character(0) else
  "C:/Users/adf44/source/r/wt-sampleprior-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm", arm, "lib:", find.package("frmtmb.sample"), "\n")
set.seed(20261022)
d <- data.frame(x = stats::rnorm(60))
d$y <- 1 + 0.3 * d$x + ((-log(stats::runif(60)))^(-0.1) - 1) / 0.1
for (route in c("formula", "fit")) {
  r <- tryCatch({
    m <- if (route == "formula") bf(y ~ x) else
      frm(bf(y ~ x), family = gen_extreme_value(), data = d)
    ds <- suppressWarnings(suppressMessages(frm_sample(
      m, data = if (route == "formula") d,
      family = if (route == "formula") gen_extreme_value(),
      chains = 1, iter = 200, refresh = 0, seed = 1)))
    sprintf("OK, xi mean %.3f", mean(ds$draws[, "xi"]))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(route, "route:", r, "\n")
}
