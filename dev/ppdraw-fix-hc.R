# Fix round: the sweep's hurdle_cumulative boundary row, z 8631: how many
# zeros the draws have against fitted()'s P(Y = 0).
# Usage: Rscript dev/ppdraw-fix-hc.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
d <- gen("hurdle_cumulative", "boundary")
fit <- suppressMessages(suppressWarnings(
  frm(bf(y ~ x), family = hurdle_cumulative(), data = d)))
set.seed(5)
p <- suppressMessages(predict(fit, ndraws = 1000, summary = FALSE))
fe <- fitted(fit)[, "Estimate", ]
cat("fitted P(Y = 0): mean", signif(mean(fe[, 1]), 4), "\n")
cat("share of draws at 0:", mean(p == 0), " cells", sum(p == 0), "of",
    length(p), "\n")
cat("hu estimate", fit$estimates$betad[["hu_(Intercept)"]], "\n")
