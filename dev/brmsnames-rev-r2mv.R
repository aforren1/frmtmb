## Reviewer: bayes_R2() on a two-response model, base and lane.
##   Rscript dev/brmsnames-rev-r2mv.R base|lane   (data of dev/brmsnames-rev-data.R, sampler seed 3)
arm <- commandArgs(trailingOnly = TRUE)[1L]
lib <- if (arm == "lane") c("C:/Users/adf44/source/r/brmsnames-lib",
                            "C:/Users/adf44/source/r/rellib-r3") else
  "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
source("dev/brmsnames-rev-data.R")
fit <- q(frm(rev_models("frmtmb")$C4$f, family = gaussian(), data = rev_dd))
ds <- q(frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 3))
cat(arm, "\n"); print(tryCatch(bayes_R2(ds), error = function(e) conditionMessage(e)))
print(tryCatch(dim(bayes_R2(ds, summary = FALSE)), error = function(e) conditionMessage(e)))
print(tryCatch(bayes_R2(ds, resp = "y2"), error = function(e) conditionMessage(e)))
