# Reviewer, lane dblbar: where did frmtmb.sample name the zero cor_ of a
# diag block on the base? variables(), VarCorr(), hypothesis(). Data
# seed 20261013; sampler seed 2, 1 chain of 400.
# Usage: Rscript ... base|lane. Output: stdout
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
fit <- frm(bf(y ~ f + diag(1 + f | g)), family = gaussian(), data = dd)
dr <- frm_sample(fit, chains = 1, iter = 400, seed = 2, refresh = 0)
cat("variables cor:", grep("cor", variables(dr), value = TRUE), "\n")
cat("VarCorr names:", names(VarCorr(dr)$g), "\n")
h <- tryCatch(hypothesis(dr, "cor_g__Intercept__fb = 0", class = NULL),
              error = function(e) paste("ERROR:", substr(conditionMessage(e),
                                                         1, 80)))
print(if (is.character(h)) h else h$hypothesis[, 1:3])
