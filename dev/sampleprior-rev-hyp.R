# Reviewer of lane sampleprior: readers of the draws that name sd_/cor_
# now that the draws carry those columns themselves.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261118)
d <- data.frame(x = rnorm(120), g = factor(rep(1:8, each = 15)))
d$y <- 1 + rnorm(8)[d$g] + (1 + rnorm(8, 0, .5)[d$g]) * d$x + rnorm(120)
ds <- suppressMessages(suppressWarnings(frm_sample(
  frm(y ~ x + (1 + x | g), data = d), chains = 2, iter = 400, refresh = 0,
  seed = 2)))
tr <- function(lab, e) cat(lab, ":", tryCatch(e, error = function(z)
  paste("ERROR", conditionMessage(z))), "\n")
tr("hypothesis sd", { h <- hypothesis(ds, "sd_g__Intercept > 0.5", class = NULL)
  sprintf("Estimate %.4f Post.Prob %.3f vs column mean %.4f",
          h$hypothesis$Estimate, h$hypothesis$Post.Prob,
          mean(ds$draws[, "sd_g__Intercept"])) })
tr("hypothesis cor", { h <- hypothesis(ds, "cor_g__Intercept__x = 0", class = NULL)
  sprintf("Estimate %.4f", h$hypothesis$Estimate) })
tr("posterior_summary rows", paste(grep("^sd_|^cor_", rownames(posterior_summary(ds)), value = TRUE), collapse = " "))
tr("fixef", paste(round(fixef(ds)[, 1], 3), collapse = " "))
tr("ranef dims", paste(dim(ranef(ds)$g), collapse = "x"))
tr("loo", round(loo::loo(ds)$estimates[1, 1], 2))
tr("bayes_R2", round(bayes_R2(ds)[1], 4))
tr("predict newdata", paste(dim(posterior_epred(ds, newdata = d[1:3, ], ndraws = 4)), collapse = "x"))
