source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages({ library(frmtmb); library(frmtmb.sample) })
ds <- structure(list(), class = "frmtmb_draws")
for (call in c("posterior_epred", "posterior_predict", "log_lik")) {
  r <- tryCatch(do.call(call, list(ds, newdata = data.frame(x = 1),
                                   allow_new_levels = TRUE)),
                error = function(e) e)
  cat(sprintf("%-20s allow_new_levels = TRUE -> %s\n", call,
              substr(gsub("[\r\n]+", " ", conditionMessage(r)), 1, 100)))
}
r <- tryCatch(fitted(structure(list(), class = "frmtmb_fit"),
                     newdata = data.frame(x = 1), allow_new_levels = TRUE),
              error = function(e) e)
cat(sprintf("%-20s allow_new_levels = TRUE -> %s\n", "fitted.frmtmb_fit",
            substr(gsub("[\r\n]+", " ", conditionMessage(r)), 1, 100)))
