source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages({ library(frmtmb); library(frmtmb.sample) })
ds <- structure(list(), class = "frmtmb_draws")
for (e in list(quote(posterior_epred(ds, sample_new_levels = "gaussian")),
               quote(posterior_predict(ds, incl_autocor = FALSE)))) {
  m <- tryCatch(eval(e), error = function(err) conditionMessage(err))
  cat(deparse(e), "\n  ", gsub("[\r\n]+", " ", m), "\n")
  cat("   lists allow_new_levels:",
      grepl("allow_new_levels", m, fixed = TRUE), "\n")
}
