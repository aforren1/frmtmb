source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(brms))
# the record's construction: rename_pars() first (dev/adefects-repro-base.R)
b1 <- brms:::rename_pars(get("brmsfit_example1", envir = asNamespace("brms")))
b_full <- b1$data[1:3, ]
b_nogrp <- b_full[, setdiff(names(b_full), "visit"), drop = FALSE]
for (fn in c("posterior_epred", "posterior_predict", "log_lik")) {
  g <- get(fn, envir = asNamespace("brms"))
  for (anl in c(TRUE, FALSE)) {
    r <- tryCatch(g(b1, newdata = b_nogrp, allow_new_levels = anl),
                  error = function(e) e)
    cat(sprintf("  %-18s allow_new_levels = %-5s -> %s\n", fn, anl,
                if (inherits(r, "condition"))
                  paste("ERR:", substr(gsub("[\r\n]+", " ",
                                            conditionMessage(r)), 1, 60))
                else paste("OK", paste(dim(r), collapse = " x "))))
  }
}
cat("  validate_newdata fills visit with:",
    paste(brms:::validate_newdata(b_nogrp, b1,
                                  allow_new_levels = TRUE)$visit,
          collapse = " "), "\n")
