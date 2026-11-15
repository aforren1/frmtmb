# Lane parity2, item 5315: what a linear predictor records about its
# intercept under `0 + Intercept`, `center = FALSE` and the default.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/parity2-data.R")
for (f in list(bf(y ~ 0 + Intercept + x, sigma ~ z),
               bf(y ~ x, sigma ~ z, center = FALSE),
               bf(y ~ x, sigma ~ 0 + Intercept + z))) {
  fit <- suppressMessages(suppressWarnings(frm(f, family = gaussian(),
                                               data = d)))
  for (nm in names(fit$frame$linpreds)) {
    lp <- fit$frame$linpreds[[nm]]
    cat(nm, " center:", deparse1(lp$center), " param_colnames:",
        lp$param_colnames, "\n")
  }
  cat("variables:", variables(fit), "\n\n")
}
