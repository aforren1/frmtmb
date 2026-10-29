## Review rescorpp: a brms 2.23.0 fit of a multivariate model with a
## non-gaussian response (rescor not allowed), and one with a random
## effect shared across responses under rescor, for the layout and the
## new-level comparisons. Rscript dev/review-rescorpp-brmsmix.R <mix|re>
w <- commandArgs(TRUE)[1]
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-data.R")
out <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-out"
dir.create(out, showWarnings = FALSE)
d <- rv_data(w)
form <- rv_brms_form(w)
t0 <- proc.time()[["elapsed"]]
b <- brm(form, data = d, chains = 2, iter = 1000, warmup = 500, seed = 21,
         refresh = 0, backend = "rstan")
cat("seconds", proc.time()[["elapsed"]] - t0, "\n")
print(summary(b))
saveRDS(b, file.path(out, paste0("brms-", w, ".rds")))
