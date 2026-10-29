## Lane rescorpp: fit brms 2.23.0 on each data set and save the fit.
##   Rscript dev/rescorpp-brmsfit.R <g3|gs|st|g2f>
w <- commandArgs(TRUE)[1]
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-data.R")
d <- rescorpp_data(rescorpp_dname(w))
fam <- get(rescorpp_family(w), asNamespace("brms"))()
out <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-brmsfits"
dir.create(out, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]
b <- brm(eval(rescorpp_forms[[w]]), data = d, family = fam,
         chains = 2, iter = 2000, warmup = 1000, seed = 21, refresh = 0,
         backend = "rstan")
cat("seconds", proc.time()[["elapsed"]] - t0, "\n")
print(summary(b))
saveRDS(b, file.path(out, paste0(w, ".rds")))
