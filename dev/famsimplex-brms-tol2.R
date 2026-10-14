# Control arm of dev/famsimplex-brms-tol.R: the same brms dirichlet call
# on the clean data and on one row off by 1e-9 and by 1e-7, in ONE
# process (one compile). Output: dev/famsimplex-brms-tol2.log
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
dd <- fs_data_dirichlet(seed = 11)
fit0 <- brm(bf(Y ~ x), data = dd, family = dirichlet(), chains = 1,
            iter = 400, seed = 1, refresh = 0, backend = "rstan")
for (eps in c(0, 1e-9, 1e-7)) {
  d <- dd
  d$Y[1, 1] <- d$Y[1, 1] + eps
  cat("\n== row 1 off by", eps, "\n")
  f <- tryCatch(update(fit0, newdata = d, recompile = FALSE, refresh = 0,
                       seed = 1),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(f)) cat(f, "\n") else {
    cat("draws:", posterior::ndraws(f), "\n")
  }
}
