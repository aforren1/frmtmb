# brms 2.23.0's simplex check is a MEAN over rows (all.equal(rowSums(Y),
# 1)), so one row off the simplex passes standata() in a data set of 300.
# What does sampling then do? Output: dev/famsimplex-brms-tol.log
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
dd <- fs_data_dirichlet(seed = 11)
dd$Y[1, 1] <- dd$Y[1, 1] + 1e-7
cat("row 1 sum - 1:", format(sum(dd$Y[1, ]) - 1, digits = 3),
    "; mean |rowSums - 1|:", format(mean(abs(rowSums(dd$Y) - 1)), digits = 3),
    "\n")
for (fam in list(dirichlet(), logistic_normal())) {
  cat("\n==", fam$family, "\n")
  r <- tryCatch({
    fit <- brm(bf(Y ~ x), data = dd, family = fam, chains = 1, iter = 400,
               seed = 1, refresh = 0, backend = "rstan")
    print(summary(fit))
    "sampled"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(r, "\n")
}
