# Reviewer of lane optima, re-check: check_laplace() on a mo() fit, on
# the arm named by the argument (base = rellib-r6, lane = wt-optima-lib).
#   Rscript dev/optima-rev2-laplace-base.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm, base = "C:/Users/adf44/source/r/rellib-r6",
               lane = c("C:/Users/adf44/source/r/wt-optima-lib",
                        "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
d <- data.frame(income, ls)
fit <- frm(bf(ls ~ mo(income)), data = d, family = gaussian())
r <- tryCatch(suppressWarnings(check_laplace(fit, chains = 2, iter = 1000,
                                             seed = 2, refresh = 0)),
              error = function(e) e)
cat("arm", arm, "\n")
if (inherits(r, "error")) cat("ERROR", conditionMessage(r), "\n") else
  print(format(r, digits = 3))
ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 300, seed = 2,
                                  refresh = 0, .diagnostic = TRUE))
cat("internal:", colnames(frmtmb.sample:::draws_internal_matrix(ds)), "\n")
cat("outer cols:", frmtmb.sample:::draws_outer_cols(ds), "\n")
cat("opt$par:", names(fit$opt$par), "\n")
