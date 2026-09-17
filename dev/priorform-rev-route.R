# Reviewer, lane wt-priorform: do ledger P15 and P27 (recorded as
# divergences because the fit route is flat) transfer on route = "sample"?
#   Rscript dev/priorform-rev-route.R     (lane build, frmtmb.sample attached)
.libPaths(c("C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
dat <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)), x = rnorm(10), g = rep(1:2, 5))
bform <- bf(mvbind(y1, y2) ~ x + (x|ID1|g)) + set_rescor(TRUE)
r <- tryCatch({p <- default_prior(bform, dat, family = gaussian(), route = "sample")
  p[p$class == "rescor", "prior"]}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("P15 rescor prior on route sample:", r, " (brms: lkj(1))\n")
dat2 <- data.frame(y = rep(c(1, 3), each = 5), off = 10)
r2 <- tryCatch({p <- default_prior(y ~ 1 + offset(off), dat2, route = "sample")
  p$prior[p$class == "Intercept"]}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("P27 Intercept prior on route sample:", r2, " (brms: student_t(3, -8, 2.5))\n")
