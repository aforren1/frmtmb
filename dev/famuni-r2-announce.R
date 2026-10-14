# The prior announcement for a scalar GEV xi after fix round 2 (review
# round 2, item 4). Output: dev/famuni-out/r2-announce.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1401); y <- (stats::rexp(40)^0.45 - 1) / -0.45
fit <- frm(y ~ 1, family = gen_extreme_value(), data = data.frame(y = y))
msgs <- character(0)
invisible(withCallingHandlers(frm_sample(fit, chains = 1, iter = 200, seed = 1, refresh = 0),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
  warning = function(w) invokeRestart("muffleWarning")))
cat(msgs, sep = "")
