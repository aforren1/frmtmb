.libPaths(c("C:/Users/adf44/source/r/wt-splitsurv-lib","C:/Users/adf44/source/r/rellib-r11","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
r <- names(frmtmb:::frmtmb_links); cat(length(r), "\n"); print(r)
cl <- frmtmb.survival:::cox_links; print(cl)
hand <- c("identity", "log", "logit", "cloglog", "probit",
              "probit_approx", "cauchit", "softit", "inverse", "logm1",
              "log1p", "tan_half", "power12", "softplus", "squareplus",
              "sqrt", "1/mu^2")
cat("roster minus hand:", setdiff(r, hand), "\n")
cat("refused by core loop:", length(setdiff(r, cl)), " by survival loop:", length(setdiff(hand, cl)), "\n")
