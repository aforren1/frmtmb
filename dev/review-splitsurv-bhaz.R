arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-splitsurv-lib","C:/Users/adf44/source/r/rellib-r11","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb)); if (arm == "lane") suppressMessages(library(frmtmb.survival))
set.seed(1); dd <- data.frame(x = rnorm(150)); dd$time <- rexp(150, exp(-0.5 + 0.7 * dd$x))
cat(arm, ":", tryCatch({frm(bf(time | bhaz(df = 6) ~ x), family = cox(), data = dd); "fit"}, error = function(e) conditionMessage(e)), "\n")
cc <- frm_compat(); cat("compat names bhaz:", any(grepl("bhaz", c(cc$feature_a, cc$feature_b))), "\n")
