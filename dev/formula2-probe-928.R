# standata:928's own case: y, x1 and x2 are the same column 1:10
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
dat <- data.frame(y = 1:10, x1 = 1:10, x2 = 1:10)
fr <- frm(y ~ ., dat, dry_run = "frame")
print(deparse1(fr$spec$responses$y$dpars$mu$fixed))
print(colnames(fr$linpreds[["y.mu"]]$X))
print(fr$linpreds[["y.mu"]]$dropped_colnames)
fr2 <- frm(y ~ x1 + x2, dat, dry_run = "frame")
print(colnames(fr2$linpreds[["y.mu"]]$X))
